#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 14:15:43 2020

@author: cjburke
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from gwcs.wcstools import wcs_from_points
from gwcs.wcs import WCS as WCS_gwcs
from gwcs.coordinate_frames import *
import os
import glob
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join( DIR,'..'))
from btjd.btjd_correction import btjd_correction
#not happy about the global variables
ephemeris_file = os.path.join( DIR,'../btjd/tess_ephem.txt')
ephemeris_data = np.genfromtxt(os.path.join(ephemeris_file))

from wcs_build.step1_get_refimg_ctrlpts import gdPRF_calc, idx_filter, ring_background

import photutils.centroids as cent

from tess_stars2px import tess_stars2px_reverse_function_entry

from astropy.modeling.core import Model
from astropy.modeling import projections
from astropy.modeling import models, fitting
from astropy import coordinates as coord

import warnings
import argparse

from statsmodels import robust


def add_btjd_info(time, ra, dec, ephem_data):
    header = fits.Header()
    #print(time, ra,dec, ephem_data[-10:, 0] )
    btjd_mid, tess_position = btjd_correction(time, ra, dec, ephem_data)
    #print(btjd_mid, tess_position)
    #this does not yet seem reliable---I think the equation and
    #geometry I'm using is wrong, can't easily reproduce spoc
    #leaving out htis keyword, for now
    #header['BTJDMID'] = (btjd_mid, 'Barycentric TJD, mid-exposure, center of frame')

    #do I know that I got X and Y right?  Possible they were mixed up?
    header['TESS_X']  = (float('{:>13.2f}'.format(tess_position[0])), 'Spacecraft X coord, J2000 (km from barycenter)')
    header['TESS_Y']  = (float('{:>13.2f}'.format(tess_position[1])), 'Spacecraft Y coord, J2000 (km from barycenter)')
    header['TESS_Z']  = (float('{:>13.2f}'.format(tess_position[2])), 'Spacecraft Z coord, J2000 (km from barycenter)')
    header['TESS_VX'] = (float('{:>6.2f}'.format(tess_position[3])), 'Spacecraft X velocity (km/s)')
    header['TESS_VY'] = (float('{:>6.2f}'.format(tess_position[4])), 'Spacecraft Y velocity (km/s)')
    header['TESS_VZ'] = (float('{:>6.2f}'.format(tess_position[5])), 'Spacecraft Z velocity (km/s)')

    return header

def binmedian(xdata, ydata, nBins=30, xmin=None, xmax=None, showDetail=False):
    if xmin == None:
        xmin = xdata.min()
    if xmax == None:
        xmax = xdata.max()
    xedges = np.linspace(xmin, xmax, nBins+1)
    midx = xedges[:-1] + np.diff(xedges)/2.0
    iargs = np.digitize(xdata, xedges)
    medata = np.zeros_like(midx)
    mndata = np.zeros_like(midx)
    stddata = np.zeros_like(midx)
    maddata = np.zeros_like(midx)
    ndata = np.zeros_like(midx)
    for i in np.arange(0,nBins):
        iuse = np.where(iargs == i+1)[0]
        if len(iuse) >= 3:
            medata[i] = np.median(ydata[iuse])
            mndata[i] = np.mean(ydata[iuse])
            stddata[i] = np.std(ydata[iuse])
            maddata[i] = robust.mad(ydata[iuse])
            ndata[i] = len(ydata[iuse])
        elif len(iuse) > 0 and len(iuse) < 3:
            medata[i] = np.median(ydata[iuse])
            mndata[i] = np.mean(ydata[iuse])
            stddata[i] = 0.0
            maddata[i] = 0.0
            ndata[i] = len(ydata[iuse])
        elif len(iuse) == 0:
            medata[i] = 0.0
            mndata[i] = 0.0
            stddata[i] = 0.0
            maddata[i] = 0.0
            ndata[i] = 0            
        
    if showDetail:
        for i in np.arange(0,nBins):
            errmn = stddata[i]/np.sqrt(ndata[i])
            sigmn = mndata[i] / errmn
            print('i: {0:d} med: {1:f} mn: {2:f} n: {3:f} errmn: {4:f} sigdif: {5:f} midx: {6:f}'.format(\
                  i, medata[i], mndata[i], ndata[i], errmn, sigmn, midx[i]))
        
        
        
    return medata, midx, mndata, stddata, maddata, iargs, ndata

# NOTE*** CJB 2021-07 Recommend removing inverse_wcs_from_points function as it is no longer
#  correct
def inverse_wcs_from_points(xy, world_coordinates, fiducial, degree=4):
    x, y = xy
    lon, lat = world_coordinates
    skyrot = models.RotateCelestial2Native(fiducial.data.lon, fiducial.data.lat, 180.0*u.deg)
    projection=projections.Sky2Pix_TAN()
    trans = (skyrot | projection)
    projection_x, projection_y = trans(lon, lat)
    poly = models.Polynomial2D(degree)
    fitter = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        poly_proj_x = fitter(poly, projection_x, projection_y, x)
        poly_proj_y = fitter(poly, projection_x, projection_y, y)
    transform = skyrot | projection | models.Mapping((0, 1, 0, 1)) | poly_proj_x & poly_proj_y
    skyframe = CelestialFrame(reference_frame=fiducial.frame)
    detector = Frame2D(name="detector")
    pipeline = [(skyframe, transform), \
                (detector, None)]
    return WCS_gwcs(pipeline)

def inverse_wcs_transform(hdr, ras, decs):
    # does the inverse transform of going from ra&dec -> col, row pixel coordinates
    #  uses the header wcs coefficients directly rather than numerical inversion.
    fitDegree = hdr['AP_ORDER']
    
    # Get the cd matrix
    # Convert columns and rows into the 
    raproj = hdr['CRVAL1']
    decproj = hdr['CRVAL2']
    c11 = hdr['CD1_1'] 
    c12 = hdr['CD1_2'] 
    c21 = hdr['CD2_1'] 
    c22 = hdr['CD2_2'] 
    proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
    # Setup a transformation pipeline that will take ra and dec
    #  do a sky rotation to be centered on CRVALs
    # do a TAN projection to the 'intermediate world coordinates' of Shupe
    # finally do matrix multiply of inverse CDmatrix to get to 
    # U,V of corrected pixel coordinates
    skyrot = models.RotateCelestial2Native(proj_point.data.lon, proj_point.data.lat, 180.0*u.deg)
    projection=projections.Sky2Pix_TAN()
    # Get the inverse of the cd matrix from wcs
    cdmat = np.array([[c11,c12],[c21,c22]])
    inv_cdmat = np.linalg.inv(cdmat)
    matrix_tran = projections.AffineTransformation2D(matrix=inv_cdmat)
    trans = (skyrot | projection | matrix_tran)
    projection_U, projection_V = trans(ras, decs)
    # now apply the SIP distortion polynomial
    cols = np.zeros_like(ras)
    rows = np.zeros_like(decs)
    cols = projection_U
    rows = projection_V
#    cnt = 0
    for i in range(1,fitDegree+1):
        j = 0
        while j+i <=fitDegree:
            Ap = hdr['AP_{0:d}_{1:d}'.format(i,j)]
            cols = cols + Ap * np.power(projection_U,i)*np.power(projection_V,j)
            Bp = hdr['BP_{0:d}_{1:d}'.format(i,j)]
            rows = rows + Bp * np.power(projection_U,i)*np.power(projection_V,j)
#            print(i,j)
            j = j + 1
#            cnt = cnt + 1
    for j in range(1,fitDegree+1):
        i = 0
#        while (i<j) and (j+i<=fitDegree) :
        Ap = hdr['AP_{0:d}_{1:d}'.format(i,j)]
        cols = cols + Ap * np.power(projection_U,i) * np.power(projection_V,j)
        Bp = hdr['BP_{0:d}_{1:d}'.format(i,j)]
        rows = rows + Bp * np.power(projection_U,i) * np.power(projection_V,j)
#        print(i,j)
        i = i + 1
#        cnt = cnt + 1
#    print(cnt)
    return cols, rows

    
def fit_wcs(ras, decs, cols, rows, tmags, \
            SECTOR_WANT,CAMERA_WANT, CCD_WANT,
            blkidxs=None, NCOL=None, \
            fitDegree=5, noClip=False, DEBUG_LEVEL=0, MAKE_FIGS=None):
    # Use gwcs to fit a WCS and look at residuals if debugging
    np.random.seed(1010101)

    REFPIXCOL = 1024.0+45.0
    REFPIXROW = 1024.0
    PIX2DEG = 21.0/3600.0 # Turn pixels to degrees roughly
        # The WCS fitter uses WCS intermediate coordinates which are in degrees
        # Relative to reference ra and dec as th eprojection zero point.
        # Trim out control points if they are offset from 
    #  wcs fit by OUTLIER_TRIM [arcsec]
    OUTLIER_TRIM = 25.0
    OUTLIER_SIGMA = 4.0
    if noClip:
        OUTLIER_SIGMA = 9.99e99
    # Record wcs fit residuals for subregions
    # Doing a mapping where targets in several subregions in the corners are
    # combined to get teh fit residual
    nResids = 4
    residSZ = 3
    #  With 5x5 grid these are corner subregions
    useRR = [[0,1,5], [3,4,9], [15,20,21], [23,24,19]]
    exResids = np.zeros((nResids,), dtype=np.float)
    
    # From reference pixel coordinates get the estimated ra and dec of this point
    raproj, decproj, scinfo = tess_stars2px_reverse_function_entry(\
                         SECTOR_WANT, CAMERA_WANT, CCD_WANT, REFPIXCOL, REFPIXROW)
    proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
    #  Reference subtracted pixel coordinates
    sclObsCols = (cols - REFPIXCOL) * PIX2DEG
    sclObsRows = (rows - REFPIXROW) * PIX2DEG
        
    xy = (sclObsCols, sclObsRows)
    radec = (ras, decs)
    gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
    gdResids = np.ones_like(cols, dtype=bool)
    # Look for outliers to trim
    gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
    deg2Rad = np.pi/180.0
    deltaRas = (gwcs_pred_ras - ras) *3600.0 * np.cos(decs*deg2Rad)
    deltaDecs = (gwcs_pred_decs - decs) *3600.0
    deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
    # First trim out ones that deviate by absolute amount
    idxbd = np.where(deltaSeps >= OUTLIER_TRIM)[0]
    gdResids[idxbd] = False
    # Next trim out ones that deviate by sigma based on tmag residual
    #  I use tmag residual in order to not remove stars based on 
    #  position that maybe didn't fit well the first time
    meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(tmags, \
                                        deltaRas)
    unarg = np.sort(np.unique(iargs))
    for curi in unarg[:-1]:
        idxcur = np.where(iargs == curi)[0]
        cury = np.abs(deltaRas[idxcur]-mndata[curi-1])
        idxclp = np.where(cury >= OUTLIER_SIGMA*maddata[curi-1])[0]
        #print('{0:d} badn: {1:d}'.format(curi, len(idxclp)))
        gdResids[idxcur[idxclp]] = False
    meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(tmags, \
                                        deltaDecs)        
    unarg = np.sort(np.unique(iargs))
    for curi in unarg[:-1]:
        idxcur = np.where(iargs == curi)[0]
        cury = np.abs(deltaDecs[idxcur]-mndata[curi-1])
        idxclp = np.where(cury >= OUTLIER_SIGMA*maddata[curi-1])[0]
        #print('{0:d} badn: {1:d}'.format(curi, len(idxclp)))
        gdResids[idxcur[idxclp]] = False
    
    idxgd = np.where(gdResids)[0]
    idxbd = np.where(np.logical_not(gdResids))[0]
    # refit
    xy = (sclObsCols[idxgd], sclObsRows[idxgd])
    radec = (ras[idxgd], decs[idxgd])
    gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
        
    # Iterate to find a better reference ra and dec
    #  such that there is no constant term in the polynomial fit
    newrefra, newrefdec = gwcs_obj(0.0, 0.0)
    if DEBUG_LEVEL>0:
        diffra = (raproj-newrefra)*3600.0
        diffdec = (decproj-newrefdec)*3600.0
        print('Old RA: {0:9.5f} NewRA: {1:9.5f} DelRa: {2:9.4f} OldDec: {3:9.5f} NewDec: {4:9.5f} DelDec: {5:9.4f}'.format(\
              raproj, newrefra, np.sqrt(diffra*diffra), decproj, newrefdec, np.sqrt(diffdec*diffdec)))
    raproj = newrefra
    decproj = newrefdec
    proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
    gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
    
    # We have the wcs in polynomial form
    # now we need to convert it to SIP form
    # Start a header
    hdr = fits.Header()
    hdr['NAXIS'] = 2
    hdr['CTYPE1'] = 'RA---TAN-SIP'
    hdr['CTYPE2'] = 'DEC--TAN-SIP'
    hdr['CRPIX1'] = REFPIXCOL
    hdr['CRPIX2'] = REFPIXROW
    hdr['CRVAL1'] = raproj
    hdr['CRVAL2'] = decproj
    hdr['A_ORDER'] = fitDegree
    hdr['B_ORDER'] = fitDegree
    hdr['AP_ORDER'] = fitDegree
    hdr['BP_ORDER'] = fitDegree
        
    # Get the transformation matrix terms
    # polynomial fit object in x (col)
    cx = gwcs_obj.forward_transform[1]
    # polynomial fit object in y (row)
    cy = gwcs_obj.forward_transform[2]
    c11 = cx.c1_0.value
    hdr['CD1_1'] = c11 * PIX2DEG
    c12 = cx.c0_1.value
    hdr['CD1_2'] = c12 * PIX2DEG
    c21 = cy.c1_0.value
    hdr['CD2_1'] = c21 * PIX2DEG
    c22 = cy.c0_1.value
    hdr['CD2_2'] = c22 * PIX2DEG
    b = c12/c11
    atmp = 1.0/c11
    c = c21/c22
    dtmp = 1.0/c22
    bc = b*c
    onembc = 1.0 - bc
    if np.isclose(onembc, 0.0, rtol=1e-10, atol=1e-12):
        print('Warning SIP polynomial conversion is ill behaved')
    ac = atmp*c
    ACoeffs = np.zeros((fitDegree+1,fitDegree+1), dtype=np.double)
    BCoeffs = np.zeros((fitDegree+1,fitDegree+1), dtype=np.double)
    # NOTE CJB in hindsight for loop over i and loop over j below
    #  sets some of the terms multiple times. For example
    #  A_1_1 is set multipletimes. This does no harm, but
    #  being inefficient and maybe confusing for someone debugging
    # down the line.
    for i in range(2,fitDegree+1):
        j = 0
        while j+i <=fitDegree:
            coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
            coeffxVal = getattr(coeffxObj, 'value');
            coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
            coeffyVal = getattr(coeffyObj, 'value');
            cura = coeffxVal * atmp
            curd = coeffyVal * dtmp
            curac = coeffxVal * ac
            BCoeffs[i,j] = (curd - curac)/ onembc 
            ACoeffs[i,j] = (cura - b*BCoeffs[i,j]) 
            hdr['A_{0:d}_{1:d}'.format(i,j)] = ACoeffs[i,j]* np.power(PIX2DEG, i+j-1)
            hdr['B_{0:d}_{1:d}'.format(i,j)] = BCoeffs[i,j]* np.power(PIX2DEG, i+j-1)
            j = j +1
    for j in range(2,fitDegree+1):
        i = 0
        while j+i <=fitDegree:
            coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
            coeffxVal = getattr(coeffxObj, 'value');
            coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
            coeffyVal = getattr(coeffyObj, 'value');
            cura = coeffxVal * atmp
            curd = coeffyVal * dtmp
            curac = coeffxVal * ac
            BCoeffs[i,j] = (curd - curac)/ onembc
            ACoeffs[i,j] = (cura - b*BCoeffs[i,j])          
            hdr['A_{0:d}_{1:d}'.format(i,j)] = ACoeffs[i,j]* np.power(PIX2DEG, i+j-1)
            hdr['B_{0:d}_{1:d}'.format(i,j)] = BCoeffs[i,j]* np.power(PIX2DEG, i+j-1)
            i = i +1
    i = 1
    j = 1
    coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
    coeffxVal = getattr(coeffxObj, 'value');
    coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
    coeffyVal = getattr(coeffyObj, 'value');
    cura = coeffxVal * atmp
    curd = coeffyVal * dtmp
    curac = coeffxVal * ac
    BCoeffs[i,j] = (curd - curac)/ onembc 
    ACoeffs[i,j] = (cura - b*BCoeffs[i,j])         
    hdr['A_{0:d}_{1:d}'.format(i,j)] = ACoeffs[i,j] * np.power(PIX2DEG, 1)
    hdr['B_{0:d}_{1:d}'.format(i,j)] = BCoeffs[i,j] * np.power(PIX2DEG, 1)

    # Do comparison to wcs SIP interpretation of fit
    my_wcs = WCS(hdr)
    pix_list = list(zip(cols[idxgd], rows[idxgd]))
    pix_coords = np.array(pix_list, dtype=np.double)
    pred_coords = my_wcs.all_pix2world(pix_coords, 1)
    predRa = np.array([x[0] for x in pred_coords], dtype=np.double)
    predDec = np.array([x[1] for x in pred_coords], dtype=np.double)
    deg2Rad = np.pi/180.0
    deltaRas = (predRa - ras[idxgd]) *3600.0 * np.cos(predDec*deg2Rad)
    deltaDecs = (predDec - decs[idxgd]) *3600.0
    deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
    if len(idxbd)>3:
        pix_listbd = list(zip(cols[idxbd], rows[idxbd]))
        pix_coordsbd = np.array(pix_listbd, dtype=np.double)
        pred_coordsbd = my_wcs.all_pix2world(pix_coordsbd, 1)
        predRabd = np.array([x[0] for x in pred_coordsbd], dtype=np.double)
        predDecbd = np.array([x[1] for x in pred_coordsbd], dtype=np.double)
        deg2Rad = np.pi/180.0
        deltaRasbd = (predRabd - ras[idxbd]) *3600.0 * np.cos(predDecbd*deg2Rad)
        deltaDecsbd = (predDecbd - decs[idxbd]) *3600.0
    c1 = 1.0/np.sqrt(2.0)
    idx = np.where(tmags[idxgd]<10.0)[0]
    std1 = np.std(deltaRas[idx])
    std2 = np.std(deltaDecs[idx])
    brightstd = np.sqrt(std1*std1+std2*std2)*c1
    idx = np.where(tmags[idxgd]>10.0)[0]
    # Protect agains too few targets
    if len(idx) < 5:
        std1 = 0.0
        std2 = 0.0
    else:
        std1 = np.std(deltaRas[idx])
        std2 = np.std(deltaDecs[idx])
    # protect against nans
    if not np.isfinite(std1):
        std1 = 0.0
    if not np.isfinite(std2):
        std2 = 0.0
    faintstd = np.sqrt(std1*std1+std2*std2)*c1
    std1 = np.std(deltaRas)
    std2 = np.std(deltaDecs)
    allstd = np.sqrt(std1*std1+std2*std2)*c1

    # Now its time to do the inverse
    #  We will following gwcs wcs_from_points function to do inverse
    #  however, the inverse needs to use the same reference points
    #  and inverse of the CDmatrix, so parameters need to be held fixed
    c11 = hdr['CD1_1'] 
    c12 = hdr['CD1_2'] 
    c21 = hdr['CD2_1'] 
    c22 = hdr['CD2_2'] 
    
    x = (cols[idxgd] - REFPIXCOL)
    y = (rows[idxgd] - REFPIXROW)
    lon, lat = radec
    skyrot = models.RotateCelestial2Native(proj_point.data.lon, proj_point.data.lat, 180.0*u.deg)
    projection=projections.Sky2Pix_TAN()
    # invert the cd matrix
    cdmat = np.array([[c11,c12],[c21,c22]])
    inv_cdmat = np.linalg.inv(cdmat)
    matrix_tran = projections.AffineTransformation2D(matrix=inv_cdmat)

    trans = (skyrot | projection | matrix_tran)
    projection_x, projection_y = trans(lon, lat)
    # Here is where we need to fix parameters for the fitting polynomial
    polyx = models.Polynomial2D(fitDegree, c0_0=0.0)
    polyx.c0_0.fixed = True
    polyy = models.Polynomial2D(fitDegree, c0_0=0.0)
    polyy.c0_0.fixed = True
    fitter = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        poly_proj_x = fitter(polyx, projection_x, projection_y, x)
        poly_proj_y = fitter(polyy, projection_x, projection_y, y)
    transform = skyrot | projection | matrix_tran | models.Mapping((0, 1, 0, 1)) | poly_proj_x & poly_proj_y
    skyframe = CelestialFrame(reference_frame=proj_point.frame)
    detector = Frame2D(name="detector")
    pipeline = [(skyframe, transform), \
                (detector, None)]
    inv_gwcs_obj = WCS_gwcs(pipeline)

    # Convert the polynomial coefficients in to the SIP wcs format
    # polynomial fit objects
    cx = inv_gwcs_obj.forward_transform[4]
    cy = inv_gwcs_obj.forward_transform[5]
    # NOTE CJB 2021-07 again like above in the loop over i and j
    # there are terms being
    # set multiple times. This does no harm other than  being
    # inefficient and confusing for future examiners of code.
    for i in range(2,fitDegree+1):
        j = 0
        while j+i <=fitDegree:
            coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
            coeffxVal = getattr(coeffxObj, 'value');
            coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
            coeffyVal = getattr(coeffyObj, 'value');
            hdr['AP_{0:d}_{1:d}'.format(i,j)] = coeffxVal
            hdr['BP_{0:d}_{1:d}'.format(i,j)] = coeffyVal
            j = j +1
    for j in range(2,fitDegree+1):
        i = 0
        while j+i <=fitDegree:
            coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
            coeffxVal = getattr(coeffxObj, 'value');
            coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
            coeffyVal = getattr(coeffyObj, 'value');
            hdr['AP_{0:d}_{1:d}'.format(i,j)] = coeffxVal
            hdr['BP_{0:d}_{1:d}'.format(i,j)] = coeffyVal
            i = i +1
    i = 1
    j = 1
    coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
    coeffxVal = getattr(coeffxObj, 'value');
    coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
    coeffyVal = getattr(coeffyObj, 'value');      
    hdr['AP_{0:d}_{1:d}'.format(i,j)] = coeffxVal
    hdr['BP_{0:d}_{1:d}'.format(i,j)] = coeffyVal
    i = 1
    j = 0
    coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
    coeffxVal = getattr(coeffxObj, 'value');
    coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
    coeffyVal = getattr(coeffyObj, 'value');      
    hdr['AP_{0:d}_{1:d}'.format(i,j)] = coeffxVal - 1.0
    hdr['BP_{0:d}_{1:d}'.format(i,j)] = coeffyVal
    i = 0
    j = 1
    coeffxObj = getattr(cx, 'c{0:d}_{1:d}'.format(i,j))
    coeffxVal = getattr(coeffxObj, 'value');
    coeffyObj = getattr(cy, 'c{0:d}_{1:d}'.format(i,j))
    coeffyVal = getattr(coeffyObj, 'value');      
    hdr['AP_{0:d}_{1:d}'.format(i,j)] = coeffxVal
    hdr['BP_{0:d}_{1:d}'.format(i,j)] = coeffyVal - 1.0

    # Do comparison to wcs SIP interpretation of fit
    my_pred_cols, my_pred_rows = inverse_wcs_transform(hdr, ras[idxgd], decs[idxgd])
    deltaCols = (my_pred_cols + REFPIXCOL - cols[idxgd])
    deltaRows = (my_pred_rows + REFPIXROW - rows[idxgd])
    deltaPixSeps = np.sqrt(deltaCols*deltaCols + deltaRows*deltaRows)
    idx = np.where(tmags[idxgd]<10.0)[0]
    std1 = np.std(deltaCols[idx])
    std2 = np.std(deltaRows[idx])
    brightstdpix = np.sqrt(std1*std1+std2*std2)*c1
    idx = np.where(tmags[idxgd]>10.0)[0]
    # Protect agains too few targets
    if len(idx) < 5:
        std1 = 0.0
        std2 = 0.0
    else:
        std1 = np.std(deltaCols[idx])
        std2 = np.std(deltaRows[idx])
    # protect against nans
    if not np.isfinite(std1):
        std1 = 0.0
    if not np.isfinite(std2):
        std2 = 0.0
    faintstdpix = np.sqrt(std1*std1+std2*std2)*c1
    std1 = np.std(deltaCols)
    std2 = np.std(deltaRows)
    allstdpix = np.sqrt(std1*std1+std2*std2)*c1
    
    
    # Also get residualsfrom targets over several regions
    tm = tmags[idxgd]
    blks = blkidxs[idxgd]
    for i in range(nResids):
        idxAll = np.array([], dtype=np.int64)
        for j in range(residSZ):
            idx = np.where((tm<10.0) & (blks == useRR[i][j]))[0]
            idxAll = np.append(idxAll, idx)
        if len(idxAll) < 5:
            mad1 = 0.0
            mad2 = 0.0
        else:
            mad1 = robust.mad(deltaRas[idxAll])
            mad2 = robust.mad(deltaDecs[idxAll])
        if np.isfinite(mad1) and np.isfinite(mad2):
            exResids[i]= np.sqrt(mad1*mad1+mad2*mad2)*c1
        else:
            exResids[i] = 0.0
    
    #if DEBUG_LEVEL>0:
    print('world_coord resid: {0:f} {1:f} {2:f}'.format(brightstd, allstd, faintstd))
    print('Pixel_coord resid:{0:f} {1:f} {2:f}'.format(brightstdpix, allstdpix, faintstdpix))
                
    if DEBUG_LEVEL>2 or (MAKE_FIGS and MAKE_FIGS.strip()):
        showDetailBool=False
        if DEBUG_LEVEL>3:
            showDetailBool=True
        # Do comparison to gwcs fit
        #gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
        #deltaRas = (gwcs_pred_ras - ras) *3600.0 * np.cos(gwcs_pred_decs*deg2Rad)
        #deltaDecs = (gwcs_pred_decs - decs)* 3600.0
        # Check gwcs inerse
        #gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
        #gwcs_pred_cols, gwcs_pred_rows = inv_gwcs_obj(gwcs_pred_ras, gwcs_pred_decs)
        #plt.plot(sclObsCols, gwcs_pred_cols - sclObsCols, '.')
        #plt.show()
        #plt.plot(sclObsRows, gwcs_pred_rows - sclObsRows, '.')
        #plt.show()
        idx = np.where(tmags[idxgd] < 16.0)[0]
        useDras = deltaRas[idx]
        useDdecs = deltaDecs[idx]
        useblkidxs = blkidxs[idx]
        uniqblk = np.unique(blkidxs)
        nblk = len(uniqblk)
        # xtmp = np.array([], dtype=np.int)
        # ytmp = np.array([], dtype=np.int)
        # delRa = np.array([], dtype=np.float)
        # delDec = np.array([], dtype=np.float)
        # for i in range(nblk):
        #     ii = i//NCOL
        #     jj = np.mod(i, NCOL)
        #     xtmp = np.append(xtmp, ii)
        #     ytmp = np.append(ytmp, jj)
        #     idx = np.where(useblkidxs == i)[0]
        #     delRa = np.append(delRa, np.mean(useDras[idx]))
        #     delDec = np.append(delDec, np.mean(useDdecs[idx]))
        # plt.quiver(xtmp, ytmp, delRa, delDec)
        # plt.xlabel('Detector Column Position; Arrow Mean RA Offset Angle on Sky [arcsec]')
        # plt.ylabel('Detector Row Position; Arrow Mean Dec Offset [arcsec]')
        # plt.ylim([-2.0, 6.0])
        # plt.xlim([-2.0, 6.0])
        # plt.show()

        plt.plot(tmags[idxgd], deltaRas, '.')
        if len(idxbd)>3:
            plt.plot(tmags[idxbd], deltaRasbd, '.k')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(tmags[idxgd], deltaRas, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Tmag')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-21.0, 21.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            print(MAKE_FIGS)
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_TmagVRaDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
            
        plt.plot(tmags[idxgd], deltaDecs, '.')
        if len(idxbd)>3:
            plt.plot(tmags[idxbd], deltaDecsbd, '.k')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(tmags[idxgd], deltaDecs, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Tmag')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-21.0, 21.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_TmagVDecDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()

        plt.plot(cols[idxgd], deltaRas, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(cols[idxgd], deltaRas, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Column [px]')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-21.0, 21.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_ColVRaDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        

        plt.plot(rows[idxgd], deltaRas, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(rows[idxgd], deltaRas, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Row [px]')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-21.0, 21.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_RowVRaDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(cols[idxgd], deltaDecs, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(cols[idxgd], deltaDecs, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Column [px]')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-21.0, 21.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_ColVDecDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(rows[idxgd], deltaDecs, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(rows[idxgd], deltaDecs, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Row [px]')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-21.0, 21.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_RowVDecDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        # now do inverse data
        plt.plot(tmags[idxgd], deltaCols, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(tmags[idxgd], deltaCols, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Tmag')
        plt.ylabel('WCS Predicted - Observed Column Position [pix]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-1.0, 1.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_TmagVColDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(tmags[idxgd], deltaRows, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(tmags[idxgd], deltaRows, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Tmag')
        plt.ylabel('WCS Predicted - Observed Row Position [pix]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-1.0, 1.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_TmagVRowDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(cols[idxgd], deltaCols, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(cols[idxgd], deltaCols, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Column [px]')
        plt.ylabel('WCS Predicted - Observed Column Position [pix]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-1.0, 1.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_ColVColDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(rows[idxgd], deltaCols, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(rows[idxgd], deltaCols, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Row [px]')
        plt.ylabel('WCS Predicted - Observed Column Position [pix]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-1.0, 1.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_RowVColDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(cols[idxgd], deltaRows, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(cols[idxgd], deltaRows, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Column [px]')
        plt.ylabel('WCS Predicted - Observed Row Position [pix]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-1.0, 1.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_ColVRowDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
        
        plt.plot(rows[idxgd], deltaRows, '.')
        meddata, midx, mndata, stddata, maddata, iargs, ndata = binmedian(rows[idxgd], deltaRows, showDetail=showDetailBool)
        plt.plot(midx, meddata, '-')
        plt.xlabel('Row [px]')
        plt.ylabel('WCS Predicted - Observed Row Position [pix]')
        plt.axhline(0.0, ls='--', color='r')
        #plt.ylim([-1.0, 1.0])
        if MAKE_FIGS and MAKE_FIGS.strip():
            plt.title(MAKE_FIGS)
            plt.savefig('{0}_RowVRowDiff.png'.format(MAKE_FIGS), dpi=300)
        if DEBUG_LEVEL <= 2:
            plt.close()
        else:
            plt.show()
    
    return hdr, allstd, brightstd, faintstd, allstdpix, brightstdpix, faintstdpix, gdResids, exResids

def fit_wcs_in_imgdir(SECTOR_WANT, CAMERA_WANT, CCD_WANT, REF_DATA, \
                    IMG_LIST_STR, outputDir, fitDegree=5, saveDiag=False,\
                    DEBUG_LEVEL=0):
    
    # ***These following parameters must match what was used in 
    #  step1 to generate the reference image control points
    # If you change these then you will want to rerun step1 
    #  to keep things consistent
    pixScl = 21.0 # arcsec
    colMin = 45
    colMax = 2092
    rowMin = 1
    rowMax = 2048
    # We want control point grid over pixels to get an even
    # distribution of stars for determining wcs
    #  These control how many regions in each dimension to find stars
    CTRL_PER_COL = 5
    CTRL_PER_ROW = 5
    # number of pixels above and below midpoint to determine centroid
    blkHlf = 5
    blkHlfCent = 2 # The analysis region for background estimate, gdPRF check, initial
                    # position offets use blkHlf. However, final centroid is calcualted
                    # with a small blkHlfCent region
    FIGOUTEVERY = 300 # Write out fit residual figures for every
                        # FIGOUTEVERY image
    
    # Load the reference position information
    fin = h5py.File(REF_DATA, 'r')    
    tics = np.array(fin['tics'])
    ras = np.array(fin['ras'])
    decs = np.array(fin['decs'])
    tmags = np.array(fin['tmags'])
    blkidxs = np.array(fin['blkidxs'])
    obscols = np.array(fin['obscols'])
    obsrows = np.array(fin['obsrows'])
    # Make sure these are sorted by tmag
    idx = np.argsort(tmags)
    tics, ras, decs, tmags, blkidxs, obscols, obsrows = idx_filter(idx, \
                tics, ras, decs, tmags, blkidxs, obscols, obsrows)


    # The following sets up the analysis subregions there are CTRL_PER_COL x CTRL_PER_Row 
    #  subregions
    spanCol = colMax - colMin
    spanRow = rowMax - rowMin
    delCol = float(spanCol) / float(CTRL_PER_COL)
    delRow = float(spanRow)/ float(CTRL_PER_ROW)
    # These are the center of each sub region for analysis
    colCtrl1D = np.linspace(delCol/2.0 + colMin, colMax - delCol/2.0, CTRL_PER_COL)
    rowCtrl1D = np.linspace(delRow/2.0 + rowMin, rowMax - delRow/2.0, CTRL_PER_ROW)
    colCtrl2D, rowCtrl2D = np.meshgrid(colCtrl1D, rowCtrl1D)    
    colCtrl2D_flat = colCtrl2D.flatten()
    rowCtrl2D_flat = rowCtrl2D.flatten()
    nCtrl = len(colCtrl2D_flat)

    # Make list of images in work dir for this camera and ccd
    inputImgList = IMG_LIST_STR#glob.glob(IMG_LIST_STR)
    # sort images
    inputImgList.sort()
#    print(inputImgList)

    # Always initialize to the reference image coordinates
    lstGdCols = np.copy(obscols)
    lstGdRows = np.copy(obsrows)
    
    # ***TESTTESTTEST
    # If you want to test on a subsample 
    #  of the image list use this line
    #inputImgList = inputImgList[0::20]
    #inputImgList = inputImgList[-13:-10]
    #inputImgList = inputImgList[3532:3538]
    # BE SURE TO COMMENT OUT IN PRODUCTION

    # The following arrays will keep the 
    #  the diagnostic information about quality of wcs fit    
    nImg = len(inputImgList)
    imgNames = np.array([], dtype=np.str)
    allStds = np.zeros((nImg,), dtype=np.float)
    brightStds = np.zeros_like(allStds)
    faintStds = np.zeros_like(allStds)
    gdFracs = np.zeros_like(allStds)
    allPixStds = np.zeros_like(allStds)
    brightPixStds = np.zeros_like(allStds)
    faintPixStds = np.zeros_like(allStds)
    ts = np.zeros_like(allStds)
    exStd0s = np.zeros_like(allStds)
    exStd1s = np.zeros_like(allStds)
    exStd2s = np.zeros_like(allStds)
    exStd3s = np.zeros_like(allStds)
    
    # get timestamp header keyword
    #  and determine whether this is QLP/TICA or SPOC data
    gotTimeStamp = False
    iImg = 0
    # Main loop over images 
    for curImg in inputImgList:
        #if DEBUG_LEVEL>0:           
#        print('{0:d} of {1:d} {2}'.format(iImg, nImg, curImg))
        # open image
        hdulistCal = fits.open(curImg)
        if len(hdulistCal) == 2:
            print('skipping,  already has WCS'.format(curImg))
            continue
        imgNames = np.append(imgNames, os.path.basename(curImg))
        if not gotTimeStamp:
            try:
                # QLP/TICA has MIDTJD
                timeVal = hdulistCal[0].header['MIDTJD']
                timeKey = 'MIDTJD'
                dataKey = 0
            except:
                # SPOC has TSTART
                timeKey = 'TSTART'
                dataKey = 1
            gotTimeStamp = True
        if 'cal' in curImg:
            dataKey = 0
        # Always initialize to the reference image coordinates
        lstGdCols = np.copy(obscols)
        lstGdRows = np.copy(obsrows)

        newCols = np.zeros_like(obscols)
        newRows = np.zeros_like(obscols)
        gdPrfs = np.zeros_like(obscols, dtype=np.int)
        # these will be used to keep track of average subregion dc offsets
        nDC = 0
        dccolsum = 0.0
        dcrowsum = 0.0
        # Loop over subregions
        # keep track of regions that cannot find any valid stars
        isGdRegion = np.zeros((len(colCtrl2D_flat),), dtype=bool)
        for ii, curCol in enumerate(colCtrl2D_flat):
            # Get the reference data corresonding with this subregion
            ia = np.where(blkidxs == ii)[0]
            curLstGdCols = lstGdCols[ia]
            curLstGdRows = lstGdRows[ia]
            curGdPrfs = np.zeros_like(curLstGdRows, dtype=np.int)
            curNewCols = np.zeros_like(curLstGdRows)
            curNewRows = np.zeros_like(curLstGdRows)
            # In this subregion first get the brightest stars
            #  to determine an initial dc offset in col and row
            #  The targets were sorted by Tmag at the begining
            #  so no need to do it here
            nTry = 5
            gotN = 0
            jj = 0
            delCols = np.array([], dtype=np.double)
            delRows = np.array([], dtype=np.double)
            while (gotN <= nTry and jj<len(curLstGdCols)):
                midCol = int(round(curLstGdCols[jj]))
                midRow = int(round(curLstGdRows[jj]))
                colX = np.arange(midCol-blkHlf, midCol+blkHlf+1)
                rowY = np.arange(midRow-blkHlf, midRow+blkHlf+1)
#                print(midRow-blkHlf-1,  midRow+blkHlf, midCol-blkHlf-1, midCol+blkHlf)
#                print(dataKey, len(hdulistCal))
#                print('type',type(hdulistCal[dataKey].data))
#                print('shape',np.shape(hdulistCal[dataKey].data))
                sciImgCal = hdulistCal[dataKey].data[midRow-blkHlf-1: midRow+blkHlf, midCol-blkHlf-1: midCol+blkHlf]
                # Check if target is isolated
                try:
                    gdPRF, contrastCol, contrastRow = gdPRF_calc(sciImgCal, blkHlf)
                except ValueError:
                    break
                    
                if gdPRF:
                    # Determine background from 2 pixel ring around box
                    bkgLev = ring_background(sciImgCal)
                    centOut = cent.centroid_com(sciImgCal-bkgLev)
                    # centroid_com produces 0-based coordinates
                    newCol = centOut[0]+colX[0]
                    newRow = centOut[1]+rowY[0]
                    delCols = np.append(delCols, newCol - curLstGdCols[jj])
                    delRows = np.append(delRows, newRow - curLstGdRows[jj])
                    gotN = gotN + 1
                jj = jj+1
            if (gotN < nTry):
                if DEBUG_LEVEL>0:
                    print('Too few targets for initial pixel position correction on Block {0:d}'.format(ii))
                delCols = np.array([0.0, 0.0])
                delRows = np.array([0.0, 0.0])
            else:
                isGdRegion[ii] = True
            corCol = np.median(delCols)
            corRow = np.median(delRows)
            if isGdRegion[ii]:
                nDC = nDC + 1
                dccolsum = dccolsum + corCol
                dcrowsum = dcrowsum + corRow
            if DEBUG_LEVEL>0:
                print('Predicted Position Tweaks Col: {0:f} Row: {1:f} Blk: {2:d}'.format(corCol, corRow, ii))
            outColPix = curLstGdCols + corCol
            outRowPix = curLstGdRows + corRow

            # Now we are ready to get all target flux weighted centroids on this subimage
            for jj, curCol in enumerate(outColPix):
                if isGdRegion[ii]: # if region found offsets measure and update if good
                    # this is predicted coordinates
                    midCol = int(round(curCol))
                    midRow = int(round(outRowPix[jj]))
                    if midCol-blkHlf >= colMin and midCol+blkHlf <= colMax and midRow-blkHlf >= rowMin and midRow+blkHlf <= rowMax:
                        colX = np.arange(midCol-blkHlf, midCol+blkHlf+1)
                        rowY = np.arange(midRow-blkHlf, midRow+blkHlf+1)
                        sciImgCal = hdulistCal[dataKey].data[midRow-blkHlf-1: midRow+blkHlf, midCol-blkHlf-1: midCol+blkHlf]
                        gdPRF, contrastCol, contrastRow = gdPRF_calc(sciImgCal, blkHlf)
                    else:
                        gdPRF = False
                else: # if region did not find any offsets dont measure and done update
                    gdPRF = False

                if gdPRF:
                    # Determine background from 2 pixel ring around box
                    bkgLev = ring_background(sciImgCal)
                    # Use small blkHlfCent region for final centroid calc
                    colX = np.arange(midCol-blkHlfCent, midCol+blkHlfCent+1)
                    rowY = np.arange(midRow-blkHlfCent, midRow+blkHlfCent+1)
                    # analysis region around target
                    sciImgCal = hdulistCal[0].data[midRow-blkHlfCent-1: midRow+blkHlfCent, midCol-blkHlfCent-1: midCol+blkHlfCent]
                    centOut = cent.centroid_com(sciImgCal-bkgLev)
                    # Once in blue moon centroids go nutsy check for it
                    idx = np.where((centOut<-1.0) | (centOut>blkHlfCent*2+2))[0]
                    if len(idx) == 0:
                        # centroid_com produces 0-based coordinates
                        newCol = centOut[0]+colX[0]
                        newRow = centOut[1]+rowY[0]
                        curGdPrfs[jj] = 1
                        curNewCols[jj] = newCol
                        curNewRows[jj] = newRow
                    else:
                        curNewCols[jj] = outColPix[jj]
                        curNewRows[jj] = outRowPix[jj]
                else:
                    #  Insert last ref position  with the dc offsets
                    curNewCols[jj] = outColPix[jj]
                    curNewRows[jj] = outRowPix[jj]
                    if not isGdRegion[ii]:
                        curGdPrfs[jj] = -2
                    
            # Done getting flux weight for this sub block record them in overall list
            newCols[ia] = curNewCols
            newRows[ia] = curNewRows
            gdPrfs[ia] = curGdPrfs
        # Done with all the blocks for this image
        # Get the overall good prf fraction
        nGd = len(np.where(gdPrfs==1)[0])
        gdFracs[iImg] = float(nGd)/float(len(gdPrfs))
        #if DEBUG_LEVEL>0:
        print('GdFrac:{0:f}'.format(gdFracs[iImg]))
        #  Now we determine wcss
        # We need to reject bad prfs from wcs fit
        idxgd = np.where(gdPrfs == 1)[0]
        # In subregion with no good prfs we need to add in a couple 
        #  targets with reference coords just to keep fit  somewhat behaved
        # first determine average dc offset
        if nDC > 0:
            dccol = dccolsum/float(nDC)
            dcrow = dcrowsum/float(nDC)
        else:
            dccol = 0.0
            dcrow = 0.0
        idx_badregion = np.where(isGdRegion == 0)[0]
        for curibad in idx_badregion:
            idx = np.where(blkidxs == curibad)[0]
            # pick three random ones in this subregion
            idxUse = np.random.choice(idx, (3,), replace=False)
            gdPrfs[idxUse] = -1 # record these random added for plotting and diagnostics
            idxgd = np.append(idxgd, idxUse)
            # These targets never had dc offset added from subregion
            #  so add the average dc offset over all subregions
            newCols[idxUse] = newCols[idxUse] + dccol
            newRows[idxUse] = newRows[idxUse] + dcrow

        # If all subregions are bad then drop down to lower fit degree
        useFitDegree = fitDegree
        useNoClipping = False
        if np.sum(isGdRegion) == 0:
            useFitDegree = np.min([4, fitDegree])
            useNoClipping = True
            print('All regions bad, using {0:d} fit degree'.format(useFitDegree))

        # Write out FIGOUTEVERY image
        FIGOUTPREFIX = None
        if saveDiag == True:
            if np.mod(iImg, FIGOUTEVERY) == 0:
                fileBase = os.path.splitext(os.path.basename(curImg))[0]

                FIGOUTPREFIX= os.path.join(outputDir, 'wcs_diags2', fileBase)
        newhdr, allStd, brightStd, faintStd, allStdPix, \
                brightStdPix, faintStdPix, gdResids, exResids = fit_wcs(ras[idxgd], decs[idxgd], \
                                                newCols[idxgd], newRows[idxgd], \
                                                tmags[idxgd], \
                                                SECTOR_WANT, CAMERA_WANT, CCD_WANT,\
                                                blkidxs[idxgd], CTRL_PER_COL,\
                                                useFitDegree, useNoClipping, DEBUG_LEVEL,\
                                                FIGOUTPREFIX)
        # Add the outliers to gdPrfs
        idx = np.where(np.logical_not(gdResids))[0]
        gdPrfs[idxgd[idx]] = -1
        # save diagnostic fit quality for image
        allStds[iImg] = allStd
        brightStds[iImg] = brightStd
        faintStds[iImg] = faintStd
        allPixStds[iImg] = allStdPix
        brightPixStds[iImg] = brightStdPix
        faintPixStds[iImg] = faintStdPix
        ts[iImg] = hdulistCal[dataKey].header[timeKey]
        exStd0s[iImg] = exResids[0]
        exStd1s[iImg] = exResids[1]
        exStd2s[iImg] = exResids[2]
        exStd3s[iImg] = exResids[3]

        #have to fix this by hand for TSO FFIs
        hdulistCal[0].header['EXPTIME'] = 475.2
        #added by Scott Flemmings request for MAST archive
        hdulistCal[0].header['EQUINOX'] = 2000.0
        hdulistCal[0].header['INSTRUME'] = "TESS Photometer"
        hdulistCal[0].header['TELESCOP'] = "TESS"
        hdulistCal[0].header['FILTER'] = "TESS" 

        try:
            hdulistCal[0].header['MJD-BEG'] = hdulistCal[0].header['STARTTJD'] +\
                                              hdulistCal[0].header['TJD_ZERO'] - 2400000.5
            hdulistCal[0].header['MJD-END'] = hdulistCal[0].header['ENDTJD'] +\
                                              hdulistCal[0].header['TJD_ZERO'] - 2400000.5
        except:
            #for spoc origin???
            hdulistCal[0].header['MJD-BEG'] = hdulistCal[0].header['TSTART'] +\
                                              hdulistCal[0].header['BJDREFI'] - 2400000.5
            hdulistCal[0].header['MJD-END'] = hdulistCal[0].header['TSTOP'] +\
                                              hdulistCal[0].header['BJDREFI'] - 2400000.5
        
            #        print(  (hdulistCal[0].header['MJD-BEG'] - hdulistCal[0].header['MJD-END']) , 600.0/86400 )
        try:
            assert (hdulistCal[0].header['MJD-END'] - hdulistCal[0].header['MJD-BEG']  -  600.0/86400 ) < 1.e-10
        except KeyError:
            #no key in SPOC, just continue for now
            pass

        # Save header and control point data to fits file
        #  separate file for now.  Not touching input data
        # Add the wcsfit diagnostics to header

        try:
            hdulistCal[0].header.extend( add_btjd_info(hdulistCal[0].header['MIDTJD'], 
                                                       newhdr['CRVAL1'], 
                                                       newhdr['CRVAL2'],
                                                       ephemeris_data)  )
        except KeyError:
            #no key in SPOC, just continue for now
            pass

        #added by Scott Flemmings request for MAST archive
        newhdr['RA_TARG']   = newhdr['CRVAL1']
        newhdr['DEC_TARG'] = newhdr['CRVAL2']


        newhdr['RMSA'] = (allStd, 'WCS fit resid all targs [arcsec]')
        newhdr['RMSB'] = (brightStd, 'WCS fit resid bright (Tmag<10) targs [arcsec]')
        #newhdr['RMSF'] = (faintStd, 'WCS fit resid faint (Tmag>10) targs [arcsec]')
        newhdr['RMSAP'] = (allStdPix, 'WCS fit resid all targs [pixel]')
        newhdr['RMSBP'] = (brightStdPix, 'WCS fit resid bright (Tmag<10) targs [pixel]')
        #newhdr['RMSBF'] = (faintStdPix, 'WCS fit resid faint (Tmag>10) targs [pixel]')
        newhdr['RMSX0'] = (exResids[0], 'WCS fit resid extra 0 [arcsec]')
        newhdr['RMSX1'] = (exResids[1], 'WCS fit resid extra 1 [arcsec]')
        newhdr['RMSX2'] = (exResids[2], 'WCS fit resid extra 2 [arcsec]')
        newhdr['RMSX3'] = (exResids[3], 'WCS fit resid extra 3 [arcsec]')
        # Skip this header element in production 
        #newhdr['TIME'] = (ts[iImg], 'Time From Original Image Header')
        newhdr['WCSGDF' ] = (gdFracs[iImg], 'Fraction of control point targs valid')
        newhdr['CTRPCOL'] = (CTRL_PER_COL, 'Subregion analysis blocks over columns')
        newhdr['CTRPROW'] = (CTRL_PER_ROW, 'Subregion analysis blocks over rows')
        newhdr['FLXWIN'] = (blkHlf*2+1, 'Width in pixels of Flux-weight centroid region')
        # Make the fits table columns
        c1 = fits.Column(name='TIC', format='K', array=tics)
        c2 = fits.Column(name='FLXCOL', format='E', unit='pix', array=newCols)
        c3 = fits.Column(name='FLXROW', format='E', unit='pix', array=newRows)
        c4 = fits.Column(name='FLXVALID', format='I', array=gdPrfs)
        # Make the extension table 
        hduex = fits.BinTableHDU.from_columns([c1, c2, c3, c4])
        #append wcs parameter to primary header            
        hdulistCal[0].header.extend(newhdr, update=True)
        # Now merge the hdus
        all_hdus = fits.HDUList([hdulistCal[0], hduex])
        # Actually write fits file
        #fileBase = os.path.splitext(os.path.basename(curImg))[0]
        #outFits = os.path.join(outputDir, fileBase+'_wcs.fits')        
        all_hdus.writeto(curImg, checksum=True, overwrite=True)

        # Show where the prfs were deemed not good
        if DEBUG_LEVEL>2:
            plt.plot(newCols, newRows, '.')
            idx = np.where(gdPrfs == 0)[0]
            plt.plot(newCols[idx], newRows[idx], '.')
            idx = np.where(gdPrfs == -1)[0]
            plt.plot(newCols[idx], newRows[idx], '.k')
            idx = np.where(gdPrfs == -2)[0]
            plt.plot(newCols[idx], newRows[idx], '.r')
            plt.show()
        iImg = iImg +1

    # Finished all images make some diagnostic figures
    if DEBUG_LEVEL > 1 or saveDiag:
        if saveDiag:
            fileoutprefix = os.path.join(outputDir,'wcs_diags2','wcs_diag_S{0:d}_{1:d}{2:d}'.format( \
                                         SECTOR_WANT, CAMERA_WANT, CCD_WANT))
            # Save wcs fit diagnostics for images
            fout = h5py.File('{0}_data.h5'.format(fileoutprefix), 'w')
            asciiList = [n.encode("ascii", "ignore") for n in imgNames]
            tmp = fout.create_dataset('imgNames', (len(asciiList),),'S100', asciiList, compression='gzip')
            tmp = fout.create_dataset('allStds', data=allStds, compression='gzip')
            tmp = fout.create_dataset('brightStds', data=brightStds, compression='gzip')
            tmp = fout.create_dataset('faintStds', data=faintStds, compression='gzip')
            tmp = fout.create_dataset('allPixStds', data=allPixStds, compression='gzip')
            tmp = fout.create_dataset('brightPixStds', data=brightPixStds, compression='gzip')
            tmp = fout.create_dataset('faintPixStds', data=faintPixStds, compression='gzip')
            tmp = fout.create_dataset('ts', data=ts, compression='gzip')
            tmp = fout.create_dataset('exStd0s', data=exStd0s, compression='gzip')
            tmp = fout.create_dataset('exStd1s', data=exStd1s, compression='gzip')
            tmp = fout.create_dataset('exStd2s', data=exStd2s, compression='gzip')
            tmp = fout.create_dataset('exStd3s', data=exStd3s, compression='gzip')
            
            fout.close()

        # This is time series for position fit residuals in arcsec
        plt.plot(ts, allStds, '.', label='Fit Resid All Targets [arcsec]')
        plt.plot(ts, brightStds, '.', label='Fit Resid Bright Targets [arcsec]')
        #plt.plot(ts, faintStds, '.', label='Fit Resid Faint Targets [arcsec]')
        plt.plot(ts, gdFracs, '.', label='Fraction References Good')
        plt.ylabel('Various')
        plt.xlabel('Time [TJD]')
        plt.legend()
        plt.savefig('{0}_arcsec.png'.format(fileoutprefix), dpi=300)
        if DEBUG_LEVEL > 1:
            plt.show()
        plt.close()

        # This is time series for position fit residuals in arcsec
        # For the extra data series corner data
        plt.plot(ts, brightStds, '.', label='Fit Resid Bright Targets')
        plt.plot(ts, exStd0s, '.', label='Fit Resid Ex0')
        plt.plot(ts, exStd1s, '.', label='Fit Resid Ex1')
        plt.plot(ts, exStd2s, '.', label='Fit Resid Ex2')
        plt.plot(ts, exStd3s, '.', label='Fit Resid Ex3')
        plt.ylabel('WCS Fit residual [arcsec]')
        plt.xlabel('Time [TJD]')
        plt.legend()
        plt.savefig('{0}_extra_arcsec.png'.format(fileoutprefix), dpi=300)
        if DEBUG_LEVEL > 1:
            plt.show()
        plt.close()

    
        plt.plot(gdFracs/10.0, '.', label='Fraction References Good/10')
        plt.plot(allPixStds, '.', label='Fit Resid All Targets [pix]')
        plt.plot(brightPixStds, '.', label='Fit Resid Bright Targets [pix]')
        #plt.plot(faintPixStds, '.', label='Fit Resid Faint Targets [pix]')
        plt.ylabel('Various')
        plt.xlabel('Time [TJD]')
        plt.legend()
        plt.savefig('{0}_pix.png'.format(fileoutprefix), dpi=300)
        if DEBUG_LEVEL > 1:
            plt.show()
        plt.close()
 
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sector", type=int, \
                        help="TESS Observing Sector Number")
    parser.add_argument("-ca", "--camera", type=int, choices=range(1,5), \
                        help="Camera Number [1-4]")
    parser.add_argument("-cd", "--ccd", type=int, choices=range(1,5), \
                        help="CCD Number [1-4]")
    parser.add_argument("-rd", "--refdata", type=argparse.FileType('rb'), \
                        help="Reference image control point data created in step1")
    parser.add_argument("-if", "--imgfiles", nargs="*",
                        help="This requires including wildcard * character in filename like"\
                             " ls that would list the single camera/ccd image files."\
                            " You MUST PUT QUOTES AROUND THIS argrument with wildcard"\
                            " or else the shell expands the argument list to all images")
    parser.add_argument("-od", "--outputdir",\
                        help="Output directory for results.  NOTE: Currently"\
                             " writes separate fits file with wcs in header"\
                             " and target information in data.  It does not"\
                             " touch the data files yet.")
    parser.add_argument("-fd", "--fitdegree", type=int,\
                        help="Degree of wcs fit")
    parser.add_argument("--savediaginfo", action="store_true",\
                        help="Save diagnostic info and figures")
    parser.add_argument("-dbg", "--debug", type=int, \
                        help="Debug level; integer higher has more output")
    args = parser.parse_args()

# DEBUG BLOCK for hard coding input parameters and testing
#     class test_arg:
#         def __init__(self):
#             self.sector=29
#             self.camera=2
#             self.ccd=1
# #            self.refdata=open('/pdo/users/cjburke/tica/wcs_build/refout/ref_S29_Orbit66_21.h5','rb')
# #            self.imgfiles='/pdo/qlp-data/orbit-66/ffi/cam2/ccd1/FITS/*'
# #            self.outputdir='fitsout/Orbit66/Cam2/Ccd1'
#             self.refdata=open('/pdo/users/cjburke/tica/wcs_build/refout/refspocunidense_S29_21.h5','rb')
#             self.imgfiles='/pdo/spoc-data/sector-029/ffis/tess*-2-1-0193-s_ffic.fits.gz'
#             self.outputdir='fitsout/SpocS29/Cam2/Ccd1'
#             self.fitdegree=6
#             self.savediaginfo=True
#             self.debug=0
#     args = test_arg()

    SECTOR_WANT = args.sector
    CAMERA_WANT = args.camera
    CCD_WANT = args.ccd
    # argparse opens file handle.  hd5 needs file name not open file handle
    args.refdata.close()
    REF_DATA = args.refdata.name
    IMG_LIST_STR = args.imgfiles
#    print(IMG_LIST_STR)
    # We need to remove quotes if they exist in IMG_LIST_STR
    # Using GNUparallel to run this command Ihad to add the quotes
    #  to avoid the files to be expanded on command line
    #IMG_LIST_STR= IMG_LIST_STR.replace('"','')
    # output directory is same as input directory
    outputDir = os.path.dirname(IMG_LIST_STR[0])
    fitDegree = args.fitdegree
    saveDiag = False
    if args.savediaginfo:
        saveDiag = True
        saveDir = os.path.join(outputDir,'wcs_diags2')
        if not os.path.isdir(saveDir):
            os.mkdir(saveDir)
    DEBUG_LEVEL = args.debug
    
    fit_wcs_in_imgdir(SECTOR_WANT, CAMERA_WANT, CCD_WANT, REF_DATA, \
                    IMG_LIST_STR, outputDir, fitDegree, saveDiag, DEBUG_LEVEL)
        
    print('Done Sector {0:d} Camera {1:d} CCD: {2:d}'.format(SECTOR_WANT, \
                    CAMERA_WANT, CCD_WANT))
