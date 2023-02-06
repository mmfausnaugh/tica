#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 18:49:13 2020

@author: cjburke
"""

import numpy as np
from tess_stars2px import tess_stars2px_reverse_function_entry
from tess_stars2px import tess_stars2px_function_entry

try:
    import pyds9 as pd
except ImportError:
    print('Warning: No pyds9 installed.  No debugging with image display available')
from astropy.io import fits
#import photutils.centroids as cent
import matplotlib.pyplot as plt
import h5py

import os
import sys

import tica
from time import time
import logging

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))

#from wcs_build.mast_filter_conesearch import mast_filter_conesearch
from wcs_build.tic_local_conesearch import tic_local_conesearch


from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gwcs.wcstools import wcs_from_points
import argparse

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list


def flux_weighted_centroid(flux_stamp):
    r,c = np.shape(flux_stamp)
    #C,R = np.meshgrid( np.r_[0:r],
    #                   np.c_[0:c])
    R,C = np.ogrid[0:r,0:c]
    norm = np.sum(flux_stamp)
    centroid_row = np.sum(flux_stamp*R)/norm
    centroid_col = np.sum(flux_stamp*C)/norm

    #centroid_row = centroid_row - np.shape(flux_stamp)[0]/2.0 + 0.5
    #centroid_col = centroid_col - np.shape(flux_stamp)[1]/2.0 + 0.5
    return np.array([centroid_col, centroid_row, norm])

def gdPRF_calc(img, blkHlf, wingFAC=0.9, contrastFAC=3.5):
    rowSum = np.sum(img, axis=1)
    colSum = np.sum(img, axis=0)
    # check to make sure middle is higher than both wings 
    gdPRF = True
    midPk = np.max(rowSum[blkHlf-2:blkHlf+3])
    lowPk = np.max(rowSum[0:3])
    hghPk = np.max(rowSum[-3:])
    delLowPk = midPk - lowPk
    delHghPk = midPk - hghPk
    maxDelPk = np.max([delLowPk, delHghPk])
    #this bounds the symmetry of the PSF:
    #the maxDelPk corresponds to the minimum of the wings, so
    #this always fails if wingFAC = 1.0
    #if wingFAC is smaller, there is space between 
    #the minimum of the wings and some fraction of the peak
    #in which the other wing can fit
    trigLev = midPk - maxDelPk*wingFAC
    
    if lowPk > trigLev or hghPk > trigLev or maxDelPk < 10.0:
        gdPRF = False
    # Check to make sure contrast of middle is higher than wings
    contrastRow = midPk/(np.max([lowPk,hghPk]))
    if contrastRow < contrastFAC:
        gdPRF = False
    midPk = np.max(colSum[blkHlf-2:blkHlf+3])
    lowPk = np.max(colSum[0:3])
    hghPk = np.max(colSum[-3:])
    delLowPk = midPk - lowPk
    delHghPk = midPk - hghPk
    maxDelPk = np.max([delLowPk, delHghPk])
    trigLev = midPk - maxDelPk*wingFAC
    if lowPk > trigLev or hghPk > trigLev or maxDelPk < 10.0:
        gdPRF = False
    # Check to make sure contrast of middle is higher than wings
    contrastCol = midPk/(np.max([lowPk,hghPk]))
    if contrastCol < contrastFAC:
        gdPRF = False
    return gdPRF, contrastCol, contrastRow

def ring_background(img):   
    tmp = np.array([], dtype=np.double)
    tmp = np.append(tmp, img[0:2,:].flatten())
    tmp = np.append(tmp, img[-2:,:].flatten())
    tmp = np.append(tmp, img[2:-2,0:2].flatten())
    tmp = np.append(tmp, img[-2:,-2:].flatten())
    return np.median(tmp)

def getRoughTranslation(col, row, cam, ccd, img, region, limits, cam_want, ccd_want, blkHlf):
    midCols = np.round(col).astype(int)
    midRows = np.round(row).astype(int)
    
    # Make sure col, row, cam ,ccd are within region of interest
    idx = np.where( (midCols > region[0]) &\
                    (midCols <= region[1]) &\
                    (midRows > region[2]) &\
                    (midRows <= region[3]) &\
                    (midCols-blkHlf >= limits[0]) &\
                    (midCols+blkHlf <= limits[1]) &\
                    (midRows-blkHlf >= limits[2]) &\
                    (midRows+blkHlf <= limits[3]) &\
                    (cam == cam_want) & (ccd == ccd_want))
    colUse = midCols[idx]
    rowUse = midRows[idx]
    track_cols = np.array([], dtype=np.int32)
    track_rows = np.array([], dtype=np.int32)
    for i in range(len(colUse)):
        curCol = colUse[i]
        curRow = rowUse[i]
        subImg = img[curRow-blkHlf-1: curRow+blkHlf, curCol-blkHlf-1: curCol+blkHlf]
        rowSum = np.sum(subImg, axis=1)
        colSum = np.sum(subImg, axis=0)
        # check to make sure middle is higher than both wings 
        midPk = np.max(rowSum[1:-2])
        bkg = np.max([rowSum[0], rowSum[-1]])
        if midPk > 2.0*bkg:
            track_rows = np.append(track_rows, np.argmax(rowSum[1:-2])+1)
        midPk = np.max(colSum[1:-2])
        bkg = np.max([colSum[0], colSum[-1]])
        if midPk > 2.0*bkg:
            track_cols = np.append(track_cols, np.argmax(colSum[1:-2])+1)
    offcol = np.median(track_cols) - blkHlf
    offrow = np.median(track_rows) - blkHlf
    #print('off column: {0:f} row: {1:f}'.format(offcol, offrow))
    #plt.hist(track_cols)
    #plt.show()
    #plt.hist(track_rows)
    #plt.show()
    
    return offcol, offrow

def apply_proper_motion(ras, decs, pmras, pmdecs, delyr):
    
    ccat = SkyCoord(ra=ras * u.deg,
                    dec=decs * u.deg,
                    pm_ra_cosdec=pmras * u.mas/u.yr,
                    pm_dec=pmdecs * u.mas/u.yr,
                    obstime=Time('J2000.0'),
                    distance=1000.0*u.pc,
                    radial_velocity = 0.0 * u.km/u.s)
    # Convert coordinates to epoch 2021.0
    curEpc = Time('J2021.0')
    ccat_cur = ccat.apply_space_motion(curEpc)
    
    return ccat_cur.ra.degree, ccat_cur.dec.degree    

def get_refimg_ctrlpts(SECTOR_WANT, CAMERA_WANT, CCD_WANT, REF_IMAGE, outputFile, 
                       DEBUG_LEVEL=0,
                       wingFAC = 0.9, contrastFAC=3.5):
        
    # Use DS9 if it exists and debug level is high
    DEBUG_DS9 = True
    if DEBUG_LEVEL>1:
        try:
            dispTargs = pd.ds9_targets()
            dispDS9 = pd.DS9(dispTargs[0].split()[1])
        except:
            DEBUG_DS9 = False
    # Open reference image
    hdulistCal = fits.open(REF_IMAGE)
    
    # Here are some hardcoded parametrers
    pixScl = 21.0 # arcsec 1 pixel pixScl arcsec
    #  The detector has 44 leading bias columns
    # The follow col and row min max parameters set the first
    #  valid science pixel.  Note this uses the ds9 wcs convention
    #  that middle of first pixel is 1 tess_stars2px.py also uses this convention
    #  for predicted locations
    colMin = 45 
    colMax = 2092
    rowMin = 1
    rowMax = 2048
    #  These control how many regions in each dimension to find stars
    #  and perform a MASt TICsearch over each one.  This is probably
    #   not needed, but it is how the code was structured originally
    #   and may save on overhead in make SQLqueries large?
    CTRL_PER_COL = 5
    CTRL_PER_ROW = 5
    # number of pixels above and below midpoint to determine centroid analysis region
    #  I wouldn't make it larger than 5 unless TESS really had a pointing excursion
    blkHlf = 5
    blkHlfCent = 2 # The analysis region for background estimate, gdPRF check, initial
                    # position offets use blkHlf. However, final centroid is calcualted
                    # with a small blkHlfCent region
    # Tmag min and max for controll point targets
    TMAG_MIN = 7.5
    TMAG_MAX = 12.0
    # Trim out control points if they are offset from 
    #  wcs fit by OUTLIER_TRIM [arcsec]
    OUTLIER_TRIM = 6.0
    # Max number of stars to keep
    #  These will get broken down by sub region
    MAX_NUSE = 1000
    MAX_SUBREGION_USE = MAX_NUSE // (CTRL_PER_COL*CTRL_PER_ROW)
    #print('MAX PER SubRegion {0:d}'.format(MAX_SUBREGION_USE))
    logging.info('MAX PER SubRegion {0:d}'.format(MAX_SUBREGION_USE))
    
    # The following sets up the analysis subregions there are CTRL_PER_COL x CTRL_PER_Row 
    #  subregions
    spanCol = colMax - colMin
    spanRow = rowMax - rowMin
    delCol = float(spanCol) / float(CTRL_PER_COL)
    delRow = float(spanRow)/ float(CTRL_PER_ROW)
    #  This sets the pixel coordinates for the center each subregion
    colCtrl1D = np.linspace(delCol/2.0 + colMin, colMax - delCol/2.0, CTRL_PER_COL)
    rowCtrl1D = np.linspace(delRow/2.0 + rowMin, rowMax - delRow/2.0, CTRL_PER_ROW)
    colCtrl2D, rowCtrl2D = np.meshgrid(colCtrl1D, rowCtrl1D)    
    colCtrl2D_flat = colCtrl2D.flatten()
    rowCtrl2D_flat = rowCtrl2D.flatten()
    nCtrl = len(colCtrl2D_flat)
    
    # this Loop will then use tess_stars2px to get a rough
    #  idea for the ra and dec for each sub region center.
    raCtrl2D_flat = np.zeros_like(colCtrl2D_flat)
    decCtrl2D_flat = np.zeros_like(rowCtrl2D_flat)
    for i,curCol in enumerate(colCtrl2D_flat):
        curRow = rowCtrl2D_flat[i]
        # on first time scinfo does not exist; it is faster to pass it in on subsequent calls
        #  so it doesn't need to be calculated each time.
        if i == 0:
            raCtrl2D_flat[i], decCtrl2D_flat[i], scinfo = tess_stars2px_reverse_function_entry(\
                         SECTOR_WANT, CAMERA_WANT, CCD_WANT, curCol, curRow)
        else:
            raCtrl2D_flat[i], decCtrl2D_flat[i], scinfo = tess_stars2px_reverse_function_entry(\
                         SECTOR_WANT, CAMERA_WANT, CCD_WANT, curCol, curRow, scinfo)
        if np.mod(i,CTRL_PER_COL) == 0 and DEBUG_LEVEL>0: # print progress with sufficient debug_level
            logging.debug('{0:d} of {1:d}'.format(i, nCtrl))
            logging.debug('{0:f} {1:f} {2:f} {3:f}'.format(curCol, curRow, raCtrl2D_flat[i], decCtrl2D_flat[i]))
            #print('{0:d} of {1:d}'.format(i, nCtrl))
            #print('{0:f} {1:f} {2:f} {3:f}'.format(curCol, curRow, raCtrl2D_flat[i], decCtrl2D_flat[i]))

    # In each control point grid keep track of how many are valid 
    numKept = np.zeros((nCtrl,), dtype=np.int32)
    # Keep track of how many were tried
    numBlk = np.zeros_like(numKept)
    
    # Now lets do some mast cone searchs
    # and find some nice clean stars to centroid
    kpTics = np.array([], dtype=np.int64)
    kpRas = np.array([], dtype=np.double)
    kpDecs = np.array([], dtype=np.double)
    kpTmags = np.array([], dtype=np.double)
    kpCtrlIdxs = np.array([], dtype=np.int32)
    kpPredCols = np.array([], dtype=np.double)
    kpPredRows = np.array([], dtype=np.double)
    kpObsCols = np.array([], dtype=np.double)
    kpObsRows = np.array([], dtype=np.double)
    kpContrastCols = np.array([], dtype=np.double)
    kpContrastRows = np.array([], dtype=np.double)
    kpApCols = np.array([], dtype=np.int32)
    kpApRows = np.array([], dtype=np.int32)

    # Specify cone search radius
    radSearch = np.max([delCol, delRow])*np.sqrt(2.0)/2.0 * pixScl
    # Here is the main work loop where over each sub image region
    #  query TIC to get stars
    #  go over stars and make sure they have well behaved PRFs for isloated stars
    for i, curRa in enumerate(raCtrl2D_flat):
        curDec = decCtrl2D_flat[i]
        # Do mast cone Search on this subimage region
        # tics, ticRas, ticDecs, ticTmags, ticKmags, ticGmags, ticpmRAs, ticpmDecs = mast_filter_conesearch(curRa, curDec, radSearch, TMAG_MIN, TMAG_MAX)
        tics, ticRas, ticDecs, ticTmags, ticKmags, ticGmags, ticpmRAs, ticpmDecs = tic_local_conesearch(curRa, curDec, radSearch, TMAG_MIN, TMAG_MAX)
        nSrch = len(tics)
        # Sort the targets with brightest first
        idx = np.argsort(ticTmags)
        tics, ticRas, ticDecs, ticTmags, ticpmRAs, ticpmDecs = idx_filter(idx, \
                        tics, ticRas, ticDecs, ticTmags, ticpmRAs, ticpmDecs)

        #print(len(tics), ' ', curRa,' ', curDec, ' ', radSearch)
        ticRas, ticDecs = apply_proper_motion(ticRas, ticDecs, ticpmRAs, ticpmDecs, 21.0)

        # Get the predicted positions for these targets found in cone search using tess_stasrs2px reverese
        outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, \
            outColPix, outRowPix, scinfo = tess_stars2px_function_entry(tics, ticRas, ticDecs, SECTOR_WANT, scinfo, aberrate=True)
        # Make sure coordinates are in the current block of interest
        colLow = int(round(colCtrl2D_flat[i]-delCol/2.0))
        colHgh = int(round(colLow + delCol))
        rowLow = int(round(rowCtrl2D_flat[i]-delRow/2.0))
        rowHgh = int(round(rowLow + delRow))
        
        # We do an initial round of peak finding for nTry targets
        #  to get a first order translation between prediction and measured
        #  This will give us better centered flux weighted centroid determinations
        #   to hopefully prevent a Tmag dependence
        cold, rowd = getRoughTranslation(outColPix, outRowPix, outCam, outCcd, \
                        hdulistCal[0].data, \
                        [colLow, colHgh, rowLow, rowHgh], [colMin, colMax, rowMin, rowMax], \
                        CAMERA_WANT, CCD_WANT, blkHlf)
        outColPix = outColPix + cold
        outRowPix = outRowPix + rowd
        if DEBUG_LEVEL > 0:
            #print('Initial Position Tweaks Col: {0:f} Row: {1:f}'.format(cold, rowd))
            logging.debug('Initial Position Tweaks Col: {0:f} Row: {1:f}'.format(cold, rowd))

        nTry = 5
        gotN = 0
        j = 0
        delCols = np.array([], dtype=np.double)
        delRows = np.array([], dtype=np.double)
        while (gotN <= nTry and j<len(outID)):
            # this is predicted coordinates
            midCol = int(round(outColPix[j]))
            midRow = int(round(outRowPix[j]))
            curCam = outCam[j]
            curCcd = outCcd[j]
            curTic = outID[j]
            # Make sure coordinates are in the current block of interest
            if midCol > colLow and midCol <= colHgh and midRow > rowLow and midRow <= rowHgh and curCam == CAMERA_WANT and curCcd == CCD_WANT:
                # Also make sure we aren't too close to edge
                if midCol-blkHlf >= colMin and midCol+blkHlf <= colMax and midRow-blkHlf >= rowMin and midRow+blkHlf <= rowMax:
                    #  pixel coordinates in x and y directions
                    colX = np.arange(midCol-blkHlf, midCol+blkHlf+1)
                    rowY = np.arange(midRow-blkHlf, midRow+blkHlf+1)
                    #  centroid analysis region
                    sciImgCal = hdulistCal[0].data[midRow-blkHlf-1: midRow+blkHlf, midCol-blkHlf-1: midCol+blkHlf]
                    
                    # Check if target is isolated
                    gdPRF, contrastCol, contrastRow = gdPRF_calc(sciImgCal, blkHlf,
                                                                 wingFAC = wingFAC, 
                                                                 contrastFAC=contrastFAC)
                    if gdPRF:                    # Determine background from 2 pixel ring around box
                        bkgLev = ring_background(sciImgCal)
                        #centOut = cent.centroid_com(sciImgCal-bkgLev)
                        centOut = flux_weighted_centroid( sciImgCal - bkgLev)
                        #print(centOut, centOut1)
                        # centroid_com produces 0-based coordinates
                        newCol = centOut[0]+colX[0]
                        newRow = centOut[1]+rowY[0]
                        delCols = np.append(delCols, newCol - outColPix[j])
                        delRows = np.append(delRows, newRow - outRowPix[j])
                        gotN = gotN + 1
            j = j+1
        if (gotN < nTry):
            print('Something went wrong determining initial pixel position correction on reference image!')
            print('Could not find nTry: {0:d} targets passing gdPRF'.format(nTry))
            logging.error('Something went wrong determining initial pixel position correction on reference image!')
            logging.error('Could not find nTry: {0:d} targets passing gdPRF'.format(nTry))
            exit()
        # From these nTry targets calculate the col and row translation from prediction
        corCol = np.median(delCols)
        corRow = np.median(delRows)
        if DEBUG_LEVEL > 0:
            #print('Predicted Position Tweaks Col: {0:f} Row: {1:f}'.format(corCol, corRow))
            logging.debug('Predicted Position Tweaks Col: {0:f} Row: {1:f}'.format(corCol, corRow))
        outColPix = outColPix + corCol
        outRowPix = outRowPix + corRow
        
        # Armed with the translation offset from prediction estimate
        #  Now loop over all potential targets to find the control points that 
        #  pass gdPRF
        for j, curCol in enumerate(outColPix):
            # this is predicted coordinates
            midCol = int(round(curCol))
            midRow = int(round(outRowPix[j]))
            curCam = outCam[j]
            curCcd = outCcd[j]
            curTic = outID[j]
            # Make sure coordinates are in the current block of interest
            if midCol > colLow and midCol <= colHgh and midRow > rowLow and midRow <= rowHgh and curCam == CAMERA_WANT and curCcd == CCD_WANT:
                # Also make sure we aren't too close to edge
                if midCol-blkHlf >= colMin and midCol+blkHlf <= colMax and midRow-blkHlf >= rowMin and midRow+blkHlf <= rowMax:
                    # Attempt to find a good star
                    numBlk[i] = numBlk[i] + 1
                    #  pixel coordinates in x and y directions
                    colX = np.arange(midCol-blkHlf, midCol+blkHlf+1)
                    rowY = np.arange(midRow-blkHlf, midRow+blkHlf+1)
                    # analysis region around target
                    sciImgCal = hdulistCal[0].data[midRow-blkHlf-1: midRow+blkHlf, midCol-blkHlf-1: midCol+blkHlf]
                    ## TEST TEST 
                    ## do an image with single pixel with counts others are zero to
                    # test zero vs one-based coordinates inject peak at col 6: row: 4 in ds9 land
                    #sciImgCal = 0*sciImgCal
                    #sciImgCal[4, 5] = 1000.0
                    
                    # Check if target is isolated
                    gdPRF, contrastCol, contrastRow = gdPRF_calc(sciImgCal, blkHlf,
                                                                 wingFAC = wingFAC, 
                                                                 contrastFAC=contrastFAC)

                    # Determine background from 2 pixel ring around box
                    bkgLev = ring_background(sciImgCal)

                    # Testing block to really see how centroiding is working            
#                    print('Col: {0:d} {1:5.3f} Row: {2:d} {3:5.3f} PixValue: {4} Tic: {5:d} Good? {6:b}'.format(\
#                          midCol, curCol, midRow, outRowPix[j], hdulistCal[0].data[midRow-1,midCol-1], curTic, gdPRF))
#                    print('Background: {0:f}'.format(bkgLev))
#                    dispDS9.set('frame 1')
#                    dispDS9.set_np2arr(sciImgCal)
#                    rowSum = np.sum(sciImgCal, axis=1)
#                    colSum = np.sum(sciImgCal, axis=0)
#
#                    centOut = cent.centroid_com(sciImgCal-bkgLev)
#                    # centroid_com produces 0-based coordinates
#                    newCol = centOut[0]+colX[0]
#                    newRow = centOut[1]+rowY[0]
#                    print('ObsCol: {0:f} {1:f} ObsRow: {2:f} {3:f}'.format(centOut[0], newCol, centOut[1], newRow))
#                    # This first display is the test case 
#                    #dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 4'.format(6, 5))
#                    dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 4 color=red'.format(outColPix[j]-colX[0]+1.0, outRowPix[j]-rowY[0]+1.0))
#                    dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 4 '.format(centOut[0]+1.0, centOut[1]+1.0))
#
#                    plt.subplot(121)
#                    plt.plot(colX, colSum, '-')
#                    plt.axvline(outColPix[j], color='r')
#                    plt.axvline(newCol, color='g')
#                    plt.xlabel('CCD Column Position [px]')
#                    plt.subplot(122)
#                    plt.plot(rowY, rowSum, '-')
#                    plt.axvline(outRowPix[j], color='r')
#                    plt.axvline(newRow, color='g')
#                    plt.xlabel('CCD Row Position [px]')
#                    plt.show()

                    if gdPRF:
                        #  pixel coordinates in x and y directions
                        # Use small blkHlfCent region for final centroid calc
                        colX = np.arange(midCol-blkHlfCent, midCol+blkHlfCent+1)
                        rowY = np.arange(midRow-blkHlfCent, midRow+blkHlfCent+1)
                        # analysis region around target
                        sciImgCal = hdulistCal[0].data[midRow-blkHlfCent-1: midRow+blkHlfCent, midCol-blkHlfCent-1: midCol+blkHlfCent]
                        #centOut = cent.centroid_com(sciImgCal-bkgLev)
                        centOut = flux_weighted_centroid( sciImgCal - bkgLev)
                        # centroid_com produces 0-based coordinates
                        newCol = centOut[0]+colX[0]
                        newRow = centOut[1]+rowY[0]
                        # Record a valid centroid
                        numKept[i] = numKept[i] + 1
                        kpTics = np.append(kpTics, outID[j])
                        jj = np.where(tics == outID[j])[0]
                        kpRas = np.append(kpRas, ticRas[jj])
                        kpDecs = np.append(kpDecs, ticDecs[jj])
                        kpTmags = np.append(kpTmags, ticTmags[jj])
                        kpCtrlIdxs = np.append(kpCtrlIdxs, i)
                        kpPredCols = np.append(kpPredCols, outColPix[j])
                        kpPredRows = np.append(kpPredRows, outRowPix[j])
                        kpObsCols = np.append(kpObsCols, newCol)
                        kpObsRows = np.append(kpObsRows, newRow)
                        kpContrastCols = np.append(kpContrastCols, contrastCol)
                        kpContrastRows = np.append(kpContrastRows, contrastRow)
                        kpApCols = np.append(kpApCols, int(midCol))
                        kpApRows = np.append(kpApRows, int(midRow))

        if DEBUG_LEVEL>0:
            cur_idx = np.where(kpCtrlIdxs == i)[0]
            #print('Image Coords Col: {0:d} {1:d} '
            #      'Row: {2:d} {3:d} NGd: {4:d} Ntry: {5:d}'.format(
            #          colLow, colHgh, rowLow, rowHgh, len(cur_idx), nSrch))
            logging.debug('Image Coords Col: {0:d} {1:d} '
                         'Row: {2:d} {3:d} NGd: {4:d} Ntry: {5:d}'.format(
                             colLow, colHgh, rowLow, rowHgh, len(cur_idx), nSrch))
                
    
                
                
        # If debugging ds9  available display image sub region with stars accepted
        #  as reference control points marked on image
        if DEBUG_LEVEL>1 and DEBUG_DS9:
            sciImgCal = hdulistCal[0].data[rowLow-1:rowHgh,colLow-1:colHgh]
            dispDS9.set('frame 1')
            dispDS9.set_np2arr(sciImgCal)
            if len(cur_idx)>0:
                corColPix = kpPredCols[cur_idx] - colLow + 1.0
                corRowPix = kpPredRows[cur_idx] - rowLow + 1.0
                corNewColPix = kpObsCols[cur_idx] - colLow + 1.0
                corNewRowPix = kpObsRows[cur_idx] - rowLow + 1.0
                for j, curCol in enumerate(corColPix):
                    curRow = corRowPix[j]
                    dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 4'.format(curCol, curRow))
                    dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 4 color=red'.format(corNewColPix[j], corNewRowPix[j]))

    # There is a trend for fields with more targets and crowding have higher
    # wcs fit residual for a fixed gdPRF contrast parameter.
    # moving towards accepting a fixed number of the highest contrast
    # Do this here
    useContrast = np.zeros_like(kpContrastCols)
    keepIt = np.zeros_like(kpCtrlIdxs)
    for i in range(len(useContrast)):
        useContrast[i] = np.min([kpContrastCols[i], kpContrastRows[i]])
    # For each subregion only use the highest contrast targets
    for i in range(len(raCtrl2D_flat)):
        cur_idx = np.where(kpCtrlIdxs == i)[0]
        curContrast = useContrast[cur_idx]
        nCur = len(cur_idx)
        if nCur > MAX_SUBREGION_USE:
            idxsort = np.argsort(curContrast)[::-1]
            idxsort = idxsort[0:MAX_SUBREGION_USE]
            keepIt[cur_idx[idxsort]] = 1
            #print('Region:{0:d} clipout {1:d}'.format(i, nCur-MAX_SUBREGION_USE))
            logging.info('Region:{0:d} clipout {1:d}'.format(i, nCur-MAX_SUBREGION_USE))
        else:
            keepIt[cur_idx] = 1
    idx = np.where(keepIt == 1)[0]
    kpTics, kpRas, kpDecs, kpTmags, kpCtrlIdxs,\
        kpPredCols, kpPredRows, kpObsCols, kpObsRows, \
        kpContrastCols, kpContrastRows, \
        kpApCols, kpApRows = idx_filter(idx, 
                                        kpTics, kpRas, kpDecs, kpTmags, kpCtrlIdxs, 
                                        kpPredCols, kpPredRows, kpObsCols, kpObsRows,
                                        kpContrastCols, kpContrastRows,
                                        kpApCols, kpApRows)


    # Fit a wcs to see the residuals and trim outliers
    fitDegree = 6
    REFPIXCOL = 1024.0+45.0
    REFPIXROW = 1024.0
    PIX2DEG = 21.0/3600.0 # Turn pixels to degrees roughly
        # The WCS fitter uses WCS intermediate coordinates which are in degrees
        # Relative to reference ra and dec as th eprojection zero point.
    # From reference pixel coordinates get the estimated ra and dec of this point
    raproj, decproj, scinfo = tess_stars2px_reverse_function_entry(\
                         SECTOR_WANT, CAMERA_WANT, CCD_WANT, REFPIXCOL, REFPIXROW)
    proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
    #  Reference subtracted pixel coordinates
    sclObsCols = (kpObsCols - REFPIXCOL) * PIX2DEG
    sclObsRows = (kpObsRows - REFPIXROW) * PIX2DEG

    xy = (sclObsCols, sclObsRows)
    radec = (kpRas, kpDecs)
    gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
    gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
    deg2Rad = np.pi/180.0
    deltaRas = (gwcs_pred_ras - kpRas) *3600.0 * np.cos(kpDecs*deg2Rad)
    deltaDecs = (gwcs_pred_decs - kpDecs) *3600.0
    deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
    idxbd = np.where(deltaSeps >= OUTLIER_TRIM)[0]
    idxgd = np.where(deltaSeps < OUTLIER_TRIM)[0]
    # Filter out the outliers 
    kpTics, kpRas, kpDecs, kpTmags, kpCtrlIdxs,\
        kpPredCols, kpPredRows, kpObsCols, kpObsRows, \
        kpContrastCols, kpContrastRows, \
        kpApCols, kpApRows  = idx_filter(idxgd, 
                                         kpTics, kpRas, kpDecs, kpTmags, kpCtrlIdxs, 
                                         kpPredCols, kpPredRows, kpObsCols, kpObsRows,
                                         kpContrastCols, kpContrastRows,
                                         kpApCols, kpApRows)
                    

    #lastly, sort by tmag
    idx_out = np.argsort(kpTmags)
    # Filter out the outliers 
    kpTics, kpRas, kpDecs, kpTmags, kpCtrlIdxs,\
        kpPredCols, kpPredRows, kpObsCols, kpObsRows, \
        kpContrastCols, kpContrastRows, \
        kpApCols, kpApRows  = idx_filter(idx_out, 
                                         kpTics, kpRas, kpDecs, kpTmags, kpCtrlIdxs, 
                                         kpPredCols, kpPredRows, kpObsCols, kpObsRows,
                                         kpContrastCols, kpContrastRows,
                                         kpApCols, kpApRows)



    # Save data reference image results
    fout = h5py.File(outputFile, 'w')


    tmp = fout.create_dataset('tics', data=kpTics, compression='gzip')
    tmp = fout.create_dataset('ras', data=kpRas, compression='gzip')
    tmp = fout.create_dataset('decs', data=kpDecs, compression='gzip')
    tmp = fout.create_dataset('tmags', data=kpTmags, compression='gzip')
    tmp = fout.create_dataset('blkidxs', data=kpCtrlIdxs, compression='gzip')
    tmp = fout.create_dataset('predcols', data=kpPredCols, compression='gzip')
    tmp = fout.create_dataset('predrows', data=kpPredRows, compression='gzip')
    tmp = fout.create_dataset('obscols', data=kpObsCols, compression='gzip')
    tmp = fout.create_dataset('obsrows', data=kpObsRows, compression='gzip')
    tmp = fout.create_dataset('contrastcols', data=kpContrastCols, compression='gzip')
    tmp = fout.create_dataset('contrastrows', data=kpContrastRows, compression='gzip')
    tmp = fout.create_dataset('aperturecols', data=kpApCols, compression='gzip')
    tmp = fout.create_dataset('aperturerows', data=kpApRows, compression='gzip')

    #example REF_IMAGE: tess2022072133152-00204040-1-crm-ffi_ccd1.cal.fits
    fout.attrs['ref_FIN'] = REF_IMAGE.split('-')[1]
    fout.attrs['wingFAC']     = wingFAC
    fout.attrs['contrastFAC'] = contrastFAC
#save the 
    fout.close()

    # Fit a wcs after trimming to report the residuals
    #  Reference subtracted pixel coordinates
    sclObsCols = (kpObsCols - REFPIXCOL) * PIX2DEG
    sclObsRows = (kpObsRows - REFPIXROW) * PIX2DEG

    xy = (sclObsCols, sclObsRows)
    radec = (kpRas, kpDecs)
    gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
    gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
    deg2Rad = np.pi/180.0
    deltaRas = (gwcs_pred_ras - kpRas) *3600.0 * np.cos(kpDecs*deg2Rad)
    deltaDecs = (gwcs_pred_decs - kpDecs) *3600.0
    deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
    c1 = 1.0/np.sqrt(2.0)
    idxb = np.where(kpTmags<10.0)[0]
    std1 = np.std(deltaRas[idxb])
    std2 = np.std(deltaDecs[idxb])
    brightstd = np.sqrt(std1*std1+std2*std2)*c1
    idx = np.where(kpTmags>10.0)[0]
    std1 = np.std(deltaRas[idx])
    std2 = np.std(deltaDecs[idx])
    faintstd = np.sqrt(std1*std1+std2*std2)*c1
    std1 = np.std(deltaRas)
    std2 = np.std(deltaDecs)
    allstd = np.sqrt(std1*std1+std2*std2)*c1
    #print('Found: {0:d} ref targets, {1:d} bright, fit resid [arcsec] bright {2:6.3f} faint {3:6.3f} nTrim: {4:d}'.format(len(kpTics), len(idxb), brightstd, faintstd, len(idxbd)))
    logging.info('Found: {0:d} ref targets, {1:d} bright, fit resid [arcsec] bright {2:6.3f} faint {3:6.3f} nTrim: {4:d}'.format(len(kpTics), len(idxb), brightstd, faintstd, len(idxbd)))
    
    if DEBUG_LEVEL>1:
        #save a bunch of plots
        outputFileRoot = os.path.splitext(outputFile)[0]


        plt.plot(kpObsCols, kpObsRows, '.')
        plt.axhline(rowMin, ls='-', color='k')
        plt.axhline(rowMax, ls='-', color='k')
        plt.axvline(colMin, ls='-', color='k')
        plt.axvline(colMax, ls='-', color='k')
        plt.xlabel('Column [pix]')
        plt.ylabel('Row [pix]')
        plt.title('Reference Target Positions on Detector')
        plt.savefig(outputFileRoot+'_ref_target_positions.png', dpi=300)
        plt.clf()

        #plt.show()
        plt.plot(kpTmags, deltaRas, '.')
        plt.xlabel('Tmag')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.title('WCS Fit RA Res Vs. Tmag')
        plt.savefig(outputFileRoot+'_tmag_resid_ra.png', dpi=300)
        plt.clf()


        #plt.show()
        plt.plot(kpTmags, deltaDecs, '.')
        plt.xlabel('Tmag')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.title('WCS FIT Dec Res Vs Tmag')
        plt.savefig(outputFileRoot+'_tmag_resid_dec.png', dpi=300)
        plt.clf()

        #plt.show()
        plt.plot(kpObsCols, deltaRas, '.')
        plt.xlabel('Column [px]')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.title('WCS Fit RA Res Vs. Col')
        plt.savefig(outputFileRoot+'_col_resid_ra.png', dpi=300)
        plt.clf()

        #plt.show()
        plt.plot(kpObsRows, deltaRas, '.')
        plt.xlabel('Row [px]')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.title('WCS Fit RA Res Vs. Row')
        plt.savefig(outputFileRoot+'_row_resid_ra.png', dpi=300)
        plt.clf()

        #plt.show()
        plt.plot(kpObsCols, deltaDecs, '.')
        plt.xlabel('Column [px]')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.title('WCS Fit Dec res Vs. Col')
        plt.savefig(outputFileRoot+'_col_resid_dec.png', dpi=300)
        plt.clf()
        #plt.show()
        plt.plot(kpObsRows, deltaDecs, '.')
        plt.xlabel('Row [px]')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.title('WCS Fit Dec Res Vs. Row')
        plt.savefig(outputFileRoot+'_row_resid_dec.png', dpi=300)
        plt.clf()
        #plt.show()

        x = np.sort(deltaSeps[idxb])
        ndata = len(x)
        y = np.arange(ndata)/ float(ndata)
        plt.plot(x, y, '-')
        plt.xlabel('RA & Dec WCS residuals quadrature add [arcsec]')
        plt.ylabel('CDF')
        plt.title('CDF of Residuals (Tmag<10)')
        plt.savefig(outputFileRoot+'_cdf_resid_bright.png', dpi=300)
        plt.clf()

    # Lets look for wcs outliers
    #  This is at very high debugging level
    if DEBUG_LEVEL>2 and DEBUG_DS9:
        ia = np.where(kpTmags < 16.0)[0]
        
        # From reference pixel coordinates get the estimated ra and dec of this point
        raproj, decproj, scinfo = tess_stars2px_reverse_function_entry(\
                             SECTOR_WANT, CAMERA_WANT, CCD_WANT, REFPIXCOL, REFPIXROW)
        proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
        #  Reference subtracted pixel coordinates
        sclObsCols = (kpObsCols - REFPIXCOL) * PIX2DEG
        sclObsRows = (kpObsRows - REFPIXROW) * PIX2DEG
    
        xy = (sclObsCols[ia], sclObsRows[ia])
        radec = (kpRas[ia], kpDecs[ia])
        print('Start computing wcs')
        gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
        # Iterate to find a better reference ra and dec
        #  such that there is no constant term in the polynomial fit
        newrefra, newrefdec = gwcs_obj(0.0, 0.0)
        diffra = (raproj-newrefra)*3600.0
        diffdec = (decproj-newrefdec)*3600.0
        print('Old RA: {0:9.5f} NewRA: {1:9.5f} DelRa: {2:9.4f} OldDec: {3:9.5f} NewDec: {4:9.5f} DelDec: {5:9.4f}'.format(\
                  raproj, newrefra, np.sqrt(diffra*diffra), decproj, newrefdec, np.sqrt(diffdec*diffdec)))
        raproj = newrefra
        decproj = newrefdec
        proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
        print('Start computing wcs')
        gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
        gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
        deg2Rad = np.pi/180.0
        deltaRas = (gwcs_pred_ras - kpRas) *3600.0 * np.cos(kpDecs*deg2Rad)
        deltaDecs = (gwcs_pred_decs - kpDecs) *3600.0
        deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
        idx = np.where(kpTmags < 9.0)[0]
        tmpCols = kpObsCols[idx]
        tmpRows = kpObsRows[idx]
        tmpSeps = deltaSeps[idx]
        tmpTics = kpTics[idx]
        tmpTmags = kpTmags[idx]
        ia = np.argsort(tmpSeps)
        ia = ia[::-1]
        for iaa in ia:
            curcol = tmpCols[iaa]
            currow = tmpRows[iaa]
            curTic = tmpTics[iaa]
            midCol = int(round(curcol))
            midRow = int(round(currow))
            colX = np.arange(midCol-blkHlf, midCol+blkHlf+1)
            rowY = np.arange(midRow-blkHlf, midRow+blkHlf+1)
            sciImgCal = hdulistCal[0].data[midRow-blkHlf-1: midRow+blkHlf, midCol-blkHlf-1: midCol+blkHlf]
            gdPRF, contrastCol, contrastRow = gdPRF_calc(sciImgCal, blkHlf,
                                                         wingFAC = wingFAC, 
                                                         contrastFAC=contrastFAC)
            print('Col: {0:d} {1:5.3f} Row: {2:d} {3:5.3f} PixValue: {4} Tic: {5:d} Good? {6:b} Tmag: {7:f}'.format(
                  midCol, curcol, midRow, currow, hdulistCal[0].data[midRow-1,midCol-1], curTic, gdPRF, tmpTmags[iaa]))
            dispDS9.set('frame 1')
            dispDS9.set_np2arr(sciImgCal)
            rowSum = np.sum(sciImgCal, axis=1)
            colSum = np.sum(sciImgCal, axis=0)

            plt.subplot(121)
            plt.plot(colX, colSum, '-')
            plt.axvline(curcol, color='r')
            plt.xlabel('CCD Column Position [px]')
            plt.subplot(122)
            plt.plot(rowY, rowSum, '-')
            plt.axvline(currow, color='r')
            plt.xlabel('CCD Row Position [px]')
            plt.show()
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sector", type=int, \
                        help="TESS Observing Sector Number")
    parser.add_argument("-ca", "--camera", type=int, choices=range(1,5), \
                        help="Camera Number [1-4]")
    parser.add_argument("-cd", "--ccd", type=int, choices=range(1,5), \
                        help="CCD Number [1-4]")
    parser.add_argument("-ri", "--refimage", type=str, \
                        help="Reference Image FITS Filename With Path")
    parser.add_argument("-o", "--outputfile", type=argparse.FileType('w'),\
                        help="Control point data storeage filename with path")
    parser.add_argument("-dbg", "--debug", type=int, \
                        help="Debug level; integer higher has more output")
    parser.add_argument("-w","--wing", type=float, default = 0.9,
                         help="For a candidate WCS star, this sets a bound on the "
                         "symmetry of the PSF; the drop in flux from the core to the "
                        "brighter PSF wing cannot be smaller than this fraction "
                        "of the drop in flux form the core to the fainter PSF wing.")
    parser.add_argument("-c","--contrast", type=float, default=3.5,
                         help = "For candidate WCS star, this sets minimum ratio of the core "
                         "of the PSF to the maximum in the wings.  Crowded stars will "
                         "tend to have a smaller ratio.")
    parser.add_argument("-l", "--log", metavar="LOG_FILE", 
                        help="save logging output to thiss file")

    

    args = parser.parse_args()

# DEBUG BLOCK for hard coding input parameters and testing
#     class test_arg:
#         def __init__(self):
#             self.sector=28
#             self.camera=4
#             self.ccd=4
# #            self.refimage=open('/pdo/qlp-data/orbit-66/ffi/cam2/ccd2/FITS/tess2020268092527-00126005-2-crm-ffi-ccd2.fits','rb')
# #            self.outputfile=open('refout/ref_S29_O66_22.h5','w')
# #            self.refimage=open('/pdo/spoc-data/sector-029/ffis/tess2020241121914-s0029-4-2-0193-s_ffic.fits.gz','rb')
# #            self.outputfile=open('refout/refspocunidense_S29_42.h5','w')
#             self.refimage=open('/pdo/spoc-data/sector-028/ffis/tess2020213183915-s0028-4-4-0190-s_ffic.fits.gz','rb')
#             self.outputfile=open('refout/refspoc_S28_44.h5','w')
#             self.debug=2
#     args = test_arg()
    
    
    
    SECTOR_WANT = args.sector
    CAMERA_WANT = args.camera
    CCD_WANT = args.ccd
    REF_IMAGE = args.refimage
    
    # argsparse actually opens a file pointer. Hd5 needs filename not file pointer
    args.outputfile.close()
    outputFile = args.outputfile.name
    DEBUG_LEVEL = args.debug


    testdir = os.path.dirname(os.path.abspath(args.log))
    if not os.path.isdir(testdir):
        os.makedirs(testdir)

    tica.setup_logging(filename=args.log)    
    info_lines = tica.platform_info()
    for info_line in info_lines:
        logging.info(info_line)
    logging.info('python environment:  {}'.format(  os.getenv('CONDA_DEFAULT_ENV') )  )
    logging.info('program: {}'.format( parser.prog ) )
    logging.info('argument parameters:')
    for key in vars(args).keys():
        logging.info('     {} = {}'.format(key, vars(args)[key]) )
    starttime = time()
    
    # Call the main functionality to get
    #  bright isolated stars on a reference image
    #  that will be used on all images for a sector camera ccd 
    #  to determine a wcs solution
    get_refimg_ctrlpts(SECTOR_WANT, CAMERA_WANT, CCD_WANT, 
                       REF_IMAGE, 
                       outputFile, DEBUG_LEVEL, 
                       args.wing, args.contrast)
    runtime = time() - starttime
    logging.info("Runtime: {0:.2f}sec".format(runtime) )
