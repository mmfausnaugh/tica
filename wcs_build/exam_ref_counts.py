#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:36:51 2020

@author: cjburke
"""


import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.wcs import WCS
from astropy.io import fits
import glob
import os
import argparse
try:
    import pyds9 as pd
except ImportError:
    print('Warning: No pyds9 installed.  No debugging with image display available')
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from gwcs.wcstools import wcs_from_points
from gwcs.wcs import WCS as WCS_gwcs
from gwcs.coordinate_frames import *
import os
import glob
from step1_get_refimg_ctrlpts import gdPRF_calc, idx_filter, ring_background, flux_weighted_centroid
#import photutils.centroids as cent
from tess_stars2px import tess_stars2px_reverse_function_entry


def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list

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
    ndata = np.zeros_like(midx)
    for i in np.arange(0,nBins):
        iuse = np.where(iargs == i)[0]
        medata[i] = np.median(ydata[iuse])
        mndata[i] = np.mean(ydata[iuse])
        stddata[i] = np.std(ydata[iuse])
        ndata[i] = len(ydata[iuse])
        
    if showDetail:
        for i in np.arange(0,nBins):
            errmn = stddata[i]/np.sqrt(ndata[i])
            sigmn = mndata[i] / errmn
            print('i: {0:d} med: {1:f} mn: {2:f} n: {3:f} errmn: {4:f} sigdif: {5:f} midx: {6:f}'.format(\
                  i, medata[i], mndata[i], ndata[i], errmn, sigmn, midx[i]))
        
        
        
    return medata, midx, mndata, stddata, ndata

if __name__ == '__main__':
    refprefix = 'refout/refspocunidense_S29_'

    nTots = np.array([], dtype=np.int)
    nBrgts = np.array([], dtype=np.int)
    noises = np.array([], dtype=np.int)
    REFPIXCOL = 1024.0+45.0
    REFPIXROW = 1024.0
    PIX2DEG = 21.0/3600.0 # Turn pixels to degrees roughly
    SECTOR_WANT=29
    fitDegree=6
    colMin = 45 
    colMax = 2092
    rowMin = 1
    rowMax = 2048
    
    for iCam in range(1,5):
        for iCcd in range(1,5):
            refh5 = '{0}{1:d}{2:d}.h5'.format(refprefix,iCam, iCcd)
            CAMERA_WANT=iCam
            CCD_WANT=iCcd
            # Load the reference position information
            fin = h5py.File(refh5, 'r')    
            tics = np.array(fin['tics'])
            ras = np.array(fin['ras'])
            decs = np.array(fin['decs'])
            tmags = np.array(fin['tmags'])
            blkidxs = np.array(fin['blkidxs'])
            obscols = np.array(fin['obscols'])
            obsrows = np.array(fin['obsrows'])
            # From reference pixel coordinates get the estimated ra and dec of this point
            raproj, decproj, scinfo = tess_stars2px_reverse_function_entry(\
                     SECTOR_WANT, CAMERA_WANT, CCD_WANT, REFPIXCOL, REFPIXROW)
            proj_point = SkyCoord(raproj, decproj, frame = 'icrs', unit=(u.deg, u.deg))
            #  Reference subtracted pixel coordinates
            sclObsCols = (obscols - REFPIXCOL) * PIX2DEG
            sclObsRows = (obsrows - REFPIXROW) * PIX2DEG
    
            xy = (sclObsCols, sclObsRows)
            radec = (ras, decs)
            gwcs_obj = wcs_from_points(xy, radec, proj_point, degree=fitDegree)
            # Look for outliers to trim
            gwcs_pred_ras, gwcs_pred_decs = gwcs_obj(sclObsCols, sclObsRows)
            deg2Rad = np.pi/180.0
            deltaRas = (gwcs_pred_ras - ras) *3600.0 * np.cos(decs*deg2Rad)
            deltaDecs = (gwcs_pred_decs - decs) *3600.0
            deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
            c1 = 1.0/np.sqrt(2.0)
            idx = np.where(tmags<10.0)[0]
            std1 = np.std(deltaRas[idx])
            std2 = np.std(deltaDecs[idx])
            brightstd = np.sqrt(std1*std1+std2*std2)*c1
            idx = np.where(tmags>10.0)[0]
            std1 = np.std(deltaRas[idx])
            std2 = np.std(deltaDecs[idx])
            faintstd = np.sqrt(std1*std1+std2*std2)*c1
            std1 = np.std(deltaRas)
            std2 = np.std(deltaDecs)
            allstd = np.sqrt(std1*std1+std2*std2)*c1
            nTots = np.append(nTots, len(tics))
            nBrgts = np.append(nBrgts, len(np.where(tmags<10.0)[0]))
            noises = np.append(noises, brightstd)
            print('Cam: {0:d} Ccd:{1:d} TotN:{2:d} BrghtN:{3:d} Nose: {4:f}'.format(\
                            iCam, iCcd, len(tics), len(np.where(tmags<10.0)[0]),\
                                brightstd))
            plt.plot(obscols, obsrows, '.')
            plt.axhline(rowMin, ls='-', color='k')
            plt.axhline(rowMax, ls='-', color='k')
            plt.axvline(colMin, ls='-', color='k')
            plt.axvline(colMax, ls='-', color='k')
            plt.xlabel('Column [pix]')
            plt.ylabel('Row [pix]')
            plt.show()            
                
    plt.plot(nTots, noises, '.')
    plt.xlabel('N Reference Targs')
    plt.ylabel('Fit noise [arcsec]')
    plt.show()
    plt.plot(nBrgts, noises, '.')
    plt.xlabel('N Bright (Tm<10) Ref Targs')
    plt.ylabel('WCS Fit Resiudal [arcsec]')
    plt.show()
    
            
            
