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
from wcs_build.step2_get_refimg_ctrlpts import calc_MAD

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--skipAll", action='store_false', \
                        help="Debug level; integer higher has more output")
    args = parser.parse_args()

#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/SpocS29/Cam1/Ccd1'
#    refh5 = 'refout/refspocunidense_S29_11.h5'
#    lookatwcs = 'tess2020240141914-s0029-1-1-0193-s_ffic.fits_wcs.fits'
#    lookatimg = '/pdo/spoc-data/sector-029/ffis/tess2020240141914-s0029-1-1-0193-s_ffic.fits.gz'
#    lookatwcs = 'tess2020240142914-s0029-1-1-0193-s_ffic.fits_wcs.fits'
#    lookatimg = '/pdo/spoc-data/sector-029/ffis/tess2020240142914-s0029-1-1-0193-s_ffic.fits.gz'


#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/SpocS29/Cam1/Ccd4'
#    refh5 = 'refout/refspoc_S29_14.h5'
#    lookatwcs = 'tess2020240001914-s0029-1-4-0193-s_ffic.fits_wcs.fits'
#    lookatimg = '/pdo/spoc-data/sector-029/ffis/tess2020240001914-s0029-1-4-0193-s_ffic.fits.gz'
#    lookatwcs = 'tess2020252001914-s0029-1-4-0193-s_ffic.fits_wcs.fits'
#    lookatimg = '/pdo/spoc-data/sector-029/ffis/tess2020252001914-s0029-1-4-0193-s_ffic.fits.gz'


#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/Orbit66/Cam2/Ccd1'
#    refh5 = 'refout/ref_S29_Orbit66_21.h5'
#    lookatwcs = 'tess2020268092527-00126867-2-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd1/FITS/tess2020268092527-00126867-2-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00126509-2-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd2/FITS/tess2020268092527-00126509-2-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00127706-2-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd1/FITS/tess2020268092527-00127706-2-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00126509-2-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd2/FITS/tess2020268092527-00126509-2-crm-ffi-ccd1.fits'


#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/Orbit66/Cam2/Ccd2'
#    refh5 = 'refout/ref_S29_Orbit66_22.h5'
#    lookatwcs = 'tess2020268092527-00126867-2-crm-ffi-ccd2_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd2/FITS/tess2020268092527-00126867-2-crm-ffi-ccd2.fits'
#    lookatwcs = 'tess2020268092527-00126509-2-crm-ffi-ccd2_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd2/FITS/tess2020268092527-00126509-2-crm-ffi-ccd2.fits'

#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/Orbit66/Cam1/Ccd1'
#    refh5 = 'refout/ref_S29_Orbit66_11.h5'
#    lookatwcs = 'tess2020268092527-00126867-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00126867-1-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00126509-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00126509-1-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00127617-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00127617-1-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00127661-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00127661-1-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00127711-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00127711-1-crm-ffi-ccd1.fits'

#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/Orbit66/Cam2/Ccd3'
#    refh5 = 'refout/ref_S29_Orbit66_23.h5'
#    lookatwcs = 'tess2020268092527-00126867-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00126867-1-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00126509-2-crm-ffi-ccd3_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd3/FITS/tess2020268092527-00126509-2-crm-ffi-ccd3.fits'
#    lookatwcs = 'tess2020268092527-00126860-2-crm-ffi-ccd3_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd3/FITS/tess2020268092527-00126860-2-crm-ffi-ccd3.fits'
#    lookatwcs = 'tess2020268092527-00126848-2-crm-ffi-ccd3_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd3/FITS/tess2020268092527-00126848-2-crm-ffi-ccd3.fits'
#    lookatwcs = 'tess2020268092527-00126849-2-crm-ffi-ccd3_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam2/ccd3/FITS/tess2020268092527-00126849-2-crm-ffi-ccd3.fits'

#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/Orbit66/Cam4/Ccd1'
#    refh5 = 'refout/refstrict_S29_Orbit66_41.h5'
#    lookatwcs = 'tess2020268092527-00126867-4-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam4/ccd1/FITS/tess2020268092527-00126867-4-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00126509-4-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam4/ccd1/FITS/tess2020268092527-00126509-4-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00127455-4-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam4/ccd1/FITS/tess2020268092527-00127455-4-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00127456-4-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam4/ccd1/FITS/tess2020268092527-00127456-4-crm-ffi-ccd1.fits'

    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/SpocS29/Cam4/Ccd1'
    refh5 = 'refout/refspocunidense_S29_41.h5'
    lookatwcs = 'tess2020259212913-s0029-4-1-0193-s_ffic.fits_wcs.fits'
    lookatimg = '/pdo/spoc-data/sector-029/ffis/tess2020259212913-s0029-4-1-0193-s_ffic.fits.gz'
    lookatwcs = 'tess2020259213913-s0029-4-1-0193-s_ffic.fits_wcs.fits'
    lookatimg = '/pdo/spoc-data/sector-029/ffis/tess2020259213913-s0029-4-1-0193-s_ffic.fits.gz'


#    wcsdir = '/pdo/users/cjburke/tica/wcs_build/fitsout/Orbit66/Cam4/Ccd3'    
#    refh5 = 'refout/ref_S29_Orbit66_43.h5'
#    lookatwcs = 'tess2020268092527-00126867-1-crm-ffi-ccd1_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam1/ccd1/FITS/tess2020268092527-00126867-1-crm-ffi-ccd1.fits'
#    lookatwcs = 'tess2020268092527-00126509-4-crm-ffi-ccd3_wcs.fits'
#    lookatimg = '/pdo/qlp-data/orbit-66/ffi/cam4/ccd3/FITS/tess2020268092527-00126509-4-crm-ffi-ccd3.fits'

    args.skipAll = True

    if args.skipAll:
        # Get the wcs fits list of images
        fitsList = glob.glob(os.path.join(wcsdir,'*.fits'))
        fitsList = np.sort(fitsList)
        
        fitsNames = np.array([])
        fitsTs = np.array([], np.float)
        fitsBrghtStds = np.array([], dtype=np.float)
        fitsFaintStds = np.array([], dtype=np.float)
        fitsFracUses = np.array([], dtype=np.float)
        fitsExStd0s = np.array([], dtype=np.float)
        fitsExStd1s = np.array([], dtype=np.float)
        fitsExStd2s = np.array([], dtype=np.float)
        fitsExStd3s = np.array([], dtype=np.float)
        
        for i, curName in enumerate(fitsList):
            hdulist=fits.open(curName)
            fitsNames = np.append(fitsNames, os.path.basename(curName))
            fitsTs = np.append(fitsTs, hdulist[0].header['TIME'])
            fitsBrghtStds = np.append(fitsBrghtStds, hdulist[0].header['FITSTDB'])
            fitsFaintStds = np.append(fitsFaintStds, hdulist[0].header['FITSTDF'])
            fitsFracUses = np.append(fitsFracUses, hdulist[0].header['WCSGDF'])
            fitsExStd0s = np.append(fitsExStd0s, hdulist[0].header['FITSTDX0'])
            fitsExStd1s = np.append(fitsExStd1s, hdulist[0].header['FITSTDX1'])
            fitsExStd2s = np.append(fitsExStd2s, hdulist[0].header['FITSTDX2'])
            fitsExStd3s = np.append(fitsExStd3s, hdulist[0].header['FITSTDX3'])

        plt.plot(fitsTs, fitsFracUses, '.')
        plt.show()
        
        plt.plot(fitsTs, fitsBrghtStds, '.')
        plt.show()
        
        plt.plot(fitsTs, fitsFaintStds, '.')
        plt.show()
        
        plt.plot(fitsTs, fitsExStd0s, '.', label='0')
        plt.plot(fitsTs, fitsExStd1s, '.', label='1')
        plt.plot(fitsTs, fitsExStd2s, '.', label='2')
        plt.plot(fitsTs, fitsExStd3s, '.', label='3')
        plt.legend()
        plt.show()

    # Load the reference position information
    fin = h5py.File(refh5, 'r')    
    tics = np.array(fin['tics'])
    ras = np.array(fin['ras'])
    decs = np.array(fin['decs'])
    tmags = np.array(fin['tmags'])
    blkidxs = np.array(fin['blkidxs'])
    obscols = np.array(fin['obscols'])
    obsrows = np.array(fin['obsrows'])
    # sort by tic number
    idx = np.argsort(tics)
    tics, ras, decs, tmags, blkidxs, obscols, obsrows = idx_filter(idx, \
                tics, ras, decs, tmags, blkidxs, obscols, obsrows)

    # Load the wcs fits info for a particular image
    hdulist = fits.open(os.path.join(wcsdir, lookatwcs))
    hdr = hdulist[0].header
    fitsTIC = hdulist[1].data['TIC']
    fitsflxcols = hdulist[1].data['FLXCOL']
    fitsflxrows = hdulist[1].data['FLXROW']
    gdPRF = hdulist[1].data['FLXVALID']
    # sort by tic number
    idx = np.argsort(fitsTIC)
    fitsTIC, fitsflxcols, fitsflxrows, gdPRF = idx_filter(idx, \
                    fitsTIC, fitsflxcols, fitsflxrows, gdPRF)
    # Double check that the data rows are lined up by tic
    checkPASS = True
    if not len(fitsTIC) == len(tics):
        checkPASS = False
    if checkPASS:
        sumd = np.sum(fitsTIC - tics)
        if sumd != 0:
            checkPASS = False
    if not checkPASS:
        print('Error: TIC arrays disagree between reference data and image')
        exit()
    idxgd = np.where(gdPRF == 1)[0]
    allBad = False
    if len(idxgd) == 0:
        # revert to looking at stars used in bad regions
        print('WARNING NO Good targets.')
        idxgd = np.where(np.abs(gdPRF) == 1)[0]
        allBad = True

    #Read WCS and see how well it recovers the input ra and decs   
    # from their flux weighted centroid measurements
    my_wcs = WCS(hdr)
    pix_list = list(zip(fitsflxcols[idxgd], fitsflxrows[idxgd]))
    pix_coords = np.array(pix_list, dtype=np.double)
    pred_coords = my_wcs.all_pix2world(pix_coords, 1)
    predRa = np.array([x[0] for x in pred_coords], dtype=np.double)
    predDec = np.array([x[1] for x in pred_coords], dtype=np.double)
    deg2Rad = np.pi/180.0
    deltaRas = (predRa - ras[idxgd]) *3600.0 * np.cos(predDec*deg2Rad)
    deltaDecs = (predDec - decs[idxgd]) *3600.0
    deltaSeps = np.sqrt(deltaRas*deltaRas + deltaDecs*deltaDecs)
    idxb = np.where(tmags[idxgd]<10.0)[0]
    brightstd = np.std(deltaSeps[idxb])
    idxf = np.where(tmags[idxgd]>10.0)[0]
    faintstd = np.std(deltaSeps[idxf])
    allstd = np.std(deltaSeps)

    # print some stats of numbers of bright and faint targets valid
    #  over the sub regions
    nCol = hdulist[0].header['CTRPCOL']
    nRow = hdulist[0].header['CTRPROW']
    nBlck = nCol*nRow
    brightns = np.zeros((nBlck,), dtype=np.int32)
    totns = np.zeros((nBlck,), dtype=np.int32)
    brightstds = np.zeros((nBlck,), dtype=np.float)
    faintstds = np.zeros((nBlck,), dtype=np.float)
    tmptmags = tmags[idxgd]
    for i in range(nBlck):
        jRow = i//nCol
        iCol = np.mod(i, nCol)
        ia = np.where(blkidxs[idxgd] == i)[0]
        if len(ia)>0:
            totns[i] = len(tmptmags[ia])
            idxt = np.where(tmptmags[ia]<10.0)[0]
            if len(idxt)>0:
                brightns[i] = len(idxt)
            else:
                brightns[i] = 0
        else:
            totns[i] = 0
            
        print('iCol: {0:d} jRow: {1:d} nTot: {2:d} nBrght {3:d}'.format(\
                    iCol, jRow, totns[i], brightns[i]))
    print('totN:{0:d} brightN:{1:d}'.format(np.sum(totns), np.sum(brightns)))
    fitsBrghtStds = hdulist[0].header['FITSTDB']
    fitsFaintStds =  hdulist[0].header['FITSTDF']
    fitsFracUses = hdulist[0].header['WCSGDF']
    print('brghtstd: {0:f} faintstd:{1:f} frac:{2:f}'.format(fitsBrghtStds, fitsFaintStds, fitsFracUses))

    # Display image with the centroids and their quality
    #  with ds9
    DEBUG_DS9 = True
    try:
        dispTargs = pd.ds9_targets()
        dispDS9 = pd.DS9(dispTargs[0].split()[1])
    except:
        DEBUG_DS9 = False
    # Open image
    hdulistCal = fits.open(lookatimg)
    # QLPor SPOC
    try:
        # QLP has MIDTJD
        timeVal = hdulistCal[0].header['MIDTJD']
        timeKey = 'MIDTJD'
        dataKey = 0
    except:
        # SPOC has TSTART
        timeKey = 'TSTART'
        dataKey = 1

    sciImgCal = hdulistCal[dataKey].data
    dispDS9.set('frame 1')
    dispDS9.set_np2arr(sciImgCal)
    for j, curCol in enumerate(fitsflxcols):
        curRow = fitsflxrows[j]
        curGd = gdPRF[j]
        if curGd == 1:
            dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 8 color=green'.format(curCol, curRow))
        if curGd == 0:
            dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 8 color=orange'.format(curCol, curRow))
        if curGd == -1:
            dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 8 color=cyan'.format(curCol, curRow))
        if curGd == -2:
            dispDS9.set('regions', 'image; point({0:7.3f},{1:7.3f}) # point=x 8 color=red'.format(curCol, curRow))


    # Show stars and their quality on detector with pyplot
    plt.plot(fitsflxcols, fitsflxrows, '.')
    idx = np.where(gdPRF == 0)[0]
    plt.plot(fitsflxcols[idx], fitsflxrows[idx], '.')
    idx = np.where(gdPRF == -1)[0]
    plt.plot(fitsflxcols[idx], fitsflxrows[idx], '.k')
    idx = np.where(gdPRF == -2)[0]
    plt.plot(fitsflxcols[idx], fitsflxrows[idx], '.r')
    plt.show()

    # Show stars broken down by bright and faint on detector
    plt.plot(fitsflxcols[idxgd][idxf], fitsflxrows[idxgd][idxf], '.')
    plt.plot(fitsflxcols[idxgd][idxb], fitsflxrows[idxgd][idxb], '.')
    plt.show()
    
    
    showDetailBool = False
    plt.plot(tmags[idxgd], deltaRas, '.')
    meddata, midx, mndata, stddata, ndata = binmedian(tmags[idxgd], deltaRas, showDetail=showDetailBool)
    plt.plot(midx, meddata, '-')
    plt.xlabel('Tmag')
    plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
    plt.axhline(0.0, ls='--', color='r')
    #plt.ylim([-21.0, 21.0])
    plt.show()
    plt.plot(tmags[idxgd], deltaDecs, '.')
    meddata, midx, mndata, stddata, ndata = binmedian(tmags[idxgd], deltaDecs, showDetail=showDetailBool)
    plt.plot(midx, meddata, '-')
    plt.xlabel('Tmag')
    plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
    plt.axhline(0.0, ls='--', color='r')
    #plt.ylim([-21.0, 21.0])
    plt.show()
    plt.plot(fitsflxcols[idxgd], deltaRas, '.')
    meddata, midx, mndata, stddata, ndata = binmedian(fitsflxcols[idxgd], deltaRas, showDetail=showDetailBool)
    plt.plot(midx, meddata, '-')
    plt.xlabel('Column [px]')
    plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
    plt.axhline(0.0, ls='--', color='r')
    #plt.ylim([-21.0, 21.0])
    plt.show()
    plt.plot(fitsflxrows[idxgd], deltaRas, '.')
    meddata, midx, mndata, stddata, ndata = binmedian(fitsflxrows[idxgd], deltaRas, showDetail=showDetailBool)
    plt.plot(midx, meddata, '-')
    plt.xlabel('Row [px]')
    plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
    plt.axhline(0.0, ls='--', color='r')
    #plt.ylim([-21.0, 21.0])
    plt.show()
    plt.plot(fitsflxcols[idxgd], deltaDecs, '.')
    meddata, midx, mndata, stddata, ndata = binmedian(fitsflxcols[idxgd], deltaDecs, showDetail=showDetailBool)
    plt.plot(midx, meddata, '-')
    plt.xlabel('Column [px]')
    plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
    plt.axhline(0.0, ls='--', color='r')
    #plt.ylim([-21.0, 21.0])
    plt.show()
    plt.plot(fitsflxrows[idxgd], deltaDecs, '.')
    meddata, midx, mndata, stddata, ndata = binmedian(fitsflxrows[idxgd], deltaDecs, showDetail=showDetailBool)
    plt.plot(midx, meddata, '-')
    plt.xlabel('Row [px]')
    plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
    plt.axhline(0.0, ls='--', color='r')
    #plt.ylim([-21.0, 21.0])
    plt.show()
    
    
        # Record wcs fit residuals for subregions
    # Doing a mapping where targets in several subregions in the corners are
    # combined to get teh fit residual
    nResids = 4
    residSZ = 3
    #  With 5x5 grid these are corner subregions
    useRR = [[0,1,5], [3,4,9], [15,20,21], [23,24,19]]
    # Also show residualsfrom targets over several regions
    tm = tmags[idxgd]
    blks = blkidxs[idxgd]
    for i in range(nResids):
        idxAll = np.array([], dtype=np.int64)
        for j in range(residSZ):
            idx = np.where((tm<10.0) & (blks == useRR[i][j]))[0]
            idxAll = np.append(idxAll, idx)
        c1 = 1.0/np.sqrt(2.0)
        std1 = np.std(deltaRas[idxAll])
        std2 = np.std(deltaDecs[idxAll])
        usestd= np.sqrt(std1*std1+std2*std2)*c1
        mad1 = calc_MAD(deltaRas[idxAll])
        mad2 = calc_MAD(deltaDecs[idxAll])
        usemad= np.sqrt(mad1*mad1+mad2*mad2)*c1

        print('Region: {0:d} std:{1:f} mad:{2:f}'.format(i, usestd, usemad))
        plt.plot(tm[idxAll], deltaRas[idxAll], '.')
        plt.xlabel('Tmag')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.show()
        plt.plot(tm[idxAll], deltaDecs[idxAll], '.')
        plt.xlabel('Tmag')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.show()
        plt.plot(fitsflxcols[idxgd[idxAll]], deltaRas[idxAll], '.')
        plt.xlabel('Column [px]')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.show()
        plt.plot(fitsflxrows[idxgd[idxAll]], deltaRas[idxAll], '.')
        plt.xlabel('Row [px]')
        plt.ylabel('GWCS Predicted - Observed RA Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.show()
        plt.plot(fitsflxcols[idxgd[idxAll]], deltaDecs[idxAll], '.')
        plt.xlabel('Column [px]')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.show()
        plt.plot(fitsflxrows[idxgd[idxAll]], deltaDecs[idxAll], '.')
        plt.xlabel('Row [px]')
        plt.ylabel('GWCS Predicted - Observed Declination Position [arcsec]')
        plt.axhline(0.0, ls='--', color='r')
        plt.show()
      

        
    
