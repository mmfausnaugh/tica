#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 16:18:07 2020

@author: cjburke
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import argparse
import os
try:
    import pyds9 as pd
except ImportError:
    print('Warning: No pyds9 installed.  No debugging with image display available')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-id", "--imagedir", \
                        help="Directory name containing images")
    parser.add_argument("-ip", "--imageprefix", \
                        help="Partial image prefix -- DO NOT INCLUDE WILDCARD; glob.glob routine will add * after prefix")
    parser.add_argument('-s', "--subsample", type=int,\
                        help="Show background every subsample image",\
                        default=10)
    parser.add_argument("-o", "--outputfile", type=argparse.FileType('w'),\
                        help="Output Image filename")
    parser.add_argument('-n', "--nimagein", type=int,\
                        help="Show the cadence numbers for nImageIn frames at beginning and end of data",\
                        default=300)
    parser.add_argument("-dbg", "--debug", type=int, \
                        help="Debug level; integer higher has more output",\
                        default=0)
    args = parser.parse_args()
    # argsparse actually opens a file pointer. Needs filename not file pointer
    args.outputfile.close()
    outputFile = args.outputfile.name

    DEBUG_LEVEL = args.debug
    DEBUG_DS9 = False
    if DEBUG_LEVEL>1:
        try:
            dispTargs = pd.ds9_targets()
            dispDS9 = pd.DS9(dispTargs[0].split()[1])
            DEBUG_DS9= True
        except:
            DEBUG_DS9 = False
            

    
    # Get the image list
    fileList = glob.glob(os.path.join(args.imagedir,\
                                      '{0}*'.format(args.imageprefix)))
    # Sort fileList to ensure images are in time order
    fileList = np.sort(fileList)
    nImgs = len(fileList)
    idxUse = np.arange(0,nImgs, args.subsample)
    nUse = len(idxUse)
    # Create storage for times and median background level
    xTime = np.zeros((nUse,), dtype=np.float)
    xCad = np.zeros((nUse,), dtype=np.int)
    yMed = np.zeros((nUse,), dtype=np.float)

    cnt = 0
    for ii, iu in enumerate(idxUse):
        curFileName = fileList[iu]    
        try:
            hdulist = fits.open(curFileName)
        except:
            print('Problem reading file {0}'.format(curFileName))
        prihdr = hdulist[0].header
        shp = hdulist[0].data.shape
        xTime[ii] = prihdr['MIDTJD']
        xCad[ii] = prihdr['CADENCE']
        
        numOrigCols = shp[1]
        oCols = np.arange(1, numOrigCols+1)
        # Below subtract 1 from pixel coordinates because python stars at zero
        # Doing Channel C
        sciColS = 1069 - 1
        sciColE = 1580 - 1
        sciRowS = 1024-1
        sciRowE = 1224-1

        oCols = oCols[sciColS:sciColE+1]
        sciImg = hdulist[0].data[sciRowS:sciRowE+1,sciColS:sciColE+1].T
        if DEBUG_DS9:
            dispDS9.set('frame 1')
            dispDS9.set_np2arr(sciImg.T)
            xtmp= np.arange(-2.0, 3.0)
            plt.plot(xtmp, xtmp*xtmp, '.')
            plt.show()

        medSci = np.median(sciImg.astype(float).flatten())
        yMed[ii] = medSci

        if np.mod(ii, 10) == 0:
            print('Done {0:d} of {1:d}'.format(ii, nUse))
    
    xTime = xTime + prihdr['TJD_ZERO']
    # Print out the filename for the nimagein at either end
    print('File {0:d} at beginning {1}'.format(args.nimagein, fileList[args.nimagein]))
    print('File {0:d} at end {1}'.format(args.nimagein, fileList[-args.nimagein]))

    plt.plot(xCad, yMed, '.')
    plt.ylim([20000, 100000])
    plt.xlabel('Cadence Number')
    plt.ylabel('Counts [e]')
    plt.savefig(outputFile)
    plt.show()
    
