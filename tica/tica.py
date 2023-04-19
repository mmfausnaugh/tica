#
#Copyright (c) 2021 Michael Fausnaugh
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the Software), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


import numpy as np
import os.path
import sys
import logging
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
plt.rcParams['xtick.labelsize']='smaller'
plt.rcParams['ytick.labelsize']='smaller'

from astropy.io import fits
from astropy.table import Table

import datetime

from .GainModel import GainModel
from .LinearityModels import LinearityModel, SpocLinearity_Updated

fstem = os.path.abspath(os.path.dirname(__file__) + '/../')
with open(os.path.join(fstem, 'VERSION'),'r') as infile:
    version = infile.read().strip()


# these are functions for fast manipulation of images
def sum_frames(framelist):
    #    for f in framelist:  f.reassemble()
    imarray = [ frame.full_image for frame in framelist ] 
    sumframe = np.sum(imarray, axis = 0)
    return FitsFrame(sumframe)

def subtract_frames(frame1, frame2):
    #    frame1.reassemble()
    #    frame2.reassemble()
    subtractframe = frame1.full_image - frame2.full_image
    return FitsFrame(subtractframe)        

def mean_frames(framelist):
    #    for f in framelist:  f.reassemble()
    imarray = [ frame.full_image for frame in framelist ] 
    meanframe = np.mean( imarray, axis=0)
    return FitsFrame(meanframe)        

def median_frames(framelist):
    #    for f in framelist:  f.reassemble()
    imarray = [ frame.full_image for frame in framelist ]  
    medianframe = np.median( imarray,axis=0)
    return FitsFrame(medianframe)        


class Sector(object):
    """
    Keep track of individual CCD output channels, i.e., sectors
    """
    def __init__(self, underclock, science, overclock):
        self.underclock = underclock.astype(float)
        self.science = science.astype(float)
        self.overclock = overclock.astype(float)

    def rebin(self, ncombine, axis, region):
        """attr is the variable for which we will rebin (along rows)"""
        ncombine = int(ncombine)
        duse = getattr(self, region)
        maxindex = duse.shape[ axis ]
        assert ncombine <= maxindex
        
        iuse = 0
        rebin = []
        while iuse < maxindex:
            if iuse + ncombine >= maxindex:
                new = duse[ iuse : : , :]
            else:
                new = duse[ iuse : iuse + ncombine , :]
                    
            rebin.append(np.mean(new,axis=axis))
            iuse += ncombine

        if np.shape(rebin)[0] == 1:
            return rebin[0]
        else:
            return rebin        


class Calibration(object):
    """Class to handle initialization of calibration models so it is
    performed just once

    """

    @staticmethod
    def _set_calibration_data(cal_file, cal_dir, logger):
        """given calibration files and directory return the FITs data
        structure or else return none if no file exists to allow for
        placeholder

        """
        fits_file = os.path.join(cal_dir, cal_file)
        if os.path.exists(fits_file):
            with fits.open(fits_file) as fits_handle:
                return fits_handle[0].data
        else:
            logger.debug("FITS calibration file not found: {}".format( fits_file ) )
            return None

    def __init__(self, calibration_dir=None):

        self.logger = logging.getLogger(__name__)
        
        if not calibration_dir:
            top_level_path=os.path.abspath(os.path.dirname(__file__))
            calibration_dir = os.path.abspath(os.path.join(top_level_path, "../calibration/"))
        
        int_time_key = calibration_dir.split('_')[-1]
        if '30min' in int_time_key:
            self.int_time = 1800
        elif '10min' in int_time_key:
            self.int_time = 600
        elif '02min' in int_time_key:
            self.int_time = 120
        elif '20sec' in int_time_key:
            self.int_time = 20
        elif '200sec' in int_time_key:
            self.int_time = 200

        twodbias_dir = os.path.join(calibration_dir, 'twodbias/')
        twodbias_files = {'cam1':{'ccd1':'cam1_ccd1_twoDcorrect.fits',
                                    'ccd2':'cam1_ccd2_twoDcorrect.fits',
                                    'ccd3':'cam1_ccd3_twoDcorrect.fits',
                                    'ccd4':'cam1_ccd4_twoDcorrect.fits'} ,
                            'cam2':{'ccd1':'cam2_ccd1_twoDcorrect.fits',
                                    'ccd2':'cam2_ccd2_twoDcorrect.fits',
                                    'ccd3':'cam2_ccd3_twoDcorrect.fits',
                                    'ccd4':'cam2_ccd4_twoDcorrect.fits'} ,
                            'cam3':{'ccd1':'cam3_ccd1_twoDcorrect.fits',
                                    'ccd2':'cam3_ccd2_twoDcorrect.fits',
                                    'ccd3':'cam3_ccd3_twoDcorrect.fits',
                                    'ccd4':'cam3_ccd4_twoDcorrect.fits'} ,
                            'cam4':{'ccd1':'cam4_ccd1_twoDcorrect.fits',
                                    'ccd2':'cam4_ccd2_twoDcorrect.fits',
                                    'ccd3':'cam4_ccd3_twoDcorrect.fits',
                                    'ccd4':'cam4_ccd4_twoDcorrect.fits'} }

        flatfields_dir = os.path.join(calibration_dir, 'flatfields/')
        flatfields_files = {'cam1':{'ccd1':'cam1_ccd1_image.fits',
                                    'ccd2':'cam1_ccd2_image.fits',
                                    'ccd3':'cam1_ccd3_image.fits',
                                    'ccd4':'cam1_ccd4_image.fits'} ,
                            'cam2':{'ccd1':'cam2_ccd1_image.fits',
                                    'ccd2':'cam2_ccd2_image.fits',
                                    'ccd3':'cam2_ccd3_image.fits',
                                    'ccd4':'cam2_ccd4_image.fits'} ,
                            'cam3':{'ccd1':'cam3_ccd1_image.fits',
                                    'ccd2':'cam3_ccd2_image.fits',
                                    'ccd3':'cam3_ccd3_image.fits',
                                    'ccd4':'cam3_ccd4_image.fits'} ,
                            'cam4':{'ccd1':'cam4_ccd1_image.fits',
                                    'ccd2':'cam4_ccd2_image.fits',
                                    'ccd3':'cam4_ccd3_image.fits',
                                    'ccd4':'cam4_ccd4_image.fits'} }
        
        # reconstruct similar dictionaries except that we use fits to load file
        self.twodbias = {camera: {ccd: self._set_calibration_data(calfile, twodbias_dir, self.logger) for ccd, calfile in outervalue.items()}
                           for camera, outervalue in twodbias_files.items()}
    
        self.flatfields = {camera: {ccd: self._set_calibration_data(calfile, flatfields_dir, self.logger) for ccd, calfile in outervalue.items()}
                           for camera, outervalue in flatfields_files.items()}

    def get_twodbias(self, camuse, ccduse):
        return self.twodbias['cam' + str(camuse)]['ccd' + str(ccduse)]

    def get_flatfields(self, camuse, ccduse):
        return self.flatfields['cam' + str(camuse)]['ccd' + str(ccduse)]
        
class CCD(object):
    """Splits a CCD into its sectors. If the CCD is passed in inverted
    (this is the case for CCDs 1 and 2), rotates it by 180 degrees
    before proceeding.

    self.sectors --- a list of Sector objects (ABCD)

    methods:
    get_image():  return concat. column with science part of sectors
    get_overclock():  return concat. column with overclock part of sectors
    get_frame():  return concat. columns of all pixels

    """
    def __init__(self, arr, calibration=None, inverted=False):
        # get logger
        self.logger = logging.getLogger(__name__)

        # if no calibration object provided, create one
        if calibration:
            self.calibration = calibration
        else:
            self.calibration = Calibration()
            
        self.ccdnum = None
        self.camnum = None

        if arr.shape != (2078, 2136):
            raise ValueError("CCD shape must be (2078, 2136)")

        if inverted:
            self.arr = arr[::-1, ::-1]
        else:
            self.arr = arr
        A = Sector(self.arr[:, 0:11],  self.arr[:, 44:556],    self.arr[:, 2092:2103])
        B = Sector(self.arr[:, 11:22], self.arr[:, 556:1068],  self.arr[:, 2103:2114])
        C = Sector(self.arr[:, 22:33], self.arr[:, 1068:1580], self.arr[:, 2114:2125])
        D = Sector(self.arr[:, 33:44], self.arr[:, 1580:2092], self.arr[:, 2125:2136])
        self.sectors = np.asarray([A, B, C, D])

        self.gains = {'A':1.0,'B':1.0,'C':1.0,'D':1.0}
        #params for linearity correction
        self.gplus  = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
        self.gminus = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}

        self.spoc_linearity = SpocLinearity_Updated()
        self.coadds = 1
        

    #===========================================
    # Convenience, return copies of CCD sections
    #===========================================
        
    def get_image(self):
        """Returns a ndarry of the science pixels
        Expected to be used when using tica in 'library' mode
        """
        return np.c_[
            self.sectors[0].science,
            self.sectors[1].science,
            self.sectors[2].science,
            self.sectors[3].science     ]

    def get_overclock(self):
        """Return a 2d array of overclock pixels
        Expected to be used when using tica in 'library' mode
        """
        return np.c_[
            self.sectors[0].overclock,
            self.sectors[1].overclock,
            self.sectors[2].overclock,
            self.sectors[3].overclock      ]

    def get_frame(self):
        """Return a 2d array of the full image
        Expected to be used when using tica in 'library' mode
        """
        out = []
        A = self.sectors[0]
        B = self.sectors[1]
        C = self.sectors[2]
        D = self.sectors[3]
        return np.c_[  A.underclock, B.underclock, C.underclock, D.underclock,
                       A.science, B.science,  C.science,  D.science,
                       A.overclock,   B.overclock,    C.overclock,   D.overclock   ]

    #==============================
    # Workhorse calibration method
    #=============================
    def calibrate(self):
        if self.ccdnum == None:
            raise Exception("Must specify CCD number to calibrate!")
        if self.camnum == None:
            raise Exception("Must specify camera number to calibrate!")

        #here, it's a CCD object
        out = CCD(self.get_frame())

        ###Apply 2D bias correction
        self.logger.debug("doing 2D correction on camera: %d and CCD: %d" % (self.camnum, self.ccdnum))
        #still a CCD object
        out = CCD(self.twod_correct(out.get_frame(), self.camnum, self.ccdnum))

        ###Apply 1D bias correction (dehoc)
        self.logger.debug("doing overclock correction on camera: %d and CCD: %d" % (self.camnum, self.ccdnum))
        #still a CCD object
        out = out.overclock_correct()


        #Convert to electrons
        #still a CCD object
        out = self.convert_to_electrons(out)
        #linearity correction
        #still a CCD object

        #version of POC model
        #out = self.linearity_correct(out)

        #in the spoc version, there is an object that stores the
        #coefficients and knows how to do the math
        out = self.spoc_linearity_correct(out)

        ###Apply smear correction
        self.logger.debug("doing smear correction on camera: %d and CCD: %d" % (self.camnum, self.ccdnum))
        #now an Ndarray
        out = self.smear_correct(out.get_frame())


        ###Apply flatfield
        self.logger.debug("doing flatfield correction on camera: %d and CCD: %d" % (self.camnum, self.ccdnum))
        #now a CCD object
        #this overwrite the gains!
        out  = CCD(self.flatfield_correct(out, self.camnum, self.ccdnum))
        out.gains = self.gains
        #print(self.gains)

        return out
    
    #===========================
    # Algorithms for calibration
    #===========================

    def twod_correct(self, out, camuse, ccduse):
        caldata = self.calibration.get_twodbias(camuse, ccduse)
        return out - caldata

    def flatfield_correct(self, out, camuse, ccduse):
        caldata = self.calibration.get_flatfields(camuse, ccduse)
        return out/caldata

    def smear_correct(self, out):
        #smear = np.mean( out[2061:2068 ,:] ,axis=0)
        smear = np.median( out[2061:2068 ,:] ,axis=0)

        #weighted mean.  remeber that error = sqrt(cts),
        #but weights = 1/error^2 = 1/cts
        #smear = 1./np.sum( (1./out[2061:2069 ,:]) ,axis=0)
        return out - smear
    
    def overclock_correct(self):
        #first few columns can be affected by scattered light
        mask_col = slice(3,  None)
        #first few hundred rows can be affected by start of 
        #frame ringing. note that the 2D bias model is slightly different than
        #the observed fram rining---overall, there is a small residual
        #at low row relative to the SPOC FFIs.
        mask_row = slice(750,None)

        cal = CCD(np.zeros((2078,2136)))
        for i in range(len(self.sectors)):
            correct = np.mean(self.sectors[i].overclock[mask_row, mask_col] )
            cal.sectors[i].underclock = self.sectors[i].underclock - correct
            cal.sectors[i].science      = self.sectors[i].science - correct
            cal.sectors[i].overclock   = self.sectors[i].overclock - correct
        return cal
    
    def convert_to_electrons(self, out):
        """default is gain == 1, unless it is CCD file with cam/ccd info in the
        header
        
        Note that the self.gains is applied on input data.  So the
        class knows about the correct value, but never applies it


        For now, this does the calculation in place!  The functional
        aspect is handled in CCD.calibrate()

        """
        guse = np.asarray([self.gains['A'], self.gains['B'], self.gains['C'], self.gains['D'] ])
        for i,s in enumerate(out.sectors):
            s.underclock *= guse[i]
            s.science *= guse[i]
            s.overclock *= guse[i]
        return out

    def linearity_correct(self, out):
        """Uses POC model, which is a kind of sigmoid.
       
        See IH, 6.4.2:

        e_out = e_in/((1 + e_in*gplus)*(1 - e_in*gminus))

        where gplus = gainGain_per_electron and gminus =
        gainLoss_per_electron are measured parameters, stored in the
        LinearityModel Class.

        e_in/e_out are in units of electrons/2 second exposure

        Do not apply to under/over clocks

        For now, this does the calculation on the original data
        structure!  The functional aspect is handled in CCD.calibrate()

        """
        gplus  = np.asarray([self.gplus['A'],  self.gplus['B'],  self.gplus['C'],  self.gplus['D'] ])
        gminus = np.asarray([self.gminus['A'], self.gminus['B'], self.gminus['C'], self.gminus['D']])
        for i,s in enumerate(out.sectors):
#            print self.camnum, self.ccdnum, i,  gplus[i],gminus[i], np.mean(( (1 + s.science/720.0*gplus[i])*(1 - s.science/720.0*gminus[i]) ))
            s.science = s.science/( (1 + s.science/720.0/0.99*gplus[i])*(1 - s.science/720.0/0.99*gminus[i]) )
        return out

    def spoc_linearity_correct(self, out):
        """Uses spoc parameterization, which is stored in a class
        LinearityModels.SpocLinearity. 

        Assign the class in this object's __init__

        Do not apply to under/over clocks

        For now, this does the calculation on the original data
        struct!  The functional aspect is handled in CCD.calibrate()

        """
        labels = ['A','B','C','D']
        for i,s in enumerate(out.sectors):
            #            print self.camnum, self.ccdnum, i,  gplus[i],gminus[i], np.mean(( (1 + s.science/720.0*gplus[i])*(1 - s.science/720.0*gminus[i]) ))

            s.science = self.spoc_linearity.linearity_correct(s.science, 
                                                              self.camnum, self.ccdnum , 
                                                              labels[i], 
                                                              self.coadds)
            
        return out

    
    #=======================================
    # Get diagnostic info and generate plots
    #=======================================
    
    def report(self, region, save=False,label=None, plot =True):
        """
        Get 3 plots per sector, organize output with one figure per CCD
        (12 plots per figure), and table of statistics
        """

        sector_labels=[ 'A', 'B', 'C', 'D' ]
        
        ndata0, ndata1 = getattr(self.sectors[0], region).shape
        collapsed_rows = self.project(ndata0, 0, region)
        collapsed_cols  = self.project(ndata1, 1, region)
        
        means = []
        stds  = []
            

        for j, sector in enumerate(self.sectors):
            duse = getattr(self.sectors[j], region)
            means.append(   duse.mean() )
            stds.append(   duse.std() )

        if plot:
            F,( axes ) = plt.subplots(4,3)
            F.set_size_inches(11,8.5)

            for j, sector in enumerate(self.sectors):
                duse = getattr(self.sectors[j], region)            
                
                axes[j][0].plot(  collapsed_rows[j] )
                axes[j][1].plot(  collapsed_cols[j]  )
                axes[j][2].hist(
                    np.ravel(duse), 50)

                #                axes[j][0].set_xlim([-0.5, 11])
                #                axes[j][1].set_xlim([-50, 2150])
                #                axes[j][0].set_xticklabels([])
                #                axes[j][1].set_xticklabels([])

                #                axes[j][0].yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
                #               axes[j][1].yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
                #              axes[j][2].yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
                
                axes[j][0].set_ylabel('Sector {sector_label}'.format(sector_label= sector_labels[j])  )


                axes[j][0].set_xlabel('col')

                #            axes[j][0].set_xticklabels(['0','5','10'])
                #            axes[j][1].set_xticklabels(['0','0','1000','2000'])

                axes[-1][1].set_xlabel('row')
                axes[-1][2].set_xlabel('counts')
                
                axes[0][0].set_title('Mean of Rows')
                axes[0][1].set_title('Mean of Columns')
            
            F.suptitle('{region} region'.format(region = region))



            if save:
                if label is not None:
                    plt.savefig(label + '_ccd' + '.pdf',fmt='pdf')
                else:
                    plt.savefig('ccd.pdf',fmt='pdf')

        table = Table(data=[sector_labels, means,stds], names=('Sector', 'Mean','Std') )
        table['Mean'].format='%.2f'
        table['Std'].format='%.2f'
        #table.pprint(max_lines=100)
        return sector_labels, means, stds


    #=========================
    # Convenience for reports
    #=========================
    def project(self, ncombine, axis, region):
        """
        Just a way of looping through all sectors, rebinning along the given axis,
        and returning in an organized list of lists
        """
        ccd_out = []
        for sector in self.sectors:                
            ccd_out.append( sector.rebin(ncombine, axis, region) )
        return np.asarray(ccd_out)


    def get_image_mode(self,bins=None):
        pixels = np.ravel( self.get_image() )
        if bins is None:
            values, bins = np.histogram( pixels, bins = 500 )
        else:
            values,bins = np.histogram(pixels,bins=bins)

        mode = bins[ values == values.max() ][0]
        return mode



    
class CCD_File(CCD):
    """
    Handles reading, writing, parsing headers, in a way unique to TESS
    frames, and assuming one CCD at a time

    'hdu' refers to astropy's 'Header Data Units'
    """
    
    def __init__(self, fname, calibration=None, outdir="."):

        self.fname = fname
        self.outdir = outdir
        
        self.hdu = fits.open(fname, mode='readonly')

        #raw TSO file
        if len(self.hdu) == 1:
            self.header = self.hdu[0].header
            image  = self.hdu[0].data
            
        #calibrated tica file
        elif len(self.hdu) == 2 and self.hdu[1].header['XTENSION'] == 'BINTABLE':
            self.header = self.hdu[0].header
            image  = self.hdu[0].data

        #handles either spoc raw or cal---maybe more checks are needed??
        #yes!, because tica calibrated FFIs now have 2 headers

        #SPOC file; raw has 2 extensions, calibrated ahs 3
        elif len(self.hdu) == 2 or len(self.hdu) == 3:
            assert self.hdu[0].data == None

            self.header = fits.Header()

            for key in self.hdu[0].header.keys():

                if key not in ['SIMPLE','BITPIX','NAXIS',
                               'EXTEND','NEXTEND',
                               'EXTNAME','EXTVER',
                               'CHECKSUM',
                               'ORIGIN','DATE','PROCVER',
                               'FILEVER','IMAGTYPE','CREATOR',
                               'CHECKSUM']:
                    self.header.set(key,
                                    self.hdu[0].header[key],
                                    self.hdu[0].header.comments[key])
                    
                if key in ['ORIGIN','DATE','PROCVER',
                           'FILEVER','IMAGTYPE','CREATOR']:
                    self.header['HISTORY'] = 'Original SPOC header {} = {}'.format(
                                        key,self.hdu[0].header[key] )

            for key in self.hdu[1].header.keys():
                if key not in ['XTENSION','BITPIX','NAXIS',
                               'NAXIS1','EXTEND','BSCALE','BZERO',
                               'INHERIT','EXTNAME','EXTVER',
                               'ORIGIN','DATE','PROCVER',
                               'FILEVER','IMAGTYPE','CREATOR',
                               'CHECKSUM']:
                    self.header.set(key,
                                    self.hdu[1].header[key],
                                    self.hdu[1].header.comments[key])

                if key in ['ORIGIN','DATE','PROCVER',
                           'FILEVER','IMAGTYPE','CREATOR']:
                    self.header['HISTORY'] = 'Original SPOC header {} = {}'.format(
                                        key,self.hdu[1].header[key] )


                    
            image  = self.hdu[1].data

        else:
            raise IOError("must be a TSO or SPOC FFI")


            
        try:
            #check that the calibration model matches the exposure time of the data.
            #SPOC keyword
            assert np.isclose(calibration.int_time,
                              self.header['EXPOSURE']*86400/0.8/0.99), \
                "It appears you are trying to calibrate {:4.0f} second data with {} second models".format(
                self.header['EXPOSURE']*86400/0.8/0.99,
                calibration.int_time)
        except KeyError:
            #TSO keyword
            assert calibration.int_time == self.header['INT_TIME'], \
                "It appears you are trying to calibrate {:4.0f} second data with {} second models".format(
                self.header['INT_TIME'],
                calibration.int_time)

        super(CCD_File, self).__init__(image, calibration=calibration)

        #SPOC uses 'CAMERA' instead of 'CAM' and 'CCDNUM' instead of 'CCD'
        #header.get returns None as a backup if it can't match the keyword
        self.ccdnum  = self.header.get('CCD')
        self.camnum = self.header.get('CAM')
        if self.camnum is None:
            self.camnum = self.header.get('CAMERA')
        if self.camnum is None:
            self.camnum = self.header.get('CAMNUM')
        if self.camnum is None:
            raise ValueError('Did not find camera number in the FITS header')

        if self.ccdnum == None:
            self.ccdnum = self.header.get('CCDNUM')


        self.gains = GainModel.gains['cam' + str(self.camnum)]['ccd' + str(self.ccdnum)]
        #for linearity
        self.gplus  = LinearityModel.gain_gain['cam' + str(self.camnum)]['ccd' + str(self.ccdnum)]
        self.gminus = LinearityModel.gain_loss['cam' + str(self.camnum)]['ccd' + str(self.ccdnum)]

        try:
            #SPOC keword
            self.coadds = self.header['NREADOUT']
        except:
            #TSO keyword
            self.coadds = self.header['INT_TIME']*0.8/2


        try:
            self.calibrated_frame = self.calibrate()
        except:
            print('calibration failed')
            raise

    def write_calibrate(self, no_gzip=False):      
        hdu_out = fits.PrimaryHDU(self.calibrated_frame.get_frame().astype(np.float32))
        for key in self.header.keys():                                                     
            #censor the standard fits headers, which have to be changed.                      
            #Let astropy handle internally                                                    
            if key not in ['SIMPLE','BITPIX','NAXIS',                                   
                           'NAXIS1','EXTEND','BSCALE','BZERO','HISTORY']:
                hdu_out.header.set(key,
                                   self.header[key],
                                   self.header.comments[key])     
        for hist_message in self.header['HISTORY']:
            hdu_out.header.set('HISTORY', hist_message)

        
        # here, the gain is the electron per ADU.
        hdu_out.header.set('GAIN_A', self.gains['A'], 'e/ADU for CCD Sector A')
        hdu_out.header.set('GAIN_B', self.gains['B'], 'e/ADU for CCD Sector B')
        hdu_out.header.set('GAIN_C', self.gains['C'], 'e/ADU for CCD Sector C')
        hdu_out.header.set('GAIN_D', self.gains['D'], 'e/ADU for CCD Sector D')
        hdu_out.header.set('UNITS',  'electrons', 'Units after calibration (ADU or electrons)')

        # indexed at one, inclusive
        hdu_out.header.set('SCIPIXS','[45:2092,1:2048]','')
        hdu_out.header.set('TICAVER',version, 'tica software version')
        hdu_out.header.set('COMMENT', 'calibration applied at {time}'.format(
            time=datetime.datetime.utcnow().isoformat())
        )

            
        hdulist = fits.HDUList([hdu_out])
        
        inpath, infilename = os.path.split(self.fname)   # get path and original filename
        inbasename, inext = os.path.splitext(infilename)  # get basename and extension 
        if inext == '.gz':
            inbasename, inext = os.path.splitext(inbasename)

        if no_gzip:
            write_gz = False
        else:
            write_gz = True

        if write_gz:
            hdupath = os.path.join(self.outdir, inbasename + '.cal' + inext + '.gz')  #generate new path
        else:
            hdupath = os.path.join(self.outdir, inbasename + '.cal' + inext)  # generate new path

            
        try:
            self.logger.info("generating output calibration file: {}".format( hdupath ) )
            hdulist.writeto(hdupath)
        except IOError:
            self.logger.info("output calibrated file already exists, skipping generation: {}".format(hdupath) )

            
        
class FFI(object):
    """
    Handles reading, writing, parsing headers, in a way unique to TESS
    frames (maybe only applicable to bias?)

    'hdu' refers to astropy's 'Header Data Units'

    self.CCDs:  list of ccds, counterclockwise order from top left

    """
    def __init__(self, full_image, calibration=None):
        self.full_image  = full_image
        # defined so ccd1 is top right, then goes counterclockwise
        # at the end, will have 4 CCDs indexed in the same convention
        ccd1 = CCD(self.full_image[2078:, 2136:], calibration=calibration,inverted=True)
        ccd1.ccdnum = 1
        ccd1.gains = {'A':1,'B':1,'C':1,'D':1}                                                                             
        ccd1.gplus  = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
        ccd1.gminus = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}


        ccd2 = CCD(self.full_image[2078:, :2136], calibration=calibration,inverted=True)
        ccd2.ccdnum = 2
        ccd2.gains = {'A':1,'B':1,'C':1,'D':1}
        ccd2.gplus  = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
        ccd2.gminus = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
             
        ccd3 = CCD(self.full_image[:2078, :2136], calibration=calibration)
        ccd3.ccdnum = 3
        ccd3.gains = {'A':1,'B':1,'C':1,'D':1}                                                                         
        ccd3.gplus  = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
        ccd3.gminus = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}

 
        ccd4 = CCD(self.full_image[:2078, 2136:], calibration=calibration)
        ccd4.ccdnum = 4
        ccd4.gains = {'A':1,'B':1,'C':1,'D':1}
        ccd4.gplus  = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
        ccd4.gminus = {'A':0.0,'B':0.0,'C':0.0,'D':0.0}
 

        self.camnum = None
        self.CCDs = np.asarray([ ccd1, ccd2, ccd3, ccd4 ])
        self.logger = logging.getLogger(__name__)


    def get_frame(self, CCDlist):
        ##  A = Sector(self.arr[:, 0:11],   self.arr[:2048, 44:556],   self.arr[:, 2092:2103])
        ##  B = Sector(self.arr[:, 11:22], self.arr[:2048, 556:1068],  self.arr[:, 2103:2114])
        ##  C = Sector(self.arr[:, 22:33], self.arr[:2048, 1068:1580], self.arr[:, 2114:2125])
        ##  D = Sector(self.arr[:, 33:44], self.arr[:2048, 1580:2092], self.arr[:, 2125:2136])
        ##  self.sectors = np.asarray([A, B, C, D])
        out = []
        for ccd in CCDlist:
            A = ccd.sectors[0]
            B = ccd.sectors[1]
            C = ccd.sectors[2]
            D = ccd.sectors[3]
            out.append( np.c_[  A.underclock, B.underclock, C.underclock, D.underclock,
                                A.science, B.science, C.science, D.science,
                                A.overclock,   B.overclock,    C.overclock,   D.overclock   ])
            
        return np.r_[
            np.c_[  out[2], out[3] ],
            np.c_[  out[1][::-1, ::-1], out[0][::-1, ::-1]  ]
        ]


    def calibrate_CCDs(self):
        calibrated_CCDs = []
        for ccd in self.CCDs:
            calibrated_CCDs.append(
                ccd.calibrate()
            )

        return calibrated_CCDs

    def write_frame(self, 
                    CCDlist, 
                    stem,
                    no_gzip=False):

        hdu_out = fits.PrimaryHDU(self.get_frame(CCDlist).astype(np.float32))
        for key in self.header.keys():                                                     
            #censor the standard fits headers, which have to be changed.                      
            #Let astropy handle internally,  Put in cam latter to be near CCD                                                    
            if key not in ['SIMPLE','BITPIX','NAXIS',                                         
                           'NAXIS1','EXTEND','BSCALE',
                           'BZERO','CAM']:
                hdu_out.header.set(key, self.header[key], self.header.comments[key])                                              


        hdu_out.header.set('CAMNUM', self.header['CAM'], 'Camera Number')
        # indexed at one, inclusive
        hdu_out.header.set('SCIPIXS','[45:2092,1:2048]','')
        hdu_out.header.set('TICAVER',version, 'tica software version')

        hdu_out.header.set('UNITS',  'electrons', 'Units (ADU or electrons)')
        hdu_out.header.set('COMMENT', 'calibration applied at {time}'.format(
            time=datetime.datetime.utcnow().isoformat())
        )

        hdulist = fits.HDUList([hdu_out])
        
        inpath, infilename = os.path.split(self.fname)   # get path and original filename
        inbasename, inext = os.path.splitext(infilename)  # get basename and extension 

        if inext == '.gz':
            inbasename, inext = os.path.splitext(inbasename)

        if no_gzip:
            write_gz = False
        else:
            write_gz = True


        if write_gz:
            hdupath = os.path.join(self.outdir, inbasename + stem + inext + '.gz')  #generate new path
        else:
            hdupath = os.path.join(self.outdir, inbasename + stem + inext)  # generate new path
                
        try:
            hdulist.writeto(hdupath)
            self.logger.info("generating output calibration file: {}".foramt( hdupath) )
        except IOError:
            self.logger.info("output calibrated file already exists, skipping generation: {}".foramt(  hdupath) )



        
        
class FFI_File(FFI):
    """
    Handles reading, writing, parsing headers, in a way unique to TESS
    frames (maybe only applicable to bias?)

    'hdu' refers to astropy's 'Header Data Units'

    self.CCDs:  list of ccds, counterclockwise order from top left

    """
    def __init__(self, fname, calibration, outdir, orbit_segment = None):
        self.fname = fname
        self.outdir = outdir
    
        self.orbit_segment = orbit_segment

        self.hdu = fits.open(fname, mode='readonly')
        self.header = self.hdu[0].header
        im  = self.hdu[0].data

    

        try:
            #SPOC keyword
            assert np.isclose(calibration.int_time,
                              self.header['EXPOSURE']*86400/0.8/0.99), \
                "It appears you are trying to calibrate {} second data with models for {} second data".format(
                self.header['EXPOSURE']*86400/0.8/0.99,
                calibration.int_time)
        except KeyError:
            #TSO keyword
            assert calibration.int_time == self.header['INT_TIME'], \
                "It appears you are trying to calibrate {} second data with models for {} second data".format(
                    self.header['INT_TIME'],
                    calibration.int_time)

        super(FFI_File, self).__init__(im, calibration=calibration)
        self.camnum = self.header['CAM']
        
        for ccd in self.CCDs:
            ccd.camnum = self.camnum
            ccd.gains = GainModel.gains['cam' + str(ccd.camnum)]['ccd' + str(ccd.ccdnum)]
            ccd.gplus  = LinearityModel.gain_gain['cam' + str(ccd.camnum)]['ccd' + str(ccd.ccdnum)]
            ccd.gminus = LinearityModel.gain_loss['cam' + str(ccd.camnum)]['ccd' + str(ccd.ccdnum)]

            ccd.coadds = self.header['INT_TIME']*0.8/2




        try:
            self.calibrated_CCDs = self.calibrate_CCDs()
        except:
            print('calibration failed')
            raise


    def write_CCD_files(self,
                        CCDlist, 
                        stem, 
                        calibrated=True,
                        no_gzip=False):

        for i,ccd in enumerate(CCDlist):      
            hdu_out = fits.PrimaryHDU(ccd.get_frame().astype(np.float32))
            for key in self.header.keys():                                                     
                #censor the standard fits headers, which have to be changed.                      
                #Let astropy handle internally,  Put in cam latter to be near CCD                                                    
                if key not in ['SIMPLE','BITPIX','NAXIS',                                         
                               'NAXIS1','EXTEND','BSCALE',
                               'BZERO','CAM','COMMENT']:
                    hdu_out.header.set(key, self.header[key], self.header.comments[key])                                              
            if self.orbit_segment is not None:
                hdu_out.header.set('ORB_SEG', self.orbit_segment,
                                   'Orbit Identifier, o1a, o1b, o2a, o2b'
                                   )


            hdu_out.header.set('CAMNUM', self.header['CAM'], 'Camera Number')
            hdu_out.header.set('CCDNUM',i+1,'CCD number')
            # indexed at one, inclusive
            hdu_out.header.set('SCIPIXS','[45:2092,1:2048]','')


            if calibrated:
                # in CCD_file, this was self.gains['A'], etc., which gives the correct value
                hdu_out.header.set('GAIN_A', ccd.gains['A'], 'e/ADU for CCD Sector A')
                hdu_out.header.set('GAIN_B', ccd.gains['B'], 'e/ADU for CCD Sector B')
                hdu_out.header.set('GAIN_C', ccd.gains['C'], 'e/ADU for CCD Sector C')
                hdu_out.header.set('GAIN_D', ccd.gains['D'], 'e/ADU for CCD Sector D')
                hdu_out.header.set('UNITS',  'electrons', 'Units (ADU or electrons)')
                hdu_out.header.set('TICAVER',version, 'tica software version')
                hdu_out.header.set('COMMENT', 'calibration applied at {time}'.format(
                    time=datetime.datetime.utcnow().isoformat())
                )


            else:
                hdu_out.header.set('UNITS',  'ADU', 'Units (ADU or electrons)')
                

            hdulist = fits.HDUList([hdu_out])
        
            inpath, infilename = os.path.split(self.fname)   # get path and original filename
            inbasename, inext = os.path.splitext(infilename)  # get basename and extension 
            if inext == '.gz' :
                inbasename, inext = os.path.splitext(inbasename)

            if no_gzip:
                write_gz = False
            else:
                write_gz = True


            if write_gz:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext + '.gz')  # generate new path
            else:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext)  # generate new path

            try:
                ccd.logger.info("generating output calibration file: {}".format(  hdupath) )
                hdulist.writeto(hdupath)
            except IOError:
                ccd.logger.info("output calibrated file already exists, skipping generation: {}".format( hdupath) )

        


    def write_trimmed_CCD_files(self,
                                CCDlist, 
                                stem, 
                                calibrated=True,
                                no_gzip=False):
        for i,ccd in enumerate(CCDlist):      
#            hdu_out = fits.PrimaryHDU(ccd.get_frame().astype(np.float32))
            hdu_out = fits.PrimaryHDU(ccd.get_image().astype(np.float32)[0:2048,0:2048])
            for key in self.header.keys():                                                     
                #censor the standard fits headers, which have to be changed.                      
                #Let astropy handle internally,  Put in cam latter to be near CCD                                                    
                if key not in ['SIMPLE','BITPIX','NAXIS',                                         
                               'NAXIS1','EXTEND','BSCALE',
                               'BZERO','CAM']:                            
                    hdu_out.header.set(key, self.header[key], self.header.comments[key])                                              
                                                       
            hdu_out.header.set('CAMNUM', self.header['CAM'], 'Camera Number')
            hdu_out.header.set('CCDNUM',i+1,'CCD number')
            # indexed at one, inclusive
            hdu_out.header.set('SCIPIXS','[1:2048,1:2048]','')

            if calibrated:
                # the gain between each sector has been normalized to the mean gain
                hdu_out.header.set('GAIN_A', ccd.gains['A'], 'e/ADU for CCD Sector A')
                hdu_out.header.set('GAIN_B', ccd.gains['B'], 'e/ADU for CCD Sector B')
                hdu_out.header.set('GAIN_C', ccd.gains['C'], 'e/ADU for CCD Sector C')
                hdu_out.header.set('GAIN_D', ccd.gains['D'], 'e/ADU for CCD Sector D')
                hdu_out.header.set('UNITS',  'electrons', 'Units  (ADU or electrons)')
                hdu_out.header.set('TICAVER',version, 'tica software version')
                hdu_out.header.set('COMMENT', 'calibration applied at {time}'.format(
                    time=datetime.datetime.utcnow().isoformat())
                )

            else:
                hdu_out.header.set('UNITS',  'ADU', 'Units (ADU or electrons)')



            hdulist = fits.HDUList([hdu_out])
        
            inpath, infilename = os.path.split(self.fname)   # get path and original filename
            inbasename, inext = os.path.splitext(infilename)  # get basename and extension 
            hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext)  # generate new path

            inpath, infilename = os.path.split(self.fname)   # get path and original filename
            inbasename, inext = os.path.splitext(infilename)  # get basename and extension
            if inext == '.gz':
                inbasename, inext = os.path.splitext(inbasename)

            if no_gzip:
                write_gz = False
            else:
                write_gz = True


            if write_gz:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext + '.gz')  # generate new path
            else:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext)  # generate new path

            try:
                hdulist.writeto(hdupath)
                ccd.logger.info("generating output calibration file: {}".format( hdupath) )
            except IOError:
                ccd.logger.info("output calibrated file already exists, skipping generation: {}".foramt(  hdupath) )



    #These handle the scripts in tica/bin
    def write_raw_CCDs(self, no_gzip = False):
        self.write_CCD_files(self.CCDs,'.raw',calibrated=False, no_gzip=no_gzip)

    def write_calibrated_CCDs(self, no_gzip = False):
        self.write_CCD_files(self.calibrated_CCDs,'.cal', no_gzip = no_gzip)

    def write_calibrated_trimmed_CCDs(self, no_gzip = False):
        self.write_trimmed_CCD_files(self.calibrated_CCDs,'.cal',
                                     no_gzip = False)
    

