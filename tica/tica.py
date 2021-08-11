# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Some code lifted from Miranda Kephart, ccd_stats.py on tessmp:/usr/local/bin

Most by MMF, Aug/Sept 2017

order is additive then multiplicative, for the expected data (I think)
"""

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

fstem = os.path.abspath(os.path.dirname(__file__) + '/../')
with open(os.path.join(fstem, 'VERSION'),'r') as infile:
    version = infile.read()


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
    Keep track of individual CCD channels, i.e., sectors
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
    """Class to handle initialization of calibrations so it is performed just once
    """

    @staticmethod
    def _set_calibration_data(cal_file, cal_dir, logger):
        """given calibration files and directory return the FITs data structure
        or else return none if no file exists to allow for placeholder"""
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
        
        # FIXME: do we need to handle all at least once in a run? or should it be semi-on-demand?
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
    """Splits a CCD into its sectors. If the CCD is passed in inverted (this
    is the case for CCDs 1 and 2), un-inverts it (i.e. rotates it 180 degrees)
    before proceeding).

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
            
        #should probably have these attributes, since all calibration is tied to them...
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
        #JPD version
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
        # FIXME: decide on subsection of virtual columns,

        # I see in the report that some CCDs have a weak trend along
        # rows.  For now, we will just take out a constant (quick look
        # pipeline, afterall), with pains to avoid start of row ringing and associated curvature.

        # each CCD has a different patter for the trend along columns.
        # It is known that the first few columns are suspect.
        # Something in the middle seems "safe", and as long as that
        # comes out ot zero after 1D correection.  some confusion
        # about what happens at the last few columns

        mask_col = slice(3,  None)
        mask_row = slice(750,None)

        cal = CCD(np.zeros((2078,2136)))
        for i in range(len(self.sectors)):
            correct = np.mean(self.sectors[i].overclock[mask_row, mask_col] )
            cal.sectors[i].underclock = self.sectors[i].underclock - correct
            cal.sectors[i].science      = self.sectors[i].science - correct
            cal.sectors[i].overclock   = self.sectors[i].overclock - correct
        return cal
    
    def convert_to_electrons(self, out):
        """default is gain == 1, less it is CCD file with cam/ccd info in the
        header
        
        Note that the self.gains is applied on input data.  So the
        class knows about the correct value, but never applies it


        For now, this does the calculation on the original data
        struct!  Orginal idea was a functional kind of thing, so that
        you can't be bit by state changes....  consider changing

        """
        guse = np.asarray([self.gains['A'], self.gains['B'], self.gains['C'], self.gains['D'] ])
        for i,s in enumerate(out.sectors):
            s.underclock *= guse[i]
            s.science *= guse[i]
            s.overclock *= guse[i]
        return out

    def linearity_correct(self, out):
        """Uses John Doty's model, which is a kind of sigmoid.
       
        See IH, 6.4.2:

        e_out = e_in/((1 + e_in*gplus)*(1 - e_in*gminus))

        where gplus = gainGain_per_electron and gminus =
        gainLoss_per_electron are measured parameters, stored in the
        LinearityModel Class.

        e_in/e_out are in units of electrons/2 second exposure

        Do not apply to under/over clocks

        For now, this does the calculation on the original data
        struct!  Orginal idea was a functional kind of thing, so that
        you can't be bit by state changes....  consider changing

        """
        gplus  = np.asarray([self.gplus['A'],  self.gplus['B'],  self.gplus['C'],  self.gplus['D'] ])
        gminus = np.asarray([self.gminus['A'], self.gminus['B'], self.gminus['C'], self.gminus['D']])
        for i,s in enumerate(out.sectors):
#            print self.camnum, self.ccdnum, i,  gplus[i],gminus[i], np.mean(( (1 + s.science/720.0*gplus[i])*(1 - s.science/720.0*gminus[i]) ))
            s.science = s.science/( (1 + s.science/720.0/0.99*gplus[i])*(1 - s.science/720.0/0.99*gminus[i]) )
        return out

    def spoc_linearity_correct(self, out):
        """Uses spoc parameterization, which is stored in a class.  Assign the
        class in this object's __init__

        Do not apply to under/over clocks

        For now, this does the calculation on the original data
        struct!  Orginal idea was a functional kind of thing, so that
        you can't be bit by state changes....  consider changing

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
        if len(self.hdu) == 1:
            self.header = self.hdu[0].header
            image  = self.hdu[0].data
            
        #handles either spoc raw or cal---maybe more checks are needed??
        elif len(self.hdu) == 2 or len(self.hdu) == 3:
            assert self.hdu[0].data == None

            self.header = fits.Header()

            for key in self.hdu[0].header.keys():
                if key not in ['SIMPLE','BITPIX','NAXIS',
                               'EXTEND','NEXTEND',
                               'EXTNAME','EXTVER',
                               'CHECKSUM']:
                    self.header.set(key,
                                    self.hdu[0].header[key],
                                    self.hdu[0].header.comments[key])
                    
            for key in self.hdu[1].header.keys():
                if key not in ['XTENSION','BITPIX','NAXIS',
                               'NAXIS1','EXTEND','BSCALE','BZERO',
                               'INHERIT','EXTNAME','EXTVER']:
                    self.header.set(key,
                                    self.hdu[1].header[key],
                                    self.hdu[1].header.comments[key])
                    
            image  = self.hdu[1].data

        else:
            raise IOError("must be a TSO or SPOC FFI")


        super(CCD_File, self).__init__(image, calibration=calibration)

        #error here on the defaul for SPOC, since 'CAM' is 'CAMERA' in their files
        self.ccdnum  = self.header.get('CCD')
        self.camnum = self.header.get('CAM')
        if self.camnum is None:
            self.camnum = self.header.get('CAMERA')
        if self.camnum is None:
            raise ValueError('Did not find camera number in the FITS header')

        self.gains = GainModel.gains['cam' + str(self.camnum)]['ccd' + str(self.ccdnum)]
        #for linearity
        self.gplus  = LinearityModel.gain_gain['cam' + str(self.camnum)]['ccd' + str(self.ccdnum)]
        self.gminus = LinearityModel.gain_loss['cam' + str(self.camnum)]['ccd' + str(self.ccdnum)]

        try:
            self.coadds = self.header['NREADOUT']
        except:
            self.coadds = self.header['INT_TIME']*0.8/2

        try:
            self.calibrated_frame = self.calibrate()
        except:
            print('calibration failed')

    def write_calibrate(self):      
        hdu_out = fits.PrimaryHDU(self.calibrated_frame.get_frame().astype(np.float32))
        for key in self.header.keys():                                                     
            #censor the standard fits headers, which have to be changed.                      
            #Let astropy handle internally                                                    
            if key not in ['SIMPLE','BITPIX','NAXIS',                                         
                           'NAXIS1','EXTEND','BSCALE','BZERO']:
                hdu_out.header.set(key,
                                   self.header[key],
                                   self.header.comments[key])     
                                               

        
        # here, the gain is the electron per ADU.
        #in the FFI version, it looks this up in the output calibrated CCD, which has this property stored as 1.0---see line 311
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
            write_gz = True
            inbasename, inext = os.path.splitext(inbasename)
        else:
            write_gz = False

        if write_gz:
            hdupath = os.path.join(self.outdir, inbasename + '.cal' + inext + '.gz')  #generate new path
        else:
            hdupath = os.path.join(self.outdir, inbasename + '.cal' + inext)  # generate new path

            
        # FIXME: should we have a forced override option?
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

    methods:
    restore (to fix if modified data are wronge, e.g., after applying calibration, etc.)
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

    def write_frame(self, CCDlist, stem):
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
            write_gz = True
            inbasename, inext = os.path.splitext(inbasename)
        else:
            write_gz = False

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

    methods:
    restore (to fix if modified data are wronge, e.g., after applying calibration, etc.)
    """
    def __init__(self, fname, calibration, outdir):
        self.fname = fname
        self.outdir = outdir
        
        self.hdu = fits.open(fname, mode='readonly')
        self.header = self.hdu[0].header
        im  = self.hdu[0].data

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


    def write_CCD_files(self,CCDlist, stem, calibrated=True):
        for i,ccd in enumerate(CCDlist):      
            hdu_out = fits.PrimaryHDU(ccd.get_frame().astype(np.float32))
            for key in self.header.keys():                                                     
                #censor the standard fits headers, which have to be changed.                      
                #Let astropy handle internally,  Put in cam latter to be near CCD                                                    
                if key not in ['SIMPLE','BITPIX','NAXIS',                                         
                               'NAXIS1','EXTEND','BSCALE',
                               'BZERO','CAM','COMMENT']:
                    hdu_out.header.set(key, self.header[key], self.header.comments[key])                                              


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
            if inext == '.gz':
                write_gz = True
                inbasename, inext = os.path.splitext(inbasename)
            else:
                write_gz = False

            if write_gz:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext + '.gz')  # generate new path
            else:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext)  # generate new path

            try:
                ccd.logger.info("generating output calibration file: {}".format(  hdupath) )
                hdulist.writeto(hdupath)
            except IOError:
                ccd.logger.info("output calibrated file already exists, skipping generation: {}".format( hdupath) )

        


    def write_trimmed_CCD_files(self,CCDlist, stem, calibrated=True):
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
                write_gz = True
                inbasename, inext = os.path.splitext(inbasename)
            else:
                write_gz = False

            if write_gz:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext + '.gz')  # generate new path
            else:
                hdupath = os.path.join(self.outdir, inbasename + '_ccd'+str(i+1) + stem + inext)  # generate new path

            try:
                hdulist.writeto(hdupath)
                ccd.logger.info("generating output calibration file: {}".format( hdupath) )
            except IOError:
                ccd.logger.info("output calibrated file already exists, skipping generation: {}".foramt(  hdupath) )



    #These hand the API of the scripts in tica/bin
    def write_raw_CCDs(self):
        self.write_CCD_files(self.CCDs,'.raw',calibrated=False)

    def write_calibrated_CCDs(self):
        self.write_CCD_files(self.calibrated_CCDs,'.cal')

    def write_calibrated_trimmed_CCDs(self):
        self.write_trimmed_CCD_files(self.calibrated_CCDs,'.cal')
    



class GainModel(object):
    """
    class properties--hard code the gain per cam/ccd/slice
    """
    gains ={'cam1':{'ccd1':{"A":5.22,
                            "B":5.21,
                            "C":5.21,
                            "D":5.26},
                    'ccd2':{"A":5.27,
                            "B":5.14,
                            "C":5.11,
                            "D":5.19},
                    'ccd3':{"A":5.32,
                            "B":5.23,
                            "C":5.2 ,
                            "D":5.24},
                    'ccd4':{"A":5.35,
                            "B":5.22, 
                            "C":5.22, 
                            "D":5.18} 
                },
            'cam2':{'ccd1':{"A":5.31,
                            "B":5.24,
                            "C":5.23,
                            "D":5.24},
                    'ccd2':{"A":5.22,
                            "B":5.28,
                            "C":5.32,
                            "D":5.22},
                    'ccd3':{"A":5.28,
                            "B":5.26,
                            "C":5.25,
                            "D":5.21},
                    'ccd4':{"A":5.33,
                            "B":5.2 ,
                            "C":5.3 ,
                            "D":5.23}
                } ,
            'cam3':{'ccd1':{"A":5.25,
                            "B":5.24,
                            "C":5.24,
                            "D":5.24},
                    'ccd2':{"A":5.36,
                            "B":5.26,
                            "C":5.36,
                            "D":5.29},
                    'ccd3':{"A":5.26,
                            "B":5.16,
                            "C":5.22,
                            "D":5.24},
                    'ccd4':{"A":5.17,
                            "B":5.15,
                            "C":5.15,
                            "D":5.17}
                } ,
            'cam4':{'ccd1':{"A":5.26,
                            "B":5.18,
                            "C":5.2,
                            "D":5.2},
                    'ccd2':{"A":5.25,
                            "B":5.17,
                            "C":5.12,
                            "D":5.2},
                    'ccd3':{"A":5.3,
                            "B":5.18,
                            "C":5.27,
                            "D":5.19,},
                    'ccd4':{"A":5.24,
                            "B":5.12,
                            "C":5.16,
                            "D":5.16}
                } 
        }
    def __init__(self,modelfile):
        pass
        
class LinearityModel(object):
    """
    class properties--hard code the gain per cam/ccd/slice
    """
    gain_loss ={'cam1':{'ccd1':{"A":5.91e-07,
                                "B":6.17e-07,
                                "C":5.72e-07,
                                "D":5.85e-07},
                        'ccd2':{"A":6.63e-07,
                                "B":6.36e-07,
                                "C":6.07e-07,
                                "D":6.19e-07},
                        'ccd3':{"A":6.41e-07,
                                "B":6.12e-07,
                                "C":6.11e-07,
                                "D":6.36e-07},
                        'ccd4':{"A":5.49e-07,
                                "B":5.82e-07, 
                                "C":5.90e-07, 
                                "D":5.97e-07} 
                    },              
                'cam2':{'ccd1':{"A":6.21e-07,
                                "B":5.82e-07,
                                "C":5.57e-07,
                                "D":6.21e-07},
                        'ccd2':{"A":6.41e-07,
                                "B":6.30e-07,
                                "C":6.34e-07,
                                "D":6.44e-07},
                        'ccd3':{"A":5.95e-07,
                                "B":5.85e-07,
                                "C":6.02e-07,
                                "D":6.15e-07},
                        'ccd4':{"A":5.72e-07,
                                "B":5.85e-07,
                                "C":5.96e-07,
                                "D":6.05e-07}
                    } ,             
                'cam3':{'ccd1':{"A":5.16e-07,
                                "B":6.17e-07,
                                "C":5.28e-07,
                                "D":5.64e-07},
                        'ccd2':{"A":5.57e-07,
                                "B":5.18e-07,
                                "C":5.98e-07,
                                "D":6.42e-07},
                        'ccd3':{"A":5.89e-07,
                                "B":5.86e-07,
                                "C":5.72e-07,
                                "D":5.96e-07},
                        'ccd4':{"A":6.03e-07,
                                "B":5.09e-07,
                                "C":5.65e-07,
                                "D":6.36e-07}
                    } ,             
                'cam4':{'ccd1':{"A":6.29e-07,
                                "B":5.87e-07,
                                "C":6.31e-07,
                                "D":6.36e-07},
                        'ccd2':{"A":5.97e-07,
                                "B":5.20e-07,
                                "C":5.86e-07,
                                "D":5.17e-07},
                        'ccd3':{"A":5.71e-07,
                                "B":6.10e-07,
                                "C":5.94e-07,
                                "D":6.21e-07},
                        'ccd4':{"A":5.78e-07,
                                "B":6.02e-07,
                                "C":5.97e-07,
                                "D":6.28e-07}
                     }
                }

    gain_gain ={'cam1':{'ccd1':{"A":6.96e-07,
                                "B":7.22e-07,
                                "C":6.65e-07,
                                "D":6.99e-07},
                        'ccd2':{"A":7.69e-07,
                                "B":7.32e-07,
                                "C":7.06e-07,
                                "D":7.31e-07},
                        'ccd3':{"A":7.69e-07,
                                "B":7.22e-07,
                                "C":7.35e-07,
                                "D":7.75e-07},
                        'ccd4':{"A":6.62e-07,
                                "B":7.06e-07, 
                                "C":7.23e-07, 
                                "D":7.25e-07} 
                    },              
                'cam2':{'ccd1':{"A":7.74e-07,
                                "B":7.01e-07,
                                "C":6.80e-07,
                                "D":7.60e-07},
                        'ccd2':{"A":7.27e-07,
                                "B":7.32e-07,
                                "C":7.47e-07,
                                "D":7.73e-07},
                        'ccd3':{"A":7.27e-07,
                                "B":7.06e-07,
                                "C":7.35e-07,
                                "D":7.51e-07},
                        'ccd4':{"A":7.14e-07,
                                "B":7.19e-07,
                                "C":7.44e-07,
                                "D":7.58e-07}
                    } ,             
                'cam3':{'ccd1':{"A":6.24e-07,
                                "B":7.56e-07,
                                "C":6.46e-07,
                                "D":6.85e-07},
                        'ccd2':{"A":6.51e-07,
                                "B":5.81e-07,
                                "C":7.08e-07,
                                "D":7.86e-07},
                        'ccd3':{"A":6.94e-07,
                                "B":6.79e-07,
                                "C":6.68e-07,
                                "D":7.10e-07},
                        'ccd4':{"A":7.16e-07,
                                "B":5.87e-07,
                                "C":6.81e-07,
                                "D":7.84e-07}
                    } ,             
                'cam4':{'ccd1':{"A":7.71e-07,
                                "B":6.86e-07,
                                "C":7.75e-07,
                                "D":7.66e-07},
                        'ccd2':{"A":7.24e-07,
                                "B":6.02e-07,
                                "C":7.04e-07,
                                "D":6.09e-07},
                        'ccd3':{"A":7.11e-07,
                                "B":7.58e-07,
                                "C":7.38e-07,
                                "D":7.78e-07,},
                        'ccd4':{"A":6.90e-07,
                                "B":7.27e-07,
                                "C":7.16e-07,
                                "D":7.75e-07}
                    }
            }
    def __init__(self,modelfile):
        pass
        
class SpocLinearity(object):
    """
    For spoc, we need the coefficients of the 2nd order polynomial,
    the scale factor for the input, and the offset of the input
    """
    coefficients ={'cam1':{'ccd1':{"A":[2.57322e-05, 823.226, 1.0,-0.0374716,0.0482595    ],
                                  "B":[2.61974e-05, 911.529, 1.0,-0.0321036,0.0472505    ],
                                  "C":[2.57796e-05, 861.398, 1.0,-0.0299381,0.0410301    ],
                                  "D":[2.60366e-05, 891.072, 1.0,-0.0359014,0.0415674    ]},
                          'ccd2':{"A":[2.65022e-05, 1003.89, 1.0,-0.0413358,0.0683761    ],
                                  "B":[2.61363e-05, 1019.84, 1.0,-0.031161,0.0534969     ],
                                  "C":[2.60988e-05, 893.095, 1.0,-0.0271916,0.0394685    ],
                                  "D":[2.60189e-05, 814.683, 1.0,-0.0377196,0.0490083    ]},
                          'ccd3':{"A":[2.59058e-05, 919.8, 1.0,-0.0393967,0.0525883      ],
                                  "B":[2.57147e-05, 976.537, 1.0,-0.0374545,0.0528702    ],
                                  "C":[2.57727e-05, 900.984, 1.0,-0.0389461,0.045786     ],
                                  "D":[2.59371e-05, 892.831, 1.0,-0.0454134,0.0524089    ]},
                          'ccd4':{"A":[2.6866e-05, 981.681, 1.0,-0.0507134,0.0562926     ],
                                  "B":[2.58883e-05, 965.747, 1.0,-0.0446359,0.0489371    ], 
                                  "C":[2.60191e-05, 775.707, 1.0,-0.0492108,0.0501099    ], 
                                  "D":[2.61775e-05, 683.515, 1.0,-0.0486886,0.0519709    ]} 
                    },                
                  'cam2':{'ccd1':{"A":[2.58758e-05, 914.283, 1.0,-0.0581077,0.0606467    ],
                                  "B":[2.59669e-05, 988.901, 1.0,-0.0439921,0.0497028    ],
                                  "C":[2.61721e-05, 952.948, 1.0,-0.0447996,0.043987     ],
                                  "D":[2.57771e-05, 969.021, 1.0,-0.0502004,0.0565481    ]},
                          'ccd2':{"A":[2.70836e-05, 1040.49, 1.0,-0.0332812,0.0576825    ],
                                  "B":[2.61238e-05, 1017.51, 1.0,-0.0369586,0.0572305    ],
                                  "C":[2.57524e-05, 904.56, 1.0,-0.0421833,0.0603724     ],
                                  "D":[2.57664e-05, 840.23, 1.0,-0.0475555,0.0606503     ]},
                          'ccd3':{"A":[2.58508e-05, 916.877, 1.0,-0.0502073,0.0547111    ],
                                  "B":[2.56261e-05, 975.949, 1.0,-0.0452499,0.05191      ],
                                  "C":[2.5947e-05, 892.329, 1.0,-0.0495091,0.0542741     ],
                                  "D":[2.58622e-05, 917.551, 1.0,-0.0508815,0.0561949    ]},
                          'ccd4':{"A":[2.60322e-05, 1008.83, 1.0,-0.0542148,0.0518747    ],
                                  "B":[2.59428e-05, 1013.74, 1.0,-0.0494182,0.0509864    ],
                                  "C":[2.57109e-05, 836.186, 1.0,-0.0563779,0.0553216    ],
                                  "D":[2.62769e-05, 765.641, 1.0,-0.0565943,0.0538907    ]}
                    } ,               
                  'cam3':{'ccd1':{"A":[2.65518e-05, 871.101, 1.0,-0.0485704,0.0507719    ],
                                  "B":[2.72475e-05, 940.163, 1.0,-0.0480001,0.0522887    ],
                                  "C":[2.57254e-05, 897.272, 1.0,-0.0524578,0.0547793    ],
                                  "D":[2.79272e-05, 916.949, 1.0,-0.0485449,0.0506356    ]},
                          'ccd2':{"A":[2.61573e-05, 983.104, 1.0,-0.0478297,0.0641121    ],
                                  "B":[2.62336e-05, 985.372, 1.0,-0.034413,0.0537814     ],
                                  "C":[2.6387e-05, 840.47, 1.0,-0.0448661,0.0600269      ],
                                  "D":[2.6297e-05, 768.779, 1.0,-0.0530794,0.0639607     ]},
                          'ccd3':{"A":[2.65843e-05, 874.064, 1.0,-0.0464842,0.0612928    ],
                                  "B":[2.69553e-05, 956.817, 1.0,-0.0368783,0.0512449    ],
                                  "C":[2.58179e-05, 858.373, 1.0,-0.038546,0.0531443     ],
                                  "D":[2.65562e-05, 874.34, 1.0,-0.0504402,0.0637462     ]},
                          'ccd4':{"A":[2.8381e-05, 1014.94, 1.0,-0.0397823,0.0471796     ],
                                  "B":[2.59052e-05, 994.259, 1.0,-0.03949,0.0520631      ],
                                  "C":[2.59421e-05, 813.502, 1.0,-0.0446166,0.0497614    ],
                                  "D":[2.578e-05, 722.269, 1.0,-0.0535503,0.0608501      ]}
                    } ,               
                  'cam4':{'ccd1':{"A":[2.57391e-05, 878.855, 1.0,-0.0496681,0.0562381    ],
                                  "B":[2.78096e-05, 952.233, 1.0,-0.0255654,0.0330627    ],
                                  "C":[2.56824e-05, 901.636, 1.0,-0.0501126,0.0565003    ],
                                  "D":[2.60466e-05, 913.999, 1.0,-0.0448427,0.0563647    ]},
                          'ccd2':{"A":[2.69813e-05, 977.398, 1.0,-0.0398781,0.0423855    ],
                                  "B":[2.73665e-05, 967.056, 1.0,-0.0313423,0.0379962    ],
                                  "C":[2.61286e-05, 847.83, 1.0,-0.0394111,0.04401       ],
                                  "D":[2.70169e-05, 749.26, 1.0,-0.0432061,0.049153      ]},
                          'ccd3':{"A":[2.63179e-05, 836.514, 1.0,-0.0545174,0.0522306    ],
                                  "B":[2.58206e-05, 918.289, 1.0,-0.052097,0.0530815     ],
                                  "C":[2.57241e-05, 819.441, 1.0,-0.0529934,0.053071     ],
                                  "D":[2.58999e-05, 850.556, 1.0,-0.0552855,0.0548873    ],},
                          'ccd4':{"A":[2.68663e-05, 984.241, 1.0,-0.0327376,0.0366624    ],
                                  "B":[2.6579e-05, 992.954, 1.0,-0.0354987,0.0387493     ],
                                  "C":[2.64267e-05, 800.128, 1.0,-0.0409345,0.0463524    ],
                                  "D":[2.58413e-05, 720.125, 1.0,-0.0516081,0.0553248    ]}
                      }
              }
    def __init__(self):
        pass

    def linearity_correct(self,e_in,cam,ccd, channel):
        """
        An FFI in units of electrons is passed in.

        Need to convert back to ADU, scale to 2 seconds, and then apply the polynomial

        """

        cuse = self.coefficients['cam'+str(cam)]['ccd'+str(ccd)][channel]
        guse = GainModel.gains['cam'+str(cam)]['ccd'+str(ccd)][channel]
        scale = cuse[0]
        offset = cuse[1]
        
        x = (e_in/guse - offset)*scale
        return e_in*(cuse[2] + cuse[3]*x/2. + cuse[4]*x**2/3.)


class SpocLinearity_Updated(object):
    """
    For spoc, we need the coefficients of the 2nd order polynomial,
    the scale factor for the input, and the offset of the input
    """
    coefficients ={'cam1':{'ccd1':{"A":[2.57322E-05, 8.23226E+02, 1.00042, -1.77135E-02, 1.60865E-02],
                                   "B":[2.61974E-05, 9.11529E+02, 1.00041, -1.49235E-02, 1.57502E-02],
                                   "C":[2.57796E-05, 8.61398E+02, 1.00035, -1.40579E-02, 1.36767E-02],
                                   "D":[2.60366E-05, 8.91072E+02, 1.00044, -1.69863E-02, 1.38558E-02]},
                           'ccd2':{"A":[2.65022E-05, 1.00389E+03, 1.00060, -1.88487E-02, 2.27920E-02],
                                   "B":[2.61363E-05, 1.01984E+03, 1.00045, -1.41545E-02, 1.78323E-02],
                                   "C":[2.60988E-05, 8.93095E+02, 1.00034, -1.26758E-02, 1.31562E-02],
                                   "D":[2.60189E-05, 8.14683E+02, 1.00042, -1.78210E-02, 1.63361E-02]},
                           'ccd3':{"A":[2.59058E-05, 9.19800E+02, 1.00050, -1.84453E-02, 1.75294E-02],
                                   "B":[2.57147E-05, 9.76537E+02, 1.00050, -1.73996E-02, 1.76234E-02],
                                   "C":[2.57727E-05, 9.00984E+02, 1.00048, -1.84099E-02, 1.52620E-02],
                                   "D":[2.59371E-05, 8.92831E+02, 1.00055, -2.14930E-02, 1.74696E-02]},
                           'ccd4':{"A":[2.68660E-05, 9.81681E+02, 1.00071, -2.38720E-02, 1.87642E-02],
                                   "B":[2.58883E-05, 9.65747E+02, 1.00059, -2.10944E-02, 1.63124E-02], 
                                   "C":[2.60191E-05, 7.75707E+02, 1.00052, -2.35940E-02, 1.67033E-02], 
                                   "D":[2.61775E-05, 6.83515E+02, 1.00045, -2.34144E-02, 1.73236E-02]} 
                       },                
                  'cam2':{'ccd1':{"A":[2.58758E-05, 9.14283E+02, 1.00072, -2.76191E-02, 2.02156E-02],
                                  "B":[2.59669E-05, 9.88901E+02, 1.00060, -2.07197E-02, 1.65676E-02],
                                  "C":[2.61721E-05, 9.52948E+02, 1.00059, -2.13027E-02, 1.46623E-02],
                                  "D":[2.57771E-05, 9.69021E+02, 1.00066, -2.36877E-02, 1.88494E-02]},
                          'ccd2':{"A":[2.70836E-05, 1.04049E+03, 1.00051, -1.50151E-02, 1.92275E-02],
                                  "B":[2.61238E-05, 1.01751E+03, 1.00053, -1.69580E-02, 1.90768E-02],
                                  "C":[2.57524E-05, 9.04560E+02, 1.00052, -1.96853E-02, 2.01241E-02],
                                  "D":[2.57664E-05, 8.40230E+02, 1.00054, -2.24647E-02, 2.02168E-02]},
                          'ccd3':{"A":[2.58508E-05, 9.16877E+02, 1.00063, -2.38069E-02, 1.82370E-02],
                                  "B":[2.56261E-05, 9.75949E+02, 1.00060, -2.13267E-02, 1.73033E-02],
                                  "C":[2.59470E-05, 8.92329E+02, 1.00060, -2.34979E-02, 1.80914E-02],
                                  "D":[2.58622E-05, 9.17551E+02, 1.00064, -2.41073E-02, 1.87316E-02]},
                          'ccd4':{"A":[2.60322E-05, 1.00883E+03, 1.00075, -2.57451E-02, 1.72916E-02],
                                  "B":[2.59428E-05, 1.01374E+03, 1.00069, -2.33682E-02, 1.69955E-02],
                                  "C":[2.57109E-05, 8.36186E+02, 1.00063, -2.69996E-02, 1.84405E-02],
                                  "D":[2.62769E-05, 7.65641E+02, 1.00059, -2.72129E-02, 1.79636E-02]}
                    } ,               
                  'cam3':{'ccd1':{"A":[2.65518E-05, 8.71101E+02, 1.00059, -2.31109E-02, 1.69240E-02],
                                  "B":[2.72475E-05, 9.40163E+02, 1.00065, -2.26606E-02, 1.74296E-02],
                                  "C":[2.57254E-05, 8.97272E+02, 1.00063, -2.49644E-02, 1.82598E-02],
                                  "D":[2.79272E-05, 9.16949E+02, 1.00065, -2.29758E-02, 1.68785E-02]},
                          'ccd2':{"A":[2.61573E-05, 9.83104E+02, 1.00066, -2.22662E-02, 2.13707E-02],
                                  "B":[2.62336E-05, 9.85372E+02, 1.00048, -1.58163E-02, 1.79271E-02],
                                  "C":[2.63870E-05, 8.40470E+02, 1.00053, -2.11018E-02, 2.00090E-02],
                                  "D":[2.62970E-05, 7.68779E+02, 1.00056, -2.52466E-02, 2.13202E-02]},
                          'ccd3':{"A":[2.65843E-05, 8.74064E+02, 1.00057, -2.18179E-02, 2.04309E-02],
                                  "B":[2.69553E-05, 9.56817E+02, 1.00051, -1.71175E-02, 1.70816E-02],
                                  "C":[2.58179E-05, 8.58373E+02, 1.00045, -1.80952E-02, 1.77148E-02],
                                  "D":[2.65562E-05, 8.74340E+02, 1.00062, -2.37400E-02, 2.12487E-02]},
                          'ccd4':{"A":[2.83810E-05, 1.01494E+03, 1.00061, -1.85321E-02, 1.57265E-02],
                                  "B":[2.59052E-05, 9.94259E+02, 1.00054, -1.84040E-02, 1.73544E-02],
                                  "C":[2.59421E-05, 8.13502E+02, 1.00049, -2.12581E-02, 1.65871E-02],
                                  "D":[2.57800E-05, 7.22269E+02, 1.00052, -2.56421E-02, 2.02834E-02]}
                    } ,               
                  'cam4':{'ccd1':{"A":[2.57391E-05, 8.78855E+02, 1.00059, -2.35619E-02, 1.87460E-02],
                                  "B":[2.78096E-05, 9.52233E+02, 1.00036, -1.19072E-02, 1.10209E-02],
                                  "C":[2.56824E-05, 9.01636E+02, 1.00061, -2.37480E-02, 1.88334E-02],
                                  "D":[2.60466E-05, 9.13999E+02, 1.00057, -2.10795E-02, 1.87882E-02]},
                          'ccd2':{"A":[2.69813E-05, 9.77398E+02, 1.00056, -1.88213E-02, 1.41285E-02],
                                  "B":[2.73665E-05, 9.67056E+02, 1.00044, -1.46656E-02, 1.26654E-02],
                                  "C":[2.61286E-05, 8.47830E+02, 1.00046, -1.87306E-02, 1.46700E-02],
                                  "D":[2.70169E-05, 7.49260E+02, 1.00046, -2.06081E-02, 1.63843E-02]},
                          'ccd3':{"A":[2.63179E-05, 8.36514E+02, 1.00063, -2.61088E-02, 1.74102E-02],
                                  "B":[2.58206E-05, 9.18289E+02, 1.00065, -2.47899E-02, 1.76938E-02],
                                  "C":[2.57241E-05, 8.19441E+02, 1.00058, -2.53780E-02, 1.76903E-02],
                                  "D":[2.58999E-05, 8.50556E+02, 1.00064, -2.64336E-02, 1.82958E-02],},
                          'ccd4':{"A":[2.68663E-05, 9.84241E+02, 1.00046, -1.53993E-02, 1.22208E-02],
                                  "B":[2.65790E-05, 9.92954E+02, 1.00050, -1.67267E-02, 1.29164E-02],
                                  "C":[2.64267E-05, 8.00128E+02, 1.00045, -1.94871E-02, 1.54508E-02],
                                  "D":[2.58413E-05, 7.20125E+02, 1.00050, -2.47745E-02, 1.84416E-02]}
                      }
              }
    def __init__(self):
        pass

    def linearity_correct(self,e_in,cam,ccd, channel, N_coadds):
        """
        An FFI in units of electrons is passed in.

        Need to convert back to ADU, scale to 2 seconds, and then apply the polynomial

        """
        cuse = self.coefficients['cam'+str(cam)]['ccd'+str(ccd)][channel]
        guse = GainModel.gains['cam'+str(cam)]['ccd'+str(ccd)][channel]
        scale = cuse[0]
        offset = cuse[1]

        x = (e_in/N_coadds/0.99/guse - offset)*scale
        return e_in*(cuse[2] + cuse[3]*x + cuse[4]*x**2)

