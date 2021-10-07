import numpy as np
from astropy.io import fits
import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))

from tica.tica import Calibration, CCD_File
from wcs_build.step2_mkwcs import fit_wcs_in_imgdir


calibration = Calibration(calibration_dir =
                          os.path.join(DIR,'../calibration_10min'))



inputs = ['input/tess2021041025908-s0035-1-1-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-1-2-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-1-3-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-1-4-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-2-1-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-2-2-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-2-3-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-2-4-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-3-1-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-3-2-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-3-3-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-3-4-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-4-1-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-4-2-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-4-3-0205-s_ffir.fits.gz',
          'input/tess2021041025908-s0035-4-4-0205-s_ffir.fits.gz']

output_checks = ['output_checks/tess2021041025908-s0035-1-1-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-1-2-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-1-3-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-1-4-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-2-1-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-2-2-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-2-3-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-2-4-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-3-1-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-3-2-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-3-3-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-3-4-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-4-1-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-4-2-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-4-3-0205-s_ffir.cal.fits.gz',
                 'output_checks/tess2021041025908-s0035-4-4-0205-s_ffir.cal.fits.gz']


def run_wcs_fit(sector, cam, ccd, refdata, imglist):

    outdir='./'
    fitdeg = 6
    save = False
    debug = 0

    #print('in func!', sector, cam, ccd, refdata, imglist)
    fit_wcs_in_imgdir(sector, cam, ccd,
                      refdata, imglist,
                      outdir, fitdeg,
                      save, debug)

def check_outputs(check_file_name):
    gold_f = fits.open(check_file_name)

    #assumes test files were written to reg_test directory, which is
    #one level above output_checks
    test_f = fits.open( os.path.basename( check_file_name ))

    
    assert (gold_f[0].data == test_f[0].data).all()

    for key in gold_f[0].header.keys():
        if key == 'COMMENT':
            if 'calibration applied at' in gold_f[0].header[key][0]:
                continue
        if 'CHECKSUM' in key:
            continue
        try:
            assert gold_f[0].header[key] == test_f[0].header[key]
        except AssertionError:
            print('gold file {} = {}'.format(key, gold_f[0].header[key]))
            print('test file {} = {}'.format(key, test_f[0].header[key]))
            raise

    #also, check that there are not new keywords in test_f
    for key in test_f[0].header.keys():
        if key == 'COMMENT':
            if 'calibration applied at' in test_f[0].header[key][0]:
                continue
        if 'CHECKSUM' in key:
            continue
        try:
            assert gold_f[0].header[key] == test_f[0].header[key]
        except AssertionError:
            print('gold file {} = {}'.format(key, gold_f[0].header[key]))
            print('test file {} = {}'.format(key, test_f[0].header[key]))
            raise

for ii in range(len(inputs)):
    ccd1 = CCD_File(inputs[ii], calibration=calibration)
    ccd1.write_calibrate()

    cam = int(ccd1.header['CAMERA'])
    ccd = int(ccd1.header['CCD'])
    ref_file = 'ref_stars/hlsp_tica_tess_ffi_s0035-cam{}-ccd{}_tess_v01_cat.h5'.format(cam,ccd)
    run_wcs_fit( 35, cam, ccd, ref_file, [os.path.basename( output_checks[ii] )] )

    check_outputs(output_checks[ii])

    os.remove(os.path.basename( output_checks[ii] ))

print('Check complete! No assertion errors')
