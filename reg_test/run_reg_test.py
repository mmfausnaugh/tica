import numpy as np
from astropy.io import fits
import os
import sys
import argparse

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))

from tica.tica import Calibration, CCD_File
from wcs_build.step2_mkwcs import fit_wcs_in_imgdir


calibration = Calibration(calibration_dir =
                          os.path.join(DIR,'../calibration_10min'))

calibration_wrong = Calibration(calibration_dir =
                                os.path.join(DIR,'../calibration_30min'))



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
    fixApertures=True
    save = False
    debug = 0
    #print('in func!', sector, cam, ccd, refdata, imglist)
    fit_wcs_in_imgdir(sector, cam, ccd,
                      refdata, imglist,
                      outdir, fitdeg,
                      fixApertures,
                      save, debug)

def check_outputs(check_file_name, verbose=False):
    gold_f = fits.open(check_file_name)

    #assumes test files were written to reg_test directory, which is
    #one level above output_checks
    test_f = fits.open( os.path.basename( check_file_name ))

    
    assert (gold_f[0].data == test_f[0].data).all()
    #print(gold_f[1].data[0:10])
    #print(test_f[1].data[0:10])
    
    #in python, np.nan == np.nan evaluts to False
    #so break this up to avoid nans
    m = gold_f[1].data == test_f[1].data
    assert all(gold_f[1].data[m] == test_f[1].data[m])
    assert all(np.isnan(test_f[1].data['FLX'][~m] ))
    assert all(np.isnan(test_f[1].data['BKG'][~m] ))
    assert all(np.isnan(gold_f[1].data['FLX'][~m] ))
    assert all(np.isnan(gold_f[1].data['BKG'][~m] ))

    assert all( gold_f[1].data['TIC'][~m]      == test_f[1].data['TIC'][~m] )
    assert all( gold_f[1].data['FLXCOL'][~m]   == test_f[1].data['FLXCOL'][~m] )
    assert all( gold_f[1].data['FLXROW'][~m]   == test_f[1].data['FLXROW'][~m] )
    assert all( gold_f[1].data['FLXVALID'][~m] == test_f[1].data['FLXVALID'][~m] )


    #header keywords can be different in the WCS for Intel vs. M1 chips,
    #as well as different versions of numpy

    #we raise a warning if they differ by less then 1 part in 10^8
    #and we raise an error if they differ by more
    
    warning_flag = 0.0
    for key in gold_f[0].header.keys():
        if key == 'COMMENT':
            if 'calibration applied at' in gold_f[0].header[key][0]:
                continue
        if 'CHECKSUM' in key:
            continue
        #if 'RMSF' in key or 'RMSBF' in key:
        #    continue
        if 'TICAVER' in key:
            #this is just too painful, no need to change the v number
            #in the test unless there are other errors

            #continue
            #assert that the major.minor versions
            #are the same.  OK if Patches are different,
            #gold_ver = gold_f[0].header[key].split('.')
            #test_ver = test_f[0].header[key].split('.')
            #assert gold_ver[0] == test_ver[0]
            #assert gold_ver[1] == test_ver[1]
            continue
        try:
            assert gold_f[0].header[key] == test_f[0].header[key]
        except AssertionError as e:
            r_diff = abs(gold_f[0].header[key] - test_f[0].header[key])/ \
                     gold_f[0].header[key]
            if verbose:
                print('WARNING!!!')
                print('Differences in gold file {}'.format(check_file_name))
                print('gold file {} = {}'.format(key, gold_f[0].header[key]))
                print('test file {} = {}'.format(key, test_f[0].header[key]))
                print('relative difference is {:3.2e}'.format(r_diff))
            if r_diff < 1.e-8:
                warning_flag = 1.0
                continue
            else:
                print(key, gold_f[0].header[key], test_f[0].header[key])
                print(e)
                raise

    #also, check that there are not new keywords in test_f
    for key in test_f[0].header.keys():
        if key == 'COMMENT':
            if 'calibration applied at' in test_f[0].header[key][0]:
                continue
        if 'CHECKSUM' in key:
            continue
        #if 'RMSF' in key or 'RMSBF' in key:
        #    continue
        if 'TICAVER' in key:
            continue
        try:
            assert gold_f[0].header[key] == test_f[0].header[key]
        except AssertionError:
            #print('WARNING')
            #print('gold file {} = {}'.format(key, gold_f[0].header[key]))
            #print('test file {} = {}'.format(key, test_f[0].header[key]))
            r_diff = abs(gold_f[0].header[key] - test_f[0].header[key])/ \
                     gold_f[0].header[key]
            if r_diff < 1.e-8:
                continue
            else:
                raise

    #check extension 1
    for key in gold_f[1].header.keys():
        if 'CHECKSUM' in key:
            continue
        try:
            assert gold_f[1].header[key] == test_f[1].header[key]
        except AssertionError as e:
            r_diff = abs(gold_f[1].header[key] - test_f[1].header[key])/ \
                     gold_f[1].header[key]
            if verbose:
                print('WARNING!!!')
                print('Differences in gold file {}'.format(check_file_name))
                print('gold file {} = {}'.format(key, gold_f[1].header[key]))
                print('test file {} = {}'.format(key, test_f[1].header[key]))
                print('relative difference is {:3.2e}'.format(r_diff))
            if r_diff < 1.e-8:
                warning_flag = 1.0
                continue
            else:
                print(key, gold_f[1].header[key], test_f[1].header[key])
                print(e)
                raise
    for key in test_f[1].header.keys():
        if 'CHECKSUM' in key:
            continue
        try:
            assert gold_f[1].header[key] == test_f[1].header[key]
        except AssertionError as e:
            r_diff = abs(gold_f[1].header[key] - test_f[1].header[key])/ \
                     gold_f[1].header[key]
            if verbose:
                print('WARNING!!!')
                print('Differences in gold file {}'.format(check_file_name))
                print('gold file {} = {}'.format(key, gold_f[1].header[key]))
                print('test file {} = {}'.format(key, test_f[1].header[key]))
                print('relative difference is {:3.2e}'.format(r_diff))
            if r_diff < 1.e-8:
                warning_flag = 1.0
                continue
            else:
                print(key, gold_f[1].header[key], test_f[1].header[key])
                print(e)
                raise


    return warning_flag

if __name__ == "__main__":

    description_string = """

    Run TICA calibration and WCS on raw SPOC files in `input`
    directory, and check for changes against standards in
    `output_checks` directory.

    """

    parser = argparse.ArgumentParser(description=description_string)
    parser.add_argument('--verbose', action='store_true',
                        help="Print warnings.  Warnings are raised for"
                        " keywords that differ by less than 1 part in 10^8."
                        " Larger differences will cause the test to fail.")

    args = parser.parse_args()
    warning_flag = 0
    for ii in range(len(inputs)):
        #check that calibration model exptime safe gaurds work
        try:
            ccd1 = CCD_File(inputs[ii], calibration=calibration_wrong)
        except AssertionError:
            ccd1 = CCD_File(inputs[ii], calibration=calibration)
        ccd1.write_calibrate()

        cam = int(ccd1.header['CAMERA'])
        ccd = int(ccd1.header['CCD'])
        ref_file = 'ref_stars/hlsp_tica_tess_ffi_s0035-cam{}-ccd{}_tess_v02_cat.h5'.format(cam,ccd)
        #ref_file = 'ref_stars_new/reftica_s35_{}-{}.h5'.format(cam,ccd)
        run_wcs_fit( 35, cam, ccd, ref_file, [os.path.basename( output_checks[ii] )] )

        warning_flag += check_outputs(output_checks[ii], verbose=args.verbose)

        os.remove(os.path.basename( output_checks[ii] ))

    if warning_flag == 0:
        print('\n\nCheck complete! No assertion errors\n\n')
    elif warning_flag > 0 and args.verbose == False:
        print('\n\nCheck complete! No assertion errors\n\n'
              'Some warnings were detected.  Warnings are raised\n'
              ' for Header Keywords that differ from the output_checks\n'
              ' by less than 1 part in 10^8.  We do not consider \n'
              ' these differences significant and so consider the\n'
              ' regression test to have passed.\n\n'
              'Consider rerunning this test with "--verbose".')
    else:
        print('\n\nCheck complete! No assertion errors\n\n'
              'Some warnings were detected.  This might be do due to\n'
              ' different versions of python or Intel vs. Apple processors.\n'
              ' We do not consider these differences significant\n'
              ' and so consider the regression test to have passed.\n\n'
              'See tica/reg_test/README for details.')
    
