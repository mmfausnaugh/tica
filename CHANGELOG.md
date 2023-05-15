# Changelog

Human-readable summary of major changes and features/data that are added or removed.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
### Added
### Removed

## [1.4.0] - 2022-5-8
Code auotomatically finds a good reference image for step1_wcs.  Other aspects of logging and checking to automate.

	
## [1.3.0] - 2022-5-8
Parallelized wcs step2, changed step1 to do local database querries
### Changed
- Invalid logging directories are made rather than crashing the scripts.
- Fixed versions of dependencies, for reproducibility
### Added
- Added `parallel` argument to step2_mkwcs.py.
- Added tic_local_conesearch.py, and changed step1 of wcs to use local server database. Now querries 16 CCDs in parallel.
- added psycopg2 (v2.9.5) dependency for local query
- Added function in tica.CCD to calculate mode of science pixels.

	
## [1.2.0] - 2022-09-09
Added 200sec models, and keyword specifying which half of an orbit FFIs came from.
### Changed
- Fixed issue with check on total integration times.
- Changed tica-calibrate-tso to deal withi new directory structure from POC.
- Changed details of output logs and production table.
### Added
- dvc file for twodbias_200s.dvc
- added relevant 200 sec exposures where needed.  
- added keyword ORB_SEG in headers for o1a, o1b, o2a, o2b, 
- added tica-cal-ffi2ccds command line argument --orb_seg to add ORB_SEG
- tica-calibrate-tso can now do orbit segments independently.  Includes catches for weird combinations (i.e., o2 is set but requesting o1a)
- added option to force .fits file, even if input is gzipped (default is to write as .gz if input is .gz)
### Removed

## [1.1.2] - 2022-06-07
Fixed some issues related to reprocessing old TICA products.
### Changed
- EXPTIME header keyword now depends on sector.
- Remved `header.extend` from wcs step2, so duplicate keywords should never be added.
### Added
### Removed

## [1.1.1] - 2022-04-27
Extended same capabilities on the production table logging.
### Changed
- Changed FLXWIN header keyword to blkHlfCent variable, as per issue 14.
- Updated tica-production-table to deal with logs that capture reruns.
- Updated tica-production-table to deal with logs from multiple orbits/segments (as in em2).	
- Updated regression test to check ext 1 data and headers.  Renamed ref star files to v02.
### Added
### Removed

## [1.1.0] - 2022-04-13
Improves configuration control, adds a file that tracks HLSP processing by sector.  Added header kewords and info in reference .h5 files, defaults to "fixed" apertures for WCS.  Adds Ref star flux and bkg estimates to ext 1 table of individual images.
### Changed
- wingFAC and contrastFAC in `step1_get_refimg_ctrlpts.py` are now command line arguments, and are stored in the .h5 files.
- Number of processers in `Pool` are command line arguments in `bin/tica-cal-ffi2ccds` and similar.
- Reg test now ignores patch number in image header `TICAVER` keyword.
- Removed np.int from `step1_get_refimg_ctrlpts.py` and `step1_mkwcs.py`, which was deprecated in numpy 1.20.
- WCS diagnostic info in the .h5 files now appends from files already with a WCS, rather than skipping those files and leaving the diagnostic info as zero.
- WCS .h5 star files now saves the Row/Col used in the reference, and forces the same row col when deriving the WCS unless the user disables this at the command line.
- Duplicates in the TIC are filtered out of the MAST query for reference stars.
- Updated regression test for fixed apertures and new ref_star files.
- Changed WCS plotting defaults.
- Changed stale SPOC header keywords to 'HISTORY' keywords.
- Added measured flux and backgrounds of WCS stars to extension 1.
### Added
- Header keywords for RMS scatter of faint stars in WCS fits.
- Script to see if the data in two reference star h5 files is the same.
- Script to make a "production table," which tracks configuration and important stats from each sector.
- Production table in Markdown format.
	
## [1.0.2] - 2022-02-04
Initial release on public github.  Primarily based on branch `faus_fix_requirements`.
### Changed
- Updated to python 3.10 and associated dependencies.
- DVC tracking of flatfields is now for individual flatfield models.
- Updated README and help menus based on MIT POC code review.
- `bin/tica-calibrate-spoc` no longer organizes outputs by camera and CCD, in order to match inputs of `bin/tica-wcs-spoc`.
- Regression tests updated to match python version 3.10 and associated dependencies.
- Regression test Header Keywords now have a 1.e-8 relative tolerance, in order to accomodate differences in Intel and M1 processors.
- Regression test directory now has a README to  log previous tests.
- Regression test now includes a verbose mode, to either display or silence warnings.
### Added
- Added MIT server to DVC remote, and set to default.
### Removed


## [1.0.1] - 2021-12-03
Refactored code to install TICA as a python package.  Added regression tests.  Primarily based on branch `faus_setup_dist`.
### Changed
- `setup.py` now resolves dependencies.
- The user can now install TICA libraries as a python package and access scripts in `tica/bin` from the command line.
- The bash scripts in `bin/` now have help functions.
### Added
- `tica/tica/tica.py` now enforces calibration of 30min/10min FFIs with the appropriate 2D bias model.
- Binary data in this repository are now tracked with DVC.
- Added the `reg_test` directory, which includes a python script to run regression test in order to test for changes in TICA data products.
### Removed
