# Changelog

Human-readable summary of major changes and features/data that are added or removed.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
### Added
### Removed

## [1.1.0] - 2022-03-30
Improves configuration control, adds a file that tracks HLSP processing by sector.  Added header kewords and info in reference .h5 files, defaults to "fixed" apertures for WCS.
### Changed
- wingFAC and contrastFAC in `step1_get_refimg_ctrlpts.py` are now command line arguments, and are stored in the .h5 files.
- Number of processers in `Pool` are command line arguments in `bin/tica-cal-ffi2ccds` and similar.
- Reg test now ignores patch number in image header `TICAVER` keyword.
- Removed np.int from `step1_get_refimg_ctrlpts.py` and `step1_mkwcs.py`, which was deprecated in numpy 1.20.
- WCS diagnostic info in the .h5 files now appends from files already with a WCS, rather than skipping those files and leaving the diagnostic info as zero.
- WCS .h5 star files now saves the Row/Col used in the reference, and forces the same row col when deriving the WCS unless the user disables this at the command line.
- Duplicates in the TIC are filtered out of the MAST query for reference stars.
- Updated regression test for fixed apertures and new ref_star files.
### Added
- Header keywords for RMS scatter of faint stars in WCS fits.
- Script to see if the data in two reference star h5 files is the same.

## [1.0.2] - 2022-02-04
Initial release on public github.  Primarily based on branch `faus_fix_requirements`.
### Changed
- Updated to python 3.10 and associated dependencies.
- DVC tracking of flatfields is now for individual flatfield models.
- Updated README and help menus based on MIT POC code review.
- `bin/tica-calibrate-spoc` no longer organizes outputs by camera and CCD, in order to match inputs of `bin/tica-wcs-spoc`..
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
