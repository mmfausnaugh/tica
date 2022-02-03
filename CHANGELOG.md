# Changelog

Human-readable summary of major changes and features/data that are added or removed.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
### Added
### Removed


## [1.0.2] - 2022-02-04
Initial release on public github.  Primarily based on branch `faus_fix_requirements`.
### Changed
- Updated to python 3.10 and associated dependencies.
- DVC tracking of flatfields is now for individual flatfield models.
- Updated README and help menus based on MIT POC code review.
- `bin/tica-calibrate-spoc` no longer organizes outputs by camera and CCD, in order to match inputs of `bin/tica-wcs-spoc`..
- Regression tests updated to match python version 3.10 and associated dependencies.
- Regression test Header Keywords now have a 1.e-8 relative tolerance, in order to accept differences in Intel and M1 processors.
- Regression test directory now has a README to  log previous tests.
### Added
### Removed


## [1.0.1] - 2021-12-03
Refactored code to install TICA as a python package.  Added regression tests.  Primarily based on branch `faus_setup_dist`.
### Changed
- `setup.py` now resolves dependencies.
- The user can now install TICA libraries a python package and access scripts in `tica/bin` from the command line.
- The bash scripts in `bin/` now have help functions.
### Added
- `tica/tica/tica.py` now Enforces calibration of 30min/10min FFIs with appropriate 2D bias model.
- Binary data in this repository are now tracked with dvc.
- Added `reg_test`, which includes a python script to run regression test in order to test for changes caused by algorithms or updates to dependencies.
### Removed
