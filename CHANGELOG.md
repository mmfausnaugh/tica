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
- Output directory structure of `bin/tica-calibrate-spoc` and `bin/tica-wcs-spoc`.
- Regression tests updated to match python version 3.10 and associated dependencies.
- Regression tests now have a 1.e-8 relative tolerance, and a README logging previous tests.
### Added
### Removed


## [1.0.1] - 2021-12-03
Refactored code to install as a python package with regression tests.  Primarily based on branch `faus_setup_dist`.
### Changed
- setup.py now resolves dependencies
- can install tica libraries and bin/ a python package
- Help functions for bash scripts in bin/
### Added
- Enforces calibration of 30min/10min FFIs with appropriate 2D bias model.
- Binary data now tracked with dvc.
- Regression test to track changes caused by algorithms or updates to dependencies
### Removed
