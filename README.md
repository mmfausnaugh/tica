# TESS Image CAlibration (TICA)

The TESS Image Calibration (TICA) is a Python module for removing instrumental effects from raw TESS images.

There are currently five steps:

 1. 2d bias correction: Remove fixed-pattern noise.
 2. 1d bias correction: Remove the baseline, time-dependent amplifier signal.
 4. Gain and Linearity: Correct non-linearity response and convert from ADU to photoelectrons.
 3. Smear correction:  Remove small contamination from shutterless transfer to the frame-store.
 4. Flat-field correction: Correct for different sensitivities of each pixel.


## Installation

Either use `pip3 install tica` or clone the repository:

  ```
  git clone https://tessgit.mit.edu/tica/tica.git
  cd tica
  python setup.py sdist
  pip3 install -e .
  ```

You can instead add the `tica` directory to you `PYTHONPATH` environment variable and `tica/bin` to your `PATH`.

2D bias and flat field calibration models are distributed via an application called [DVC](htpps://www.dvc.org).  To retrieve the model, install DVC on your system, and use 

```
cd calibration_models
dvc pull flat_combine
dvc pull twodbias_<exptime>
```

where `<exptime>` corresponds to whatever exposure your FFIs are, `30min` for Sectors 1-26, `10min` for Sectors 37-52.

## Quick Start

The work horse script for calibrating reaw TESS data is `bin/tica-cal-ccd2ccd`.  This script is installed by default, and can be run with `--help` to  see an explanation of the options available.  

An example bash script to run tica on raw FFIs downloaded from MAST is in `bin/tica-calibrate-spoc`.  A help option is also available, but fundamentally the user inputs the location of the FFIs and the location of the calilbration models:

```
mkdir tica_outputs
cd tica_outputs
tica-calibrate-spoc /absolutepath/to/raw/data ~/python/tica/calibration_10min
```

The script will make directories that organize the calibrated files by cam and ccd in the current directory.

`bin/tica-calibrate-spoc` is also installed in the users PATH by default, and likely covers 99% of use cases.  Users can also write their own scripts to call `tica-cal-ccd2ccd` or run `tica-cal-ccd2ccd` directly from the comamnd line.

### Setting the calibration directory

***At the moment, users point to the appropriate calibration directory.  Flatfields are always the same, and don't matter.  2Dbiases just need to be rescaled for exptime.  We could make TICA automatically detect this (at some load on Users who must then store all files), or have the user set an environement variable, or some combination of the three***


## To Do

1. Document the remaining scripts---what should the public actually see given that they can't use everythign?
2. Check install and that it runs.
3. Do more details need to be provided about DVC/how to use DVC?
4.  Need some method for time dependent twodbias corrections (as of April 2021)
5.  Would like some tutorial on tica data structures, which are useful for analyzing and manipulation TESS FFIs. (similar to my mapspec module)
