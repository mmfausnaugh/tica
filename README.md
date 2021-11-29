# TESS Image CAlibration (TICA)

The TESS Image Calibration (TICA) module is a python module for removing instrumental effects from raw TESS images.

There are currently six steps:

 1. 2d bias correction: Remove fixed-pattern noise.
 2. Scalar bias correction: Remove the time-dependent amplifier pedestal.
 3. Gain: Convert from ADU to photoelectrons.
 4. Linearity: Correct the non-linear response.
 5. Smear correction:  Remove small contamination from the shutterless transfer to the frame-store.
 6. Flat-field correction: Correct for the non-uniform response of each pixel.


## Installation

You will need to clone this repository to get the calibration files, so you should install from here:

  ```
  git clone https://tessgit.mit.edu/tica/tica.git
  cd tica
  python setup.py sdist
  pip3 install -e .
  ```

Setting up a virtual environement is considered best practice.

Instead, you can  add the `tica` directory to you `PYTHONPATH` environment variable and `tica/bin` to your `PATH`.

2D bias and flat field calibration models are distributed via an application called [DVC](https://www.dvc.org) (Data Version Control).  DVC uses metadata files in the TICA git repository to track different versions of the calibration models.  The calibration models themselves live in the cloud, currently a GDrive account associated with MIT.  

First, install DVC following the directions on [the DVC website](https://dvc.org/doc/install).  Then issue the following commands and follow the prompts to retrieve the calibration models:

```
cd calibration_models
dvc pull flat_combine
dvc pull twodbias_<exptime>
```

where `<exptime>` corresponds to whatever exposure your FFIs are, `30min` for Sectors 1-26, `10min` for Sectors 37-55.

## Updating Calibration Files

With DVC, you will be able to easily update to the latest calibration models whenever they are available.  You can also revert to old versions if you need to reproduce prior results.  

To check for updates, run `git pull.`  If new calibration models are availabel, `.dvc` files in `calibration_models` will be updated.  If this is the case, move to the `calibration_models` directory, and run `dvc pull`.  

Calibration models will be updated no more than once per year, and we expect that updates will be much less frequent than that.

Reverting to an old model is slightly more complicated.  Roughly speaking, this consists of (1) using `git checkout` to get the version of the `.dvc` file that corresponds to the model that you need, and then (2) running `dvc checkout`.  DVC handles access to the cloud, and caches the downloads so that you can easily switch back and forth between different versions of the models.  See the [DVC docs](https://dvc.org/doc/start) for more information and tutorials.

There are directories at the top level of the repo for each calibration mode (30 minute or 10 minute data).  Those directories link back to the data in `calibration_models`; when you run TICA, you specify a single calibration directory, which will load in the correct models, assuming that you have the desired version checked out in `calibration_models`.

Calibration models also exist for 20 second and 2 minute data, but calibrating raw pixels from SPOC Target Pixel Files is not supported at this time.

## Quick Start

The workhorse script for calibrating raw TESS data is `bin/tica-cal-ccd2ccd`.  This script is installed by default, and can be run with `--help` to  see an explanation of the options available.  

An example bash script to run tica on raw FFIs downloaded from MAST is in `bin/tica-calibrate-spoc`.  A help option is also available for this script (run with `--help` or with no arguments). 

The user inputs the location of the FFIs and the location of the calilbration models:

```
mkdir tica_outputs
cd tica_outputs
tica-calibrate-spoc input_dir=/absolutepath/to/raw/data CALDIR=~/python/tica/calibration_10min
```

The script will make directories that organize the calibrated files by cam and ccd in the current directory.  The "raw data" are SPOC `*ffir*` files, available on MAST (https://archive.stsci.edu/tess/bulk_downloads.html).

`bin/tica-calibrate-spoc` is installed in the user's PATH by default, and likely covers 99% of use cases.  Users can also write their own scripts to call `tica-cal-ccd2ccd` or run `tica-cal-ccd2ccd` directly from the comamnd line.

Note that several other scripts are also installed by default---these are used to calibrate FFIs for MIT's Quick Look Pipeline.  The format of these FFIs are different than archival data products available at MAST, and so these scripts are not useful for the public at this time.

### Setting the calibration directory

***Users must point to the appropriate calibration directory when running `tica-calibrate-spoc`.   The code will stop you if the exposure time in the data files does not match the calibration directory that you specified.***

## Data Structures and Calibration Algorithms

Algorithms used in the calibration itself can be found in `tica/tica.py`.  In particular, `CCD.calibrate` actually produces the calibrated data.  Ambituous users may find the data structures in `tica.py` useful for their own scripts.  We would recommend that they use `bin/tica-cal-ccd2ccd` or `bin/tica-cal-ffi2ccds` as a model.


## To Do

1. Check install and that it runs.
3. Need some method for time dependent twodbias corrections (as of April 2021)
4. Would like some tutorial on tica data structures, which are useful for analyzing and manipulation TESS FFIs. (similar to my mapspec module)
5. Configuration control, for example, when modifying goodPRF parameters in wcs_step2.