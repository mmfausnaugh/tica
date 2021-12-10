# TESS Image CAlibration (TICA)

The TESS Image Calibration (TICA) module is a python module for removing instrumental effects from raw TESS images.

There are currently six steps:

 1. 2d bias correction: Remove fixed-pattern noise.
 2. Scalar bias correction: Remove the time-dependent amplifier pedestal.
 3. Gain: Convert from ADU to photoelectrons.
 4. Linearity: Correct the non-linear detector response.
 5. Smear correction:  Remove contamination from the shutterless transfer to the frame-store.
 6. Flat-field correction: Correct for the non-uniform response of each pixel.

TICA also includes a module to find bright, relatively isolated stars in TESS FFIs, measure their centroids, and derive World Coordinate Solutions.  

## Installation

First, setting up a virtual environment for the TICA installation is considered best practice.

Next, you will need to clone this repository to get the calibration files, so you should install from here:

  ```
  git clone https://tessgit.mit.edu/tica/tica.git
  cd tica
  python setup.py sdist
  pip3 install -e .
  ```

Instead, you can  add the `tica` directory to you `PYTHONPATH` environment variable and `tica/bin` to your `PATH`.

See the `requirements.txt` file for a list of required packages.  Besides typical scientific computing packages (`numpy`, `scipy`, `matplotlib`, and `astropy`), TICA requires `gwcs` (https://github.com/spacetelescope/gwcs) and `tess-point` (https://github.com/christopherburke/tess-point).  The `setup.py` file should install any missing packages.

2D bias and flat field calibration models are distributed by an application called [DVC](https://www.dvc.org) (Data Version Control).  DVC uses metadata files in the TICA git repository to track different versions of the calibration models.  The calibration models themselves are currently stored in a GDrive account associated with MIT.  

To retrieve the calibration models, first install DVC using the directions on [the DVC website](https://dvc.org/doc/install).  Then run the following commands and follow the prompts:

```
cd calibration_models
dvc pull flat_combine
dvc pull twodbias_<exptime>
```

where `<exptime>` corresponds to whatever exposure your FFIs are, `30min` for Sectors 1-26, `10min` for Sectors 37-55.

***If you install `dvc` using `python` or `pip`, you also need to install the packages necessary for accessing GDrive.  This can be accomplished as follows:***

```
pip3 install dvc
pip3 install dvc[gdrive]
```

## Quick Start

The workhorse script for calibrating raw TESS data is `bin/tica-cal-ccd2ccd`.  This script is installed in your working environment and can be run with `--help` to  see an explanation of the available options.  

An example bash script to run TICA on raw FFIs downloaded from MAST is `bin/tica-calibrate-spoc`.  This script is also installed in your working environment and a help option is  available (run with `--help` or with no arguments). 

The user gives the bash script the location of the FFIs and the location of the calibration models:

```
mkdir tica_outputs
cd tica_outputs
tica-calibrate-spoc input_dir=/absolute/path/to/raw/data CALDIR=~/python/tica/calibration_10min
```

The script will make directories that organize the calibrated files by camera and CCD in the current directory.  The "raw data" are SPOC `*ffir*` files, available on MAST (https://archive.stsci.edu/tess/bulk_downloads.html).

`bin/tica-calibrate-spoc` likely covers 99% of use cases.  Users can also write their own scripts to call `tica-cal-ccd2ccd`, or they can run `tica-cal-ccd2ccd` directly from the command line.

Several other scripts are also installed by default.  These scripts are used to calibrate FFIs for MIT's [Quick Look Pipeline](https://archive.stsci.edu/hlsp/qlp) at the Payload Operations Center/TESS Science Office (POC/TSO).  Note that the POC/TSO FFIs are formatted differently than archival data products available at MAST.

### Setting the calibration directory

***Users must point to the appropriate calibration directory when running `tica-calibrate-spoc` or `tica-cal-ccd2ccd`.   The code will raise an error if the exposure time in the data files does not match the specified calibration directory.***

## Updating and Reverting Calibration Models

With DVC, you can easily update to the latest calibration models.  You can also revert to old versions of the models if you need to reproduce prior results.  

To check for updates, run `git pull` in the top-level `tica` directory.  If new calibration models are available, then the `calibration_models/*.dvc` files will be updated.  If this is the case, `cd` to the `calibration_models` directory and run `dvc pull`.  

Calibration models will be updated no more than once per year, and probably much less frequently.

To revert to an old model, first use `git checkout` to get the version of the `.dvc` file that corresponds to the model that you need.  Then run `dvc checkout`.  DVC handles access to the cloud, and caches the downloads so that you can easily switch between different versions of the models.  See the [DVC docs](https://dvc.org/doc/start) for more information and tutorials.

There are directories at the top level of the repository for each calibration model (30 minute or 10 minute data).  Those directories link to the data in `calibration_models`. When you run TICA, you specify a single calibration directory, which will load all of the models that you need to run the code (based on what is checked out in `calibration_models`).

Calibration models also exist for 20-second and 2-minute data, but calibrating raw pixels from SPOC Target Pixel Files is not supported at this time.

## Data Structures and Calibration Algorithms

Algorithms used for pixel calibrations can be found in `tica/tica/tica.py`.  In particular, `CCD.calibrate` is the function that actually produces the calibrated data.  

Some users may find the data structures in `tica.py` useful for their own scripts.  We recommend using `bin/tica-cal-ccd2ccd` or `bin/tica-cal-ffi2ccds` as a guide for scripting with TICA data structures.


## Regression Tests

There is a regression test in `tica/reg_test` that can be used to test any changes to the TICA code.  To run the regression test:

```
cd tica/reg_test
python run_reg_test.py
```

This script will take SPOC `*ffir*` files in `reg_test/input` and apply calibrations and WCSs.  The results of those calibrations are checked against the contents of `reg_test/output_checks`; the test will fail if there are any changes.

To run `run_reg_test.py`, it is necessary to pull the input data, output checks, and tables of WCS stars from the TICA GDrive storage.  DVC is used to retrieve these data in the same way as for the calibration models.  The regression test uses data from Sector 35 and 10 minute calibration models.  

In very rare cases, `output_checks` will be updated if we implement changes that improve the TICA calibration or WCSs.  Changes to the output checks are under version control with the `reg_test/output_checks.dvc` file, and can be updated in the same way as the calibration models.

Note that the regression test does not check `tica/wcs_build/step1_get_refimg_ctrlpts.py`, which generates the table of stars used for WCS fitting.  

## World Coordinate Solutions

There are two WCSs available to end-users: SPOC WCSs in the `*ffir*`/`*ffic*` files, and TICA WCSs in the [High Level Science Products on MAST](https://archive.stsci.edu/hlsp/tica).  There are small systematic differences between the SPOC and TICA WCSs that may be important for users with high astrometric requirements.

The WCS differences can be analyzed by comparing the transforms of stars using the two WCSs.  For TICA, astrometric residuals from the WCSs of the fitted stars are also given in a binary table extension for every FFI.

We recommend that users employ the existing WCSs to analyze these astrometric differences.  However, as a scripting example, we have also provided a bash script that refits the TICA WCSs and overwrites the WCS keywords in the SPOC `*ffir*` files.  This script is `bin/tica-wcs-spoc`.  This bash script is installed in the user's environment and has a help menu (use `--help` or call with no arguments).


## To Do

1. Need some method for time dependent twodbias corrections (as of April 2021)
2. Configuration control, for example, when modifying goodPRF parameters in wcs_step2.