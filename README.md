# TESS Image CAlibration (TICA)

The TESS Image Calibration (TICA) is a Python module for removing instrumental effects from raw TESS images.

There are currently five steps:

 1. 2d bias correction: Remove fixed-pattern noise (a property of the electronics), using models from the POC
 2. 1d bias correction: Remove the base level (time-dependent) amplifier signal, using the mean value of virtual columns.
 4. Gain and Linearity: Correct non-linearity response and convert from ADU to photoelectrons.
 3. Smear correction:  Remove small contamination from shutterless transfer to the frame-store.
 4. Flat-field correction: Correct for different sensitivities of each pixel, using models from the POC.

Steps 1 is a single model, steps 2, 3 and 4 are array operations calculated for each image, and step 5 is a single model.

## Installation

Either use `pip3 install tica` or clone the repository:

  ```
  git clone https://tessgit.mit.edu/tica/tica.git
  cd tica
  python setup.py sdist
  pip3 install -e .
  ```

You can instead add the `tica` directory to you `PYTHONPATH` environment variable and `tica/bin` to your path***true???**  

2D bias and flat field calibration models are distributed by (DVC)[htpps://www.dvc.org].  To retrieve the model, install DVC on your system, and use 

```
cd calibration_models
dvc pull calibration_<exptime>
```

where `<exptime>` corresponds to whatever exposure your FFIs are (30 minutes in Sectors 1-26, 10 minutes in Sectors 37-52).

## Quick Start

The work horse script for calibrating reaw TESS data is `bin/tica-cal-ccd2ccd`.  This script is installed by default, and can be run with `--help` to  see an explanation of the options available.  

An example bash script to run tica on raw FFIs downloaded from MAST is in `bin/tica-calibrate-spoc`.  A help option is also available, but the user inputs the location of the FFIs as argument one and the location of the calilbration models as argument two.  The script will mkdir directories that organize the calibrated files by cam and ccd in the current directory.

`bin/tica-calibrate-spoc` is also installed in the users PATH by default, and likely covers 99% of usecases.

### Setting the calibration directory

***At the moment, users point to the appropriate calibration directory.  Flatfields are always the same, and don't matter.  2Dbiases just need to be rescaled for exptime.  We could make TICA automatically detect this (at some load on Users who must then store all files), or have the user set an environement variable, or some combination of the three***


## To Do

1. Document the remaining scripts---what should the public actually see given that they can't use everythign?
2. Check install and that it runs.
3. Do more details need to be provided about DVC/how to use DVC?
4. Documentation of the algorithms and models?  Probably best to put that in Fausnaugh et al.	