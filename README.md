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

2D bias and flat field calibration models are distributed by (dvc.org)[DVC].  To retrieve the model, install DVC on your system, and use 

```
cd calibration_models
dvc pusll calibration_<exptime>
```

where `<exptime>` corresponds to whatever exposure your FFIs are (30 minutes in Sectors 1--26, 10 minutes in Sectors 37--52).

## Quick Start

The work horse script for calibrating reaw TESS data is `bin/tica-cal-ccd2ccd`.  This script is installed by default, and can be run with `--help` to  see an explanation of the options available.  

An example bash script to run tica on raw FFIs downloaded from MAST is in `bin/tica-calibrate-spoc`.  A help option is also available, but the user inputs the location of the FFIs as argument one and the location of the calilbration models as argument two.  The script will mkdir directories that organize the calibrated files by cam and ccd in the current directory.

#### calibrating a set of files

Use `cal-ccd2ccd` (in `bin` directory) to run example set of files:

     cd examples/
     ../bin/cal-ccd2ccd testcal.txt

You can also put the output in a different directory

     cd examples/
     mkdir output
     ../bin/cal-ccd2ccd --outdir output testcal.txt

## Support

Please submit bug reports and feature requests

    https://tessgit.mit.edu/tess/tica/issues

## Development

The code for TICA is at

    https://tessgit.mit.edu/tess/tica.git

## Copyright and License

TICA is Copyright(c) 2017, Massachusetts Institute of Technology (MIT)

TICA is distributed under the terms of GPLv3


