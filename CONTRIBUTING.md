
## Maintenance Plan

-The primary developers for TICA are Michael Fausnaugh (MF) and Christopher Burke (CB).  We are research scientists at the MIT Kavli Institute.  We work directly with the TESS Payload Operation Center (POC), Science Processing Operation Center (SPOC), and TESS Science Office (TSO).
-We plan to respond to bug reports and any issues that cause the core functionality of TICA to fail.
--We will also check for changes based on updates to python and core packages such as numpy, astropy, etc.  In some cases, we may pin TICA to old versions of these packages, rather than make updates to the TICA code.
-MF is responsible for the calibration algorithms and general structure/invocation of the package and scripts.  CB is responsible for the WCS implementation.
-Tica uses semantic versioning, with the initial release being v1.0.2.  
--The first minor version increment (to 1.1.0) will be an update to the cloud storage of the calibration models, so tht a Google Account is not required to run TICA.  This change is expected in the first half of 2022.
-- The second minor release will include calibration models for 200 second FFIs.  This is expected to be in the second half of 2022. 
-Issues/questions should be reported to tesshelp@bigbang.gsfc.nasa.gov, which automatically logs the email (and responses) to the TESSHELP JIRA project.
-TESSHELP is hosted by MIT, but team members at GSFC and MAST have access.
--We will ask the GI office and MAST to advertise tesshelp as the avenue for questions, bug reports, and feedback.
--MF and CB will identify TESSHELP issues that require code changes/development, and open an issue on the TICA github to track progress.
-Development is on a best effort basis, and so we cannot guarantee a time line for responding to TESSHELP questions or closing issues on github.  That said, for urgent issues we will endeavor to provide a response and prototype version of a merge request within 1 month of opening an issue on github.

## How to Contribute

-Users are free to clone, fork, and track this repository.
-Users may submit pull requests, as long as they provide
--An explanation of the proposed changes.
--Evidence that the changes make substantial differences to the TICA products and would be valuable for the scientific community.
-We will review the code for any merge requests, and require that the regression tests pass on the branch submitted in the pull request.
