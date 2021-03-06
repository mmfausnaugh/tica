## Maintenance Plan

- The primary developers for TICA are Michael Fausnaugh (MF) and Christopher Burke (CB).  
    - We are research scientists at the MIT Kavli Institute.  We work directly with the TESS Payload Operation Center (POC), Science Processing Operation Center (SPOC), and TESS Science Office (TSO).
- We plan to respond to bug reports and any issues that cause the core functionality of TICA to fail.
    - We will also check for changes based on updates to python and core packages such as `numpy`, `astropy`, `tess-point`, etc.  In most cases, we will pin TICA to old versions of these packages rather than make updates to the TICA code.
- MF is responsible for the calibration algorithms and structure of the python package and scripts.  CB is responsible for the WCS implementation.
- TICA uses semantic versioning, with the initial release being v1.0.2.  
- Questions, issues, and bug reports should be sent to tesshelp@bigbang.gsfc.nasa.gov.
- TESSHELP automatically logs the email in the TESSHELP JIRA project, along with responses from TESS team members.
- TESSHELP is hosted by MIT, but team members at GSFC and MAST also have access.
   - We will ask the GI office and MAST to advertise the tesshelp email address as the avenue for questions, bug reports, and feedback.
    - MF and CB will identify TESSHELP issues that require code changes or development effort, and open an issue on the TICA github to track progress.
- Development is on a best effort basis, and so we cannot guarantee a timeline for responding to TESSHELP questions or closing issues on github.  That said, for urgent issues we will endeavor to provide a response and prototype merge request within 1 month of opening an issue on github.

## How to Contribute

- Users are free to clone, fork, and track this repository.
- Typically, we will not approve pull request from outside developers.  This is because we have designed TICA from the perspective of the instrument team and we do not wish to introduce changes unrelated to effects in the detector.  
- In exceptional cases, a pull request may be appropriate.  If you feel this is the case, you can submit an email to tesshelp@bigbang.gsfc.nasa.gov, with the following details:
    - An explanation of the proposed changes.
    - Evidence that the changes make a substantial differences to the TICA data products or the way in which users interact with the TICA code. 
    - An explanation of why the change would be valuable for the scientific community as a whole.
- If we feel that the proposed changes are exceptionally useful or warranted, we will invite you to submit a pull request.
