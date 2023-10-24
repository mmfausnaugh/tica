from setuptools import setup, find_packages
with open('VERSION','r') as infile:
    version = infile.read().strip()
    
setup(
    name='tica',
    version = version,
    author='Michael Fausnaugh',
    author_email='faus@mit.edu',
    url='https://tessgit.mit.edu/tica/tica',
    license='MIT',
    packages=find_packages(include=['tica',
                                    'btjd']),
    scripts=['bin/tica-cal-ccd2ccd', 
             'bin/tica-cal-ffi2ccds',
             'bin/tica-cal-ffi2ffi', 
             'bin/tica-ccd-report',  
             'bin/tica-production-table',  
             'bin/tica-ffi2ccds',
             'bin/tica-ccds2ffi',
             'bin/tica-trim-ccds2ffi',

             'bin/tica-calibrate-spoc',
             'bin/tica-calibrate-tso',
             'bin/tica-wcs-step1-allcam-allccd',
             'bin/tica-wcs-step2-allcam-allccd',
             'bin/tica-wcs-spoc',

             'wcs_build/step1_get_refimg_ctrlpts.py',
             'wcs_build/step2_mkwcs.py',

             'bin/tica-check-delivery',
             'bin/tica-stage-delivery',

             'bin/tica-darktime',
             'bin/tica-doall-automate'

         ],
    install_requires=[
        'numpy==1.22.1',
        'scipy==1.7.3',
        'astropy==5.0',
        'matplotlib==3.5.1',
        'h5py==3.6.0',
        'tess-point>=0.8.0',
        'gwcs==0.15.0',
        'psycopg2-binary==2.9.5'
        ]

)    

