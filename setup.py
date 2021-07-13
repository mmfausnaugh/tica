from setuptools import setup, find_packages
with open('VERSION','r') as infile:
    version = infile.read()
    
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
             'bin/tica-ffi2ccds',
             'bin/tica-ccds2ffi',
             'bin/tica-trim-ccds2ffi',

             'bin/tica-calibrate-spoc',
             'bin/tica-calibrate-tso',
             'bin/tica-wcs-step1-allcam-allccd',
             'bin/tica-wcs-step2-allcam-allccd',

             'wcs_build/step1_get_refimg_ctrlpts.py',
             'wcs_build/step2_mkwcs.py',

         ],
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'matplotlib',
        'logging',
        'h5py',
        'photutils',
        'statsmodels',
        'gwcs'
        ]

)    

