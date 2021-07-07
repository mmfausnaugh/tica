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
             'bin/tica-trim_ccds2ffi',

             'bin/tica-calibrate_spoc',
             'bin/tica-calibrate_tso',
             'bin/tica-step1_allcam_allccd.sh',
             'bin/tica-step2_allcam_allccd.sh'
         ],
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'matplotlib',
        'os',
        'sys',
        'glob',
        'logging',
        'multiprocessing',
        'argparse',
        'datetime',
        'time',

)    


