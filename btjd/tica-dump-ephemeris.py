import numpy as np
import spiceypy as s
import os
import sys
import argparse
from datetime import datetime


#Designed in a modular way, so that all the spice stuff is separate from
#the TICA and WCS stuff.
#
#output of this is tess_ephem.txt. Requires inputting a NAIF ephemeris file for SPICE
#uses spiceypy to get TESS position relative to barycent
#
#writes or appends to tess_ephem.txt.  That file has X,Y,Z, VX,VY, VZ every 15 minutes
#
#tes_ephem.txt is getting long!  what is a better option? 
#well, only 200,000 rows so far, and that is 5 years. 17MB, probably this is OK until 2030, at least

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--ephemeris", required=True, type=str, nargs='+',
      help="One or more ephemeris files to load into spice")
    parser.add_argument("--outdir", default="./",type=str)

    return parser.parse_args()


tofffin = 577364091.1874759      # Calculated using s.str2et("2018-04-18 22:53:42.001873 UTC"), which is in turn from poctime -t  0 -i SPM3 -o UTC
tofftjd = 1227.45475911# Calculated using poctime -i SPM3 -o TJD -t 0

def utc2tjd(utc_str):


    #astart is in UTC...
    t = s.str2et(utc_str)
    
    #print('raw from spiceypy',t)
    #print('check RAW', t - tofffin , (t - tofffin )/86400.0)

    #just lines up the RAW time in Space with the TJD for FFI_INDEX = 0
    return (t - tofffin )/86400.0 + tofftjd

def tjd2et(time_in):
    return (time_in - tofftjd)*86400.0  +  tofffin 

    
def main():
    args = parse_args()

    s.kclear()
    s.boddef("TESS",-95)

    for fn in args.ephemeris:
        s.furnsh(fn)

    outfile = os.path.join(args.outdir,'tess_ephem.txt')

    if not os.path.isfile(outfile):
        t0 = 1250.00080076  #TJD of launch + 23 days, to make spice happy
    else:
        t = np.genfromtxt(outfile,usecols=(0))
        t0 = t[-1]


    d = datetime.utcnow()
    utc_t1 = d.strftime('%Y-%m-%d-%H:%M:%S')
    print('in UTC',utc_t1)
    t1 = utc2tjd(utc_t1)
    print('convert tjd',t1)

    #spacing is about 15 minutes
    time_in = np.r_[t0:t1:0.010412]
    if len(time_in) < 2:
        sys.exit()
    t_spiceypy = tjd2et(time_in)
    #t_spiceypy = s.str2et(utc)

    #s.spkezr returns 6-element array of position and velocity 
    #and a scalar light travel time
    print(t_spiceypy)
    tess_position = np.array([s.spkezr('TESS', tm, 'J2000', 'NONE', '0')[0] for tm in t_spiceypy ])


    
    if not os.path.isfile(outfile):
        with open(outfile,'w') as fout:
            for ii,row in enumerate(tess_position):
                fout.write('{:15.6f} {:15e} {:15e} {:15e} {:8.4f} {:8.4f} {:8.4f}\n'.format(time_in[ii], row[0],row[1],row[2],row[3],row[4],row[5]))
    else:
        with open(outfile,'a') as fout:
            for ii,row in enumerate(tess_position):
                fout.write('{:15.6f} {:15e} {:15e} {:15e} {:8.4f} {:8.4f} {:8.4f}\n'.format(time_in[ii], row[0],row[1],row[2],row[3],row[4],row[5]))

if __name__ == "__main__":
    main()
