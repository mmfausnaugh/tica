#!/usr/bin/env python
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import sys
import glob


def btjd_correction(time_in,ra,dec, ephemeris_data):
    c = 2.99792458e5 #km/s
 
#pass this in only once for tica production   
#    duse = os.path.dirname(os.path.abspath(__file__))
#    t,tess_x,tess_y,tess_z = np.genfromtxt(os.path.join(duse,'tess_ephem.txt'),
#                                           unpack=1,usecols=(0,1,2,3))
    t, tess_x,tess_y,tess_z = ephemeris_data[:,0], ephemeris_data[:,1], \
                              ephemeris_data[:,2], ephemeris_data[:,3]
    tess_vx,tess_vy,tess_vz = ephemeris_data[:,4], ephemeris_data[:,5], \
                              ephemeris_data[:,6]

    interp_x = interp1d(t,tess_x)
    interp_y = interp1d(t,tess_y)
    interp_z = interp1d(t,tess_z)

    ix = interp_x(time_in)
    iy = interp_y(time_in)
    iz = interp_z(time_in)

    #output for image headers
    interp_vx = interp1d(t,tess_vx)
    interp_vy = interp1d(t,tess_vy)
    interp_vz = interp1d(t,tess_vz)
    ivx = interp_vx(time_in)
    ivy = interp_vy(time_in)
    ivz = interp_vz(time_in)

#    tess_position = np.transpose([ix,iy,iz])
    tess_position = np.transpose([iy,ix, iz])

    #ra/dec in radians
    ra_use  = ra*np.pi/180.0
    dec_use = dec*np.pi/180.0
    #turns out this will already be normalized as a unit vector
    source_vector = np.array([np.cos(dec_use)*np.cos(ra_use), 
                              np.cos(dec_use)*np.sin(ra_use), 
                              np.sin(dec_use)])    

    #assuming tess_position is in kilometers, this gives the correction in days
    #seems as if the default tess_ephem file is rotated 90 degress or so?
    dot_prod = np.dot(tess_position, source_vector)
    norm2 = np.sqrt(np.sum(tess_position**2))
    #print('dot product',dot_prod/norm2)
    #print('arccos',np.arccos(dot_prod/norm2)*180./np.pi)

    dtime = np.dot(tess_position, source_vector)/c/86400.0

#    print tess_position
#    print dtime
#    print dtime*24*60
    return time_in + dtime, [ix,iy,iz,ivx,ivy,ivz]
