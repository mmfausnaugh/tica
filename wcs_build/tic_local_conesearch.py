#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Created on Mon Dec 26 13:40 2022

Uses psycopg to query the local postgres tic, for production on MIT
serveres

author:  Michael Fausnaugh

"""

import numpy as np
import os
import sys
import json
import time

import psycopg2

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))

from wcs_build.tic_postgres_param import  tic_host, tic_db, tic_user, tic_port, password
    

def tic_local_conesearch(starRa, starDec, radius, minTmag, maxTmag):
    #dop  a cone search around the ra and dec to get the nearby stars
    #depends on configuration of TIC and postgres on MIT servers

    con = psycopg2.connect(
        "host='{}' dbname='{}' user='{}' port='{}' password='{}'".format(
            tic_host, tic_db, tic_user, tic_port, password
        ))

    cur = con.cursor()

    #example query string from /pcds/ops/pochome/orbits/orbit-129/guides/log

    # Executing query "select id,ra,dec,pmra,pmdec,tmag,e_pmra,e_pmdec from ticentries where spoint(radians(ra),radians(dec)) @ scircle '< (2.24024786248, -0.0240763683039), 0.296192195877 >' and tmag between -5 and 21.0;"

    rad_ra = starRa*np.pi/180.
    rad_dec = starDec*np.pi/180.
    #radius is passed in as arcsec, but submit to postgres as radius
    ##querystr="select id,ra,dec,tmag,kmag,gaiamag,pmra,pmdec,disposition"\
    ##    " from ticentries where"\
    ##    " spoint(radians(ra),radians(dec)) @ scircle '< ({},{}), {} >'"\
    ##    " and tmag between {} and {}"\
    ##    " and pmra between -100.0 and 100.0"\
    ##    " and pmdec between -100.0 and 100.0 ; ".format(rad_ra, rad_dec, 
    ##                                                radius/3600*np.pi/180.,
    ##                                                minTmag, maxTmag)
    querystr="select id,ra,dec,tmag,kmag,gaiamag,pmra,pmdec,disposition"\
        " from ticentries where"\
        " q3c_radial_query(ra,dec, {},{}, {}) "\
        " and tmag between {} and {}"\
        " and pmra between -100.0 and 100.0"\
        " and pmdec between -100.0 and 100.0 ; ".format(starRa, starDec,
                                                    radius/3600.,
                                                    minTmag, maxTmag)

    cur.execute(querystr)
    output = cur.fetchall()
    output = np.array(output)


    ##if disposition is not set to None, then remove the sources
    m = output[:,8] == None
    output = output[m]
    idx = np.argsort(output[:,0])
    output = output[idx]

    ticList  = output[:,0].astype(np.int64)
    ticRas   = output[:,1].astype(np.float64)
    ticDecs  = output[:,2].astype(np.float64)
    ticTmags = output[:,3].astype(np.float64)
    ticKmags = output[:,4].astype(np.float64)
    ticGmags = output[:,5].astype(np.float64)
    ticpmRA  = output[:,6].astype(np.float64)
    ticpmDec = output[:,7].astype(np.float64)

    return ticList, ticRas, ticDecs, ticTmags, ticKmags, ticGmags, ticpmRA, ticpmDec

if __name__ == "__main__":
    out = tic_local_conesearch( 180.0 , 45.0, 8640.0, 7.5, 10)
    for o in out:
        print(o)

