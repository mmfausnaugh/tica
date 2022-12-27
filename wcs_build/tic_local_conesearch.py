#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Created on Mon Dec 26 13:40 2022

Uses psycopg to query the local postgres tic, for production on MIT
serveres

author:  Michael Fausnaugh

"""

import numpy as np
import sys
import json
import time

import psycopg2

from tic_postgres_param import  tic_host, tic_db, tic_user, tic_port
    

def tic_local_conesearch(starRa, starDec, radius, minTmag, maxTmag):
    #dop  a cone search around the ra and dec to get the nearby stars
    #depends on configuration of TIC and postgres on MIT servers

    con = psycopg2.connect(
        "host='{} dbname='{}' user='{}' port='{}'".format(
            tic_host, tic_db, tic_user, tic_port
        ))

    cur = con.cursor()

    #example query string from /pcds/ops/pochome/orbits/orbit-129/guides/log

    # Executing query "select id,ra,dec,pmra,pmdec,tmag,e_pmra,e_pmdec from ticentries where spoint(radians(ra),radians(dec)) @ scircle '< (2.24024786248, -0.0240763683039), 0.296192195877 >' and tmag between -5 and 21.0;"

    rad_ra = starRa*np.pi/180.
    rad_dec = starDec*np.pi/180.
    querystr="selectid,ra,dec,tmag,kmag,gaiamag,pmra,pmdec,disposition"\
        " from ticentries where"\
        "spoint(radians(ra),radians(dec)) @ scircle '< ({},{}), {} >'"\
        "and tmag between {} and {}"\
        "and pmra between -100.0 and 100.0"\
        "and pmdec between -100.0 and 100.0".format(rad_ra, read_dec, 
                                                    radius/3600.*np.pi/180.,
                                                    minTmag, maxTmag)

    cur.execute(querystr)
    

    output = cur.fetchall()

    print(output)
    print(np.shape(output))

    #while True:    
    #    headers, outString = mastQuery(request)
    #    try:
    #        outObject = json.loads(outString)
    #        if outObject['status'] != 'EXECUTING':
    #            break
    #    except:
    #        print('Problem at MAST. Resting and trying again')
    #        time.sleep(10)
    #    if time.time() - startTime > 30:
    #            print('Working...')
    #            startTime = time.time()
    #    time.sleep(5)
    #
    #try:
    #    outObject2 = []
    #    for x in outObject['data']:
    #        if x['disposition'] == 'DUPLICATE' or x['disposition'] == 'SPLIT':
    #            continue
    #        else:
    #            outObject2.append(x)
    #    ticList = np.array([x['ID'] for x in outObject2], dtype=np.int64)
    #    ticRas = np.array([x['ra'] for x in outObject2], dtype=np.float)
    #    ticDecs = np.array([x['dec'] for x in outObject2], dtype=np.float)
    #    ticTmags = np.array([x['Tmag'] for x in outObject2], dtype=np.float)
    #    ticKmags = np.array([x['Kmag'] for x in outObject2], dtype=np.float)
    #    ticGmags = np.array([x['GAIAmag'] for x in outObject2], dtype=np.float)
    #    ticpmRA = np.array([x['pmRA'] for x in outObject2], dtype=np.float)
    #    ticpmDec = np.array([x['pmDEC'] for x in outObject2], dtype=np.float)
    #    
    #except:
    #    # Try rerunning search
    #    while True:    
    #        headers, outString = mastQuery(request)
    #        try:
    #            outObject = json.loads(outString)
    #            if outObject['status'] != 'EXECUTING':
    #                break
    #        except:
    #            print('Problem at MAST. Resting and trying again')
    #            time.sleep(20)
    #        if time.time() - startTime > 30:
    #                print('Working...')
    #                startTime = time.time()
    #        time.sleep(5)
    #        
    #    try:
    #        ticList = np.array([x['ID'] for x in outObject['data']], dtype=np.int64)
    #        ticRas = np.array([x['ra'] for x in outObject['data']], dtype=np.float)
    #        ticDecs = np.array([x['dec'] for x in outObject['data']], dtype=np.float)
    #        ticTmags = np.array([x['Tmag'] for x in outObject['data']], dtype=np.float)
    #        ticKmags = np.array([x['Kmag'] for x in outObject['data']], dtype=np.float)
    #        ticGmags = np.array([x['GAIAmag'] for x in outObject['data']], dtype=np.float)
    #        ticpmRA = np.array([x['pmRA'] for x in outObject['data']], dtype=np.float)
    #        ticpmDec = np.array([x['pmDEC'] for x in outObject['data']], dtype=np.float)
    #    except:
    #        print('Tried MAST cone search twice and failed. Exiting')
    #        exit()
    #    
    #
    #return ticList, ticRas, ticDecs, ticTmags, ticKmags, ticGmags, ticpmRA, ticpmDec

if __name__ == "__main__":
    tic_local_conesearch( 180.0 , 45.0, 2.39, 7.5, 10):    
