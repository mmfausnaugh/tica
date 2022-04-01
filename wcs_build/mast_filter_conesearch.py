
import numpy as np
import sys
import json
import time
try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
    from urllib.parse import urlencode as dict_urlencode
    from urllib.request import urlopen
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve
    from urllib import urlencode as dict_urlencode
    from urllib import urlopen
try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib
    

## [Mast Query]
def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content
## [Mast Query]

def mast_filter_conesearch(starRa, starDec, radius, minTmag, maxTmag):
    # Do a MAST cone search around the ra and dec to get the nearby stars
    startTime = time.time()
    request = {'service':'Mast.Catalogs.Filtered.Tic.Position.Rows', 
               'params':{'columns':'ID,ra,dec,Tmag,Kmag,GAIAmag,pmRA,pmDEC,disposition', 
                         'filters':[ 
                             {'paramName':'Tmag',
                              'values':[{'min':minTmag, 'max':maxTmag}]},
                             {'paramName':'pmRA',
                              'values':[{'min':-100.0, 'max':100.0}]}, 
                             {'paramName':'pmDEC',
                              'values':[{'min':-100.0, 'max':100.0}]}],
                         'ra':'{:10.5f}'.format(starRa),
                         'dec':'{:10.5f}'.format(starDec),
                         'radius':'{:10.7f}'.format(radius/3600.0) },
               'format':'json', 'removenullcolumns':False}
    while True:    
        headers, outString = mastQuery(request)
        try:
            outObject = json.loads(outString)
            if outObject['status'] != 'EXECUTING':
                break
        except:
            print('Problem at MAST. Resting and trying again')
            time.sleep(10)
        if time.time() - startTime > 30:
                print('Working...')
                startTime = time.time()
        time.sleep(5)

    try:
        outObject2 = []
        for x in outObject['data']:
            if x['disposition'] == 'DUPLICATE' or x['disposition'] == 'SPLIT':
                continue
            else:
                outObject2.append(x)
        ticList = np.array([x['ID'] for x in outObject2], dtype=np.int64)
        ticRas = np.array([x['ra'] for x in outObject2], dtype=np.float)
        ticDecs = np.array([x['dec'] for x in outObject2], dtype=np.float)
        ticTmags = np.array([x['Tmag'] for x in outObject2], dtype=np.float)
        ticKmags = np.array([x['Kmag'] for x in outObject2], dtype=np.float)
        ticGmags = np.array([x['GAIAmag'] for x in outObject2], dtype=np.float)
        ticpmRA = np.array([x['pmRA'] for x in outObject2], dtype=np.float)
        ticpmDec = np.array([x['pmDEC'] for x in outObject2], dtype=np.float)
        
    except:
        # Try rerunning search
        while True:    
            headers, outString = mastQuery(request)
            try:
                outObject = json.loads(outString)
                if outObject['status'] != 'EXECUTING':
                    break
            except:
                print('Problem at MAST. Resting and trying again')
                time.sleep(20)
            if time.time() - startTime > 30:
                    print('Working...')
                    startTime = time.time()
            time.sleep(5)
            
        try:
            ticList = np.array([x['ID'] for x in outObject['data']], dtype=np.int64)
            ticRas = np.array([x['ra'] for x in outObject['data']], dtype=np.float)
            ticDecs = np.array([x['dec'] for x in outObject['data']], dtype=np.float)
            ticTmags = np.array([x['Tmag'] for x in outObject['data']], dtype=np.float)
            ticKmags = np.array([x['Kmag'] for x in outObject['data']], dtype=np.float)
            ticGmags = np.array([x['GAIAmag'] for x in outObject['data']], dtype=np.float)
            ticpmRA = np.array([x['pmRA'] for x in outObject['data']], dtype=np.float)
            ticpmDec = np.array([x['pmDEC'] for x in outObject['data']], dtype=np.float)
        except:
            print('Tried MAST cone search twice and failed. Exiting')
            exit()
        

    return ticList, ticRas, ticDecs, ticTmags, ticKmags, ticGmags, ticpmRA, ticpmDec

