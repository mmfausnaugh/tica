import numpy as np
import h5py

import sys


fname1 = sys.argv[1]
fname2 = sys.argv[2]

f1 = h5py.File(fname1,'r')
f2 = h5py.File(fname2,'r')

#check they have the same data sets
flag = True
for k in f1.keys():
    try:
        assert k in f2
    except AssertionError:
        flag = False
        print('Key "{}" foun in {} but not in {}'.format(
            k,fname1, fname2))

for k in f2.keys():
    try:
        assert k in f1
    except AssertionError:
        flag = False
        print('Key "{}" foun in {} but not in {}'.format(
            k,fname2, fname1))


#and the same attriburtes
#for k in f1.attrs.keys():
#    assert k in f2.attrs.keys()

#for k in f2.attrs.keys():
#    assert k in f1.attrs.keys()

d1 = f1['tics'][:]
d2 = f2['tics'][:]
print('N tics in {} = {}'.format(f1,len(d1)))
print('N tics in {} = {}'.format(f2,len(d2)))

m1 = np.in1d(f1['tics'][:], f2['tics'][:])
idx1 = np.argsort(f1['tics'][:][m1])
m2 = np.in1d(f2['tics'][:], f1['tics'][:][m1])
idx2 = np.argsort(f2['tics'][:][m2])

d1 = d1[m1][idx1]
d2 = d2[m2][idx2]

mask = abs(d1-d2) > 0

for k in f1.keys():
    try:
        d1 = f1[k][:][m1][idx1]
    except KeyError:
        print('Key "{}" not in {}'.format(k,fname1))
        flag = False
        continue
    try:
        d2 = f2[k][:][m2][idx2]
    except KeyError:
        print('Key "{}" not in {}'.format(k,fname2))
        flag = False
        continue

   # print(np.c_[d1[mask],d2[mask]])
   # print(np.c_[d1,d2])
    #print(np.array_equal(d1,d2))
    try:
        assert np.array_equal(d1,d2)
    except AssertionError:
        flag = False
        print('Data not the same for key "{}"'.format(k))
        print('data that is not identical:')
        print('tics1','d1','diff','d2','tics2')
        print(np.c_[ f1['tics'][:][m1][idx1][ abs(d1-d2) > 0], 
                     d1[ abs(d1-d2) > 0],
                     d1[ abs(d1-d2)  > 0] - d2[ abs(d1-d2) > 0],
                     d2[ abs(d1-d2) > 0],
                     f2['tics'][:][m2][idx2][ abs(d1-d2) > 0] ])

    

if flag == True:
    print("no assertion errors, all data is the same")
