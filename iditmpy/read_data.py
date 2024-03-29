#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 4 10:45:03 2020

@author: jtu
"""
import os
import numpy as np
import h5py

# input idimt3 arrays
def input_data(fname):
    ##### open and read hdf5 file, and then close the file ########
    print('Reading HDF5 data file ...')
    if (os.path.exists(fname)):
        hf=h5py.File(fname, 'r')
    else:
        print('HDF5 file '+fname+' does not exist')
        return [-1],-1,-1

    key, iditm_arr=next(iter(hf.items()))
    iditm_arr=np.array(iditm_arr)   #array of type iditm_arr[][][][]

    print(hf[key])
    dset = hf[key]

    attr_value=[]
    i=0
    for item in dset.attrs.keys():
        attr_value.append(dset.attrs[item])

    param=[]
    for i in range(len(attr_value[1])):
        param.append(attr_value[1][i])
    for i in range(len(attr_value[2])):
        param.append(attr_value[2][i])

    hf.close()

    print(param)

    fname_B=os.path.dirname(fname)+'/B0.h5'

    hf=h5py.File(fname_B, 'r')

    B0=hf.get('B_BGD')
    B0=np.array(B0)

    print(B0.shape)
    hf.close()

    return iditm_arr, param, B0
    ###############################################################

def input_grids(fpath, a1, a2, a3):
    ################ read altitude, latitude & longitude #################
    print('Reading grid system ...')
    frd=open(fpath+'/grids.dat',"r")

    #skip first two lines
    frd.readline()
    frd.readline()

    alt=np.zeros(a1)
    lat=np.zeros(a2)
    lon=np.zeros(a3)
    cr=np.zeros(a1)
    rdth=np.zeros(a1)
    cos_rsinth=np.zeros((a2,a1))
    rsinth_dph=np.zeros((a2,a1))

    #read altitudes
    for i in range(a1):
        lineR=frd.readline()
        astr=lineR.split()
        alt[i]=float(astr[3])

    frd.readline()

    #read latitudes
    for j in range(a2):
        lineR=frd.readline()
        astr=lineR.split()
        lat[j]=90.0-float(astr[1])

    #reverse latitude to seuquence from south to north
    lat=lat[::-1]

    frd.readline()

    #read longitutdes
    for k in range(a3):
        lineR=frd.readline()
        astr=lineR.split()
        lon[k]=float(astr[1])

    #read cr, r_dth
    """
    for i in range(a1):
        lineR=frd.readline()
        astr=lineR.split()
        cr[i]=float(astr[1])
        rdth[i]=float(astr[2])

    #read cos_rsinth & rsinth_dph
    for j in range(a2):
        for i in range(a1):
            lineR=frd.readline()
            astr=lineR.split()

            cos_rsinth[j][i]=float(astr[2])
            rsinth_dph[j][i]=float(astr[3])
    """

    frd.close()

    return alt, lat, lon, cr, rdth, cos_rsinth, rsinth_dph
    ######################################################################
