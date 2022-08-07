#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue March 9 12:24 2021

@author: jtu
"""

#####################################################################################
#####################################################################################
import os
import numpy as np
import math

import matplotlib.pyplot as plt
#from matplotlib import rcParams
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from scipy import stats
from scipy import interpolate

if __name__ == '__main__':
    rad=math.pi/180.0

    fid=open('/home/jtu/mycode/imit3d/conveloc.dat', 'r')
    #for i in range(10):
    #    fid.readline()

    ff=fid.readline()
    a=ff.split()
    mnstr=a[2]
    if len(mnstr) ==1:
        mnstr='0'+mnstr
    daystr=a[3]
    if len(daystr) == 1:
        daystr='0'+daystr
    utstr=a[5]
    if len(utstr)==1:
        uststr='0'+utstr
    minstr=a[6]
    if len(minstr)==1:
        minstr='0'+minstr
    secstr=a[7]
    if len(secstr)==1:
        secstr='0'+secstr

    DateUT=a[1]+'-'+mnstr+'-'+daystr+'  UT '+utstr+':'+minstr+':'+secstr

    ff=fid.readline()
    a=ff.split()
    imf='IMF Bx: '+a[1]+' By: '+a[3]+' Bz: '+a[5]

    ff=fid.readline()
    a=ff.split()
    sw=r'SW N$_p$: '+a[2]+' Vx: '+a[5]+' Vy: '+a[8]+' Vz: '+a[11]

    ff=fid.readline()
    a=ff.split()
    UseAL=a[2]

    ff=fid.readline()
    a=ff.split()
    ALindex=float(a[2])

    if UseAL == 'Yes':
        imf=imf+' AL '+str(int(ALindex))

    fid.readline()

    fl=fid.readlines()
    fid.close()

    ndata=len(fl)

    mLT=[]
    mlon=[]

    mlat=[]
    pot=[]
    vth=[]
    vph=[]

    data_arr=[]

    k=1
    for i in range(ndata):
        a=fl[i].split()

        LT=float(a[1])

        if i == 0:
            mLT.append(LT)
            mlon.append(15.0*LT)

        if i > 0 and LT != mLT[k-1]:
            mLT.append(LT)
            mlon.append(15.0*LT)

            data_arr.append([mlat, pot, vth, vph])

            k=k+1

            mlat=[]
            pot=[]
            vth=[]
            vph=[]

            mlat.append(float(a[0]))
            pot.append(float(a[3]))
            vth.append(float(a[5]))
            vph.append(float(a[6]))
        else:
            mlat.append(float(a[0]))
            pot.append(float(a[3]))
            vth.append(float(a[5]))
            vph.append(float(a[6]))

        if i == ndata-1:
            data_arr.append([mlat, pot, vth, vph])

    nx=len(data_arr)+1
    ny=int(len(data_arr[0][0])/2)

    xn=np.zeros((nx,ny))
    yn=np.zeros((nx,ny))
    xs=np.zeros((nx,ny))
    ys=np.zeros((nx,ny))

    potn=np.zeros((nx,ny))
    pots=np.zeros((nx,ny))

    vxn=np.zeros((nx,ny))
    vyn=np.zeros((nx,ny))

    vxs=np.zeros((nx,ny))
    vys=np.zeros((nx,ny))

    mlon.append(mlon[0])

    for k in range(nx):
        theta=(mlon[k]-90.0)*rad

        for j in range(ny):
            if k < nx-1:
                r=90.0-data_arr[k][0][j]
            else:
                r=90.0-data_arr[0][0][j]

            xn[k][j]=r*math.cos(theta)
            yn[k][j]=r*math.sin(theta)

            xs[k][j]=r*math.cos(theta)
            ys[k][j]=r*math.sin(theta)

            if k < nx-1:
                potn[k][j]=data_arr[k][1][j]

                vxn[k][j]=( data_arr[k][2][j]*math.cos(theta) \
                           -data_arr[k][3][j]*math.sin(theta))
                vyn[k][j]=( data_arr[k][2][j]*math.sin(theta) \
                           +data_arr[k][3][j]*math.cos(theta))

                pots[k][j]=data_arr[k][1][j+ny]

                vxs[k][j]=( data_arr[k][2][j+ny]*math.cos(theta) \
                           -data_arr[k][3][j+ny]*math.sin(theta))
                vys[k][j]=( data_arr[k][2][j+ny]*math.sin(theta) \
                           +data_arr[k][3][j+ny]*math.cos(theta))

    # periodic values
    for j in range(ny):
        potn[nx-1][j]=potn[0][j]
        pots[nx-1][j]=pots[0][j]

        vxn[nx-1][j]=vxn[0][j]
        vyn[nx-1][j]=vyn[0][j]

        vxs[nx-1][j]=vxs[0][j]
        vys[nx-1][j]=vys[0][j]

    #find maximum magnitude of the velocity vectors
    vmax_k=np.zeros(nx)
    vmaxn=np.zeros(ny)
    vmaxs=np.zeros(ny)
    for k in range(nx):
        for j in range(ny):
            vmaxn[j]=math.sqrt(vxn[k][j]**2+vyn[k][j]**2)
            vmaxs[j]=math.sqrt(vxs[k][j]**2+vys[k][j]**2)

        vmax_k[k]=max(max(vmaxn), max(vmaxs))

    vmax=max(vmax_k)
    print(vmax, 'm/s')

    xL=[0.05,0.51,0.05,0.51]
    yB=[0.05,0.05,0.46,0.46]
    wh=0.43
    ht=0.39

    fig=plt.figure(figsize=(1.6, 1.7), dpi=400)

    axes=[]
    for i in range(4):
        axes.append(fig.add_axes([xL[i], yB[i], wh, ht]))
        fig.add_axes(axes[i])

        axes[i].axis('off')

    ring=[10.,20.,30.,40.]
    rlen=len(ring)
    clen=len(mlon)
    xr=np.zeros((rlen, clen))
    yr=np.zeros((rlen, clen))

    for j in range(rlen):
        for k in range(clen):
            xr[j][k]=ring[j]*math.cos(mlon[k]*rad)
            yr[j][k]=ring[j]*math.sin(mlon[k]*rad)

        if j < rlen-1:
            axes[0].plot(xr[j], yr[j], lw=0.2, ls='dashed', color='black')
            axes[1].plot(xr[j], yr[j], lw=0.2, ls='dashed', color='black')
            axes[2].plot(xr[j], yr[j], lw=0.2, ls='dashed', color='black')
            axes[3].plot(xr[j], yr[j], lw=0.2, ls='dashed', color='black')
        else:
            axes[0].plot(xr[j], yr[j], lw=0.2, ls='solid', color='black')
            axes[1].plot(xr[j], yr[j], lw=0.2, ls='solid', color='black')
            axes[2].plot(xr[j], yr[j], lw=0.2, ls='solid', color='black')
            axes[3].plot(xr[j], yr[j], lw=0.2, ls='solid', color='black')

    mring=max(ring)
    for i in range(4):
        axes[i].plot([-mring-1.2, -mring], [0.0,0.0], lw=0.2, color='black')
        axes[i].plot([mring, mring+1.1], [0.0,0.0], lw=0.2, color='black')
        axes[i].plot([0.0,0.0], [-mring-1.2, -mring], lw=0.2, color='black')
        axes[i].plot([0.0,0.0], [mring, mring+1.1], lw=0.2, color='black')

        axes[i].text(-4, -mring-5.3, '0 MLT', fontsize=2)
        axes[i].text(mring+3, -1, '6', fontsize=2)
        axes[i].text(-2,mring+2.5, '12', fontsize=2)
        axes[i].text(-mring-6.8, -1, '18', fontsize=2)

        if i < 2:
            axes[i].text(-10., -ring[3]-3, '50', fontsize=2)
            axes[i].text(-8.5, -ring[2]-3, '60', fontsize=2)
            axes[i].text(-7.0, -ring[1]-3, '70', fontsize=2)
            axes[i].text(-5.5, -ring[0]-3, '80', fontsize=2)

    axes[2].text(0, mring+8, 'Northern Hemisphere', fontsize=2.3, ha='center')
    axes[3].text(0, mring+8, 'Southern Hemisphere', fontsize=2.3, ha='center')

    axes[2].contour(xn, yn, potn, levels=20, colors='tab:gray', linewidths=0.2)
    axes[3].contour(xs, ys, pots, levels=20, colors='tab:gray', linewidths=0.2)

    axes[2].quiver(xn, yn, vxn, vyn, color='blue')
    axes[3].quiver(xs, ys, vxs, vys, color='blue')

    axes[0].contourf(xn, yn, potn, levels=20, cmap='rainbow')
    axes[1].contourf(xs, ys, pots, levels=20, cmap='rainbow')

    axes[2].text(-mring,mring+26, DateUT, fontsize=2.3, ha='left')
    axes[2].text(-mring,mring+20, imf, fontsize=2.3, ha='left')
    axes[2].text(-mring,mring+14, sw, fontsize=2.3, ha='left')

    plt.show()
