#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:20:06 2019

@author: jtu
"""
import os
import numpy as np
import math
from scipy import stats
import tkinter as tk

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as figtk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTk as figtk

try:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as navitk
except ImportError:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as navitk

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#implement th edefault matplotlib key bindings
import matplotlib.ticker as ticker
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#####################################################################################
#####################################################################################
def latitude_altitude_contour(wd, canv, fig, toolbar, latf, altf, gridNumEnt, varib, \
    ptype, var_index, add_text, fname, paramList):

    j0=int(gridNumEnt[1].get())
    j1=int(gridNumEnt[4].get())
    i0=int(gridNumEnt[0].get())
    i1=int(gridNumEnt[3].get())

    lat=latf[j0:j1]
    alt=altf[i0:i1]

    print("\nPlot contour of "+paramList[var_index]+" ...")

    colormap=plt.cm.get_cmap('rainbow')

#------------ make plots on the tkinter window -----------------------------
    #clear figure if it exists
    if fig != None:
        fig.clf()
        
    #clear canvas if it exists
    if canv != None:
        canv.get_tk_widget().destroy()

    fig=Figure(figsize=(2,1.5),dpi=400)

    canv = figtk(fig, master=wd)
    canv.draw()
    canv.get_tk_widget().place(x=0,y=60)

    if toolbar != 0:
        toolbar.destroy()

    toolbar=navitk(canv, wd)
    toolbar.update()

    def on_key_press(event):
        key_press_handler(event, canv, toolbar)
        canv.mpl_connect("key_press_event", on_key_press)

    #set size, resolution, and position of axes
    xL=0.12
    yB=0.12
    wh=0.87
    ht=0.76

    axes = fig.add_axes([xL, yB, wh, ht])
    fig.add_axes(axes)

    #set line width of axes
    for ax in ['left','bottom','right','top']:
        axes.spines[ax].set_linewidth(0.3)

# set axis, labels, color bar etc for left three panels
    axes.tick_params(axis='both',which='major',length=1.3,width=0.3,pad=1, \
                    labelsize=2.7,left=True,right=True,bottom=True,top=True, \
                    direction='out')
    axes.tick_params(axis='both',which='minor',length=0.8,width=0.3,pad=1, \
                    labelsize=0,left=True,right=True,bottom=True,top=True, \
                    direction='out')
    axes.xaxis.set_minor_locator(AutoMinorLocator())
    axes.yaxis.set_minor_locator(AutoMinorLocator())
    #axes.set_xticks(np.arange(-90,90, 15))
    #axes.set_xlabel(['-90','-75','-60','-45','-30','-15','0','15','30','45','60','75','90'])

    axes.set_ylabel('Altitude (km)', fontsize=2.8, labelpad=0.7)
    axes.set_xlabel ('Latitude (degree)', fontsize=2.8, labelpad=0.7)
    axes.text(0.01,1.02,add_text,ha='left',va='bottom',fontsize=2.6, \
        transform=axes.transAxes)

    if var_index < 7 or (var_index>=20 and var_index<=26) or (var_index>=40 and var_index <= 43):
        axes.text(1.02,0.97,'Log Scale',ha='left',va='bottom',fontsize=2.5,transform=axes.transAxes)

    #--------------------- plot contour -------------------------------
    sc=axes.contourf(lat, alt, varib, levels=30, cmap=colormap)

    # ---- vertical color bar ----
    cbar=fig.colorbar(sc, ax=axes, shrink=0.9, pad=0.03, aspect=18)
    cbar.ax.tick_params(labelsize=2.4, length=0.85, width=0.4, pad=0.4)
    cbar.ax.set_ylabel(paramList[var_index], fontsize=2.5, labelpad=1.1)
    cbar.outline.set_visible(False)

    print("Contour plot completed!")

    return canv, fig, toolbar