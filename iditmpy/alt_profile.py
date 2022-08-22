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
def plot_altitude_profile(wd, canv, fig, toolbar, altf, gridNumEnt, varib, ptype, \
                          var_index, add_text, fname, paramList):

    i0=int(gridNumEnt[0].get())
    i1=int(gridNumEnt[3].get())

    alt=altf[i0:i1]

    print("\nPlot altitude profile of "+paramList[var_index]+" ...")

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
    yB=0.10
    wh=0.80
    ht=0.76

    axes = fig.add_axes([xL, yB, wh, ht])
    fig.add_axes(axes)

    #set line width of axes
    for ax in ['left','bottom','right','top']:
        axes.spines[ax].set_linewidth(0.3)

# set axis, labels, color bar etc for left three panels
    axes.tick_params(axis='both',which='major',length=1.3,width=0.3,pad=1, \
                    labelsize=2.7,left=True,right=True,bottom=True,top=True, \
                    direction='in')
    axes.tick_params(axis='both',which='minor',length=0.8,width=0.3,pad=1, \
                    labelsize=0,left=True,right=True,bottom=True,top=True, \
                    direction='in')
    axes.xaxis.set_minor_locator(AutoMinorLocator())
    axes.yaxis.set_minor_locator(AutoMinorLocator())

    if (var_index  < 7) or (var_index>=12 and var_index<=18) or (var_index>=32 and var_index<=35):
        axes.set_xscale('log')

    axes.set_ylabel('Altitude (km)', fontsize=3, labelpad=0.7)
    axes.set_xlabel (paramList[var_index], fontsize=3, labelpad=0.7)
    axes.text(0.01,1.01,add_text,ha='left',va='bottom',fontsize=3, \
        transform=axes.transAxes)

    #----- plot  ------
    if var_index < 34 or var_index > 36:
        axes.plot(varib, alt, lw=0.3, color='blue')
        axes.scatter(varib, alt, s=0.01, color='black')

    elif var_index == 34:
        colors=['black','red','blue','purple','yellow','green','cyan']
        labeltxt=['O+','H+','He+','O2+','N2+','NO+','N+']
        #axes.set_xlim([10,1.0e4])
        for s in range(7):
            axes.plot(varib[s][0:], alt, lw=0.3, color=colors[s], label=labeltxt[s])

        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.93,0.99), \
            frameon=False)

    elif var_index == 35:
        colors=['black','red','blue','purple','yellow','green','cyan']
        labeltxt=['O','H','He','O2','N2','NO','N']
        axes.set_xlim([0.1,1.0e17])
        for l in range(7):
            axes.plot(varib[l][0:], alt, lw=0.3, color=colors[l], label=labeltxt[l])
    
        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.93,0.99), \
            frameon=False)

    elif var_index == 36:
        colors=['black','red','blue']
        labeltxt=[r'B$_r$',r'B$_\theta$',r'B$_\phi$']
        for l in range(3):
            axes.plot(varib[l][0:], alt, lw=0.3, color=colors[l], label=labeltxt[l])
    
        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.93,0.99), \
            frameon=False)

    elif (var_index == 70):
        colors=['black','red','blue','purple','yellow']
        labeltxt=['O','H','He','O2','N2']
        for l in range(5):
            axes.plot(varib[l][0:], alt, lw=0.3, color=colors[l], label=labeltxt[l])
    
        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.10,0.99), \
            frameon=False)

    elif (var_index == 71 or var_index==73):
        axes.plot(varib, alt, lw=0.3, color='blue')

    elif (var_index == 72):
        colors=['black','red','blue','purple','yellow','green']
        labeltxt=['O+','H+','He+','O2+','N2+','NO+']
        #axes.set_xlim([0.1,1.0e17])
        for l in range(6):
            axes.plot(varib[l][0:], alt, lw=0.3, color=colors[l], label=labeltxt[l])
    
        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.93,0.99), \
            frameon=False)

    elif (var_index == 74):
        colors=['black','red','blue','purple','yellow','green','cyan']
        labeltxt=['e - O+','e - O2+','e - N2+','e - H+','e - He+','e - NO+', 'e - neutrals']
        for l in range(7):
            axes.plot(varib[l][0:], alt, lw=0.3, color=colors[l], label=labeltxt[l])
    
        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.90,0.99), \
            frameon=False)

    elif (var_index >= 75 and var_index<=80):
        colors=['black','red','blue','purple','yellow','green','cyan']
        if var_index==75:
            labeltxt=['O+ - e','O+ - O2+','O+ - N2+','O+ - H+','O+ - He+','O+ - NO+', 'O+ - neutrals']
        elif var_index==76:
            labeltxt=['O2+ - e','O2+ - O+','O2+ - N2+','O2+ - H+','O2+ - He+','O2+ - NO+', 'O2+ - neutrals']
        elif var_index==77:
            labeltxt=['N2+ - e','N2+ - O+','N2+ - O2+','N2+ - H+','N2+ - He+','N2+ - NO+', 'N2+ - neutrals']
        elif var_index==78:
            labeltxt=['H+ - e','H+ - O+','H+ - O2+','H+ - N2+','H+ - He+','H+ - NO+', 'H+ - neutrals']
        elif var_index==79:
            labeltxt=['He+ - e','He+ - O+','He+ - O2+','He+ - N2+','He+ - H+','He+ - NO+', 'He+ - neutrals']
        else:
            labeltxt=['NO+ - e','NO+ - O+','NO+ - O2+','NO+ - N2+','NO+ - H+','NO+ - He+', 'NO+ - neutrals']

        for l in range(7):
            axes.plot(varib[l][0:], alt, lw=0.3, color=colors[l], label=labeltxt[l])
    
        axes.legend(loc='upper center', fontsize=2, bbox_to_anchor=(0.86,0.99), \
            frameon=False)

    print("Plot completed!")

    return canv, fig, toolbar