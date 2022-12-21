import os
import math
import numpy as np
import h5py
import tkinter as tk

#import plot_panels
import def_classes
import read_data
import alt_profile
import lati_height_cont
import latitude_profile
import long_profile
import long_lat_cont
import divB
import efield
import ele_velocity
import fast_wave_speed
import parameters
import shapiro
import day_number
import subsolar_colat_lon

# subsolar point longitude and then magnetic local time string
def get_MLT_str(year, days, UT, long):
    # subsolar point longitude in deg
    _, subsolar_lon = subsolar_colat_lon.subsolar_colat_lon(year, days, UT)

    # magnetic local time in hours
    MLT = (long-subsolar_lon)/15.0+12.0
    if MLT < 0.0:
        MLT = MLT+24.0
    elif MLT >= 24.0:
        MLT = MLT-24.0

    # convert to string in hr:mm:ss
    MLT_hr=int(MLT)
    MLT_mi=int((MLT-MLT_hr)*60.0)
    MLT_se=int((MLT-MLT_hr-MLT_mi/60.0)*3600.0)
    MLT_str=str(MLT_hr)+':'+str(MLT_mi)+':'+str(MLT_se)

    return MLT_str

class plot_button():
    def __init__(self, rwd, getfile, varbox, plottpbox, varList, gridNumEnt, file_sel):
        self.rwd=rwd
        self.file_sel=file_sel
        self.getfile=getfile
        self.varbox=varbox
        self.plottpbox=plottpbox
        self.varList=varList
        self.gridNumEnt=gridNumEnt

        self.canv=None
        self.toolbar=0
        self.fig=None
        self.data_inputted=0
        self.iditm_arr=None
        self.params=None
        self.oldfilename=None
        self.alt=None
        self.lat=None
        self.lon=None
        self.cr=None
        self.rdth=None
        self.cos_rsinth=None
        self.rsinth_dph=None
        self.B0=None

    def wplot(self):
        #atomic mass for O, O2, N2, H, He, NO, N
        #ms=[16.0,32.0,28.0,1.0,4.0,30.0,14.0]

        fullname=self.getfile.get_filename()
        if fullname == '':
            def_classes.print_message('File not selected')
            return -1
        else:
            if self.oldfilename != fullname:
                self.oldfilename=fullname
                self.data_inputted=0

        wd=self.rwd
        fig=self.fig
        canv=self.canv
        toolbar=self.toolbar
        varbox=self.varbox
        plottpbox=self.plottpbox
        varList=self.varList
        gridNumEnt=self.gridNumEnt

        for var_index in range(len(varList)):
            if varbox.get() == varList[var_index]:
                break

        fpath=os.path.dirname(fullname)
        fname=self.file_sel.get()
        if fname == "":
            print('file name not provided!')
            return -1

        else:
            fullname=fpath+'/'+fname
            if self.oldfilename != fullname:
                self.oldfilename=fullname
                self.data_inputted=0

                self.getfile.set_filename(fullname)

        #fname=os.path.basename(fullname)
        paramList=[ \
            r'O$^+$ Density (cm$^{-3}$)',r'H$^+$ Density (cm$^{-3}$)',r'H$_e^+$ Density (cm$^{-3}$)', \
            r'O$_2^+$ Density (cm$^{-3}$)',r'N$_2^+$ Density (cm$^{-3}$)', r'NO$^+$ Density (cm$^{-3}$)', \
            r'N$^+$ Density (cm$^{-3}$)', \
            r'V_${O+,r}$ (m/s)',r'V$_{O+,\theta}$ (m/s)',r'V$_{O+,\phi}$ (m/s)',
            r'V_${H+,r}$ (m/s)',r'V$_{H+,\theta}$ (m/s)',r'V$_{H+,\phi}$ (m/s)',
            r'V_${He+,r}$ (m/s)',r'V$_{He+,\theta}$ (m/s)',r'V$_{He+,\phi}$ (m/s)',
            r'T$_{O+}$ (K)', r'T$_{H+}$ (K)', r'T$_{He+}$ (K)', r'T$_e$ (K)', \
            r'O Density (cm$^{-3}$)',r'H Density (cm$^{-3}$)',r'H$_e$ Density (cm$^{-3}$)', \
            r'O$_2$ Density (cm$^{-3}$)',r'N$_2$ Density (cm$^{-3}$)',r'NO Density (cm$^{-3}$)', \
            r'N Density (cm$^{-3}$)', \
            r'V$_{n,r}$ (m/s)',r'V$_{n,\theta}$ (m/s)',r'V$_{n,\phi}$ (m/s)',r'T$_n$ (K)', \
            r'$\delta$B$_r$ (nT)',r'$\delta$B$_\theta$ (nT)',r'$\delta$B$_\phi$ (nT)', \
            r'E$_r$ (mV/m)',r'E$_\theta$ (mV/m)',r'E$_\phi$ (mV/m)', \
            r'V$_{e,r}$',r'V$_{e,\theta}$',r'V$_{e,\phi}$', \
            r'Electron Density (cm$^{-3}$)', r'Total Neutral Density (cm$^{-3}$)', \
            r'All Ion Density (cm$^{-3}$)',r'All Neutral Density (cm$^{-3}$)', \
            r'$\delta$B Components (nT)', r'Divergence B (nT/km)', \
            r'V$_{f,r}$ (km/s)',r'V$_{f,\theta}$ (km/s)',r'V$_{f,\phi}$ (km/s)', \
            r'C$_s$ (km/s)','Max Time Step (s)', \
            r'Electron Thermal Conductivity (Joule m$^{-1}$ s$^{-1}$ K$^{-1}$)' ,\
            r'Ion Thermal Conductivity (Joule m$^{-1}$ s$^{-1}$ K$^{-1}$)', \
            r'Neutral Thermal Conductivity (Joule m$^{-1}$ s$^{-1}$ K$^{-1}$)', \
            r'Electron Collision Frequency (s$^{-1}$)', \
            r'O$^+$ Collision Frequency (s$^{-1}$)', \
            r'H$^+_{2}$ Collision Frequency (s$^{-1}$)', \
            r'He$^+_{2}$ Collision Frequency (s$^{-1}$)', \
            r'O$^+_{2}$ Collision Frequency (s$^{-1}$)', \
            r'N$^+_{2}$ Collision Frequency (s$^{-1}$)', \
            r'NO$^+_{2}$ Collision Frequency (s$^{-1}$)', \
            r'N$^+_{2}$ Collision Frequency (s$^{-1}$)']

####### Display data file name in an entry field
        wd=self.rwd
        #tk.Label(self.rwd, text="Data File: "+fname, bg='white').place(x=580,y=45)
        #self.file_sel.delete(0, 'end')
        #self.file_sel.insert(0, fname)
        
        ptype=plottpbox.get()

        #params = read_data.input_param(fpath)

        if (self.data_inputted == 0):
            print(fullname)
            self.iditm_arr, self.params, self.B0 = read_data.input_data(fullname)
            if len(self.iditm_arr) == 1:
                return -1

            a1=self.params[0]
            a2=self.params[1]
            a3=self.params[2]
            self.alt, self.lat, self.lon, self.cr, self.rdth, self.cos_rsinth, \
                self.rsinth_dph = read_data.input_grids(fpath, a1, a2, a3)

            self.data_inputted=1

        a1=self.params[0]
        a2=self.params[1]
        a3=self.params[2]
        sl=self.params[6]
        sm=self.params[7]
        r0=1.0e6

        date_str=str(self.params[3])+'-'+str(self.params[4])+'-'+str(self.params[5])
        UT_hr=int(self.params[8]/3600.0)
        UT_mi=int((self.params[8]-UT_hr*3600.0)/60.0)
        UT_se=self.params[8]-UT_hr*3600.0-UT_mi*60.0
        UT_str='UT='+str(UT_hr)+':'+str(UT_mi)+":"+"{:5.3f}".format(UT_se)

        dayno=day_number.dayofyear()
        days = dayno.dayOfYear(date_str)

        print("\nVariable index: ", var_index)

################## make altitude profile plot #########################################
        if ptype == 'Altitude_Profile':
            k=int(gridNumEnt[2].get())
            j=int(gridNumEnt[1].get())
            i0=int(gridNumEnt[0].get())
            i1=int(gridNumEnt[3].get())

            #print(gridNumEnt[0].get(),gridNumEnt[1].get(),gridNumEnt[2].get(),gridNumEnt[3].get(), \
            #    gridNumEnt[4].get(),gridNumEnt[5].get(),j,a2)
            if (i0 >= i1):
                print("Start i0 must be less than i1")
                return -1
            if (i0 >= a1 or i1 > a1):
                print("Start i0 must be < a1 and end i1 must be <= a1")
                return -1
            if (k >= a3):
                print("Start Index k must be < a3")
                return -1
            if (j >= a2):
                print("Start index j must be < a2")
                return -1

            MLT_str = get_MLT_str(self.params[3], days, self.params[8], self.lon[k])

            add_text=date_str+' '+UT_str+' MLat=' +'{:.2f}'.format(-self.lat[j])+' MLT=' +MLT_str

            varib=np.zeros(i1-i0)
            if var_index==42 or var_index==43:
                varib=np.zeros((sl, i1-i0))
            elif var_index==44:
                varib=np.zeros((3, i1-i0))

            for i in range(i0, i1):
                if var_index < 37:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 7:
                    #ion density in cm^-3
                    varib[i-i0]=self.iditm_arr[k][j][i][var_index]
 
                elif var_index >=7 and var_index <=15:
                    #ion velocity in m/s
                    varib[i-i0]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=16 and var_index <=19:
                    #ion and electron temperature in K
                    varib[i-i0]= self.iditm_arr[k][j][i][var_index]

                elif var_index >= 20 and var_index <=26:
                    #neutral density in cm^-3
                    varib[i-i0]= self.iditm_arr[k][j][i][var_index]

                elif var_index >=27 and var_index <=29:
                    #neutral velocity in m/s
                    varib[i-i0]=self.iditm_arr[k][j][i][var_index]

                elif var_index == 30:
                    #neutral temperature in K
                    varib[i-i0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=31 and var_index <=33:
                    #perturbation magnetic field in nT
                    varib[i-i0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=34 and var_index <=36: #electric field in mV/m
                    varib[i-i0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=37 and var_index <=39: #electron velocity in m/s
                    evel=ele_velocity.evelocity(self.iditm_arr, i, j, k)

                    for s in range(3):
                        if (math.isnan(evel[s])):
                            print('(i,j,k,l)=(',i,j,k,s,'), evel-component Nan')
                            return -1
                        elif(math.isinf(evel[s])):
                            print('(i,j,k,l)=(',i,j,k,s,'), evel-component Inf')
                            return -1

                    varib[i-i0]=evel[var_index-37]  #evel in m/s

                elif var_index == 40: #electron density in cm^-3
                    varib[i-i0]=0.0
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1
                        else:
                            varib[i-i0]=varib[i-i0]+self.iditm_arr[k][j][i][s]

                elif var_index == 41: #total neutral density in cm^-3
                    varib[i-i0]=0.0
                    for s in range(sm):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1
                        else:
                            varib[i-i0]=varib[i-i0]+self.iditm_arr[k][j][i][s20]

                elif var_index == 42: # every ion density
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1

                        varib[s][i-i0]=self.iditm_arr[k][j][i][s]

                elif var_index == 43: # every neutral density
                    for s in range(sm):
                        if (math.isnan(self.iditm_arr[k][j][i][20+s])):
                            print('(i,j,k,s)=(',i,j,k,20+s,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][20+s])):
                            print('(i,j,k,s)=(',i,j,k,20+s,'), Neutral Density Inf')
                            return -1

                        varib[s][i-i0]=self.iditm_arr[k][j][i][20+s]

                elif var_index == 44: # all 3 perturbation mfd components in nT
                    for s in range(3):
                        if (math.isnan(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s+31,'), B-field Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s+31,'), B-field Inf')
                            return -1

                        varib[s][i-i0]=self.iditm_arr[k][j][i][s+31]

                elif var_index == 45: #div-B in T/m
                    for s in range(3):
                        if (math.isnan(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s+31,'), B-component Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s+31,'), B-component Inf')
                            return -1

                    varib[i-i0]=divB.divB(self.iditm_arr, i, j, k, self.alt, self.cr, \
                        self.rdth, self.cos_rsinth, self.rsinth_dph, a1, a2, a3)

                elif var_index >= 46 and var_index <= 48: #fast mode speed in km/s
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                    if var_index == 46:
                        varib[i-i0] = vf_r 
                    elif var_index == 47:
                        varib[i-i0] = vf_t
                    elif var_index == 48:
                        varib[i-i0] = vf_f

                    varib[i-i0]=varib[i-i0]*1.0e-3 #fast mode speed in km/s

                elif var_index == 49: #sound speed in m/s 
                    _, _, _, varib[i-i0] = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                elif var_index == 50: #maximum time step allowed
                    r0=1.0e6
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr,self.B0,i,j,k,sl)

                    if i==0:
                        dr=self.alt[i+1]-self.alt[0]
                    else:
                        dr=self.alt[i]-self.alt[i-1]
   
                    varib[i-i0]=min(dr/vf_r, self.rdth[i]/vf_t, self.rsinth_dph[j][i]/vf_f)*r0

                elif var_index == 51: # electron thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[i-i0]=lamda[0] #in J /(m s K)

                elif var_index == 52: # ion thermal conductivities
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[i-i0]=0.0
                    for s in range(sl):
                        varib[i-i0]=varib[i-i0]+lamda[s+1] #ion thermal conductivity in J /(m s K)

                elif var_index == 53: # neutral thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[i-i0]=lamda[8] #in J /(m s K)

                elif var_index == 54: # electron collision frequencies
                    nuet, _, _=parameters.paramaters(self.iditm_arr, i, j, k)

                    for s in range(8):
                        varib[i-i0]=varib[i-i0]+nuet[s]

                elif var_index >= 55 and var_index<=167: # ion collision frequencies
                    _, nust, _=parameters.paramaters(self.iditm_arr, i, j, k)

                    for s in range(7):
                        varib[i-i0]=varib[i-i0]+nust[(var_index-47)*11+s]

            if varib.ndim ==1:
                print("min value:", min(varib), "; max value: ",max(varib))
                print("Index of min value:", np.argmin(varib), "; Index of max value: ", \
                    np.argmax(varib))

            self.canv, self.fig, self.toolbar = \
            alt_profile.plot_altitude_profile(wd, canv, fig, toolbar, \
                self.alt, gridNumEnt, varib, ptype, var_index, add_text, fname, paramList)

############# completed altitude profile plot #########################################

################## make latitudinal profile plot #########################################
        if ptype == 'Latitude_Profile':
            k=int(gridNumEnt[2].get())
            i=int(gridNumEnt[0].get())
            j0=int(gridNumEnt[1].get()) #colatitude start index
            j1=int(gridNumEnt[4].get()) #colatitude end index

            if (j0 >= j1):
                print("start j0 must be less than end j1")
                return -1
            if (j0 >= a2 or j1 > a2):
                print("Start j0 must be < a2 and end j1 must be <= a2")
                return -1
            if (k >= a3):
                print("Start Index k must be < a3")
                return -1
            if (i >= a1):
                print("Start index i must be <= a1")
                return -1

            MLT_str = get_MLT_str(self.params[3], days, self.params[8], self.lon[k])

            add_text=date_str+' '+UT_str+' Height='+'{:.2f}'.format(self.alt[i])+' MLT='+MLT_str

            varib=np.zeros(j1-j0)
            if var_index==42 or var_index==43:
                varib=np.zeros((sl, j1-j0))
            elif var_index==44:
                varib=np.zeros((3, j1-j0))

            for j in range(j0,j1):
                if var_index < 37:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 7:
                    #ion density in cm^-3
                    varib[j1-1-j]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=7 and var_index <=15:
                    #ion velocity in m/s
                    varib[j1-1-j]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=16 and var_index <=19:
                    #ion and electron temperature in K
                    varib[j1-1-j]= self.iditm_arr[k][j][i][var_index]

                elif var_index >= 20 and var_index <=26:
                    #neutral density in cm^-3
                    varib[j1-1-j]= self.iditm_arr[k][j][i][var_index]

                elif var_index >=27 and var_index <=29:
                    #neutral velocity in m/s
                    varib[j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index == 30:
                    #neutral temperature in K
                    varib[j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=31 and var_index <=33:
                    #perturbation magnetic field in nT
                    varib[j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=34 and var_index <=36: #electron veclocity in m/s
                    varib[j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=37 and var_index <=39: #electron veclocity in m/s
                    evel=ele_velocity.evelocity(self.iditm_arr, i, j, k)

                    for s in range(3):
                        if (math.isnan(evel[s])):
                            print('(i,j,k,l)=(',i,j,k,s,'), evel-component Nan')
                            return -1
                        elif(math.isinf(evel[s])):
                            print('(i,j,k,l)=(',i,j,k,s,'), evel-component Inf')
                            return -1

                    varib[j1-1-j]=evel[var_index-37]  #evel in m/s

                elif var_index == 40: #electron density in cm^-3
                    varib[j1-1-j]=0.0
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1
                        else:
                            varib[j1-1-j]=varib[j1-1-j]+self.iditm_arr[k][j][i][s]

                elif var_index == 41: #total neutral density in cm^-3
                    varib[j1-1-j]=0.0
                    for s in range(sm):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1
                        else:
                            varib[j1-1-j]=varib[j1-1-j]+self.iditm_arr[k][j][i][s20]

                    #varib[j1-1-j]=varib[j1-1-j]*1.0e-6

                elif var_index == 42: # every ion density
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1

                        varib[s][j1-1-j]=self.iditm_arr[k][j][i][s]

                elif var_index == 43: # every neutral density
                    for s in range(sm):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1

                        varib[s][j1-1-j]=self.iditm_arr[k][j][i][s20]

                elif var_index == 44: # every perturbation mfd components in nT
                    for s in range(3):
                        if (math.isnan(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-field Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-field Inf')
                            return -1

                        varib[s][j1-1-j]=self.iditm_arr[k][j][i][s+31]

                elif var_index == 45: #div-B in T/m
                    for s in range(3):
                        if (math.isnan(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-component Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-component Inf')
                            return -1

                    varib[j1-1-j]=divB.divB(self.iditm_arr, i, j, k, self.alt, self.cr, \
                        self.rdth, self.cos_rsinth, self.rsinth_dph, a1, a2, a3)

                elif var_index >= 46 and var_index <= 48: #fast mode speed in km/s
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr,self.B0,i,j,k,sl)

                    if var_index == 46:
                        varib[j1-1-j] = vf_r 
                    elif var_index == 47:
                        varib[j1-1-j] = vf_t
                    elif var_index == 48:
                        varib[j1-1-j] = vf_f

                    varib[j1-1-j]=varib[j1-1-j]*1.0e-3 #fast mode speed in km/s

                elif var_index == 49: #sound speed in m/s 
                    _, _, _, varib[j1-1-j] = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                elif var_index == 50: #maximum time step allowed
                    r0=1.0e6
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                    if i==0:
                        dr=self.alt[i+1]-self.alt[0]
                    else:
                        dr=self.alt[i]-self.alt[i-1]
                    
                    varib[j1-1-j]=min(dr/vf_r, self.rdth[i]/vf_t, self.rsinth_dph[j][i]/vf_f)*r0

                elif var_index == 41: # electron thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[j1-1-j]=lamda[0] #in J /(m s K)

                elif var_index == 52: # ion thermal conductivities
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[j1-1-j]=0.0
                    for s in range(sl):
                        varib[j1-1-j]=varib[j1-1-j]+lamda[s+1] #in J /(m s K)

                elif var_index == 53: # neutral thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[j1-1-j]=lamda[8] #in J /(m s K)

            if varib.ndim ==1:
              print("min value:", min(varib), "; max value: ",max(varib))
              print("Index of min value:", np.argmin(varib), "; Index of max value: ", \
                  np.argmax(varib))

            self.canv, self.fig, self.toolbar = \
            latitude_profile.latitude_profile(wd, canv, fig, toolbar, \
                self.lat, gridNumEnt, varib, ptype, var_index, add_text, fname, paramList)
############# completed latitudinal profile plot #########################################

################## make longitude profile plot #########################################
        if ptype == 'Longitude_Profile':
            k0=int(gridNumEnt[2].get())
            k1=int(gridNumEnt[5].get())
            i=int(gridNumEnt[0].get())
            j=int(gridNumEnt[1].get())

            if (k0 >= k1):
                print("start k0 must be less than end k1")
                return -1
            if (k0 >= a3 or k1 > a3):
                print("Start k0 must be < a3 and end k1 must be <= a3")
                return -1
            if (i >= a1):
                print("Start Index i must be less than a1")
                return -1
            if (j >= a2):
                print("Start index j must be less than a2")
                return -1

            add_text=date_str+' '+UT_str +' Height='+'{:.2f}'.format(self.alt[i]) \
                +' Latitude='+'{:.2f}'.format(-self.lat[j])

            varib=np.zeros(k1-k0)
            if var_index==42 or var_index==43:
                varib=np.zeros((sl, k1-k0))
            elif var_index==44:
                varib=np.zeros((3, k1-k0))

            for k in range(k0, k1):
                if var_index < 37:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 7:
                    #ion density in cm^-3
                    varib[k-k0]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=7 and var_index <=15:
                    #ion velocity in m/s
                    varib[k-k0]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=16 and var_index <=19:
                    #ion and electron temperature in K
                    varib[k-k0]= self.iditm_arr[k][j][i][var_index]

                elif var_index >= 20 and var_index <=26:
                    #neutral density in cm^-3
                    varib[k-k0]= self.iditm_arr[k][j][i][var_index]

                elif var_index >=27 and var_index <=29:
                    #neutral velocity in m/s
                    varib[k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index == 30:
                    #neutral temperature in K
                    varib[k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=31 and var_index <=33:
                    #perturbation magnetic field in nT
                    varib[k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=34 and var_index <=36:
                    #e-field in mV/m
                    varib[k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=37 and var_index <=29: #electron veclocity in m/s
                    evel=ele_velocity.evelocity(self.iditm_arr, i, j, k)

                    for s in range(3):
                        if (math.isnan(evel[s])):
                            print('(i,j,k,l)=(',i,j,k,s,'), evel-component Nan')
                            return -1
                        elif(math.isinf(evel[s])):
                            print('(i,j,k,l)=(',i,j,k,s,'), evel-component Inf')
                            return -1

                    varib[k-k0]=evel[var_index-37]  #evel in m/s

                elif var_index == 40: #electron density in cm^-3
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1
                        else:
                            varib[k-k0]=varib[k-k0]+self.iditm_arr[k][j][i][s]

                    #electron density in cm^-3
                    #varib[k-k0]=varib[k-k0]*1.0e-6

                elif var_index == 41: #total neutral density in cm^-3
                    for s in range(sm):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1
                        else:
                            varib[k-k0]=varib[k-k0]+self.iditm_arr[k][j][i][s20]

                    #total neutral density in cm^-3
                    #varib[k-k0]=varib[k-k0]*1.0e-6

                elif var_index == 42: # every ion density
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1

                        varib[s][k-k0]=self.iditm_arr[k][j][i][s]

                elif var_index == 43: # every neutral density
                    for s in range(7):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1

                        varib[s][k-k0]=self.iditm_arr[k][j][i][s20]

                elif var_index == 44: # every perturbation mfd components in nT
                    for s in range(3):
                        if (math.isnan(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-field Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-field Inf')
                            return -1

                        varib[s][k-k0]=self.iditm_arr[k][j][i][s+31]

                elif var_index == 45: #div-B in T/m
                    for s in range(3):
                        if (math.isnan(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-field Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s+31])):
                            print('(i,j,k,s)=(',i,j,k,s,'), B-field Inf')
                            return -1

                    varib[k-k0]=divB.divB(self.iditm_arr, i, j, k, self.alt, self.cr, \
                        self.rdth, self.cos_rsinth, self.rsinth_dph, a1, a2, a3)

                elif var_index >= 46 and var_index <= 48: #fast mode speed in km/s
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                    if var_index == 46:
                        varib[k-k0] = vf_r 
                    elif var_index == 47:
                        varib[k-k0] = vf_t
                    elif var_index == 48:
                        varib[k-k0] = vf_f

                    varib[k-k0]=varib[k-k0]*1.0e-3 #fast mode speed in km/s

                elif var_index == 49: #sound speed in m/s 
                    _, _, _, varib[k-k0] = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                elif var_index == 50: #maximum time step allowed
                    r0=1.0e6
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                    if i==0:
                        dr=self.alt[i+1]-self.alt[0]
                    else:
                        dr=self.alt[i]-self.alt[i-1]
                    
                    varib[k-k0]=min(dr/vf_r, self.rdth[i]/vf_t, self.rsinth_dph[j][i]/vf_f)*r0

                elif var_index == 51: # electron thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[k-k0]=lamda[0] #electron thermal conductivity in J/(m s K)

                elif var_index == 52: # ion thermal conductivities
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[k-k0]=0.0
                    for s in range(sl):
                        varib[k-k0]=varib[k-k0]+lamda[s+1] #ion thermal conductivity in J/(m s K)

                elif var_index == 53: # neutral thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[k-k0]=lamda[8] #in J /(m s K)

            if varib.ndim ==1:
              print("min value:", min(varib), "; max value: ",max(varib))
              print("Index of min value:", np.argmin(varib), "; Index of max value: ", \
                  np.argmax(varib))

            self.canv, self.fig, self.toolbar = \
                long_profile.longitude_profile(wd, canv, fig, toolbar, \
                self.lon, gridNumEnt, varib, ptype, var_index, add_text, fname, paramList)
############# completed longitude profile plot #########################################

############# make Latitude - Height contour plot #####################################
        elif (ptype == 'Lat-Height_Contour'):
            if (var_index >= 42 and var_index <= 44):
                print("Cann't make contour plot for multiple variables!")
                return -1

            k=int(gridNumEnt[2].get())
            j0=int(gridNumEnt[1].get())
            j1=int(gridNumEnt[4].get())
            i0=int(gridNumEnt[0].get())
            i1=int(gridNumEnt[3].get())

            if(j0 >= j1):
                print("start j0 must be less than end j1")
                return -1
            if (i0 >= i1):
                print("start i0 must be less than end i1")
                return -1
            if (j0 >= a2 or j1 > a2):
                print("Start j0 must be < a2 and end j1 must be <= a2")
                return -1
            if (i0 >= a1 or i1 > a1):
                print("Start i0 must be < a1 and end i1 must be <= a1")
                return -1
            if (k >= a3):
                print("Start index k must be less than a3")
                return -1

            MLT_str = get_MLT_str(self.params[3], days, self.params[8], self.lon[k])

            add_text=date_str+' '+UT_str+" MLT="+MLT_str

            varib=np.zeros((i1-i0, j1-j0), dtype=float)

            for j in range(j0, j1):
              for i in range(i0, i1):
                if var_index < 37:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 7:
                    #ion density in cm^-3
                    varib[i-i0][j1-1-j]=np.log10(self.iditm_arr[k][j][i][var_index])

                elif var_index >=7 and var_index <=15:
                    #ion velocity in m/s
                    varib[i-i0][j1-1-j]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=16 and var_index <=19:
                    #ion and electron temperature in K
                    varib[i-i0][j1-1-j]=self.iditm_arr[k][j][i][var_index]

                elif var_index >= 20 and var_index <=26:
                    #neutral density in cm^-3
                    varib[i-i0][j1-1-j]=np.log10(self.iditm_arr[k][j][i][var_index])

                elif var_index >=27 and var_index <=29:
                    #neutral velocity in m/s
                    varib[i-i0][j1-1-j]= self.iditm_arr[k][j][i][var_index]

                elif var_index == 30:
                    #neutral temperature in K
                    varib[i-i0][j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=31 and var_index <=33:
                    #perturbation magnetic field in nT
                    varib[i-i0][j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=34 and var_index <=36:
                     #e-field in mV/m
                    varib[i-i0][j1-1-j] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=37 and var_index <=29: #electron veclocity in m/s
                    evel=ele_velocity.evelocity(self.iditm_arr, i, j, k)
                    varib[i-i0][j1-1-j]=evel[var_index-37] #evel in m/s

                elif var_index == 40: #electron density in cm^-3
                    varib[i-i0][j1-1-j]=0.0
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1

                        varib[i-i0][j1-1-j]=varib[i-i0][j1-1-j]+self.iditm_arr[k][j][i][s]

                    #electron density in cm^-3
                    varib[i-i0][j1-1-j]=np.log10(varib[i-i0][j1-1-j])

                elif var_index == 41: #total neutral density in cm^-3
                    varib[i-i0][j1-1-j]=0.0
                    for s in range(sm):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1

                        varib[i-i0][j1-1-j]=varib[i-i0][j1-1-j]+self.iditm_arr[k][j][i][s20]

                    #neutral density in cm^-3
                    varib[i-i0][j1-1-j]=np.log10(varib[i-i0][j1-1-j])

                elif var_index == 45: #div-B in T/m (x 10^12)
                    varib[i-i0][j1-1-j]=divB.divB(self.iditm_arr, i, j, k, \
                        self.alt, self.cr, self.rdth, self.cos_rsinth, \
                        self.rsinth_dph, a1, a2, a3)*1.0e12

                elif var_index >= 46 and var_index <= 48: #fast mode speed in km/s
                    if var_index == 46:
                        varib[i-i0][j1-1-j], _, _, _ = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)
                    elif var_index == 47:
                        _, varib[i-i0][j1-1-j], _, _ = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)
                    elif var_index == 48:
                        _, _, varib[i-i0][j1-1-j], _ = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)

                    #fast mode speed in km/s
                    varib[i-i0][j1-1-j]=varib[i-i0][j1-1-j]*1.0e-3

                elif var_index == 49: #sound speed in m/s 
                    _, _, _, varib[i-i0][j1-1-j] = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)

                elif var_index == 50: #maximum time step allowed
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                    if i==0:
                        dr=self.alt[i+1]-self.alt[0]
                    else:
                        dr=self.alt[i]-self.alt[i-1]
                    
                    varib[i-i0][j1-1-j]=min(dr/vf_r, self.rdth[i]/vf_t, \
                        self.rsinth_dph[j][i]/vf_f)*r0
                    varib[i-i0][j1-1-j]=np.log10(varib[i-i0][j1-1-j])

                elif var_index == 51: # electron thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[i-i0][j1-1-j]=math.log10(lamda[0]) #in J /(m s K) logscale

                elif var_index == 52: # ion thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[i-i0][j1-1-j]=0.0
                    for s in range(sl):
                        varib[i-i0][j1-1-j] = varib[i-i0][j1-1-j]+lamda[s+1]

                    varib[i-i0][j1-1-j]=math.log10(varib[i-i0][j1-1-j]) #in J /(m s K) logscale

                elif var_index == 53: # neutral thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[i-i0][j1-1-j]=math.log10(lamda[8]) #in J /(m s K) logscale

            minval=np.min(varib)
            maxval=np.max(varib)
            print("min value:", minval, "; max value: ",maxval)
            print("Index of min value:", np.argwhere(varib == minval), \
                "; Index of max value: ", np.argwhere(varib == maxval))

            self.canv, self.fig, self.toolbar = \
                lati_height_cont.latitude_altitude_contour(wd, canv, fig, toolbar, \
                    self.lat, self.alt, gridNumEnt, varib, ptype, var_index, \
                    add_text, fname, paramList)
############# completed Latitude - Height contour plot #################################

############# make Longitude - Latitude contour plot #####################################
        elif (ptype == 'Long-Lat_Contour'):
            if (var_index >= 42 and var_index <= 44):
                print("Cann't make contour plot for multiple variables!")
                return -1

            k0=int(gridNumEnt[2].get())
            k1=int(gridNumEnt[5].get())
            j0=int(gridNumEnt[1].get())
            j1=int(gridNumEnt[4].get())
            i =int(gridNumEnt[0].get())

            if(j0 >= j1):
                print("start j0 must be less than end j1")
                return -1
            if (k0 >= k1):
                print("start k0 must be less than end k1")
                return -1
            if (j0 >= a2 or j1 > a2):
                print("Start j0 must be < a2 and end j1 must be <= a2")
                return -1
            if (k0 >= a3 or k1 > a3):
                print("Start k0 must be < a3 and end k1 must be <= a3")
                return -1
            if (i >= a1):
                print("Start index i must be less than a1")
                return -1

            add_text=date_str+' '+UT_str+' Altitude='+'{:.2f}'.format(self.alt[i])+' km'

            varib=np.zeros((j1-j0, k1-k0), dtype=float)

            for k in range(k0, k1):
              for j in range(j0, j1):
                if var_index < 37:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 7:
                    #ion density in cm^-3
                    varib[j1-1-j][k-k0]=math.log10(self.iditm_arr[k][j][i][var_index])

                elif var_index >=7 and var_index <=15:
                    #ion velocity in m/s
                    varib[j1-1-j][k-k0]=self.iditm_arr[k][j][i][var_index]

                elif var_index >=16 and var_index <=19:
                    #ion and electron temperature in K
                    varib[j1-1-j][k-k0]=self.iditm_arr[k][j][i][var_index]

                elif var_index >= 20 and var_index <=26:
                    #neutral density in cm^-3
                    varib[j1-1-j][k-k0]=np.log10(self.iditm_arr[k][j][i][var_index])

                elif var_index >=27 and var_index <=29:
                    #neutral velocity in m/s
                    varib[j1-1-j][k-k0]= self.iditm_arr[k][j][i][var_index]

                elif var_index == 30:
                    #neutral temperature in K
                    varib[j1-1-j][k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=31 and var_index <=33:
                    #perturbation magnetic field in nT
                    varib[j1-1-j][k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=34 and var_index <=36:
                    #electric field in mV/m
                    varib[j1-1-j][k-k0] = self.iditm_arr[k][j][i][var_index]

                elif var_index >=37 and var_index <=39: #electron veclocity in m/s
                    evel=ele_velocity.evelocity(self.iditm_arr, i, j, k)
                    varib[j1-1-j][k-k0]=evel[var_index-37] #evel in m/s

                elif var_index == 40: #electron density in cm^-3
                    varib[j1-1-j][k-k0]=0.0
                    for s in range(sl):
                        if (math.isnan(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s])):
                            print('(i,j,k,s)=(',i,j,k,s,'), Ion Density Inf')
                            return -1

                        varib[j1-1-j][k-k0]=varib[j1-1-j][k-k0]+self.iditm_arr[k][j][i][s]

                    #electron density in cm^-3
                    varib[j1-1-j][k-k0]=np.log10(varib[j1-1-j][k-k0])

                elif var_index == 41: #total neutral density in cm^-3
                    varib[j1-1-j][k-k0]=0.0
                    for s in range(sm):
                        s20=20+s
                        if (math.isnan(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s20])):
                            print('(i,j,k,s)=(',i,j,k,s20,'), Neutral Density Inf')
                            return -1

                        varib[j1-1-j][k-k0]=varib[j1-1-j][k-k0]+self.iditm_arr[k][j][i][s20]

                    #neutral density in cm^-3
                    varib[j1-1-j][k-k0]=np.log10(varib[j1-1-j][k-k0])

                elif var_index == 45: #div-B in T/m (x 10^12)
                    varib[j1-1-j][k-k0]=divB.divB(self.iditm_arr, i, j, k, \
                        self.alt, self.cr, self.rdth, self.cos_rsinth, \
                        self.rsinth_dph, a1, a2, a3)*1.0e12

                elif var_index >= 46 and var_index <= 48: #fast mode speed in km/s
                    if var_index == 46:
                        varib[j1-1-j][k-k0], _, _, _ = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)
                    elif var_index == 47:
                        _, varib[j1-1-j][k-k0], _, _ = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)
                    elif var_index == 48:
                        _, _, varib[j1-1-j][k-k0], _ = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)

                    #fast mode speed in km/s
                    varib[j1-1-j][k-k0]=varib[j1-1-j][k-k0]*1.0e-3

                elif var_index == 49: #sound speed in m/s 
                    _, _, _, varib[j1-1-j][k-k0] = fast_wave_speed.fast_speed \
                            (self.iditm_arr,self.B0, i, j, k, sl)

                elif var_index == 50: #maximum time step allowed
                    vf_r, vf_t, vf_f, _ = fast_wave_speed.fast_speed(self.iditm_arr, \
                        self.B0,i,j,k,sl)

                    if i==0:
                        dr=self.alt[i+1]-self.alt[0]
                    else:
                        dr=self.alt[i]-self.alt[i-1]
                    
                    varib[j1-1-j][k-k0]=min(dr/vf_r, self.rdth[i]/vf_t, \
                        self.rsinth_dph[j][i]/vf_f)*r0
                    varib[j1-1-j][k-k0]=np.log10(varib[i-i0][j1-1-j])

                elif var_index == 51: # electron thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[j1-1-j][k-k0]=math.log10(lamda[0]) #in J /(m s K) logscale

                elif var_index == 52: # ion thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[j1-1-j][k-k0]=0.0
                    for s in range(sl):
                        varib[j1-1-j][k-k0] = varib[j1-1-j][k-k0]+lamda[s+1]

                    varib[j1-1-j][k-k0]=math.log10(varib[j1-1-j][k-k0]) #in J /(m s K) logscale

                elif var_index == 53: # neutral thermal conductivity
                    _, _, lamda=parameters.paramaters(self.iditm_arr, i, j, k)

                    varib[j1-1-j][k-k0]=math.log10(lamda[8]) #in J /(m s K) logscale

            minval=np.min(varib)
            maxval=np.max(varib)
            print("min value:", minval, "; max value: ",maxval)
            print("Index of min value:", np.argwhere(varib == minval), \
                "; Index of max value: ", np.argwhere(varib == maxval))

            self.canv, self.fig, self.toolbar = \
                long_lat_cont.lon_lat_contour(wd, canv, fig, toolbar, \
                    self.lon, self.lat, gridNumEnt, varib, ptype, var_index, \
                    add_text, fname)
############# completed Latitude - Height contour plot #################################

        return 0
