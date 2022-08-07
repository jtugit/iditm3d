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
        paramList=[r'$\nu_{e,O+}$ (s$^{-1}$)',r'$\nu_{e,O2+}$ (s$^{-1}$)', \
            r'$\nu_{e,N2+}$ (s$^{-1}$)',r'$\nu_{e,H+}$ (s$^{-1}$)', \
            r'$\nu_{e,He+}$ (s$^{-1}$)',r'$\nu_{e,NO+}$ (s$^{-1}$)', \
            r'$\nu_{O+,e+}$ (s$^{-1}$)',r'$\nu_{O+,O2+}$ (s$^{-1}$)', \
            r'$\nu_{O+,N2+}$ (s$^{-1}$)',r'$\nu_{O+,H+}$ (s$^{-1}$)', \
            r'$\nu_{O+,He+}$ (s$^{-1}$)',r'$\nu_{O+,NO+}$ (s$^{-1}$)', \
            r'$\nu_{O+,O}$ (s$^{-1}$)',r'$\nu_{O+,O2}$ (s$^{-1}$)', \
            r'$\nu_{O+,N2}$ (s$^{-1}$)',r'$\nu_{O+,H}$ (s$^{-1}$)', \
            r'$\nu_{O+,He}$ (s$^{-1}$)', \
            r'$\nu_{O2+,e+}$ (s$^{-1}$)',r'$\nu_{O2+,O+}$ (s$^{-1}$)', \
            r'$\nu_{O2+,N2+}$ (s$^{-1}$)',r'$\nu_{O2+,H+}$ (s$^{-1}$)', \
            r'$\nu_{O2+,He+}$ (s$^{-1}$)',r'$\nu_{O2+,NO+}$ (s$^{-1}$)', \
            r'$\nu_{O2+,O}$ (s$^{-1}$)',r'$\nu_{O2+,O2}$ (s$^{-1}$)', \
            r'$\nu_{O2+,N2}$ (s$^{-1}$)',r'$\nu_{O2+,H}$ (s$^{-1}$)', \
            r'$\nu_{O2+,He}$ (s$^{-1}$)', \
            r'$\nu_{N2+,e+}$ (s$^{-1}$)',r'$\nu_{N2+,O+}$ (s$^{-1}$)', \
            r'$\nu_{N2+,O2+}$ (s$^{-1}$)',r'$\nu_{N2+,H+}$ (s$^{-1}$)', \
            r'$\nu_{N2+,He+}$ (s$^{-1}$)',r'$\nu_{N2+,NO+}$ (s$^{-1}$)', \
            r'$\nu_{N2+,O}$ (s$^{-1}$)',r'$\nu_{N2+,O2}$ (s$^{-1}$)', \
            r'$\nu_{N2+,N2}$ (s$^{-1}$)',r'$\nu_{N2+,H}$ (s$^{-1}$)', \
            r'$\nu_{N2+,He}$ (s$^{-1}$)', \
            r'$\nu_{H+,e+}$ (s$^{-1}$)',r'$\nu_{H+,O+}$ (s$^{-1}$)', \
            r'$\nu_{H+,O2+}$ (s$^{-1}$)',r'$\nu_{H+,N2+}$ (s$^{-1}$)', \
            r'$\nu_{H+,He+}$ (s$^{-1}$)',r'$\nu_{H+,NO+}$ (s$^{-1}$)', \
            r'$\nu_{H+,O}$ (s$^{-1}$)',r'$\nu_{H+,O2}$ (s$^{-1}$)', \
            r'$\nu_{H+,N2}$ (s$^{-1}$)',r'$\nu_{H+,H}$ (s$^{-1}$)', \
            r'$\nu_{H+,He}$ (s$^{-1}$)', \
            r'$\nu_{He+,e+}$ (s$^{-1}$)',r'$\nu_{He+,O+}$ (s$^{-1}$)', \
            r'$\nu_{He+,O2+}$ (s$^{-1}$)',r'$\nu_{He+,N2+}$ (s$^{-1}$)', \
            r'$\nu_{He+,H+}$ (s$^{-1}$)',r'$\nu_{He+,NO+}$ (s$^{-1}$)', \
            r'$\nu_{He+,O}$ (s$^{-1}$)',r'$\nu_{He+,O2}$ (s$^{-1}$)', \
            r'$\nu_{He+,N2}$ (s$^{-1}$)',r'$\nu_{He+,H}$ (s$^{-1}$)', \
            r'$\nu_{He+,He}$ (s$^{-1}$)', \
            r'$\nu_{NO+,e+}$ (s$^{-1}$)',r'$\nu_{NO+,O+}$ (s$^{-1}$)', \
            r'$\nu_{NO+,O2+}$ (s$^{-1}$)',r'$\nu_{NO+,N2+}$ (s$^{-1}$)', \
            r'$\nu_{NO+,H+}$ (s$^{-1}$)',r'$\nu_{NO+,He+}$ (s$^{-1}$)', \
            r'$\nu_{NO+,O}$ (s$^{-1}$)',r'$\nu_{NO+,O2}$ (s$^{-1}$)', \
            r'$\nu_{NO+,N2}$ (s$^{-1}$)',r'$\nu_{NO+,H}$ (s$^{-1}$)', \
            r'$\nu_{NO+,He}$ (s$^{-1}$)', \
            r'$\nu_{O,O2}$ (s$^{-1}$)',r'$\nu_{O,N2}$ (s$^{-1}$)', \
            r'$\nu_{O2,N2|$ (s$^{-1}$)',r'$\nu_{H,O}$ (s$^{-1}$)', \
            r'$\nu_{H,O2}$ (s$^{-1}$)',r'$\nu_{H,N2}$ (s$^{-1}$)', \
            r'$\nu_{H,He}$ (s$^{-1}$)',r'$\nu_{He,O}$ (s$^{-1}$)', \
            r'$\nu_{He,O2}$ (s$^{-1}$)',r'$\nu_{He,N2}$ (s$^{-1}$)', \
            r'$\nu_{e,O}$ (s$^{-1}$)',r'$\nu_{e,O2}$ (s$^{-1}$)', \
            r'$\nu_{e,N2}$ (s$^{-1}$)', r'$\nu_{e,H}$ (s$^{-1}$)', \
            r'$\nu_{e,He}$ (s$^{-1}$)', \
            r'Photoelectron Heating Rate (eV cm$^{-3}$ s$^{-1}$)', \
            r'Electron Cooling Rate (eV cm$^{-3}$ s$^{-1}$)', \
            r'Neutral Heating Rate (eV cm$^{-3}$ s$^{-1}$)', \
            r'Neutral Cooling Rate (eV cm$^{-3}$ s$^{-1}$)', \
            r'Electron Thermal Conductivity  (eV cm$^{-3}$ s$^{-1}$)', \
            r'O+ Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'O2+ Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'N2+ Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'H+ Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'He+ Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'NO+ Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'Neutral Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)', \
            r'Electron Collision Frequency  (s$^{-1}$)', \
            r'O+ Collision Frequency  (s$^{-1}$)', \
            r'O2+ Collision Frequency  (s$^{-1}$)', \
            r'N2+ Collision Frequency  (s$^{-1}$)', \
            r'H+ Collision Frequency  (s$^{-1}$)', \
            r'He+ Collision Frequency  (s$^{-1}$)', \
            r'NO+ Collision Frequency  (s$^{-1}$)', \
            r'Ion Thermal Conductivity (eV cm$^{-1}$ s$^{-1}$)']

        print('Variable index:', var_index)

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
        UT_str='/UT='+str(UT_hr)+':'+str(UT_mi)+":"+"{:5.3f}".format(UT_se)

        #height and cof for reducing electron thermal conductivity
        ne_pk=400.0
        ind=0.0 #4.5

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

            LT=self.params[8]/3600.0+self.lon[k]/15.0
            if (LT >= 24.0):
                LT=LT-24.0
            LT_str='{:5.2f}'.format(LT)
            add_text=date_str+" LT="+LT_str+UT_str +' Latitude=' \
                +'{:.2f}'.format(-self.lat[j])+' Longitude=' \
                +'{:.2f}'.format(self.lon[k]) 

            varib=np.zeros(i1-i0)
            for i in range(i0, i1):
                if var_index <99:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 87:
                    #collision frequencies
                    varib[i-i0]=self.iditm_arr[k][j][i][var_index]

                elif (var_index >= 87 and var_index<=90):
                    #heating and cooling rates in eV cm^-3 s^-1
                    varib[i-i0]=self.iditm_arr[k][j][i][var_index]*1.0e13/1.6022

                elif var_index>=91 and var_index<=98:
                    #thermal conductivities in eV cm^-1 s^-1
                    varib[i-i0]= self.iditm_arr[k][j][i][var_index]*1.0e17/1.6022

            if var_index==99:
                varib=np.zeros((7,i1-i0))

                for i in range(i0,i1):
                    for l in range(6):
                        varib[l][i-i0]=self.iditm_arr[k][j][i][l]

                    varib[6][i-i0]=0.0
                    for l in range(82,87):
                        varib[6][i-i0]=varib[6][i-i0]+self.iditm_arr[k][j][i][l]

            elif var_index>=100 and var_index<=105:
                varib=np.zeros((7,i1-i0))

                for i in range(i0,i1):
                    for l in range(6):
                        varib[l][i-i0]=self.iditm_arr[k][j][i][l+6+11*(var_index-100)]

                    for l in range(5):
                        varib[6][i-i0]= varib[6][i-i0] \
                                       +self.iditm_arr[k][j][i][l+12+11*(var_index-100)]

            elif var_index==106:
                varib=np.zeros((6,i1-i0))

                for i in range(i0,i1):
                    for l in range(6):
                        varib[l][i-i0]=self.iditm_arr[k][j][i][l+92]*1.0e17/1.6022

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

            LT=self.params[8]/3600.0+self.lon[k]/15.0
            if (LT >= 24.0):
                LT=LT-24.0
            LT_str='{:5.2f}'.format(LT)
            add_text=date_str+" LT="+LT_str+UT_str +' Height=' \
                +'{:.2f}'.format(self.alt[i])+' Longitude='+'{:.2f}'.format(self.lon[k]) 

            varib=np.zeros(j1-j0)
            for j in range(j0,j1):
                if var_index <99:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 87:
                    #collision frequencies
                    varib[j1-1-j]=self.iditm_arr[k][j][i][var_index]

                elif var_index >= 87 and var_index<=90:
                    #heating and cooling rates in eV cm^-3 s^-1
                    varib[j1-1-j]=self.iditm_arr[k][j][i][var_index]*1.0e13/1.6022

                elif var_index>=91 and var_index<98:
                    #thermal conductivities in eV cm^-1 s^-1
                    varib[j1-1-j]=self.iditm_arr[k][j][i][var_index]*1.0e17/1.6022

            if var_index==99:
                varib=np.zeros((7,j1-j0))

                for i in range(j0,j1):
                    for l in range(6):
                        varib[l][j1-i-j]=self.iditm_arr[k][j][i][l]

                    varib[6][j1-i-j]=0.0
                    for l in range(82,87):
                        varib[6][j1-i-j]=varib[6][j1-i-j]+self.iditm_arr[k][j][i][l]

            elif var_index>=100 and var_index<=105:
                varib=np.zeros((7,j1-j0))

                for i in range(j0,j1):
                    for l in range(6):
                        varib[l][j1-i-j]=self.iditm_arr[k][j][i][l+6+11*(var_index-100)]

                    for l in range(5):
                        varib[6][j1-i-j]= varib[6][j1-i-j] \
                                         +self.iditm_arr[k][j][i][l+12+11*(var_index-100)]

            elif var_index==106:
                varib=np.zeros((6,j1-j0))

                for i in range(j0,j1):
                    for l in range(6):
                        varib[l][j1-i-j]=self.iditm_arr[k][j][i][l+92]*1.0e17/1.6022

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

            add_text=date_str+UT_str +' Height='+'{:.2f}'.format(self.alt[i]) \
                +' Latitude='+'{:.2f}'.format(-self.lat[j]) 

            varib=np.zeros(k1-k0)
            for k in range(k0, k1):
                if var_index <99:
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 87:
                    #collision frequencies
                    varib[k-k0]=self.iditm_arr[k][j][i][var_index]

                elif var_index >= 87 and var_index<=90:
                    #heating and cooling rates in eV cm^-3 s^-1
                    varib[k-k0]=self.iditm_arr[k][j][i][var_index]*1.0e13/1.6022

                elif var_index>=91 and var_index<=98:
                    #thermal conductivities in eV cm^-1 s^-1
                    varib[k-k0]= self.iditm_arr[k][j][i][var_index]*1.0e17/1.6022

            if var_index==99:
                varib=np.zeros((7,k1-k0))

                for i in range(k0,k1):
                    for l in range(6):
                        varib[l][k-k0]=self.iditm_arr[k][j][i][l]

                    varib[6][k-k0]=0.0
                    for l in range(82,87):
                        varib[6][k-k0]=varib[6][k-k0]+self.iditm_arr[k][j][i][l]

            elif var_index>=100 and var_index<=105:
                varib=np.zeros((7,k1-k0))

                for i in range(k0,k1):
                    for l in range(6):
                        varib[l][k-k0]=self.iditm_arr[k][j][i][l+6+11*(var_index-100)]

                    for l in range(5):
                        varib[6][k-k0]= varib[6][k-k0] \
                                       +self.iditm_arr[k][j][i][l+12+11*(var_index-100)]

            elif var_index==106:
                varib=np.zeros((6,k1-k0))

                for i in range(k0,k1):
                    for l in range(6):
                        varib[l][k-k0]=self.iditm_arr[k][j][i][l+92]*1.0e17/1.6022

            self.canv, self.fig, self.toolbar = \
                long_profile.longitude_profile(wd, canv, fig, toolbar, \
                self.lon, gridNumEnt, varib, ptype, var_index, add_text, fname, paramList)
############# completed longitude profile plot #########################################

############# make Latitude - Height contour plot #####################################
        elif (ptype == 'Lat-Height_Contour'):
            if var_index > 98:
                print("cann't make contour plot for multiple variables!")
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

            LT=self.params[8]/3600.0+self.lon[k]/15.0
            if (LT >= 24.0):
                LT=LT-24.0
            LT_str='{:5.2f}'.format(LT)
            add_text=' '+date_str+" LT="+LT_str+UT_str+' Longitude=' \
                     +'{:.2f}'.format(self.lon[k])

            varib=np.zeros((i1-i0, j1-j0), dtype=float)

            for j in range(j0, j1):
                for i in range(i0, i1):
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if var_index < 87:
                    #collision frequencies in logacle 
                    varib[i-i0][j1-1-j]=math.log10(self.iditm_arr[k][j][i][var_index])

                elif var_index >= 87 and var_index <=90:
                    #heating and cooling rates in eV cm^-3 s^-1 in logscale 
                    varib[i-i0][j1-1-j]=math.log10(self.iditm_arr[k][j][i][var_index]
                                                   *1.0e13/1.6022)

                elif var_index>=91 and var_index<=98:
                    #thermal conductivities in eV cm^-1 s^-1 in logscale 
                    varib[k-k0]= self.iditm_arr[k][j][i][var_index]*1.0e17/1.6022
                    varib[i-i0][j1-1-j]=math.log10(self.iditm_arr[k][j][i][var_index] \
                                                   *1.0e17/1.6022)

            self.canv, self.fig, self.toolbar = \
                lati_height_cont.latitude_altitude_contour(wd, canv, fig, toolbar, \
                    self.lat, self.alt, gridNumEnt, varib, ptype, var_index, \
                    add_text, fname)
############# completed Latitude - Height contour plot #################################

############# make Longitude - Latitude contour plot #####################################
        elif (ptype == 'Long-Lat_Contour'):
            if (var_index > 48 and var_index != 50 and var_index !=51):
                print("var_index < 48 or = 50 or = 51 for contour plots!")
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

            add_text=' '+date_str+UT_str+' Altitude='+'{:.2f}'.format(self.alt[i]) \
                +' km'

            varib=np.zeros((j1-j0, k1-k0), dtype=float)

            for k in range(k0, k1):
              for j in range(j0, j1):
                if (var_index < 45):
                    if (math.isnan(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable Nan')
                        return -1
                    elif(math.isinf(self.iditm_arr[k][j][i][var_index])):
                        print('(i,j,k,s)=(',i,j,k,var_index,'), Variable inf')
                        return -1

                if (var_index < 3):
                    #magnetic field in nT
                    varib[j1-1-j][k-k0]=self.iditm_arr[k][j][i][var_index]*1.0e9

                elif (var_index == 3):
                    #electron temperature in K
                    varib[j1-1-j][k-k0]=self.iditm_arr[k][j][i][var_index]

                elif (((var_index-4) % 5 ==0 and var_index < 34) or \
                      (var_index>=34 and var_index<=40)):
                    #density in cm^-3
                    varib[j1-1-j][k-k0]=np.log10(self.iditm_arr[k][j][i][var_index]*1.0e-6)

                elif ((var_index > 3 and var_index < 34) and (var_index-8) % 5 ==0):
                    #temperature in K
                    varib[j1-1-j][k-k0]= self.iditm_arr[k][j][i][var_index]

                elif (var_index == 48):
                    #neutral temperature in K
                    varib[j1-1-j][k-k0] = self.iditm_arr[k][j][i][var_index]

                elif (var_index == 50):
                    varib[j1-1-j][k-k0]=0.0
                    for l in range(sl):
                        s4=4+5*l
                        if (math.isnan(self.iditm_arr[k][j][i][s4])):
                            print('(i,j,k,s)=(',i,j,k,s4,'), Ion Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s4])):
                            print('(i,j,k,s)=(',i,j,k,s4,'), Ion Density Inf')
                            return -1

                        varib[j1-1-j][k-k0]=varib[j1-1-j][k-k0]+self.iditm_arr[k][j][i][s4]

                    #electron density in cm^-3
                    varib[j1-1-j][k-k0]=np.log10(varib[j1-1-j][k-k0]*1.0e-6)

                elif (var_index == 51):
                    varib[j1-1-j][k-k0]=0.0
                    for l in range(sm):
                        s34=34+l
                        if (math.isnan(self.iditm_arr[k][j][i][s34])):
                            print('(i,j,k,s)=(',i,j,k,s34,'), Neutral Density Nan')
                            return -1
                        elif(math.isinf(self.iditm_arr[k][j][i][s34])):
                            print('(i,j,k,s)=(',i,j,k,s34,'), Neutral Density Inf')
                            return -1

                        varib[j1-1-j][k-k0]=varib[j1-1-j][k-k0]+self.iditm_arr[k][j][i][s34]

                    #neutral density in cm^-3
                    varib[j1-1-j][k-k0]=np.log10(varib[j1-1-j][k-k0]*1.0e-6)

                elif ((((var_index-5) % 5 ==0 or (var_index-6) % 5 ==0 \
                      or (var_index-7) % 5 ==0) and var_index <34) \
                      or (var_index >=41 and var_index <=47)):
                    #velocity in m/s
                    varib[j1-1-j][k-k0]=self.iditm_arr[k][j][i][var_index]

            self.canv, self.fig, self.toolbar = \
                long_lat_cont.lon_lat_contour(wd, canv, fig, toolbar, \
                    self.lon, self.lat, gridNumEnt, varib, ptype, var_index, \
                    add_text, fname)
############# completed Latitude - Height contour plot #################################

        return 0
