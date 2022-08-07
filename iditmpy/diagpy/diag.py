#!/usr/bin/python3
import os
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

import def_classes
import class_plot

def main():
    ################## Read configureation file ############################
    conf_file=os.getcwd()
    if conf_file[len(conf_file)-1] == '/':
        conf_file=conf_file.strip()+'init_dir.cfg'
    else:
        conf_file=conf_file.strip()+'/init_dir.cfg'

    if not os.path.exists(conf_file.strip()):
        mtext='No such file: '+conf_file+'. Use default directory'
        def_classes.print_message(mtext)
    else:
        frd=open(conf_file,"r")

        fname=frd.readline()
        fname=fname.strip()
        fbase=os.path.basename(fname)

        variable=frd.readline()
        variable=variable.strip()

        ptype=frd.readline()
        ptype=ptype.strip()

        line=frd.readline()
        astr=line.split()
        alt_num=[astr[0].strip(), astr[1].strip()]

        line=frd.readline()
        astr=line.split()
        lat_num=[astr[0].strip(), astr[1].strip()]

        line=frd.readline()
        astr=line.split()
        lon_num=[astr[0].strip(), astr[1].strip()]

        frd.close()

    ################## End Read configureation file ############################

    rwd=tk.Tk()
    rwd.title('IDITM Simulation Results Analysis')
    rwd.geometry('830x704')

    ################## Drop boxs and entry fields ##########################
    varList=['e_O+ coll','e_O2+ coll','e_N2+ coll','e_H+ coll', \
            'e_He+ coll','e_NO+ coll', \
            'O+_e coll','O+_O2+ coll','O+_N2+ coll','O+_H+ coll','O+_He+ coll','O+_NO+ coll', \
            'O+_O coll','O+_O2 coll','O+_N2 coll','O+_H coll','O+_He coll', \
            'O2+_e coll','O2+_O+ coll','O2+_N2+ coll','O2+_H+ coll','O2+_He+ coll','O2+_NO+ coll', \
            'O2+_O coll','O2+_O2 coll','O2+_N2 coll','O2+_H coll','O2+_He coll', \
            'N2+_e coll','N2+_O+ coll','N2+_O2+ coll','N2+_H+ coll','N2+_He+ coll','N2+_NO+ coll', \
            'N2+_O coll','N2+_O2 coll','N2+_N2 coll','N2+_H coll','N2+_He coll', \
            'H+_e coll','H+_O+ coll','H+_O2+ coll','H+_N2+ coll','H+_He+ coll','H+_NO+ coll', \
            'H+_O coll','H+_O2 coll','H+_N2 coll','H+_H coll','H+_He coll', \
            'He+_e coll','He+_O+ coll','He+_O2+ coll','He+_N2+ coll','He+_H+ coll','He+_NO+ coll', \
            'He+_O coll','He+_O2 coll','He+_N2 coll','He+_H coll','He+_He coll', \
            'NO+_e coll','NO+_O+ coll','NO+_O2+ coll','NO+_N2+ coll','NO+_H+ coll','NO+_He+ coll', \
            'NO+_O coll','NO+_O2 coll','NO+_N2 coll','NO+_H coll','NO+_He coll', \
            'O_O2 coll','O_N2 coll','O2_N2 coll','H_O coll','H_O2 coll','H_N2 coll', \
            'H_He coll','He_O coll','He_O2 coll','He_N2 coll', \
            'e_O coll','e_O2 coll', 'e_N2 coll', 'e_H coll','e_He coll', \
            'e_Heating_rate','e_Cooling_rate','Neu_Heating_rate','Neu_Cooling_rate', \
            'e_therm-cond','O+_therm-cond','O2+_therm-cond', \
            'N2+_therm-cod','H+_therm-cond','He+_therm-cond', \
            'NO+_therm-cond','Neu_therm-cond', 'All e coll','All O+ coll', \
            'All O2+ coll','All N2+ coll','All H+ coll','All He+ coll','All NO+ coll', \
            'Ion_themal-cond']

    for i in range(len(varList)):
        if variable == varList[i]:
            break

    tk.Label(rwd,text='Var to Plot',bg='white').place(x=15,y=5)
    varbox=ttk.Combobox(rwd,values=varList,width=10)
    varbox.place(x=10,y=25)
    varbox.current(i)

    ptypeList=['Altitude_Profile','Latitude_Profile','Longitude_Profile', \
        'Lat-Height_Contour','Long-Lat_Contour','Long-Height_Cont']
    for i in range(6):
        if ptype == ptypeList[i]:
            break

    tk.Label(rwd,text='Plot Type',bg='white').place(x=130,y=5)
    plottpbox=ttk.Combobox(rwd,values=ptypeList,width=14)
    plottpbox.place(x=110,y=25)
    plottpbox.current(i)

    gridNumEnt=[]*6
    tk.Label(rwd,text='AltIndex',bg='white').place(x=234,y=5)
    tk.Label(rwd,text='Start',bg='white').place(x=245,y=45)
    gridNumEnt.append(tk.Entry(rwd))
    gridNumEnt[0].place(x=240,y=25,width=40)
    gridNumEnt[0].delete(0, 'end')
    gridNumEnt[0].insert(0, alt_num[0])

    tk.Label(rwd,text='CoLatIndex',bg='white').place(x=287,y=5)
    tk.Label(rwd,text='Start',bg='white').place(x=301,y=45)
    gridNumEnt.append(tk.Entry(rwd))
    gridNumEnt[1].place(x=296,y=25,width=40)
    gridNumEnt[1].delete(0, 'end')
    gridNumEnt[1].insert(0, lat_num[0])

    tk.Label(rwd,text='LonIndex',bg='white').place(x=349,y=5)
    tk.Label(rwd,text='Start',bg='white').place(x=356,y=45)
    gridNumEnt.append(tk.Entry(rwd))
    gridNumEnt[2].place(x=352,y=25,width=40)
    gridNumEnt[2].delete(0, 'end')
    gridNumEnt[2].insert(0, lon_num[0])

    tk.Label(rwd,text='AltIndex',bg='white').place(x=411,y=5)
    tk.Label(rwd,text='End',bg='white').place(x=425,y=45)
    gridNumEnt.append(tk.Entry(rwd))
    gridNumEnt[3].place(x=418,y=25,width=40)
    gridNumEnt[3].delete(0, 'end')
    gridNumEnt[3].insert(0, alt_num[1])

    tk.Label(rwd,text='CoLatIndex',bg='white').place(x=464,y=5)
    tk.Label(rwd,text='End',bg='white').place(x=480,y=45)
    gridNumEnt.append(tk.Entry(rwd))
    gridNumEnt[4].place(x=473,y=25,width=40)
    gridNumEnt[4].delete(0, 'end')
    gridNumEnt[4].insert(0, lat_num[1])

    tk.Label(rwd,text='LonIndex',bg='white').place(x=527,y=5)
    tk.Label(rwd,text='End',bg='white').place(x=535,y=45)
    gridNumEnt.append(tk.Entry(rwd))
    gridNumEnt[5].place(x=528,y=25,width=40)
    gridNumEnt[5].delete(0, 'end')
    gridNumEnt[5].insert(0, lon_num[1])

    tk.Label(rwd,text='Data File Selected',bg='white').place(x=600,y=5)
    file_sel=tk.Entry(rwd)
    file_sel.place(x=580,y=25,width=140)
    file_sel.delete(0,'end')
    file_sel.insert(0, fbase)

    ################# End Drop boxs and entry fields #########################

    ################### Menubar and related actions ##########################
    menubar=tk.Menu(rwd)

    # create a pulldown menu, and add it to the menu bar
    filemenu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="File", menu=filemenu)

    helpmenu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Help", menu=helpmenu)
    helpmenu.add_command(label="About", command=def_classes.about_message)

    rwd.configure(menu=menubar,background='white')

    #opendir=def_classes.directory_dialog(fdir)
    #filemenu.add_command(label="Select Directory",command=opendir.select_dir)

    getfile=def_classes.file_dialog(fname, file_sel)
    filemenu.add_command(label="Select A hdf5 file",command=getfile.select_file)

    saveconfig=def_classes.save_cfg(getfile, conf_file, varbox, plottpbox, \
        gridNumEnt)
    filemenu.add_command(label="Save Config", command=saveconfig.save_config)

    exitmenu = def_classes.close_button(rwd, getfile, conf_file, varbox, \
        plottpbox, gridNumEnt)
    filemenu.add_command(label="Exit", command=exitmenu.clean_close)

    ################### End menubar and related actions ########################

#plot button for making plots
    wplotbt=class_plot.plot_button(rwd, getfile, varbox, plottpbox, varList, \
        gridNumEnt, file_sel)
    pbtn=tk.Button(rwd,text='Plot',width=2,height=2,bg='grey',font=("Courier", 10), \
                   command=wplotbt.wplot)
    pbtn.place(x=735,y=5)

#close button for close the GUI
    wclosebt = def_classes.close_button(rwd, getfile, conf_file, varbox, \
        plottpbox, gridNumEnt)
    pbtn=tk.Button(rwd,text='Close',width=2,height=2,bg='grey', \
        font=("Courier", 10), command=wclosebt.clean_close)
    pbtn.place(x=780,y=5)

    rwd.mainloop()

if __name__== "__main__":
    main()
