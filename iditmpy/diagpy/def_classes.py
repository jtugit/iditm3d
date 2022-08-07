import os
import tkinter as tk
from tkinter import filedialog

#save configure parameters to a file
def save_dir(getfile, conf_file, varbox, plottpbox, gridNumEnt):
    fname=getfile.get_filename()
    if len(fname) < 2:
        return -1

    fwd=open(conf_file,"w")
    fwd.write(fname+'\n')
    fwd.write(varbox.get()+'\n')
    fwd.write(plottpbox.get()+'\n')
    fwd.write(gridNumEnt[0].get()+'  '+gridNumEnt[3].get()+'\n')
    fwd.write(gridNumEnt[1].get()+'  '+gridNumEnt[4].get()+'\n')
    fwd.write(gridNumEnt[2].get()+'  '+gridNumEnt[5].get()+'\n')
    fwd.close()

#close botton action function
class close_button():
    def __init__(self, wd, getfile, conf_file, varbox, plottpbox, gridNumEnt):
        self.wd=wd
        self.getfile=getfile
        self.conf_file=conf_file
        self.varbox=varbox
        self.plottpbox=plottpbox
        self.gridNumEnt=gridNumEnt

    def clean_close(self):
        save_dir(self.getfile, self.conf_file, self.varbox, self.plottpbox, \
            self.gridNumEnt)

        self.wd.destroy()

class save_cfg():
    def __init__(self, getfile, conf_file, varbox, plottpbox, gridNumEnt):
        self.getfile=getfile
        self.conf_file=conf_file
        self.varbox=varbox
        self.plottpbox=plottpbox
        self.gridNumEnt=gridNumEnt

    def save_config(self):
        save_dir(self.getfile, self.conf_file, self.varbox, self.plottpbox, \
            self.gridNumEnt)

class file_dialog():
    def __init__(self, filen, file_sel):
        self.filen=filen
        self.init_dir=os.path.dirname(self.filen)
        self.file_sel=file_sel

    def select_file(self):
        fname=filedialog.askopenfilename(initialdir=self.init_dir.strip(), \
              filetypes =(("hdf5 file", "*.h5"),("All Files","*.*")), \
              title = "Choose a file.")

        self.filen=fname

        fbasen=os.path.basename(fname)

        self.file_sel.delete(0, 'end')
        self.file_sel.insert(0, fbasen)

    def get_filename(self):
        return self.filen

    def set_filename(self, filen):
        self.filen=filen

def print_message(mtext):
    popwd=tk.Tk()

    popwd.title('Error!')

    tk.Label(popwd,text=mtext,bg='white').pack()
    tk.Button(popwd,text='Ok',width=1,height=1,bg='grey',font=("Courier", 10), \
              command=popwd.destroy).pack()

    popwd.configure(background='white')

    popwd.mainloop()

def about_message():
    popwd=tk.Tk()
    popwd.geometry('270x230')
    popwd.title('About')
    tk.Frame(popwd,width=270,height=230,bg="",bd=0,highlightbackground="blue", \
             highlightcolor="blue",highlightthickness=5).place(x=0,y=0)

    tk.Label(popwd, text='',bg='white').pack()
    tk.Label(popwd, text='IDIMT Simulation Results Analysis',bg='white').pack()
    tk.Label(popwd, text='Version 0.0',bg='white').pack()
    tk.Label(popwd, text='',bg='white').pack()
    tk.Label(popwd, text='Author: Jiannan Tu',bg='white').pack()
    tk.Label(popwd, text='2020.11',bg='white').pack()
    tk.Label(popwd, text='',bg='white').pack()
    tk.Label(popwd, text='Space Science Laboratory',bg='white').pack()
    tk.Label(popwd, text='University of Massachusetts Lowell',bg='white').pack()
    tk.Label(popwd, text='',bg='white').pack()

    btn=tk.Button(popwd,text='Ok',width=1,height=1,bg='grey',font=("Courier", 10), \
                  command=popwd.destroy)
    btn.place(x=125,y=192)

    popwd.configure(background='white')

    popwd.mainloop()

class directory_dialog():
    def __init__(self,fdir):
        self.init_dir=fdir
        self.fpath=fdir

    def select_dir(self):
        self.fpath=filedialog.askdirectory(initialdir=self.init_dir, \
                                           title="Choose a director.")
        if len(self.fpath) == 0:
            print_message("No file is choosen or exists")
            return -1

    def get_dirname(self):
        return self.fpath
########################################################################
