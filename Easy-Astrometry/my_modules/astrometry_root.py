import numpy as np

import matplotlib.pyplot as plt

from matplotlib.figure import Figure

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk

import math
from functools import partial
import os

from pathlib import Path

import traceback


import ipywidgets as widgets

from IPython.display import display


from tkinter import *

from tkinter import filedialog

from tkinter import messagebox


from my_modules.tooltip import CreateToolTip, _Tooltip_strings #not my own module, from the internet
from my_modules.astrometry import astrometry



#plt.ion()


class astrometry_root():


    def __init__(self):


        self.astrometry_plot_class=None

        self.easy_astrometry_GUI=Tk(className="Easy-Astrometry")

        width, height = self.easy_astrometry_GUI.winfo_screenwidth(), self.easy_astrometry_GUI.winfo_screenheight()
        self.easy_astrometry_GUI.geometry(str(width)+"x"+str(height)+"+0+0")

        #self.easy_image_identifier_GUI.geometry("1900x1000+10+20")
        self.easy_astrometry_GUI.iconbitmap('_images/icons/astrometry_favicon.ico') #@georgios this line can also be uncommentend


        self.Menubar=Menu(self.easy_astrometry_GUI)

        filemenu = Menu(self.Menubar, tearoff=0)

        filemenu.add_command(label="Open file", command=self.open_file)
        filemenu.add_separator()

        filemenu.add_command(label="Exit", command=self.easy_astrometry_GUI.quit)

        self.Menubar.add_cascade(label="File", menu=filemenu)
        self.easy_astrometry_GUI.config(menu=self.Menubar)
        self.easy_astrometry_GUI.mainloop()




    def open_file(self):


        try:

            lights_path=filedialog.askopenfilenames(initialdir =" ", title = "Select light files",filetypes = (("fits files","*.fits"),("newly solved files",".new"),("fit files","*.fit")))       

            self.astrometry_plot_class=astrometry(lights_path,self.easy_astrometry_GUI,self.Menubar)




        except Exception:

            traceback.print_exc()

            messagebox.showerror("Error", "something went wrong importing light files")



    



