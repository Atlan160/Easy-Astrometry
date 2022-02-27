
import numpy as np
import matplotlib.pyplot as plt
import math
from functools import partial
import os
from pathlib import Path

import ipywidgets as widgets
from IPython.display import display

from my_modules.astrometry import astrometry
from my_modules.tooltip import CreateToolTip, _Tooltip_strings #not my own module, from the internet

import photutils
from photutils import datasets
from photutils import DAOStarFinder
from photutils import CircularAperture

import astropy
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import PercentileInterval
from astropy.visualization import simple_norm
from astropy.table import Table

from tkinter import *
from tkinter import filedialog
from tkinter import messagebox




class astrometry_root():

    def __init__(self):
        self.ametry=astrometry()
        Tooltip_strings=_Tooltip_strings()
        self.easy_astrometry_root=Tk(className=" Easy-Astrometry")
        self.easy_astrometry_root.geometry("300x200+400+200")
        self.easy_astrometry_root.iconbitmap('_images/icons/astrometry_favicon.ico')

        menubar=Menu(self.easy_astrometry_root)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open lights", command=self.open_lights)
        filemenu.add_command(label="Open darks", command=self.open_dark)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.easy_astrometry_root.quit)
        menubar.add_cascade(label="File", menu=filemenu)


        runmenu=Menu(menubar,tearoff=0)
        runmenu.add_command(label="Set parameters", command=self.set_settings_tab)
        runmenu.add_command(label="Plot all imported lights", command=self.ametry.plot_all_sources)
        runmenu.add_command(label="Search for stars and write to file", command=self.ametry.find_sources)
        runmenu.add_command(label="Search for same stars and write to file",command=self.ametry.search_find)
        runmenu.add_command(label="Search for moving targets", command=self.ametry.search_for_moving_stars)
        runmenu.add_command(label="Search for moving targets relative", command=self.ametry.search_without_platesolve)
        menubar.add_cascade(label="Astrometry", menu=runmenu)

        self.easy_astrometry_root.config(menu=menubar)
        self.easy_astrometry_root.mainloop()

    def update_parameters(self,_fwhm, _sigma, _threshold, _distance_tolerance_arcsec, _moving_star_tolerance,_number): 
        self.ametry.set_parameters(float(_fwhm.get()),float(_sigma.get()),float(_threshold.get()),float(_distance_tolerance_arcsec.get()),float(_moving_star_tolerance.get()))
        n=int(_number.get())
        self.ametry.test_plot(n)



    def open_lights(self):

        try:
            lights_path=filedialog.askopenfilenames(initialdir =" ", title = "Select light files",filetypes = (("newly solved files",".new"),("fit files","*.fit"),("fits files","*.fits")))       
            self.ametry.import_lights(lights_path)
            self.ametry.same_stars_searched=False
            self.ametry.moving_stars_searched=False
            #messagebox.showinfo("success","imported "+str(len(lights_path))+" files.")
            #root.mainloop()
        except:
            messagebox.showerror("Error", "something went wrong importing light files")
            #root.mainloop()

    def open_dark(self):
        try:
            dark_path=filedialog.askopenfilenames(initialdir= " ",title ="Select dark file", filetypes = (("fits files","*.fits"),("fit file","*.fit")) )        
            if dark_path == "":
                messagebox.showinfo(" ","imported no files, median correction will then be done")
                self.ametry.perform_median_correction()
            else:
                self.ametry.import_dark(dark_path)
                self.ametry.perform_dark_correction()
                messagebox.showinfo("success","imported dark files")
        except:
            self.ametry.perform_median_correction()
            messagebox.showerror("Error", "something went wrong importing dark files, median correction will then be done")
            
    

    def set_settings_tab(self):


        ###############
        # w = widgets.IntSlider()
        # w.value=5
        # w.min=0
        # w.max=10
        # display(w)



        ###############


        t_fwhm=StringVar() #textvariables, will be updated each time
        t_fwhm.set(str(self.ametry.fwhm))
        t_sigma=StringVar()
        t_sigma.set(str(self.ametry.sigma))
        t_thr=StringVar()
        t_thr.set(str(self.ametry.threshold))
        t_dist_tolerance=StringVar()
        t_dist_tolerance.set(str(self.ametry.distance_tolerance_arcsec))
        t_number=StringVar()
        t_number.set(str(self.ametry.get_test_file_number()))
        t_moving_star_tolerance=StringVar()
        t_moving_star_tolerance.set(str(self.ametry.moving_star_tolerance))
        
        _update_parameters=partial(self.update_parameters,t_fwhm,t_sigma,t_thr,t_dist_tolerance,t_moving_star_tolerance,t_number)
        _update_parameters()
        lfwhm=Label(self.easy_astrometry_root, text="fwhm in pixel")
        lfwhm.grid(row=0,column=0)
        Entry(self.easy_astrometry_root, textvariable=t_fwhm).grid(row=0,column=1)

        lsigma=Label(self.easy_astrometry_root,text="sigma")
        lsigma.grid(row=1,column=0)
        Entry(self.easy_astrometry_root, textvariable=t_sigma).grid(row=1,column=1)

        lthreshold=Label(self.easy_astrometry_root,text="threshold")
        lthreshold.grid(row=2,column=0)
        Entry(self.easy_astrometry_root, textvariable=t_thr).grid(row=2,column=1)
        

        ldist=Label(self.easy_astrometry_root,text="distance tolerance of stars in arcsec")
        ldist.grid(row=3,column=0)
        Entry(self.easy_astrometry_root, textvariable=t_dist_tolerance).grid(row=3,column=1)    

        lmstartol=Label(self.easy_astrometry_root,text="moving star tolerance")
        lmstartol.grid(row=5,column=0)
        Entry(self.easy_astrometry_root, textvariable=t_moving_star_tolerance).grid(row=5,column=1)

        limnum=Label(self.easy_astrometry_root,text="Image number for test")
        limnum.grid(row=4,column=0)
        Scale(self.easy_astrometry_root,from_=1, to=self.ametry.get_number_of_lights(), variable=t_number, orient=HORIZONTAL).grid(row=4,column=1)


        Button(self.easy_astrometry_root, text="set and test settings",command=_update_parameters).grid(row=6,column=0)

        __=CreateToolTip(lfwhm, self.Tooltip_strings.tooltip_fwhm)
        __=CreateToolTip(lsigma, self.Tooltip_strings.tooltip_sigma)
        __=CreateToolTip(lmstartol,self.Tooltip_strings.tooltip_tolerance)
        self.easy_astrometry_root.mainloop()

    def hello(self): # just for testing
        print("Hello world")
