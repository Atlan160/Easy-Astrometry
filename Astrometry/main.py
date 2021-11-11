
import numpy as np
import matplotlib.pyplot as plt
import math
from functools import partial
import os
from pathlib import Path

import ipywidgets as widgets
from IPython.display import display

import astromath #my own module
import astrometry #my own class
import spectroscopy # my own class
import save #my own module
import star #my own module
import tooltip #not my own module, from the internet

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

ametry=astrometry.astrometry()
scopy=spectroscopy.spectroscopy()

def update_parameters(_fwhm, _sigma, _threshold, _distance_tolerance_arcsec, _moving_star_tolerance,_number): 
    ametry.set_parameters(float(_fwhm.get()),float(_sigma.get()),float(_threshold.get()),float(_distance_tolerance_arcsec.get()),float(_moving_star_tolerance.get()))
    n=int(_number.get())
    ametry.test_plot(n)



def open_lights():

    try:
        lights_path=filedialog.askopenfilenames(initialdir =" ", title = "Select light files",filetypes = (("newly solved files",".new"),("fit files","*.fit"),("fits files","*.fits")))       
        ametry.import_lights(lights_path)
        scopy.import_lights(lights_path)
        messagebox.showinfo("success","imported "+str(len(lights_path))+" files.")
        #root.mainloop()
    except:
        messagebox.showerror("Error", "something went wrong importing light files")
        #root.mainloop()

def open_dark():
    try:
        dark_path=filedialog.askopenfilenames(initialdir= " ",title ="Select dark file", filetypes = (("fits files","*.fits"),("fit file","*.fit")) )        
        if dark_path == "":
            messagebox.showinfo(" ","imported no files, median correction will then be done")
            ametry.perform_median_correction()
        else:
            ametry.import_dark(dark_path)
            ametry.perform_dark_correction()
            messagebox.showinfo("success","imported dark files")
    except:
        ametry.perform_median_correction()
        messagebox.showerror("Error", "something went wrong importing dark files, median correction will then be done")
        
  

def set_settings_tab():


    ###############
    # w = widgets.IntSlider()
    # w.value=5
    # w.min=0
    # w.max=10
    # display(w)



    ###############


    t_fwhm=StringVar() #textvariables, will be updated each time
    t_fwhm.set(str(ametry.fwhm))
    t_sigma=StringVar()
    t_sigma.set(str(ametry.sigma))
    t_thr=StringVar()
    t_thr.set(str(ametry.threshold))
    t_dist_tolerance=StringVar()
    t_dist_tolerance.set(str(ametry.distance_tolerance_arcsec))
    t_number=StringVar()
    t_number.set(str(ametry.get_test_file_number()))
    t_moving_star_tolerance=StringVar()
    t_moving_star_tolerance.set(str(ametry.moving_star_tolerance))
    
    _update_parameters=partial(update_parameters,t_fwhm,t_sigma,t_thr,t_dist_tolerance,t_moving_star_tolerance,t_number)
    _update_parameters()
    lfwhm=Label(root, text="fwhm in pixel")
    lfwhm.grid(row=0,column=0)
    Entry(root, textvariable=t_fwhm).grid(row=0,column=1)

    lsigma=Label(root,text="sigma")
    lsigma.grid(row=1,column=0)
    Entry(root, textvariable=t_sigma).grid(row=1,column=1)

    lthreshold=Label(root,text="threshold")
    lthreshold.grid(row=2,column=0)
    Entry(root, textvariable=t_thr).grid(row=2,column=1)
    

    ldist=Label(root,text="distance tolerance of stars in arcsec")
    ldist.grid(row=3,column=0)
    Entry(root, textvariable=t_dist_tolerance).grid(row=3,column=1)    

    lmstartol=Label(root,text="moving star tolerance")
    lmstartol.grid(row=5,column=0)
    Entry(root, textvariable=t_moving_star_tolerance).grid(row=5,column=1)

    limnum=Label(root,text="Image number for test")
    limnum.grid(row=4,column=0)
    Scale(root,from_=1, to=ametry.get_number_of_lights(), variable=t_number, orient=HORIZONTAL).grid(row=4,column=1)


    Button(root, text="set and test settings",command=_update_parameters).grid(row=6,column=0)

    __=tooltip.CreateToolTip(lfwhm, tooltip.tooltip_fwhm)
    __=tooltip.CreateToolTip(lsigma, tooltip.tooltip_sigma)
    __=tooltip.CreateToolTip(lmstartol,tooltip.tooltip_tolerance)
    root.mainloop()

def hello():
    print("Hello world")

root=Tk()
menubar=Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Open lights", command=open_lights)
filemenu.add_command(label="open darks", command=open_dark)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)

runmenu=Menu(root)
runmenu=Menu(menubar,tearoff=0)
runmenu.add_command(label="set parameters", command=set_settings_tab)
runmenu.add_command(label="plot all imported lights", command=ametry.plot_all_sources)
runmenu.add_command(label="search for stars and write to file", command=ametry.find_sources)
runmenu.add_command(label="search for same stars and write to file",command=ametry.search_find)
runmenu.add_command(label="search for moving stars", command=ametry.search_for_moving_stars)
menubar.add_cascade(label="astrometry", menu=runmenu)

runmenu=Menu(root)
runmenu=Menu(menubar,tearoff=0)
runmenu.add_command(label="show spectra", command=scopy.show_spectra)
runmenu.add_command(label="plot spectra", command=scopy.plot_spectra)
menubar.add_cascade(label="spectroscopy", menu=runmenu)

# helpmenu = Menu(menubar, tearoff=0)
# helpmenu.add_command(label="About", command=hello)
# menubar.add_cascade(label="Help", menu=helpmenu)
root.config(menu=menubar)
root.mainloop()