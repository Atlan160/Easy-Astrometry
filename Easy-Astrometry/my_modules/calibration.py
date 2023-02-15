
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, ion, show
import math
from functools import partial
import os
from pathlib import Path


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

class calibration:

    def __init__(self):
        self.light_imported=False
        self.dark_imported=False
        self.dark_corrected=False
        self.dark=[] #is imported 
        self.lights=[] #is imported
        self.headers=[] # is imported
        self.lights_path=""
        self.data_bin2=[]


    def print_hello(self):
        print("hello world")


    def import_lights(self,path):
        _lights=[]
        _headers=[]
        self.lights_path=path #import paths
        self.light_imported=True

        for fi in path:
            file=fits.open(fi)
            _lights.append(file[0].data)
            _headers.append(file[0].header)
        self.lights=_lights
        self.headers=_headers
        print("imported "+str(len(path))+" files")

    def generate_data(self):

        for i in range(len(self.lights)):
            pass
            #self.data_bin2.append(self.data_raw[::2,::2])


    def get_number_of_lights(self):
        return len(self.lights)

    def import_dark(self,path):
        _darks=[]
        for fi in path:
            file=fits.open(fi)
            _darks.append(file[0].data)
        self.darks=_darks
        self.dark_imported=True


    def perform_dark_correction(self):
        n=len(self.darks)
        master_dark=self.darks[0] # init with first dark

        #generate master_dark by stacking all darks
        for i in range(1,len(self.darks)): 
            master_dark+=self.darks[i]
        self.dark=master_dark/n

        for i in range(len(self.lights)):
            self.lights[i]=self.lights[i]-self.dark
        self.dark_corrected=True
        print("dark correction done")
            

    def perform_median_correction(self):
        for i,li in enumerate(self.lights):
            _mean, median, std = sigma_clipped_stats(li, sigma=self.sigma) #get data statistic
            self.lights[i]=self.lights[i]-median

