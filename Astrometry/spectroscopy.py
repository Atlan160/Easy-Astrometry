
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.ndimage
from functools import partial
import os
from pathlib import Path

from scipy.ndimage.interpolation import rotate

import astromath #my own module
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

class spectroscopy:

    def __init__(self):
        self.dark_path="" #is imported 
        self.lights_path="" #is imported
        self.data=[]


    def set_light_path(self,path):
        self.lights_path=path

    def set_dark_path(self,path):
        self.dark_path=path
        self.perform_dark_correction=True

    def dark_correction(self, data,dark):
        return data-dark


    def show_spectra(self):

        light=fits.open(self.lights_path)
        self.data=light[0].data
        print("shape:",np.shape(self.data))
        plt.figure(1)
        plt.title(os.path.basename(self.lights_path))
        plt.imshow(self.data, cmap='Greys', origin='upper', interpolation='none')
        plt.show()


    def plot_spectra(self):
        _lights=fits.open(self.lights_path)
        light=_lights[0]
        self.data=light.data


        x2=403
        y2=545

        x1=1939
        y1=3098

        # line = [(881, 1266), (352, 384)]
        dx=x2-x1
        dy=y2-y1

        print("dx,dy:",dx,dy)
        d=math.sqrt(dx**2+dy**2)

        angle=math.degrees(math.asin(dx/d))
        print("angle",angle)


        self.data=rotate(self.data,angle)
        plt.plot([x1,x2],[y1,y2],color="red",linewidth=1)
        plt.imshow(self.data, cmap='Greys', origin='upper', interpolation='none')
        plt.show()

        row=1855

        # self.data=self.data/np.max(self.data)
        # self.data+=0.5
        spec=np.zeros(len(self.data[:,0]),dtype=float)
        for i in range(-3,4):
            spec+=self.data[:,row+i]/7

            plt.title("single spectra")
            plt.plot(self.data[:,row+i])
            plt.show()
        plt.title("summated spectra")
        plt.plot(spec)
        plt.show()


        # Coordinates of the line we'd like to sample along


        # Generate some data...
        # x=np.arange(0,len(self.data[1]),1)
        # y=np.arange(0,len(self.data[0]),1)
        # print(x,y)
        # x,y=np.meshgrid(x,y)



        # # Coordinates of the line we'd like to sample along
        #

        # # Convert the line to pixel/index coordinates
        # x_world, y_world = np.array(list(zip(*line)))
        # col = self.data.shape[1] * (x_world - x.min()) / x.ptp()
        # row = self.data.shape[0] * (y_world - y.min()) / y.ptp()

        # # Interpolate the line at "num" points...
        # num = 1000
        # row, col = [np.linspace(item[0], item[1], num) for item in [row, col]]

        # # Extract the values along the line, using nearest-neighbour
        # zi = self.data[row.astype(int), col.astype(int)]


        # plt.plot(zi)
        # plt.show()