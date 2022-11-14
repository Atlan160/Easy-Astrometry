from tkinter import LEFT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.pyplot import axes, plot, ion, show
import matplotlib.lines as lines
import math
from functools import partial
import os
import time
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk

import my_modules.astromath as astromath #my own module
from my_modules.calibration import calibration # my own module
import my_modules.save as save #my own module
import my_modules.star as star #my own module
from my_modules.tooltip import CreateToolTip, _Tooltip_strings #not my own module, from the internet


from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import LEFT

import astropy
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import PercentileInterval
from astropy.visualization import simple_norm
from astropy.table import Table

import photutils
from photutils import datasets
from photutils import DAOStarFinder
from photutils import CircularAperture

import ipywidgets as widgets
from IPython.display import display


class image_container(calibration):

    def __init__(self,lights_path,GUI,menubar):
        super().__init__()
        self.GUI=GUI #TK interface
        ### GUI elements and figure ###
        self.GUI_Frame=Frame(GUI)
        self.Image_Frame=Frame(self.GUI_Frame)
        self.Elements_Frame=Frame(self.GUI_Frame)
        self.menubar=menubar
        self.figure=plt.figure(figsize=(20,12),dpi=100,frameon=True,constrained_layout=True)
        self.axes=plt.axes()

        ###  getting image data, inherited from calibration ###
        self.import_lights(lights_path)
        self.data=self.lights[0] #inherited from calibration
        self.header=self.headers[0] #same


        self.canvas=FigureCanvasTkAgg(self.figure, self.Image_Frame)
        self.canvas.get_tk_widget().grid(row=0,column=0)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events() # pauses event loop until next event is triggered

        #connect mpls
        self.figure.canvas.mpl_connect('button_press_event', self.onclick)
        self.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.figure.canvas.mpl_connect('motion_notify_event', self.mouse_moved)
        self.figure.canvas.mpl_connect('scroll_event',self.mousewheel_moved)


    def onclick(self):
        pass

    def onrelease(self):
        pass

    def mouse_moved(self):
        pass

    def mousewheel_moved(self):
        pass
