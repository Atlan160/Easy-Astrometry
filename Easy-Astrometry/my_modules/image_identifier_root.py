
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
from my_modules.image_identifier import image_identifier

#plt.ion()

class image_identifier_root():

    def __init__(self):

        self.image_identifier_plot_class=None

        self.easy_image_identifier_GUI=Tk(className="Easy-Image Identifier")
        width, height = self.easy_image_identifier_GUI.winfo_screenwidth(), self.easy_image_identifier_GUI.winfo_screenheight()
        self.easy_image_identifier_GUI.geometry(str(width)+"x"+str(height)+"+0+0")
        #self.easy_image_identifier_GUI.geometry("1900x1000+10+20")
        self.easy_image_identifier_GUI.iconbitmap('_images/icons/image_identifier_favicon.ico')

        self.Menubar=Menu(self.easy_image_identifier_GUI)
        filemenu = Menu(self.Menubar, tearoff=0)
        filemenu.add_command(label="Open file", command=self.open_file)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.easy_image_identifier_GUI.quit)
        self.Menubar.add_cascade(label="File", menu=filemenu)
        self.easy_image_identifier_GUI.config(menu=self.Menubar)
        self.easy_image_identifier_GUI.mainloop()



    def open_file(self):

        try:
            lights_path=filedialog.askopenfilenames(initialdir =" ", title = "Select one image file",filetypes = (("fits files","*.fits"),("newly solved files",".new"),("fit files","*.fit")))       
            self.image_identifier_plot_class=image_identifier(lights_path,self.easy_image_identifier_GUI,self.Menubar)



        except Exception:
            traceback.print_exc()
            messagebox.showerror("Error", "something went wrong importing light files")



    


