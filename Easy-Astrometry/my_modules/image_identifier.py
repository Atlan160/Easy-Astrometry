from msilib.schema import CheckBox
from tkinter import LEFT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import axes, plot, ion, show
import matplotlib.lines as lines
import math
from functools import partial
import os
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk

import my_modules.astromath as astromath #my own module
from my_modules.calibration import calibration # my own module
import my_modules.save as save #my own module
import my_modules.star as star #my own module
from my_modules.tooltip import CreateToolTip #not my own module, from the internet


import photutils
from photutils import datasets
from photutils import DAOStarFinder
from photutils import CircularAperture

import ipywidgets as widgets
from IPython.display import display
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


class image_identifier(calibration):

    def __init__(self,lights_path,GUI,menubar):
        super().__init__()
        self.GUI=GUI #TK interface
        self.menubar=menubar
        self.import_lights(lights_path)
        self.data=self.lights[0] #inherited from calibration
        self.header=self.headers[0] #same

        self.figure=plt.figure(figsize=(20,12),dpi=100)
        self.axes=plt.axes()
        self.axes.imshow(np.log(self.data)/np.log(20),cmap="Greys_r")
        ##############################
        self.canvas=FigureCanvasTkAgg(self.figure, self.GUI)
        self.canvas.get_tk_widget().grid(row=0,column=10,rowspan=50,sticky=E)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events() # pauses event loop until next event is triggered

        self.line_container=[]
        self.circle_container=[]
        self.text_container=[]

        #[x,y]
        self.clicked=[0,0]
        self.released=[0,0]
        self.current_position=[0,0]


        self.clickstate_draw=0 #1 is click#1 and 2 is click#2

        ### States ####
        self.draw_line_state=False
        self.draw_circles_state=False
        self.draw_point_state=False
        self.snap_on_stars=False
        self.mouse_on_hold=False

        ### image properties ###
        self.image_scale_arcsec_per_pixel=self.header["scale"]

        ### Textvariables ###
        self.stringvar_RA=StringVar(self.GUI)
        self.stringvar_RA.set(str(0))
        self.stringvar_DEC=StringVar(self.GUI)
        self.stringvar_DEC.set(str(0))
        self.stringvar_current_coordinates_RA=StringVar(self.GUI)
        self.stringvar_current_coordinates_RA.set(str(0))
        self.stringvar_current_coordinates_DEC=StringVar(self.GUI)
        self.stringvar_current_coordinates_DEC.set(str(0))


        # setup GUI
        self.init_GUI()

        self.figure.canvas.mpl_connect('button_press_event', self.onclick)
        self.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.figure.canvas.mpl_connect('motion_notify_event', self.mouse_moved)
        self.figure.canvas.mpl_connect('scroll_event',self.mousewheel_moved)

    #only for test
    def hello(self):
        #for test
        print("Hello World")

    def find_nearby_star(self,x,y,tol=20):
        """
        Finds nearby stars in proximity of coordinates x,y.
        returns the coordinates of the proximate star.


        Args:
            x (float): [description]
            y (float): [description]
            tol (int,float ): tolerance, possible distance of star to coordinates x,y. Defaults to 15.
        """



        _mean, median, std = sigma_clipped_stats(self.data, sigma=4.0) #get data statistic

        daofind = DAOStarFinder(fwhm=4, threshold=9*std)
        sources=daofind(self.data)
        index=-1
        distance=100000 #should be more than enough for any image

        for i in range(len(sources["id"])):

            if astromath.return_distance_pixel(x,y,sources["xcentroid"][i],sources["ycentroid"][i])<distance:
                distance=astromath.return_distance_pixel(x,y,sources["xcentroid"][i],sources["ycentroid"][i])
                index=i
        
        if distance<tol:
            return sources["xcentroid"][index],sources["ycentroid"][index]
        else:
            messagebox.showinfo("Error","couldnt find a star at the given position")
            return -1,-1


    def onclick(self,event):
        """
        triggered when mouse is clicked

        Args:
            event ([type]): event
        """

        print("clicked")
        #print("dtype",type(event.xdata))

        if event.xdata is not None and event.ydata is not None:
            if event.xdata>=0 and event.ydata>=0:
            
                print("after event_check")
                self.clicked[0]=event.xdata 
                self.clicked[1]=event.ydata

                if self.draw_line_state==True and self.clickstate_draw==0:
                    if self.snap_on_stars==True:
                        x,y=self.find_nearby_star(self.clicked[0],self.clicked[1])

                        if x>=0 and y>=0:                            
                            self.draw_line(self.clicked[0],self.clicked[1])
                            self.clickstate_draw+=1
                    else:   
                        self.clickstate_draw+=1
                        self.draw_line(self.clicked[0],self.clicked[1])

                else:
                    self.clickstate_draw=0
                    self.draw_line_state=False

                if self.draw_point_state==True:

                    if self.snap_on_stars==True:
                        x,y=self.find_nearby_star(self.clicked[0],self.clicked[1])
                        if x>=0 and y>=0:                            
                            self.draw_point(x,y)

                    else:
                        self.draw_point(event.xdata,event.ydata)

                    self.draw_point_state=False
        
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()






        

    def onrelease(self,event):
        pass
        #check if click is in data

    def mouse_moved(self,event):
        """triggered when mouse is moved

        Args:
            event (event): event
        """

        if event.xdata is not None and event.ydata is not None:
            if event.xdata>=0 and event.ydata>=0:
    
                self.current_position[0]=event.xdata 
                self.current_position[1]=event.ydata
                RA,DEC=astromath.return_coordinates_RA_DEC(self.header,event.xdata,event.ydata)
                self.stringvar_current_coordinates_RA.set(str(np.round(RA,5)))
                self.stringvar_current_coordinates_DEC.set(str(np.round(DEC,5)))                

                self.update_plot()
        #check if movement is in data

    def mousewheel_moved(self,event):
        print("mousewheelevent",event.button)
        if event.button=="down":
            #zoom out
            xlim=self.axes.get_xlim()
            ylim=self.axes.get_ylim()
            self.axes.set_xlim(xlim[0]*0.9,xlim[1]*0.9)
            self.axes.set_ylim(ylim[0]*0.9,ylim[1]*0.9)


        if event.button=="down":
            #zoom in
            xlim=self.axes.get_xlim()
            ylim=self.axes.get_ylim()
            self.axes.set_xlim(xlim[0]*1/0.9,xlim[1]*1/0.9)
            self.axes.set_ylim(ylim[0]*1/0.9,ylim[1]*1/0.9)



    def update_plot(self):
        """
        updates the current line so that it is drawn from first click position to current mouse position. 
        is updated whenever mouse is moved
        """

        if self.draw_line_state==True and self.clickstate_draw==1 and len(self.line_container)>0 and self.snap_on_stars==False:
            x_current_line=self.line_container[-1].get_xdata(orig=False)
            y_current_line=self.line_container[-1].get_ydata(orig=False)
            length_pixel=astromath.return_distance_pixel(x_current_line[0],y_current_line[0],self.current_position[0],self.current_position[1])
            length_arcsec=length_pixel*self.image_scale_arcsec_per_pixel
            self.line_container[-1].set_data([x_current_line[0],self.current_position[0]],[y_current_line[0],self.current_position[1]])
            self.text_container[-1].set_position((x_current_line[0]+(self.current_position[0]-x_current_line[0])/2,y_current_line[0]+(self.current_position[1]-y_current_line[0])/2))
            self.text_container[-1].set_text("length in arcmin "+str(length_arcsec/60),color="deepskyblue")

        elif self.draw_line_state==True and self.clickstate_draw==1 and len(self.line_container)>0 and self.snap_on_stars==True:
            x_current_line=self.line_container[-1].get_xdata(orig=False)
            y_current_line=self.line_container[-1].get_ydata(orig=False)
            x,y=self.find_nearby_star(self.clicked[0],self.clicked[1],tol=100000)
            if x>=0 and y>=0:                            

                length_pixel=astromath.return_distance_pixel(x_current_line[0],y_current_line[0],x,x)
                length_arcsec=length_pixel*self.image_scale_arcsec_per_pixel
                self.line_container[-1].set_data([x_current_line[0],x],[y_current_line[0],y])
                #TODO needs a rework of position
                self.text_container[-1].set_position((x_current_line[0]+(x-x_current_line[0])/2,y_current_line[0]+(y-y_current_line[0])/2))
                self.text_container[-1].set_text("length in arcmin "+str(length_arcsec/60),color="deepskyblue")
            
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()


    def draw_line(self,x,y):
        """
        adds a new line to current axes with length 0. line is updated when mouse is moved
        """
        print("drawing line")

        l=lines.Line2D([x,x],[y,y],color="deepskyblue")
        self.axes.add_line(l)
        self.line_container.append(l)
        self.text_container.append(self.axes.text(0,0,""))

    def find_draw_coordinates(self):
        """
        draws the coordinates typed in the textfield RA,Dec on the image as a circle
        
        """

        RA=float(self.stringvar_RA.get())
        DEC=float(self.stringvar_DEC.get())
        X,Y=astromath.return_X_Y_coordinates(self.header,RA,DEC)
        #TODO X,Y coordinates are always the same...
        self.axes.scatter(X,Y,color="None",edgecolors="orangered",alpha=0.5)
        text="RA:"+str(np.round(RA,6))+" DEC:"+str(np.round(DEC,6))
        self.axes.text(X+10,Y,text,color="orangered")
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
        #self.figure.canvas.draw()


    def set_drawline(self):
        self.draw_line_state=True

    def set_drawpoint(self):
        self.draw_point_state=True
        
    def draw_point(self,x,y):
        self.axes.scatter(x,y,color="None",edgecolors="springgreen",alpha=0.5)
        RA,DEC=astromath.return_coordinates_RA_DEC(self.header,x,y)
        text="RA:"+str(np.round(RA,6))+" DEC:"+str(np.round(DEC,6))
        self.axes.text(x+10,y,text,color="springgreen")
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()



    def set_drawcircle(self):
        self.draw_circles_state=True

    def invert_snap_on_stars(self):
        self.snap_on_stars= not self.snap_on_stars
        print("snap on stars is",self.snap_on_stars)

    def init_GUI(self):
        """
        generates all the graphical elements when an image is loaded

        """

        calibration_menu = Menu(self.menubar, tearoff=0)
        calibration_menu.add_command(label="Show image information", command=self.hello)
        calibration_menu.add_separator()
        calibration_menu.add_command(label="Set image object distance", command=self.hello)
        self.menubar.add_cascade(label="Image calibration", menu=calibration_menu)
        self.GUI.config(menu=self.menubar)


        Label_drawing=Label(self.GUI,text="Drawing")
        Label_drawing.grid(row=0,column=0)
        Line_draw_button=Button(self.GUI,text="Line",command=self.set_drawline)
        Line_draw_button.grid(row=1,column=0)
        Point_draw_button=Button(self.GUI,text="Point",command=self.set_drawpoint)
        Point_draw_button.grid(row=1,column=1)
        Circle_draw_button=Button(self.GUI,text="Circle",command=self.set_drawcircle)
        Circle_draw_button.grid(row=1,column=2)
        snap_stars_checkbox=Checkbutton(self.GUI,text="snap on stars",onvalue=1, offvalue=0,command=self.invert_snap_on_stars)
        snap_stars_checkbox.deselect()
        snap_stars_checkbox.grid(row=1,column=3)

        Label_coordinates=Label(self.GUI,text="Find Coordinates")
        Label_coordinates.grid(row=2,column=0)
        Label_coordinates_RA=Label(self.GUI,text="RA in decimals:")
        Label_coordinates_RA.grid(row=3,column=0)
        Label_coordinates_DEC=Label(self.GUI,text="DEc in decimals")
        Label_coordinates_DEC.grid(row=4,column=0)

        textfield_RA=Entry(self.GUI, textvariable=self.stringvar_RA)
        textfield_RA.grid(row=3,column=1)
        textfield_DEC=Entry(self.GUI, textvariable=self.stringvar_DEC)
        textfield_DEC.grid(row=4,column=1)
        Find_coord_button=Button(self.GUI,text="Find coordinates",command=self.find_draw_coordinates)
        Find_coord_button.grid(row=5,column=1)
        print("entry is",self.stringvar_DEC.get())

        Label_show_coords=Label(self.GUI,text="Current mouse coordinates in degree")
        Label_show_coords.grid(row=6,column=0)
        Label_show_coords_RA=Label(self.GUI,textvariable=self.stringvar_current_coordinates_RA)
        Label_show_coords_RA.grid(row=7,column=0)
        Label_show_coords_DEC=Label(self.GUI,textvariable=self.stringvar_current_coordinates_DEC)
        Label_show_coords_DEC.grid(row=7,column=1)
        # dont use this function otherwise program is traped in this mainloop
        #self.GUI.mainloop()