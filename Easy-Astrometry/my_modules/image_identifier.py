import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.pyplot import axes, plot, ion, show
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import matplotlib.lines as lines

import math
from functools import partial
import os
import time
from pathlib import Path


import my_modules.astromath as astromath #my own module
from my_modules.calibration import calibration # my own module
import my_modules.save as save #my own module
import my_modules.star as star #my own module
from my_modules.image_container import image_container
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

class image_identifier(image_container):

    def __init__(self,lights_path,GUI,menubar):

        ### getting calibration constructor and setting up GUI ###
        super().__init__(lights_path,GUI,menubar)

        self.Tooltip_strings=_Tooltip_strings() # all the text for the tooltips are stored in a class

        ### mouse positions [x,y] ###
        self.clicked=[0,0]
        self.released=[0,0]
        self.current_position=[0,0]

        self.clickstate_draw_line=0 #1 is click#1 and 2 is click#2 TODO check if one can change that to a boolean variable
        self.clickstate_draw_circle=0 #1 is click#1 and 2 is click#2

        ### States ####
        self.draw_line_state=False
        self.draw_circle_state=False
        self.draw_point_state=False
        self.snap_on_stars=False
        self.mouse_on_hold=False
        self.show_magnitude_state=False
        self.show_orientation_state=True
        self.calibrate_magnitude_state=False


        ### image properties ###
        self.header=self.headers[0] #import header of first image (only the first image is plotted)
        self.data=self.lights[0].copy() #import data
        self.image_scale_arcsec_per_pixel=self.header["scale"]
        self.sources=self.find_stars() # stars in image


        ### Textvariables ###
        self.stringvar_RA=StringVar(self.GUI,value=0) #TODO why not use DoubleVar ?
        self.stringvar_DEC=StringVar(self.GUI,value="0")
        self.stringvar_current_coordinates_RA=StringVar(self.GUI,value="0")
        self.stringvar_current_coordinates_DEC=StringVar(self.GUI,value="0")
        self.stringvar_calibrate_magnitude=StringVar(self.GUI,value="0")
        self.doublevar_gamma=DoubleVar(self.GUI,value=1)


        # the actual plot, backend iherited from image_container
        im1=self.axes.imshow(self.data**-self.doublevar_gamma.get(),cmap="Greys",interpolation="none")

        self.line_container=[]
        self.ellipse_container=[]
        self.point_container=[]
        self.text_container=[]
        self.text_magnitude_container=[]
        self.arrow_container=[]


        # setup GUI
        self.init_GUI()


    #only for test
    def hello(self):
        #for test
        print("Hello World")


    def find_stars(self):
        """
        Find all stars in image with given parameters and return them
        """

        _mean, median, std = sigma_clipped_stats(self.data, sigma=4.0) #get data statistic
        daofind = DAOStarFinder(fwhm=4, threshold=9*std)
        return daofind(self.data)



    def find_nearby_star(self,x,y,tol=20):
        """
        Finds nearby stars in proximity of coordinates x,y.
        returns the coordinates of the proximate star.


        Args:
            x (float): [description]
            y (float): [description]
            tol (int,float ): tolerance in pixel, possible distance of star to coordinates x,y. Defaults to 15.
        """


        index=-1
        distance=100000 #should be more than enough for any image

        for i in range(len(self.sources["id"])):

            if astromath.return_distance_pixel(x,y,self.sources["xcentroid"][i],self.sources["ycentroid"][i])<distance:
                distance=astromath.return_distance_pixel(x,y,self.sources["xcentroid"][i],self.sources["ycentroid"][i])
                index=i
        
        if distance<tol:
            return self.sources["xcentroid"][index],self.sources["ycentroid"][index], index
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

        if event.xdata is not None and event.ydata is not None: # just check if coordinates are valid
            if event.xdata>=0 and event.ydata>=0:
            

                self.clicked[0]=event.xdata 
                self.clicked[1]=event.ydata

                if self.draw_line_state==True and self.clickstate_draw_line==0: 
                    """
                    if lines shall be drawn and it is first point
                    """
                    if self.snap_on_stars==True:
                        x,y,_=self.find_nearby_star(self.clicked[0],self.clicked[1])

                        if x>=0 and y>=0:        
                            self.clickstate_draw_line+=1               
                            self.draw_line(x,y)
                            
                    else:   
                        self.clickstate_draw_line+=1
                        self.draw_line(self.clicked[0],self.clicked[1])

                else:
                    self.clickstate_draw_line=0
                    self.draw_line_state=False

                if self.draw_point_state==True: 
                    """
                    if point shall be drawn  
                    """

                    if self.snap_on_stars==True:
                        x,y,_=self.find_nearby_star(self.clicked[0],self.clicked[1])
                        if x>=0 and y>=0:                            
                            self.draw_point(x,y)

                    else:
                        self.draw_point(self.clicked[0],self.clicked[1])

                    self.draw_point_state=False
                
                if self.draw_circle_state==True and self.clickstate_draw_circle==0: 
                    """
                    if a circle shall be drawn
                    """

                    if self.snap_on_stars==True:
                        x,y,_=self.find_nearby_star(self.clicked[0],self.clicked[1])
                        if x>=0 and y>=0:
                            print("in circle snap on stars")
                            print("xy",x,y)
                            self.clickstate_draw_circle+=1                            
                            self.draw_circle(x,y,0)
                    else:
                        self.clickstate_draw_circle+=1
                        self.draw_circle(self.clicked[0],self.clicked[1],0)
                else:
                    self.clickstate_draw_circle=0
                    self.draw_circle_state=False
                
                if self.calibrate_magnitude_state==True:
                    """
                    Stars magnitude shall be calibrated
                    """
                    x,y,index=self.find_nearby_star(self.clicked[0],self.clicked[1])
                    if x>=0 and y>=0:
                        print("calibrating star")
                        mag=float(self.stringvar_calibrate_magnitude.get())
                        self.calibrate_magnitude(mag, index)

                    self.calibrate_magnitude_state=False

        
        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()

    def show_RA_DEC_orientation(self):

        self.show_orientation_state = not self.show_orientation_state

        if self.show_orientation_state==True:
            dx_ra,dy_ra=astromath.rotate_RA_DEC_vector(self.header, 1, 0)
            dx_dec,dy_dec=astromath.rotate_RA_DEC_vector(self.header, 0, 1)

            #normalize vector, since image also applies stretching
            length_ra,length_dec=np.sqrt(dx_ra**2+dy_ra**2),np.sqrt(dx_dec**2+dy_dec**2)
            dx_ra/=length_ra
            dy_ra/=length_ra 
            dx_dec/=length_dec
            dy_dec/=length_dec

            image_x=float(self.header["NAXIS1"]) #get image size x,y
            image_y=float(self.header["NAXIS2"])
            dx_ra*=image_y/10
            dy_ra*=image_y/10
            dx_dec*=image_y/10
            dy_dec*=image_y/10

            print("dx ra arrow",dx_ra)
            print("dy ra arrow",dy_ra)

            self.arrow_container.append(self.axes.arrow(image_x/10,image_y/10,dx_ra,dy_ra,color="red",alpha=0.8)) #red seems ok, green not
            self.arrow_container.append(self.axes.arrow(image_x/10,image_y/10,dx_dec,dy_dec,color="green",alpha=0.8))


        if self.show_orientation_state==False:
            while len(self.arrow_container)>0:
                self.arrow_container.remove(self.arrow_container[-1])

        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()
        

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
                ra_h,ra_min,ra_sec=astromath.decimal_rec_to_hours(RA)
                dec_deg,dec_min,dec_sec=astromath.decimal_dec_to_hours(DEC)

                self.stringvar_current_coordinates_RA.set("RA: "+str(ra_h)+"h "+str(ra_min)+"min "+str(np.round(ra_sec,1))+"sec ")
                self.stringvar_current_coordinates_DEC.set("DEC: "+str(dec_deg)+"° "+str(dec_min)+"min "+str(np.round(dec_sec,1))+"sec ")                

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


        if event.button=="up":
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

        if self.draw_line_state==True and self.clickstate_draw_line==1 and len(self.line_container)>0 and self.snap_on_stars==False:
            x_current_line=self.line_container[-1].get_xdata(orig=False)
            y_current_line=self.line_container[-1].get_ydata(orig=False)
            length_pixel=astromath.return_distance_pixel(x_current_line[0],y_current_line[0],self.current_position[0],self.current_position[1])
            length_arcsec=length_pixel*self.image_scale_arcsec_per_pixel
            self.line_container[-1].set_data([x_current_line[0],self.current_position[0]],[y_current_line[0],self.current_position[1]])
            self.text_container[-1].set_position((x_current_line[0]+(self.current_position[0]-x_current_line[0])/2,y_current_line[0]+(self.current_position[1]-y_current_line[0])/2))
            self.text_container[-1].set_text("length in arcmin "+str(np.round(length_arcsec/60,6)))


        elif self.draw_line_state==True and self.clickstate_draw_line==1 and len(self.line_container)>0 and self.snap_on_stars==True:
            x_current_line=self.line_container[-1].get_xdata(orig=False)
            y_current_line=self.line_container[-1].get_ydata(orig=False)
            x,y,_=self.find_nearby_star(self.current_position[0],self.current_position[1],tol=100000)
            if x>=0 and y>=0:                            

                length_arcsec=astromath.return_distance_pixel_scaled(x_current_line[0],y_current_line[0],x,y,self.image_scale_arcsec_per_pixel)
                self.line_container[-1].set_data([x_current_line[0],x],[y_current_line[0],y])
                self.text_container[-1].set_position((x_current_line[0]+(x-x_current_line[0])/2,y_current_line[0]+(y-y_current_line[0])/2))
                self.text_container[-1].set_text("length in arcmin "+str(np.round(length_arcsec/60,6)))


        if self.draw_circle_state==True and self.clickstate_draw_circle==1 and len(self.ellipse_container)>0:
            x,y=self.ellipse_container[-1].get_center()
            radius=astromath.return_distance_pixel(x,y,self.current_position[0],self.current_position[1])
            length_arcsec=self.image_scale_arcsec_per_pixel*2*radius
            self.ellipse_container[-1].width=2*radius
            self.ellipse_container[-1].height=2*radius
            self.text_container[-1].set_text("diameter in arcmin "+str(np.round(length_arcsec/60,6)))
            
        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()


    def draw_line(self,x,y):
        """
        adds a new line to current axes with length 0. line is updated when mouse is moved
        """
        print("drawing line")

        l=lines.Line2D([x,x],[y,y],color="magenta")
        self.axes.add_line(l)
        self.line_container.append(l)
        self.text_container.append(self.axes.text(x,y,"",color="magenta"))

    def find_draw_coordinates(self):
        """
        draws the coordinates typed in the textfield RA,Dec on the image as a circle
        
        """

        RA=float(self.stringvar_RA.get())
        DEC=float(self.stringvar_DEC.get())
        X,Y=astromath.return_X_Y_coordinates(self.header,RA,DEC)
        self.point_container.append(self.axes.scatter(X,Y,color="None",edgecolors="orangered",alpha=0.5))
        ra_h,ra_min,ra_sec=astromath.decimal_rec_to_hours(RA)
        dec_deg,dec_min,dec_sec=astromath.decimal_dec_to_hours(DEC)
        text="RA: "+str(ra_h)+"h "+str(ra_min)+"min "+str(np.round(ra_sec,1))+"sec "+"DEC: "+str(dec_deg)+"° "+str(dec_min)+"min "+str(np.round(dec_sec,1))+"sec "
        self.text_container.append(self.axes.text(X+10,Y,text,color="orangered"))
        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()



    def set_drawline(self):
        self.draw_line_state=True

    def set_drawpoint(self):
        self.draw_point_state=True
        
    def draw_point(self,x,y):
        self.point_container.append(self.axes.scatter(x,y,color="None",edgecolors="springgreen",alpha=0.5))
        RA,DEC=astromath.return_coordinates_RA_DEC(self.header,x,y)
        ra_h,ra_min,ra_sec=astromath.decimal_rec_to_hours(RA)
        dec_deg,dec_min,dec_sec=astromath.decimal_dec_to_hours(DEC)
        text="RA: "+str(ra_h)+"h "+str(ra_min)+"min "+str(np.round(ra_sec,1))+"sec "+"DEC: "+str(dec_deg)+"° "+str(dec_min)+"min "+str(np.round(dec_sec,1))+"sec "
        self.text_container.append(self.axes.text(x+10,y,text,color="springgreen"))
        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()

    def set_drawcircle(self):
        self.draw_circle_state=True

    def draw_circle(self,x,y,r):

        ellipse=Ellipse((x,y),r,r,0,fill=False,edgecolor="red")
        self.ellipse_container.append(ellipse)
        self.axes.add_patch(ellipse)
        self.text_container.append(self.axes.text(x,y,"",color="red"))
        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()

    def show_magnitude(self):
        self.show_magnitude_state= not self.show_magnitude_state
        if self.show_magnitude_state ==True:
            positions = np.transpose((self.sources['xcentroid'], self.sources['ycentroid']))
            #apertures = CircularAperture(positions, r=3.5)
            #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))

            for ii in range(len(self.sources['xcentroid'])):

                self.text_magnitude_container.append(self.axes.text(self.sources['xcentroid'][ii]+7, self.sources['ycentroid'][ii],np.round(self.sources['mag'][ii],2),color="cornflowerblue",alpha=0.5))
            #apertures.plot(color='blue', lw=1.5, alpha=0.5)
        

        if self.show_magnitude_state==False:
            while len(self.text_magnitude_container)>0:
                self.text_magnitude_container.remove(self.text_magnitude_container[-1])

        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()


    def invert_snap_on_stars(self):
        self.snap_on_stars= not self.snap_on_stars
        print("snap on stars is",self.snap_on_stars)

    def reset_all_drawings(self):


        while len(self.line_container)>0:
            print("line_container length: ",len(self.line_container))
            self.line_container.remove(self.line_container[-1])
            print("line length after remove: ",len(self.line_container))
        while len(self.text_container)>0:
            self.text_container.remove(self.text_container[-1])
        while len(self.ellipse_container)>0:
            self.ellipse_container.remove(self.ellipse_container[-1])
        while len(self.text_magnitude_container)>0:
            self.text_magnitude_container.remove(self.text_magnitude_container[-1])
        while len(self.point_container)>0:
            self.point_container.remove(self.point_container[-1])
        while len(self.arrow_container)>0:
            self.arrow_container.remove(self.arrow_container[-1])

        self.axes.clear()
        self.axes.imshow(self.data**-self.doublevar_gamma.get(),cmap="Greys",interpolation="none")

        self.Checkbox_show_magnitude.deselect()
        self.Checkbox_show_orientation.deselect()
        self.show_magnitude_state=False
        self.show_orientation_state=False

        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()

    def set_calibrate_magnitude(self):
        self.calibrate_magnitude_state=True

    def calibrate_magnitude(self,mag_calibrated_reference_star,index):
        mag_instrument_reference_star=self.sources["mag"][index]
        
        for i in range(len(self.sources["mag"])):
            mag_instrument_star=self.sources["mag"][i]
            self.sources["mag"][i]=mag_calibrated_reference_star+(mag_instrument_star-mag_instrument_reference_star)

    def update_gamma_value(self,event):
        #self.axes.clear()
        self.axes.imshow(self.data**-self.doublevar_gamma.get(),cmap="Greys",interpolation="none")

        print("new gamma value",self.doublevar_gamma.get())
        self.figure.canvas.draw()
        #self.figure.canvas.flush_events()




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


        Frame_drawing=Frame(self.Elements_Frame,borderwidth=1,relief=RIDGE)
        Label_drawing=Label(Frame_drawing,text="Drawing",font=("Arial", 12))
        Label_drawing.grid(row=0,column=0)
        Line_draw_button=Button(Frame_drawing,text="Line",command=self.set_drawline)
        Line_draw_button.grid(row=1,column=0)
        Point_draw_button=Button(Frame_drawing,text="Point",command=self.set_drawpoint)
        Point_draw_button.grid(row=1,column=1)
        Circle_draw_button=Button(Frame_drawing,text="Circle",command=self.set_drawcircle)
        Circle_draw_button.grid(row=1,column=2)
        snap_stars_checkbox=Checkbutton(Frame_drawing,text="snap on stars",onvalue=1, offvalue=0,command=self.invert_snap_on_stars)
        snap_stars_checkbox.deselect()
        snap_stars_checkbox.grid(row=2,column=1)
        #__=CreateToolTip(Frame_drawing, self.Tooltip_strings.tooltip_click_draw)
        #__=CreateToolTip(snap_stars_checkbox,self.Tooltip_strings.tooltip_snap_on_stars)



        Frame_find_coordinates=Frame(self.Elements_Frame,borderwidth=1,relief=RIDGE)
        Label_coordinates=Label(Frame_find_coordinates,text="Find Coordinates",font=("Arial", 12))
        Label_coordinates.grid(row=0,column=0)
        Label_coordinates_RA=Label(Frame_find_coordinates,text="RA in decimals:")
        Label_coordinates_RA.grid(row=1,column=0)
        Label_coordinates_DEC=Label(Frame_find_coordinates,text="DEC in decimals:")
        Label_coordinates_DEC.grid(row=2,column=0)

        Textfield_RA=Entry(Frame_find_coordinates, textvariable=self.stringvar_RA)
        Textfield_RA.grid(row=1,column=1)
        Textfield_DEC=Entry(Frame_find_coordinates, textvariable=self.stringvar_DEC)
        Textfield_DEC.grid(row=2,column=1)
        Button_Find_coord=Button(Frame_find_coordinates,text="Find entered coordinates",command=self.find_draw_coordinates)
        Button_Find_coord.grid(row=3,column=1)
        print("DEC entry in decimals is",self.stringvar_DEC.get())
        print("RA entry in decimals is",self.stringvar_RA.get())

        #currently deactivated
        self.Checkbox_show_orientation=Checkbutton(Frame_find_coordinates,text="Show axes orientation",onvalue=1,offvalue=0,command=self.show_RA_DEC_orientation)
        self.Checkbox_show_orientation.select()
        self.show_RA_DEC_orientation() #do initial plot
        self.Checkbox_show_orientation.grid(row=4,column=0)
        #__=CreateToolTip(self.Checkbox_show_orientation, self.Tooltip_strings.tooltip_show_orientation)
        

        Frame_current_mouse_position=Frame(self.Elements_Frame,borderwidth=1,relief=RIDGE)
        Label_show_coords=Label(Frame_current_mouse_position,text="Current mouse coordinates",font=("Arial", 12))
        Label_show_coords.grid(row=0,column=0,columnspan=2)
        Label_show_coords_RA=Label(Frame_current_mouse_position,textvariable=self.stringvar_current_coordinates_RA)
        Label_show_coords_RA.grid(row=1,column=0)
        Label_show_coords_DEC=Label(Frame_current_mouse_position,textvariable=self.stringvar_current_coordinates_DEC)
        Label_show_coords_DEC.grid(row=1,column=1)

        Frame_calibrate_magnitude=Frame(self.Elements_Frame,border=1,relief=RIDGE)
        Label_magnitude=Label(Frame_calibrate_magnitude,text="Calibrate Stars Magnitude",font=("Arial", 12))
        Label_magnitude.grid(row=0,column=0)
        Textfield_calibrate_magnitude=Entry(Frame_calibrate_magnitude, textvariable=self.stringvar_calibrate_magnitude)
        Textfield_calibrate_magnitude.grid(row=1,column=0)
        Button_calibrate_magnitude=Button(Frame_calibrate_magnitude,text="calibrate magnitude",command=self.set_calibrate_magnitude)
        Button_calibrate_magnitude.grid(row=1,column=1)
        self.Checkbox_show_magnitude=Checkbutton(Frame_calibrate_magnitude,text="show stars relative magnitude",onvalue=1, offvalue=0, command=self.show_magnitude)
        self.Checkbox_show_magnitude.deselect()
        self.Checkbox_show_magnitude.grid(row=2,column=0)
        #__=CreateToolTip(self.Checkbox_show_magnitude, self.Tooltip_strings.tooltip_show_stars_magnitude)
        #__=CreateToolTip(Button_calibrate_magnitude, self.Tooltip_strings.tooltip_calibrate_stars_magnitude)



        Frame_drawing.grid(row=0,column=0,sticky="NW",pady=20)
        Frame_current_mouse_position.grid(row=1,column=0,pady=20,sticky="NW")
        Frame_find_coordinates.grid(row=2,column=0,sticky="NW",pady=20)
        Frame_calibrate_magnitude.grid(row=4,column=0,sticky="NW",pady=20)


        Slider_gamma=Scale(self.Elements_Frame,from_=0.1,to=10,resolution=0.1,length=300,orient="horizontal",variable=self.doublevar_gamma, command=self.update_gamma_value)
        Slider_gamma.grid(row=5,column=0,sticky="NW")
        Label_gamma_slider=Label(self.Elements_Frame,text="Change Gamma value of image")
        Label_gamma_slider.grid(row=6,column=0,sticky="NW",pady=10)
        Button_reset=Button(self.Elements_Frame,text="reset all drawings",command=self.reset_all_drawings)
        Button_reset.grid(row=7,column=0,pady=20)

        self.Elements_Frame.grid(row=0,column=0,sticky="NW")
        self.Image_Frame.grid(row=0,column=1,sticky="NN")

        self.GUI_Frame.pack()

        
        # dont use this function otherwise program is traped in this mainloop
        #self.GUI.mainloop()