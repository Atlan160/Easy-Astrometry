
import numpy as np
import matplotlib.pyplot as plt
import math
from functools import partial
import os
from pathlib import Path

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

# TODO implement dark correction correctly
# TODO create new folder when data_stars folder doesnt exist
# TODO take middle light as reference frame (regarding time), so that moving targets can be easier spottet


plot_all=True
perform_dark_correction=False
perform_median_correction=True #irrelevant if dark_correction is true
sigma=4.0
threshold=6.0
fwhm=4.5
moving_star_tolerance=0.6 #in arcsec
distance_tolerance_arcsec=5
distance_tolerance_pixel=100
dark_path=" " #is imported 
lights_path=" " #is imported

star_list=[]
sources_list=[]

def dark_correction(data,dark):
    return data-dark

def median_correction(data,median):
    return data-median

def plot_all_sources():

    for i,fi in enumerate(lights_path):
        light=fits.open(fi)
        data=light[0].data
        _mean, median, std = sigma_clipped_stats(data, sigma=sigma) #get data statistic

        #calibration
        if perform_dark_correction:
            dark=fits.open(dark_path)
            dark_data=dark[0].data
            data=dark_correction(data,dark_data)
        else:
             data=median_correction(data,median)

        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std)
        sources=daofind(data)


        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))
        plt.figure(i)
        plt.title(os.path.basename(fi))
        plt.imshow(np.log(data), cmap='Greys', origin='lower', interpolation='none')
        for ii in range(len(sources['xcentroid'])):
            plt.text(sources['xcentroid'][ii]+10, sources['ycentroid'][ii],sources['id'][ii])

        apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()
            
def test_plot():


    light=fits.open(lights_path[0])
    data=light[0].data
    _mean, median, std = sigma_clipped_stats(data, sigma=sigma) #get data statistic
    #calibration
    if perform_dark_correction:
        dark=fits.open(dark_path)
        dark_data=dark[0].data
        data=dark_correction(data,dark_data)
    else:
         data=median_correction(data,median)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std)
    sources=daofind(data)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.)
    #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))
    plt.figure(1)
    plt.title("first lightfile, found stars: "+str(len(sources['id'])))
    plt.imshow(np.log(data), cmap='Greys', origin='lower', interpolation='none')
    for ii in range(len(sources['xcentroid'])):
        plt.text(sources['xcentroid'][ii]+10, sources['ycentroid'][ii],sources['id'][ii])
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()


def find_sources(save_files=True):
    sources_list=[]
    _found=False
    _number=0
    while _found==False: #define filename for output
        try:
            os.mkdir(os.path.dirname(lights_path[0])+"/data_stars"+str(_number))
            _found=True
        except FileExistsError:
            _number+=1

    for fi in lights_path:
        #print("os.path.dirname(fi)", os.path.dirname(fi))


        #####import files  ############
        light=fits.open(fi)
        data = light[0].data
        _mean, median, std = sigma_clipped_stats(data, sigma=sigma) #get data statistic

        #calibration
        if perform_dark_correction:
            dark=fits.open(dark)
            dark_data=dark[0].data
            data=dark_correction(data,dark_data)
        else:
            data=median_correction(data,median)



        ####find stars in file################

        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std)
        sources=daofind(data) #sources is a astropy.table type, sources contains all the found stars by daofind algorithm
        # print("Found "+str(len(sources['id']))+" Stars in file number "+str(i))
        # if plot_all==True:
        #     positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        #     apertures = CircularAperture(positions, r=4.)
        #     #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))
        #     plt.figure(i)
        #     plt.imshow(np.log(data), cmap='Greys', origin='lower', interpolation='none')
        #     for ii in range(len(sources['xcentroid'])):
        #         plt.text(sources['xcentroid'][ii]+10, sources['ycentroid'][ii],sources['id'][ii])

        #     apertures.plot(color='blue', lw=1.5, alpha=0.5)
        #     plt.show()
            

        RA,DEC =astromath.return_coordinates_RA_DEC(light[0].header, sources['xcentroid'], sources['ycentroid'])
        sources['RA']=RA
        sources['DEC']=DEC
        sources.sort('id')
        del sources['sky']; del sources['roundness1']
        del sources['roundness2']; del sources['npix']; del sources['sharpness']

        sources_list.append(sources) #write into an array
        
        if save_files==True:
      
            sources.write(os.path.dirname(fi)+"/data_stars"+str(_number)+"/"+Path(fi).stem+".dat",format='ascii',overwrite=True)
        
    return sources_list
            




    ###preprocessing
    # x_ref=sources_list[0]['xcentroid'][0]
    # y_ref=sources_list[0]['ycentroid'][0]
    # x=sources_list[1]['xcentroid'][0]
    # y=sources_list[1]['ycentroid'][0]

    # print("distance",astromath.return_pixel_distance(x_ref,y_ref,x,y))

    # print("done")
    #find same stars

def search_for_moving_stars():
    print("in search for moving stars")
    #star_list[star][file] type: star.star object
    global star_list

    for i, star in enumerate(star_list): #star i
        print("i",i)

        dev_dec=[]
        dev_ra=[]
        for s in star:  
            dev_dec.append(s.get_dev_dec()*3600)
            dev_ra.append(s.get_dev_ra()*3600)
        try:
            m_dec,_=np.polyfit(range(len(dev_dec)),dev_dec,1)
            m_ra,_ =np.polyfit(range(len(dev_ra)),dev_ra,1)
        except:
            m_dec=0
            m_ra=0

        print("m_dec ", m_dec)
        print("m_ra ",m_ra)
        if i<=5 or m_dec>moving_star_tolerance or m_ra>moving_star_tolerance:
            plt.figure(i) #too many
            plt.xlabel("dev RA in arcsec")
            plt.ylabel("dev DEC in arcsec")
            plt.title("m_dec "+str(np.round(m_dec,4))+" m_ra"+str(np.round(m_ra,4)))
            plt.plot(dev_ra,dev_dec,'r+')

    plt.show()





def search_same_stars(sources_list):
    global star_list

    star_list=[] #reset list
    _number=0
    _found=False

    while _found==False: #define filename for output
        try:
            os.mkdir("data_stars"+str(_number))
            _found=True
        except FileExistsError:
            _number+=1
    
    ref_i=math.floor(len(sources_list)/2.0) #last entry 

    for i in range(len(sources_list[ref_i]['id'])): #reference image, loop over stars i in reference image

        print("at star number",i)
        ra_ref=sources_list[ref_i]['RA'][i] #RA of reference star nr. i
        dec_ref=sources_list[ref_i]['DEC'][i] #DEC of reference star nr. i
        flux_ref=sources_list[ref_i]['flux'][i] # flux of reference star nr. i
        mag_ref=sources_list[ref_i]['mag'][i]
        id_ref=sources_list[ref_i]['id'][i]
        y_ref=sources_list[ref_i]['ycentroid'][i]
        x_ref=sources_list[ref_i]['xcentroid'][i]
        file_ref=os.path.basename(lights_path[ref_i])
        ref_star=star.star(star_ID=id_ref,ra=ra_ref,dec=dec_ref,mag=mag_ref,flux=flux_ref,x=x_ref,y=y_ref, filename=file_ref,dev_dec=0, dev_ra=0)
        star_files_list=[] #every reference star has an entry in here. All same stars from other files are added here
        star_files_list.append(ref_star)



        for j, star_file in enumerate(sources_list): #loop over files j, exclude reference image

            
            if j!=ref_i:
                for k in range(len(star_file['id'])): #loop over stars k in file j, !exclude reference image! maybe exclude distance=0 for that?
                
                    #print("star number",k)
                    ra=star_file['RA'][k]
                    dec=star_file['DEC'][k]
                    x=star_file['xcentroid'][k]
                    y=star_file['ycentroid'][k]

                    if astromath.return_pixel_distance(x,y,x_ref,y_ref)<=distance_tolerance_pixel:
                        #print("detected")
                        s=star.star(star_ID=int(star_file['id'][k]),x=x,y=y,ra=ra, dec=dec, flux=star_file['flux'][k],mag=star_file['mag'][k],filename=os.path.basename(lights_path[j]), dev_dec=dec-dec_ref,dev_ra=ra-ra_ref)
                        star_files_list.append(s)  #if star is found in file j, add to list, next file then
              

                    #if astromath.return_distance_arsec(ra,dec,ra_ref,dec_ref)<=distance_tolerance_arcsec:
                    #    #print("detected")
                    #    s=star.star(star_ID=int(star_file['id'][k]),ra=ra, dec=dec, flux=star_file['flux'][k],mag=star_file['mag'][k],filename=os.path.basename(lights_path[j]), dev_dec=dec-dec_ref,dev_ra=ra-ra_ref)
                    #    star_files_list.append(s)  #if star is found in file j, add to list, next file then
                    

        # now print star_help_list to file
        ra_help=[]
        dec_help=[]
        flux_help= []
        id_help=[]
        mag_help=[]
        filename_help=[]
        dev_ra_help=[]
        dev_dec_help=[]
        rms_help=[]
        x_help=[]
        y_help=[]


        for s in star_files_list: #better call s.get()
            ra_help.append(s.ra)
            dec_help.append(s.dec)
            flux_help.append(s.flux)
            id_help.append(s.star_ID)
            mag_help.append(s.mag)
            x_help.append(s.xpos)
            y_help.append(s.ypos)
            filename_help.append(s.filename) # I dont know yet how to get this to file
            dev_ra_help.append(s.dev_ra*3600)
            dev_dec_help.append(s.dev_dec*3600)
            rms_help.append(np.sqrt(s.dev_dec**2+s.dev_ra**2)*3600)


        header="  ID    	x	      y	     flux	     mag 	  dev_ra in arces  dev_dec in arcsec  RMS in arcsec "
        #this is done for every star i  (s) 
        #save.save_to_file_5D(id_help,ra_help,dec_help, flux_help,mag_help, filename="data_stars"+str(_number)+"/id_"+str(sources_list[0]['id'][i]), acuracy=6,header=header)
        output=np.array([id_help,x_help,y_help, flux_help,mag_help,dev_ra_help, dev_dec_help, rms_help])
        save.save_to_file_ND(output, filename="data_stars"+str(_number)+"/id_"+str(sources_list[ref_i]['id'][i]), acuracy=6,header=header)
                
        star_list.append(star_files_list)

def search_find():
    help1=find_sources(False)
    search_same_stars(help1)

def open_lights():
    global lights_path
    try:
        lights_path=filedialog.askopenfilenames(initialdir =" ", title = "Select light files",filetypes = (("fits files","*.fits"),("newly solved files",".new")))
        #print(lights_path)
        messagebox.showinfo("success","imported "+str(len(lights_path))+" files")
        #root.mainloop()
    except:
        messagebox.showerror("Error", "something went wrong importing light files")
        #root.mainloop()
def open_dark():
    global dark_path
    try:
        dark_path=filedialog.askopenfile(initialdir= " ", title ="Select dark file", filetypes = ("fits file","*.fits") )
        perform_dark_correction=True
        if dark_path == " ":
            perform_dark_correction=False
            messagebox.showinfo(" ","imported no files, median correction will then be done")
        else:
            messagebox.showinfo("success","imported dark files")
    except:
        perform_dark_correction=False
        messagebox.showerror("Error", "something went wrong importing light files, median correction will then be done")

def set_parameters(_fwhm, _sigma, _threshold, _distance_tolerance_arcsec): 
    global fwhm
    global sigma
    global threshold
    global distance_tolerance_arcsec

    fwhm=float(_fwhm.get())
    sigma=float(_sigma.get())
    threshold=float(_threshold.get())
    distance_tolerance_arcsec=float(_distance_tolerance_arcsec.get())

    

def set_settings():
    t_fwhm=StringVar()
    t_fwhm.set(str(fwhm))
    t_sigma=StringVar()
    t_sigma.set(str(sigma))
    t_thr=StringVar()
    t_thr.set(str(threshold))
    t_dist_tolerance=StringVar()
    t_dist_tolerance.set(str(distance_tolerance_arcsec))
    
    _set=partial(set_parameters,t_fwhm,t_sigma,t_thr,t_dist_tolerance)


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

    Button(root,text="test settings", command=test_plot).grid(row=4,column=0)
    Button(root, text="set settings",command=_set).grid(row=4,column=1)

    __=tooltip.CreateToolTip(lfwhm, tooltip.tooltip_fwhm)
    __=tooltip.CreateToolTip(lsigma, tooltip.tooltip_sigma)
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
runmenu.add_command(label="set parameters", command=set_settings)
runmenu.add_command(label="plot all imported lights", command=plot_all_sources)
runmenu.add_command(label="search for stars and write to file", command=find_sources)
runmenu.add_command(label="search for same stars and write to file",command=search_find)
runmenu.add_command(label="search for moving stars", command=search_for_moving_stars)
menubar.add_cascade(label="Run", menu=runmenu)

# helpmenu = Menu(menubar, tearoff=0)
# helpmenu.add_command(label="About", command=hello)
# menubar.add_cascade(label="Help", menu=helpmenu)
root.config(menu=menubar)
root.mainloop()