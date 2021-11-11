
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, ion, show
import math
from functools import partial
import os
from pathlib import Path

import astromath #my own module
from calibration import calibration # my own module
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

class astrometry(calibration):

    def __init__(self):
        super().__init__()
        self.perform_median_correction=True #irrelevant if dark_correction is true
        self.sigma=4.0
        self.threshold=9.0
        self.fwhm=4
        self.test_file_number=1
        self.output_folder_number=0
        self.moving_star_tolerance=0.1 #in arcsec
        self.distance_tolerance_arcsec=8
        self.distance_tolerance_pixel=3

        self.ref_i=0
        self.same_stars_searched=False
        self.moving_stars_searched=False


        self.star_list=[]
        self.sources_list=[]

        ion()
        



    
    def get_test_file_number(self):
        return self.test_file_number
    
    def set_parameters(self,fwhm,sigma,threshold,distance_tolerance_arcsec,moving_star_tolerance):
        self.fwhm=fwhm
        self.sigma=sigma
        self.threshold=threshold
        self.distance_tolerance_arcsec=distance_tolerance_arcsec
        self.moving_star_tolerance=moving_star_tolerance

    def plot_all_sources(self):
        """
        This method plots all images that have been imported with its found stars
        
        """
        if self.light_imported==False:
            messagebox.showerror("Error", "import lights first")
        
        else:
            plt.clf()
            for i,data in enumerate(self.lights):
                fi=self.lights_path[i]

                _mean, median, std = sigma_clipped_stats(data, sigma=self.sigma) #get data statistic


                daofind = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold*std)
                sources=daofind(data)


                positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
                apertures = CircularAperture(positions, r=4.)
                #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))
                plt.figure(i+1)
                plt.title(os.path.basename(fi))
                plt.imshow(np.log(data), cmap='Greys', origin='lower', interpolation='none')
                for ii in range(len(sources['xcentroid'])):
                    plt.text(sources['xcentroid'][ii]+10, sources['ycentroid'][ii],sources['id'][ii])

                apertures.plot(color='blue', lw=1.5, alpha=0.5)
            plt.show()
                
    def test_plot(self,file_number):
        """Plot only the image with the given filenumber

        Args:
            file_number (int): number of file to plot (starting with 1)
        """
        if self.light_imported==False:
            messagebox.showerror("Error", "import lights first")
        
        else:
            data=self.lights[file_number-1]

            _mean, median, std = sigma_clipped_stats(data, sigma=self.sigma) #get data statistic
            #calibration

            daofind = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold*std)
            sources=daofind(data)
            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
            apertures = CircularAperture(positions, r=4.)
            #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))
            plt.clf()
            plt.figure(1)
            plt.title("File number:"+str(file_number)+", found stars: "+str(len(sources['id'])))
            plt.imshow(np.log(data), cmap='Greys', origin='lower', interpolation='none')
            for ii in range(len(sources['xcentroid'])):
                plt.text(sources['xcentroid'][ii]+10, sources['ycentroid'][ii],sources['id'][ii])
            apertures.plot(color='blue', lw=1.5, alpha=0.5)
            plt.show()


    def find_sources(self,save_files=True): 
        """searches for stars in the loaded images. The sources are searched via DAOStarFinder.
           stars are saved in sources_list. Beforehands a calibration is done if activated.
           if save_files is True then the stars with its properties are also saved to a file.

        Args:
            save_files (bool): Decides if found stars are also written to file. Defaults to True.

        Returns:
            sources_list [list]: list of found stars. structure is sources_list[filenumber][starnumber]
        """
        if self.light_imported==False:
            messagebox.showerror("Error", "import lights first")
        
        else:
            sources_list=[]
            _found=False
            while _found==False: #define filename for output and find folder number
                try:
                    os.mkdir(os.path.dirname(self.lights_path[0])+"/data_stars"+str(self.output_folder_number))
                    _found=True
                except FileExistsError:
                    self.output_folder_number+=1

            for i,data in enumerate(self.lights):

                _mean, median, std = sigma_clipped_stats(data, sigma=self.sigma) #get data statistic



                ####find stars in file################

                daofind = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold*std)
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
                    

                RA,DEC =astromath.return_coordinates_RA_DEC(self.headers[i], sources['xcentroid'], sources['ycentroid'])
                sources['RA']=RA
                sources['DEC']=DEC
                sources.sort('id')
                del sources['sky']; del sources['roundness1']
                del sources['roundness2']; del sources['npix']; del sources['sharpness']

                sources_list.append(sources) #write found stars of file into list
                
                if save_files==True:
            
                    sources.write(os.path.dirname(self.lights_path[i])+"/data_stars"+str(self.output_folder_number)+"/"+Path(self.lights_path[i]).stem+".dat",format='ascii',overwrite=True)
                
            self.sources_list=sources_list
                    




    def search_for_moving_stars(self):
        """
        currently under development. not finished yet.

        """
        if self.light_imported==False:
            messagebox.showerror("Error", "import lights first")
        
        else:
                
            if self.same_stars_searched==False:
                self.search_find()
            print("in search for moving stars")
            #star_list[star][file] type: star.star object
            ref_star_ID=-1 #error if not found

            #star_list[found stars in reference image][filenumber of related stars]  type=object star
            for i, star in enumerate(self.star_list): 
                print("i",i)

                dev_dec=[]
                dev_ra=[]
                rms=[]

                for s in star:  #browse through files in which star was found
                    dev_dec.append(s.get_dev_dec()*3600)
                    dev_ra.append(s.get_dev_ra()*3600)
                    rms.append(np.sqrt(s.get_dev_dec()**2+s.get_dev_ra()**2)*3600)
                    if s.get_dev_dec()==0.0:
                        ref_star_ID=s.get_ID()
                try:
                    m_dec,_=np.polyfit(range(len(dev_dec)),dev_dec,1)
                    m_ra,_ =np.polyfit(range(len(dev_ra)),dev_ra,1)
                    m_rms,_=np.polyfit(range(len(rms)),rms,1)
                except:
                    m_dec=0
                    m_ra=0
                    m_rms=0

                print("m_dec ", m_dec)
                print("m_ra ",m_ra)
                if (i<=2 or m_rms>self.moving_star_tolerance):
                    plt.figure(i) #too many
                    plt.xlabel("dev RA in arcsec")
                    plt.ylabel("dev DEC in arcsec")
                    plt.title("Reference image: "+str(self.ref_i)+" star id: "+str(ref_star_ID)+" m_dec "+str(np.round(m_dec,4))+" m_ra"+str(np.round(m_ra,4)))
                    plt.plot(dev_ra,dev_dec,linestyle='-',marker='+')
            self.test_plot(self.ref_i)
            plt.show()





    def search_same_stars(self):
        """
        This method searches for the same stars in different files and relates them.
        
        """
        if self.light_imported==False:
            messagebox.showerror("Error", "import lights first")
        
        else:
                
            self.star_list=[] #reset list
            # _number=0
            # _found=False

            # while _found==False: #define filename for output
            #     try:
            #         os.mkdir("data_stars"+str(_number))
            #         _found=True
            #     except FileExistsError:
            #         _number+=1
            
            #self.ref_i=math.floor(len(sources_list)/2.0) #middle entry 
            self.ref_i=astromath.get_image_with_highest_index(self.sources_list) # entry with highest number of stars
            print("reference image is image file number:",self.ref_i)

            for i in range(len(self.sources_list[self.ref_i]['id'])): #reference image, loop over stars i in reference image

                print("at star number",i)
                ra_ref=self.sources_list[self.ref_i]['RA'][i] #RA of reference star nr. i
                dec_ref=self.sources_list[self.ref_i]['DEC'][i] #DEC of reference star nr. i
                flux_ref=self.sources_list[self.ref_i]['flux'][i] # flux of reference star nr. i
                mag_ref=self.sources_list[self.ref_i]['mag'][i]
                id_ref=self.sources_list[self.ref_i]['id'][i]
                y_ref=self.sources_list[self.ref_i]['ycentroid'][i]
                x_ref=self.sources_list[self.ref_i]['xcentroid'][i]
                file_ref=os.path.basename(self.lights_path[self.ref_i])
                ref_star=star.star(star_ID=id_ref,ra=ra_ref,dec=dec_ref,mag=mag_ref,flux=flux_ref,x=x_ref,y=y_ref, filename=file_ref,dev_dec=0, dev_ra=0)
                star_files_list=[] #every reference star has an entry in here. All same stars from other files are added here
                



                for j, star_file in enumerate(self.sources_list): #loop over files j, exclude reference image

                    
                    if j!=self.ref_i:
                        for k in range(len(star_file['id'])): #loop over stars k in file j, !exclude reference image! maybe exclude distance=0 for that?
                        
                            #print("star number",k)
                            ra=star_file['RA'][k]
                            dec=star_file['DEC'][k]
                            x=star_file['xcentroid'][k]
                            y=star_file['ycentroid'][k]



                            if astromath.return_distance_arsec(ra,dec,ra_ref,dec_ref)<=self.distance_tolerance_arcsec:
                                #print("detected")
                                s=star.star(star_ID=int(star_file['id'][k]),ra=ra, dec=dec, flux=star_file['flux'][k],mag=star_file['mag'][k],filename=os.path.basename(self.lights_path[j]), dev_dec=dec-dec_ref,dev_ra=ra-ra_ref)
                                star_files_list.append(s)  #if star is found in file j, add to list, next file then
                    else:
                        star_files_list.append(ref_star) #so that reference frame is at correct position
                            


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
                    filename_help.append(s.filename) # I dont know yet how to get this to file
                    dev_ra_help.append(s.dev_ra*3600)
                    dev_dec_help.append(s.dev_dec*3600)
                    rms_help.append(np.sqrt(s.dev_dec**2+s.dev_ra**2)*3600)


                header="  ID    	RA/deg	      DEC/deg	     flux	     mag 	  dev_ra in arces  dev_dec in arcsec  RMS in arcsec"
                #this is done for every star i  (s) 
                #save.save_to_file_5D(id_help,ra_help,dec_help, flux_help,mag_help, filename="data_stars"+str(_number)+"/id_"+str(sources_list[0]['id'][i]), acuracy=6,header=header)
                output=np.array([id_help,ra_help,dec_help, flux_help,mag_help,dev_ra_help, dev_dec_help, rms_help])
                #errors can occur here due to index overflow!
                save.save_to_file_ND(output, filename=os.path.dirname(self.lights_path[0])+"/data_stars"+str(self.output_folder_number)+"/id_"+str(self.sources_list[self.ref_i]['id'][i]), acuracy=6,header=header)
                        
                self.star_list.append(star_files_list)
            print("done finding same stars")

    def search_find(self):
        if self.light_imported==False:
            messagebox.showerror("Error", "import lights first")
        
        else:
            self.find_sources(False)
            self.search_same_stars()
