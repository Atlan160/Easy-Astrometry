import photutils
import astropy
import numpy as np
import math
import save
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import PercentileInterval
from astropy.visualization import simple_norm
from astropy.table import Table
import astroquery
from astroquery.astrometry_net import AstrometryNet

# ast = AstrometryNet()
# ast.api_key = 'dgbjzdkqcoeynxxq'

#dark=fits.open("master_dark.fits")
plot_all=True

def dark_correction(data,dark):
    return data-dark

def median_correction(data,median):
    return data-median
number_of_files=1
for i in range(number_of_files):
    #####import files and plot the data ############
    light=fits.open("Light_5_secs_002.new")
    data = light[0].data
    save.save_text(np.str(light[0].header),filename="header_test")
    #dark_data=dark[0].data
    norm=simple_norm(data, stretch='sqrt')
    mean, median, std = sigma_clipped_stats(data, sigma=4.0) #get data statistic
    print(data.shape)
    #print(data) 

    if plot_all==True:
        plt.subplot(121)
        plt.title("wo median clipping")
        plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='none')
        
        plt.subplot(122)
        plt.title("with median clipping")
        plt.imshow(median_correction(data,median), cmap='Greys', origin='lower', norm=norm, interpolation='none')
        plt.show() 

    data=median_correction(data,median)


    ####find stars in file################

    #print(mean, median, std)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources=daofind(data) #sources is a astropy.table type

    ############################
    #print("Number of stars:", len(sources['id']))

    # ##print data of stars###
    # for col in sources.colnames: 
    #     sources[col].info.format = '%.7g'
    # print(sources)
    #print("sources colnames",sources.colnames)

    if plot_all==True:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        #norm = ImageNormalize(stretch=SqrtStretch()+PercentileInterval(70.))
        plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='none')

        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.show()

    
    save.save_text(np.str(sources.colnames)+"\n"+np.str(sources.as_array()),filename='data'+str(i+1))

    # f=open('data'+str(i+1)+'.txt', 'w')
    # f.write(np.str(sources.as_array()))
    # f.close()


#ast.api_key = 'XXXXXXXX'

# sources.sort('FLUX')
# # Reverse to get descending order
# sources.reverse()

# image_width=np.size(data,axis=1)
# image_height=np.size(data,axis=0)
# print("widht, height", image_height, image_height)
# print("solving")
# wcs_header=ast.solve_from_source_list(sources['X_IMAGE'], sources['Y_IMAGE'], image_width, image_height, solve_timeout=120)
# print("wcs_header")
# print("done")

