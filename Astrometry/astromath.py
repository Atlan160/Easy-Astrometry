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


def return_coordinates_RA_DEC(header, x,y):

    try:
        print("image has a renewed plate solve")
        c11=header['CD1_1']
        c12=header['CD1_2']
        c21=header['CD2_1']
        c22=header['CD2_2']
    except:
        print("original plate solve")
        cdelt1=header['CDELT1']
        cdelt2=header['CDELT2']
        crota2=header['CROTA1']
        c11=cdelt1*math.cos(crota2)
        c12=-cdelt2*math.sin(crota2)
        c21=cdelt1*math.sin(crota2)
        c22=cdelt2*math.sin(crota2)
    
    try:
        xref=header['CRPIX1']
        yref=header['CRPIX2']
        RAref=header['CRVAL1']
        DECref=header['CRVAL2']
    except:
        print("something went wrong with extracting reference points")
        return (0,0)

    RAnew=np.zeros_like(x)
    DECnew=np.zeros_like(y)

    for i in range(len(RAnew)):

        RAnew[i]=c11*(x[i]-xref)+c12*(y[i]-yref)+RAref
        DECnew[i]=c21*(x[i]-xref)+c22*(y[i]-yref)+DECref

    return (RAnew,DECnew)
    



def return_distance_arsec(ra1, dec1, ra2,dec2):
    return math.sqrt((ra1-ra2)**2+(dec1-dec2)**2)*3600

def return_pixel_distance(x1,y1,x2,y2):
    return math.sqrt((x1-x2)**2+(y1-y2)**2)