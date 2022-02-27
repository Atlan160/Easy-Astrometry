import photutils
import astropy
import numpy as np
import math
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
    """[summary]

    Args:
        header (header): the header of the current fitsfile
        x ([type]): x coordinates which should be transformed
        y ([type]): y coordinates which should be transformed

    Returns:
        touple of size 2: first entry is ra, second is dec
    """

    try:
        #print("image has a renewed plate solve")
        c11=header['CD1_1']
        c12=header['CD1_2']
        c21=header['CD2_1']
        c22=header['CD2_2']
    except:
        #print("original plate solve")
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

    if (type(x)==float and type(y)==float) or (type(x)==np.float64 and type(y)==np.float64):
        RAnew=c11*(x-xref)+c12*(y-yref)+RAref
        DECnew=c21*(x-xref)+c22*(y-yref)+DECref

        return (RAnew,DECnew)


    elif (type(x)==list and type(y)==list) or (type(x)==np.ndarray and type(y)==np.ndarray):
        RAnew=np.zeros_like(x)
        DECnew=np.zeros_like(y)

        for i in range(len(RAnew)): #ranew and decnew are arrays with length of number of stars

            RAnew[i]=c11*(x[i]-xref)+c12*(y[i]-yref)+RAref
            DECnew[i]=c21*(x[i]-xref)+c22*(y[i]-yref)+DECref

        return (RAnew,DECnew)
    else:
        print("in else statement")
        print("type of RA,Dec",type(x))
        return (0,0)

def return_X_Y_coordinates(header,RA,DEC):
    """[summary]

    Args:
        header (header): the header of the current fitsfile
        RA ([type]): RA coordinates which should be transformed
        DEC ([type]): DEC coordinates which should be transformed

    Returns:
        touple of size 2: first entry is x, second is y
    """
    try:
        #print("image has a renewed plate solve")
        c11=header['CD1_1']
        c12=header['CD1_2']
        c21=header['CD2_1']
        c22=header['CD2_2']
    except:
        #print("original plate solve")
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


    if (type(RA)==float and type(DEC)==float) or (type(RA)==np.float64 and type(DEC)==np.float64) or ((type(RA)==int and type(DEC)==int)):
        print("in single number")

        determinante=(c11*c22-c21*c12)

        Xnew=(RA-RAref)*c22/determinante-(DEC-DECref)*c12/determinante+xref
        Ynew=-(RA-RAref)*c21/determinante+(DEC-DECref)*c11/determinante+yref


        return (Xnew,Ynew)
     
    elif (type(RA)==list and type(DEC)==list) or (type(RA)==np.ndarray and type(DEC)==np.ndarray):
        Xnew=np.zeros_like(RA)
        Ynew=np.zeros_like(DEC)

        determinante=(c11*c22-c21*c12)
        for i in range(len(Xnew)): #Xnew and Ynew are arrays with length of number of stars
            Xnew[i]=(RA[i]-RAref)*c22/determinante-(DEC[i]-DECref)*c12/determinante+xref
            Ynew[i]=-(RA[i]-RAref)*c21/determinante+(DEC[i]-DECref)*c11/determinante+yref


        return (Xnew,Ynew)
    else:
        print("in else statement")
        print("type of RA,Dec",type(RA))
        return(0.0,0.0)
        



def return_distance_arsec(ra1, dec1, ra2,dec2):
    return math.sqrt((ra1-ra2)**2+(dec1-dec2)**2)*3600

def return_distance_pixel(x1,y1,x2,y2):
    return math.sqrt((x1-x2)**2+(y1-y2)**2)


def return_distance_pixel_scaled(x1,y1,x2,y2,scale):
    return scale*math.sqrt((x1-x2)**2+(y1-y2)**2)


def get_image_with_highest_index(sources):
    max=0
    index=0
    for i in range(len(sources)):
        number=len(sources[i]['id'])
        if number>max:
            index=i
    
    return index