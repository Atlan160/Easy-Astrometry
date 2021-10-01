import math

class star():
    "This is a star"
    def __init__(self,mag=0, x=0, y=0, star_ID=-1, filename="",obs_date=0, ra=0, dec=0, flux=0, dev_x=0, dev_y=0, dev_ra=0, dev_dec=0, dev_mag=0):

        self.xpos=x
        self.ypos=y
        self.mag=mag
        self.ra=ra
        self.dec=dec
        self.flux=flux
        self.dev_x=dev_x
        self.dev_y=dev_y
        self.dev_ra=dev_ra
        self.dev_dec=dev_dec
        self.dev_mag=dev_mag
        self.star_ID=star_ID
        self.filename=filename
        self.obs_date=obs_date
        
    
    def get_mag(self):
        return float(self.mag)
    
    def get_xpos(self):
        return float(self.xpos)

    def get_ypos(self):
        return float(self.ypos)
    def get_ra(self):
        return float(self.ra)
    def get_dec(self):
        return float(self.dec)
    def get_flux(self):
        return float(self.flux)
    def get_dev_x(self):
        return float(self.dev_x)
    def get_dev_y(self):
        return float(self.dev_y)
    def get_dev_ra(self):
        return float(self.dev_ra)
    def get_dev_dec(self):
        return float(self.dev_dec)
    def get_dev_mag(self):
        return float(self.dev_mag)   
    def get_obs_date(self):
        return float(self.obs_date)
    def get_ID(self):
        return float(self.star_ID)
    def get_filename(self):
        return self.filename
    def set_filename(self,_):
        self.filename=_
    def set_ID(self,_):
        self.star_ID=_
    
    def distance(self,xposition, yposition):
        return math.sqrt((float(self.xpos)-xposition)**2+(float(self.ypos)-yposition)**2)

    def set_mag(self,mag):
        self.magnitude=mag


#def distance(s1,s2):
#    return (float(s1.get_xpos())- float(s2.get_xpos()))+(float(s1.get_ypos())-float(s2.get_ypos()))