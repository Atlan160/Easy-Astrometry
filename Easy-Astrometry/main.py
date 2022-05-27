from tkinter import *
from tkinter import filedialog
from tkinter import messagebox

from my_modules.astrometry_root import astrometry_root
from my_modules.image_identifier_root import image_identifier_root


easy_astronomy_root=Tk(className=" Easy-Astronomy")
width, height = easy_astronomy_root.winfo_screenwidth(), easy_astronomy_root.winfo_screenheight()
easy_astronomy_root.geometry("200x200+100+100")
easy_astronomy_root.iconbitmap('_images/icons/astrometry_favicon.ico')


def start_astrometry():
    _astrometry_class=astrometry_root()

def start_spectroscopy():
    pass

def start_image_identifier():
    _image_identifier_class=image_identifier_root()


Button(easy_astronomy_root, text="Easy Astrometry",command=start_astrometry).grid(row=1,column=1)
Button(easy_astronomy_root, text="Easy Spectroscopy",command=start_spectroscopy).grid(row=2,column=1)
Button(easy_astronomy_root, text="Easy Image Identifier",command=start_image_identifier).grid(row=3,column=1)



easy_astronomy_root.mainloop()

