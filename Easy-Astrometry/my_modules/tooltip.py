import tkinter as tk

class CreateToolTip(object):
    '''
    create a tooltip for a given widget
    '''

    def __init__(self, widget, text='widget info'):
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.close)


    def enter(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background='yellow', relief='solid', borderwidth=1,
                       font=("times", "10", "normal"))
        label.pack(ipadx=1)

    def close(self, event=None):
        if self.tw:
            self.tw.destroy()


class _Tooltip_strings():

    def __init__(self):
      
        self.tooltip_fwhm="adjust the size of the stars in pixel. Defocussed stars need bigger values. 3-4 is a good value"
        self.tooltip_sigma="clipping value to estimate background. Definitely has to be lower than threshold. 3-4 is a good value"
        self.tooltip_tolerance="tolerance stars can move without being falsly detected in arcsecond"

        self.tooltip_show_orientation="Is supposed to show the orientation of RA,DEC axis. \n This is not properly working yet, the orientation seems not to be correct"
        self.tooltip_click_draw="click on one of the buttons and then click on a point in the image. you need to click on the button again for a new drawing"
        self.tooltip_snap_on_stars="when this checkbox is activated, the programm will automatically check for stars nearby and then snap the point onto them. \n \
            If this takes longer than some seconds and returns an error, then you need to click closer to a recognized star."
        self.tooltip_show_stars_magnitude="show the instrument magnitude of the recognized stars. You need to calibrate with a reference star first to get correct values."
        self.tooltip_calibrate_stars_magnitude="enter the magnitude of a star you know, then click on the button and the star which belongs to the entered magnitude.\n \
             If the instrument magnitude is already shown, you have to redraw everything."