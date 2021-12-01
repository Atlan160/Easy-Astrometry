import spectroscopy
import numpy as np
import matplotlib.pyplot as plt
import math



spec=spectroscopy.spectroscopy()
spec.set_light_path("spectra/aldebaran_neu.fits")
spec.show_spectra()
spec.plot_spectra()