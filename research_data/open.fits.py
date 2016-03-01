##=============== OPENING FITS TABLES =============
import numpy as np
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt


##=============== settting environment ============
fermi = fits.open('3C454.3_604800.lc')

data = fermi[1].data

print data
