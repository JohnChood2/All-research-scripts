##=============== OPENING FITS TABLES =============
import matplotlib
matplotlib.use( 'Agg' )
import numpy as np
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
import math
import datetime as dt


##=============== settting environment ============


fermi = fits.open('PMNJ2345-1555_86400.lc')

tbdata = fermi[1].data

#print tbdata[0]

Var_Names = tbdata.names
print Var_Names

##def MET_to_MJD(x_axe):
##    #convert MET to days since 2001
##    mjd = ((x_axe)/(3.15e7))
##
##    
##    return mjd 
##mjd = ((x_axe)/(3.15e7))

##===== using data ==========
x_axis = 'START'
y_axis = 'FLUX_1000_300000'
y_axis_err = 'ERROR_1000_300000'
UL = 'UL_1000_300000'
Plot_color = 'k'

x_axe = tbdata[x_axis]
y_axe = tbdata[y_axis]
y_axe_err = tbdata[y_axis_err]
UL_axes   = tbdata[UL]
print len(y_axe)

#mod_x = MET_to_MJD(x_axe)

def MET_to_MJD(x_axe):
    #convert MET to days since 2001
    MJD = ((((x_axe)/(3.15e7))*(365))+ (51910.00074287037))
   
    return MJD

MJD = ((((x_axe)/(3.15e7))*(365))+ (51910.00074287037))
x_axe = MJD 
print type(x_axe)
print type(MJD)
#print MJD
#print UL_axes

##============ data filter =================
if y_axis == 'FLUX_1000_300000':
    y_id = np.where(y_axe >= 0)[0]
    x_axe_c = x_axe        [y_id]
    y_axe_c = y_axe        [y_id]
    y_axe_err_c = y_axe_err[y_id]
   
print len(x_axe),len(x_axe_c)
print len(y_id), max(y_id), min(y_id)
IDX = np.where(UL_axes == 'F')[0]
print len(IDX), max(IDX)
#print x_axe[2672]

t_N = x_axe[IDX]
F_N = y_axe[IDX]
E_N = y_axe_err[IDX]
print F_N
arrow_dat = np.where(UL_axes == 'T')[0]
# Directions
u = np.zeros(arrow_dat.size)
v = np.ones(arrow_dat.size)*-500.

##============ plotting =====================

##plt.plot(x_axis, y_axis)
##plt.xlabel('MJD (d)')
##plt.ylabel('Flux_300_1000')
##plt.title('Source 3C 454.3')
##plt.show()

fig = plt.figure()
ax = fig.add_subplot( 111, axisbg='white')
ax.set_xlabel('MJD (d)', fontsize=18)
ax.set_ylabel('FLUX_1000_300000', fontsize=18)
ax.set_ylim(0, .0000005)
ax.set_title( 'Source = PMNJ2345-1555 DURATION = 86400')
ax.errorbar(t_N, F_N, yerr=E_N, 
    color=Plot_color, markersize=3, fmt = 'o', linewidth=1, 
    ecolor=Plot_color, elinewidth=1, capthick = 1 
    )#, label=r'\boldmath $\overline{\mathrm{Cens}}$' )
ax.quiver(x_axe[arrow_dat], y_axe[arrow_dat],u,v, color='r', zorder=5) 
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.tight_layout()
fname = '#sPMNJ2345-1555.png'#.format(x_axis, y_axis)
plt.savefig( fname, bbox_inches='tight' )
print( '    Figure: {0}'.format( fname ))
plt.show()

