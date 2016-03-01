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

y_axis1 = 'FLUX_100_300000'
y_axis_err1 = 'ERROR_100_300000'
UL1 = 'UL_100_300000'

y_axis2 = 'FLUX_1000_300000'
y_axis_err2 = 'ERROR_1000_300000'
UL2 = 'UL_1000_300000'

y_axis3 = 'FLUX_300_1000'
y_axis_err3 = 'ERROR_300_1000'
UL3 = 'UL_300_1000'

Plot_color = 'k'

x_axe = tbdata[x_axis]

##============== ALL FLUX DATA ===========
y_axe1 = tbdata[y_axis1]
y_axe_err1 = tbdata[y_axis_err1]
UL_axes1   = tbdata[UL1]

y_axe2 = tbdata[y_axis2]
y_axe_err2 = tbdata[y_axis_err2]
UL_axes2   = tbdata[UL2]

y_axe3 = tbdata[y_axis3]
y_axe_err3 = tbdata[y_axis_err3]
UL_axes3   = tbdata[UL3]

#mod_x = MET_to_MJD(x_axe)

def MET_to_MJD(x_axe):
    #convert MET to days since 2001
    MJD = ((((x_axe)/(3.15e7))*(365))+ (51910.00074287037))
   
    return MJD

MJD = ((((x_axe)/(3.15e7))*(365))+ (51910.00074287037))

x_axe = MJD

#print type(x_axe)
#print type(MJD)
#print MJD
#print UL_axes1
#print UL_axes2
#print UL_axes3
##============ data filter =================
if y_axis1 == 'FLUX_100_300000':
    y_id1 = np.where(y_axe1 >= 0)[0]   
    x_axe_c = x_axe            [y_id1]   
    y_axe1_c = y_axe1          [y_id1]
    y_axe_err1_c = y_axe_err1  [y_id1]
    
elif y_axis2 == 'FLUX_1000_300000':   
    x_axe_c = x_axe            [y_id2]
    y_id2 = np.where(y_axe2 >= 0)[0]
    y_axe2_c = y_axe2          [y_id2]
    y_axe_err2_c = y_axe_err2  [y_id2]
    
elif y_axis3 == 'FLUX_300_1000':   
    x_axe_c = x_axe            [y_id3]
    y_id3 = np.where(y_axe3 >= 0)[0]
    y_axe3_c = y_axe3          [y_id3]
    y_axe_err3_c = y_axe_err3  [y_id3]
    
IDX1 = np.where(UL_axes1 == 'F')[0]
IDX2 = np.where(UL_axes2 == 'F')[0]
IDX3 = np.where(UL_axes3 == 'F')[0]

t_N1 = x_axe       [IDX1]
F_N1 = y_axe1      [IDX1]
E_N1 = y_axe_err1  [IDX1]

t_N2 = x_axe       [IDX2]
F_N2 = y_axe2      [IDX2]
E_N2 = y_axe_err2  [IDX2]

t_N3 = x_axe       [IDX3]
F_N3 = y_axe3      [IDX3]
E_N3 = y_axe_err3  [IDX3]

arrow_dat1 = np.where(UL_axes1 == 'T')[0]
arrow_dat2 = np.where(UL_axes2 == 'T')[0]
arrow_dat3 = np.where(UL_axes3 == 'T')[0]

# Directions

u1 = np.zeros(arrow_dat1.size)
v1 = np.ones(arrow_dat1.size)*-500.

u2 = np.zeros(arrow_dat2.size)
v2 = np.ones(arrow_dat2.size)*-500.

u3 = np.zeros(arrow_dat3.size)
v3 = np.ones(arrow_dat3.size)*-500.
print len(t_N2)
print len(F_N2)
print len(E_N2)
print len(IDX2)

##============ plotting =====================

##plt.plot(x_axis, y_axis)
##plt.xlabel('MJD (d)')
##plt.ylabel('Flux_300_1000')
##plt.title('Source 3C 454.3')
##plt.show()

fig = plt.figure(figsize = (12,12))
ax1 = fig.add_subplot( 221, axisbg='white')
ax1.set_xlabel('MJD (d)', fontsize=18)
ax1.set_ylabel('FLUX_100_300000', fontsize=18)
ax1.set_ylim(0, .000004)
ax1.set_title( 'Source = PMNJ2345-1555 DURA = 86400')
ax1.errorbar(t_N1, F_N1, yerr=E_N1, 
    color=Plot_color, markersize=2, fmt = 'o', linewidth=1, 
    ecolor=Plot_color, elinewidth=1, capthick = 1 
    )#, label=r'\boldmath $\overline{\mathrm{Cens}}$' )
ax1.quiver(x_axe[arrow_dat1], y_axe1[arrow_dat1],u1,v1, color='r', zorder=5)

ax2 = fig.add_subplot( 222, axisbg='white')
ax2.set_xlabel('MJD (d)', fontsize=18)
ax2.set_ylabel('FLUX_1000_300000', fontsize=18)
ax2.set_ylim(0, .0000005)
ax2.set_title( 'Source = PMNJ2345-1555 DURA = 86400')
ax2.errorbar(t_N2, F_N2, yerr=E_N2, 
    color=Plot_color, markersize=2, fmt = 'o', linewidth=1, 
    ecolor=Plot_color, elinewidth=1, capthick = 1 
    )#, label=r'\boldmath $\overline{\mathrm{Cens}}$' )
ax2.quiver(x_axe[arrow_dat2], y_axe2[arrow_dat2],u2,v2, color='r', zorder=5)

ax3 = fig.add_subplot( 223, axisbg='white')
ax3.set_xlabel('MJD (d)', fontsize=18)
ax3.set_ylabel('FLUX_300_1000', fontsize=18)
ax3.set_ylim(0, .000004)
ax3.set_title( 'Source = PMNJ2345-1555 DURA = 86400')
ax3.errorbar(t_N1, F_N1, yerr=E_N1, 
    color=Plot_color, markersize=2, fmt = 'o', linewidth=1, 
    ecolor=Plot_color, elinewidth=1, capthick = 1 
    )#, label=r'\boldmath $\overline{\mathrm{Cens}}$' )
ax3.quiver(x_axe[arrow_dat3], y_axe3[arrow_dat3],u3,v3, color='r', zorder=5)

#ax1.set_yscale('log')
#ax2.set
#ax3.set_yscale('log')
#ax.set_xscale('log')
plt.tight_layout()
fname = '#PMNJ2345-1555.png'#.format(x_axis, y_axis)
plt.savefig( fname, bbox_inches='tight' )
print( '    Figure: {0}'.format( fname ))
plt.show()
