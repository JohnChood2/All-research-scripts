import numpy as np
import pylab as P
import urllib
import matplotlib.pyplot as plt

jdstart,bjdstart,bmag,berr,vjdstart,vmag,verr,rjdstart,rmag,rerr,jjdstart,jmag,jerr,kjdstart,kmag,kerr = np.genfromtxt('http://www.astro.yale.edu/smarts/glast/tables/3C279.tab', unpack=True,dtype=float,skip_header=3)
#np.float
#print len()

#jdstart = jdstart[595:]
#bmag = bmag[595:]
#jmag = jmag[595:]
yarray = (jdstart,bjdstart,bmag,berr,vjdstart,vmag,verr,rjdstart,rmag,rerr,jjdstart,jmag,jerr,kjdstart,kmag,kerr)
yarray2 = (bmag,jmag,rmag,vmag,kmag,)

gdats = []

for i in yarray:
    gdats.append(np.where((i!=999.)&(i!=-999.)))

##gdat = np.where((yarray!=999.0) & (yarray!=-999.0))
##gdat2 = np.where((yarray2!=999.0) & (yarray2!=-999.0))
##
##
##plt.scatter(jdstart,bmag-jmag)
##plt.title("b-j vs. time")
##plt.xlabel("J")
##plt.ylabel("B-J")
##plt.ylim([2000,-2000])
##plt.xlim()
##plt.show()

plt.scatter(jdstart[gdats[2]], bmag[gdats[2]]-jmag[gdats[2]])
plt.show()
