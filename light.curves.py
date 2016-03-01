import numpy as np
import pylab as P
import urllib
import matplotlib.pyplot as plt

jdstart,bjdstart,bmag,berr,vjdstart,vmag,verr,rjdstart,rmag,rerr,jjdstart,jmag,jerr,kjdstart,kmag,kerr = np.genfromtxt('http://www.astro.yale.edu/smarts/glast/tables/3C279.tab', unpack=True,dtype=float,skip_header=3)
#np.float
#print len() 

plt.scatter(jdstart,kmag)
plt.title("Flux vs. time")
plt.xlabel("Time")
plt.ylabel("B")
plt.ylim([8,14])
plt.show()
