import numpy as np
import pylab as P
import urllib2
import matplotlib.pyplot as plt

#url='http://www.astro.yale.edu/smarts/glast/tables/3C279.tab'
#open("3c279.dat",'wb').write(urllib2.urlopen(url).read())



jdstart,bjdstart,bmag,berr,vjdstart,vmag,verr,rjdstart,rmag,rerr,jjdstart,jmag,jerr,kjdstart,kmag,kerr = np.genfromtxt('3c279.dat', unpack=True,dtype=float,skip_header=3)


print jdstart.dtype
print bjdstart.dtype
print vjdstart.dtype
print rjdstart.dtype
print jjdstart.dtype
print kjdstart.dtype
print bmag.dtype
print vmag.dtype
print rmag.dtype
print jmag.dtype
print kmag.dtype
print berr.dtype
print verr.dtype
print rerr.dtype
print jerr.dtype
print kerr.dtype


#print bjdstart
#bjdstart = "999."

#if bjdstart == "999.":
    #print('ummmmm')

#print bjdstart    
# all objects are arrays (numpy.ndarray)


print type(jdstart)
print type(bjdstart)
print type(vjdstart)
print type(rjdstart)
print type(jjdstart)
print type(kjdstart)
print type(bmag)
print type(vmag)
print type(rmag)
print type(jmag)
print type(kmag)
print type(berr)
print type(verr)
print type(rerr)
print type(jerr)
print type(kerr)

#plotting error bars
#plt.errorbar(jdstart,bmag-jmag, bmagerr=berr, jmagerr = jerr)
#bad = '999.0'
#bad2 = '-999.0'

#if bad or bad2 in jdstart or bjdstart or bmag or berr or vjdstart or vmag or verr or rjdstart or rmag or rerr or jjdstart or jmag or jerr or kjdstart or kmag or kerr:
    #print "bad"
#else:
    #print bmag

#np.float
#print len()

#jdstart = jdstart[595:]
#bmag = bmag[595:]
#jmag = jmag[595:]




#plt.scatter(jdstart,bmag-jmag)
#plt.title("b-j vs. time")
#plt.xlabel("J")
#plt.ylabel("B-J")
#plt.ylim([2.5,5])
#plt.xlim()
#plt.show()



