#from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


##============= variables in cgs =====================
h = 6.6260755e-34         #J*s
c = 2.99792458e8          #m/s
T_pot =   533.15          #Kelvin
kb = 1.380658e-23         #boltzmann
v = 6e-14              #

#============== Planck ========================
def Planck(temp, lamb):
    global h, c, kb
#    print h, c, kb
    lamb = lamb*10**(-6)         #from micro to meters
#    print 'lamb:',lamb
    exterm = ((h*c)/(lamb*kb*temp))
#   print 'exterm:', exterm
    numerator = (2*h*(c**2))
#    print 'numerator:', numerator
    denominator = (lamb**5)*(np.exp(exterm)-1.0)
#    print 'denominator:', denominator
    I = numerator/denominator

    return I    

def Plancks(temp, v):
    global h, c, kb
             #Angstroms to meters
    exterm = (h*v)/(kb*temp)
    numerator = 2*h*(v**3)
    denominator = (c**2)*(np.exp(exterm)-1.0)
    R = numerator/denominator

    return R

I = Planck(310, (9.35))
print 'I = ',I

R = Plancks(500, v)
print 'R = ', R



wave = np.arange(1, 40, .01)           #creates an array (initial index, final index, increment jump)

Intense = Planck(T_pot, wave)

free = np.logspace(12, 14.5, 1000)

I_V = Plancks(533.13, free)    
    
plt.plot(wave, Intense)
plt.title("Hot Potato 533.15K")
plt.ylabel("Intensity B_v")
plt.xlabel("Wavelength (nm)")
plt.savefig("Hot_Potato_BB2.pdf")
plt.show()
plt.plot(free, I_V)
plt.title("Hot Potato 533.15K")
plt.ylabel("Intensity B_v")
plt.xlabel("Frequancy (Hz)")
plt.savefig("Hot_Potato_BB.pdf")
plt.show()
