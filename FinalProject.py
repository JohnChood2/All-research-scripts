import numpy as np
import matplotlib.pyplot as plt

#=============== IMPORTANT VALUES FOR THE MODEL =============================================

X = 0.7                             #SOLAR COMPOSITION
Y = 0.28                            #SOLAR COMPOSITION
Z = 0.02                            #METALLICITY
pi = np.pi                          #PI
Rsun = 6.955e10                     #SOLAR RADIUS IN CGS
G = 6.67259e-8                      #GRAVITATIONAL CONSTANT
k = 1.380658e-16                    #BOLTZMAN CONSTANT
mp = 1.672631e-24                   #PROTON MASS
Pc = 10.0**(16.675)                 #CENTRAL BOUNDRY PRESSURE IN LOG SCALE
Tc = 10.0**(7.46)                   #CENTRAL BOUNDRY TEMPERATURE LOG SCALE IN KELVIN
M_sun = 1.99e33                     #MASS OS THE SUN IN CGS UNITS
M1 = 7.08*M_sun                     #SOLAR MASS FOR PART 1
M2 = 2.82*M_sun                     #SOLAR MASS FOR PART 2
a = 7.5e-15                         #RADIATION CONSTANT
N = 6.022e23                        #AVACADOS NUMBER
me = 9.1094e-28                     #ELECTRON MASSS
mu_e = 1.176                        #MOLECULAR WEIGHT OF GAS NOT IONIZED
mu_i = 1.298                        #MOLECULAR WEIGHT OF GAS WHEN IONIZED
mu = .613                           #MOLECULAR WEIGHT OF GAS TOTAL
c = 2.99792458e10                   #SPEED OF LIGHT IN CGS
L_sun = 3.9e33                      #SOLAR LUMINOSITY IN CGS
h = 6.6260755e-27                   #PLANCK CONSTANT IN CGS
dm1 = M1*1.0e-4                     #DM
Ae = 5.865e29                       #USED FOR MEMAN MOLECULAR WEIGHT
d = (k**4)/c                        #USED IN DT/DM RAD IN PLACE OF A
m = 0.0                             #SOME VALUE AT THE CENTER
m1_surf = M1/2                      #SOME SURFACE VALUE
m2_surf = M2/2                      #SOME SURFACE VALUE
r_center = 0                        #radius at center
#======================= solving for central variables =======================================

def find_rho(P, T):
    global mu, a, N, k
##    print mu
##    print a
##    print N
##    print k
##    print P
##    print T
    rho = (mu/(N*k*T))*(P-(1./3.)*a*(T**4))
    
    return rho
#rho = find_rho(Pc, Tc)
#r = ((3*dm/4*pi*rho)**(1./3.))

def dr_dm(r,rho):
    global pi
   # print "drdm", (1./(4*pi*r**2*rho))
    return (1./(4*pi*r**2*rho))

#dm = dr_dm
#print "radius =", r
#print "rho =", rho
#======================= Calcuating Central Pressure Variables ===================================================

def findP(rho, T):
    global mu_i, mu_e, pi, a, N, mu, k, h, c, me
    
    Prad = (a*(T**4))/3

    Pend = (N*k/mu_e)*rho*T

    Ped = (3./(8.0*pi))**(2./3)*((h**2)/(5.0*me))*(N/mu_e)**(5./3)*(rho)**(5./3)

    Pion = (N*k/mu_i)*rho*T

    Pe = Pend

    Pgas = Pe + Pion

    Ptot = Pion + Pe + Prad

    return Ptot, Pgas

#Ptot, Pgas = findP(rho, Tc)
#print "total pressure =", Ptot
#print "total gas pressure =", Pgas

#====================================== calculating Epsilon =====================================

def findE(rho, T):                             
    global X, Y, Z, Tc

    T9 = T/1e9
    
    Epp = (2.4e4*((rho*(X**2))/((T9)**(2./3.))))*np.exp((-3.38)/((T9)**(1./3.)))

    Ecno = (4.4e25*((rho*X*Z)/((T9)**(2./3.))))*np.exp((-15.228)/((T9)**(1./3.)))

    E_3a = (5.0e8*((rho**2*(Y**3))/((T9)**(3))))*np.exp((-4.4)/(T9))
 
    Et = Epp + Ecno + E_3a

    return Et
             
#Et = findE(11.41, (10**(7.45)))
#print "total energy =", Et

#================================== Calculating Luminosity =====================================
def Luminosity_center(Et, m):
    Lum_c = (Et*m)
    return (Et*m)
#Lum_c = (Et*m)
#Lc = Lum_c
#print "central luminosity =", Lc

def dL_dm(Et, dm):
    L = Et*dm
    return L
#L = Et*dm
#print "Luminocity =", L

#================================= Calculating Gamma ===========================================
def beta(Pgas, Ptot):
    beta2 = (Pgas/Ptot)
    return (Pgas/Ptot)
#beta2 = (Pgas/Ptot)
#B = beta2
#print "Beta initial =", B

def findGamma(B):
    Gamma = ((32 - 24*(B) - 3*(B**2))/(24 - 18*(B) - 3*(B**2)))
    return ((32 - 24*(B) - 3*(B**2))/(24 - 18*(B) - 3*(B**2)))

#print "Gamma_2 =", Ga

#====================================== Calculating Temperature =================================

def dt_dm_rad(k,L,pi,a,c,r,Tc):
    Trad = ((-3.0*k*L)/(64.0*(pi**2)*a*c*(r**4)*(Tc**3)))
    return Trad
#Trad = (3*k*L/64*(pi**2)*a*c*(Tc**3))

def dt_dm_conv(Ga, G, m, Tc, pi, r, Ptot):

    Tcon =  -(1 - 1/Ga)*((G*m*Tc)/(4*pi*(r**4)*Ptot))
   # print "dtdm", Tcon
    return Tcon
#Tcon =  -(1 - 1/Ga)*(G*m*Tc/4*pi*(r**4)*Ptot)

#print "radiation temp =", Trad
#print "convection temp =", Tcon
#=================================== Calculating Pressure ======================================

def dP_dm(m, r):
    global G, pi
    
    Press = -(G*m/(4*pi*(r**4)))
 #   print "dpdm", Press
    return Press

#Press = dP_dm(G, M1)

#print "Pressure =", Press

def Pressure(T, rho):
    global mu, k, N, a
    Pre = ((1./3.)*a*T**4) + ((N*k/mu)*rho*T)
    return Pre

#P = Pressure
#======================================= Calculating Opacity ===================================

name1 = "/Users/johnhood/Desktop/StellarData/OPAL_X0.7_Y0.28.dat"
data = np.genfromtxt(name1, unpack = True)

def f_x(x):
    return (x**7)/(np.e**(x)*(1.0-np.e**(-x))**3)

def linear_interpol(x1,x2,y1,y2,X):
    m = ((y2-y1)/(x2-x1))
    return m*X-m*x1+y1

def get_opacity(rho, T, data):

    R = np.log10((rho*1e18)/(T**3))
    temp = np.log10(T)
    if R < data[1][0] or R > data[-1][0] or temp < data[0][1] or temp > data[0][-1]:
        kf = approx_opacity(temp,R)
        return kf
    r1 = 0
    r2 = 0

    data2 = np.where((data[0] > temp)&(data[0]!=9999))
    
    t1 = data[0][data2[0][0]-1]
    t2 = data[0][data2][0]
    
    k1 = 0
    k2 = 0
    k3 = 0
    k4 = 0

    for j in range(len(data)-1):
        j += 1
        if data[j-1][0] < R:
            r1 = data[j-1][0]
            r2 = data[j][0]
            k1 = data[j-1][np.where(data[0] == t1)][0]
            k2 = data[j][np.where(data[0] == t1)][0]
            k3 = data[j-1][np.where(data[0] ==t2)][0]
            k4 = data[j][np.where(data[0] == t2)][0]

    if -9999 in [k1, k2, k3, k4]:
        kf = approx_opacity(temp, R)
        return kf

    
    kr = linear_interpol(r1,r2,k1,k2,R)
    kr2 = linear_interpol(r1,r2,k1,k2,R)
    kf = linear_interpol(t1,t2,kr,kr2,temp)
    return kf

def approx_opacity(T,r):
    global X,Y,Z
    R = 10**(r)
    T = 10**(T)
    rho = R*(T**3)/1e18
    n,s = 0.5, 9.
    k_h = (2.5e-31)*(Z/0.02)*(rho)*(T**(s))
    k_e = 0.2*(1.+X)
    n,s = 1, -3.5
    k_ff = (4.e22)*(X+Y)*(1+X)*(rho)*(T**(s))
    n,s = 1, -3.5
    k_bf = (4.e25)*Z*(1+X)*(rho)*(T**(s))
    return np.log10(1./((1./k_h) + 1./(k_e + k_ff + k_bf)))

#Go = 10**(get_opacity(rho, Tc, data))
#print "Opcaity =", Go

#=============================radiation and convection check =====================================

#T_r = Tc-((Go*rho**2*Et/8*a*c*Tc**3)*r**2)                      #radiative

#T_r2 = Tc-(1-1./Ga)*((2*pi*G*(rho**2)*Tc/Pc)*r**2)              #convective

def check_conv(L, Go, Ga, T, m, P, pi):
    global a, c, G
    A = (16.0*pi*a*c*G)/(3.0*Go)
    B = (1.0 - (1.0/Ga))
    C = ((T**4)*m)/P

    if L > (A*B*C):
        return True, L/(A*B*C)
    else: return False, L/(A*B*C)
#============================= CALCULATING CENTRAL VALUES AND FIRST STEP =========================
mass = np.linspace(0, M1, 10**4)
dm = mass[1]
#print dm
#### CENTRAL DENSITY
rho = find_rho(Pc, Tc)
#print rho
#### BETA
Ptot, Pgas = findP(rho, Tc)
beta2 = (Pgas/Ptot)
B = beta2
#print B
#### CENTRAL GAMMA
Gamma = ((32 - 24*(B) - 3*(B**2))/(24 - 18*(B) - 3*(B**2)))
Ga = Gamma
#print Ga
#### CENTRAL OPACITY
Go = 10**(get_opacity(rho, Tc, data))
#print Go
#### RADIUS AT FIRST STEP
r = ((3*dm/(4*pi*rho))**(1./3.))
#print r
#### ENERGY
Et = findE(11.41, (10**(7.45)))
#print Et
#### LIMINOCITY AT FIRST STEP
L = Et*dm

#### PRESSURE AT FIRST STEP
P = Pc - (2./3.)*pi*G*(rho**2)*r**2

#### TEMPERATURE AT FIRST STEP
O, S = check_conv(L,Go,Ga,Tc,m,Pc, pi)
if O:
    T = Tc-(1-1./Ga)*((2*pi*G*(rho**2)*Tc/Pc)*r**2)
else:
    T = Tc-((Go*rho**2*Et/8*a*c*Tc**3)*r**2)

#### ENERGY
#Et = findE(11.41, (10**(7.45)))



#==================================== LOOP OF DOOM ===============================================

#================ Lists needed ============================

R = [0.0,r]                                  #radius
Rho = [rho]                                  #density
Pressu = [Pc,P]                              #Pressure
E = [Et]                                     #epsilon
Gamma = [Ga]                                 #Gamma
Kappa = [Go]                                 #Opacity
Temp = [Tc,T]                                #Temperature
Lum = [0.0,L]                                #Luminocity
Rrag = [0]
##print R
##print Rho
##print Pressu
##print E
##print Gamma
##print Kappa
##print Temp
##print Lum

##def dr_dm(rho):
##    global pi, r
##    return (1./(4*pi*r**2*rho))
##dr = dr_dm(R[-1],rho[-1]*dm + R[-1])
##
#mass = np.linspace(0, M1, 10**4)
#dm = mass[1]
for i in range(len(mass)-2):
    i += 2
    rho = find_rho(Pressu[-1], Temp[-1])
    
    if rho < 0:
        print "rho broke"
        print rho
        break
    ptot, pgas = findP(rho,Temp[-1])
    
    B = beta(Pgas, Ptot)
    
    Gammas = findGamma(B)
    
    Kappas = 10**(get_opacity(rho, Temp[-1], data))
    
    radius = R[-1] + dr_dm(R[-1],rho)*dm
    
    Energy = findE(rho, Temp[-1])
    
    Lumin = dL_dm(Energy, dm) + Lum[-1]
    
    Pres = dP_dm(mass[i], R[-1])*dm + Pressu[-1]
    
    if Pres < 0:
        print "presure broke"
        break
    con, spock = check_conv(Lum[-1], Kappa[-1], Gamma[-1], Temp[-1], mass[i], Pressu[-1], pi)
    if con:
       #print "con"
       T = Temp[-1]+dt_dm_conv(Gamma[-1], G, mass[i], Temp[-1], pi, R[-1],Pressu[-1])*dm
    else:
       #print "rad"
       T = Temp[-1]+dt_dm_rad(Kappa[-1], Lum[-1], pi, a, c, R[-1], Temp[-1])*dm
    if T < 0:
        print " Temp broke"
        break
    R.append(radius)
    Rho.append(rho)
    Pressu.append(Pres)
    E.append(Energy)
    Gamma.append(Gammas)
    Kappa.append(Kappas)
    Temp.append(T)
    Lum.append(Lumin)
    Rrag.append(spock)

print 'solar Luminosity =', Lum[-1]/L_sun
print 'solar Radius =', R[-1]/Rsun
print 'solar Mass =', mass[-1]/M_sun


print "done"    

#=================================== PLOTTING GRAPHS ===================================================

plt.subplot(4, 2, 1)
plt.plot(np.array(mass)/M_sun, np.log10(Pressu))
plt.title('Pressure VS. Mass')
plt.ylabel('presure')
plt.xlabel('mass')

plt.subplot(4, 2, 2)
plt.plot(np.array(mass)/M_sun, np.log10(Temp))
plt.title('Temperature VS. Mass')
plt.ylabel('Temperature')
plt.xlabel('mass')

plt.subplot(4, 2, 3)
plt.plot(np.array(mass[:-1])/M_sun, np.log10(Rho))
plt.title('Density VS. Mass')
plt.ylabel('density')
plt.xlabel('mass')

plt.subplot(4, 2, 4)
plt.plot(np.array(mass)/M_sun, np.array(R)/Rsun)
plt.title('radius VS. mass')
plt.ylabel('radius')
plt.xlabel('mass')

plt.subplot(4, 2, 5)
plt.plot(np.array(mass)/M_sun, np.array(Lum)/L_sun)
plt.title('Luminosity VS. Mass')
plt.ylabel('luminosity')
plt.xlabel('mass')

plt.subplot(4, 2, 6)
plt.plot(np.array(mass[:-1])/M_sun, Kappa)
plt.title('Opacity VS. Mass')
plt.ylabel('opacity')
plt.xlabel('mass')
plt.ylim(0,10)

plt.subplot(4, 2, 7)
plt.plot(np.array(mass[:-1])/M_sun, E)
plt.title('Energy VS. Mass')
plt.ylabel('energy')
plt.xlabel('mass')

plt.subplot(4, 2, 8)
plt.plot(np.array(mass[:-1])/M_sun, Rrag)
plt.title('ratio VS. mass')
plt.ylabel('ratio')
plt.xlabel('mass')
plt.ylim(0,10)


plt.tight_layout()
plt.show()

