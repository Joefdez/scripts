#!/usr/bin/python
from mpi4py import MPI
from numpy import *
from time import sleep
from numpy.random import seed, uniform
from scipy.integrate import quad

comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank   = comm.Get_rank()

# Constants
Msun = 1.989e30      # Solar mass, kg
GG   = 6.674e-11     # Newton's gravitational constant, SI units
cc   = 2.997e8       # Speed of light in vacuum, m/s
pc   = 3.086e16      # 1 parsec in metres
au   = 1.496e11
# Problem data
MM = 4.0e6*Msun # In solar masses
m1 = 30.0*Msun  # In solar masses
m2 = 30.0*Msun  # In solar masses
# Scharzschild radii
RS = 2.*GG*MM/(cc**2.)
rs = 2.*GG*m1/(cc**2.)

# Other parameters
beta = 64./5.*(m1*m2*(m1+m2)*  GG**3)/(cc**5)
Mf1  = (MM/(m1+m2))**(1./3.)
Rs   = 2.*GG*MM/(cc**2)       # EH radius
esma = 1.*pc                  # Semi-major axis of elliptic orbit
kepT = 2.*pi*1./sqrt(GG*MM) * esma**(3./2.)


def mulQuad(ff, xx0, xx):
    #   Wrapper to use quad method for integrating area under curve on several arguments
    res = zeros(shape(xx))
    for ii in range(len(xx)):
       res[ii] = quad(ff, xx0[ii], xx[ii])[0]
       #print ii
    return res


def GWtimearg(ee):
    return ee**(29./19.) *(1.+121./304. *ee**2. )**(1181./2299.)/((1.-ee**2.)**(3./2.))


def GWtime(ee, cc0, func):
    # Calculate GW merger time for an array of initial eccentricities
    f1 = 12./19. * cc0**4/beta               # Numerical factor
    ze = zeros(len(ee))
    times = f1 * mulQuad(func, ze, ee)
    return times

print rank, "Reading data"
predata = loadtxt("survD075.dat")               # All cores loading is faster than sharing$
data    = predata[:,[0, 8, 9, 10, 14, 7]]       # Pull out D, ee0, ee, aa, chieff, Lz0
data[:,3] = data[:,3]*data[:,0]                 # Change aa for aa in units of initial aa0
data[data[:,5]>0][:,5] = 1        # Flag progrades
data[data[:,5]<0][:,5] = 0        # Flag retrogrades
#for jj in range(10):
nn = shape(data)[0]
res = zeros([nn,9])                         # Array for results   save D, e, sma pre encount$
res[:,0] = data[:,0]                        # Save D
res[:,1] = data[:,1]                        # save pre-encounter eccentricity
res[:,2] = data[:,2]                        # save post-encounter eccentricity
sleep(rank+5)                               # Put the process to sleep for a few seconds (so$
seed()                                      # Re-seed PRNG
aa0 = uniform(log10(0.01*au), log10(1.*au), nn)  # Uniform distribution in log space of ini$
#aa0 = uniform(0.01,1.,nn)        # log(aa0) uniformly distributed between 0.01 and 1. au
aa0 = (10**aa0)
res[:,3] = aa0                              # Rescaled initial separation
res[:,4] = aa0*data[:,3]                    # Rescaled final separation
cc0B  = aa0*(1. -res[:,1]**2.)*res[:,1]**(-12./19                    # cc0 before the encounter
                    ) * (1.+121./304. * res[:,1]**2.)**(-870./2299.)
cc0A  = res[:,4]*(1. -res[:,2]**2.)*res[:,2]**(-12./19               # cc0 after the encounter
                    ) * (1.+121./304. * res[:,2]**2.)**(-870./2299.)

res[:,5] = GWtime(res[:,1], cc0B, GWtimearg)                            # Save pre-encounter merger time
res[:,6] = GWtime(res[:,2], cc0A, GWtimearg)                            # Save post-encounter merger time
res[:,7] = data[:,-2]                                                     # Save effective spins
res[:,8] = data[:,-1]
print rank, "finished"
filname = "pyoutput/mergTimes_" + str(rank) + "_" + ".dat"
savetxt(filname, res)                                                   # Save file
