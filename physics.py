from toolBox import *
from constants import *

def cc0Factor(aa0, ee0):
    # Integration constant for Peter's equations
    c1  = aa0*(1.-ee0**2.)*ee0**(-12./19.)
    c2 = 1. + (121./304.)*(ee0**2.)

    return c1*(c2**(-870./2299.))

def betaFactor(m1, m2):
    # Beta factor of binary of componentmasses m1, m2

    return (64./5.)*GG**3. * m1*m2*(m1+m2)/(cc**5)

    return

def mergerTimeInt(ee):
    # Integral of merger time formula of a compact binary object due to GW radiation according to Peters equations:
    # "Gravitational Radiation and the Motion of Two Point Masses", Phys. Rev. 136 4b, 1964
    # ee: eccentricity

    num = ee**(29./19.) * (1. + (121./304.)*ee**2.)**(1181./2299.)
    den = (1. - ee**2)**(3./2.)

    return num/den

def mergerTimeecc(aa0s, ee0s, beta):
    # Merger time formula of an eccentric compact binary object due to GW radiation according to Peters equations:
    # "Gravitational Radiation and the Motion of Two Point Masses", Phys. Rev. 136 4b, 1964
    # Assumes equal mass binaries (i.e single value of beta factor)
    # ee0s: array of initial eccentricities
    # aa0s: array of initial semi-major axes

    numS = shape(aa0s)
    cc0s = cc0Factor(aa0s, ee0s)

    prefactor = (12./19.) * cc0s**4. / beta

    ints = quadArray(mergerTimeInt, zeros(numS), ee0s)

    return prefactor*ints

def mergerTimecirc(aa0s, beta):
    # Merger time of a circular binary ue to GW radiation according to Peters equations:
    # "Gravitational Radiation and the Motion of Two Point Masses", Phys. Rev. 136 4b, 1964
    # Assumes equal mass binaries (i.e single value of beta factor)
    # aa0s: array of initial semi-major axes

    return aa0s**4./beta

def minSMA(ees, tt, beta):
    # semi-major axis at which given an eccentricity ee the binary with beta factor beta will merge in a time tt
    # Inputs are asrrays, at least ee.

    numS = shape(ees)
    ints = quadArray(mergerTimeInt, zeros(numS), ees)

    prefactor = ees**(12./19.)/(1.-ees**2.) * (1. + (121./304.)*(ees**2.))**(870./2299.)

    return prefactor*(19./12. * tt*beta/ints)**(1./4.)


def deltaE_encounter(MM, mm1, mm2, rm):
    # Calculate the change in energy of a mass upon a close encounter with another mass
    # on a parabolic orbit (e=1)

    # rm -> periastron distance
    # MM -> MBH mass, mm1, mm2 -> binary member masses

    ge = (85./12.) * pi
    totM =  MM + mm1 + mm2
    mf = sqrt(totM) * MM**2. * (mm1+mm2)**2.

    delE = (GG**(7./2.)/cc**5.) * mf/(rm**(7./2.)) * ge


    return delE

def binary_energy(mm1, mm2, aa):
    # Binding energy of a binary
    # It is negativ (returned with positive sign though)
    return GG*(mm1*mm2)/(2*aa)
