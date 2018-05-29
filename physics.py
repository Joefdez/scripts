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

    return aa0**4./beta
