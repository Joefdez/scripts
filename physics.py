from toolBox import *
from constants import *
from scipy.special import jv

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


def orbFreqEcc(ee, ff0, ee0):
    # Compute orbital frequency as a function of eccentricity, given the initial orbital binary eccentricity and frequency arXvi:1805.06194 eq 5.

    return ff0*((1.-ee0**2.)/(1.-ee**2.) * (ee/ee0)**(12./19.
            ) * (1.+121./304. * ee**2.)/(1.+121./304. * ee0**2.
            )**(870./229))**(-3./2.)

def peakFreqEcc(ee, ff0, ee0):
    # Compute the peak gravitational wave frequency given the initial eccenricity and initial orbital frequency
    orbF = orbFreqEcc(ee, ff0, ee0)

    return 2.*orbF*(1-ee)**(-3./2.)

def chirpMass(m1, m2):

    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

def eccF(ee):

    return (1. + (73./24.) * ee**2. + (37./96.) * ee**4.)/((1.-ee**2.)**(7./2.))

def harmG(nn, ee):

    ne = nn*ee
    cc1 = jv(nn-2, ne) - 2.*ee*jv(nn-1,ne) + 2./nn * jv(nn, ne) + 2.*ee*jv(nn+1, ne) - jv(nn+2, ne)
    cc2 = (1.-ee**2.) * (jv(nn-2,ne) - 2.*jv(nn, ne) + jv(nn+2, ne))**2.
    cc3 = 4./(3.*nn**2.) * jv(nn,ne)**2.

    return (nn**4.)/32. * (cc1**2. + cc2 + cc3)


def dEdf(cMass, ff, zz, ee, nn):

    gg = harmG(nn, ee)
    enhancement = eccF(ee)
    preFac =  (GG*cMass)**(5./3.)/(3.*pi**(1./3.)*(1.+ zz)**(1./3.)*ff**(1./3.))

    return preFac*(2./nn)**(2./3.) * gg/enhancement


def eccStrainNN(cMass, ff, ee, zz, nn, DD):

    return 1./(pi*DD) * sqrt(2.*GG/cc**3. * dEdf(cMass, ff, zz, ee, nn))


def circStrain(cMass, ff, DD):

    return 1./DD * 2.*(4.*pi)**(1./3.)*GG**(5./3.)/cc**4. * ff**(2./3.) * cMass**(5./3.)
