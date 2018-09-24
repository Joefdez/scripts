
from physics import *
from matplotlib.pylab import *
from matplotlib.colors import *
ion()

#rcParams['font.family'] = 'serif'
#rcParams['font.serif']  = 'dejavusans'
#rcParams['mathtext.fontset'] = 'dejavuserif'
#rcParams['mathtext.bf'] = 'dejavuserif'



nn = 3
data1 = loadtxt("/home/arijfern//Desktop/lsmc/ud_final.dat")
print shape(data1)[0]
data1 = data1[data1[:,-1]==1]
data1 = data1[data1[:,10]>0]
data1 = data1[data1[:,8]<0.951]
numS1 = shape(data1)[0]
initMT1 = zeros(nn*numS1)
finalMT1 = zeros(nn*numS1)
initA1   = zeros(nn*numS1)
finalA1  = zeros(nn*numS1)
#ratE1    = zeros(3*numS1)
rps1     = zeros(nn*numS1)
ts1      = zeros(nn*numS1)
ls1      = zeros(nn*numS1)
l12      = zeros(nn*numS1) # Angular momentum squared

data2 = loadtxt("/home/arijfern/Desktop/lsmc/pln1_final.dat")
print shape(data2)[0]
data2 = data2[data2[:,-1]==1]
data2 = data2[data2[:,10]>0]
data2 = data2[data2[:,8]<0.951]
numS2 = shape(data2)[0]
initMT2 = zeros(nn*numS2)
finalMT2 = zeros(nn*numS2)
initA2   = zeros(nn*numS2)
finalA2  = zeros(nn*numS2)
#ratE2    = zeros(3*numS2)
rps2     = zeros(nn*numS2)
ts2      = zeros(nn*numS2)
ls2      = zeros(nn*numS2)
l22      = zeros(nn*numS2) # Angular momentum squared


data3 = loadtxt("/home/arijfern/Desktop/lsmc/pln05_final.dat")
print shape(data3)[0]
data3 = data3[data3[:,-1]==1]
data3 = data3[data3[:,10]>0]
data3 = data3[data3[:,8]<0.951]
numS3 = shape(data3)[0]
initMT3 = zeros(nn*numS3)
finalMT3 = zeros(nn*numS3)
initA3   = zeros(nn*numS3)
finalA3  = zeros(nn*numS3)
#ratE3    = zeros(3*numS3)
rps3     = zeros(nn*numS3)
ts3      = zeros(nn*numS3)
ls3      = zeros(nn*numS3)
l32      = zeros(nn*numS3) # Angular momentum squared


m1 = 30.*Msun
m2 = 30.*Msun
beta = betaFactor(m1, m2)

MM = 4.e6*Msun

MMm = (MM/(m1+m2))**(1./3.)

minA = 1.5*au
maxA1 = (3.e2)*au
maxA2 = (3.e2)*au


print "Generating initial separations and initial and final merger times."
for ii in range(nn):

    initA1[ii*numS1:(ii+1)*numS1]     = uniformLog(minA, maxA1, numS1)
    initA2[ii*numS2:(ii+1)*numS2]     = uniformLog(minA, maxA2, numS2)
    initA3[ii*numS3:(ii+1)*numS3]     = uniformLog(minA, maxA1, numS3)
    finalA1[ii*numS1:(ii+1)*numS1]    = data1[:,0]*data1[:,10]*initA1[ii*numS1:(ii+1)*numS1]
    finalA2[ii*numS2:(ii+1)*numS2]    = data2[:,0]*data2[:,10]*initA2[ii*numS2:(ii+1)*numS2]
    finalA3[ii*numS3:(ii+1)*numS3]    = data3[:,0]*data3[:,10]*initA3[ii*numS3:(ii+1)*numS3]

    initMT1[ii*numS1:(ii+1)*numS1]  = mergerTimeecc(initA1[ii*numS1:(ii+1)*numS1], data1[:,8], beta)
    initMT2[ii*numS2:(ii+1)*numS2]  = mergerTimeecc(initA2[ii*numS2:(ii+1)*numS2], data2[:,8], beta)
    initMT3[ii*numS3:(ii+1)*numS3]  = mergerTimeecc(initA3[ii*numS3:(ii+1)*numS3], data3[:,8], beta)
    finalMT1[ii*numS1:(ii+1)*numS1] = mergerTimeecc(finalA1[ii*numS1:(ii+1)*numS1], data1[:,9], beta)
    finalMT2[ii*numS2:(ii+1)*numS2] = mergerTimeecc(finalA2[ii*numS2:(ii+1)*numS2], data2[:,9], beta)
    finalMT3[ii*numS3:(ii+1)*numS3] = mergerTimeecc(finalA3[ii*numS3:(ii+1)*numS3], data3[:,9], beta)

initMT1, initMT2, initMT3 = initMT1/sTy, initMT2/sTy, initMT3/sTy
finalMT1, finalMT2, finalMT3 = finalMT1/sTy, finalMT2/sTy, finalMT3/sTy

fil1 = finalMT1<10**10.
fil2 = finalMT2<10**10.
fil3 = finalMT3<10**10.

print "Calculating tidal radii"
rt1 = MMm*initA1
rt2 = MMm*initA2
rt3 = MMm*initA3

print "Calculating periastron distances, length scales, time scales."
for ii in range(nn):
    rps1[ii*numS1:(ii+1)*numS1] = data1[:,0]*rt1[ii*numS1:(ii+1)*numS1]
    rps2[ii*numS2:(ii+1)*numS2] = data2[:,0]*rt2[ii*numS2:(ii+1)*numS2]
    rps3[ii*numS3:(ii+1)*numS3] = data3[:,0]*rt3[ii*numS3:(ii+1)*numS3]
    ls1[ii*numS1:(ii+1)*numS1] =  rps1[ii*numS1:(ii+1)*numS1]/MMm
    ls2[ii*numS2:(ii+1)*numS2] =  rps2[ii*numS2:(ii+1)*numS2]/MMm
    ls3[ii*numS3:(ii+1)*numS3] =  rps3[ii*numS3:(ii+1)*numS3]/MMm
    ts1[ii*numS1:(ii+1)*numS1] =  sqrt((rps1[ii*numS1:(ii+1)*numS1])**3./(GG*MM))
    ts2[ii*numS2:(ii+1)*numS2] =  sqrt((rps2[ii*numS2:(ii+1)*numS2])**3./(GG*MM))
    ts3[ii*numS3:(ii+1)*numS3] =  sqrt((rps3[ii*numS3:(ii+1)*numS3])**3./(GG*MM))


#initMT  = initMT/sTy
#finalMT = mergerTimeecc(data[:,0]*data[:,10]*initA, data[:,9], beta)
#finalMT = finalMT/sTy


# data -> array object from data file

fe, ae     = subplots(1,2)
#fa, aa     = subplots(1,2)
ft, at     = subplots(1,2)
ftt1, att1   = subplots(1,2)
ftt2, att2   = subplots(1,2)
ftt3, att3   = subplots(1,2)
fdelGE, adelGE = subplots(1,2)
fdelE0, adelE0 = subplots(1,2)
fE0GE, aE0GE   = subplots(1,2)
fef, aef   = subplots(1,2)
#fee, aee   = subplots(1,2)
#fte, ate   = subplots(1,2)
#fre, are   = subplots(1,2)

#fe, ae1     = subplots()
#ae2 = ae1.twinx()
#fa, aa1     = subplots()
#aa2 = aa1.twinx()
#ft, at1     = subplots()
#at2 = at1.twinx()
#ftt1, att11   = subplots()
#att12 = att11.twinx()
#ftt2, att21   = subplots()
#att22 = att21.twinx()
#ftt3, att31   = subplots()
#att32 = att31.twinx()
#fef, aef1   = subplots()
#aef2 = aef1.twinx()
#fee, aee1   = subplots()
#aee2 = aee1.twinx()
#fte, ate1   = subplots()
# ate2 = ate1.twinx()
#fre, are1   = subplots()
#are2 = are1.twinx()



# Calculate eccentricity distribution
ecc1 = concatenate([data1[:,9], data1[:,9], data1[:,9]])
ecc1 = 1. - ecc1[fil1]
vals1, bins1 = histogram(ecc1, bins=100, density=True)
points1 = 0.5*(bins1[:-1]+bins1[1:])
cumDist1 = cumsum(vals1*diff(bins1))*100

ecc2 = concatenate([data2[:,9], data2[:,9], data2[:,9]])
ecc2 = 1. - ecc2[fil2]
vals2, bins2 = histogram(ecc2, bins=100, density=True)
points2 = 0.5*(bins2[:-1]+bins2[1:])
cumDist2 = cumsum(vals2*diff(bins2))*100

ecc3 = concatenate([data3[:,9], data3[:,9], data3[:,9]])
ecc3 = 1. - ecc3[fil3]
vals3, bins3 = histogram(ecc3, bins=100, density=True)
points3 = 0.5*(bins3[:-1]+bins3[1:])
cumDist3 = cumsum(vals3*diff(bins3))*100
# Plot eccentricity distribution
ae[0].step(points1, vals1, linewidth=2, label=r'$p(e_{0}) \propto const.$')
ae[0].step(points2, vals2, linewidth=2, label=r'$p(e_{0}) \propto e_{0}$')
ae[0].step(points3, vals3, linewidth=2, label=r'$p(e_{0}) \propto e_{0}^{-1/2}$')
ae[0].set_xscale("log")
ae[1].loglog(points1, cumDist1)
ae[1].loglog(points2, cumDist2)
ae[1].loglog(points3, cumDist3)
#ae[1].step(points2, vals2, linewidth=2)
#ae[1].step(points2, vals2, linewidth=2)
#ae[1].loglog(points, cumDist, linewidth=2)
ae[0].set_xlabel(r"$1-e$", fontsize=24)
ae[0].set_ylabel(r"$p(1-e)$", fontsize=24)
ae[1].set_xlabel(r"$1-e$", fontsize=24)
ae[1].set_ylabel(r"Cumulative fraction $[\%]$", family='sans-serif', fontsize=24)
"""
ae1.step(points, vals, linewidth=2)
ae2.loglog(points, cumDist, linewidth=2, color='red')
ae1.set_xlabel(r"$1-e$", fontsize=24)
ae1.set_ylabel(r"$p(1-e)$", fontsize=24)
ae2.set_xlabel(r"$1-e$", fontsize=24)
ae2.set_ylabel(r"Cumulative fraction $[\%]$", family='sans-serif', fontsize=24)
ae1.set_xlim([0.005,1])
ae1.grid()
ae2.grid()
"""

"""
# Calculate semi-major axis distribution
aas = data[:,10]*data[:,0]
binse = logspace(log10(aas.min()), log10(aas.max()), num=100, endpoint=True, base=10.0)
vals, bins = histogram(aas, bins=binse, density=True)
points = 0.5*(bins[:-1]+bins[1:])
cumDist = cumsum(vals*diff(bins))*100
# Plot semi-major axis distribution
#aa[0].step(points, vals, linewidth=2)
#aa[1].loglog(points, cumDist, linewidth=2)
#aa[0].set_xlabel(r"$a/a_{0}$", fontsize=24)
#aa[0].set_ylabel(r"$p(a/a_{0})$", fontsize=24)
#aa[1].set_xlabel(r"$a/a_{0}$", fontsize=24)
#aa[1].set_ylabel(r"Cumulative fraction $[\%]$", fontsize=24)
aa1.step(points, vals, linewidth=2)
aa2.loglog(points, cumDist, linewidth=2, color='red')
aa1.set_xlabel(r"$a/a_{0}$", fontsize=24)
aa1.set_ylabel(r"$p(a/a_{0})$", fontsize=24)
aa2.set_xlabel(r"$a/a_{0}$", fontsize=24)
aa2.set_ylabel(r"Cumulative fraction $[\%]$", fontsize=24)
aa1.set_xlim([10**(-1.5) ,10**(1.5)])
aa1.set_ylim([0.,1.53])
aa2.set_ylim([0.01,100.])
aa1.grid()
aa2.grid()
"""
# Calculate merger time ratio distribution
#mt = finalMT/initMT
#binse = logspace(log10(mt.min()), log10(mt.max()), num=100, endpoint=True, base=10.0)
#vals, bins = histogram(mt, bins=binse, density=True)
#points = 0.5*(bins[:-1]+ bins[1:])
#cumDist = cumsum(vals*diff(bins))*100
#vals, bins = histogram(log10(mt), bins=100, density=True)
#points2 = 0.5*(bins[:-1]+ bins[1:])
#cumDist = cumsum(vals*diff(bins))*100
# Plot merger time ratio distribution
#at[0].step(points2, vals, linewidth=2)
#at[1].loglog(points, cumDist, linewidth=2)
#at[0].set_xlabel(r"$\log_{10}(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)
#at[0].set_ylabel(r"$p(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)
#at[1].set_xlabel(r"$\tau_{GW, out}/\tau_{GW, in}$", fontsize=24)
#at[1].set_ylabel(r"Cumulative fraction $[\%]$", fontsize=24)
"""
at1.step(points2, vals, linewidth=2)
at2.semilogy(points2, cumDist, linewidth=2, color='red')
at1.set_xlabel(r"$\log_{10}(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)
at1.set_ylabel(r"$p(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)
at2.set_xlabel(r"$\tau_{GW, out}/\tau_{GW, in}$", fontsize=24)
at2.set_ylabel(r"Cumulative fraction $[\%]$", fontsize=24)
at1.set_xlim([-10.,10.])
at1.set_ylim([0.,0.64])
at2.set_ylim([0.,100.])
at1.grid()
at2.grid()
"""


# Initial and final merger time distributions
valsb1, binsb1 = histogram(log10(initMT1/(10.**10)), bins=100, density=True)
pointsb1 = 0.5*(binsb1[:-1]+ binsb1[1:])
cumdistB1 = cumsum(valsb1*diff(binsb1))*100
att1[1].plot(pointsb1, cumdistB1, color='blue', linewidth=2)
valsb1, pointsb1 = append(zeros(1), valsb1), append(pointsb1[0:1], pointsb1)
valsb1, pointsb1 = append(valsb1, zeros(1)), append(pointsb1, pointsb1[-1:])

valsa1, binsa1 = histogram(log10(finalMT1/(10.**10)), bins=100, density=True)
pointsa1 = 0.5*(binsa1[:-1]+ binsa1[1:])
cumdistA1 = cumsum(valsa1*diff(binsa1))*100
att1[1].semilogy(pointsa1, cumdistA1, color='red', linewidth=2)
att1[1].grid()
valsa1, pointsa1 = append(zeros(1), valsa1), append(pointsa1[0:1], pointsa1)
valsa1, pointsa1 = append(valsa1, zeros(1)), append(pointsa1, pointsa1[-1:])

valsb2, binsb2 = histogram(log10(initMT2/(10.**10)), bins=100, density=True)
pointsb2 = 0.5*(binsb2[:-1]+ binsb2[1:])
cumdistB2 = cumsum(valsb2*diff(binsb2))*100
att2[1].plot(pointsb2, cumdistB2, color='blue', linewidth=2)
valsb2, pointsb2 = append(zeros(1), valsb2), append(pointsb2[0:1], pointsb2)
valsb2, pointsb2 = append(valsb2, zeros(1)), append(pointsb2, pointsb2[-1:])

valsa2, binsa2 = histogram(log10(finalMT2/(10.**10)), bins=100, density=True)
pointsa2 = 0.5*(binsa2[:-1]+ binsa2[1:])
cumdistA2 = cumsum(valsa2*diff(binsa2))*100
att2[1].semilogy(pointsa2, cumdistA2, color='red', linewidth=2)
att2[1].grid()
valsa2, pointsa2 = append(zeros(1), valsa2), append(pointsa2[0:1], pointsa2)
valsa2, pointsa2 = append(valsa2, zeros(1)), append(pointsa2, pointsa2[-1:])

valsb3, binsb3 = histogram(log10(initMT3/(10.**10)), bins=100, density=True)
pointsb3 = 0.5*(binsb3[:-1]+ binsb3[1:])
cumdistB3 = cumsum(valsb3*diff(binsb3))*100
att3[1].plot(pointsb3, cumdistB3, color='blue', linewidth=2)
valsb3, pointsb3 = append(zeros(1), valsb3), append(pointsb3[0:1], pointsb3)
valsb3, pointsb3 = append(valsb3, zeros(1)), append(pointsb3, pointsb3[-1:])

valsa3, binsa3 = histogram(log10(finalMT3/(10.**10)), bins=100, density=True)
pointsa3 = 0.5*(binsa3[:-1]+ binsa3[1:])
cumdistA3 = cumsum(valsa3*diff(binsa3))*100
att3[1].semilogy(pointsa3, cumdistA3, color='red', linewidth=2)
att3[1].grid()
valsa3, pointsa3 = append(zeros(1), valsa3), append(pointsa3[0:1], pointsa3)
valsa3, pointsa3 = append(valsa3, zeros(1)), append(pointsa3, pointsa3[-1:])

# Plot merger time distributions
att1[0].step(pointsb1, valsb1, linewidth=2, color='blue', label="Pre-encounter")
att1[0].step(pointsa1, valsa1, linewidth=2, color='red', label="Post-encounter")
att1[0].set_xlabel(r"$\log_{10}(\tau_{GW}/10\ Gyr)$", fontsize=24)
att1[0].set_ylabel(r"$p(\tau_{GW})$", fontsize=24)
att1[1].set_xlabel(r"$\log_{10}(\tau_{GW}/10\ Gyr)$", fontsize=24)
att1[1].set_ylabel(r"Cumulative fration [%]", fontsize=24)
att1[0].grid()
att1[1].grid()

# Plot merger time distributions
att2[0].step(pointsb2, valsb2, linewidth=2, color='blue', label="Pre-encounter")
att2[0].step(pointsa2, valsa2, linewidth=2, color='red', label="Post-encounter")
att2[0].set_xlabel(r"$\log_{10}(\tau_{GW}/10\ Gyr)$", fontsize=24)
att2[0].set_ylabel(r"$p(\tau_{GW})$", fontsize=24)
att2[1].set_xlabel(r"$\log_{10}(\tau_{GW}/10\ Gyr)$", fontsize=24)
att2[1].set_ylabel(r"Cumulative fration [%]", fontsize=24)
att2[0].grid()
att2[1].grid()

# Plot merger time distributions
att3[0].step(pointsb3, valsb3, linewidth=2, color='blue', label="Pre-encounter")
att3[0].step(pointsa3, valsa3, linewidth=2, color='red', label="Post-encounter")
att3[0].set_xlabel(r"$\log_{10}(\tau_{GW}/10\ Gyr)$", fontsize=24)
att3[0].set_ylabel(r"$p(\tau_{GW})$", fontsize=24)
att3[1].set_xlabel(r"$\log_{10}(\tau_{GW}/10\ Gyr)$", fontsize=24)
att3[1].set_ylabel(r"Cumulative fration [%]", fontsize=24)
att3[0].grid()
att3[1].grid()


"""
#Eccentricity change distribution
vals, bins = histogram((data[:,9]-data[:,8]), bins=100, density=True)
points = 0.5*(bins[:-1]+bins[1:])
points = append(points,array([points[-1]]))
points = append(array([points[0]]), points)
vals = append(vals, zeros(1))
vals = append(zeros(1),vals)
# Plot the distributions
aee[0].step(points, vals, linewidth=2)
aee[1].scatter(data[:,8], (data[:,9]-data[:,8]), s=0.1)
aee[0].set_xlabel(r"$e-e_{0}$", fontsize=24)
aee[0].set_ylabel(r"$p(e-e_{0})$", fontsize=24)
aee[0].set_xlabel(r"$e-e_{0}$", fontsize=24)
aee[0].set_ylabel(r"$p(e-e_{0})$", fontsize=24)
"""
"""
# Ratio of initial and final eccentricities
re = data[:,9]/data[:,8]
binse = logspace(log10(re.min()), log10(re.max()), num=100, endpoint=True, base=10.0)
vals, bins = histogram(re, bins=binse, density=True)
#vals, bins = histogram(log10(data[:,9]/data[:,8]), bins=100, density=True)
points = 0.5*(bins[:-1]+bins[1:])
are[0].hist2d(log10(re), log10(mt), bins=100, norm=LogNorm())
are[1].loglog(points, cumsum(vals*diff(bins))*100, linewidth=2)
are[0].set_xlabel(r"$\log_{10}(e/e_{0})$", fontsize=24)
are[0].set_ylabel(r"$\log_{10}(\tau_{GW, out}/\log_{10}(\tau_{GW, in})$", fontsize=24)
are[1].set_xlabel(r"$\log_{10}(e/e_{0})$", fontsize=24)
are[1].set_ylabel(r"Cumulative fraction [%]", fontsize=24)
are[1].set_xlim([0.005, 21.])
are[1].set_ylim([3.e-4, 99])
are[1].grid()
# Effective spins
fil1 = log10(finalMT)<10
fil2 = finalMT[fil1]<initMT[fil1]
vals, bins = histogram(data[:,-10], bins=100, density=True)
valsQ, binsQ = histogram(data[fil1][fil2][:,-10], bins=100, density=True)
points = 0.5*(bins[:-1]+bins[1:])
pointsQ = 0.5*(binsQ[:-1]+binsQ[1:])
cumDist = cumsum(vals*diff(bins))*100
cumDistQ = cumsum(valsQ*diff(binsQ))*100
aef1.step(points, vals, linewidth=2, color='blue', label='Total population')
aef1.step(pointsQ, valsQ, linewidth=2, color='red', label=r'$\tau_{GW,out}< 10\ Gyr$')
aef2.semilogy(points, cumDist, linewidth=2, color='blue', linestyle='--')
aef2.semilogy(pointsQ, cumDistQ, linewidth=2, color='red', linestyle='--')
aef1.set_xlabel(r"$\chi_{eff}$", fontsize=24)
aef1.set_ylabel(r"$p(\chi_{eff})$", fontsize=24)
aef2.set_ylabel("Cumulative fraction [%]", fontsize=24)
aef1.set_xlim([-0.99,0.99])
aef1.set_ylim([0., 5.4])
aef2.set_ylim([0., 100.])
aef1.grid()
aef2.grid()
"""

"""
att3[1].set_ylim([0.215,100])
att3[0].set_xlim([-3,8])
att3[1].set_xlim([-3,8])
att3[1].set_ylim([0.049,100])
att3[1].grid()
att1[1].set_ylim([0.08,100])
att1[1].set_xlim([-5,10.1])
att1[0].set_xlim([-5,10.1])
att1[1].grid()
att1[0].set_ylim([0.,0.408])
att2[0].set_ylim([0.,0.304])
att2[0].set_xlim([-7.9,10])
att2[1].set_xlim([-7.9,10])
att2[1].set_ylim([0.005,100])
att2[1].grid()
"""
"""

"""
# Change in orbital energy
ratE1 = zeros(nn*numS1)
ratE2 = zeros(nn*numS2)
ratE3 = zeros(nn*numS3)
fil1 = log10(finalMT1)<= 10
fil2 = log10(finalMT2)<= 10
fil3 = log10(finalMT3)<= 10

for ii in range(nn):
    ratE1[ii*numS1:(ii+1)*numS1] = (initA1[ii*numS1:(ii+1)*numS1]/finalA1[ii*numS1:(ii+1)*numS1] - 1.) * 100.
    ratE2[ii*numS2:(ii+1)*numS2] = (initA2[ii*numS2:(ii+1)*numS2]/finalA2[ii*numS2:(ii+1)*numS2] - 1.) * 100.
    ratE3[ii*numS3:(ii+1)*numS3] = (initA3[ii*numS3:(ii+1)*numS3]/finalA3[ii*numS3:(ii+1)*numS3] - 1.) * 100.

ratE1 = ratE1[fil1]
ratE2 = ratE2[fil2]
ratE3 = ratE3[fil3]


vals1, bins1 = histogram(ratE1[ratE1<300], bins=40, density=True)
points1 = 0.5*(bins1[:-1]+bins1[1:])
cumDist1 = cumsum(vals1*diff(bins1)) * 100

vals2, bins2 = histogram(ratE2[ratE2<200], bins=40, density=True)
points2 = 0.5*(bins2[:-1]+bins2[1:])
cumDist2 = cumsum(vals2*diff(bins2)) * 100

vals3, bins3 = histogram(ratE3[ratE3<200], bins=40, density=True)
points3 = 0.5*(bins3[:-1]+bins3[1:])
cumDist3 = cumsum(vals3*diff(bins3)) * 100


adelE0[0].step(points1, vals1, linewidth=2, label=r'$p(e_{0}) \propto const.$')
adelE0[0].step(points2, vals2, linewidth=2, label=r'$p(e_{0}) \propto e_{0}$')
adelE0[0].step(points3, vals3, linewidth=2, label=r'$p(e_{0}) \propto e_{0}^{-1/2}$')
adelE0[1].plot(points1, cumDist1)
adelE0[1].plot(points2, cumDist2)
adelE0[1].plot(points3, cumDist3)


# Energy radiated due to gravitational waves
# deltaE_encounter(MM, m1, m2, rm)
binEn01 = binary_energy(m1, m2, initA1[fil1])
binEn1  = binary_energy(m1, m2, finalA1[fil1])
delGW1 = deltaE_encounter(MM, m1, m2, rps1[fil1])
vals1, bins1 = histogram(log10(delGW1/binEn01), bins=40, density=True)
points1 = 0.5*(bins1[:-1]+bins1[1:])
cumDist1 = cumsum(vals1*diff(bins1)) * 100

binEn02 = binary_energy(m1, m2, initA2[fil2])
binEn2  = binary_energy(m1, m2, finalA2[fil2])
delGW2 = deltaE_encounter(MM, m1, m2, rps2[fil2])
vals2, bins2 = histogram(log10(delGW2/binEn02), bins=40, density=True)
points2 = 0.5*(bins2[:-1]+bins2[1:])
cumDist2 = cumsum(vals2*diff(bins2)) * 100

binEn03 = binary_energy(m1, m2, initA3[fil3])
binEn3  = binary_energy(m1, m2, finalA3[fil3])
delGW3 = deltaE_encounter(MM, m1, m2, rps3[fil3])
vals3, bins3 = histogram(log10(delGW3/binEn03), bins=40, density=True)
points3 = 0.5*(bins3[:-1]+bins3[1:])
cumDist3 = cumsum(vals3*diff(bins3)) * 100

adelGE[0].step(points1, vals1, linewidth=2, label=r'$p(e_{0}) \propto const.$')
adelGE[0].step(points2, vals2, linewidth=2, label=r'$p(e_{0}) \propto e_{0}$')
adelGE[0].step(points3, vals3, linewidth=2, label=r'$p(e_{0}) \propto e_{0}^{-1/2}$')
adelGE[1].plot(points1, cumDist1)
adelGE[1].plot(points2, cumDist2)
adelGE[1].plot(points3, cumDist3)



ratr1 = delGW1/(binEn1-binEn01)
ratr1 = ratr1[ratr1<.25]
ratr1 = ratr1[ratr1>-.25]
ratr2 = delGW2/(binEn2-binEn02)
ratr2 = ratr2[ratr2<.25]
ratr1 = ratr2[ratr2>-.25]
ratr3 = delGW3/(binEn3-binEn03)
ratr3 = ratr3[ratr3<.25]
ratr1 = ratr3[ratr3>-.25]

vals1, bins1 = histogram(ratr1, bins=40)
points1 = 0.5*(bins1[:-1]+bins1[1:])
cumDist1 = cumsum(vals1*diff(bins1)) * 100
vals2, bins2 = histogram(ratr2, bins=40)
points2 = 0.5*(bins2[:-1]+bins2[1:])
cumDist2 = cumsum(vals2*diff(bins2)) * 100
vals3, bins3 = histogram(ratr3, bins=40)
points3 = 0.5*(bins3[:-1]+bins3[1:])
cumDist3 = cumsum(vals3*diff(bins3)) * 100


aE0GE[0].step(points1, vals1, linewidth=2, label=r'$p(e_{0}) \propto const.$')
aE0GE[0].step(points2, vals2, linewidth=2, label=r'$p(e_{0}) \propto e_{0}$')
aE0GE[0].step(points3, vals3, linewidth=2, label=r'$p(e_{0}) \propto e_{0}^{-1/2}$')
aE0GE[1].plot(points1, cumDist1)
aE0GE[1].plot(points2, cumDist2)
aE0GE[1].plot(points3, cumDist3)

# Calculate post-encoutner COM-MBH orbital eccentricity distribution

cosf1 = cos(data1[:,22])
sinf1 = sin(data1[:,22])
cosf2 = cos(data2[:,22])
sinf2 = sin(data2[:,22])
cosf3 = cos(data3[:,22])
sinf3 = sin(data3[:,22])

binEn01 = binary_energy(m1, m2, initA1)
binEn1  = binary_energy(m1, m2, finalA1)
binEn02 = binary_energy(m1, m2, initA2)
binEn2  = binary_energy(m1, m2, finalA2)
binEn03 = binary_energy(m1, m2, initA3)
binEn3  = binary_energy(m1, m2, finalA3)

for ii in range(nn):

    rmod1 = rps1[ii*numS1:(ii+1)*numS1] /(1. + cosf1)
    vmod1 = rps1[ii*numS1:(ii+1)*numS1] * sqrt(2.)/4. * (1. + cosf1)**2. / ts1[ii*numS1:(ii+1)*numS1]

    rmod2 = rps2[ii*numS2:(ii+1)*numS2] /(1. + cosf2)
    vmod2 = rps2[ii*numS2:(ii+1)*numS2] * sqrt(2.)/4. * (1. + cosf2)**2. / ts2[ii*numS2:(ii+1)*numS2]

    rmod3 = rps3[ii*numS3:(ii+1)*numS3] /(1. + cosf3)
    vmod3 = rps3[ii*numS3:(ii+1)*numS3] * sqrt(2.)/4. * (1. + cosf3)**2. / ts3[ii*numS3:(ii+1)*numS3]

    l12[ii*numS1:(ii+1)*numS1] = rmod1**2. * vmod1**2 * (
                                     cosf1/(1.+cosf1) + cosf1*sinf1/((1+cosf1)**2.))**2.

    l22[ii*numS2:(ii+1)*numS2] = rmod2**2. * vmod2**2 * (
                                     cosf2/(1.+cosf2) + cosf2*sinf2/((1+cosf2)**2.))**2.

    l32[ii*numS3:(ii+1)*numS3] = rmod3**2. * vmod3**2 * (
                                     cosf3/(1.+cosf3) + cosf3*sinf3/((1+cosf3)**2.))**2.



ee_t1 = sqrt(1. -((binEn1-binEn01)*l12)/(GG**2. * MM**2. * (m1+m2)) )
ee_t2 = sqrt(1. -((binEn2-binEn02)*l22)/(GG**2. * MM**2. * (m1+m2)) )
ee_t3 = sqrt(1. -((binEn3-binEn03)*l32)/(GG**2. * MM**2. * (m1+m2)) )


print "Dist 1., percentage of mergers", shape(finalMT1[log10(finalMT1)<10]
                                            )[0]/float(shape(finalMT1)[0]) * 100

print "Dist 2., percentage of mergers", shape(finalMT2[log10(finalMT2)<10]
                                            )[0]/float(shape(finalMT2)[0]) * 100

print "Dist 3., percentage of mergers", shape(finalMT3[log10(finalMT3)<10]
                                            )[0]/float(shape(finalMT3)[0]) * 100
