from physics import *
from matplotlib.pylab import *
ion()

data = loadtxt("/home/joe/Documents/PhD_datafiles/eccentricity/cEudD/survE01.dat")
print shape(data)[0]
data = data[data[:,10]>0]
print shape(data)[0]

m1 = 30*Msun
m2 = 30*Msun
beta = betaFactor(m1, m2)

MM = 4.e6*Msun

minA = 1.e-2*au
maxA = (1.)*au
numS = shape(data)[0]



initA   = uniformLog(minA, maxA, numS)
initMT  = mergerTimeecc(initA, data[:,8], beta)
initMT  = initMT/sTy
finalMT = mergerTimeecc(data[:,0]*data[:,10]*initA, data[:,9], beta)
finalMT = finalMT/sTy


def plots(data, initMT, finalMT):
    # data -> array object from data file

    fe, ae     = subplots(1,2)
    fa, aa     = subplots(1,2)
    ft, at     = subplots(1,2)
    ftt, att   = subplots(1,2)
    #fef, aef   = subplots(1,2)
    fee, aee   = subplots(1,2)
    #fte, ate   = subplots(1,2)
    fre, are   = subplots(1,2)

    # Calculate eccentricity distribution
    vals, bins = histogram(1.- data[:,9], bins=40, density=True)
    points = 0.5*(bins[:-1]+bins[1:])
    cumDist = cumsum(vals*diff(bins))*100
    # Plot eccentricity distribution
    ae[0].step(points, vals, linewidth=2)
    ae[1].loglog(points, cumDist, linewidth=2)
    ae[0].set_xlabel(r"$1-e$", fontsize=24)
    ae[0].set_ylabel(r"$p(1-e)$", fontsize=24)
    ae[1].set_xlabel(r"$1-e$", fontsize=24)
    ae[1].set_ylabel(r"Cumulative fraction $[\%]$", family='sans-serif', fontsize=24)

    # Calculate semi-major axis distribution
    aas = data[:,10]*data[:,0]
    binse = logspace(log10(aas.min()), log10(aas.max()), num=40, endpoint=True, base=10.0)
    vals, bins = histogram(aas, bins=binse, density=True)
    points = 0.5*(bins[:-1]+bins[1:])
    cumDist = cumsum(vals*diff(bins))*100
    # Plot semi-major axis distribution
    aa[0].step(points, vals, linewidth=2)
    aa[1].loglog(points, cumDist, linewidth=2)
    aa[0].set_xlabel(r"$a/a_{0}$", fontsize=24)
    aa[0].set_ylabel(r"$p(a/a_{0})$", fontsize=24)
    aa[1].set_xlabel(r"$a/a_{0}$", fontsize=24)
    aa[1].set_ylabel(r"Cumulative fraction $[\%]$", fontsize=24)

    # Calculate merger time ratio distribution
    mt = finalMT/initMT
    binse = logspace(log10(mt.min()), log10(mt.max()), num=100, endpoint=True, base=10.0)
    vals, bins = histogram(mt, bins=binse, density=True)
    points = 0.5*(bins[:-1]+ bins[1:])
    cumDist = cumsum(vals*diff(bins))*100
    vals, bins = histogram(log10(mt), bins=100, density=True)
    points2 = 0.5*(bins[:-1]+ bins[1:])
    # Plot merger time ratio distribution
    at[0].step(points2, vals, linewidth=2)
    at[1].loglog(points, cumDist, linewidth=2)
    at[0].set_xlabel(r"$\log_{10}(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)
    at[0].set_ylabel(r"$p(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)
    at[1].set_xlabel(r"$\tau_{GW, out}/\tau_{GW, in}$", fontsize=24)
    at[1].set_ylabel(r"Cumulative fraction $[\%]$", fontsize=24)

    # Initial and final merger time distributions
    valsb, binsb = histogram(log10(initMT), bins=100, density=True)
    pointsb = 0.5*(binsb[:-1]+ binsb[1:])
    valsb, pointsb = append(zeros(1), valsb), append(pointsb[0:1], pointsb)
    valsb, pointsb = append(valsb, zeros(1)), append(pointsb, pointsb[-1:])
    valsa, binsa = histogram(log10(finalMT), bins=100, density=True)
    pointsa = 0.5*(binsa[:-1]+ binsa[1:])
    valsa, pointsa = append(zeros(1), valsa), append(pointsa[0:1], pointsa)
    valsa, pointsa = append(valsa, zeros(1)), append(pointsa, pointsa[-1:])
    # Plot merger time distributions
    att[0].step(pointsb, valsb, linewidth=2, color='blue', label="Pre-encounter")
    att[0].step(pointsa, valsa, linewidth=2, color='red', label="Post-encounter")
    att[1].scatter(log10(initMT), log10(mt), s=0.1)
    att[0].set_xlabel(r"$\tau_{GW}$", fontsize=24)
    att[0].set_ylabel(r"$p(\tau_{GW})$", fontsize=24)
    att[1].set_xlabel(r"$\log_{10}(\tau_{GW, in})$", fontsize=24)
    att[1].set_ylabel(r"$\log_{10}(\tau_{GW, out}/\tau_{GW, in})$", fontsize=24)


    #Eccentricity change distribution
    vals, bins = histogram((data[:,9]-data[:,8]), bins=40, density=True)
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

    # Ratio of initial and final eccentricities
    vals, bins = histogram((data[:,9]/data[:,8]), bins=40, density=True)
    points = 0.5*(bins[:-1]+bins[1:])
    points = append(points,array([points[-1]]))
    points = append(array([points[0]]), points)
    vals = append(vals, zeros(1))
    vals = append(zeros(1),vals)
    are[0].scatter(log10(data[:,9]/data[:,8]), log10(finalMT), s=0.1)
    are[1].loglog(points, cumsum(vals)*100*(points[2]-points[0]), linewidth=2)
    are[0].set_xlabel(r"$\log_{10}(e/e_{0})$", fontsize=24)
    are[0].set_ylabel(r"$\log_{10}(\tau_{GW, out})$", fontsize=24)
    are[1].set_xlabel(r"$\log_{10}(e/e_{0})$", fontsize=24)
    are[1].set_ylabel(r"Cumulative fraction [%]", fontsize=24)
