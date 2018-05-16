from physics import *
from mpi4py import MPI

#Init()

comm   = MPI.COMM_WORLD
rank   = comm.Get_rank()
nprocs = comm.Get_size()
a0_min = 0.01*au
a0_max = 1.0*au
nS = 2                         # Number of rescalings

Mmbh = 4.e6*Msun
m1   = 30.*Msun
m2   = 30.*Msun
mm   = m1 + m2
Mm  = (Mmbh/mm)**(1./3.)
rK  = 1.*pc
orbT = 2.*pi/sqrt(GG*Mmbh) * rK**(3./2.)                    # Full elliptical orbit time
beta = betaFactor(m1, m2)

# master proc load and share data with other procs
if rank == 0:
    print "Proc 0 reading data file"
    dataT = loadtxt("/home/joe/Documents/PhD_datafiles/UD_10rt_Dea.dat")     # Change data.dat for an adequate directory
    dataT = dataT[dataT[:,2]>0]
    nE = shape(dataT)[0]            # Number of binary encounters
else:
    nE = 0

nE = comm.bcast(nE, root=0)

comm.Barrier()
modN = nE%nprocs
if modN == 0:
    #print "case 1"
    nEP = nE/nprocs
    #print rank, nEP
    #print "rank" + str(rank) + str(nE) + str(nEP)

else:
    if rank < modN:
        nEP = nE/nprocs + 1
        #print "rank" + str(rank) + str(nE) + str(nEP)
    else:
        nEP = nE/nprocs
dataP = zeros([nEP,4])              # Prepare empty arrays to recieve data
#print rank, nEP, shape(dataP)


if rank == 0:
    if nE%nprocs != 0:
        dataP = dataT[rank*(nE/nprocs + 1):(rank+1)*(nE/nprocs + 1),:]
        for ii in range(1, nprocs-1):
            if ii < modN:
            #print dataT[rank*(nE/nprocs + 1):(rank+1)*(nE/nprocs + 1),:]
            #print shape(dataT[(rank+1)*(nE/nprocs + 1):(rank+2)*(nE/nprocs + 1),:])
                comm.Send([dataT[(ii)*(nE/nprocs + 1):(ii+1)*(nE/nprocs + 1),:
                                                                ], MPI.DOUBLE], dest=ii, tag=(77+ii))
            else:
                comm.Send([dataT[(ii*(nE/nprocs)+nE%nprocs):((ii+1)*(nE/nprocs)+nE%nprocs),:
                                                                ], MPI.DOUBLE], dest=ii, tag=(77+ii))

            comm.Send([dataT[(nprocs-1)*(nE/nprocs) + modN:,:], MPI.DOUBLE], dest=(nprocs-1), tag=(77+nprocs-1))
    else:
        dataP = dataT[:nEP,:]
        for ii in range(1, nprocs-1):
            comm.Send([dataT[ii*nEP:(ii+1)*nEP,:], MPI.DOUBLE], dest=ii, tag=(77+ii))

        comm.Send([dataT[nE-nEP:,:], MPI.DOUBLE], dest=(nprocs-1), tag=(77+nprocs-1))

else:
    comm.Recv([dataP, MPI.DOUBLE], source=0, tag=(rank+77))
    print "Proc " + str(rank) + " recieved data."
#print rank, nEP, dataP[1:10,1]
#comm.Barrier()
#if rank == 0: del dataT
# Generate binary distribution at Apoastron

# Pre-calcutiations
print "Proc. "+ str(rank) + " calculating integrl terms of merger times"
mtInt= quadArray(mergeTimeInt, zeros(nEP), dataP[:,1])         # Integral part of merger time, only depends on eccentricity. This line reduces the number of integrals to be calculated
ff0  = arccos(dataP[:,0]/5. - 1)                               # ff0 -> doesn't depend on scaled values
tf0  = tan(ff0/2.)                               # tan(ff0)
tt0  = (sqrt(2.)/3.)*tf0*(3. + tf0**2.)          # Scaled initial time. Also time at which the binary reaches 10 Rt after the encounter
nB    = nS*nEP                                   # Number of trajectories this process will generate
dataP[:,2] = dataP[:,0]*dataP[:,2]                  # Post-encounter semi-major axes in units of initial semi-major axis
results = zeros([nB, 6])

print "Proc. "+ str(rank) + " sampling binary separations"
for ii in range(nS):
    dists = zeros(nEP)
    tmf   = zeros(nEP)
    mt    = zeros(nEP)
    aa0s = uniformLog(a0_min, a0_max, nEP)       # Binary separations at Apoastron
    Rt   = Mm*aa0s                               # Tidal radius
    Rp   = dataP[:,0]*Rt                         # Periastron distance: closest approach
    Ls   = Rp/Mmbh                               # Problem length scale
    Ts   = sqrt(Rp**3./(GG*Mmbh))                 # Problem time scale
    print rank, shape(Ts)
#    aa0s = (aa0s**4. - 4*beta*orbT/2.)**(1./4.)  # Hardening ue to GW radiation on transit from apoastron to periastron
    aa   = aa0s*dataP[:,2]                        # Post-encounter separation. Evaluated at 10Rt

    cc0s = cc0Factor(aa, dataP[:,1])               # Integration constant from Peter's equations, depends on pre-merger e and sma
    mt   = (12./19.) * cc0s**4. / beta * mtInt   # Merger time for eccentric binaries, as given by Peter's equations
    tt   = mt + Ts*tt0                           # Total time after encounter at which the binary will have merged.
    polCoefs = zeros([nEP, 4])
    polCoefs[:,0] = 1.
    polCoefs[:,2] = 3.
    polCoefs[:,3] = -1.*3./(Ts*sqrt(2.)) * tt
    for jj in range(nEP):
        tmfU = roots(polCoefs[jj,:])
        tmf[jj] = tmfU.real[abs(tmfU.imag)<1e-6]
    mf  = 2.*arctan(tmf)
    dists = 2*Rp/(1.+cos(mf))

    results[ii*nEP:(ii+1)*nEP,0] = dataP[:,0]     # D
    results[ii*nEP:(ii+1)*nEP,1] = aa0s           # Pre-encounter sma
    results[ii*nEP:(ii+1)*nEP,2] = aa             # Post-encounter SMA
    results[ii*nEP:(ii+1)*nEP,3] = dataP[:,1]     # Eccentricity
    results[ii*nEP:(ii+1)*nEP,4] = mt             # Merger time
    results[ii*nEP:(ii+1)*nEP,5] = dists[:]       # Distance of merger

fname = "outputDist/results_" + str(rank) + ".dat"
savetxt(fname, results)
