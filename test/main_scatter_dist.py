#! /usr/bin/python

import sys
import math
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import Integrator
import Detector
import MultipleScatter
import Params

#Detector.LoadBField("bfield/bfield.txt")
Params.BFieldType = 'none'
Params.Q = 1.0
Params.MSCtype = 'PDG'
Params.MatSetup = 'iron'

prefix = "data/scatter_dist/"

dt = 0.2
nsteps = 100

#pvals = np.array([3000,5000,7000,10000])
#pvals = np.array([3000])
pvals = np.arange(3000,20001, 1000.)

L = 2.
Nsamp = 1000

rvals={}
thvals={}

for p in pvals:

    print p

    p0 = [p,0,0]
    x0 = np.array([0,0,0]+p0)

    rvals[p] = []
    thvals[p] = []
    
    ntuple = ROOT.TNtuple("p"+str(int(p))+"_data","title","r:theta")

    for i in range(Nsamp):

        traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=L, cutoffaxis=0)
        
        # get the coordinates as it leaves material
        k = (L-traj[0,-2])/(traj[0,-1]-traj[0,-2])
        cross = traj[:3,-2]+k*(traj[:3,-1]-traj[:3,-2])

        rvals[p].append(np.sqrt(cross[1]**2 + cross[2]**2)*100)
        
        # get the projected angle as it exits material
        projth = np.arctan2(traj[4,-1],traj[3,-1])
        thvals[p].append(projth*180/np.pi)

        ntuple.Fill(rvals[p][-1],thvals[p][-1])

    # rfile = ROOT.TFile(prefix+"p"+str(int(p))+".root","RECREATE")
    # ntuple.Write()
    # rfile.Close()

plt.figure(1)

plt.hist(rvals[3000], bins=70, range=(0,25), histtype='step', color='r', label="3 GeV")
plt.hist(rvals[5000], bins=70, range=(0,25), histtype='step', color='g', label="5 GeV")
plt.hist(rvals[7000], bins=70, range=(0,25), histtype='step', color='b', label="7 GeV")
plt.hist(rvals[10000],bins=70, range=(0,25), histtype='step', color='c', label="10 GeV")

plt.legend()
plt.xlabel("Deviation (cm)")
plt.title("Spatial deviation through 2 m of iron")

plt.figure(2)

plt.hist(thvals[3000], bins=70, range=(-10,10), histtype='step', color='r', label="3 GeV")
plt.hist(thvals[5000], bins=70, range=(-10,10), histtype='step', color='g', label="5 GeV")
plt.hist(thvals[7000], bins=70, range=(-10,10), histtype='step', color='b', label="7 GeV")
plt.hist(thvals[10000],bins=70, range=(-10,10), histtype='step', color='c', label="10 GeV")

plt.legend()
plt.xlabel("Angular deflection (deg)")
plt.title("Angular deflection through 2 m of iron")

rmeans = []
rstd = []
thstd = []
for p in pvals:
    rmeans.append(np.mean(rvals[p]))
    rstd.append(np.std(rvals[p])/np.sqrt(Nsamp))
    thstd.append(np.std(thvals[p]))

plt.figure(3)
plt.errorbar(pvals/1000.,rmeans, yerr=rstd, fmt='-or')

plt.xlabel("p (GeV)")
plt.ylabel("Deviation (cm)")
plt.title("Mean Deviation through 2m of iron")
plt.gca().set_ylim(ymin=0)

plt.figure(4)
plt.plot(pvals/1000.,thstd, '-or')

plt.xlabel("p (GeV)")
plt.ylabel("Width (deg)")
plt.title("Width of theta distribution")
plt.gca().set_ylim(ymin=0)

plt.show()
