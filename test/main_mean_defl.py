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
Params.BFieldOn = False
Params.Q = 0.5

isBarrell = True

prefix = "data/barrell_Boff_qhalf/"

dt = 0.2
nsteps = 5000

pvals = np.arange(2000,20001,1000)

meanrvals = []
sdevrvals = []

Nsamp = 100

#function to get r value of a point
getr = lambda x: np.sqrt(x[0]**2+x[1]**2)

if isBarrell:
    co = 8
    coa = 3
else:
    co = 14
    coa = 2

for p in pvals:

    print p

    if isBarrell:
        p0 = [p,0,0]
    else:
        th = np.arctan(1./4)
        p0 = [p*np.sin(th)*np.cos(np.pi/4),p*np.sin(th)*np.sin(np.pi/4),p*np.cos(th)]
    
    x0 = np.array([0,0,0]+p0)

    rvals = []
    thvals = []

    ## get the trajectory & end coords with no multiple scattering
    Params.MSCtype = 'none'
    traj_noMSC = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=co, cutoffaxis=coa)
    if isBarrell:
        k = (8-getr(traj_noMSC[:,-2]))/(getr(traj_noMSC[:,-1])-getr(traj_noMSC[:,-2]))
    else:
        k = (14-traj_noMSC[2,-2])/(traj_noMSC[2,-1]-traj_noMSC[2,-2])
    cross_noMSC = traj_noMSC[:3,-2]+k*(traj_noMSC[:3,-1]-traj_noMSC[:3,-2])
    projth_noMSC = np.arctan2(traj_noMSC[1,-1],traj_noMSC[0,-1])
    
    Params.MSCtype = 'kuhn'

    ntuple = ROOT.TNtuple("p"+str(p)+"_data","title","r:theta")

    for i in range(Nsamp):

        traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=co, cutoffaxis=coa)
        
        # get the coordinates as it leaves material
        if isBarrell:
            k = (8-getr(traj[:,-2]))/(getr(traj[:,-1])-getr(traj[:,-2]))
        else:
            k = (14-traj_noMSC[2,-2])/(traj_noMSC[2,-1]-traj_noMSC[2,-2])
        cross = traj[:3,-2]+k*(traj[:3,-1]-traj[:3,-2])

        if isBarrell:
            rvals.append(np.sqrt((cross[1]-cross_noMSC[1])**2 + (cross[2]-cross_noMSC[2])**2))
        else:
            rvals.append(np.sqrt((cross[1]-cross_noMSC[1])**2 + (cross[0]-cross_noMSC[0])**2))
        
        # get the projected angle as it exits material
        projth = np.arctan2(traj[1,-1],traj[0,-1])
        thvals.append(abs(projth-projth_noMSC))

        ntuple.Fill(rvals[-1],thvals[-1])

    fid = open(prefix+"p"+str(p),"w")
    for i in range(Nsamp):
        fid.write("{0}\t{1}\n".format(rvals[i],thvals[i]))
    fid.close()

    rfile = ROOT.TFile(prefix+"p"+str(p)+".root","RECREATE")
    ntuple.Write()
    rfile.Close()

    meanrvals.append(np.mean(rvals))
    sdevrvals.append(np.std(rvals))


meanrvals = np.array(meanrvals)
sdevrvals = np.array(sdevrvals)

fid = open(prefix+"mean_defl","w")
for i in range(len(pvals)):
    fid.write("{0}\t{1}\t{2}\n".format(pvals[i],meanrvals[i],sdevrvals[i]))
fid.close()

#plt.errorbar(pvals/1000, meanrvals, yerr=sdevrvals/np.sqrt(Nsamp), fmt='ok')
#plt.errorbar(pvals/1000, meanyvals, yerr=sdevyvals/np.sqrt(Nsamp), fmt='or')

#plt.show()
