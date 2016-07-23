#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import Params
import Integrator
import Detector

Detector.LoadCoarseBField("../bfield/bfield_coarse.pkl")

Params.BFieldOn = True
Params.BFieldType = 'cms'
Params.MSCtype = 'none'
Params.Q = 1.0
Params.m = 105.

dt = 0.1
nsteps = 5000
thvals = []
bvals = []

pvals = np.arange(2000.,20001.,100.)

nth = 30.

for p in pvals:
    thvals.append(0)
    bvals.append(0)
    for thi in np.arange(0,360-360/nth+1,360/nth):
        p0 = [p*np.cos(thi),p*np.sin(thi),0]
        x0 = np.array([0,0,0]+p0)
        traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=9, cutoffaxis=3)
        
        pf = traj[3:,-1]
        xf = traj[:3,-1]
        
        t = -np.dot(pf,xf)/np.dot(pf,pf)
        b = np.linalg.norm(xf+t*pf)
        bvals[-1] += b
        
        th = abs( np.arccos(np.dot(pf,p0)/np.linalg.norm(p0)/np.linalg.norm(pf)) * 180/np.pi )
        thvals[-1] += th
    thvals[-1] /= nth
    bvals[-1] /= nth

print thvals[-1]
        
plt.figure(1)
plt.plot(pvals/1000,thvals)
plt.xlabel('p (GeV)')
plt.ylabel(r'$\Delta\theta$ (deg)')
plt.title(r'$\Delta\theta$ vs. initial $p$, $\eta=0$, q=1.0e')

plt.figure(2)
plt.plot(pvals/1000,bvals)
plt.xlabel('p (GeV)')
plt.ylabel('Impact parameter (m)')
plt.title(r'Impact parameter vs. initial $p$, $\eta=0$, q=1.0e')

plt.show()
