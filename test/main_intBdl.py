#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
import Detector
import Params
import Drawing

Detector.LoadCoarseBField("../bfield/bfield_coarse.pkl")

Params.BFieldType = 'cms'
Params.Q = 1.0
Params.MSCtype = 'none'
Params.EnergyLossOn = False
Params.Interpolate = True
Params.UseFineBField = False

z = 0
phi = 0
Bdl=0
Bdlvals = [0]
dr = 0.5
rvals = np.arange(0,901,dr)
for r in rvals:

    igd = Detector.getBField(r/100.0,0,0)[2]*dr/100 * 0.3/20 * 180/np.pi
    Bdl += igd
    Bdlvals.append(Bdl)

print "\nPredicted deflection for 20 GeV particle:",Bdl,"deg\n"

plt.plot(rvals,Bdlvals[:-1],'-b', linewidth=2, label=r"$\int$ $Bdl$ straight line")


x0 = np.array([0,0,0,20000,0,0])
dt = 0.2
nsteps = 1000
traj = Integrator.rk4(x0,Integrator.traverseBField,dt,nsteps,cutoff=9,cutoffaxis=3)

#trajRvals = traj[0,:]*100
trajRvals = np.linalg.norm(traj[:3,:],axis=0)*100
trajTvals = np.absolute(np.arctan(traj[4,:]/traj[3,:]) * 180/np.pi)

Bdlpath = 0
Bdlpathvals = [0]
for i in range(traj.shape[1]-1):
    r1 = traj[:3,i]
    r2 = traj[:3,i+1]
    rav = (r1+r2)/2
    dr = np.linalg.norm(r2-r1)
    igd = Detector.getBField(rav[0],rav[1],rav[2])[2] * dr * 0.3/20 * 180/np.pi
    Bdlpath += igd
    Bdlpathvals.append(Bdlpath)

plt.plot(trajRvals,Bdlpathvals, '-g', linewidth=2, label=r'$\int$ $Bdl$ along path')
plt.plot(trajRvals,trajTvals,'--r', linewidth=2, label="simulated")


plt.legend(loc='lower right')
plt.xlabel('r (cm)')
plt.ylabel(r'$\Delta\theta$')
plt.title(r'$\Delta\theta$ vs. r')

plt.show()
