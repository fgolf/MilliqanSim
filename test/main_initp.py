#! /usr/bin/env python

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
scatterType = 'PDG'
Params.MSCtype = scatterType
Params.MatSetup = 'cms'
Params.EnergyLossOn=True

#print "\nDone loading B-field...\nB-field at center:", Detector.getBField(0,0,0), '\n'


# define the initial momenta (in MeV)
init_p = []

init_p.append([2000,0,0])
init_p.append([3000,0,0])
init_p.append([5000,0,0])
init_p.append([10000,0,0])
init_p.append([20000,0,0])

colors = ['r', 'g', 'b', 'c', 'm']

print 'Initial Momenta (colors r,g,b,c):'
for i in range(len(init_p)):
    print ' -',round(np.linalg.norm(init_p[i])/1000, 1), "GeV"

dt = 0.1
nsteps = 700

trajs = [None for i in init_p]

#compute trajectories
for i in range(len(init_p)):
    x0 = np.array([0,0,0]+init_p[i])
    trajs[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)

fig = plt.figure(2, figsize=(14.5,5.5))
## xz slice
x = np.arange(-Params.RMAX, Params.RMAX+1e-10, Params.DR)/100
z = np.arange(Params.ZMIN, Params.ZMAX+1e-10, Params.DZ)/100
    
Z,X = np.meshgrid(z,x)

plt.subplot2grid((1,5),(0,0),colspan=3)

# draw mag field

if Params.BFieldLoaded:
    mag = np.append(Params.Bmag[::-1,:,0],Params.Bmag[1:,:,180/Params.DPHI],0)
    bmplot = plt.pcolor(Z,X,mag,cmap='afmhot',vmax=5.0)
    #bmcb = plt.colorbar(bmplot, orientation='horizontal')
    

sl = Params.solLength
sr = Params.solRad

# draw trajectory
for i in range(len(init_p)):
    plt.plot(trajs[i][2,:],trajs[i][0,:],'-', linewidth=2, color=colors[i])

plt.axis([-15,15,-9,9])
plt.xlabel("z (m)")
plt.ylabel("x (m)")

## xy slice

plt.subplot2grid((1,5),(0,3), colspan=2)

Drawing.DrawXYslice(trajs,ax=plt.gca())

plt.axis([-9,9,-9,9])
plt.xlabel("x (m)")
plt.ylabel("y (m)")

## 3d view

fig = plt.figure(num=3, figsize=(8, 8))

Drawing.Draw3Dtrajs(trajs)


plt.show()



