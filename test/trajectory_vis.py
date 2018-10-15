#! /usr/bin/env python

## Computes trajectories with a few initial conditions
## and makes plots for visualization. Makes 3 plots:
##   -"XZ" slice - looking at CMS from the side. CMS magnetic field overlayed
##   -"XY" slice - looking at CMS from the endcap. Get classic S-shape trajectories
##   - 3D view - a rotatable 3D visualization of the trajectories

import sys
try:
    sys.path.remove('/home/users/bemarsh/.local/lib/python2.7/site-packages/matplotlib-1.4.3-py2.7-linux-x86_64.egg')
except:
    pass
import numpy as np
import matplotlib.pyplot as plt
import math
import Params
import Integrator
import Detector
import Drawing

## this must be present at the beginning of the main script
## it unpickles the object and stores the array in memory
Detector.LoadCoarseBField("../bfield/bfield_coarse.pkl")

# this tells it to use the CMS magnetic field
Params.BFieldType = 'cms'
# charge of the particle
Params.Q = 1.0
# use PDG scattering. This should probably always be the case
scatterType = 'PDG'
Params.MSCtype = scatterType
# use our simple model of CMS for the material
Params.MatSetup = 'cms'
# turn on dE/dx energy loss
Params.EnergyLossOn=True

# define the initial momenta (in MeV)
# just define 5 for now, all in the +x direction
init_p = []
init_p.append([2000* np.cos( 0*np.pi/180),0,2000* np.sin( 0*np.pi/180)])
init_p.append([3000* np.cos(10*np.pi/180),0,3000* np.sin(10*np.pi/180)])
init_p.append([5000* np.cos(20*np.pi/180),0,5000* np.sin(20*np.pi/180)])
init_p.append([10000*np.cos(30*np.pi/180),0,10000*np.sin(30*np.pi/180)])
init_p.append([20000*np.cos(40*np.pi/180),0,20000*np.sin(40*np.pi/180)])

colors = ['r', 'g', 'b', 'c', 'm']
print 'Initial Momenta (colors r,g,b,c):'
for i in range(len(init_p)):
    print ' -',round(np.linalg.norm(init_p[i])/1000, 1), "GeV"

# set the timestep and number of steps to simulate
dt = 0.1
nsteps = 700

# compute trajectories. Each trajectory is a 2D numpy array,
# where the first dimension contains (x,y,z,px,py,pz),
# and the second dimension is the "timestep"
# e.g. trajs[0][2,5] is the z coordinate at the 6th step of the 1st trajectory
trajs = [None for i in init_p]
for i in range(len(init_p)):
    x0 = np.array([0,0,0]+init_p[i])
    trajs[i] = Integrator.rk4(x0, dt, nsteps)

fig = plt.figure(1, figsize=(14.5,5.5))

## Do the drawing. the Drawing module contains methods to do this automatically

## XZ slice
plt.subplot2grid((1,5),(0,0),colspan=3)
Drawing.DrawXZslice(trajs, colors=colors, ax=plt.gca(), drawBField=True)

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
fig = plt.figure(num=2, figsize=(8, 8))
Drawing.Draw3Dtrajs(trajs)


plt.show()



