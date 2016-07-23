# Integrator.py
# contains routines to perform the numeric integration

import numpy as np
import Detector
from MultipleScatter import *
import Params

def doEnergyLoss(x, dt):
    ## returns new x after losing proper amount of energy according to Bethe-Bloch

    p = x[3:]
    magp = np.linalg.norm(p)
    E = np.sqrt(magp**2+Params.m**2)
    gamma = E/Params.m
    beta = magp/E;
    me = 0.511;  #electron mass in MeV
    
    Wmax = 2*me*beta**2*gamma**2/(1+2*gamma*me/Params.m + (me/Params.m)**2)
    K = 0.307075  # in MeV cm^2/mol

    mat = Detector.getMaterial(x[0],x[1],x[2])
    Z,A,rho,X0 = Params.materials[mat]
    I,a,k,x0,x1,Cbar,delta0 = Params.dEdx_params[mat]

    I = I/1e6  ## convert from eV to MeV

    xp = np.log10(magp/Params.m)
    if xp>=x1:
        delta = 2*np.log(10)*xp - Cbar
    elif xp>=x0:
        delta = 2*np.log(10)*xp - Cbar + a*(x1-xp)**k
    else:
        delta = delta0*10**(2*(xp-x0))

    # mean energy loss in MeV/cm
    dEdx = K*rho*Params.Q**2*Z/A/beta**2*(0.5*np.log(2*me*beta**2*gamma**2*Wmax/I**2) - beta**2 - delta/2)

    dE = dEdx * beta*2.9979e1 * dt

    if dE>(E-Params.m):
        return np.zeros(6)

    newmagp = np.sqrt((E-dE)**2-Params.m**2)
    x[3:] = p*newmagp/magp

    return x

def traverseBField(t, x):
    # x is a 6-element vector (x,y,z,px,py,pz)
    # returns dx/dt
    #
    # if B is in Tesla, dt is in ns, p is in units of MeV/c,  then the basic eq is
    # dp/dt = (89.8755) Qv x B,
    
    
    dxdt = np.zeros(6)

    p = x[3:]
    magp = np.linalg.norm(p)
    E = np.sqrt(magp**2 + Params.m**2)
    v = p/E
    dxdt[:3] = v * 2.9979e-1

    B = Detector.getBField(x[0],x[1],x[2])

    dxdt[3:] = (89.8755) * Params.Q * np.cross(v,B)

    return dxdt


# 4th order runge-kutta integrator
def rk4(x0, update_func, dt, nsteps, cutoff=None, cutoffaxis=None):
    # x0 is a vector of initial values e.g. (x0,y0,z0,px0,py0,pz0)
    # update func is as in dx/dt = update_func(x,t)
    # return value is an N by nsteps+1 array, where N is the size of x0
    # each column is x at the next time step
    #
    # option to cutoff integration once particle reaches certain coordinate
    # along the cutoff axis. 0,1,2 correspond to x,y,z axes. 3 is radial.

    if cutoff!=None and cutoffaxis==None:
        print "Warning: cutoff axis not specified! Not using cutoff"
        cutoff=None

    x0 = np.array(x0)
    x = np.zeros((x0.size, nsteps+1))
    x[:,0] = x0
    t = 0

    # perform the runge-kutta integration
    for i in range(nsteps):
        k1 = update_func(t, x[:,i])
        k2 = update_func(t+dt/2., x[:,i]+dt*k1/2.)
        k3 = update_func(t+dt/2., x[:,i]+dt*k2/2.)
        k4 = update_func(t+dt, x[:,i]+dt*k3)
        dx_Bfield = dt/6. * (k1 + 2*k2 + 2*k3 + k4)
        #dx_Bfield = dt * k1

        # add on the effect of MSC if desired
        dx_MS = np.zeros(x0.size)
        if Params.MSCtype.lower()=='pdg':
            dx_MS = multipleScatterPDG(x[:,i], dt)
        elif Params.MSCtype.lower()=='kuhn':
            dx_MS = multipleScatterKuhn(x[:,i], dt)


        x[:,i+1] = x[:,i] + dx_Bfield + dx_MS

        if Params.EnergyLossOn:
            x[:,i+1] = doEnergyLoss(x[:,i+1], dt)
        if not Params.SuppressStoppedWarning and np.all(x[:,i+1]==0):
            print "Warning: stopped particle! (initial p =", round(np.linalg.norm(x[3:,0])/1000,2), "GeV)"

        # check if particle has stopped
        if np.all(x[:,i+1]==0):
            return x[:,:i+1]
        
        if cutoff!=None:
            if cutoffaxis==3 and x[0,i+1]**2+x[1,i+1]**2>=cutoff**2:
                return x[:,:i+2]
            elif 0<=cutoffaxis<=2 and x[cutoffaxis,i+1]>=cutoff:
                return x[:,:i+2]

        t += dt
        
    if cutoff!=None and x[cutoffaxis,-1]<cutoff:
        print "Warning: cutoff not reached!"

    return x
