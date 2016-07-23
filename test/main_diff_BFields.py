#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import Integrator
import Detector
import Params
import Drawing

#Detector.LoadCoarseBField("bfield/bfield_coarse.pkl")
#print "Loaded coarse field"
#Detector.LoadFineBField("bfield/bfield_x.pkl","bfield/bfield_y.pkl","bfield/bfield_z.pkl")
#print "Loaded fine field"

Params.BFieldType = 'cms'
Params.Q = 1.0
Params.MSCtype = 'none'
Params.EnergyLossOn=False
Params.Interpolate=True

## y vs x, eta 0
def yvsx():
    dt = 0.05
    nsteps = 800
    
    p0 = [20000, 0, 0]
    x0 = np.array([0,0,0]+p0)
    
    detectorDict = { "norm": np.array([1,0,0]),
                     "dist": 0,
                     "v": np.array([0,1,0]),
                     "w": np.array([0,0,1]),
                     "width": 20,
                     "height": 20  }
    
    Params.UseFineBField = False
    traj_coarse = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    Params.UseFineBField = True
    traj_fine = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    
    yfine = []
    ycoarse = []
    zfine = []
    zcoarse = []
    rvals = np.linspace(0,8.5,100)
    for r in rvals:
        detectorDict["dist"] = r
        
        intersect_coarse = Detector.FindIntersection(traj_coarse,detectorDict)
        ycoarse.append(100*intersect_coarse[0][1])
        zcoarse.append(100*intersect_coarse[0][2])
        intersect_fine = Detector.FindIntersection(traj_fine,detectorDict)
        yfine.append(100*intersect_fine[0][1])
        zfine.append(100*intersect_fine[0][2])
        
    yfine = np.array(yfine)
    ycoarse = np.array(ycoarse)
    zfine = np.array(zfine)
    zcoarse = np.array(zcoarse)
        
    plt.figure(1)

    plt.subplot(1,1,1)
    ax1 = plt.gca()
    fplot, = ax1.plot(rvals, yfine, '-b', linewidth=2, label='1 cm field')
    cplot, = ax1.plot(rvals, ycoarse, '--r', linewidth=2, label='10 cm field')
    ax2 = ax1.twinx()
    dplot, = ax2.plot(rvals,abs(yfine-ycoarse), '-', color='0.75', label='Difference')
    ax1.set_ylabel('y (cm)')
    ax2.set_ylabel('diff (cm)')
    
    # plt.subplot(2,1,2)
    # ax1 = plt.gca()
    # ax1.plot(rvals, zfine, '-b')
    # ax1.plot(rvals, zcoarse, '--r')
    # ax2 = ax1.twinx()
    # ax2.plot(rvals,abs(zfine-zcoarse), '-g')
    # ax1.set_ylabel('y (cm)')
    # ax2.set_ylabel('diff (cm)')
    
    ax1.set_xlabel('x (m)')
    ax1.set_title(r'20 GeV muon trajectory ($\eta=0$)')
    
    plt.legend(handles=[fplot,cplot,dplot],loc='lower left')
    
    plt.savefig("/home/users/bemarsh/public_html/backup/B-integrator/fine_vs_coarse/r_vs_x_eta0_interp.png")
    
    plt.show()

## deltaTheta vs theta
def dthvsth():
    thvals = np.linspace(0,np.pi,200)
    Bdlvvals_fine = []
    Bdlvvals_coarse = []
    Bdlyvals_fine = []
    Bdlyvals_coarse = []
    
    for th in thvals:
        rhat = np.array([np.sin(th),0,np.cos(th)])
        yhat = np.array([0,1,0])
        vhat = np.cross(rhat,yhat)
        dr = 1.0
        ivvals_fine = []
        ivvals_coarse = []
        iyvals_fine = []
        iyvals_coarse = []
        rvals = np.arange(0,1500,dr)
        for r in rvals:
            rvec = r*rhat/100.0
    
            Params.UseFineBField = False
            B = Detector.getBField(rvec[0],rvec[1],rvec[2])
            Bv = np.dot(vhat,B)
            ivvals_coarse.append(Bv)
            iyvals_coarse.append(B[1])
    
            Params.UseFineBField = True
            B = Detector.getBField(rvec[0],rvec[1],rvec[2])
            Bv = np.dot(vhat,B)
            ivvals_fine.append(Bv)
            iyvals_fine.append(B[1])
    
        intv_coarse = 2*np.sum(ivvals_coarse)-ivvals_coarse[0]-ivvals_coarse[-1]
        Bdlvvals_coarse.append(intv_coarse * dr/100/2)
        intv_fine = 2*np.sum(ivvals_fine)-ivvals_fine[0]-ivvals_fine[-1]
        Bdlvvals_fine.append(intv_fine * dr/100/2)
        inty_coarse = 2*np.sum(iyvals_coarse)-iyvals_coarse[0]-iyvals_coarse[-1]
        Bdlyvals_coarse.append(inty_coarse * dr/100/2)
        inty_fine = 2*np.sum(iyvals_fine)-iyvals_fine[0]-iyvals_fine[-1]
        Bdlyvals_fine.append(inty_fine * dr/100/2)
    
    plt.figure(2)
    plt.plot(thvals*180/np.pi,Bdlvvals_fine, '-b', linewidth=2, label='1 cm field')
    plt.plot(thvals*180/np.pi,Bdlvvals_coarse, '--r', linewidth=2, label='10 cm field')
    plt.plot(thvals*180/np.pi,Bdlyvals_fine, '-g', linewidth=2)
    plt.plot(thvals*180/np.pi,Bdlyvals_coarse, '--y', linewidth=2)
    plt.xlabel(r'$\theta$ (deg)')
    plt.title(r'$\int$ $B dl$ ($\phi=0$ plane)')
    plt.legend(fontsize='small')
    plt.text(60,0.35,"y direction")
    plt.text(60,9.2,"xz direction")
    
    plt.savefig('/home/users/bemarsh/public_html/backup/B-integrator/fine_vs_coarse/intBdl_allTheta_interp.png')

    plt.show()


## y at 15 m away
def yat15():

    dt = 0.05
    nsteps = 1200
    
    detectorDict = { "norm": np.array([1,0,0]),
                     "dist": 15,
                     "v": np.array([0,1,0]),
                     "w": np.array([0,0,1]),
                     "width": 20,
                     "height": 20  }
    
    thvals = np.linspace(0,np.pi,200)
    yfine = []
    ycoarse = []
    
    for th in thvals:
        detectorDict["norm"] = np.array([np.sin(th), 0, np.cos(th)])
        detectorDict["w"] = np.cross(detectorDict["norm"],detectorDict["v"])
        p0 = 20000 * detectorDict["norm"]
        x0 = np.append([0,0,0],p0)
    
        Params.UseFineBField = False
        traj_coarse = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
        Params.UseFineBField = True
        traj_fine = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    
        intersect_coarse = Detector.FindIntersection(traj_coarse,detectorDict)
        ycoarse.append(100*intersect_coarse[0][1])
        intersect_fine = Detector.FindIntersection(traj_fine,detectorDict)
        yfine.append(100*intersect_fine[0][1])
    
    yfine = np.array(yfine)
    ycoarse = np.array(ycoarse)
    thvals *= 180/np.pi
    
    plt.figure(3)
    
    plt.subplot(1,1,1)
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax2.plot(thvals,abs(yfine-ycoarse), '-', color='0.75')
    ax1.plot(thvals, yfine, '-b', linewidth=2)
    ax1.plot(thvals, ycoarse, '--r', linewidth=2)
    ax1.set_ylabel('y (cm)')
    ax2.set_ylabel('diff (cm)')
    
    ax1.set_xlabel(r'$\theta$ (deg)')
    ax1.set_title(r'y at r=15m, 20 GeV muon')
    
    plt.savefig("/home/users/bemarsh/public_html/backup/B-integrator/fine_vs_coarse/deltaR_allTheta_nointerp.png")
    

    plt.show()


def compareTimeSteps(dt1, dt2):
    
    Params.UseFineBField=True

    total_time = 35
    nsteps1 = int(round(total_time/dt1))
    nsteps2 = int(round(total_time/dt2))

    p0 = [20000, 0, 0]
    x0 = np.array([0,0,0]+p0)
    
    detectorDict = { "norm": np.array([1,0,0]),
                     "dist": 0,
                     "v": np.array([0,1,0]),
                     "w": np.array([0,0,1]),
                     "width": 20,
                     "height": 20  }
    
    traj1 = Integrator.rk4(x0, Integrator.traverseBField, dt1, nsteps1)
    traj2 = Integrator.rk4(x0, Integrator.traverseBField, dt2, nsteps2)
    
    y1 = []
    y2 = []
    rvals = np.linspace(0,8.5,100)
    for r in rvals:
        detectorDict["dist"] = r
        
        intersect1 = Detector.FindIntersection(traj1,detectorDict)
        y1.append(100*intersect1[0][1])
        intersect2 = Detector.FindIntersection(traj2,detectorDict)
        y2.append(100*intersect2[0][1])
        
    y1 = np.array(y1)
    y2 = np.array(y2)
        
    plt.figure(1)

    plt.subplot(1,1,1)
    ax1 = plt.gca()
    fplot, = ax1.plot(rvals, y1, '-b', linewidth=2, label='dt='+str(dt1)+' ns')
    cplot, = ax1.plot(rvals, y2, '--r', linewidth=2, label='dt='+str(dt2)+' ns')
    ax2 = ax1.twinx()
    dplot, = ax2.plot(rvals,abs(y1-y2), '-', color='0.75', label='Difference')
    ax1.set_ylabel('y (cm)')
    ax2.set_ylabel('diff (cm)')
        
    ax1.set_xlabel('x (m)')
    ax1.set_title(r'20 GeV muon trajectory ($\eta=0$)')
    
    plt.legend(handles=[fplot,cplot,dplot],loc='lower left')
    
    plt.savefig("/home/users/bemarsh/public_html/backup/B-integrator/fine_vs_coarse/comp_dt_eta0_interp.png")
    
    plt.show()
    



