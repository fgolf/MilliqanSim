#! /usr/bin/env python

import math
import os.path
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import Params
import Integrator
import Detector
import Drawing

# VIS or STATS
mode = "STATS"
visWithStats = True

if mode=="VIS":
    ntrajs = 15
    trajs = []
if mode=="STATS":
    ntrajs = 5
    trajs = []

outname = "data/test"

Params.BFieldType = 'cms'
Params.MSCtype = 'pdg'
Params.EnergyLossOn = True
Params.UseFineBField = False
Params.Q = 1.0
Params.m = 105.
Params.SuppressStoppedWarning = True

Detector.LoadCoarseBField("../B-integrator/bfield/bfield_coarse.pkl")

# make sure numbers are new each run
ROOT.gRandom.SetSeed(0)

rootfile = ROOT.TFile("p_eta_dists/DY_CutPtSpect_v2.root")
pt_dist = rootfile.Get("pt")

for i in range(1,pt_dist.GetNbinsX()+1):
        if pt_dist.GetBinContent(i) < 0:
            pt_dist.SetBinContent(i,0)

dt = 0.2
nsteps = 700

# set up detector plane

# normal to the plane (x,0,z)
normToDetect = np.array([9.0,0,0])
# distance from origin to plane
distToDetect = np.linalg.norm(normToDetect)
normToDetect = normToDetect/np.linalg.norm(normToDetect)

#center point
center = normToDetect*distToDetect
# y-axis (in plane)
detV = np.array([0,1,0])
# another orthogonal vector to norm in plane ((normToDetect, detV, detW) form an ON basis)
detW = np.cross(normToDetect,detV)

detWidth = 3.0
detHeight = 3.0

detectorDict = {"norm":normToDetect, "dist":distToDetect, "v":detV, 
                "w":detW, "width":detWidth, "height":detHeight}

# the four corners (for drawing)
c1 = center + detW*detWidth/2 + detV*detHeight/2
c2 = center + detW*detWidth/2 - detV*detHeight/2
c3 = center - detW*detWidth/2 - detV*detHeight/2
c4 = center - detW*detWidth/2 + detV*detHeight/2

intersects = []

stats = ROOT.TNtuple("stats","stats","q:p:pT:eta:phi:theta:thW:thV:w:v:pInt")
ntotaltrajs = 0

if mode=="STATS":
    if os.path.isfile(outname+".txt"):
        ow = 'q'
        while ow not in 'yYnN':
            ow = raw_input("Overwrite file? (y/n) ")
        if ow in 'yY':
            txtfile = open(outname+".txt",'w')
        else:
            txtfile = open(outname+".txt",'a')
    else:
        txtfile = open(outname+".txt",'w')
    txtfile.close()

while len(trajs)<ntrajs:
    magp = ROOT.Double(-1)
    eta = ROOT.Double(-1)

    etalow =  -0.1
    etahigh =  0.1

    while magp<10:
        magp = pt_dist.GetRandom()

    eta = np.random.rand()*(etahigh-etalow) + etalow

    th = 2*np.arctan(np.exp(-eta))
    magp = magp/np.sin(th)
    phimin, phimax =  -0.2,0.2
    phi = np.random.rand() * (phimax-phimin) + phimin
    Params.Q *= np.random.randint(2)*2 - 1
    phi *= Params.Q/abs(Params.Q)

    p = 1000*magp * np.array([np.sin(th)*np.cos(phi),np.sin(th)*np.sin(phi),np.cos(th)])
    x0 = np.array([0,0,0,p[0],p[1],p[2]])

    traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    ntotaltrajs += 1
    if mode=="VIS":
        trajs.append(traj)

    intersection, theta, thW, thV, pInt = Detector.FindIntersection(traj, detectorDict)
    if intersection is not None:
        intersects.append(intersection)
        print len(trajs), ": p =",magp, ", eta =", eta, ", phi =", phi, ", eff =", float(len(intersects))/ntotaltrajs
        if mode=="VIS":
            pass
        if mode=="STATS":
            if visWithStats:
                trajs.append(traj)
            else:
                trajs.append(0)
            w = np.dot(intersection, detW)
            v = np.dot(intersection, detV)
            stats.Fill(Params.Q,magp,magp*np.sin(th),eta,phi,theta, thW, thV, w, v, pInt)
            txtfile = open(outname+".txt",'a')
            txtfile.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:f}\t{9:f}\t{10:f}\n".format(Params.Q,magp,magp*np.sin(th),eta,phi,theta,thW,thV,w,v,pInt))

if mode=="STATS":
    fid = ROOT.TFile(outname+".root","RECREATE")
    stats.Write()
    fid.Close()

print "Efficiency:", float(len(intersects))/ntotaltrajs

if mode=="VIS" or visWithStats:
    plt.figure(num=1, figsize=(15,7))

    Drawing.Draw3Dtrajs(trajs, subplot=121)

    Drawing.DrawLine(c1,c2,is3d=True)
    Drawing.DrawLine(c2,c3,is3d=True)
    Drawing.DrawLine(c3,c4,is3d=True)
    Drawing.DrawLine(c4,c1,is3d=True)

    for i in range(len(intersects)):
        Drawing.DrawLine(intersects[i],intersects[i],is3d=True,linestyle='None',marker='o',color='r')

    Drawing.DrawXYslice(trajs, subplot=122)

    plt.figure(num=2)
    Drawing.DrawXZslice(trajs)

    plt.show()
