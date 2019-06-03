#1! /usr/bin/env python

## Test script to show how to generate particles from a random distribution,
## propagate throught the magnetic field, and simulate their incidence
## on an external Milliqan detector.
##
## For the purposes of visualization, the external detector has been moved
## much closer and enlarged to a 3m x 3m square.
##
## use mode "VIS" to generate a few trajectries and then visualize with matplotlib
## use mode "STATS" to collect statistics on particle hits on the detector, and output
## to a text file
##
## The output format is (Q,m,p,pT,eta,phi,theta,thetaW,thetaV,w,v,pInt), one per column
##  -Q is the charge (program randomly chooses each to be pos/neg)
##  -m is the particle mass
##  -p/pT are the initial (transverse) momentum in MeV
##  -eta/phi are the initial eta/phi of the particle
##  -theta is the of incidence on the detector plane, w.r.t. the normal
##  -w,v are coordinates in the detector plane, where (0,0) is the center.
##  -thetaW, thetaV are the angles of incidence projected on the w,v directions
##  -pInt is the momentum magnitude upon incidence with the detector plane
##
## The above variables are somewhat obscure and were used for testing/debugging.
## To get more useful variables, run the tools/formatOutput.py script.

import math
import time
import os.path
import sys
try:
    sys.path.remove('/home/users/bemarsh/.local/lib/python2.7/site-packages/matplotlib-1.4.3-py2.7-linux-x86_64.egg')
except:
    pass
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import Params
import Integrator
import Detector
import Drawing
from MilliTree import MilliTree
import run_params as rp

suffix = sys.argv[1]

outname = "../output_{0}".format(suffix)

# must run at the beginning of main script. Loads the B field map into memory
Detector.LoadCoarseBField("../bfield/bfield_coarse.pkl")

# turn on CMS magnetic field and PDG multiple scattering
Params.BFieldType = rp.BFieldType
Params.MSCtype = 'pdg'
Params.MatSetup='cms'
# turn on dE/dx energy loss (Bethe-Bloch)
Params.EnergyLossOn = True
# charge and mass of the particle. Q in units of electric charge, m in MeV
Params.Q = rp.particleQ
Params.m = rp.particleM
#suppress annoying warnings
Params.SuppressStoppedWarning = False
Params.RockBegins = rp.RockBegins

if rp.useCustomMaterialFunction:
    Params.matFunction = rp.matFunction

# make sure numbers are new each run
ROOT.gRandom.SetSeed(0)

#rootfile = ROOT.TFile(rp.pt_spect_filename)
# this is a 1D pT distribution (taken from small-eta events)
#pt_dist = rootfile.Get("pt")

dt = rp.dt
nsteps = rp.max_nsteps    

center = rp.centerOfDetector
distToDetect = np.linalg.norm(center)
normToDetect = center/distToDetect

detV = np.array([0., 1., 0.])
detW = np.cross(normToDetect, detV)

detWidth = rp.detWidth
detHeight = rp.detHeight
detDepth = rp.detDepth

detectorDict = {"norm":normToDetect, "dist":distToDetect, "v":detV, 
            "w":detW, "width":detWidth, "height":detHeight, "depth":detDepth}

# the four corners (only for drawing)
c1,c2,c3,c4 = Detector.getDetectorCorners(detectorDict)

intersects = []
ntotaltrajs = 0

mt = MilliTree()

#out_rootfile = ROOT.TFile(outname+'.root','recreate')
#out_tree = ROOT.TTree('Events','tree with particles propagated to face of detector')
#out_mother_id = array('i', [ 0 ] )
#out_id = array('i', [ 0 ] )
#out_pt = array('f', [ 0. ] )
#out_eta = array('f', [ 0. ] )
#out_phi = array('f', [ 0. ] )
#out_pthat = array('f', [ 0. ] )
#out_weight = array('f', [ 0. ] )
#out_charge = array('f', [ 0. ] )
#out_event_num = array('f', [ 0. ] )
#out_opt = array('f', [ 0. ] )
#out_oeta = array('f', [ 0. ] )
#out_ophi = array('f', [ 0. ] )
#out_tree.Branch('pt', d, 'out_pt/F')
#out_tree.Branch('eta', d, 'out_eta/F')
#out_tree.Branch('phi', d, 'out_phi/F')

starttime = time.time()

rootfile = ROOT.TFile(rp.infile)
tree = rootfile.Get('EventTree')
nevents = tree.GetEntries()
for event in tree:
    
    # convert to cartesian momentum
    p = event.pT * np.array([np.cos(event.phi),np.sin(event.phi),np.sinh(event.eta)])
    x0 = np.array([0,0,0,p[0],p[1],p[2]])
    
    # simulate until nsteps steps is reached, or the particle passes x=10
    traj,tvec = Integrator.rk4(x0, dt, nsteps, cutoff=rp.cutoff, cutoffaxis=3, use_var_dt=rp.use_var_dt)

    # compute the intersection. Will return None if no intersection
    findIntersection = Detector.FindIntersection
    if rp.useCustomIntersectionFunction:
        findIntersection = rp.intersectFunction
    intersection, t, theta, thW, thV, pInt = findIntersection(traj, tvec, detectorDict)
    if intersection is not None:
        intersects.append(intersection)
        print len(trajs), ": p =",magp, ", eta =", eta, ", phi =", phi, ", eff =", float(len(intersects))/ntotaltrajs
        w = np.dot(intersection, detectorDict['w'])
        v = np.dot(intersection, detectorDict['v'])
        magpint = np.linalg.norm(pInt)
        mt.SetValues(intersection, pInt)
        mt.Fill()
    else:
        intersection = np.array([-99,-99,-99])
        pInt = np.array([-99,-99,-99])
        mt.SetValues(intersection, pInt)
        mt.Fill()

endtime = time.time()

print "Efficiency:", float(len(intersects))/nevents
print "Total time: {0:.2f} sec".format(endtime-starttime)
print "Time/Hit: {0:.2f} sec".format((endtime-starttime)/nevents)

mt.Write("../output_{0}.root".format(suffix))

#fid = ROOT.TFile("../output_{0}.root".format(suffix), "UPDATE")
#
#hhits = ROOT.TH1F("hhits","",1,0,2)
#hsims = ROOT.TH1F("hsims","",1,0,2)
#hhits.Fill(1, ntrajs)
#hsims.Fill(1, ntotaltrajs)
#hhits.Write()
#hsims.Write()
#
#fid.Close()

