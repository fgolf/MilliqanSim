import numpy as np
import Custom

mode = "STATS"
ntrajs = 10
pt_spect_filename = "../p_eta_dist/combined_PtSpect_Eta0p16.root"
dt = 0.05   #timestep in ns
max_nsteps = 10000
cutoff = 37  # stop simulation when dist(IP) > cutoff

particleQ = 1.0  # in electron charge units
particleM = 105. # in MEV
RockBegins = 9999999.

distToDetector = 33.
eta = 0.16
## don't touch
theta = np.pi/2 - 2*np.arctan(np.exp(-eta))
x = distToDetector*np.cos(theta)
z = distToDetector*np.sin(theta)
centerOfDetector = np.array([x, 0., z])
##

detWidth = 0.11
detHeight = 0.11

etabounds = (eta-0.05, eta+0.05)
ptCut = 20.
phibounds = (0.00, 0.22)

useCustomMaterialFunction = True
matFunction = Custom.getMaterial

useCustomOutput = True
outputFunction = Custom.trajOutput
