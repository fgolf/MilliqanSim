import numpy as np

mode = "VIS"
ntrajs = 1
pt_spect_filename = "../p_eta_dist/combined_PtSpect_Eta0p16.root"
dt = 0.1   #timestep in ns
max_nsteps = 5000
cutoff = 34.
use_var_dt = False
BFieldType = "cms"

#particleQ = 1.0  # in electron charge units
#particleM = 105. # in MEV
particleQ = 0.01  # in electron charge units
particleM = 50. # in MEV

distToDetector = 33.
eta = 0.16
## don't touch
theta = np.pi/2 - 2*np.arctan(np.exp(-eta))
x = distToDetector*np.cos(theta)
z = distToDetector*np.sin(theta)
centerOfDetector = np.array([x, 0., z])
##
RockBegins = distToDetector - 17.

detWidth = 0.5
detHeight = 0.5
detDepth = 1.0

etabounds = (eta-0.08, eta+0.08)
ptCut = 17.
phibounds = (0.00, 0.22)

useCustomMaterialFunction = False
useCustomIntersectionFunction = False
useCustomOutput = False

infile = "/Users/fgolf/Dropbox/milliQan/sim/Milliqan-Pythia-Simulation/out.root"
