# Params.py
# control all global variables

import Detector
import numpy as np

## NOTE: all numbers for the below two lists are found at http://pdg.lbl.gov/2015/AtomicNuclearProperties/

## materials definition. (atomic num, atomic weight, density (g/cm3), radiation length (m))
materials = { "fe"   : (26, 55.845, 7.874, .01757), 
              "si"   : (14, 28.0855, 2.329, .0937),
              "pb"   : (82, 207.2, 11.35, .005612),
              "air"  : (7.34, 14.719, 1.205e-3, 3.04e2),  
              "pbwo4": (31.3, 75.8, 8.3, 0.008903),
              "bc408": (3.37, 6.23, 1.032, 0.4254),
              "concrete": (11.10, 22.08, 2.3, .1155),
              "rock": (11.0, 22.0, 2.65, .1002)  }

## parameters for dEdx (I, a, k, x0, x1, Cbar, delta0)
dEdx_params = { "fe"   : (286.0, 0.14680, 2.9632, -0.0012, 3.1531, 4.2911, 0.12),
                "si"   : (173.0, 0.14921, 3.2546, 0.2015, 2.8716, 4.4355, 0.14),
                "pb"   : (823.0, 0.09359, 3.1608, 0.3776, 3.8073, 6.2018, 0.14),
                "air"  : (85.7, 0.10914, 3.3994, 1.7418, 4.2759, 10.5961, 0.0),
                "pbwo4": (600.7, 0.22758, 3.0, 0.4068, 3.0023, 5.8528, 0.0),
                "bc408": (64.7, 0.16101, 3.2393, 0.1464, 2.4855, 3.1997, 0.0),
                "concrete": (135.2, .07515, 3.5467, .1301, 3.0466, 3.9464, 0.0),
                "rock": (136.4, .08301, 3.4210, .0492, 3.0549, 3.7738, 0.0) }

## these may be updated in main program
Q = 1  ## in units of e
m = 105.658  ## in MeV 
solRad = 3.6  ## in m
solLength = 15.0   ## in m
MSCtype = 'PDG'
EnergyLossOn = False
SuppressStoppedWarning = True
BFieldType = 'CMS'
BFieldUsePickle = True
UseFineBField = False
MatSetup = 'cms'
RockBegins = 999999.  #past this distance from IP, solid concrete
Interpolate = True
matFunction = Detector.getMaterial

## internal parameters. don't touch
MSCWarning = False
BFieldLoaded = False

## parameters to load bfield
ZMIN = -1500
ZMAX = 1500
DZ = 10
RMIN = 0
RMAX = 900
DR = 10
PHIMIN = 0
PHIMAX = 355
DPHI = 5
Bx = np.array([])
By = np.array([])
Bz = np.array([])
Bmag = np.array([])

## parameters to load fine bfield
ZMINf = -1500
ZMAXf = 1500
DZf = 1
RMINf = 0
RMAXf = 900
DRf = 1
PHIMINf = 0
PHIMAXf = 355
DPHIf = 5
Bxf = np.array([])
Byf = np.array([])
Bzf = np.array([])
Bmagf = np.array([])
