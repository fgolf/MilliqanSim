## Detector.py
## methods relating to detector and environment properties

import cPickle as pickle
import numpy as np
import Params

def LoadCoarseBField(fname):
    
    if Params.BFieldUsePickle:
        Params.Bx,Params.By,Params.Bz,Params.Bmag = pickle.load(open(fname,"rb"))
        Params.BFieldLoaded = True
        return

    NR = Params.RMAX/Params.DR + 1
    NZ = (Params.ZMAX-Params.ZMIN)/Params.DZ + 1
    NPHI = (Params.PHIMAX-Params.PHIMIN)/Params.DPHI + 1

    Params.Bx   = np.zeros((NR,NZ,NPHI))
    Params.By   = np.zeros((NR,NZ,NPHI))
    Params.Bz   = np.zeros((NR,NZ,NPHI))
    Params.Bmag = np.zeros((NR,NZ,NPHI))
    
    with open(fname,'r') as fid:
        for line in fid:
            sp = line.strip().split()
            if len(sp)!=4:
                continue

            r = float(sp[0])
            z = float(sp[1])
            phi = float(sp[2])
            B = tuple([float (a) for a in sp[3].strip("()").split(",")])

            iz = int((z-Params.ZMIN)/Params.DZ)
            ir = int(r/Params.DR)
            iphi = int((phi-Params.PHIMIN)/Params.DPHI)

            Params.Bmag[ir,iz,iphi] = np.linalg.norm(B)
            Params.Bx[ir,iz,iphi] = B[0]
            Params.By[ir,iz,iphi] = B[1]
            Params.Bz[ir,iz,iphi] = B[2]

    Params.BFieldLoaded = True


def LoadFineBField(fnamex, fnamey=None, fnamez=None):
    
    if fnamey==None:
        Params.Bx,Params.By,Params.Bz = pickle.load(open(fname,"rb"))
        Params.BFieldLoaded = True
        return

    Params.Bxf = pickle.load(open(fnamex,"rb"))
    Params.Byf = pickle.load(open(fnamey,"rb"))
    Params.Bzf = pickle.load(open(fnamez,"rb"))

    return


def getMaterial(x,y,z):

    if Params.MatSetup == 'iron':
        return 'fe'

    if Params.MatSetup == 'sife':
        if x<4:
            return 'si'
        else:
            return 'fe'

    withinLength = -Params.solLength/2 < z < Params.solLength/2
    r = np.sqrt(x**2+y**2)

    if not withinLength:
        return 'air'
    
    if r < 1.29:
        mat = 'air'
    elif r < 1.8:
        mat = 'pbwo4'
    elif r < 2.95:
        mat = 'fe'
    elif r < 4.0:
        mat = 'fe'
    elif r < 7.0:
        mat = 'fe'
    else:
        mat = 'air'

    if x>16.5:
        mat = 'concrete'

    return mat
    

def getBField(x,y,z):

    if Params.BFieldType.lower() == 'none':
        return np.zeros(3)

    if Params.BFieldType.lower() == 'uniform':
        return np.array([0.,0.,1.])

    if Params.BFieldType.lower() == 'updown':
        r = np.sqrt(x**2+y**2)
        if r<4:
            return np.array([0,0,3])
        else:
            return np.array([0,0,-3])/r

    ## correct for cm usage in bfield file
    x *= 100
    y *= 100
    z *= 100
    r = np.sqrt(x**2+y**2)

    if Params.UseFineBField:
        ZMIN,ZMAX,DZ,RMIN,RMAX,DR,PHIMIN,PHIMAX,DPHI = \
            Params.ZMINf, Params.ZMAXf, Params.DZf, \
            Params.RMINf, Params.RMAXf, Params.DRf, \
            Params.PHIMINf, Params.PHIMAXf, Params.DPHIf
    else:
        ZMIN,ZMAX,DZ,RMIN,RMAX,DR,PHIMIN,PHIMAX,DPHI = \
            Params.ZMIN, Params.ZMAX, Params.DZ, \
            Params.RMIN, Params.RMAX, Params.DR, \
            Params.PHIMIN, Params.PHIMAX, Params.DPHI

    if z>ZMIN and z<ZMAX and r<RMAX:

        r = np.sqrt(x**2+y**2)
        phi = np.arctan2(y,x) * 180/np.pi

        if phi<0:
            phi += 360
        
        # nearR = int(DR*round(r/DR))
        nearZ = int(DZ*round(z/DZ))
        nearPHI = int(DPHI*round(phi/DPHI))
        
        if nearPHI==360:
            nearPHI = 0

        iz = (nearZ-ZMIN)/DZ
        iphi = (nearPHI-PHIMIN)/DPHI

        ir = r/DR
        irlow = int(np.floor(ir))
        irhigh = int(np.ceil(ir))
        irfrac = ir-irlow

        if not Params.Interpolate:
            irfrac = 1
            irhigh = int(np.round(ir))
        
        if Params.UseFineBField:
            Bx = irfrac*Params.Bxf[irhigh,iz,iphi]+(1-irfrac)*Params.Bxf[irlow,iz,iphi]
            By = irfrac*Params.Byf[irhigh,iz,iphi]+(1-irfrac)*Params.Byf[irlow,iz,iphi]
            Bz = irfrac*Params.Bzf[irhigh,iz,iphi]+(1-irfrac)*Params.Bzf[irlow,iz,iphi]
        else:
            Bx = irfrac*Params.Bx[irhigh,iz,iphi]+(1-irfrac)*Params.Bx[irlow,iz,iphi]
            By = irfrac*Params.By[irhigh,iz,iphi]+(1-irfrac)*Params.By[irlow,iz,iphi]
            Bz = irfrac*Params.Bz[irhigh,iz,iphi]+(1-irfrac)*Params.Bz[irlow,iz,iphi]
            
        return np.array([Bx,By,Bz])

    else:
        return np.zeros(3)

def FindIntersection(traj, detectorDict):
    # find the intersection with a plane with normal norm
    # and distance to origin dist. returns None if no intersection

    norm = detectorDict["norm"]
    dist = detectorDict["dist"]

    for i in range(traj.shape[1]-1):
        p1 = traj[:3,i]
        p2 = traj[:3,i+1]
        
        proj2 = np.dot(p2,norm)

        if proj2>=dist:
            proj1 = np.dot(p1,norm)
            intersect = p1+(dist-proj1)/(proj2-proj1)*(p2-p1)

            vHat = detectorDict["v"]
            wHat = detectorDict["w"]
            center = norm*dist

            w = np.dot(intersect,wHat)
            v = np.dot(intersect,vHat)
            
            if abs(w) < detectorDict["width"]/2 and \
               abs(v) < detectorDict["height"]/2:
                unit = (p2-p1)/np.linalg.norm(p2-p1)
                theta = np.arccos(np.dot(unit,norm))
                
                projW = np.dot(unit,wHat)
                projV = np.dot(unit,vHat)

                thW = np.arcsin(projW/np.linalg.norm(unit-projV*vHat))
                thV = np.arcsin(projV/np.linalg.norm(unit-projW*wHat))

                # momentum when it hits detector
                pInt = np.linalg.norm(traj[3:,i])

                return intersect,theta,thW,thV,pInt

            break

    return None, None, None, None, None












