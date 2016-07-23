#! /usr/bin/python

## compare scatter angle distributions for 2 different MSC algorithms

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
import Detector
import MultipleScatter
import Params

# modify globals from their default values
Params.BFieldOn = False
Params.MatSetup = "sife"

dt = 0.2

p0 = [10000,0,0]
x0 = np.array([0,0,0]+p0)

thValsKuhn = []
thValsPDG = []

for i in range(200000):
    th = MultipleScatter.getScatterAngleKuhn(x0, dt) * 180/np.pi
    thValsKuhn.append(th)
    # sys.stdout.write("\b\b\b" + "{0:3d}%".format(int(i/200000.*100)))

Params.MSCtype = 'PDG'
for i in range(200000):
    thx,thy,yx,yy = MultipleScatter.getScatterAnglePDG(x0, dt)
    th = np.sqrt(thx**2+thy**2) * 180/np.pi
    thValsPDG.append(th)

plt.figure(1)

#w = 1./200000*np.ones(200000)
plt.hist(thValsKuhn, range=(0,.5), bins=90, log=False, alpha=0.4, histtype='stepfilled', label='Kuhn')
plt.hist(thValsPDG,  range=(0,.5), bins=90, log=False, alpha=0.4, histtype='stepfilled', label='PDG')
#plt.axis([0,.5,5e-1,1e5])
plt.xlabel('Angle (deg)')
plt.title('Scattering Angle through 6 cm of Si')
plt.legend()

print "Finished dtheta distribution"

plt.show()

