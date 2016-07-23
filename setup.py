#! /usr/bin/env python

import os
import sys
import pipes

print "Checking for required modules"
try:
    import numpy
    print "Successfully imported numpy"
except:
    print "ERROR: could not find the numpy module!"
    sys.exit(1)

try:
    import matplotlib.pyplot
    print "Successfully imported matplotlib.pyplot"
except:
    print "WARNING: coule not find matplotlib. Simulation will work, but you will not be able to plot or use the Drawing module."
