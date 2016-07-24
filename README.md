# MilliqanSim

## What is it?

A general-purpose particle trajectory simulator that can propagate a particle through a magnetic field, implement
basic multiple scattering and Bethe-Bloch energy loss, and generate statistics on the intersections of these trajectories
with an exernal plane (e.g., the Milliqan detector).

## Dependencies

The package is entirely written in Python (v 2.7). The only absolute requirement is the [numpy](http://www.numpy.org) 
library. Optional but useful are the [matplotlib](http://matplotlib.org/) module to make plots, and pyROOT to get input from or output to a ROOT file.

## Code structure

The `src` directory contains 5 main modules, which are imported into a user-written main script.

1. **Params.py** - contains all of the configurable parameters, material definitions, etc.
2. **Integrator.py** -  routines to propagate the particle through a magnetic field
3. **Detector.py** - Functions related to detector and environment setup - magnetic field, material, etc.
4. **MultipleScatter.py** - methods implementing multiple scattering and Bethe-Bloch energy loss
5. **Drawing.py** - various methods to make plots of trajectories

The `test` directory contains a few well-commented example scripts that should give a good idea of how to use the program. In particular, `milliqan_test.py` is a program to simulate particle hits on a pseudo Milliqan detector.
