#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from parmed import *
from parmed.amber import AmberParm

import mdtraj
import mdtraj.reporters
from parmed.openmm.reporters import RestartReporter 

def run_md(usr_defined_platform, iter, steps, stepsize, vdw_cutoff, temp, output_freq):
    """
    Use openmm for the production phase
    """
    prmtop = AmberPrmtopFile('solvated_perturbed.prmtop')
    restart = AmberInpcrdFile('rst.%02d'%(iter-1))

    # Convert the stepsize from fs to ps
    stepsize = stepsize * 0.001

    # cutoff distance is hand coded here, will change in the future
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=vdw_cutoff*angstroms,
         constraints=HBonds)

    integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, stepsize*picoseconds)

    # Monte Carlo barostat is hard coded here, will change in the future
    system.addForce(MonteCarloBarostat(1*atmospheres, temp*kelvin))

    if usr_defined_platform == 'CUDA':
      platform = Platform.getPlatformByName('CUDA')
      properties = {'CudaPrecision': 'mixed'}
    elif usr_defined_platform == 'OpenCL':
      platform = Platform.getPlatformByName('OpenCL')
      properties = {'OpenCLPrecision': 'mixed'}
    elif usr_defined_platform == 'CPU':
      platform = Platform.getPlatformByName('CPU')
    else:
      platform = Platform.getPlatformByName('Reference')
   
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

    simulation.context.setPositions(restart.positions)
    if restart.boxVectors is not None:
      simulation.context.setPeriodicBoxVectors(*restart.boxVectors)

    simulation.context.setVelocities(restart.velocities) 
    #simulation.context.setVelocitiesToTemperature(temp*kelvin)

    # Simulation.reporters.append(PDBReporter('traj.%02d.pdb'%(iter), 500))
    simulation.reporters.append(StateDataReporter('traj.%02d.out'%(iter), output_freq, step=True,
        potentialEnergy=True, temperature=True, volume=True, density=True))
    simulation.reporters.append(mdtraj.reporters.HDF5Reporter('traj.%02d.h5'%(iter), output_freq))
    restrt = RestartReporter('rst.%02d'%(iter), steps, write_multiple=False,netcdf=False)
 
    simulation.step(steps)            
    state = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    restrt.report(simulation, state)

