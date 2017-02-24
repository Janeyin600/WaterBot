#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys


def write_min_in():
    """
     Write an Amber input file for minimization
    """
    if os.path.isfile('min.in'):
       sp.call('rm min.in', shell = True)
    
    min_file = open('min.in', 'a')
    min_in = '''Minimizing.
    &cntrl
     imin = 1,
     ntx = 1, 
     ntpr = 50,
     maxcyc = 50000,
     ncyc = 5000,
     ntxo = 1,
     irest = 0,
     ntf = 1,
     ntc = 1,
     ntb = 1,
     cut = 9.0,
    /\n'''
    min_file.write(min_in)
    min_file.close()

def write_therm1_in():
    """
     Write an Amber input file for thermalizing
    """
    if os.path.isfile('therm1.in'):
       sp.call('rm therm1.in', shell = True)

    therm1_file = open('therm1.in', 'a')
    therm1_in = '''Thermalizing, NVT.
    &cntrl
     imin = 0,
     ntx = 1,
     irest = 0,
     ntpr = 500,
     ntwr = 500,
     ntxo = 1,
     ntwx = 0,
     ioutfm = 1,
     iwrap = 1,
     ntf = 2,
     ntc = 2,
     cut = 9.0,
     ntt = 3,
     tempi = 10.0,
     temp0 = 10.0,
     gamma_ln = 1.0,
     ig = -1,
     ntp = 0,
     nstlim = 500,
     dt = 0.002,
    /\n'''
    therm1_file.write(therm1_in)
    therm1_file.close()


def write_therm2_in(temperature):
    """
     Write an Amber input file for the heating process
    """
    if os.path.isfile('therm2.in'):
       sp.call('rm therm2.in', shell = True)
    
    therm2_file = open('therm2.in', 'a')
    therm2_in = '''Thermalizing, NVT.
    &cntrl
     imin = 0,
     ntx = 5,
     irest = 1,
     ntpr = 500,
     ntwr = 50000,
     ntwx = 0,
     ntxo = 1,
     ioutfm = 1,
     iwrap = 1,
     ntf = 2,
     ntc = 2,
     cut = 9.0,
     ntt = 3,
     gamma_ln = 1.0,
     ig = -1,
     ntp = 0,
     nstlim = 50000,
     dt = 0.002,
    /\n''' 
           
    therm2_file.write(therm2_in)
    therm2_file.write('\n   &wt type=\'TEMP0\', istep1=0, istep2=25000, value1=10.0, value2=%5.2f,/'%(temperature))
    therm2_file.write('\n   &wt type = \'END\', /')
    therm2_file.close()

def write_eqnpt_in(temperature):
    """
     Write an Amber input file for the npt equilibration
    """
    if os.path.isfile('eqnpt.in'):
       sp.call('rm eqnpt.in', shell = True)

    eqnpt_file = open('eqnpt.in', 'a')
    eqnpt_in = '''Equilibrate, NPT.
    &cntrl
     imin = 0,
     ntx = 5,
     irest = 1,
     ntpr = 500,
     ntwr = 2500,
     ntwx = 0,
     ntxo = 1,
     ioutfm = 1,
     iwrap = 1,
     ntf = 2,
     ntc = 2,
     cut = 9.0,
     ntt = 3,
     gamma_ln = 1.0,
     ig = -1,
     ntp = 1,
     barostat = 2,
     nstlim = 2500,
     dt = 0.002,
     temp0 = %5.2f
    /\n'''%(temperature)
    eqnpt_file.write(eqnpt_in)
    eqnpt_file.close()


