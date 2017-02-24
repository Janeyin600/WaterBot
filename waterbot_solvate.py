#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys

def setup_solvate(water_model, waters, ion_list, isAmber16):
    """
    Solvate the system to a certain number of waters.
    :param target_waters: the target number of waters in the window
    :param water_box: type of water box (must be a string that tleap recognizes)
    :param ions: type of counterions to be added for neutralization
    :return:
    """
    print('Solvating ...')

    water_box = water_model.upper()+'BOX'
    # Update target_waters: the pdb contains one water already,so we add one less water than the specified number 
    target_waters = waters - 1 

    # Initialize a log file to record the thickness of the water buffer
    # How much we should adjust the buffer in each iteration, trying to match the correct number of waters
    buffer_adjustment = 0.1  # Angstroms
    buffer_thickness = 9.0  # Angstroms

    # Call tleap for an initial estimate
    write_leap_input(water_model, buffer_thickness, isAmber16)
    current_waters = run_leap_and_read_output()

    count = 0
    max_count = 100
    refinement_threshold = 10

    while current_waters != target_waters:
	count += 1
        if count > max_count:
	# Try a larger refinement threshold
             refinement_threshold += 1
        # If we are really close, redefine the target_adjustment to make smaller steps.
        manual_removal = None
        if abs(current_waters - target_waters) < 20:
            buffer_adjustment = 0.01
        # If we are really really close and we have more waters than necessary, we can manually remove a couple
        if current_waters > target_waters and (current_waters - target_waters) < refinement_threshold:
            difference = current_waters - target_waters
            manual_removal = [target_waters + 1 + i for i in range(difference)]
            write_leap_input(water_model, buffer_thickness, isAmber16, manual_removal)
            run_leap_and_read_output()
            break
        # If we have more waters than we want, reduce the buffer thickness
        if current_waters > target_waters:
            buffer_thickness -= buffer_adjustment
        # Otherwise, we need more waters, increase the buffer thickness
        if current_waters < target_waters:
            buffer_thickness += buffer_adjustment
        write_leap_input(water_model, buffer_thickness,isAmber16)
        current_waters = run_leap_and_read_output()	

    # Compose the final leap input file from the solvate leap input stub and the final number of waters
    shutil.copy('tleap.in', 'solvate_tleap.in')
    sol_file = open('solvate_tleap.in', 'a')
    if isAmber16!='yes' and water_model!='TIP3P':
        sol_file.write('# load the water parameters\n')
        sol_file.write('loadamberparams frcmod.%s\n'%(water_model.lower()))
    elif isAmber16=='yes':
        sol_file.write('# load the water parameters\n')        
        sol_file.write('source leaprc.water.%s\n'%(water_model.lower()))

    sol_file.write('solvatebox WAT ' + water_box + ' {10.0 10.0 ' + str(buffer_thickness) + '}\n\n')
    if manual_removal is not None:
        for water in manual_removal:
            sol_file.write('remove WAT WAT.%s\n' % water)
        sol_file.write('\n')
    # Make the system neutral by putting a final '0'
    sol_file.write('\n')
    sol_file.write('desc model\n')
    sol_file.write('savepdb WAT solvated.pdb\n')
    sol_file.write('saveamberparm WAT solvated.prmtop solvated.inpcrd\n')
    sol_file.write('quit')
    sol_file.close()
    p = sp.call('tleap -s -f solvate_tleap.in > solvation.log', shell=True)
    
    f = open('solvation.log', 'r')
    added_waters = target_waters
    for line in f:
        if "Could not open file" in line:
           print 'WARNING!!!'
           print line
	   sys.exit(1)
        if "WARNING: The unperturbed charge of the unit:" in line:
           print line
	   # Program not terminated in case the user actually prefers a non-neutralized system,for whatever reason.	
           print ('The system is not neutralized properly. Are you sure you want to continue?')
           print('Maybe you should check the solvation.log file, and ion types specified in the APR input file.')
        if "addIonsRand: Argument #2 is type String must be of type: [unit]" in line:
	   print('Aborted.The ion types specified in the APR input file could be wrong.')	           
           print('Please check the solvation.log file, and the ion types specified in the APR input file.\n')
           sys.exit(1)
    f.close()

    # Clean up
    sp.call('rm tmp.log', shell=True)
    sp.call('rm tmp_tleap.in', shell=True)
    sp.call('rm leap.log', shell=True)


def write_leap_input(water_model, local_buff, isAmber16, manual_removal=None):
    """
    Write a tleap input file that solvates according to the water type and buffer size in Angstroms.
    :param water_box: type of water box (must be a string that tleap recognizes)
    :param local_buff: size of z-axis buffer in Angstroms
    :param manual_removal: can be a list to tleap to move a few, specific waters.
    :return:
    """

    # Copy tleap.in to tmp_tleap.in
    shutil.copy('tleap.in', 'tmp_tleap.in')
    tmp_file = open('tmp_tleap.in', 'a')
    if isAmber16!='yes' and water_model!='TIP3P':
        tmp_file.write('loadamberparams frcmod.%s\n'%(water_model.lower()))
    elif isAmber16=='yes':
        tmp_file.write('source leaprc.water.%s\n'%(water_model.lower()))
    tmp_file.write('solvatebox WAT ' + water_model+'BOX' + ' {10.0 10.0 ' + str(local_buff) + '}\n')
    if manual_removal is not None:
        for water in manual_removal:
            tmp_file.write('remove model model.%s' % water)
        tmp_file.write('\n')
    tmp_file.write('quit')
    tmp_file.close()


def run_leap_and_read_output():
    """
    Run tleap and get the number of waters for a given buffer padding.
    :return:
    """
    p = sp.call('tleap -s -f tmp_tleap.in > tmp.log', shell=True)
    # Search the key word "Added" in tmp.log
    current_waters = None
    f = open('tmp.log', 'r')
    for line in f:
        if "Added" in line:
            current_waters = line[line.find('Added') + len('Added'):].split()
            current_waters = int(current_waters[0])
    f.close()
    if current_waters is None:
        print('It appears there was a problem running tleap.')
  	print('Maybe check your AMBER version setting in your APR input file.\n')
        sys.exit(1)
    return current_waters



