#!/usr/bin/env python2
import os as os
import sys as sys
import glob as glob
import shutil as shutil
import datetime as dt
import signal as signal
import subprocess as sp
import re as re

import simtk.openmm.version as simtk_openmm_version
import waterbot_solvate as waterbot_solvate
import waterbot_eq as waterbot_eq
import waterbot_mdtraj as waterbot_mdtraj
import waterbot_parmed as waterbot_parmed
import waterbot_openmm as waterbot_openmm

class WaterBot:
    def __init__(self):
        """
        Default user and system parameters specification. These are now overwritten with the values found in WaterBot.in.
        More information on these variables can be found in the input file and the online tutorial.
        """
        self.amber16 = 'yes'
        self.exe_path = 'pmemd.cuda'
        self.output_file = 'waterbot.out'
        self.input_file = 'None'
        self.param_file = 'None'
        self.param_list = [[],[],[],[]]  # a list of lists
        self.num_water_model = 1
        self.tip3p_param = [1.7683, 0.1520, -0.834, 0.417]

        # A list of bulk properties of liquid water: density, heat of vaporization, dielectric constant, isothermal_compressibility
        # and thermal expansion coefficient [property, experimental value, unit, weight (for calculating the water score) 
        self.water_properties = [['density',0.997,'g/cm^3',0],['heat of vaporization',10.52, 'kcal/mol', 0], 
        ['dielectric constant',78.4, 'n/a', 0],
        ['isothermal compressibility',45.3,'10^(-6)bar^-1', 0],['thermal expansion coefficient', 2.56, '10^(-4)K^-1', 0]]
        self.ions =[]         
       
    def check_arguments(self):
        """
        The location/name of the input file is required; output file is optional  
        :return:
        """
        if len(sys.argv) == 1:
            help_message()
            # Return without an error
            sys.exit(0)
        elif (len(sys.argv) == 2) and ('-h' in sys.argv[1].lower()):
            help_message()
            # Return without an error
            sys.exit(0)
        elif (len(sys.argv) == 7) and ('-i' in sys.argv) and ('-p' in sys.argv):
            print('Running WaterBot ...\n')
        elif (len(sys.argv) == 5) and ('-i' in sys.argv) and ('-p' in sys.argv):
            print('Running WaterBot ...') 
        elif (len(sys.argv) < 5):
            help_message()
            print('Too few arguments. You need to provide an input file with the -i and -p flags.\n')
            sys.exit(1)
        else:
            help_message()
            print('Not the right command.You need the provide an input file with the -i and -p flags.\n')
            sys.exit(1)

        # Loop over the list of command line arguments and look for the `-p`  and `-i` flags
        for i, argv in enumerate(sys.argv[:-1]):
          if argv.lower() == '-i':           
            self.input_file = sys.argv[i+1]
            if not os.path.isfile(self.input_file):
              print('The input file %s does not exist.' % self.input_file)
              sys.exit(1)
          elif argv.lower() == '-p':
            self.param_file = sys.argv[i+1]
            if not os.path.isfile(self.param_file):
              print('The water parameter file %s does not exist.' % self.param_file)
              sys.exit(1)
          elif argv.lower() == '-o':
            self.output_file = sys.argv[i+1]
          
        if self.input_file == 'None':
            # if the input file is not provided
            help_message()
            print('Please provide an input file!')
            sys.exit(1)
        
        if self.param_file == 'None':
            # if the parameter file is not provided
            help_message()
            print('Please provide a parameter file for your water model!')
            sys.exit(1)
        

    def process_input_file(self):
        """
        Process the WaterBot input file
        """
        with open(self.input_file) as f_in:
            # Delete all spaces/tabs at both ends       
            lines = (line.strip(' \t\n\r') for line in f_in)
            lines = list(line for line in lines if line)  # Non-blank lines in a list

        for i in range(0, len(lines)):
            # Combine the lines that belong to the same entry
            if not lines[i][0] == ';':
                lines[i] = lines[i].split(';')[0].split('=')
                if len(lines[i]) == 1:
                    j = i
                    while True:
                        if lines[j - 1] != ';':
                            lines[j - 1][1] += lines[i][0]
                            lines[i] = ';'
                            break
                        j -= 1
        
        for i in range(0, len(lines)):
            if not lines[i][0] == ';':
                lines[i][0] = lines[i][0].strip().lower()
                lines[i][1] = lines[i][1].strip()
                if lines[i][0] == 'amber16':
                    if lines[i][1].lower() == 'yes':
                        self.amber16 = 'yes'
                    elif lines[i][1].lower() == 'no':
                        self.amber16 = 'no'
                    else:
                        print('Wrong input! Please use yes or no to indicate whether it is AMBER16 or an older version.')
                        sys.exit(1)
                elif lines[i][0] == 'quiet':
                    if lines[i][1].lower() == 'yes':
                        self.quiet = 'yes'
                    elif lines[i][1].lower() == 'no':
                        self.quiet = 'no'
                    else:
                        print('Wrong input! Please use yes or no to indicate whether you need the bulk properties printed after each iteration:')
                        print('yes means only a summary output will be generated in the end; no means log files will be generated for each water model')
                        print('after every iteration. Choosing to be quiet will save a lot of analyzing time. However, being noisy is helpful especially')
                        print('when you try to determine how long it takes for each property to reach good convergence.\n') 
                        sys.exit(1)
                elif lines[i][0] == 'platform':
                    if lines[i][1].lower() == 'cuda':
                        self.platform = 'CUDA'
                    elif lines[i][1].lower() == 'opencl':
       	       	       	self.platform =	'OpenCL'
                    elif lines[i][1].lower() == 'cpu':
                        self.platform = 'CPU'
                    elif lines[i][1].lower() == 'reference':
                        self.platform = 'Reference'
                    else:
                        print('Wrong input! Only provide CUDA, OpenCL, CPU or Reference as the OpenMM platform options.')
                        print('CUDA and OpenCL will be my first and second best choices becuase they are significantly faster.')
                        sys.exit(1)  
                elif lines[i][0] == 'exe_path':
                    self.exe_path = lines[i][1].strip('\'\"-,.:;][')
                elif lines[i][0] == 'water_name':
                    self.water_name = lines[i][1].upper()
                elif lines[i][0] == 'water_site':
                    water_site = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                    if water_site == 3:
                       self.water_model = 'TIP3P'
                    else:
                       print('We are not supporting 4 or 5 site water model at this point.\n')
                       sys.exit(0)      
                elif lines[i][0] == 'waters':
                    self.waters = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'temperature':
                    self.temp = ismyinstance('float', lines[i][1], self.input_file, lines[i][0]) 
                elif lines[i][0] == 'dt':
                    self.stepsize = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'nstlim':
                    self.steps = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'vdw_cutoff':
                    self.vdw_cutoff = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'output_freq':
                    self.output_freq = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'cycle':
                    self.cycle = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif 'density' in lines[i][0]:
                    self.water_properties[0][3] = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif 'heat' in lines[i][0] and 'vaporization' in lines[i][0]:
                    self.water_properties[1][3] = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif 'dielectric' in lines[i][0] and 'constant' in lines[i][0]:
                    self.water_properties[2][3] = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif 'isothermal' in lines[i][0] and 'compressibility' in lines[i][0]:
                    self.water_properties[3][3] = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif 'thermal' in lines[i][0] and 'expansion' in lines[i][0]:
                    self.water_properties[4][3] = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                else:
                    print('Wrong entry name in %s! ' % self.input_file)
                    print(lines[i][0] + '\n')
                    print('Please use the same keywords as in the template input file. Aborted.')
                    sys.exit(1)

        if sum ([float(property[3]) for property in self.water_properties])!= 1:
          print('The sum of weights from all properties has to be equal to 1')  
          sys.exit(1)
        else:
          print('The sum of weights from all properties is equal to 1. Perfect!')

        for property in self.water_properties:
          print '  The weight of %s (exp val: %.3f %s) is %.1f percent'%(property[0],property[1],property[2], property[3]*100)  

    def parsing_params(self):
        """
        Parse the parameter file for the water model(s)
        """
        with open(self.param_file) as f:
          lines = f.read().splitlines()
           
        # Processing the first line, which is supposed to indicate the types of the parameters
        # listed in this file
        self.param_type = lines[0].split()
        
        # Check if every parameter type starts with a symbol '#'
        for type in self.param_type:
          if type[0]!='#':
            print('Use the symbol\'#\' to indicate what those parameters are for each column')
            sys.exit(1)

        # Get the number of water models that will be evaluate
        lines.pop(0)
        # Remove the blank lines
        lines = (line.strip(' \t\n\r') for line in lines)
        lines = list(line for line in lines if line)
        self.num_water_models = len(lines)
        if self.num_water_models == 1:
          print ('\nYour parameter file contains only one water model.')
        else:
          print ('\nYour parameter file contains  %d water models.'%(self.num_water_models))

        # Check to see if the parameter type is one of those four: radius, epsilon, charge_oxy and charge_hyd

        # A dictionary data strucutre is used, the value of each type is corresponding to the column in the parameter data file
        supported_types = {'ow_radius':999,'ow_epsilon':999,'charge_ow':999,'charge_hw':999}

        for i, type in enumerate(self.param_type):
           for key, value in supported_types.items():
              if key == type[1:].lower():
                supported_types[key] = i
    
        for key, value in sorted(supported_types.items(),reverse=True):
           if value!=999:
             print ('Reading parameters for %s, column %d ...'%(key, value))

        for i, line in enumerate(lines):
          splitline = line.split()
          # Check if the number of parameters per row matches the header
          if len(splitline) != len(self.param_type):
            print ('\nProblem detected in water model %d'%(i+1))
            print ('Expecting %d parameters, got %d'%(len(self.param_type), len(splitline)))
            sys.exit(1)

          for j in range(len(splitline)):
             self.param_list[j].append(float(splitline[j]))
     
    def solvate_system(self):
        """
        Solvate the system based on the given numebr of water molecules
        """
        ###Build a folder called simulation_tryXXX in the current directory, and then create subfolders for each water model

        if self.num_water_model > 999:
          print ('Please try less than 1000 water models at one time.')
          system.exit(1)

        if os.path.exists('simulation_try999'):
          print('Folder simulation_try999 was detected. Please clean up your directory and try again!\n')

        # Analyze the exsiting simulation folders
        simulation_folders = []

        for i in range(1,1000):
          if os.path.exists('simulation_try%03d'%(i)):
            simulation_folders.append('simulation_try%03d'%(i))

        # It looks like the folders are automately sorted, so we just need to analyze the last element/folder
        if not simulation_folders:
        # for an empty list
           curr_folder_index = 0
        else:
          curr_folder_index =  int(simulation_folders[-1][-3:])

        destination = 'simulation_try%03d'%(curr_folder_index+1)
        self.destination = destination

        print ('The simulation data will be stored in folder %s'%(destination))
        os.makedirs('%s'%(destination))

        os.chdir(self.destination)

        # Write water PDB file
        TIP3P_PDB_file = open('TIP3P.pdb','w')
        TIP3P_PDB_file.write('HETATM    1  O   WAT A 1       -21.776  35.513  28.223  1.00 59.96           O')
        TIP3P_PDB_file.close()

        if os.path.isfile('tleap.in'):
              sp.call('rm tleap.in', shell = True)

        header_file = open('tleap.in', 'a')
        if self.amber16!='yes':
          header_file.write('source leaprc.ff12SB\n')
        elif self.amber16=='yes':
          header_file.write('source leaprc.protein.ff14SB\n')
        header_file.write('source leaprc.gaff\n')
     
        if self.water_model == 'TIP3P':
          header_file.write('WAT = loadpdb TIP3P.pdb\n')
        else:
          head_file.write('WAT = loadpdb tip4p.pdb\n')
        header_file.close()

        waterbot_solvate.setup_solvate(self.water_model, self.waters, self.ions, self.amber16)
        sp.call('rm tleap.in', shell=True)

    def equilibrate_system(self):
        """
        Equilibrate the system using AMBER
        """ 
        # Walk into each subfolder of water models
        print ('### Equilibration Starts ###')

        for model in range(1,self.num_water_models+1):
          os.chdir('waterModel%03d'%(model)) 
          print('Water model %03d'%(model))
          waterbot_eq.write_min_in()
          waterbot_eq.write_therm1_in()        
          waterbot_eq.write_therm2_in(self.temp)
          waterbot_eq.write_eqnpt_in(self.temp)        
          print('Minimizing ...')
          sp.call('%s -O -i min.in -p solvated_perturbed.prmtop -c solvated.inpcrd -o min.out -r solvated.rst7 '
            '-inf mdinfo -ref solvated.inpcrd'%(self.exe_path), shell=True)
          print('Thermalizing ...')
          sp.call('%s -O -i therm1.in -p solvated_perturbed.prmtop -c solvated.rst7 -o therm1.out -r therm1.rst7 '
            '-x therm1.nc -inf mdinfo -ref solvated.inpcrd'%(self.exe_path), shell=True)
          print('Heating up ...')
          sp.call('%s -O -i therm2.in -p solvated_perturbed.prmtop -c therm1.rst7 -o therm2.out -r therm2.rst7 '
            '-x therm2.nc -inf mdinfo -ref solvated.inpcrd'%(self.exe_path), shell=True)
          print('Equilibrating using NPT ...')
          sp.call('%s -O -i eqnpt.in -p solvated_perturbed.prmtop -c therm2.rst7 -o eqnpt.out -r eqnpt01.rst7 -x eqnpt.nc '
            '-inf mdinfo -ref solvated.inpcrd'%(self.exe_path), shell=True)
       
          # Do 50 iterations of short quilibration           
          for run in range(1, 50):
            sp.call('%s -O -i eqnpt.in -p solvated_perturbed.prmtop -c eqnpt%02d.rst7 -o eqnpt.out -r eqnpt%02d.rst7 '
                '-x eqnpt.nc -inf mdinfo -ref solvated.inpcrd' % (self.exe_path, run, run + 1), shell=True)
            os.remove('eqnpt%02d.rst7' % run)
          os.chdir('../')        

    def run_md(self):
        """
        Run MD simulations for the following bulk property analysis 
        """
        # Walk into each subfolder of water models

        print ('### Production Starts ###')

        for model in range(1,self.num_water_models+1):
          os.chdir('waterModel%03d'%(model))
          print('Water model %03d'%(model))
          shutil.copyfile('eqnpt50.rst7','rst.00')
          for iter in range(self.cycle):
            print('  iteration %d'%(iter+1))
            waterbot_openmm.run_md(self.platform, iter+1, self.steps, self.stepsize, self.vdw_cutoff, self.temp, self.output_freq)
          os.chdir('../')

    def analyze_trajectory(self):
        """
        Analyze the h5 trajectory using mdtraj program
        """
        # Write a summary dat file for all water models
        """
        if os.path.isfile(self.output_file):
          shutil.copyfile(self.output_file, self.output_file+'.backup')
          os.remove(self.output_file)        
        """ 
        summary_file = open(self.output_file, 'w')
        summary_file.write('#Water Model #Temperature (K)  #Density (g/cm^3)  #Dielectric Constant  #Isothermal compressibility (10^(-6) bar^(-1)) ...\n')
        summary_file.write('... #Thermal Expansion Coefficient (10^(-4)K^(-1)) #Heat of Vaporization (kcal/mol) #Water Score\n')
        summary_file.close()
        for model in range(1,self.num_water_models+1):
          os.chdir('waterModel%03d'%(model))
          # Write a local log file for water properties
          local_output = 'water_model%03d.log'%(model)
          local_file = open(local_output, 'w')
          waterbot_mdtraj.mdtraj_analysis(self.quiet, self.water_properties, model, self.destination+'/'+self.output_file, local_output,\
                           self.cycle, self.waters, self.temp, self.output_freq, self.steps, self.stepsize)
          os.chdir('../')

    def evaluate_all_water_models(self):
        
        for model in range(1,self.num_water_models+1):
          # Start from 1 (not 0)
          os.makedirs('waterModel%03d'%(model))                 
          # Copy the topology and coordinate files to each water subfolder:
          shutil.copy('solvated.prmtop','waterModel%03d'%(model))   
          shutil.copy('solvated.inpcrd','waterModel%03d'%(model))
          os.chdir('waterModel%03d'%(model))          

          # Using the ParmEd API to edit solvated.prmtop
          waterbot_parmed.parmed_topology('solvated.prmtop', self.param_type, self.param_list, self.tip3p_param, model-1)                    
          os.chdir('../')

def ismyinstance(variable_type, parameter_value, filename, parameter):
    """
    Check the data type of the user inputs
    """ 
    if not parameter_value:
    # If the value of entry is left as blank    
         if variable_type == 'string':
             return 'None'
         elif variable_type == 'list':
             return []
         elif variable_type == 'float':
             return 0.0
         elif variable_type == 'int':
             return 0
    
    # if input was provided
    if variable_type == 'float':
        try:
            float(parameter_value)
        except ValueError:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter a float value for %s.' % parameter)
            sys.exit()

        if float(parameter_value) < 0:
            print('Wrong input in %s:' % (filename))
            print(parameter_value)
            print('Please enter a non-negative value for %s.' % parameter)
            sys.exit()
        else:
            return float(parameter_value)

    elif variable_type == 'int':
        try:
            int(parameter_value)
        except ValueError:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter an int value for %s.' % parameter)
            sys.exit()

        if int(parameter_value) < 0:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter a non-negative value for %s.' % parameter)
            sys.exit()
        else:
            return int(parameter_value)

    elif (variable_type == 'list') or (variable_type == 'string'):
        return parameter_value

def welcome_message():
    print('*************************************************************************')
    print(' Welcome to WaterBot:')
    print('     a program that calculates the bulk properties of new water models.')
    print('                    a beta version \n')
    print(' Written by: Jian (Jane) Yin')
    print(' Copyright (c) 2017-, University of California, San Diego')
    print('*************************************************************************')
    
def help_message():
    print('Usage:')
    print('-i  WaterBot input file\n')
    print('-p  WaterBot parameter file\n')
    print('-o  WaterBot output file (optional; default file name will be provided)\n')
    print('For example:')
    print('  python2 waterbot.py -i waterbot.in -p water_param.dat -o janewaterbot.out')
    print('or')
    print('  python2 waterbot.py -i waterbot.in -p water_param.dat')

def check_versions():
    """
    Check the versions of python and openmm
    """
    latest_openmm_version = '7.1'

    print 'Checking the version of OpenMM installed ...'
    openmm_version = simtk_openmm_version.short_version.split('.')
    if openmm_version[0] == latest_openmm_version.split('.')[0] and openmm_version[1] == latest_openmm_version.split('.')[1]:
      print ('You are using OpenMM %s.%s, which is the latest version.\n'%(openmm_version[0],openmm_version[1]))
    else:
      print ('You are using OpenMM %s.%s.'%(openmm_version[0],openmm_version[1]))
      print ('You can now download the lastest version of OpenMM (%s).\n'%(latest_openmm_version))
              
def main():
    welcome_message()
    check_versions()
    this = WaterBot()
    this.check_arguments()
    this.process_input_file()
    this.parsing_params()    
    this.solvate_system()    
    this.evaluate_all_water_models()
    this.equilibrate_system()
    this.run_md()  
    this.analyze_trajectory()

if __name__ == "__main__":
    main()   
