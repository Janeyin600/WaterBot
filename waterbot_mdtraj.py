#!/usr/bin/env python2
import glob as glob
import os
import sys
import subprocess as sp
import mdtraj as md
import math
import numpy as np
import mdtraj.utils.unit.unit_definitions as u

convert_pressure =   1.45837843401E-5 # This factor converts the unit atm to kcal/mol/Angstrom^3
GAS_CONSTANT     =   8.314459848
kj2kcal          =   0.239006 # This factor converts kj/mol to kcal/mol
nm3toA3          =   1000
j2kj             =   0.001    

# Correction terms for computing enthalpy of vaporization, obtained from the classic TIP3P-EW paper:Horn et al. J.Chem. Phys. 120, 9665 (2004)
Cvib             =  -0.0651    # in Kcal/mol
Cni              =  -0.0048    # in Kcal/mol
Cx               =   0         # This correction term is neglectable    

# Partial charge of the TIP3P water models. hand coded for now, will be removed soon for charge pertubing functionality!
hw_charge        =   0.417
ow_charge        = - 0.834 

def mdtraj_analysis(quiet, bulk, water_index, summary_file, local_file, cycle, waters, temp, output_freq, steps, stepsize):
    """
    Using mdtraj to analyze the Amber trajectory
    """
    print 'Generating', local_file, '...'
    # Append the local report for each water model
    if quiet == 'no':
      f_local = open(local_file,'a')    
      f_local.write('###### Temperature (unit: K) ######\n')  

    # Insert the final results to the summary file
    f_summary = open(os.path.dirname(__file__) + '/%s'%(summary_file),'a')
    f_summary.write('\n%-8d'%(water_index))

    frames_per_iteration = int(steps/output_freq)
    
    temp_list    = []
    energy_list  = []
    density_list = []
    volume_list  = []

    # Read openmm output files
    for iteration in range(cycle): 
      with open('traj.%02d.out'%(iteration+1)) as f:
        lines = f.read().splitlines()
        lines.pop(0)
        new_temp = [float(line.split(',')[2]) for line in lines]
        new_energy = [float(line.split(',')[1]) for line in lines]
        new_volume = [float(line.split(',')[3]) for line in lines]
        temp_list += new_temp
        energy_list += new_energy
        volume_list += new_volume
        if quiet == 'no':
          average_temp = sum(temp_list)/len(temp_list)
          f_local.write('After iteration %03d: %10.4f\n'%(iteration+1, average_temp))

    average_temp2 = np.array(temp_list, float)
    f_summary.write('%6.2f'%(np.mean(average_temp2)))
   
    ##############################################################
    if quiet == 'no':
      f_local.write('###### Density (unit: g/cm^3) ######\n') 

    for iteration in range(cycle):
        with open('traj.%02d.out'%(iteration+1)) as f:
          lines = f.read().splitlines()
          lines.pop(0)
       	  new_density = [float(line.split(',')[4]) for line in lines]
       	  density_list +=	new_density
          average_density = sum(density_list)/len(density_list)
          if quiet == 'no':
            f_local.write('After iteration %03d: %10.4f\n'%(iteration+1, average_density))
    
    final_density = np.mean(density_list)     
    f_summary.write('%10.3f'%(final_density))
      
    ###########################################################
    if quiet == 'no':
      f_local.write('###### Dielectric Constant ######\n')
    
    charges = []
    for i in range(0,waters):
      charges.append(ow_charge)
      charges.append(hw_charge)
      charges.append(hw_charge) 

    # Compute dielectric constant
    if os.path.isfile('traj.01.h5'):
      t = md.load('traj.01.h5')

    if quiet == 'no':
      dielectric_constant = md.static_dielectric(t, charges, temp)
      f_local.write('After iteration 001: %10.4f\n'%(dielectric_constant))

    for iteration in range(1,cycle):
      new_t = md.load('traj.%02d.h5'%(iteration+1))
      t = t.join(new_t,True)
      if quiet == 'no':
        dielectric_constant = md.static_dielectric(t, charges, temp)
        f_local.write('After iteration %03d: %10.4f\n'%(iteration+1, dielectric_constant))

    if quiet == 'yes':
      dielectric_constant = md.static_dielectric(t, charges, temp)

    f_summary.write('%10.3f'%(dielectric_constant))
      
    ###########################################################
    if quiet == 'no':
      f_local.write('### Isothermal compressibility (unit: 10^(-6) bar^(-1)) ###\n')
    # Compute isothermal compressibility
    if os.path.isfile('traj.01.h5'):
      t = md.load('traj.01.h5')
    if quiet == 'no':
      ic = md.isothermal_compressability_kappa_T(t, temp)*1000000
      f_local.write('After iteration 001: %10.4f\n'%(ic))
    
    for iteration in range(1,cycle):
      new_t = md.load('traj.%02d.h5'%(iteration+1))
      t = t.join(new_t,True)
      if quiet == 'no':
        ic = md.isothermal_compressability_kappa_T(t, temp)*1000000
        f_local.write('After iteration %03d: %10.4f\n'%(iteration+1, ic))

    if quiet == 'yes':
      ic = md.isothermal_compressability_kappa_T(t, temp)*1000000

    f_summary.write('%10.3f'%(ic))
  
    ##########################################################################
    if quiet == 'no':
      f_local.write('### Thermal Expansion Coefficient (unit:10^(-4)K^-(1))###\n')

    # Compute thermal expansion coefficient
    if os.path.isfile('traj.01.h5'):
      t = md.load('traj.01.h5')
      if quiet == 'no':
        thermal_expansion = 10000.0*thermal_expansion_alpha_P(t, temp, energy_list[0:frames_per_iteration])
        f_local.write('After iteration 001: %10.4f\n'%(thermal_expansion))

    for iteration in range(1,cycle):
      new_t = md.load('traj.%02d.h5'%(iteration+1))
      t = t.join(new_t,True)
      if quiet == 'no':
        thermal_expansion = 10000*thermal_expansion_alpha_P(t, temp,  energy_list[0: (iteration+1)*frames_per_iteration])  
        f_local.write('After iteration %03d: %10.4f\n'%(iteration+1, thermal_expansion))
    
    if quiet == 'yes':
      thermal_expansion = 10000*thermal_expansion_alpha_P(t, temp,  energy_list)

    f_summary.write('%10.3f'%(thermal_expansion))          

    ###################################################################
    if quiet == 'no':
      f_local.write('###### Heat of Vaporization (unit: kcal/mol)######\n')
    # Compute enthalpy of vaporization
    average_energy = np.array(energy_list, float)
    average_volume = np.array(volume_list, float)

    if quiet == 'no':
      for iteration in range(cycle):
        accumulative_energy = np.array(energy_list[0:(iteration+1)*frames_per_iteration],float)
        accumulative_volume = np.array(volume_list[0:(iteration+1)*frames_per_iteration],float)
        
        # See the TIP4P-EW paper by Horn et al., or the SI of OPC paper by Izadi and Onufriev
        heat = - np.mean(accumulative_energy) * kj2kcal / float(waters) + GAS_CONSTANT * temp * kj2kcal * j2kj\
             - np.mean(accumulative_volume) * nm3toA3 * convert_pressure/float(waters) + Cvib + Cni + Cx # Epol is not counted for now    
        f_local.write('After iteration %03d: %10.4f\n'%(iteration+1, heat))                  
    
    final_heat = - np.mean(average_energy) * kj2kcal / float(waters) + GAS_CONSTANT * temp * kj2kcal * j2kj\
             - np.mean(average_volume) * nm3toA3 * convert_pressure/float(waters) + Cvib + Cni + Cx
    # The Cpol term is not included here because:
    # "There is a bit of a philosophical question about whether to include it. -- Niel Henriksen"
    f_summary.write('%10.3f'%(final_heat)) 
    
    ### compute a water score
    water_score = float(bulk[0][3])*abs(final_density - bulk[0][1])/bulk[0][1] + float(bulk[1][3])*abs(final_heat - bulk[1][1])/bulk[1][1]\
                  + float(bulk[2][3])*abs(dielectric_constant - bulk[2][1])/bulk[2][1] + float(bulk[3][3])*abs(ic - bulk[3][1])/bulk[3][1]\
                  + float(bulk[4][3])*abs(thermal_expansion - bulk[4][1])/bulk[4][1]  
    f_summary.write('%10.1f'%(10/water_score))
     
    f_summary.close()
    
    if quiet == 'no':
      f_local.close()

def thermal_expansion_alpha_P(traj, temperature, energies):
    """
    Obtained from the source code of MDtraj
    """    
    gas_constant = GAS_CONSTANT * u.joule / u.kelvin / u.mole
    temperature = temperature * u.kelvin  
    mean_volume = traj.unitcell_volumes.mean()

    alpha = np.cov(traj.unitcell_volumes, energies)[0, 1]  # <HV> - <H><V> = cov(H, V)
    alpha /= mean_volume
    alpha *= u.kilojoules_per_mole
    alpha /= (gas_constant * temperature ** 2)

    return alpha * u.kelvin 
    
