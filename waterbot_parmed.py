#!/usr/bin/env python2
import sys
import os
from parmed.tools import changeLJSingleType
from parmed.tools import parmout
from parmed.tools import change
from parmed.amber import AmberParm
from parmed.structure import Structure

def parmed_topology(top_file, param_type, param_list,TIP3P_param_list, index):
    """
    Use ParmEd to edit water parameters
    """
    parm = AmberParm(top_file)    
    local_param_list = TIP3P_param_list
    ow_charge_column = hw_charge_column = ow_rad_column = ow_eps_column = 999

    # Write a log file for parmEd
    f = open('waterbot_parmed.log','w')
    logfile = open('parmed.log','w')
    for i, nonbonded_term in enumerate(param_type):
      if 'charge' in nonbonded_term.lower() and 'ow' in nonbonded_term.lower(): 
        ow_charge_column = i
        f.write('The partial charge of water oxygen will be perturbed; ')
        f.write('the value of the new parameter is %.4f.\n'%(param_list[i][index]))
        local_param_list[2] = param_list[i][index] 
      elif 'charge' in nonbonded_term.lower() and 'hw' in nonbonded_term.lower():
        hw_charge_column = i
        f.write('The partial charge of water hydrogen will be perturbed; ')
        f.write('the value of the new parameter is %.4f.\n'%(param_list[i][index]))
        local_param_list[3] = param_list[i][index]
      elif 'radius' in nonbonded_term.lower():
        ow_rad_column = i
        f.write('The radii parameter of water hydrogen will be perturbed; ')
        f.write('the value of the new parameter is %.4f.\n'%(param_list[i][index]))
        local_param_list[0] = param_list[i][index]
      elif 'epsilon' in nonbonded_term.lower():
        ow_eps_column = i 
        f.write('The epsilon parameter of water hydrogen will be perturbed; ')
        f.write('the value of the new parameter is %.4f.\n'%(param_list[i][index]))
        local_param_list[1] = param_list[i][index]   

    # if only one of the charges were provided:
    if ow_charge_column == 999 and hw_charge_column!=999:
      f.write('\nThe charge of water oxygen was not provided.')
      local_param_list[2] = -2.0 * local_param_list[3] 
      f.write('A value of %.4f was computed and used to make sure this water model is neutral.'%(local_param_list[2]))
    elif hw_charge_column == 999 and ow_charge_column!=999:
      f.write('\nThe charge of water hydrogen was not provided.')
      local_param_list[3] = -local_param_list[2]/2.0
      f.write('A value of %.4f was computed and used to make sure this water model is neutral.'%(local_param_list[3]))
               
    # Check whether the water model has a neutral charge
    if (local_param_list[2] + 2*local_param_list[3]) != 0:
       f.write('\nAborted.The new water model is not neutral!!!\n')
       sys.exit(1)          
     
    if ow_rad_column!=999 or ow_eps_column!=999:
      # It looks like there is no way to only change radius or epsilon
      action = changeLJSingleType(parm, "@%OW", local_param_list[0], local_param_list[1])
      action.execute()
      logfile.write(('%s\n' % action))
    if ow_charge_column!=999 or hw_charge_column!=999:
      action = change(parm, 'CHARGE', "@%OW", local_param_list[2])
      action.execute()
      logfile.write(('%s\n' % action))              
    Structure.save(parm,'solvated_perturbed.prmtop') 


