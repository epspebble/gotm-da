# This code is written in Python 3, but for Python 2.6 or above, we need...
from __future__ import print_function
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# Use tabulate to beatify output if possible
try:
    from tabulate import tabulate
except ImportError:
    pass # silently ignore if it is not installed.
    
# neccesary library
import os, sys
from datetime import datetime
from netCDF4 import Dataset
from numpy import pi, cos, sin

from os import getenv, mkdir
from os.path import isfile, isdir, join

### Global settings for the p_sossta project medsea runs.
userhome = os.getenv('HOME')

## For general GOTM setup

# Global "input" variables.
project_name = 'medsea_GOTM'
GOTM_version = 'latest'
# Some essential GOTM settings overrider
#timestep = 30
#max_depth = 150
#nsave = 120
# For setting unit of combined result files.
epoch = datetime(1981,1,1)

def setup(**overrides):
    """
    Set up other global variables including: 

       project_folder, GOTM_executable, GOTM_nml_path, GOTM_nml_list

    using preset module global variables including:

       project_name, GOTM_version, epoch

    This function is automatically invoked with default values at loading. Pass a dictionary to override any of 
    the preset global variables and recompute / check the rest.

    Version date: 2017-11-16.
    """
    
    # Global "output" variables computed from the "input" global preset module variables.
    global project_folder, GOTM_executable, GOTM_nml_path, GOTM_nml_list

    # Hardcoded lists
    preset = ['project_name', 'GOTM_version', 'epoch']
    output = ['project_folder', 'GOTM_executable', 'GOTM_nml_path', 'GOTM_nml_list']

    # Split user overrides by name.
    preset_overrides = {name: overrides[name] for name in preset if name in overrides}
    output_overrides = {name: overrides[name] for name in output if name in overrides}

    # Update the preset globals.
    current_settings = globals()
    current_settings.update(preset_overrides)

    # Now compute the output global variables one by one.

    # 1.
    if 'project_folder' in output_overrides:
        project_folder = output_overrides['project_folder']
    else:
        project_folder = join(userhome,project_name) # notebooks, jobs, scripts, config templates

    # 2.
    if 'GOTM_executable' in output_overrides:
        GOTM_executable = output_overrides['GOTM_executable']
    else:
        GOTM_executable = 'gotm'+'-'+GOTM_version # Just the name.
    # The actual file is assumed to reside in join(project_folder,'bin'). Check now.
    gotm_fullfn = join(project_folder,'bin',GOTM_executable)
    if not(isfile(gotm_fullfn))
        raise FileNotFoundError("The GOTM executable cannot be found at " + gotm_fullfn)

    # 3.
    if 'GOTM_nml_path' in output_overrides:
        GOTM_nml_path = output_overrides['GOTM_nml_path']
    else:
        GOTM_nml_path = join(project_folder,'config')
    # Check.
    if not isdir(GOTM_nml_path):
        raise FileNotFoundError('The GOTM namelist path {!s} is invalid.'.format(GOTM_nml_path))

    # 4.
    if 'GOTM_nml_list' in output_overrides:
        GOTM_nml_list = output_overrides['GOTM_nml_list']
    else:
        GOTM_nml_list = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp']
    for nml in GOTM_nml_list:
        GOTM_nml_template = join(GOTM_nml_path,nml)
        if not(isfile(GOTM_nml_template)):
            raise FileNotFoundError("The GOTM config namelist file: " + GOTM_nml_template + " is invalid.")

    print('Current pyGOTM settings:\n')
    if 'tabulate' not in sys.modules.keys():
        for name in preset:
            print('\t',name,':\t', current_settings[name])
            for name in output:
                print('\t',name,':\t', current_settings[name])
    else:
        print('Settings:')
        print(tabulate([[name,current_settings[name]] for name in preset]))
        print('Paths:')
        print(tabulate([[name,current_settings[name]] for name in output]))
        
# Set default values.
setup()
