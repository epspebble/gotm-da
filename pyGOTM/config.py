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

# check whether are in interactive shell like ipython or jupyter notebook, which
# does not have the '__file__' attribute. Maybe all main() based runs has '__file__'?
def is_interactive():
    import __main__ as main
    return not hasattr(main,'__file__')

# Set some default flags for interactive use of the module.
verbose = True if is_interactive() else False
    
# neccesary library
import os, sys
from datetime import datetime
from netCDF4 import Dataset
from numpy import pi, cos, sin

from os import getenv, mkdir
from os.path import isfile, isdir, join

### Global settings for the project runs for the medsea.
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

def setup(verbose=verbose,**overrides):
    """
    Set up other global variables including: 

       project_folder, GOTM_executable, GOTM_executable_path, GOTM_nml, GOTM_nml_path

    using preset module global variables including:

       project_name, GOTM_version, epoch

    This function is automatically invoked with default values at loading. Pass a dictionary to override any of 
    the preset global variables and recompute / check the rest.

    Version date: 2017-11-16.
    """
    
    # Global "paths" variables computed from the "input" global preset module variables.
    global project_folder, GOTM_executable, GOTM_executable_path, GOTM_nml_path, GOTM_nml

    # Hardcoded defaults
    preset = ['project_name', 'GOTM_version', 'epoch']
    paths = ['project_folder', 'GOTM_executable', 'GOTM_executable_path', 'GOTM_nml', 'GOTM_nml_path']

    # Split user overrides by name.
    preset_overrides = {name: overrides[name] for name in preset if name in overrides}
    paths_overrides = {name: overrides[name] for name in paths if name in overrides}

    # Update the preset globals.
    current_settings = globals()
    current_settings.update(preset_overrides)

    # Now compute the output global variables one by one.

    # 1. Project folder
    if 'project_folder' in paths_overrides:
        project_folder = paths_overrides['project_folder']
    else:
        project_folder = join(userhome,project_name) # notebooks, jobs, scripts, config templates

    # 2. GOTM program
    if 'GOTM_executable' in paths_overrides:
        GOTM_executable = paths_overrides['GOTM_executable']
    else:
        GOTM_executable = 'gotm'+'-'+GOTM_version # Just the name.
    if 'GOTM_executable_path' in paths_overrides:
        GOTM_executable_path = paths_overrides['GOTM_executable_path']
    else:
        # The actual file is assumed to reside in join(project_folder,'bin'). Check now.
        GOTM_executable_path = join(project_folder,'bin')
    if not(isfile(join(GOTM_executable_path,GOTM_executable))):
        raise FileNotFoundError("The GOTM executable {!s} cannot be found at {!s} ".format(GOTM_executable, GOTM_executable_path))

    # 3. GOTM namelists
    if 'GOTM_nml_path' in paths_overrides:
        GOTM_nml_path = paths_overrides['GOTM_nml_path']
    else:
        GOTM_nml_path = join(project_folder,'nml')
    if not isdir(GOTM_nml_path):
        raise FileNotFoundError('The GOTM namelist templates path {!s} is invalid.'.format(GOTM_nml_path))
    if 'GOTM_nml' in paths_overrides:
        GOTM_nml = paths_overrides['GOTM_nml']
    else:
        GOTM_nml = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp']
    for nml in GOTM_nml:
        GOTM_nml_template = join(GOTM_nml_path,nml)
        if not(isfile(GOTM_nml_template)):
            raise FileNotFoundError("The GOTM config namelist file: " + GOTM_nml_template + " is invalid.")

    # Print config.
    if verbose:
        print('Current pyGOTM settings:')
        if 'tabulate' not in sys.modules.keys():
            for name in preset:
                print('\t',name,':\t', current_settings[name])
                for name in paths:
                    print('\t',name,':\t', current_settings[name])
        else:
            print(tabulate([[name,current_settings[name]] for name in preset]))
            print('\nPaths and filenames:')
            print(tabulate([[name,current_settings[name]] for name in paths]))
        
# Set default values.
setup()
