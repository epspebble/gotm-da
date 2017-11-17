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
        project_folder = os.path.join(userhome,project_name) # notebooks, jobs, scripts, config templates

    # 2.
    if 'GOTM_executable' in output_overrides:
        GOTM_executable = output_overrides['GOTM_executable']
    else:
        GOTM_executable = os.path.join(project_folder,'bin','gotm'+'-'+GOTM_version)
    # Check. 
    if not(os.path.isfile(GOTM_executable)):
        raise FileNotFoundError("The GOTM executable not found at " + GOTM_executable)

    # 3.
    if 'GOTM_nml_path' in output_overrides:
        GOTM_nml_path = output_overrides['GOTM_nml_path']
    else:
        GOTM_nml_path = os.path.join(project_folder,'config')
    # Check.
    if not os.path.isdir(GOTM_nml_path):
        raise FileNotFoundError('The GOTM namelist path {!s} is invalid.'.format(GOTM_nml_path))

    # 4.
    if 'GOTM_nml_list' in output_overrides:
        GOTM_nml_list = output_overrides['GOTM_nml_list']
    else:
        GOTM_nml_list = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp']
    for nml in GOTM_nml_list:
        GOTM_nml_template = os.path.join(GOTM_nml_path,nml)
        if not(os.path.isfile(GOTM_nml_template)):
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
        
## Global config for medsea simulations 
# The following are not meant to be changed interactively, as functions in other modules in thisos.listdir(ms.base_folder)
# package will not see the changes. 
#
# Configs that are meant to be changed in runtime should go into medsea.py.
#

# List of run_profiles done in the past.

# TODO: We should specify the corresponding GOTM code version in the profiles as well.
run_profiles = {
    # The V2 runs. They were run at 150m deep with 150 levels.
    'ASM0': dict(assimilation_type=0, extinct_method=9),
    'ASM1': dict(assimilation_type=0, extinct_method=12),
    'ASM2': dict(assimilation_type=2, assim_window=1, extinct_method=9),
    'ASM3': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat',),

    # Alternative vertical grid runs.
    # Maybe we should dynamically calculate 'nlev' and 'depth' from number of lines in the 'grid.dat' files.
    'ASM3-100m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat',
                      depth = 99.282236525788903, nlev = 132,
                      grid_method = 2, grid_file = 'grid_100m.dat'), 
    'ASM3-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),
    'ASM3-MFC-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
                         depth = 74.389762997627258, nlev = 15,
                         grid_method = 2, grid_file = 'grid_MFC_75m.dat'),
    
    # Added for completion and comparison to ASM3-75m
    'ASM0-75m': dict(assimilation_type=0, extinct_method=9, 
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),
    'ASM1-75m': dict(assimilation_type=0, extinct_method=12,
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),
    'ASM2-75m': dict(assimilation_type=2, assim_window=1, extinct_method=9, 
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),
    # 20171031 Previously ASM4-75m uses extinct_method=13, but the GOTM code
    # has not incorporated this case and it fell through to the default case
    # which is extinct_method=1, Jerlev (1976) water type 1.
    'ASM4-75m': dict(assimilation_type=2, assim_window=1,
                     extinct_method=1,  
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),
    # 20170912 Hybrid grid
    'ASM3-HYB-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
                         depth = 74.389762997627258, nlev = 41,
                         grid_method = 2, grid_file = 'grid_HYB_75m.dat'),
    # 20171031 Two levels as per email with Sam Pimentel, 2017-10-27.
    # This is identical to the previous ASM4-75m, now should properly use iop.dat
    'ASM5-75m': dict(assimilation_type=2, assim_window=1,
                     extinct_method=13, extinct_file='iop.dat', 
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),
    # This is identical to ASM3-75m, the difference lies in the GOTM code. The fixed GOTM
    # code has the albedo issue fixed, so that Payne's albedo is not applied to Ohlmann-Siegel (2000)
    # formulas which implicitly include albedo.
    'ASM3.1-75m': dict(assimilation_type=2, assim_window=1, 
                       extinct_method=12, extinct_file='chlo.dat', 
                       depth = 74.539324233308434, nlev = 122,
                       grid_method = 2, grid_file = 'grid_75m.dat'),

    # 20171106, same as the previous, except with a higher resolution output.
    'ASM5-75m_4t': dict(assimilation_type=2, assim_window=1,
                        extinct_method=13, extinct_file='iop.dat', 
                        depth = 74.539324233308434, nlev = 122,
                        grid_method = 2, grid_file = 'grid_75m.dat',
                        nsave=30),
    'ASM3.1-75m_4t': dict(assimilation_type=2, assim_window=1,
                          extinct_method=12, extinct_file='chlo.dat', 
                          depth = 74.539324233308434, nlev = 122,
                          grid_method = 2, grid_file = 'grid_75m.dat',
                          nsave=30),
    # 20171113, same as the previous, except with a higher resolution output.
    'ASM5-75m_no_salt': dict(assimilation_type=2, assim_window=1,
                               s_prof_method=0, # skip s_prof 
                               extinct_method=13, extinct_file='iop.dat', 
                               depth = 74.539324233308434, nlev = 122,
                               grid_method = 2, grid_file = 'grid_75m.dat'),
    'ASM3.1-75m_no_salt': dict(assimilation_type=2, assim_window=1,
                                 s_prof_method=0, # skip s_prof 
                                 extinct_method=12, extinct_file='chlo.dat', 
                                 depth = 74.539324233308434, nlev = 122,
                                 grid_method = 2, grid_file = 'grid_75m.dat'),
    'ASM2-75m_no_salt': dict(assimilation_type=2, assim_window=1, extinct_method=9, 
                               depth = 74.539324233308434, nlev = 122,
                               grid_method = 2, grid_file = 'grid_75m.dat'),

}
