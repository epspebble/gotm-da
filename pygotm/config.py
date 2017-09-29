# This code is written in Python 3, but for Python 2.6 or above, we need...
from __future__ import print_function
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

### Global settings for the p_sossta project medsea runs.
import os
from datetime import datetime
from netCDF4 import Dataset
from numpy import pi, cos, sin

## For general GOTM setup

# Set the default GOTM executable and namelist locations.
userhome = os.getenv('HOME')
project_folder = os.path.join(userhome,'medsea_GOTM') # notebooks, jobs, scripts, config templates
GOTM_executable = os.path.join(project_folder,'bin','gotm')
if not(os.path.isfile(GOTM_executable)):
    raise FileNotFoundError("The GOTM executable not found at " + GOTM_executable)
GOTM_nml_path = os.path.join(project_folder,'config')
GOTM_nml_list = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp']
for nml in GOTM_nml_list:
    GOTM_nml_template = os.path.join(GOTM_nml_path,nml)
    if not(os.path.isfile(GOTM_nml_template)):
           raise FileNotFoundError("The GOTM config namelist " + GOTM_nml_template + " is invalid.")

# GOTM namelist values, copied here manually as global variables.
timestep = 30
#max_depth = 150
            
## Global config for medsea simulations 
# The following are not meant to be changed interactively, as functions in other modules in thisos.listdir(ms.base_folder)
# package will not see the changes. 
#
# Configs that are meant to be changed in runtime should go into medsea.py.
#

# List of run_profiles done in the past.
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

    # 20170727 New test runs
    'ASM4-75m': dict(assimilation_type=2, assim_window=1,
                     extinct_method=13, extinct_file='iop.dat', 
                     depth = 74.539324233308434, nlev = 122,
                     grid_method = 2, grid_file = 'grid_75m.dat'),

    # 20170912 Hybrid grid
    'ASM3-HYB-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
                         depth = 74.389762997627258, nlev = 41,
                         grid_method = 2, grid_file = 'grid_HYB_75m.dat'),
}

epoch = datetime(1981,1,1)
