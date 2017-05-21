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

# For our netCDF4 files
epoch = datetime(1981,1,1)

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
            
## For medsea simulations

# Top-level project folders
data_folder = os.path.join('/global/scratch',os.getenv('USER'),'medsea_data')
while not(os.path.isdir(data_folder)):
    print('The data folder ' + data_folder + 'is either not accessible or created.')
    data_folder = input("Enter new data folder location.")
    
base_folder = os.path.join(data_folder,'medsea_GOTM')
while not(os.path.isdir(base_folder)):
#    raise IOError('The base folder: ' + base_folder + ' is either not accessible or created.')
    print('The base folder: ' + base_folder + ' is either not accessible or created.')
    base_folder = input("Enter new folder location.")
    
# run_folder = os.path.join(base_folder,'run')
# if not(os.path.isdir(run_folder)):
#     os.mkdir(run_folder)

# Ocean and Satellite products datasets source folders.
p_sossta_folder = os.path.join(data_folder,'p_sossta')
ERA_folder = os.path.join(p_sossta_folder,'medsea_ERA-INTERIM','3-hourly')
rea_folder = os.path.join(p_sossta_folder,'medsea_rea')

# GOTM dat files' netCDF reformatted dataset sources.
def data_sources(year=None, month=None, mode='r', dat=['heat','met','tprof','sprof','chlo'], region='medsea', grid='144x'):
    """ 
    Return the netCDF4 Dataset (or MFDataset) handles for the data source. 
    Calling data_sources() returns MFDataset of all available data for dat = ['heat','met','tprof','sprof']
    by default in read-only mode.
    
    """
    from netCDF4 import Dataset, MFDataset
    import os

    if (year is None) or (month is None):
        MF = True 
        # Need to use MFDataset
        NCDataset = MFDataset
    else:
        # When both year and month is given, we can specifically return handle to our reformatted datasets 
        # organized by months.
        NCDataset = Dataset

    if isinstance(dat,str):
        # if only one type is requested, still make it into a list so that list comprehension still works.
        dat = [dat]

    # nc_dict = dict(heat = MFDataset(os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_*.nc')), 
    #                met = MFDataset(os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_*.nc')), 
    #                tprof = MFDataset(os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_*.nc')), 
    #                sprof = MFDataset(os.path.join(data_folder,'medsea_rea','medsea_rea_vosaline_*.nc')))

    suffix = '_'
    if (year is not None):
        suffix += '{:d}'.format(year)
    else:
        suffix += '*'
    if (month is not None):
        suffix += '{:02d}'.format(month)
    else:
        suffix += '*'
    suffix += '.nc'

    fn_dict = {'heat' : os.path.join(data_folder,grid,region+'_ERA-INTERIM',region+'_ERA_heat' + suffix),
               'met'  : os.path.join(data_folder,grid,region+'_ERA-INTERIM',region+'_ERA_met' + suffix),
               'tprof': os.path.join(data_folder,grid,region+'_rea',region+'_rea_votemper' + suffix),
               'sprof': os.path.join(data_folder,grid,region+'_rea',region+'_rea_vosaline' + suffix),
               'sst'  : os.path.join(data_folder,grid,region+'_OSTIA',region+'_OSTIA_sst' + suffix),
               'chlo' : os.path.join(data_folder,grid,region+'_MODIS',region+'_MODIS_chlor_a' + suffix)}

    assert all([each in fn_dict.keys() for each in dat]) # Check that the function is called correctly.

    for each in fn_dict.keys():
        try:
            ds_dict = {each : NCDataset(fn_dict[each],mode) for each in dat}
        except OSError:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            raise
        except:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            raise

    if len(ds_dict.keys()) == 1:
        # If only one dataset requested, return the netcdf dataset unwrapped from the dict. 
        return ds_dict[dat[0]] 
    else:
        # Else return a dictionary of datasets.
        return ds_dict

# Global setting for the core_dat() routines (and possibly the ERA routines as well)
overwrite=True
