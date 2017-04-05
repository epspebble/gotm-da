# This code is written in Python 3, but for Python 2.6 or above, we need...
from __future__ import print_function
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

### Global settings
import os
from datetime import datetime
from netCDF4 import Dataset
from numpy import pi, cos, sin

# For our netCDF4 files
epoch = datetime(1981,1,1)

## For general GOTM setup

# Set the default GOTM executable and namelist locations.
userhome = os.getenv('HOME')
project_folder = os.path.join(userhome,'medsea_GOTM')
GOTM_executable = os.path.join(project_folder,'bin','gotm')
if not(os.path.isfile(GOTM_executable)):
    raise FileNotFoundError("The GOTM executable not found at " + GOTM_executable)
GOTM_nml_path = os.path.join(project_folder,'config')
GOTM_nml_list = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp']
for nml in GOTM_nml_list:
    GOTM_nml_template = os.path.join(GOTM_nml_path,nml)
    if not(os.path.isfile(GOTM_nml_template)):
           raise FileNotFoundError("The GOTM config namelist " + GOTM_nml_template + " is invalid.")

# GOTM namelist values
timestep = 30
max_depth = 150
            
## For medsea simulations

# Top-level project folders
data_folder = os.path.join('/global/scratch',os.getenv('USER'))
base_folder = os.path.join(data_folder,'medsea_GOTM')

if not(os.path.isdir(base_folder)):
    raise IOError('The base folder: ' + base_folder + ' is either not accessible or created.')
run_folder = os.path.join(base_folder,'run')
if not(os.path.isdir(run_folder)):
    os.mkdir(run_folder)

# Ocean and Satellite products datasets source folders.
p_sossta_folder = os.path.join(data_folder,'p_sossta')
ERA_folder = os.path.join(p_sossta_folder,'medsea_ERA-INTERIM','3-hourly')
rea_folder = os.path.join(p_sossta_folder,'medsea_rea')

# GOTM dat files' netCDF reformatted dataset sources.
def data_sources(year=None, month=None, mode='r', dat=['heat','met','tprof','sprof']):
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
        NCDataset = Dataset

    if isinstance(dat,str):
        dat = [dat] # So that list comprehension still works.

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

    # print(suffix) # debug
    fn_dict = {'heat' : os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_heat' + suffix),
               'met'  : os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_met' + suffix),
               'tprof': os.path.join(data_folder,'medsea_rea','medsea_rea_votemper' + suffix),
               'sprof': os.path.join(data_folder,'medsea_rea','medsea_rea_vosaline' + suffix),
               'sst'  : os.path.join(data_folder,'medsea_OSTIA','medsea_OSTIA_sst' + suffix),
               'chlo' : os.path.join(data_folder,'medsea_MODIS','medsea_MODIS_chlor_a' + suffix)}

    assert all([each in fn_dict.keys() for each in dat]) # Check that the function is called correctly.

    for each in fn_dict.keys():
        try:
            ds_dict = {each : NCDataset(fn_dict[each],mode) for each in dat}
        except OSError:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            raise
        except:
            raise

    if len(ds_dict.keys()) == 1:
        # If only one dataset requested, return the netcdf dataset unwrapped from the dict. 
        return ds_dict[dat[0]] 
    else:
        # Else return a dictionary of datasets.
        return ds_dict

# Global setting for the core_dat() routines (and possibly the ERA routines as well)
overwrite=True

# Our grid points. Maybe we can reduce dependence on numpy by just using a Python array.
medsea_lats = tuple(30.75+0.75*i for i in range(21))
medsea_lons = tuple(-6.0+0.75*i for i in range(57))

# The corresponding index ranges in medsea_ERA-INTERIM datasets.
ERA_lat_ind = slice(-8,4,-1)
ERA_lon_ind = slice(12,-18,1)

# The corresponding index ranges in medsea_rea datasets.
rea_lat_ind = slice(9,250,12)
rea_lon_ind = slice(0,673,12)
rea_depth_ind = slice(0,18)

# Enumerate the grid points 
import itertools
mm, nn = zip(*itertools.product(range(21),range(57)))
# Make use of a medsea_rea dataset to infer sea, shallow and land locations.
with data_sources(2014,1,dat='tprof') as rea_ds:
    votemper = rea_ds['votemper']
    # Preallocate
    is_sea = list(None for i in range(21*57))
    is_shallow = list(None for i in range(21*57))
    is_land = list(None for i in range(21*57))
    for i in range(21*57):
        # Since fill value is 1e20, and sea water should not be boiling...
        is_sea[i] = (votemper[0,-1,mm[i],nn[i]]<max_depth) # deepest location in our data should be about 100m.
        is_land[i] = (votemper[0,0,mm[i],nn[i]]>max_depth) # shallowest data
        is_shallow[i] = \
            (votemper[0,0,mm[i],nn[i]]<max_depth) and \
            (votemper[0,-1,mm[i],nn[i]]>max_depth)

# Check that there are no logical loopholes.
assert sum(is_sea) + sum(is_land) + sum(is_shallow) == 21*57

# Return the counters i for lat/lon index arrays mm and nn.
sea_i = tuple(itertools.compress(range(21*57),is_sea))
land_i = tuple(itertools.compress(range(21*57),is_land))
shallow_i = tuple(itertools.compress(range(21*57),is_shallow))

# Return the actual lat/lon index pairs (m,n)
sea_mn = tuple((mm[i],nn[i]) for i in sea_i)
shallow_mn = tuple((mm[i],nn[i]) for i in shallow_i)
land_mn = tuple((mm[i],nn[i]) for i in land_i)

# Return the actual (lat,lon) coorindates as well
sea_locations = tuple((medsea_lats[m],medsea_lons[n]) for (m,n) in sea_mn)
shallow_locations = tuple((medsea_lats[m],medsea_lons[n]) for (m,n) in shallow_mn)
land_locations = tuple((medsea_lats[m],medsea_lons[n]) for (m,n) in land_mn)
