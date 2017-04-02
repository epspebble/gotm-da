# This code is written in Python 3, but for Python 2.6 or above, we need...
from __future__ import print_function
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

### Global settings
import os
from datetime import datetime
from numpy import pi, cos, sin

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
def data_sources(year, month, mode = 'r', dat = ['heat','met','tprof','sprof']):
    from netCDF4 import Dataset
    import os

    if isinstance(dat,str):
        dat = [dat] # So that list comprehension still works.

    fn_dict = {'heat' : os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_{:d}{:02d}.nc'.format(year,month)),
               'met'  : os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_{:d}{:02d}.nc'.format(year,month)),
               'tprof': os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_{:d}{:02d}.nc'.format(year,month)),
               'sprof': os.path.join(data_folder,'medsea_rea','medsea_rea_vosaline_{:d}{:02d}.nc'.format(year,month)),
               'sst'  : os.path.join(data_folder,'medsea_OSTIA','medsea_OSTIA_sst_{:d}{:02d}.nc'.format(year,month)),
               'chlo' : os.path.join(data_folder,'medsea_MODIS','medsea_MODIS_chlor_a_{:d}{:02d}.nc'.format(year,month))}
    assert all([each in fn_dict.keys() for each in dat])

    for each in fn_dict.keys():
        try:
            ds_dict = {each : Dataset(fn_dict[each],mode) for each in dat}
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
        is_sea[i] = (votemper[0,-1,mm[i],nn[i]]<100) # deepest location in our data should be about 100m.
        is_land[i] = (votemper[0,0,mm[i],nn[i]]>100) # shallowest data
        is_shallow[i] = \
            (votemper[0,0,mm[i],nn[i]]<100) and \
            (votemper[0,-1,mm[i],nn[i]]>100)

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