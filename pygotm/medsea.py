from .config import *
from .gotmks import *

import numpy as np

## Global settings and initializations.

# If a change to the following is desired, do so right after loading the module and then call set_grid() and set_folders() in this order.
ASM = 3 # This selects the GOTM extra parameters profile
max_depth = 75 
grid = '144x'

# This will also be in the name (as a suffix before the .nc extension) in the nc file per grid point.
run = 'ASM3-75m-144x' # Will be reset when set_grid() is called.

region = 'medsea' # practically just a prefix of filenames for now, the code is strongly tied to this assumption
overwrite = True # True means running at the same grid point will overwrite files if already present (notably the *.inp etc...)

# Routines to set global values in this module. Can be used in interactive session to change config.

def set_folders():
    global scratch_folder, data_folder, base_folder, run_folder, p_sossta_folder, ERA_folder, rea_folder, cache_folder
    # Top-level project folders
    scratch_folder = os.path.join(userhome,'scratch')
    # the grid subfolder is now part of the data_folder
#    data_folder = os.path.join(userhome,'medsea_data', grid)
    data_folder = os.path.join(scratch_folder,'medsea_data')
    while not(os.path.isdir(data_folder)):
        print('The data folder ' + data_folder + ' is either not accessible or created.')
        data_folder = input("Enter new data folder location.")
    base_folder = os.path.join(scratch_folder,'medsea_GOTM')
    while not(os.path.isdir(base_folder)):
        #    raise IOError('The base folder: ' + base_folder + ' is either not accessible or created.')
        print('The base folder: ' + base_folder + ' is either not accessible or created.')
        base_folder = input("Enter new folder location.")
    #run_folder = os.path.join(base_folder,run)
    run_folder = base_folder
    if not(os.path.isdir(run_folder)):
        print('Run folder: {:s} not found. Creating it now.'.format(run_folder))
        os.mkdir(run_folder)

    # Ocean and Satellite products datasets source folders.
    p_sossta_folder = os.path.join(scratch_folder,'p_sossta')
    ERA_folder = os.path.join(p_sossta_folder,'medsea_ERA-INTERIM','3-hourly')
    rea_folder = os.path.join(p_sossta_folder,'medsea_rea')
    cache_folder = '/dev/shm'
    return scratch_folder, data_folder, base_folder, run_folder, p_sossta_folder, ERA_folder, rea_folder, cache_folder

def get_folders():
    global scratch_folder, data_folder, base_folder, run_folder, p_sossta_folder, ERA_folder, rea_folder, cache_folder
    return scratch_folder, data_folder, base_folder, run_folder, p_sossta_folder, ERA_folder, rea_folder, cache_folder

# A CRUCIAL routine for parallelizing over grid points.
set_folders()
def set_grid(new_grid=grid,
             new_max_depth=max_depth, # These names just need to be different... Because we cannot declare an input name global below...
             subindices=None,
             plot = False, stat = False,
             ):
    """ Obtain the 1/16 degree grid used in medsea_rea, and classify each grid points according to 'max_depth'. 
        A sub-grid is set to the global variables in the module, and also returned by specifying 'subindices':
            ** 9x test grid (the only one with both lat and lon being multiples of 0.25): 
                    subindices=(slice(1,None,4), slice(0,None,4))
            ** 9x grids (16 of them, including the 9x test grid above): 
                    subindices=(slice(i,None,4), slice(j,None,4)) 
               for (i,j) in {0,1,2,3} x {0,1,2,3}
            ** 1x grid that is co-locational with ERA data grid:
                    subindices=(slice(9,None,12), slice(None,None,12))
            ** a mini grid with only 23 sea locations with depth >= 75 that can be used for testing:
                    subindices=(slice(1,None,48), slice(None,None,48))
               NOTE: In medsea_ERA dataset, the latitudes are arranged, exceptionally, in descending order.
 
       Returns three tuples:
            subgrid = (grid_lats, grid_lons, medsea_flags, max_depth)
            rea_indices = (medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth)
            grid_indices = (M, N, sea_mn, sea_m, sea_n) 
    """
    print('Initializing grid...')
    # Declaring global is necessary for modifying them interactively after importing this module.
    global run, grid, max_depth, ASM 
    global grid_lats, grid_lons, medsea_flags, max_depth
    global medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth 
    global M, N, sea_mn, sea_m, sea_n

    # Override the global variables using values passed from function call.
    grid = new_grid
    max_depth = new_max_depth

    # Update the run name
    run = 'ASM{:d}-{:d}m-{:s}'.format(ASM,max_depth,grid)

    # Set up the slices for the subgrid using the name.
    if new_grid in ['9x_{:d}'.format(i) for i in range(16)]:
        k = int(new_grid[3:])
        #print('k=',k)
        if subindices is not None:
            raise Exception("Mistakes? 'subindices' need not be given if a known 'grid' is provided.")
        else:
            i = int(k/4)
            j = k%4
            subindices = (slice(i,None,4), slice(j,None,4))

    if new_grid == '144x':
        subindices = (slice(None,None,None),slice(None,None,None))
        
    if new_grid == '1x':
        subindices = (slice(9,None,12), slice(None,None,12))

    if new_grid == 'mini':
        subindices=(slice(1,None,48), slice(None,None,48))
        
    medsea_rea_lat_ind = subindices[0]
    medsea_rea_lon_ind = subindices[1]
    #print(medsea_rea_lat_ind)
    #print(medsea_rea_lon_ind)

    # Load the rea grid, and a sample set of data for its masks.
    with Dataset(os.path.join(p_sossta_folder, 'medsea_rea/2013/20130101_TEMP_re-fv6.nc'),'r') as ds:
        lat_rea = ds['lat'][:]
        lon_rea = ds['lon'][:]
        temp_rea = ds['votemper'][:]
        depth_rea = ds['depth'][:]

    ndepth = sum(depth_rea<max_depth)+1
 
    # 2 means deeper than max_depth, 1 means less than max_depth, 0 means land
    loc_type = 0 + ~temp_rea[0,0,:].mask + ~temp_rea[0,ndepth,:].mask
        
    # Setting values to global names.
    grid_lats = lat_rea[subindices[0]]
    grid_lons = lon_rea[subindices[1]]
    medsea_flags = loc_type[subindices[0],subindices[1]]
    #print(grid_lats,grid_lons)
    assert grid_lats.shape, grid_lons.shape == medsea_flags.shape

    M, N = medsea_flags.shape
    assert M == grid_lats.size and N == grid_lons.size

    sea_m, sea_n = np.where(medsea_flags==2)
    sea_mn = [(sea_m[i],sea_n[i]) for i in range(sea_m.size)]
    assert len(sea_mn) == sea_m.size

    #print(grid_lats.size,grid_lons.size,grid_lats.min(),grid_lats.max(),grid_lons.min(),grid_lons.max())
    print('Finished setting up a subgrid of shape {!s} x {!s} with {!s} <= latitude <= {!s}, {!s} <= longitude <= {!s}.'.format(\
            grid_lats.size,grid_lons.size,grid_lats.min(),grid_lats.max(),grid_lons.min(),grid_lons.max()))
    if stat:
        # The following are for the current subgrid.
        def print_stat(bl_array):  
            drei = bl_array.size
            zwei = (bl_array == 2).sum()
            ein = (bl_array == 1).sum()
            null = (bl_array == 0).sum()
            assert null+ein+zwei == drei
            print('\t' + 'sea (>{!s}m) locations: {!s}'.format(max_depth,zwei) + ', {:4.2f}% out of {!s}'.format(zwei*100/drei,drei))
            print('\t' + 'sea (<{!s}m) locations: {!s}'.format(max_depth,ein) + ', {:4.2f}% out of {!s}'.format(ein*100/drei,drei))
            print('\t' + 'land locations: {!s}'.format(null) + ', {:4.2f}% out of {!s}'.format(null*100/drei,drei))        

        print('In the current grid:')
        print_stat(medsea_flags)
        print('In the full medsea_rea grid:')
        print_stat(loc_type)

    if plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2,1,figsize=(14,14))
        ax1, ax2 = axes
        ax1.imshow(loc_type,origin='lower',extent=[lon_rea.min(),lon_rea.max(),lat_rea.min(),lat_rea.max()],cmap=cm.Blues)
        ax1.set_title('full medsea_rea grid')
        ax2.imshow(medsea_flags,origin='lower',extent=[lon_rea.min(),lon_rea.max(),lat_rea.min(),lat_rea.max()],cmap=cm.Blues)
        ax2.set_title('current subgrid')
        for ax in axes:
            ax.set_xlabel('longitude')
            ax.set_ylabel('latitude')
          
    # These values have been written directly to global variables as well.
    subgrid = (grid_lats, grid_lons, medsea_flags, max_depth)
    rea_indices = (medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth)
    grid_indices = (M, N, sea_mn, sea_m, sea_n) 
    return subgrid, rea_indices, grid_indices

def get_grid():
    """ 
    A simple getter for the global variables, returning three tuples.

    subgrid = (grid_lats, grid_lons, medsea_flags, max_depth)
    rea_indices = (medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth)
    grid_indices = (M, N, sea_mn, sea_m, sea_n) 
    """

    global grid_lats, grid_lons, medsea_flags, max_depth
    global medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth 
    global M, N, sea_mn, sea_m, sea_n
    return (grid_lats, grid_lons, medsea_flags, max_depth), \
        (medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth), \
        (M, N, sea_mn, sea_m, sea_n)

# Set default global values at the loading of this module.

# For the grid info. Load from the file if it exists.
if not(os.path.isfile(os.path.join(project_folder,'medsea_grid_data.npy'))):
    subgrid_data = set_grid()
    np.save(os.path.join(project_folder,'medsea_grid_data.npy'),subgrid_data)
else:
    subgrid, rea_indices, grid_indices = np.load(os.path.join(project_folder,'medsea_grid_data.npy'))
    grid_lats, grid_lons, medsea_flags, max_depth = subgrid
    medsea_rea_lat_ind, medsea_rea_lon_ind, ndepth = rea_indices
    M, N, sea_mn, sea_m, sea_n = grid_indices

# For the folders.
set_folders()


# GOTM dat files' netCDF reformatted dataset sources.
def data_sources(year=None, month=None, mode='r', dat=['heat','met','tprof','sprof','chlo']):
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

    fn_dict = {'heat' : os.path.join(data_folder,region+'_ERA-INTERIM',region+'_ERA_heat' + suffix),
               'met'  : os.path.join(data_folder,region+'_ERA-INTERIM',region+'_ERA_met' + suffix),
               'tprof': os.path.join(data_folder,region+'_rea',region+'_rea_votemper' + suffix),
               'sprof': os.path.join(data_folder,region+'_rea',region+'_rea_vosaline' + suffix),
               'sst'  : os.path.join(data_folder,region+'_OSTIA',region+'_OSTIA_sst' + suffix)}

    if year is None:
        MODIS_suffix = '*.nc'
    else:
        MODIS_suffix = '_{:d}.nc'

    # Temporary special treatment just to make it work for new MODIS data.
    fn_dict.update(chlo = os.path.join(data_folder,region+'_MODIS','8days',region+'_MODIS_chlor_a_8D' + suffix)) 
    fn_dict.update(iop = os.path.join(data_folder,region+'_MODIS','8days',region+'_MODIS_IOP_8D' + suffix))
    
    assert all([each in fn_dict.keys() for each in dat]) # Check that the function is called correctly.

    for each in fn_dict.keys():
        try:
            ds_dict = {each : NCDataset(fn_dict[each],mode) for each in dat}
        except OSError:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            raise
        except:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            print('Requested list of dat files: {!s}.'.format(dat))
            print('Available built-in data sources: {!s}.'.format(fn_dict))
            raise

    if len(ds_dict.keys()) == 1:
        # If only one dataset requested, return the netcdf dataset unwrapped from the dict. 
        return ds_dict[dat[0]] 
    else:
        # Else return a dictionary of datasets.
        return ds_dict

# Global setting for the core_dat() routines (and possibly the ERA routines as well)



# 2017-05-20 First time doing this, let's be safe.
#assert all(grid_lats == np.arange(30.25,45.75+0.25,0.25))
#assert all(grid_lons == np.arange(-6.0,36.25+0.25,0.25))

# 20170521 Following no longer needed.
## The global variables that are still needed by some code.
#M = grid_lats.size
#N = grid_lons.size
#
#sea_m, sea_n = np.where(medsea_flags == 2)
#sea_mn = [(sea_m[i],sea_n[i]) for i in range(sea_m.size)]

## Helper functions
def timestr(nctime,i):
    " Return a formatted time string from a nc time variable at index i."
    from netCDF4 import datetime, num2date
    try:
        ts = datetime.strftime(num2date(nctime[i],nctime.units),'%Y-%m-%d %H:%M:%S')
    except:
        print("Converting datetime to string failed!")
        print('i,len(nctime),nctime[i],nctime.units')
        print(i,len(nctime),nctime[i],nctime.units)
        raise
    return ts

def change_base(new_base_folder):
    global base_folder, run_folder
    if not(os.path.isdir(new_base_folder)):
        raise IOError('The base folder: ' + new_base_folder + ' is either not accessible or created.')
    base_folder = new_base_folder    
    run_folder=os.path.join(base_folder,run)
    if not(os.path.isdir(run_folder)):
        os.mkdir(run_folder)
       
def print_lat_lon(lat,lon,fmt_str='g'):
    "Helper function for printing (lat,lon) as 10.5N2.1E etc. "

    # if not(isinstance(lat,float)):
    #     raise Exception("`lat` is of type " + str(type(lat)))
    # if not(isinstance(lon,float)):
    #     raise Exception("`lon` is of type " + str(type(lon)))
    lat = float(lat)
    lon = float(lon)
    
    template = '{:' + fmt_str + '}'
    lat_str = template.format(lat) + 'N' if lat>=0 else template.format(-lat) + 'S'
    lon_str = template.format(lon) + 'E' if lon>=0 else template.format(-lon) + 'W'
    #return lat_str + ' ' + lon_str
    return lat_str + lon_str

def chunk(c,k,fun,*args,**kwargs):
    """
    calls fun() for the c-th chunk with k grid points (c-1)*10, (c-1)*10+1, ... (c-1)*10+9. 
    This partitions the 390 grid points on sea into 30 chunks.
    """


def get_m_n(lat,lon):
    "Return the grid index (m,n) given latlong. Assume uniformly-spaced grid."
    spacing = grid_lats[1]-grid_lats[0]
    m = (lat-grid_lats[0])/spacing
    n = (lon-grid_lons[0])/spacing
    if m%1 != 0. or n%1 != 0.:
        print('Warning: given latlong ({!s},{!s}) does not correspond exactly to a point in the {:s} grid.'.format(lat,lon,grid))
        m = int(round(m))
        n = int(round(n))
        print('Returning the indices to the closest grid point instead: ({!s},{!s}).'.format(grid_lats[m],grid_lons[n]))
    return int(m), int(n)

def get_lat_lon(m,n):
    "Return the latlong given the grid index (m,n)."
    return grid_lats[m],grid_lons[n]

def create_dimensions(nc, lat=grid_lats, lon=grid_lons):
    " Declaring dimensions and creating the coordinate variables for each dimension."

    # Dimensions
    nc.createDimension('time') # unlimited
    nc.createDimension('lat', size = len(lat))
    nc.createDimension('lon', size = len(lon))
    
    # Dimension variables.
    nctime = nc.createVariable('time','i4',dimensions=('time',))
    nctime.units = 'hours since ' + str(epoch) # epoch is to be set in global config, and loaded by importing .config
    nclat = nc.createVariable('lat','f4',dimensions=('lat',))
    nclat.units = 'degrees north'
    nclon = nc.createVariable('lon','f4',dimensions=('lon',))
    nclon.units = 'degrees east'
    nclat[:] = lat
    nclon[:] = lon
    #print('Done initializing dimensions.')
    return nctime, nclat, nclon

def create_variable(nc,varname,datatype,dimensions=('time','lat','lon'),zlib=True, fill_value=1e+20):
    " Default settings applied to create a netCDF variable. 'fill_value' of the rea dataset is used here."
    ncvar = nc.createVariable(varname,datatype,dimensions=dimensions,zlib=zlib,fill_value=fill_value)
    #print('Done initializing variables.')
    return ncvar

def prepare_engine():
    " Prepare ipyparallel engines by importing settings and dependencies. "
    from ipyparallel import Client
    rc = Client()
    dv = rc[:]
    lv = rc.load_balanced_view()

    with dv.sync_imports():
        import os, sys
    dv.execute("userhome = os.getenv('HOME')")
    dv.execute('from gotm import *')

    # Push all possibly changed folder locations from local namesapce to each engine. 
    dv.execute("project_folder=" + project_folder)
    dv.push(dict(project_folder = project_folder,
                 data_folder = data_folder,
                 run_folder = run_folder,
                 base_folder = base_folder,
                 p_sossta_folder = p_sossta_folder,
                 ERA_folder = ERA_folder,
                 rea_folder = rea_folder))
    dv.apply(change_base,base_folder)
    dv.execute('os.chdir("{}")'.format(base_folder))
    return rc, lv

## Medsea serial / parallel run toolbox

#def get_local_folder(m=None,n=None,new_run=None,i=None,new_lats=None,new_lons=None,create=False):
def get_local_folder(*args, new_lats=None,new_lons=None,create=False):
    """
    Return the corresponding local folder for given grid point indices (m,n) or linear index i of the 
    index arrays sea_m and sea_n.
 
    Defaults to the run_folder set by the module 'run' global variable, and 'run_folder' generated by
    set_folders(). If not desired, pass a new_run keyword argument to get a different run's local folder.
    
    TODO: It will be convenient that this function be overloaded by argument type and number (single/multiple
    dispatch). We want get_local_folder(i) and get_local_folder(m,n) or get_local_folder((m,n)) to both work. 
    Oh well, that's an overkill. Just use *args and count the number for now, that single dispatch thing is 
    better if we need to distinguish by the type of argument.
    """
    # # Temporary hack, be forgiving if the provided lat, lon are actually indices of our medsea grid.
    # if isinstance(lat_or_m,int):
    #     lat = grid_lats[lat_or_m]
    # else:
    #     lat = lat_or_m
    # if isinstance(lon_or_n,int):
    #     lon = grid_lons[lon_or_n]
    # else:
    #     lon = lat_or_m

    ## Setting default argument values.

    # use module defaults if not provided:
    if new_lats is None:
        new_lats = grid_lats
    if new_lons is None:
        new_lons = grid_lons

    if args is None:
        print('No grid indices given, defaulting to the favourite grid point (lat, lon) = (36.00, 25.50) near buoy 61277...')
        m, n = get_m_n(36.,25.5) 
        print('m =',m,'n =',n)
    elif len(args) == 1:
        i = args[0]
        assert isinstance(i,int), str(i)
        m = sea_m[i]
        n = sea_n[i]
        print('Using linear index of sea locations, i = {!s}, (m,n) = ({!s},{!s})...'.format(i,m,n))
    elif len(args) == 2:
        m,n = args
        if not(int(m) == m and int(n) == n):
            raise Exception('Given (m, n) = ({!s}, {!s}) do not seem to be grid indices.'.format(m,n))
    else:
        print('Position arguments given: ', args)
        raise Exception('Wrong number of positional arguments given.')

    ## Set up local variables
    lat  = new_lats[m]
    lon = new_lons[n]
    latlong = print_lat_lon(lat,lon)

    # 20170524 Need to be reviewed after the new folder structure considerations. 
    #local_folder = os.path.join(base_folder,new_run,latlong)
    local_folder = os.path.join(base_folder,latlong)

    if not(os.path.isdir(local_folder)):
        if create:
            print('The folder {:s} is not found, creating it.'.format(local_folder))
            os.mkdir(local_folder)
        else:
            raise IOError("The local folder {:s} is not found. Have you run local_dat()?".format(local_folder))
    return local_folder

def get_core_folder(year,month,lat,lon):
    """ Create folder structure initially or mid-way (if not yet done). 
        The innermost subfolder, which is called 'core_folder' is 
        where a GOTM run is executed for one grid point per time period.  
        It contains settings (*.inp), input data (*.dat) and output data (*.nc) """
    # Temporary hack, be forgiving if the provided lat, lon are actually indices of our medsea grid.
    if isinstance(lat,int):
        lat = grid_lats[lat]
    if isinstance(lon,int):
        lon = grid_lons[lon]
        
    monthly_folder = os.path.join(run_folder,'{:d}{:02d}'.format(year,month))
    if not(os.path.isdir(monthly_folder)):
        os.mkdir(monthly_folder)
    latlong = print_lat_lon(lat,lon)
    core_folder = os.path.join(monthly_folder,latlong)
    if not(os.path.isdir(core_folder)):
        os.mkdir(core_folder)
    return core_folder 

def write_dat(m,n,dat_fn,nc,outdir):
    " Write dat files for each lat/lon in grid_lats/grid_lons from a given netCDF Dataset or MFDataset."

    from numpy.ma import is_masked
    from netCDF4 import Dataset, MFDataset

    if outdir is None:
        outdir = get_local_folder(m,n)
    if isinstance(nc,Dataset) or isinstance(nc,MFDataset):
        time = nc['time']
    elif os.path.isfile(nc):
        nc = Dataset(nc,'r')
        time = nc['time']
    fn = os.path.join(outdir,dat_fn+'.dat')
    if os.path.isfile(fn) and not overwrite:
        print(fn + " exists, skipping.\n")
        return    

    with open(fn,'w') as f:
        # Recipes for each type of dat file.
        # print(fn) # Print the filename for debug.
        if dat_fn == 'tprof':
            for i in range(len(time)):
                nc_depth = nc['depth']
                f.write(timestr(time,i) + ' {0:d} 2\n'.format(len(nc_depth))) # Always two columns.
                for j in range(len(nc_depth)):
                    line = ('{0:g} {1:g}\n').format(-nc['depth'][j],nc['votemper'][i,j,m,n])
                    f.write(line)
            done = True
        elif dat_fn == 'sprof':
            for i in range(len(time)):
                nc_depth = nc['depth']
                f.write(timestr(time,i) + ' {0:d} 2\n'.format(len(nc_depth))) # Always two columns.
                for j in range(len(nc_depth)):
                    line = ('{0:g} {1:g}\n').format(-nc['depth'][j],nc['vosaline'][i,j,m,n])
                    f.write(line)
            done = True
        elif dat_fn == 'heat':
            col = [None for i in range(4)]
            # Temporary hack #1: repeat the first record if it starts at 03:00:00 so that the 
            # simulation can start at midnight instead. Otherwise, GOTM will generate a
            # plethora of nan values.
            if timestr(time,0)[-8:] == '03:00:00':
                col[0] = timestr(time,0)[:-8] + '00:00:00'
                I_0_obs = nc['swrd'][0,m,n]
                col[1] = I_0_obs

                # Let's not use cloud_factor, some values got accidentally masked.
                if 'swrd_clear_sky' not in nc.variables:
                    col[2] = 1
                else:
                    I_0_calc = nc['swrd_clear_sky'][0,m,n]
                    if I_0_obs < 1 or I_0_calc < 1:
                        col[2] = 0
                    else:
                        col[2] = I_0_obs / I_0_calc
                #col[2] = 1 if 'cloud_factor' not in nc.variables else nc['cloud_factor'][0,m,n]

                col[3] = nc['lwrd'][0,m,n]
                try:
                    line = ('{:s}' + ' {:g}'*3 + '\n').format(*col)
                    f.write(line)
                except Exception:
                    print('col',col)
                    print('m,n',m,n)
                    raise  

            # The following loop is hacked temporarily to accomodate GOTM Fortran code assumption,
            # reinterpreting the timsteamp to mean the beginning of 3-hourly periods.
            for i in range(len(time)-1): #  Last record is not used.
                #col[0] = timestr(time,i)
                col[0] = timestr(time,i)

                I_0_obs = nc['swrd'][i+1,m,n]
                col[1] = I_0_obs
                
                # Let's not use cloud_factor, some values got accidentally masked.
                if 'swrd_clear_sky' not in nc.variables:
                    col[2] = 1
                else:
                    I_0_calc = nc['swrd_clear_sky'][i+1,m,n]
                    if I_0_obs < 1 or I_0_calc < 1:
                        col[2] = 0
                    else:
                        col[2] = I_0_obs / I_0_calc

                #col[2] = 1 if 'cloud_factor' not in nc.variables else nc['cloud_factor'][i+1,m,n]

                col[3] = nc['lwrd'][i+1,m,n]
                try:
                    line = ('{:s}' + ' {:g}'*3 + '\n').format(*col)
                    f.write(line)
                except Exception:
                    print('col',col)
                    print('m,n',m,n)
                    raise
        
            # Temporary hack #3: include the last line being first day next month midnight, whose
            # value should not be used because of the new interpretation of timing. However, for some 
            # reason GOTM halted at the last hour of time. So let's just repeat the value 3 hours earlier
            # at 21:00:00 last day of month.
            assert timestr(time,-1)[-8:] == '00:00:00' # The last record is at midnight.
            col[0] = timestr(time,-1)

            I_0_obs = nc['swrd'][-1,m,n]
            col[1] = I_0_obs
            
            # Let's not use cloud_factor, some values got accidentally masked.
            if 'swrd_clear_sky' not in nc.variables:
                col[2] = 1
            else:
                I_0_calc = nc['swrd_clear_sky'][-1,m,n]
                if I_0_obs < 1 or I_0_calc < 1:
                    col[2] = 0
                else:
                    col[2] = I_0_obs / I_0_calc

            #col[2] = 1 if 'cloud_factor' not in nc.variables else nc['cloud_factor'][-1,m,n]
            col[3] = nc['lwrd'][-1,m,n]
            try:
                line = ('{:s}' + ' {:g}'*3 + '\n').format(*col)
                f.write(line)
            # debug
            except Exception:
                print('col',col)
                print('m,n',m,n)
                raise

            # Complicated write. Finally done.
            done = True

        elif dat_fn == 'met':
            col = [None for i in range(9)]
            # Temporary hack #1 (with the one for 'heat'). Just repeat the first value but use it for midnight.
            if timestr(time,0)[-8:] == '03:00:00':
                col[0] = timestr(time,0)[:-8] + '00:00:00'
                col[1] = nc['u10m'][0,m,n]
                col[2] = nc['v10m'][0,m,n]
                col[3] = nc['sp'][0,m,n]/100 # surface pressure, convert from Pa to hPa
                col[4] = nc['t2m'][0,m,n] - 273.15 # a0r temperature at 2m, convert to Celsius
                col[5] = nc['q2m'][0,m,n] # specific humidity at 2m
                col[6] = 0 # "cloud" value?
                col[7] = nc['precip'][0,m,n]
                col[8] = nc['snow'][0,m,n]
                line = ('{:s}'+' {:g}'*8 + '\n').format(*col)
                f.write(line)
                done = True
            for i in range(len(time)):
                col[0] = timestr(time,i)
                col[1] = nc['u10m'][i,m,n]
                col[2] = nc['v10m'][i,m,n]
                col[3] = nc['sp'][i,m,n]/100 # surface pressure, convert from Pa to hPa
                col[4] = nc['t2m'][i,m,n] - 273.15 # air temperature at 2m, convert to Celsius
                col[5] = nc['q2m'][i,m,n] # specific humidity at 2m
                col[6] = 0 # "cloud" value?
                col[7] = nc['precip'][i,m,n]
                col[8] = nc['snow'][i,m,n]
                line = ('{:s}'+' {:g}'*8 + '\n').format(*col)
                f.write(line)
            done = True
        elif dat_fn == 'sst':
            for i in range(len(time)):
                line = '{:s} {:g}\n'.format(timestr(time,i),nc['analysed_sst'][i,m,n])
                f.write(line)
            done = True
        elif dat_fn == 'chlo':
            count = 0
            for i in range(len(time)):
                if is_masked(nc['chlor_a'][i,m,n]):
                    #print('i,m,n')
                    #print(i,m,n)
                    #print('time[i]')
                    #print(timestr(time,i))
                    count +=1
                    continue
                line = '{:s} {:g}\n'.format(timestr(time,i),nc['chlor_a'][i,m,n])
                f.write(line)
            if count > 4:
                print('WARNING: {:d} consecutive nan values.'.format(count))
                #raise Exception("3 consecutive nan values.")
                done = False
            else:
                done = True
                
        elif dat_fn == 'iop':
            count = 0
            for i in range(len(time)):
                if is_masked(nc['a_488_giop'][i,m,n]) or is_masked(nc['bb_488_giop'][i,m,n]):
                    #print('i,m,n')
                    #print(i,m,n)
                    #print('time[i]')
                    #print(timestr(time,i))
                    count +=1
                    continue
                line = '{:s} {:g} {:g}\n'.format(timestr(time,i),nc['a_488_giop'][i,m,n],nc['bb_488_giop'][i,m,n])
                f.write(line)
            if count > 4:
                print('WARNING: {:d} consecutive nan values.'.format(count))
                #raise Exception("3 consecutive nan values.")
                done = False
            else:
                done = True
        else:
            raise Exception("Requested {}.dat has no recipes defined in core_dat()".format(dat_fn))

    if done:
        print('Done writing {}.\n'.format(fn))
    else:
        if os.path.isfile(fn):
            print('Writing failed. Deleting {:s}... '.format(fn))
            os.remove(fn)

def local_dat(mm,nn,dat=['heat','met','tprof','sprof','chlo','iop']):
    """
    Generate *.dat files from all available data. See core_dat() for other explanations.
    mm, nn can be a sequence of m,n's
    """

    from netCDF4 import MFDataset
    from tempfile import mkdtemp
    from shutil import copyfile
    from glob import glob
    
    if isinstance(dat,str):
        dat = [dat]

# 20170621, we no longer create separate folders for different runs, but share the same set of subfolders named
# by print_lat_lon(), to avoid creating too many files and draining disk quota too fast.

#    run_folder = os.path.join(base_folder,run)
    run_folder = base_folder

#    if not(os.path.isdir(run_folder)):
#        os.mkdir(run_folder)

    # temp_folder = mkdtemp(prefix=temp_base_folder)
    # ERA_files = glob(os.path.join(data_folder,'medsea_ERA-INTERIM','*.nc')) 
    # rea_files = glob(os.path.join(data_folder,'medsea_rea','*.nc'))
    # ERA_tempfiles = [copyfile(fn,os.path.join(temp_folder,os.path.basename(fn))) for fn in ERA_files]
    # rea_tempfiles = [copyfile(fn,os.path.join(temp_folder,os.path.basename(fn))) for fn in rea_files]

    print("Looking for data sources from " + data_folder + "...")
    nc_dict = data_sources(dat=dat)
    if not(isinstance(nc_dict, dict)):
        nc_dict = {dat[0]: nc_dict} 

    for dat_fn, nc in nc_dict.items():
        ## Assume m, n are iterable and of the same length:
        if isinstance(mm,int) and isinstance(nn,int):
            mm = [mm]
            nn = [nn]
        else:
            assert len(m) == len(n)
        try:
            for m, n in zip(mm, nn):
                lat = grid_lats[m]
                lon = grid_lons[n]

                # Create the local folder if necessary.
                local_folder = os.path.join(run_folder,print_lat_lon(lat,lon))
                if not(os.path.isdir(local_folder)):
                    os.mkdir(local_folder)
                    
                # Actually write.
                write_dat(m,n,dat_fn,nc,local_folder)
        except Exception:
                print('mm',mm)
                print('nn',nn)
        
        ## When m,n are integers, exception occurs at zip() instead, and it not a TypeError.
        # except TypeError as te:
        #     if isinstance(mm,int) and isinstance(nn,int):
        #         m = mm
        #         n = nn
        #         # So the provided sequences are just a single pair of indices. Just write.
        #         write_dat(m,n,dat_fn,nc,local_folder)
        #    else:

#        write_dat(m,n,dat_fn,nc,local_folder)
        
    for nc in nc_dict.values():
        nc.close()

def core_dat(year,month,m,n,**nc_dict):
    """
    Generate *.dat files for each core folder, by months. 
    The critical keyword argument 'nc_dict' provides dat filename to nc Dataset
    handle correspondance, e.g. {'heat': heat_nc} where 
                   heat_nc = Dataset('ERA_heat_yyyymm.nc','r') 
    has been run before calling this function. The file 'heat.dat' would be generated by
    using information from 'ERA_heat_yyyymm.nc', and recipe defined in an inner function.
    """

    # Get the location of the core_folder.
    lat = grid_lats[m]
    lon = grid_lons[n]
    # latlong = print_lat_lon(lat,lon)
    core_folder = get_core_folder(year,month,lat,lon) 
    
    for dat_fn, nc in nc_dict.items():
        # print(dat_fn)
        write_dat(m,n,dat_fn,nc,core_folder)
    return

def prepare_run(start,stop,run_folder,out_dir='.',out_fn='results',m=None,n=None,lat=None,lon=None, **gotm_user_args):
    "Transfer config files and GOTM executable to the folder in which GOTM will be run."
    import shutil
    # Determine the grid point location.
    if (m is None) and (n is None):
        (m,n) = get_m_n(lat,lon)
    if (lat is None) and (lon is None):
        (lat,lon) = get_lat_lon(m,n)
    assert not(m is None or n is None or lat is None or lon is None), 'Either a (m,n) or (lat,lon) should be specified.'

    run_name = 'medsea_GOTM, {:s} grid, #(m,n)=({:d},{:d})'.format(grid,m,n)    

    # Set up GOTM arguments.
    gotm_args = dict(name = run_name,
                     start = str(start), stop = str(stop), 
                     latitude = float(lat), longitude = float(lon), 
                     out_dir = str(out_dir), out_fn = out_fn)
    gotm_args.update(**gotm_user_args) 

    # Get the externally specified sea level widths if method is 2.
    #print(gotm_user_args.keys())
    #print(gotm_user_args['grid_method'])
    #print(gotm_user_args['grid_file'])
    if 'grid_method' in gotm_user_args.keys() and gotm_user_args['grid_method'] == 2:
        assert 'grid_file' in gotm_user_args.keys()
        grid_file = gotm_user_args['grid_file']

    # Copy over GOTM executable and config files, the latter to be updated when gotm() is called.
    if not(os.path.exists(GOTM_executable)):
        os.symlink(GOTM_executable,os.path.join(folder,'gotm'))    
    for each in GOTM_nml_list:
        # All config files are overwritten every time GOTM is run.
        src = os.path.join(GOTM_nml_path,each)
        dst = os.path.join(run_folder,each)
        # print('Copying {:s} to {:s}'.format(src,dst))
        shutil.copyfile(src,dst)       
    
    # Temporary hack. Also copy the grid data file.
    src = os.path.join(GOTM_nml_path,grid_file)
    dst = os.path.join(run_folder,grid_file)
    # print('Copying {:s} to {:s}'.format(src,dst))
    shutil.copyfile(src,dst)

    updatecfg(path=run_folder, **gotm_args)
    return gotm_args

def local_run(year,month,m,n,run,start=None,stop=None,create=False,verbose=False,**gotm_user_args):
    """ 
    
    Generate GOTM results for the (m,n)-th grid point at the specified month. Only *.dat files are expected to be 
    in the local folder. The config file, GOTM run time will be generated or copied over. 

    Temporary hack: pass None to month to run for the whole year.

    """
    from datetime import datetime
    from netCDF4 import Dataset, num2date
    import shutil
    
    lat, lon = get_lat_lon(m,n)
#    local_folder = get_local_folder(lat,lon,run)
    local_folder = get_local_folder(m,n,create=create)

    # Argument handling
    if year is None and month is None: # Run for a specific period.
        assert start is not None
        assert stop is not None
        def check(string):
            if isinstance(string,str):
                assert len(string) == 4+3+3+ 3+3+3 # 2013-01-01 00:00:00
            else:
                assert isinstance(string, datetime)
        check(start)
        check(stop)
        startdayofyear = (start-datetime(start.year-1,12,31)).days # E.g. start is 2013-01-07, then it is 7.
        stopdayofyear = (stop-datetime(stop.year-1,12,31)).days # E.g. stop is 2013-12-31, then it is 365, if stop is 2014-01-01, then it is 1.
        suffix = '_' + \
                 '{:04d}{:04d}{:04d{:04d}'.format(start.year,startdayofyear,stop.year,stopdayofyear) + \
                 '_' + run 
        
    elif month is None: # Run for a year.
        start = datetime(year,1,1)
        stop = datetime(year+1,1,1)
        suffix = '_' + run + '_' + '{:04d}'.format(start.year)
        
    else: # Run for a month
        assert year is not None
        start = datetime(year,month,1);
        stop = datetime(year,month+1,1) if month < 12 else datetime(year+1,1,1)
        suffix = '_' + run + '_' + '{:04d}{:02d}'.format(start.year, start.month)

    # Should GOTM write to the local folder or a cached folder?
    out_dir = local_folder if cache_folder is None else cache_folder
    out_fn = 'results' + suffix

    if run in run_profiles.keys():
        # Should subclass an Exception to tell people what happened.
        if gotm_user_args != {}:
            print(('{:s} = {!s}' * len(gotm_user_args)).format(*(gotm_user_args.items())))
            raise Exception("A recorded run profile {:s} is specified, rejecting all user arguments for GOTM.\n")
        print('Using pre-defined profile: {:s}...'.format(run))
        gotm_args = prepare_run(start,stop,local_folder,lat=lat,lon=lon,
                                out_fn=out_fn,
                                out_dir=out_dir,
                                daily_stat_fn='daily_stat'+suffix+'.dat',
                                **run_profiles[run])
    else:
        print('Running without a pre-defined profile and creating new folder structures for {:s}...'.format(run))
        gotm_args = prepare_run(start,stop,local_folder,lat=lat,lon=lon,
                                out_dir=out_dir,
                                out_fn=out_fn,
                                daily_stat_fn='daily_stat'+suffix+'.dat',
                                **gotm_user_args)

    os.chdir(local_folder)
    stat = dict()
    try:
        tic()
        logfn = 'GOTM_' + print_ctime(sep='_') + '.log'
        print('GOTM running... ')
        print('  Working from: {:s}...'.format(local_folder))
        gotm(verbose=verbose,logfn=logfn)
        print('  Results written to {:s}.nc...'.format(os.path.join(out_dir,out_fn)))

        if cache_folder is not None:
            src = os.path.join(out_dir,out_fn+'.nc')
            dst = os.path.join(local_folder,out_fn+'.nc')
            try:
                print('Moving {:s} to {:s}...'.format(src,dst))
                shutil.move(src,dst)
            except:
                raise
            
        stat['elapsed'] = toc()
        # statfn = 'stat_{:d}{:02d}.dat'.format(year,month)
        statfn = 'run_stat_' + print_ctime(sep='_') + '.log'
        with open(statfn,'w') as f:
            print('Writing diagnostic statistics to {0}...\n'.format(statfn))
            f.write('--------------------------------------------------------------\n')
            f.write('Run parameters:\n')
            for key, val in gotm_args.items():
                f.write('    {:s} = {!s}\n'.format(key,val))
            f.write('--------------------------------------------------------------\n')
            f.write('Run statistics:\n')
            f.write('Elapsed: {:.2f} seconds.\n'.format(stat['elapsed']))
            f.write('--------------------------------------------------------------\n')
            f.write('Data statisitics:\n')
            with Dataset(os.path.join(local_folder,gotm_args['out_fn']+'.nc'),'r') as ds:
                sst = ds['sst'][:]
                time = ds['time']
                stat.update(sst_mean = sst.mean(),
                            sst_max = sst.max(),
                            sst_time_max = num2date(time[sst.argmax()],time.units),
                            sst_min = sst.min(),
                            sst_time_min = num2date(time[sst.argmin()],time.units))
                f.write('   SST:\n')
                f.write('      Mean: {sst_mean:.4g}\n'.format(**stat))
                f.write('      Max: {sst_max:.4g} at {sst_time_max!s}\n'.format(**stat))
                f.write('      Min: {sst_min:.4g} at {sst_time_min!s}\n'.format(**stat))
            f.write('--------------------------------------------------------------\n')
    except:
        raise
    return stat

def core_run(year,month,m=None,n=None,lat=None,lon=None,verbose=False,**gotm_user_args):
    """ Generate GOTM results for the (m,n)-th grid point for a specified month in the MxN (lat,long) medsea grid. 
    
        All necessary files (*.inp, *.dat) are assumed to be present in `core_folder` (see get_core_folder(year,month,m,n)), 
        except the GOTM executable. The program changes directory into the core folder to run and the log and results are 
        both saved in core_folder. """
    from datetime import datetime
    import os, shutil

    ## Setup GOTM arguments for this run.
    start = datetime(year,month,1,0,0)
    stop = datetime(year,month+1,1,0,0) if month < 12 else datetime(year+1,1,1,0,0)
    if not m is None:
        lat = grid_lats[m]
    if not n is None:
        lon = grid_lons[n]
    latlong = print_lat_lon(lat,lon)
    run_name = 'medsea_GOTM, #(m,n)=({1:d},{2:d})'.format(latlong,m,n)
    core_folder = get_core_folder(year,month,lat,lon)    
    
    # NOTE 1: The default values for out_fn, t_prof_file, s_prof_file, heatflux_file, meteo_file, extinct_file, sst_file
    # are set in the template namelist files in the base_folder, but not here. Same for the values of heights of measurements
    # wind_h, rh_h, airt_h adapted for use of our medsea dataset.
    # NOTE 2: Explicit cast types to avoid ValueError in f90nml, which only supports the fundmental data types.
    gotm_args = dict(name = run_name,
                     start = str(start), stop = str(stop), 
                     latitude = float(lat), longitude = float(lon), out_dir = str(core_folder))
    
    gotm_args.update(**gotm_user_args) 
    #print(gotm_args) # Debug
    
    ## Prepare the core folder.
    # Symlink the executable.
    if not(os.path.exists(GOTM_executable)):
        os.symlink(GOTM_executable,os.path.join(core_folder,'gotm'))    
    for each in GOTM_nml_list:
        # All config files are overwritten every time GOTM is run.
        shutil.copyfile(os.path.join(GOTM_nml_path,each),os.path.join(core_folder,each))
    # NOTE: The actual updates of the namelists are currently done in the gotm() call.

    # Actual GOTM run.
    os.chdir(core_folder)   
    try:
        print('GOTM run: ' + run_name + '...')
        logfn = 'GOTM_' + print_ctime(sep='_') + '.log'
        gotm(verbose=verbose, logfn=logfn, run_folder = core_folder, varsout = {}, **gotm_args)
    except:
        os.chdir(base_folder)
        raise # Maybe we should define some sort of exception if GOTM fails.
        
    os.chdir(base_folder)

## Deprecated
def combine_run(year, month, run,
                varsout = None, # The subset of variables to output, defaults to all available variables
                format = 'NETCDF3_CLASSIC', # Do not store in HDF5 format unless we know our collaborators use tools that are compatible.
                cleanup = False): # If number or sizes of files generated is a concern... maybe True here.
    " Combine GOTM results nc files from each grid point into a single monthly nc file. "
    from netCDF4 import Dataset
    from numpy.ma import masked_all, zeros
    
    if month is None:
        start = datetime(year,1,1)
        stop = datetime(year+1,1,1)
        print('Combining medsea GOTM results for {:s}-{:d}...'.format(year,run))
        outfn = 'medsea_GOTM_{:s}-{:d}.nc'.format(run,year)

    else:
        start = datetime(year,month,1);
        stop = datetime(year,month+1,1) if month < 12 else datetime(year+1,1,1)
        print('Combining medsea GOTM results for {:s}-{:d}-{:02d}...'.format(run,year,month))
        # Try a monthly file first.
        outfn = 'medsea_GOTM_{:s}-{:d}{:02d}.nc'.format(run,year,month)
        if not(os.path.isfile(os.path.join(get_local_folder(sea_m[0],sea_n[0]),outfn))):
            # Assume a yearly run was done instead.
            outfn = 'medsea_GOTM_{:s}-{:d}.nc'.format(run,year)
            assert os.path.isfile(os.path.join(get_local_folder(sea_m[0],sea_n[0]),outfn)), print(outfn) 

    print('Writing dimensions and metadata...')
    elapsed = 0
    tic()

    with Dataset(outfn,'w',format=format) as nc:
        # Default dimensions for medsea
        nctime, nclat, nclon = create_dimensions(nc, *grid)
        
        # Having precomputed the sea locations, and use the first point (30.75N,18.75E) in our 21 x 57 medsea grid. Still works
        # for finer grids as long as the sea_locaions global variable is updated.
        fn = os.path.join(base_folder,print_lat_lon(*sea_locations[0]),'results-{0:d}{1:02d}.nc'.format(year,month))

        # Transfer units, dimensions and create the nc variables.
        with Dataset(fn,'r') as first:
            # Is there any need for these? We did put units when calling create_dimesions()
            #nclat.units = first['lat'].units 
            #nclon.units = first['lon'].units
            
            # Create the depth dimension, save the depth levels by reversing order and sign of the GOTM z variable.
            nz = len(first['z'])
            nc.createDimension('depth', size = nz)
            ncdepth = nc.createVariable('depth','f4',dimensions=('depth',))
            ncdepth.units = first['z'].units
            ncdepth[:] = -first['z'][::-1]

            nctime[:] = first['time'][:] 
            nctime.units = first['time'].units # Small danger, here we are overwriting the units set in create_dimensions()

            # Test outputs
            # for var in ds.variables:
            #     if len(ds[var].dimensions) > 1:
            #         print(var,ds[var].units,ds[var].dimensions)

            # DEPRECATED. Create adaptively using the first sea location results nc file instead...
            # Create nc variables for each GOTM output variable specified.
            # ncvar3d = {name: create_variable(nc,name,'f8', dimensions=('time','lat','lon')) for name in var3dnames}
            # ncvar4d = {name: create_variable(nc,name,'f8',dimensions=('time','depth','lat','lon')) for name in var4dnames}
            # for name in var3dnames:
            #     ncvar3d[name].units = first[name].units
            # for name in var4dnames:
            #     ncvar4d[name].units = first[name].units
            
            # New verison. 2017-04-15
            var3d_nc = dict()
            var4d_nc = dict()
            for var in first.variables:
                var_dim = first[var].dimensions
                if var_dim == ('time','lat','lon'): # Make very sure we're talking about the same things.
                    var3d_nc[var] = create_variable(nc,var,'f8', dimensions=('time','lat','lon'))
                    var3d_nc[var].units = first[var].units
                if var_dim == ('time','z','lat','lon'): # For this we need to replace z by depth.
                    var4d_nc[var] = create_variable(nc,var,'f8',dimensions=('time','depth','lat','lon'))
                    var4d_nc[var].units = first[var].units
        elapsed += toc()
        
        print('Begin reading data into memory...')
        tic()

        # Initialize temp arrays to store data.
        var3d_data = dict()
        var4d_data = dict()
        if month == 12:
            num_hr = (datetime(year+1,1,1)-datetime(year,month,1)).days*24;
        else:
            num_hr = (datetime(year,month+1,1)-datetime(year,month,1)).days*24;
        for var in var3d_nc.keys():
            var3d_data[var] = masked_array(zeros((num_hr,M,N)),mask=True)
        for var in var4d_nc.keys():
            var4d_data[var] = masked_array(zeros((num_hr,nz,M,N)),mask=True) # Danger, here the dimensions depends on a preselected grid.
            
        # Now proceed to read from each GOTM result nc file.
        for m, n in indices:
            fn = os.path.join(base_folder,run,print_lat_lon(*get_lat_lon(m,n)),'results-{0:d}{1:02d}.nc'.format(year,month))
            if not(os.path.isfile(fn)):
                print(fn, ' not found. Skipping this grid point...')
                continue
            with Dataset(fn,'r') as each:
                for var in var3d_nc.keys():
                    # print(each[var][:,0,0].shape)
                    # print(var3d_data[var][:,m,n].shape)
                    var3d_data[var][:,m,n] = each[var][:,0,0]
                for var in var4d_nc.keys():
                    # Make sure the depth axis is reversed.
                    var4d_data[var][:,::-1,m,n] = each[var][:,:,0,0] 
            if cleanup:
                os.remove(fn)        
        elapsed += toc()

        print('Begin writing to {}'.format(outfn))
        tic()
        for var in var3d_nc.keys():
            var3d_nc[var][:] = var3d_data[var]
        for var in var4d_nc.keys():
            var4d_nc[var][:] = var4d_data[var]
        elapsed += toc()

        print('Finished combining GOTM results after {0:.0f}s'.format(elapsed))
    return outfn

## This function has not been updated after folder structure change (the 'run' string does not name a folder but a substring in output filenames)

# def combine_stat(run,year,month,format='NETCDF3_CLASSIC', indices=sea_mn):

#     from numpy import loadtxt
#     from numpy import ma

#     print('Combining GOTM daily statistics for {:d}-{:02d}...'.format(year,month))
#     outfn = os.path.join(base_folder,run,'medsea_GOTM_daily_stat_{:d}{:02d}.nc'.format(year,month))

#     print('Writing dimensions and metadata...')
#     elapsed = 0
#     tic()

#     with Dataset(outfn,'w',format=format) as nc:

#         # Dimensions
#         nc.createDimension('day_of_year')
#         nc.createDimension('lat', size = len(grid_lats))
#         nc.createDimension('lon', size = len(grid_lons))
        
#         # Dimension variables.
#         ncday = nc.createVariable('day_of_year','i4',dimensions=('day_of_year',))
#         ncday.units = "number of days since last new year's eve"
#         nclat = nc.createVariable('lat','f4',dimensions=('lat',))
#         nclat.units = 'degrees north'
#         nclon = nc.createVariable('lon','f4',dimensions=('lon',))
#         nclon.units = 'degrees east'
#         nclat[:] = grid_lats
#         nclon[:] = grid_lons
#         print('Done initializing dimensions.')
        
#         nc_assim_time = nc.createVariable('assim_time','i4',dimensions=('day_of_year','lat','lon'))
#         nc_assim_time.units = 'number of seconds since midnight in local time'
#         nc_SST_max = nc.createVariable('SST_max','f4',dimensions=('day_of_year','lat','lon'))
#         nc_SST_max.units = 'degree Celsius'
#         nc_SST_min_day = nc.createVariable('SST_min_day','f4',dimensions=('day_of_year','lat','lon'))
#         nc_SST_min_day.units = 'degree Celsius'
#         nc_SST_min_night = nc.createVariable('SST_min_night','f4',dimensions=('day_of_year','lat','lon'))
#         nc_SST_min_night.units = 'degree Celsius'
#         nc_SST_max_time = nc.createVariable('SST_max_time','i4',dimensions=('day_of_year','lat','lon'))
#         nc_SST_max_time.units = 'number of seconds since midnight in local time'
#         nc_SST_min_day_time = nc.createVariable('SST_min_day_time','i4',dimensions=('day_of_year','lat','lon'))
#         nc_SST_min_day_time.units = 'number of seconds since midnight in local time'
#         nc_SST_min_night_time = nc.createVariable('SST_min_night_time','i4',dimensions=('day_of_year','lat','lon'))
#         nc_SST_min_night_time.units = 'number of seconds since midnight in local time'
#         print('Done creating variables')

#         for m,n in indices:
#             stat_fn = os.path.join(get_local_folder(m,n),'daily_stat-{:d}{:02d}.dat'.format(year,month))
#             assert os.path.isfile(stat_fn)
#             try:
#                 tmp = loadtxt(stat_fn)
#             except ValueError:
#                 print('Problematic location: ' + print_lat_lon(*get_lat_lon(m,n)))
#                 raise
            
#             ncday[:] = tmp[:,0]
#             nc_assim_time[:,m,n] = ma.masked_equal(tmp[:,1],0)
#             nc_SST_min_day_time[:,m,n] = ma.masked_equal(tmp[:,2],0)
#             nc_SST_min_day[:,m,n] = ma.masked_equal(tmp[:,3],99.)
#             nc_SST_max_time[:,m,n] = ma.masked_equal(tmp[:,4],0)
#             nc_SST_max[:,m,n] = ma.masked_equal(tmp[:,5],-99.)
#             nc_SST_min_night_time[:,m,n] = ma.masked_equal(tmp[:,6],0)
#             nc_SST_min_night[:,m,n] = ma.masked_equal(tmp[:,7],99.)
          
#         print("Done creating " + outfn)
#         elapsed += toc()        

## Medsea results visualization toolbox

# Import matplotlib colormap for assigning default values to functions.
from matplotlib import cm
from numpy.ma import mean

def medsea_heatmap(data, # The only necessary argument: a netCDF Dataset to be opened by user.
                   varnames=['sst'],fun = lambda varargin: mean(varargin[0][:],axis=0),
                   ax=None, draw_colorbar=True, vlim = None,cmap=cm.coolwarm):
    # Expect a netCDF4 Data
    from netCDF4 import Dataset
    from matplotlib.pyplot import subplots
    if not(isinstance(data,Dataset)):
        raise Exception('First argument must be a netCDF4 Dataset.')
    
    # Load essential data
    lat = data['lat'][:]
    lon = data['lon'][:]
    #epoch = data['time'].units
    #time = data['time'][:]
    #z = data['z'][:]

    # Evaluate data to be plotted.
    varargin = [data[each] for each in varnames]
    val = fun(varargin)
    
    # Create figures / axis if not provided.
    if ax is None:
        fig, ax = subplots(figsize=(11,8.5)) # Full letter landscape size.
    fig = ax.get_figure()
    
    # Generate heatmap
    if vlim is None:
        im = ax.imshow(val, extent=(lon.min(), lon.max(), lat.min(), lat.max()),
                        interpolation='nearest', origin='lower', cmap=cmap)
    else:
        im = ax.imshow(val, extent=(lon.min(), lon.max(), lat.min(), lat.max()),
                        vmin = vlim[0], vmax = vlim[1],
                        interpolation='nearest', origin='lower', cmap=cmap)
        
    if draw_colorbar:
        cbar = fig.colorbar(im)
    
    ax.grid('on')
    ax.axis('tight')
    
    return fig, ax, im # Return them for further annotations.

def medsea_plot_mean(varname='sst', ax=None, cmap = cm.coolwarm,
                     fig_fn_prefix = 'medsea_ASM0_avg_sst', show_fig = True,
                     year=2014, month=1,
                     nc_folder='/home/simon/gotm-dst/medsea_GOTM/results',nc_fn_prefix='no_assim/p_sossta/medsea_GOTM'):

    from numpy.ma import mean
    data = medsea_data(year=year,month=month,nc_folder=nc_folder,nc_fn_prefix=nc_fn_prefix)
    fig, ax = medsea_heatmap(data, ax=ax, cmap=cmap,
                             varnames=[varname],fun=lambda varargin: mean(varargin[0][:],axis=0))
    # Annotations
    ax.set_title('mean {} for {:d}-{:02d}'.format(varname,year,month))

    # Save the figure
    fig_fn = '{}_{}_mean_{:d}{:02d}.png'.format(fig_fn_prefix,varname,year,month)
    print('Saving {}...'.format(fig_fn))
    fig.savefig(fig_fn)                   
    
    if not(show_fig):
        close(fig)
    
    # Remember to close the file.
    data.close()

def sst_range_ASM0(time,lon,sst):
    nday = int(time[-1]/86400)
    sst_max = ones((nday,M,N))*NaN
    sst_min = ones((nday,M,N))*NaN
    sst_range = ones((nday,M,N))*NaN
    for day in range(nday):
        for k in range(M*N):
            m = int(k/N)
            n = mod(k,N)
            offset = int(ceil(lon[n]/15))

            h1 = 11 - offset + day*24
            h2 = 21 - offset + day*24
            sst_max[day,m,n] = max(sst[h1:h2,m,n])
            
            
            h3 = 3 - offset + day*24 # We increased the starting hour from 1 (HCMR ppt) to 3 to avoid going to the previous day.
            h4 = 12 - offset + day*24
            #sst_min[day,m,n] = min(sst[h1:h2,m,n]) # OH NO! TYPO! The min was taken from the same range as in taking max!!! That means we're measuring the amount of cooling, but it's not all the way to the coolest point yet...
            sst_min[day,m,n] = min(sst[h3:h4,m,n])
    sst_range = sst_max - sst_min
    return sst_range
    
def sst_range_ASM2(time,lon,sst,swr):
    import numpy
    nday = int(time[-1]/86400)
    sst_max = ones((nday,M,N))*NaN
    sst_min = ones((nday,M,N))*NaN
    sst_range = ones((nday,M,N))*NaN
    for day in range(nday):
        for k in range(M*N):
            m = int(k/N)
            n = mod(k,N)
            offset = int(ceil(lon[n]/15))

            h1 = 11 - offset + day*24
            h2 = 21 - offset + day*24
            sst_max[day,m,n] = max(sst[h1:h2,m,n])
            
            swr_hourly_today = swr[day*24:(day+1)*24,m,n] # This could be a masked array! So find could return empty array.
            if size(find(swr_hourly_today)) == 0:
                swr_first_nonzero = -1000 # So that max(3-offset,swr_first_nonzero)=3-offset
            else:
                swr_first_nonzero = find(swr_hourly_today)[0]
            # Without taking the max of swr_first_nonzero, it could cover a jagged portion of sunrise assim.
            h3 = max([3-offset,swr_first_nonzero]) + day*24 
            h4 = 12 - offset + day*24
            sst_min[day,m,n] = min(sst[h3:h4,m,n])
    sst_range = sst_max - sst_min
    return sst_range

def medsea_monthly_mean_heatmaps(assim = 0, show_fig = False, save_fig = True, 
                                 months = [(2013,each) for each in range(1,13)] + [(2014,each) for each in range(1,13)],
                                 fig_folder = 'fig/monthly_means'):
    """ assim = 0 (none), 1 (local midnight), 2 (swr onset after midnight, i.e. sunrise) """
    
    # Choose the right subfolder of nc files by the assim argument.
    if assim == 0:
        nc_folder = 'results/p_sossta/no_assim'
    elif assim == 1:
        nc_folder = 'results/p_sossta/assim_midnight'
    elif assim == 2:
        nc_folder = 'results/p_sossta/assim_sunrise'
    else:
        raise Exception('assim = 0/1/2')
        
    for year,month in months:
        data = medsea_data(year=year,month=month,nc_folder=nc_folder)

        # A template plotter for monthly means.
        def do_heatmap(fig_title_prefix, fig_filename_prefix, **kwargs):
            fig, ax = medsea_heatmap(data,**kwargs)
            ax.set_title(fig_title_prefix + ' for {:d}-{:02d}'.format(year,month))
            if save_fig:
                fig_fn = os.path.join(base_folder,fig_folder,fig_filename_prefix+'_ASM{:d}_{:d}{:02d}.png'.format(assim,year,month))
                print('Saving ' + format(fig_fn) + '...')
                fig.savefig(fig_fn)
            if not(show_fig):
                close(fig) 
                
        # Mean values of "sea surface variables", i.e. those with dimensions time, lat, lon but not depth.
        for varname in data.variables:
            if data[varname].dimensions == ('time','lat','lon'):
                do_heatmap('mean ' + varname,
                           'medsea_mean_' + varname,
                           varnames=[varname],
                           fun=lambda varargin: mean(varargin[0][:],axis=0),
                           cmap=cm.coolwarm)
        # Mean value of cooling effect
        do_heatmap('mean skin cooling effect',
                   'medsea_mean_cooling',
                   varnames=['sst','skint'],
                   fun=lambda args: mean(args[0][:]-args[1][:],axis=0),
                   cmap=cm.cool)
        
        # Mean diurnal warming amount
        if assim == 0:
            varnames = ['time','lon','sst']
            fun = lambda args: mean(sst_range_ASM0(args[0][:],args[1][:],args[2][:]),axis=0)
        elif assim == 2:
            varnames = ['time','lon','sst','swr']
            fun = lambda args: mean(sst_range_ASM2(args[0][:],args[1][:],args[2][:],args[3][:]),axis=0)
        else:
            raise Exception('Assim = 1 not implemented.')
            
        do_heatmap('mean diural warming amount in SST',
                   'medsea_mean_SST_range',
                   varnames=varnames, fun=fun,
                   cmap=cm.hot)
        data.close()
        
def buoy_comparisons(ax=None,year=2014,month=4,station='61277',plot_GOTM_SST=False):
    import os
    home = os.getcwd()
    os.chdir(base_folder)
    
    if ax is None:
        fig, ax = subplots()

    # Buoy
    if station == '61280' or station == '61281' or station == '61430':
        buoy_fn = 'buoys/{2}/IR_{0:d}{1:02d}_TS_MO_{2}.nc'.format(year,month,station)
    elif station == '61277' or station == '68422' or station == 'SARON':
        buoy_fn = 'buoys/{2}/MO_{0:d}{1:02d}_TS_MO_{2}.nc'.format(year,month,station)
    if station == '61277' or station == '68422':
        ind = 1
    else:
        ind = 0

    with MFDataset('buoys/{0:}/*{0:}.nc'.format(station),'r') as buoy_data:
        buoy_lat = (float(buoy_data.geospatial_lat_max) + float(buoy_data.geospatial_lat_min))/2
        buoy_lon = (float(buoy_data.geospatial_lon_max) + float(buoy_data.geospatial_lon_min))/2
        depth_buoy = buoy_data['DEPH'][0,ind]
        
    data = list()
    if os.path.isfile(buoy_fn):
        with Dataset(buoy_fn,'r') as buoy_data:
            print('Buoy {0} at {1} found for {2:d}-{3:02d}'.format(station,print_lat_lon(buoy_lat,buoy_lon),year,month))
            temp_buoy = buoy_data['TEMP'][:,ind]
            temp_buoy[temp_buoy==5.0] = NaN # Spurious 5.0 degree readings.
            time_buoy = buoy_data['TIME']
            date_buoy = num2date(time_buoy[:],time_buoy.units)
            ax.plot(date_buoy,temp_buoy,'red',
                    label='{} at {} ({:d}m)'.format(station,print_lat_lon(buoy_lat,buoy_lon),int(depth_buoy)))    
            data.append((date_buoy,temp_buoy))

    # GOTM - ASM0
    with Dataset(os.path.join('results/p_sossta/no_assim','medsea_GOTM_{:d}{:02d}.nc'.format(year,month)),'r') as medsea_data:
        medsea_lat = medsea_data['lat'][:]
        medsea_lon = medsea_data['lon'][:]
        m = argmin(abs(medsea_lat-buoy_lat),axis=0)
        n = argmin(abs(medsea_lon-buoy_lon),axis=0)
        print('Closest medsea_GOTM grid point at (m,n) = ({},{})'.format(m,n), print_lat_lon(medsea_lat[m],medsea_lon[n]))
        #medsea_sst = medsea_data['sst']
        time_GOTM = medsea_data['time']
        date_GOTM = num2date(time_GOTM[:],time_GOTM.units)
        
        # Plot GOTM SST
        if plot_GOTM_SST:
            sst_GOTM = medsea_data['temp'][:,0,m,n]
            ax.plot(date_GOTM, sst_GOTM,'gray',linestyle='dashed',label='GOTM ASM0 SST at ' + \
                    print_lat_lon(medsea_lat[m],medsea_lon[n]))
            data.append((date_GOTM,sst_GOTM))       
        
        # 2016-12-19, temporary solution for plotting GOTM TEMP at buoy depth.
        from scipy.interpolate import RectBivariateSpline as interp
        depth_GOTM = medsea_data['depth'][:]
        temp_GOTM = medsea_data['temp'][:,:,m,n]
        f_GOTM = interp(time_GOTM,depth_GOTM,temp_GOTM,kx=3,ky=3, # bicubic interpolation
                        bbox=[0,time_GOTM[-1],0,max(depth_GOTM)])
        temp_GOTM_at_buoy_depth = f_GOTM(time_GOTM,depth_buoy)
        ax.plot(date_GOTM, temp_GOTM_at_buoy_depth,'gray', linestyle='solid',
                label='GOTM ASM0 TEMP ({0:d}m) '.format(int(depth_buoy)) + \
                print_lat_lon(medsea_lat[m],medsea_lon[n]))
        data.append((date_GOTM,temp_GOTM_at_buoy_depth))                            
            
    # GOTM - ASM2 
    with Dataset(os.path.join('results/p_sossta/assim_sunrise','medsea_GOTM_{:d}{:02d}.nc'.format(year,month)),'r') as medsea_data:
        medsea_lat = medsea_data['lat'][:]
        medsea_lon = medsea_data['lon'][:]
        m = argmin(abs(medsea_lat-buoy_lat),axis=0)
        n = argmin(abs(medsea_lon-buoy_lon),axis=0)
        print('Closest medsea_GOTM grid point at (m,n) = ({},{})'.format(m,n), print_lat_lon(medsea_lat[m],medsea_lon[n]))
        #medsea_sst = medsea_data['sst']
        time_GOTM = medsea_data['time']
        date_GOTM = num2date(time_GOTM[:],time_GOTM.units)
        
        # Plot GOTM SST
        if plot_GOTM_SST:
            sst_GOTM = medsea_data['temp'][:,0,m,n]
            ax.plot(date_GOTM, sst_GOTM,'blue',linestyle='dashed',label='GOTM ASM2 SST at ' + \
                    print_lat_lon(medsea_lat[m],medsea_lon[n]))
            data.append((date_GOTM,sst_GOTM))
        
        # 2016-12-19,  temporary solution for plotting GOTM TEMP at buoy depth.
        from scipy.interpolate import RectBivariateSpline as interp
        depth_GOTM = medsea_data['depth'][:]
        temp_GOTM = medsea_data['temp'][:,:,m,n]
        f_GOTM = interp(time_GOTM,depth_GOTM,temp_GOTM,kx=3,ky=3, # bicubic interpolation
                        bbox=[0,time_GOTM[-1],0,max(depth_GOTM)])
        temp_GOTM_at_buoy_depth = f_GOTM(time_GOTM,depth_buoy)
        ax.plot(date_GOTM, temp_GOTM_at_buoy_depth,'blue', linestyle='solid',
                label='GOTM ASM2 TEMP ({0:d}m) '.format(int(depth_buoy)) + \
                print_lat_lon(medsea_lat[m],medsea_lon[n]))
        data.append((date_GOTM,temp_GOTM_at_buoy_depth))
        
    # REA
    with Dataset(os.path.join('profiles','medsea_rea_votemper_{:d}{:02d}.nc'.format(year,month)),'r') as rea_data:
        votemper = rea_data['votemper']
        time_rea = rea_data['time']
        date_rea = num2date(time_rea[:],time_rea.units)
        temp_REA = votemper[:,0,m,n]
        if votemper[0,0,m,n] > 1e10: # the first (1st index) temperature record at shallowest depth (2nd index)
            print('Possible land location!')
        elif votemper[0,-1,m,n] > 1e10: # the first (1st index) temperature record at deepest depth (2nd index)
            print('Possible shallow water location!')
        
        depth_REA = rea_data['depth'][0]
        ax.plot(date_rea,temp_REA,'green', label='REA ({:.2f}m) at '.format(depth_REA) + \
                print_lat_lon(rea_data['lat'][m],rea_data['lon'][n]))
        data.append((date_rea,temp_REA))

    os.chdir(home)
    return ax.get_figure(), ax, data

def buoy_comparisons_full_year(year=2014,station='61277',showfig=True,fig_folder='fig/buoys'):
    fig, axes = subplots(ncols=1,nrows=12,figsize=(17,22))
    for i,ax in enumerate(axes):
        buoy_comparisons(ax=ax,year=year,month=i+1,station=station)
    
    # 1.4721m is the shallowest depth in REA data
    suptitle('GOTM-ASM0 SST (gray), GOTM-ASM2 SST (blue), Buoy data (red) vs REA at 1.4721m (green) for near {}'.format(station),fontsize=16,y=0.92)
    fig_fn = os.path.join(base_folder,fig_folder,'SST_GOTM_vs_REA_near_{}_{:d}.png'.format(station,year))
    print('Saving ' + fig_fn + '...')
    fig.savefig(fig_fn)
    if not(showfig):
        close(fig)

def SST_monthly_comparisons(year=2013,month=2,m=7,n=42,runs=['ASM0','ASM2'],colors=None,
                            ax=None,pretty=True,output_data=True,fig_fn=None):
    from netCDF4 import Dataset, num2date
    import matplotlib.cm as cm
    from matplotlib.pyplot import subplots
    from numpy import linspace

    if ax is None:
        fig, ax = subplots(figsize=(16,4))
        
    if colors is None:
        cfun = cm.rainbow
        ax.set_prop_cycle('color',[cfun(x) for x in linspace(0,1,len(runs))]) # Extra "run" for rea data
        colors = [cfun(i/len(runs)) for i in range(len(runs))]
        
    data = list()
    # GOTM runs
    for i, run in enumerate(runs):
        with Dataset(os.path.join(base_folder,run,'medsea_GOTM_{:d}{:02d}.nc'.format(year,month)),'r') as medsea_data:
            time_GOTM = medsea_data['time']
            date_GOTM = num2date(time_GOTM[:],time_GOTM.units)
            sst_GOTM = medsea_data['sst'][:,m,n]
            ax.plot(date_GOTM, sst_GOTM, color=colors[i], label=run)
            if output_data:
                data.append({run: (date_GOTM,sst_GOTM)})

    # REA
    with Dataset(os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_{:d}{:02d}.nc'.format(year,month)),'r') as rea_data:
        votemper = rea_data['votemper']
        time_rea = rea_data['time']
        date_rea = num2date(time_rea[:],time_rea.units)
        temp_REA = votemper[:,0,m,n]
        
        if votemper[0,0,m,n] > 1e10: # the first (1st index) temperature record at shallowest depth (2nd index)
            print('Possible land location!')
        elif votemper[0,-1,m,n] > 1e10: # the first (1st index) temperature record at deepest depth (2nd index)
            print('Possible shallow water location!')

        depth_REA = rea_data['depth'][0]
        ax.plot(date_rea,temp_REA,color='black',linewidth=2,label='REA'.format(depth_REA))
        ax.set_xlim(left=date_rea[0],right=date_rea[-1])
        ax.set_xticks(date_rea,minor=False) 
        ax.set_xticklabels([date.day for date in date_rea])
        if output_data:
            data.append({'rea': (date_rea,temp_REA)})

    fig = ax.get_figure()
    if pretty:
        # fig.autofmt_xdate()            
        ax.set_title('SST comparisons at {} for {:d}-{:02d}'.format(print_lat_lon(*get_lat_lon(m,n)),year,month))
        ax.grid('on')
        ax.legend()
    if output_data:
        return fig, ax, data 
    else:
        return fig, ax

def SST_yearly_comparisons(year,m,n,runs=['ASM0','ASM2'],showfig=True,
                           fig_subfolder = 'fig/spotchk'):
    from matplotlib.pyplot import subplots, suptitle, legend, close

    fig, axes = subplots(ncols=1,nrows=12,figsize=(17,22),sharex=False)
    for i,ax in enumerate(axes):
        SST_monthly_comparisons(year=year,month=i+1,m=m,n=n,runs=['ASM0','ASM2'],ax=ax,pretty=False,output_data=False)
        ax.grid('on')
        ax.set_xticklabels([])
        if i == 0:
            ax.legend()
            lines = ax.get_lines()

    latlon = print_lat_lon(*get_lat_lon(m,n))
    # 1.4721m is the shallowest depth in REA data
    suptitle('SST comparisons for the year {:d} at {:s}'.format(year, latlon),
             fontsize=16,y=0.92)
    labels = (*runs, 'REA at 1.4721m')
    # legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0,-0.1,1,1),
    #        bbox_transform = fig.transFigure )
    fig_fn = os.path.join(base_folder,fig_subfolder,'SST_{}_GOTM_ASM0_vs_ASM2_vs_REA_{:d}.png'.format(latlon, year))
    print('Saving ' + fig_fn + '...')
    fig.savefig(fig_fn)
    if not(showfig):
        close(fig)
