# __all__ = ['medsea_ERA_reformat',
#            'medsea_ECMWF_reformat',
#            'medsea_MFC_midnights_eformat',
#            'medsea_MFC_sunrise_reformat']

from pyGOTM.gotmks import subindices # Use this to subgrid a finer datset.
from pyGOTM import medsea
from os.path import isfile, isdir, join
from os import mkdir
from numpy import array, array_equal
from datetime import date, datetime, timedelta

p_sossta_folder = medsea.p_sossta_folder
data_folder = medsea.data_folder
dst_lats = medsea.grid_lats
dst_lons = medsea.grid_lons

# Precomputed indices map MFC reanalysis data to medsea 1x grid up to 75m deep.
nlev = 16

MFC_lev = \
          [1.4721018075942993,
           4.587478160858154,
           7.944124221801758,
           11.558627128601074,
           15.448707580566406,
           19.633302688598633,
           24.132646560668945,
           28.968355178833008,
           34.16352844238281,
           39.742835998535156,
           45.73264694213867,
           52.1611213684082,
           59.05834197998047,
           66.4564437866211,
           74.38976287841797,
           82.89495086669922]
dst_lev = MFC_lev

    
### ERA-INTERIM
# This improves current code in pyGOTM/ncdf_reformat.py, version on 2017-07-13
def make_abs(rel_path, side_effect = 'raise'):
    """
    Returns the absolute path by prepending user home folder. 
    Options:
    1. side_effect = 'raise' 
       Raise OSError if the directory does not exist.
    2. side_effect = 'mkdir'
       Assumes it is a folder and create it if it does not exist.
    Always raise an OSError if the path leads to an existing file.
    """
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join

    if isabs(rel_path):
        abs_path = rel_path
    else:
        abs_path = join(getenv('HOME'),rel_path)
    if isfile(abs_path):
        raise OSError('The path to a file is given, expecting a directory.')
    if not isdir(abs_path):
        # The mkdir() command could itself raise an OSError.
        if side_effect == 'mkdir':
            try:
                mkdir(abs_path)
            except OSError as oe:
                raise oe('Creating the directory {!s} failed. '.format(abs_path))
        elif side_effect == 'raise':
            raise OSError('The directory {!s} does not exist.'.format(abs_path))
    return abs_path

def get_subgrid_data(ds,varname,lat_slice,lon_slice,depth_slice=None):
    """ 
    Pass lat_slice=None, lon_slice=None to compute the slices, e.g. for the first time. 
    depth_slice is to be precomputed for now.
    """
    # First, compute the right indices if necessary, or confirm they're identical.
    if lat_slice is None or lon_slice is None:
        lat_slice = slice(*subindices(ds['lat'][:],dst_lats))
        lon_slice = slice(*subindices(ds['lon'][:],dst_lons))
        
        
    # Now read.
    assert array_equal(ds['lat'][lat_slice], dst_lats)
    assert array_equal(ds['lon'][lon_slice], dst_lons)
    if depth_slice is not None:
        assert array_equal(ds['depth'][depth_slice], dst_lev)
        data = ds[varname][:,depth_slice,lat_slice,lon_slice]
    else:
        # Now, read the data unpacked from the netCDF Dataset object.
        data = ds[varname][:,lat_slice,lon_slice]
    
    return data, lat_slice, lon_slice # Also return the lat_slice, lon_slice for reuse 

def get_interpolated_data(ds,varname):
    """ 
    Use ds[varname][i,:,:] to interpolate linearly to obtain data at i-th time record.
    The destination grid defined by 'dst_lats', 'dst_lons' globally in this module need
    to be strictly within the target Dataset's lat lon limits.
    """
    from scipy.interpolate import RectBivariateSpline
    from numpy import array

    lon = ds['lon'][:]
    lat = ds['lat'][:]
    
    lon_dir = 1 if (lon[1]-lon[0] > 0) else -1
    lat_dir = 1 if (lat[1]-lat[0] > 0) else -1

    lon_min, lon_max = (lon[0], lon[-1])[::lon_dir]
    lat_min, lat_max = (lat[0], lat[-1])[::lat_dir]

    #print(lon_min,lon_max,lat_min,lat_max)
    
    if any((lon_min > dst_lons[0],
           lon_max < dst_lons[-1],
           lat_min > dst_lats[0],
           lat_max < dst_lats[-1])):
        raise ValueError("The grid with 'lat' \in [{:g},{:g}], 'lon' \in [{:g},{:g}] do not cover the destination grid. \n Check metadata: \n{!s}\n{!s}".format(lat_min, lat_max, lon_min, lon_max, ds['lat'], ds['lon'] ))
    
    # Reorder in strictly ascending order for RectBivariateSpline to work.
    lon = ds['lon'][::lon_dir]
    lat = ds['lat'][::lat_dir]
    val = ds[varname][:,::lat_dir,::lon_dir]
    val_intp = list()
    for i in range(val.shape[0]):
        interp = RectBivariateSpline(lat,lon,val[i,:])
        val_intp.append(interp(dst_lats,dst_lons))
    return array(val_intp)

def get_src_dst_folders(dataset):
    assert isinstance(dataset,str), 'Check argument'
    
    if dataset == 'ERA-INTERIM':
        src_folder = join(p_sossta_folder,'medsea_ERA-INTERIM','3-hourly')        
        dst_folder = join(data_folder,'medsea_ERA-INTERIM')
    elif dataset == 'ECMWF':
        src_folder = join(p_sossta_folder,'medsea_ECMWF','3-HOURLY')
        dst_folder = join(data_folder,'medsea_ECMWF')
    elif dataset == 'MFC_midnights':
        src_folder = join(p_sossta_folder,'medsea_rea') # but stored in yearly subfolders
        dst_folder = join(data_folder, 'medsea_MFC_midnights')
    elif dataset == 'MFC_sunrise':
        src_folder = join(p_sossta_folder, 'medsea_rea','sunrise_nrt')
        dst_folder = join(data_folder,'medsea_MFC_sunrise')
    else:
        raise NotImplementedError('Unknown dataset: ' + dataset)

    if not isdir(src_folder):
        raise OSError('Invalid path: {!s}'.format(src_folder))
    
    if not isdir(dst_folder):
        try:
            mkdir(dst_folder)
        except:
            raise
    return src_folder, dst_folder

def medsea_ERA_reformat(year, grid='1x'):
    """ Combine variables from yearly ERA data into intermediate netCDF4 files, grouped by 
    whether they are to be eventually written to heat.dat or met.dat.
    """
    ## Preparations.
    
    # Local imports.
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join
    from time import time
    from datetime import datetime, timedelta
    from netCDF4 import Dataset, date2num
    from numpy import array_equal
    from pyGOTM import medsea
    from pyGOTM.gotmks import tic, toc
    from pyGOTM.config import epoch # Should be datetime(1981,1,1,0,0,0)    
    
    # # Find absolute path for the src_folder and dst_folder
    # if not isabs(src_folder):
    #     src_folder = join(getenv('HOME'),src_folder)
    # if not isabs(dst_folder):
    #     dst_folder = join(getenv('HOME'),dst_folder)
    #     if not isdir(dst_folder):
    #         mkdir(dst_folder)
    
    dataset = 'ERA-INTERIM'
    src_folder, dst_folder = get_src_dst_folders(dataset)
    
    # Filename keyword 'name' / internal nc variable name 'alias' (when they don't agree)
    names = dict(met = ['u10m','v10m','sp','t2m','q2m','precip'],
                 heat = ['lwrd','swrd'])
    aliases = dict(precip='var144', snow='var228',
                   lwrd='var175', swrd='var169')

    # Convenience function to get one variable's yearly data.
    lat_slice, lon_slice, do_interp = None, None, False
    def get_data(name,src_fn):
        nonlocal lat_slice, lon_slice, do_interp
        with Dataset(join(src_folder,src_fn),'r') as ds:
            # Switch to the correct netCDF4 variable name if necessary
            if name in aliases.keys():
                varname = aliases[name]
            else:
                varname = name
            if do_interp:
                data = get_interpolated_data(ds,varname)
            else:
                try: 
                    data = get_subgrid_data(ds,varname,lat_slice,lon_slice)
                except ValueError as ve:
                    print('ValueError occurred using get_subgrid_data():')
                    print(ve)
                    print('Attempting to interpolate instead...')
                    do_interp = True
                    data = get_interpolated_data(ds,varname)
        return data

    # if grid == '1x':
    #     from pyGOTM import medsea
    #     medsea.set_grid('1x')
    #     # Precomputed indices to map ERA-INTERIM to GOTM 1x grid 
    #     src_lat_idx = slice(-8,4,-1)
    #     src_lon_idx = slice(12,-18,1)
    #     dst_lats = medsea.grid_lats
    #     dst_lons = medsea.grid_lons
    # else:
    #     raise NotImplementedError('Me no do grid = {!s} (yet).'.format(grid))

    ## Creating intermediate file on 1x grid for creating heat.dat / met.dat later.
    tic = time()
    elapsed = 0
    for dat_type in ['heat','met']:
        dst_fn = 'medsea_ERA_{:s}_{:d}.nc'.format(dat_type, year)
        print('Writing {:s}...'.format(join(dst_folder,dst_fn)))
        # Create the dimensions
        with Dataset(join(dst_folder,dst_fn),'w',format="NETCDF3_CLASSIC") as ds:
            ds.createDimension('time') # unlimited
            ds.createDimension('lat', size=len(dst_lats))
            ds.createDimension('lon', size=len(dst_lons))

            ds_time = ds.createVariable('time','i4',dimensions=('time',))
            ds_time.units = 'hours since {!s}'.format(epoch) # Data given in 3-hourly periods since UTC midnight.
            ds_lat = ds.createVariable('lat','f4',dimensions=('lat',))
            ds_lon = ds.createVariable('lon','f4',dimensions=('lon',))

            # Write the values of the lat/lon dimensions using user-specified grid when loading the medsea module.
            ds_lat[:] = dst_lats
            ds_lon[:] = dst_lons

            # Write the values of the record times BY MAKING ASSUMPTIONS.
            # CHECK: The first ERA record is at new year day, 03:00:00, last record at next new year day midnight.
            if dat_type == 'met':
                # For met.dat records, we use the same timestamp, indicating the preceding interval:
                first_record = datetime(year,1,1,3,0,0)
            elif dat_type == 'heat':
                # For heat.dat records, we move the timestamp 3-hours earlier to indicate the subsequent interval:
                first_record = datetime(year,1,1,0,0,0)
            hour_offset = (first_record - epoch).total_seconds() // 3600
            yrdays = 366 if year%4 == 0 else 365
            ds_time[:] = [hour_offset+3*i for i in range(yrdays*8)] # 8 records per day.
        toc = time()-tic
        elapsed += toc
        print('Finished writing dimensions. Elapsed {!s}s'.format(toc))
    
        # Append the data one variable after another.
        for name in names[dat_type]:
            tic = time()
            src_fn = 'MEDSEA_ERA-INT_{:s}_y{:d}.nc'.format(name,year)
            # Append data for one variable for one year.
            with Dataset(join(dst_folder,dst_fn),'a') as ds:
                var = ds.createVariable(name,'f8',dimensions=('time','lat','lon'))
                # The lat_slice and lon_slice is computed in the first run, and shold
                # stay the same after every function call if assumption that
                # the lat/lon grid of the source dataset is unique and identical
                # across nc files.
                val = get_data(name,src_fn) 
                if name == 't2m':
                    var[:] = val - 273.15 # Convert from degrees K to degrees C
                elif name == 'sp':
                    var[:] = val / 100 # Convert from Pa to hPa
                else:
                    var[:] = val
            toc = time()-tic
            elapsed += toc
            print('Finished appending data for {:s}. Elapsed {!s}s'.format(name,toc))

        print('Total time: {!s}s'.format(elapsed))

### ECMWF
# This improves current code in pyGOTM/ncdf_reformat.py, version on 2017-07-13
def medsea_ECMWF_reformat(year,month,grid='1x'):
#                          src_folder='p_sossta/medsea_ECMWF/3-HOURLY',
#                          dst_folder='medsea_data/medsea_ECMWF'):    
    """ Combine variables from monthly ECMWF data into intermediate netCDF4 files, grouped by 
    whether they are to be eventually written to heat.dat or met.dat.
    """
    ## Preparations.
    
    # Local imports.
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join
    from time import time
    from datetime import datetime, timedelta
    from netCDF4 import Dataset, date2num
    from numpy import array_equal    
    from pyGOTM.gotmks import tic, toc
    from pyGOTM.config import epoch # Should be datetime(1981,1,1,0,0,0)    
    
    # # Find absolute path for the src_folder and dst_folder
    # if not isabs(src_folder):
    #     src_folder = join(getenv('HOME'),src_folder)
    # if not isabs(dst_folder):
    #     dst_folder = join(getenv('HOME'),dst_folder)
    #     if not isdir(dst_folder):
    #         mkdir(dst_folder)

    src_folder = join(p_sossta_folder,'medsea_ECMWF','3-HOURLY')
    dst_folder = join(data_folder,'medsea_ECMWF')
    if not isdir(dst_folder):
        mkdir(dst_folder)

    
    # Filename keyword 'name' / internal nc variable name 'alias' (when they don't agree)
    names = dict(met = ['u10m','v10m','sp','t2m','q2m','precip'],
                 heat = ['lwrd','swrd'])
    aliases = dict() # Checked: they are identical.

    ## No longer needed when we can compute on the fly.
    # if grid == '1x':
    #     from pyGOTM import medsea
    #     medsea.set_grid(grid)
    #     # Precomputed indices to map ECMWF data to medsea 1x grid.
    #     src_lat_idx = slice(154,33,-6) 
    #     src_lon_idx = slice(72,409,6)
    #     dst_lats = medsea.grid_lats
    #     dst_lons = medsea.grid_lons
    # else:
    #     raise NotImplementedError('Me no do grid = {!s} (yet).'.format(grid))

    # Convenience function to get one variable's yearly data.
    lat_slice, lon_slice, do_interp = None, None, False
    def get_data(name,src_fn):
        nonlocal lat_slice, lon_slice, do_interp
        # Switch to the correct netCDF4 variable name if necessary
        if name in aliases.keys():
            varname = aliases[name]
        else:
            varname = name
        # Now read in the data by appealing to subgridding or interpolation.
        with Dataset(join(src_folder,src_fn),'r') as ds:
            if do_interp:
                data = get_interpolated_data(ds,varname)
            else:
                try: 
                    data, lat_slice, lon_slince = get_subgrid_data(ds,varname,lat_slice,lon_slice)
                except ValueError as ve:
                    print('ValueError occurred using get_subgrid_data():')
                    print(ve)
                    print('Attempting to interpolate instead...')
                    do_interp = True
                    data = get_interpolated_data(ds,varname)
        return data

    # # Convenience function to get one variable's yearly data.
    # lat_slice, lon_slice = None, None
    # def get_data(name,src_fn):
    #     nonlocal lat_slice, lon_slice
    #     with Dataset(join(src_folder,src_fn),'r') as ds:
            
    #         # Compute the right indices for the first time, or confirm they're correct.
    #         if lat_slice is None or lon_slice is None:
    #             lat_slice = slice(*subindices(ds['lat'][:],dst_lats))
    #             lon_slice = slice(*subindices(ds['lon'][:],dst_lons))
    #         assert array_equal(ds['lat'][lat_slice], dst_lats)
    #         assert array_equal(ds['lon'][lon_slice], dst_lons)

    #         # Switch to the correct netCDF4 variable name if necessary
    #         if name in aliases.keys():
    #             alias = aliases[name]
    #         else:
    #             alias = name

    #         # Then return the data unpacked from the netCDF object.
    #         data = ds[alias][:,lat_slice,lon_slice]

    #     return data, lat_slice, lon_slice # Also return the lat_slice, lon_slice for use next time.
    
    ## Creating intermediate file on 1x grid for creating heat.dat / met.dat later.
    tic = time()
    elapsed = 0
    for dat_type in ['heat','met']:        
        dst_fn = 'medsea_ECMWF_{:s}_{:d}{:02d}.nc'.format(dat_type,year,month)
        print('Writing {:s}...'.format(join(dst_folder,dst_fn)))
        # Create the dimensions
        with Dataset(join(dst_folder,dst_fn),'w',format="NETCDF3_CLASSIC") as ds:
            ds.createDimension('time') # unlimited
            ds.createDimension('lat', size=len(dst_lats))
            ds.createDimension('lon', size=len(dst_lons))

            ds_time = ds.createVariable('time','i4',dimensions=('time',))
            ds_time.units = 'hours since {!s}'.format(epoch) # Data given in 3-hourly periods since UTC midnight.
            ds_lat = ds.createVariable('lat','f4',dimensions=('lat',))
            ds_lon = ds.createVariable('lon','f4',dimensions=('lon',))

            # Write the values of the lat/lon dimensions using user-specified grid when loading the medsea module.
            ds_lat[:] = dst_lats
            ds_lon[:] = dst_lons

            # Write the values of the record times BY MAKING ASSUMPTIONS.
            # CHECK: The first ECMWF record is at new year day, 03:00:00, last record at next new year day midnight.
            if dat_type == 'met':
                # For met.dat records, we use the same timestamp, indicating the preceding interval:
                first_record = datetime(year,month,1,3,0,0)
            elif dat_type == 'heat':
                # For heat.dat records, we move the timestamp 3-hours earlier to indicate the subsequent interval:
                first_record = datetime(year,month,1,0,0,0)
            hour_offset = (first_record - epoch).total_seconds() // 3600
#             yrdays = 366 if year%4 == 0 else 365
            next_month = datetime(year,month+1,1) if month < 12 else datetime(year+1,1,1)
            mthdays = (next_month-datetime(year,month,1)).days
            ds_time[:] = [hour_offset+3*i for i in range(mthdays*8)] # 8 records per day.
            
        toc = time()-tic
        elapsed += toc
        print('Finished writing dimensions. Elapsed {!s}s'.format(toc))
    

        # Append the data one variable after another.
        lat_slice, lon_slice = None, None # To instruct get_data to calculate them.
        for name in names[dat_type]:
            tic = time()
            src_fn = 'ECMWF_{:s}_y{:d}m{:02d}.nc'.format(name,year,month)
            # Append data for one variable for one year.
            with Dataset(join(dst_folder,dst_fn),'a') as ds:
                var = ds.createVariable(name,'f8',dimensions=('time','lat','lon'))
                val = get_data(name,src_fn)
                if name == 't2m':
                    var[:] = val - 273.15 # Convert from degrees K to degrees C
                elif name == 'sp':
                    var[:] = val / 100 # Convert from Pa to hPa
                else:
                    var[:] = val
            toc = time()-tic
            elapsed += toc
            print('Finished appending data for {:s}. Elapsed {!s}s'.format(name,toc))

        print('Total time: {!s}s'.format(elapsed))

### MFC_midnights
def medsea_MFC_midnights_reformat(year,month=None,grid='1x'):
    """ Combine daily MFC midnightly mean TEMP / PSAL data into intermediate yearly or monthly netCDF4 files, 
    to be eventually written to tprof.dat / sprof.dat.
    """
    ## Preparations.
    
    # Local imports.
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join
    from time import time
    from datetime import datetime, timedelta
    from netCDF4 import Dataset, date2num, MFDataset
    from numpy import array_equal, array    
    from pyGOTM.gotmks import tic, toc
    from pyGOTM.config import epoch # Should be datetime(1981,1,1,0,0,0)    


    src_folder, dst_folder = get_src_dst_folders('MFC_midnights')
    # MFC reanalysis data stored in yearly subfolders.
    src_folder = join(src_folder,str(year))
    
    # # Find absolute path for the src_folder and dst_folder
    # if not isabs(src_folder):
    #     src_folder = join(getenv('HOME'),src_folder)
    # if not isabs(dst_folder):
    #     dst_folder = join(getenv('HOME'),dst_folder)
    #     if not isdir(dst_folder):
    #         mkdir(dst_folder)
            
    ## Hard-coding some information for the subfunction get_ERA_yearly_data()

    # Filename keyword 'name' / internal nc variable 'alias' (when they don't agree)
    names = dict(tprof = 'TEMP',
                 sprof = 'PSAL')
    src_varname = dict(TEMP = 'votemper',
                       PSAL = 'vosaline')
    dst_varname = dict(TEMP = 'temp',
                       PSAL = 'salt')

    depth_slice = slice(None,nlev,None)
    lat_slice, lon_slice = None, None

    # Convenience function to get one variable's yearly data.
    def get_data(name,src_fn):
        nonlocal lat_slice, lon_slice, depth_slice
        # Use MFDataset instead of Dataset, src_fn is a glob-able pattern instead of a specific filename. 
        print('Reading from {!s}...'.format(join(src_folder,src_fn)))
        # Switch to the correct netCDF4 variable name if necessary
        if name in src_varname.keys():
            varname = src_varname[name]
        else:
            varname = name
        with MFDataset(join(src_folder,src_fn),'r') as ds:
            data, lat_slice, lon_slice = get_subgrid_data(ds,varname,lat_slice,lon_slice,depth_slice)
        return data
    
    ## Creating intermediate file on 1x grid for creating heat.dat / met.dat later.
    tic = time()
    elapsed = 0
    for dat_type in ['tprof','sprof']:
        dst_fn = 'medsea_MFC_midnights_' + dat_type + '_' 
        dst_fn += '{:d}{:02d}.nc'.format(year,month) if month is not None else '{:d}.nc'.format(year)
        
        print('Writing {:s}...'.format(join(dst_folder,dst_fn)))
        # Create the dimensions
        with Dataset(join(dst_folder,dst_fn),'w',format="NETCDF3_CLASSIC") as ds:
            ds.createDimension('time') # unlimited
            ds.createDimension('depth', size=nlev)             
            ds.createDimension('lat', size=len(dst_lats))
            ds.createDimension('lon', size=len(dst_lons))

            ds_time = ds.createVariable('time','i4',dimensions=('time',))
            ds_time.units = 'hours since {!s}'.format(epoch) # Data given in 3-hourly periods since UTC midnight.
            ds_depth = ds.createVariable('depth','f4',dimensions=('depth',))
            ds_lat = ds.createVariable('lat','f4',dimensions=('lat',))
            ds_lon = ds.createVariable('lon','f4',dimensions=('lon',))

            # Write the depth values using precomputed values.
            ds_depth[:] = dst_lev

            # Write the values of the lat/lon dimensions using user-specified grid when loading the medsea module.
            ds_lat[:] = dst_lats
            ds_lon[:] = dst_lons

            # Write the values of the record times BY MAKING ASSUMPTIONS.
            # CHECK: The first MFC record is at new year day, 00:00:00, last record at new year eve midnight.
            first_record = datetime(year,1,1,0,0,0)
            hour_offset = (first_record - epoch).total_seconds() // 3600
            if month is not None:
                next_month = date(year,month+1,1) if month < 12 else date(year+1,month,1)
                mthdays = (next_month - date(year,month,1)).days
                ds_time[:] = [hour_offset+24*i for i in range(mthdays)]
            else:
                yrdays = 366 if year%4 == 0 else 365
                ds_time[:] = [hour_offset+24*i for i in range(yrdays)]
            
        toc = time()-tic
        elapsed += toc
        print('Finished writing dimensions. Elapsed {!s}s'.format(toc))
    
        # Append the data.
        tic = time()
        with Dataset(join(dst_folder,dst_fn),'a') as ds:
            name = names[dat_type] # Just TEMP / PSAL
            
            # Filename pattern instead of specific filename
            src_fn = '{:d}{:02d}??'.format(year,month) if month else '{:d}????'.format(year)
            src_fn += '_' + name + '_re-fv6.nc'
            var = ds.createVariable(dst_varname[name],'f8',dimensions=('time','depth','lat','lon'))
            var[:] = get_data(name,src_fn)
            
        toc = time()-tic
        elapsed += toc
        print('Finished appending data for {:s}. Elapsed {!s}s'.format(dat_type,toc))

        print('Total time: {!s}s'.format(elapsed))


### MFC_sunrise
def medsea_MFC_sunrise_reformat(start_day,stop_day,grid='1x',
                                src_folder='p_sossta/medsea_rea/sunrise_nrt',
                                dst_folder='medsea_data/medsea_MFC_sunrise'):
    """ Combine hourly MFC TEMP / PSAL profile chosen closest to sunrise into intermediate yearly or monthly netCDF4 files, 
    to be eventually written to tprof.dat / sprof.dat.
    
    Version date: 2017-11-09.
    """
    ## Preparations.
    
    # Local imports.
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join
    from time import time
    from datetime import datetime, timedelta, date
    from netCDF4 import Dataset, date2num, MFDataset
    from numpy import array_equal, array
    from numpy.ma import masked_all
    from glob import glob
    from pyGOTM.gotmks import tic, toc
    from pyGOTM.swr import solar_times
    from pyGOTM.config import epoch # Should be datetime(1981,1,1,0,0,0)    

    src_folder, dst_folder = get_src_dst_folders('MFC_midnights')
    
    # # MFC reanalysis data stored in yearly subfolders.
    # #     src_folder = join(src_folder,str(year))
    
    # # Find absolute path for the src_folder and dst_folder
    # if not isabs(src_folder):
    #     src_folder = join(getenv('HOME'),src_folder)
    # if not isabs(dst_folder):
    #     dst_folder = join(getenv('HOME'),dst_folder)
    #     if not isdir(dst_folder):
    #         mkdir(dst_folder)
            
    # # User need to make sure the medsea module is loaded and set up.
    # import sys
    # assert 'medsea' not in sys.modules, "'import pyGOTM.medsea as medsea' first!'"
    
    # Restrict the date object type for clarity.
    assert isinstance(start_day,date)
    assert isinstance(stop_day,date)


    ## Hard-coding some information for the subfunction get_ERA_yearly_data()

    # Filename keyword 'name' / internal nc variable 'alias' (when they don't agree)
    names = dict(tprof = 'TEMP',
                 sprof = 'PSAL')
    src_varname = dict(TEMP = 'votemper',
                       PSAL = 'vosaline')
    dst_varname = dict(TEMP = 'temp',
                       PSAL = 'salt')
    if grid == '1x':
        from pyGOTM import medsea
        medsea.set_grid(grid)
        # Precomputed indices map MFC reanalysis data to medsea 1x grid up to 75m deep.
        nlev = 16
        src_lev = \
                  [1.4721018075942993,
                   4.587478160858154,
                   7.944124221801758,
                   11.558627128601074,
                   15.448707580566406,
                   19.633302688598633,
                   24.132646560668945,
                   28.968355178833008,
                   34.16352844238281,
                   39.742835998535156,
                   45.73264694213867,
                   52.1611213684082,
                   59.05834197998047,
                   66.4564437866211,
                   74.38976287841797,
                   82.89495086669922]
        
        assert len(src_lev) == nlev
        src_lat_idx = slice(9,None,12)
        # This is different from MFC_midnights because the lon starts from -15 instead of -6, so there is 9 / 0.0625 = 144 index shift
        src_lon_idx = slice(144,None,12) 
        
        # Target grid
        dst_lats = medsea.grid_lats
        dst_lons = medsea.grid_lons
        dst_lev = src_lev
       
    else:
        raise NotImplementedError('Me no do grid = {!s} (yet).'.format(grid))

    # Convenience function to get one variable's data.
    def get_data(name,src_fn):
#         print('Reading from {!s}...'.format(join(src_folder,src_fn)))
        with Dataset(join(src_folder,src_fn),'r') as ds:
            # Unfortunately, all the dimensional data were missing, and the dimensions
            # for lon is different!
            #assert array_equal(ds['lat'][src_lat_idx], dst_lats)
            #assert array_equal(ds['lon'][src_lon_idx], dst_lons)
            #assert array_equal(ds['depth'][:nlev],src_lev), ds['depth'][:nlev]
            # Switch to the correct netCDF4 variable name if necessary
            if name in src_varname.keys():
                alias = src_varname[name]
            else:
                alias = name
            # The data lacks the time dimension.
            # Then return the data unpacked from the netCDF object.
            data = ds[alias][:nlev,src_lat_idx,src_lon_idx]
        return data
    
    ## Creating intermediate file on 1x grid for creating heat.dat / met.dat later.
    tic = time()
    elapsed = 0
    for dat_type in ['tprof','sprof']:
        dst_fn = 'medsea_MFC_sunrise_' + dat_type + '_' 
        dst_fn += '{0.year:02d}{0.month:02d}{0.day:02d}'.format(start_day)
        dst_fn += '{0.year:02d}{0.month:02d}{0.day:02d}'.format(stop_day)
        dst_fn += '.nc'
        
        print('Writing {:s}...'.format(join(dst_folder,dst_fn)))
        # Create the dimensions
        with Dataset(join(dst_folder,dst_fn),'w',format="NETCDF3_CLASSIC") as ds:
            ds.createDimension('day') # aggregating dimension, UTC day.
            ds.createDimension('depth', size=nlev)             
            ds.createDimension('lat', size=len(dst_lats))
            ds.createDimension('lon', size=len(dst_lons))

            # Just the days from 0, 1, 2, ...
            ds_day = ds.createVariable('day','i4',dimensions=('day',))
            ds_day.units = "days since {!s}".format(epoch) # When the data has holes, this is safest.
            ndays = (stop_day-start_day).days
            ds_day[:] = range(0,ndays)
            
            # Write the depth values using precomputed values.
            ds_depth = ds.createVariable('depth','f4',dimensions=('depth',))
            ds_depth[:] = src_lev

            # Write the values of the lat/lon dimensions using user-specified grid when loading the medsea module.
            ds_lat = ds.createVariable('lat','f4',dimensions=('lat',))
            ds_lon = ds.createVariable('lon','f4',dimensions=('lon',))
            ds_lat[:] = dst_lats
            ds_lon[:] = dst_lons
            
        toc = time()-tic
        elapsed += toc
        print('Finished writing dimensions. Elapsed {!s}s'.format(toc))

        # Read the data and build up the tprof data array and sunrise times to be written.
        tic = time()
        
        data = masked_all((ndays,nlev,dst_lats.size,dst_lons.size))
        seconds_of_sunrise = masked_all((ndays,dst_lats.size,dst_lons.size))
        assert stop_day == start_day + timedelta(days=ndays)
        
        name = names[dat_type] # Just TEMP / PSAL            
        for i in range(ndays):
            current_day = start_day+timedelta(days=i)
            src_fn = 'sunrise_{0:s}_{1.year:d}{1.month:02d}{1.day:02d}.nc'.format(name,current_day)
            data[i,:,:,:] = get_data(name,src_fn)
            for m,n in medsea.sea_mn:
                lat, lon = dst_lats[m], dst_lons[n]
                seconds_of_sunrise[i,m,n] = solar_times(current_day,lat,lon,events='sunrise')*60
        toc = time()-tic
        elapsed += toc
        print('Finished reading data and calculating sunrise times. Elapsed {!s}s'.format(toc))

        # Append the data.
        tic = time()
        with Dataset(join(dst_folder,dst_fn),'a') as ds:
            ds_sec = ds.createVariable('seconds','i4',dimensions=('day','lat','lon'))
            ds_sec.units = "seconds since UTC midnight of the day"   
            ds_sec[:] = seconds_of_sunrise
            # Filename pattern instead of specific filename
            var = ds.createVariable(dst_varname[name],'f4',dimensions=('day','depth','lat','lon'))
            var[:] = data
            
        toc = time()-tic
        elapsed += toc
        print('Finished appending data for {:s}. Elapsed {!s}s'.format(dat_type,toc))

        print('Total time: {!s}s'.format(elapsed))

if __name__ == '__main__':
    import sys
    from datetime import date, datetime, timedelta
    from pyGOTM import medsea
    grid = '1x' # Default testing grid.

    if len(sys.argv) == 1:
        grid = '1x'
        datasets = ['ERA','ECMWF','MFC_midnights','MFC_sunrise']
    elif len(sys.argv) == 2:
        grid = sys.argv[1]
        datasets = ['ERA','ECMWF','MFC_midnights','MFC_sunrise']
    elif len(sys.argv) == 3:
        grid = sys.argv[1]
        datasets = [sys.argv[2]]
    else:
        print(
        """
        Usage:
              python reformat.py {grid} {ERA|ECMWF|MFC_midnights|MFC_sunrise}
        """)
        raise RuntimeError('Wrong number of arguments.')

    print('Attempting to regenerate reformatted data for medsea_{!s}...'.format(grid))
    input('Press any key to continue...')
    medsea.set_grid(grid)
            
    ## ERA & MFC midnightly means for years 2013 and 2014
    medsea_ERA_reformat(2013)
    medsea_ERA_reformat(2014)
    medsea_MFC_midnights_reformat(2013)
    medsea_MFC_midnights_reformat(2014)

    ## ECMWF and MFC sunrise profiles for 2016-04 to 2017-07
    for i in range(4,13):
        medsea_ECMWF_reformat(2016,i)
    for i in range(1,8):
        medsea_ECMWF_reformat(2017,i)

    # Date range #1
    start_day = date(2016,4,5)
    stop_day = date(2016,4,11) 
    medsea_MFC_sunrise_reformat(start_day,stop_day)

    # Date range #2
    start_day = date(2016,4,19)
    stop_day = date(2016,5,2) 
    medsea_MFC_sunrise_reformat(start_day,stop_day)

    # Date range #3
    start_day = date(2016,5,10)
    stop_day = date(2016,11,17) 
    medsea_MFC_sunrise_reformat(start_day,stop_day)

    # Date range #4
    start_day = date(2016,11,19)
    stop_day = date(2017,7,24) 
    medsea_MFC_sunrise_reformat(start_day,stop_day)
