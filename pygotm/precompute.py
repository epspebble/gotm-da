def date2ndaysnsecs(dt):
    """
    Convert a datetime object to (ndays, nsecs), number of days 
    since last new year eve and number of seconds since midnight.
    """
    from datetime import datetime,date
    
    ndays = (date(dt.year,dt.month,dt.day) - date(dt.year-1,12,31)).days 
    nsecs = int((dt-datetime(dt.year,dt.month,dt.day,0,0,0)).total_seconds())
    return ndays, nsecs

def medsea_ERA_append_cloud_factor(year,,grid='1x', # Grid info used to get sea_mn
                                   dst_folder='medsea_data/medsea_ERA-INTERIM'):
    """
    Append downward short wave radiation 3-hourly means to the nc variable 'swrd_cs' in each file.
    
    Version date: 2017-11-11
    """
    
    ## Preparations.
    
    # Local imports.
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join
    from time import time
    from netCDF4 import Dataset, date2num
    from pygotm.swr import swr_3hourly_mean

    # Assumed filename pattern.
    dst_fn = 'medsea_ERA_heat_{:d}.nc'.format(year)
    
    # Find absolute path for the src_folder and dst_folder
    if not isabs(dst_folder):
        dst_folder = join(getenv('HOME'),dst_folder)
    if not isdir(dst_folder) or not isfile(join(dst_folder,dst_fn)):
        raise(FileNotFoundError('The target directory must contain reformatted ERA data.'))
    
    print('Appending swrd_cs and cloud_factor to {:s}...'.format(join(dst_folder,dst_fn)))
    with Dataset(join(dst_folder,dst_fn),'r') as ds:
        full_dates = num2date(ds['time'][:],ds['time'].units)
        swrd = ds['swrd'][:]
        
    # Create a masked array and compute only for sea locations to save time.
    swrd_cs = masked_all(swrd.shape)
    cloud_factor = masked_all(swrd.shape)
    elapsed = 0

    # Grab some indices and lat/lon info about the medsea...
    # ... coz they are NOT available from ANY sort atmospheric forecast data.
    from pygotm import medsea
    medsea.set_grid('1x')
    for m,n in medsea.sea_mn:
        tic = time()
        lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]
        print('Calculating for lat,lon =',lat,',',lon)
        for i,dt in enumerate(full_dates):
            ndays, nsecs = date2ndaysnsecs(dt)
            swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
            if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
        elapsed += time() - tic
        print('Elapsed {:.2f}s'.format(elapsed))
    
    # Append the data
    tic = time()
    print('Writing to ',join(dst_folder,dst_fn))
    with Dataset(join(dst_folder,dst_fn),'a') as ds:
        var_sc = ds.createVariable('swrd_cs','f8',dimensions=('time','lat','lon'))
        var_sc[:] = swrd_cs
        var_cf = ds.createVariable('cloud_factor','f8',dimensions=('time','lat','lon'))
        var_cf[:] = cloud_factor
    elapsed = time() - tic
    
    print('Total time: {:.2f}s'.format(elapsed))

def medsea_ECMWF_append_cloud_factor(year,month,grid='1x', # Grid info used to get sea_mn
                                     dst_folder='medsea_data/medsea_ECMWF'):
    """
    Append downward short wave radiation 3-hourly means to the nc variable 'swrd_cs' in each file.
    
    Version date: 2017-11-11
    """
    
    ## Preparations.
    
    # Local imports.
    from os import getenv, mkdir
    from os.path import isfile, isdir, isabs, join
    from time import time
    from netCDF4 import Dataset, num2date
    from pygotm.swr import swr_3hourly_mean
    from numpy.ma import masked_all

    # Assumed filename pattern.
    dst_fn = 'medsea_ECMWF_heat_{:d}{:02d}.nc'.format(year,month)    
    
    # Find absolute path for the src_folder and dst_folder
    if not isabs(dst_folder):
        dst_folder = join(getenv('HOME'),dst_folder)
    if not isdir(dst_folder) or not isfile(join(dst_folder,dst_fn)):
        raise(FileNotFoundError('The target directory must contain reformatted ECMWF data: {!s}'.format(join(dst_folder,dst_fn))))
    
    print('Appending swrd_cs and cloud_factor to {:s}...'.format(join(dst_folder,dst_fn)))
    with Dataset(join(dst_folder,dst_fn),'r') as ds:
        full_dates = num2date(ds['time'][:],ds['time'].units)
        swrd = ds['swrd'][:]
        
    # Create a masked array and compute only for sea locations to save time.
    swrd_cs = masked_all(swrd.shape)
    cloud_factor = masked_all(swrd.shape)
    elapsed = 0

    # Grab some indices and lat/lon info about the medsea...
    # ... coz they are NOT available from ANY sort atmospheric forecast data.
    from pygotm import medsea
    medsea.set_grid('1x')
    for m,n in medsea.sea_mn:
        tic = time()
        lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]
        print('Calculating for lat,lon =',lat,',',lon)
        for i,dt in enumerate(full_dates):
            ndays, nsecs = date2ndaysnsecs(dt)
            swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
            if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
        elapsed += time() - tic
        print('Elapsed {:.2f}s'.format(elapsed))
    
    # Append the data
    tic = time()
    print('Writing to ',join(dst_folder,dst_fn))
    with Dataset(join(dst_folder,dst_fn),'a') as ds:
        var_sc = ds.createVariable('swrd_cs','f8',dimensions=('time','lat','lon'))
        var_sc[:] = swrd_cs
        var_cf = ds.createVariable('cloud_factor','f8',dimensions=('time','lat','lon'))
        var_cf[:] = cloud_factor
    elapsed = time() - tic
    
    print('Total time: {:.2f}s'.format(elapsed))
