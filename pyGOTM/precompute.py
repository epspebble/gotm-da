# Local imports.
from os import getenv, mkdir
from os.path import isfile, isdir, isabs, join
from time import time
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime, date, timedelta
from pyGOTM.swr import swr_3hourly_mean
from numpy.ma import masked_all

# Hopefully grid has been set correctly...
from pyGOTM import medsea
grid = medsea.grid
data_folder = medsea.data_folder

def date2ndaysnsecs(dt):
    """
    Convert a datetime object to (ndays, nsecs), number of days 
    since last new year eve and number of seconds since midnight.
    """
    
    ndays = (date(dt.year,dt.month,dt.day) - date(dt.year-1,12,31)).days 
    nsecs = int((dt-datetime(dt.year,dt.month,dt.day,0,0,0)).total_seconds())
    return ndays, nsecs

def medsea_append_cloud_factor(year,month=None,thrifty=True,data_source='ECMWF'): # Grid info used to get sea_mn
    """
    Append downward short wave radiation 3-hourly means to the nc variable 'swrd_cs' in each file.
    
    Version date: 2017-11-19
    """
    
    ## Preparations.    

    # Assumed filename pattern.
    if month is None:
        dst_fn = 'medsea_{:s}_heat_{:d}.nc'.format(data_source,year)
    else:
        dst_fn = 'medsea_{:s}_heat_{:d}{:02d}.nc'.format(data_source,year,month)
    dst_folder = join(data_folder,'medsea_{:s}'.format(data_source))
    if not isdir(dst_folder) or not isfile(join(dst_folder,dst_fn)):
        raise(FileNotFoundError('The target directory must contain reformatted {:s} data.'.format(data_source)))
    
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
    tic = time()
    if thrifty:
        print('Precomputing swr_cs and cloud_factor values.')
        for j,(m,n) in enumerate(medsea.sea_mn):
            lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]
            #print('Calculating for lat,lon =',lat,',',lon)
            for i,dt in enumerate(full_dates):
                ndays, nsecs = date2ndaysnsecs(dt)
                swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
                if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                    cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
            elapsed = time() - tic
            avg_time = elapsed/(j+1)
            J = medsea.sea_m.size
            print('\t {:d}/{:d} completed. About {:.2f}s to go.'.format(j+1,medsea.sea_m.size,avg_time*(J-j-1)),end='\r')
        print('\nElapsed {:.2f}s.'.format(elapsed))
    else:
        print('Precomputing swr_cs and cloud_factor values.')        
        j=0
        J = medsea.M*medsea.N
        for m in range(medsea.M):
            for n in range(medsea.N):
                lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]                
                for i,dt in enumerate(full_dates):
                    ndays, nsecs = date2ndaysnsecs(dt)
                    swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
                    if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                        cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
                elapsed = time() - tic
                j += 1
                avg_time = elapsed/j
                print('\t {:d}/{:d} completed. About {:.2f}s to go.'.format(j,J,avg_time*(J-j)),end='\r')
        print('\nElapsed {:.2f}s.'.format(elapsed))
    
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

def medsea_ERA_append_cloud_factor(year,month=None,thrifty=True): # Grid info used to get sea_mn
    """
    Append downward short wave radiation 3-hourly means to the nc variable 'swrd_cs' in each file.
    
    Version date: 2017-11-19
    """
    
    ## Preparations.    

    # Assumed filename pattern.
    if month is None:
        dst_fn = 'medsea_ERA-INTERIM_heat_{:d}.nc'.format(year)
    else:
        dst_fn = 'medsea_ERA-INTERIM_heat_{:d}{:02d}.nc'.format(year,month)
    dst_folder = join(data_folder,'medsea_ERA-INTERIM')
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
    tic = time()
    if thrifty:
        print('Precomputing swr_cs and cloud_factor values.')
        for j,(m,n) in enumerate(medsea.sea_mn):
            lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]
            #print('Calculating for lat,lon =',lat,',',lon)
            for i,dt in enumerate(full_dates):
                ndays, nsecs = date2ndaysnsecs(dt)
                swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
                if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                    cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
            elapsed = time() - tic
            avg_time = elapsed/(j+1)
            J = medsea.sea_m.size
            print('\t {:d}/{:d} completed. About {:.2f}s to go.'.format(j+1,medsea.sea_m.size,avg_time*(J-j-1)),end='\r')
        print('\nElapsed {:.2f}s.'.format(elapsed))
    else:
        print('Precomputing swr_cs and cloud_factor values.')        
        j=0
        J = medsea.M*medsea.N
        for m in range(medsea.M):
            for n in range(medsea.N):
                lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]                
                for i,dt in enumerate(full_dates):
                    ndays, nsecs = date2ndaysnsecs(dt)
                    swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
                    if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                        cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
                elapsed = time() - tic
                j += 1
                avg_time = elapsed/j
                print('\t {:d}/{:d} completed. About {:.2f}s to go.'.format(j,J,avg_time*(J-j)),end='\r')
        print('\nElapsed {:.2f}s.'.format(elapsed))
    
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

def medsea_ECMWF_append_cloud_factor(year,month):
    """
    Append downward short wave radiation 3-hourly means to the nc variable 'swrd_cs' in each file.
    
    Version date: 2017-11-19
    """

    print('WARNING!!! DEPRECATED!!!')
    import sys
    sys.exit(1)
    
    # Assumed filename pattern.
    dst_fn = 'medsea_ECMWF_heat_{:d}{:02d}.nc'.format(year,month)    
    dst_folder = join(data_folder,'medsea_ECMWF')
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
    tic = time()
    print('Precomputing swr_cs and cloud_factor values.')
    for j,(m,n) in enumerate(medsea.sea_mn):
        lat, lon = medsea.grid_lats[m],medsea.grid_lons[n]
        # print('Calculating for lat,lon =',lat,',',lon)

        for i,dt in enumerate(full_dates):
            ndays, nsecs = date2ndaysnsecs(dt)
            swrd_cs[i,m,n] = swr_3hourly_mean(ndays,nsecs,lat,lon)
            if swrd_cs[i,m,n] > 0: # This implies fill_value for night time on sea locations. 
                cloud_factor[i,m,n] = swrd[i,m,n]/swrd_cs[i,m,n] 
        elapsed = time() - tic
        avg_time = elapsed/(j+1)
        J = medsea.sea_m.size
        print('\t {:d}/{:d} completed. About {:.2f}s to go.'.format(j+1,medsea.sea_m.size,avg_time*(J-j-1)),end='\r')
    print('\nElapsed {:.2f}s.'.format(elapsed))
    
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
