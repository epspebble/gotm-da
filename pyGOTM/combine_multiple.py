from pyGOTM import medsea
from pyGOTM.combine import reload, hour_range, outfn

def read(varnames,fn=None,hr=None,year=None,month=None,failed_list=None):
    from time import time
    from os.path import join
    from netCDF4 import Dataset
    from numpy.ma import masked_all, masked_invalid
    
    # Handle arguments
    if fn is None:
        fn = outfn(year,month)
           
    # If not provided, compute 'hr' from one of the nc files.
    # If a specific date range was used to generate GOTM output, it is better to recalculate hr outside of this function.
    if hr is None:

        assert year is not None, 'At least the year should be specified .'
        # Use the first readable file's metadata to determine start time, assuming GOTM's unit for 'time' in the ncdf file is the start time.
        i=0
        start = None
        from netCDF4 import num2date
        while start is None:
            m,n = medsea.sea_m[i], medsea.sea_n[i]            
            try:
                with Dataset(join(medsea.get_local_folder(m,n),fn),'r') as ds:
                    start = num2date(0,ds['time'].units)
            except Exception as e:
                if i == medsea.sea_m.size-1:
                    raise e
                else:
                    i += 1
                    pass
        hr = hour_range(year,month,start=start) # Here month might be None.

    # Use medsea global settings and 'hr' to determine array dimensions.
    assert not isinstance(varnames,str), print('Must be a list of variables to combine.')
    
    # Do not allow these two 4D variables to be lumped in one file.
    assert 'temp' not in varnames
    assert 'salt' not in varnames
    
    # Use a common dimension for every 3D variable in the 'varnames' list.
    dims = (len(hr),medsea.M,medsea.N)

    # Define how to read from each grid point.
    def read_each(varname,ds):
        assert isinstance(ds,Dataset)
        nonlocal err_count, m, n, fullfn 
        return masked_invalid(ds[varname][hr,0,0])
        
    # Reading in from scattered netCDF3 files.
    grid_size = medsea.sea_m.size
    tic = time()
    err_count = 0    
    print('Initializing {:d} masked arrays of dimension {!s}...'.format(len(varnames),dims))
    data_dict = {varname: masked_all(dims) for varname in varnames}
    print('Elapsed {:.2f} s.'.format(time()-tic))
    
    print('Reading {:s} for these variables for {:d} grid points, serially:'.format(fn,grid_size))
    print(varnames)
    
    tic = time()
    for i in range(grid_size):
        m,n = medsea.sea_m[i], medsea.sea_n[i]
        fullfn = join(medsea.get_local_folder(m,n),fn)
        try:
            with Dataset(fullfn,'r') as ds:
                for varname in varnames:
                    data_dict[varname][:,m,n] = read_each(varname,ds)
        except Exception as e:
            if failed_list is not None:
                if (m,n) not in failed_list:
                    failed_list.append((m,n))
                    err_count += 1
                print('Error #{:d} occurred at (m,n) = ({:d},{:d}), for the file {:s}: \n {!s} \n'.format(err_count,m,n,fullfn,e))
                
        elapsed = time()-tic
        togo = elapsed / (i+1) * (grid_size-i)
        progress_msg = 'Elapsed {:02d}:{:02d}:{:02d} /// {:d}/{:d} grid points processed. /// About {:02d}:{:02d}:{:02d} to go... ' 
        print(progress_msg.format(
            int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60),
            i,grid_size,
            int(togo/3600),int(togo%3600/60),int(togo%60)),end='\r')

    elapsed = time()-tic
    print(' '*len(progress_msg),end='\r')
    print('Elapsed {:02d}:{:02d}:{:02d}.'.format(int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60)))
    
    if len(varnames) == 1:
        return data_dict[varnames[0]]
    else:
        return data_dict

def write(data_dict,varnames,year,month=None,fn=None,hr=None,epoch=None,grid_fn='grid_75m.dat'):
    from time import time
    import os
    from os.path import join
    from netCDF4 import Dataset, num2date

    assert not isinstance(varnames,str)
    assert 'temp' not in varnames
    assert 'salt' not in varnames
    
    if fn is None:        
        dr = join(medsea.project_folder,'results','{0.run}_{0.grid}'.format(medsea))
        if not(os.path.isdir(dr)):
            os.mkdir(dr)
            
        # Need to agree with that from outfn()    
        if month is None:
            fn = join(dr,'medsea_GOTM_{:s}_{:s}_{:d}.nc'.format(medsea.run,medsea.grid,year))
        else:
            fn = join(dr,'medsea_GOTM_{:s}_{:s}_{:d}{:02d}.nc'.format(medsea.run,medsea.grid,year,month))
    tic = time()
    print('Writing to {:s}...'.format(fn))

    def write_each(varname,data):
        tic = time()
        ncvar = medsea.create_variable(ds,varname,'f4',dimensions=('time','lat','lon'))
        # Actual write                    
        ncvar[:] = data[:]
        return time()-tic
        
    with Dataset(fn,'w',format='NETCDF3_CLASSIC') as ds:
        nctime, nclat, nclon = medsea.create_dimensions(ds,lat=medsea.grid_lats,lon=medsea.grid_lons)
        # Setting time units only, write the time variable after all data is written.
        if epoch is None:
            nctime.units = 'hours since {:d}-01-01 00:00:00'.format(year)
            #nctime[:] = hour_range(year,month=month)
        else:
            nctime.units = 'hours since {!s}'.format(epoch)
            #nctime[:] = hour_range(year,month=month,start=num2date(0,'hours since {:s}'.format(epoch)))

        # Use the override hr value if supplied, else resort to a default.
        if hr is None:
            hr = hour_range(year,month=month,start=epoch)

        # Actual write
        for varname in varnames:
            print('Begin writing {:s}...'.format(varname))
            assert data_dict[varname].shape[0] == len(hr) # check the data lengths 
            elapsed = write_each(varname,data_dict[varname])
            print('Elapsed {:d} min {:d} sec'.format(int(elapsed/60),int(elapsed%60)))

        # Python indices begin 0, 1, 2 but actual GOTM hours are 1, 2, 3, ...
        true_hr = [hr[i]+1 for i in range(len(hr))]
        nctime[:] = true_hr 
        
    print('Finished writing to {:s}.'.format(fn)) 
    elapsed = time()-tic
    print('Elapsed {:d} min {:d} sec'.format(int(elapsed/60),int(elapsed%60)))
