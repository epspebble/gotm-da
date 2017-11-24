from pyGOTM import medsea
from time import time
from os.path import join
from netCDF4 import Dataset
from numpy.ma import masked_all, masked_invalid
    
# def reload(grid='144x',run='ASM3-75m',nlev=122):
#     from importlib import reload
#     reload(medsea)
    
#     medsea.set_grid(grid)
#     medsea.run = run
#     medsea.nlev = nlev
    
# def hour_range(year,month=None,start=None,avoid_year_end=True):
#     """
#     Returns the GOTM time index range for the given month of year, 
#     i.e. the number of hours since the beginning to the simulation to the beginning of the month / beginning of the next month.
#     """ 
#     from datetime import date, datetime

#     # Assume the epoch is the beginning of year.
#     if start is None:
#         t0 = date(year,1,1)
#     else:
#         start_date = date(start.year,start.month,start.day) # should work for both date, datetime instances 
#         if isinstance(start,datetime):
#             assert start.hour + start.minute + start.second == 0, "We always assume GOTM begins at UTC midnight"
#         t0 = start_date
        
#     if month is None: # A year-long record.
#         t1 = date(year,1,1)
#         t2 = date(year,12,31) if avoid_year_end else date(year+1,1,1)
#     else: # A monthly record
#         t1 = date(year,month,1)
#         if month == 12:
#             t2 = date(year,12,31) if avoid_year_end else date(year+1,1,1)
#         else:
#             t2 = date(year,month+1,1)
        
#     offset = (t1-t0).days*24
#     nrec = (t2-t1).days*24
#     return range(offset,offset+nrec)

# # Checked this to agree with the definition of 'suffix' in medsea.local_run() as of 20171118
# def outfn(year=None, month=None, start=None, stop=None):
#     #         if month is None:
#     #             tag += '-{:d}'.format(year)
#     #         else:
#     #             tag += '-{:d}{:02d}'.format(year,month)
#     from datetime import date
#     if year is not None:
#         if month is None:
#             return 'results_{0.run}_{1:d}.nc'.format(medsea,year)
#         else:
#             return 'results_{0.run}_{1:d}{2:02d}.nc'.format(medsea,year,month)
#     else:
#         assert month is None
#         assert start is not None
#         assert stop is not None
        
#         ndays1 = (start - date(start.year,1,1)).days+1
#         ndays2 = (stop - date(stop.year,1,1)).days+1
#         return 'results_{:d}{:03d}{:d}{:03d}_{.run}.nc'.format(start.year,ndays1,stop.year,ndays2,medsea)    

def read(varnames,fn=None,year=None,month=None,use_ipp=False,failed_list=None):
   
    # Handle arguments
    if fn is None:
        #fn = outfn(year,month)
        fn = 'results_{0.run}_{0.GOTM_version}.nc'.format(medsea)
           
    # # If not provided, compute 'hr' from one of the nc files.
    # if hr is None:
    #     # Use the first readable file's metadata to determine start time, assuming GOTM's unit for 'time' in the ncdf file is the start time.        
    #     i=0
    #     start = None
    #     from datetime import date
    #     from netCDF4 import num2date
    #     while start is None:
    #         m,n = medsea.sea_m[i], medsea.sea_n[i]
    #         try:
    #             with Dataset(join(medsea.get_local_folder(m,n),fn),'r') as ds:
    #                 start = num2date(0,ds['time'].units).date()
    #         except:
    #             if i == medsea.sea_m.size-1:
    #                 raise OSError('Result file {!s} not found in all grid point folders for {:s}.'.format(fn,medsea.grid))
    #             else:
    #                 i += 1
    #                 pass
    #     hr = hour_range(year,month,start=start) 


    
    
    # Use medsea global settings and 'hr' to determine array dimensions.
    if isinstance(varnames,str):
        varname = varnames
        if varname == 'temp' or varname == 'salt':
            dims = (len(hr),medsea.nlev,medsea.M,medsea.N)
        else:
            dims = (len(hr),medsea.M,medsea.N)
        varnames = [varname]
    else:
        # Do not allow these two 4D variables to be lumped in one file.
        assert 'temp' not in varnames
        assert 'salt' not in varnames
        
        # Common dimensions for every other 3D variables.
        dims = (len(hr),medsea.M,medsea.N)

    # Allocate memory for a temporary array.
    def tmp_array():
        tic = time()
        print('Initializing a masked array of shape {!s} for temporary RAM storage...'.format(dims))
        tmp = masked_all(dims)
        print('Elapsed: {!s}'.format(time()-tic))
        return tmp

    # Define how to read from each grid point.
    def read_each(i):
        nonlocal err_count
        m,n = medsea.sea_m[i], medsea.sea_n[i]
        fullfn = join(medsea.get_local_folder(m,n),fn) 
        try:
            with Dataset(fullfn,'r') as ds:
                if len(dims) == 3:
                    return masked_invalid(ds[varname][hr,0,0])
                elif len(dims) == 4:
                    # Reversing the z dimension to agree with "depth"
                    return masked_invalid(ds[varname][hr,::-1,0,0])
                else:
                    assert False, "Problematic array dimensions: {!s}.".format(dims)
        except Exception as e:
            if failed_list is not None:
                if (m,n) not in failed_list:
                    failed_list.append((m,n))
                    err_count += 1
                    print('Error #{:d} occurred at (m,n) = ({:d},{:d}), for the file {:s}: \n {!s} \n'.format(err_count,m,n,fullfn,e))
            
            return masked_all(dims[:-2]) # Returning a masked array without the last two dimensions (lat, lon).
        
    # Reading in from scattered netCDF3 files.
    grid_size = medsea.sea_m.size
    
    data_dict = dict()
    for varname in varnames:
        tic = time()
        err_count = 0
        data = tmp_array()

        if use_ipp:
            print('Reading {:s} for {:s} data, {:d} grid points, parallelly...'.format(fn,varname,grid_size))
            from ipyparallel import Client
            rc = Client()
            lv = rc.load_balanced_view()
            for i, data_each in enumerate(lv.map(read_each,range(grid_size))):
                m,n = medsea.sea_m[i], medsea.sea_n[i]
                if len(dims) == 3:
                    data[:,m,n] = data_each
                elif len(dims) == 4:
                    data[:,:,m,n] = data_each
                else:
                    assert False
                    elapsed = time()-tic
                    print('Elapsed {:02d}:{:02d}:{:02d}.'.format(int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60)))

        else:
            print('Reading {:s} for {:s} data, {:d} grid points, serially...'.format(fn,varname,grid_size))
            for i in range(grid_size):
                m,n = medsea.sea_m[i], medsea.sea_n[i]
                if len(dims) == 3:
                    data[:,m,n] = read_each(i)
                elif len(dims) == 4:
                    data[:,:,m,n] = read_each(i)
                else:
                    assert False
                    
                elapsed = time()-tic
                togo = elapsed / (i+1) * (grid_size-i)
                progress_msg = 'Elapsed {:02d}:{:02d}:{:02d} /// {:d}/{:d} grid points processed. /// About {:02d}:{:02d}:{:02d} to go... ' 
                print(progress_msg.format(
                    int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60),
                    i+1,grid_size,
                    int(togo/3600),int(togo%3600/60),int(togo%60)),end='\r')
        data_dict[varname] = data.copy()

    elapsed = time()-tic
    print(' '*len(progress_msg),end='\r')
    print('Elapsed {:02d}:{:02d}:{:02d}'.format(int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60)))
    
    if len(varnames) == 1:
        return data_dict[varnames[0]]
    else:
        return data_dict

def write(data,varname,year,month=None,hr=None,fn=None,epoch=None,grid_fn='grid_75m.dat'):
    from time import time
    from os.path import join
    from netCDF4 import Dataset, num2date
    import os

    if fn is None:
        dr = join(medsea.project_folder,'results','{0.run}_{0.grid}'.format(medsea))
        if not os.path.isdir(dr):
            os.mkdir(dr)
            
        # Temporary add 75m after grid name e.g. '144x', we should later update the grid name to just '144x75m'
        if month is None:
            fn = join(dr,'medsea_GOTM_{:s}_{:s}_{:s}_{:d}.nc'.format(varname,
                                                                     medsea.run,
                                                                     medsea.grid,
                                                                     year))
        else:
            fn = join(dr,'medsea_GOTM_{:s}_{:s}_{:s}_{:d}{:02d}.nc'.format(varname,
                                                                           medsea.run,
                                                                           medsea.grid,
                                                                           year,month))
    tic = time()
    print('Writing to {:s}...'.format(fn))
    with Dataset(fn,'w',format='NETCDF3_CLASSIC') as ds:
        nctime, nclat, nclon = medsea.create_dimensions(ds,lat=medsea.grid_lats,lon=medsea.grid_lons)
        if varname == 'temp' or varname == 'salt':
            ds.createDimension('depth', size = medsea.nlev)
            depth = ds.createVariable('depth', 'f4', dimensions=('depth',))
            from numpy import loadtxt, cumsum
            depth[:] = cumsum(loadtxt(join(medsea.project_folder,'nml',grid_fn),skiprows=1))
            ncvar = medsea.create_variable(ds,varname,'f4',dimensions=('time','depth','lat','lon'))
        else:
            ncvar = medsea.create_variable(ds,varname,'f4',dimensions=('time','lat','lon'))

        # Setting the time units only
        if epoch is None:
            nctime.units = 'hours since {:d}-01-01 00:00:00'.format(year)
            #nctime[:] = hour_range(year,month) # A "None" could be passed to hour_range.
        else:
            nctime.units = 'hours since {!s}'.format(epoch)
            #nctime[:] = hour_range(year,month=month,start=num2date(0,'hours since {:s}'.format(epoch)))

        # Calculate time variable.
        if hr is None:
            hr = hour_range(year,month=month,start=epoch)

        assert data.shape[0] == len(hr), "Number of hours specified is NOT equal to the time dimension of {!s}!".format(varname)

        # Python indices begin 0, 1, 2 but actual GOTM hours are 1, 2, 3, ...
        true_hr = [hr[i]+1 for i in range(len(hr))]
        nctime[:] = true_hr 
        # Actual write
        ncvar[:] = data[:]
        
    elapsed = time()-tic
    print('Elapsed {:d} min {:d} sec'.format(int(elapsed/60),int(elapsed%60)))
