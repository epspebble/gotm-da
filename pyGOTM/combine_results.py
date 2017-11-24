#!/home/simontse/miniconda3/bin/python
from time import time
import os,sys
from os.path import join
from netCDF4 import Dataset, MFDataset, num2date, date2num
from numpy.ma import masked_all, masked_invalid
from datetime import date, datetime, timedelta

def hour_range(year,month=None,skip_year_end=False):
    if month is None:
        first = datetime(year,1,1,1,0,0)
        last = datetime(year+1,1,1,0,0,0)
        if skip_year_end:
            last -= timedelta(days=1)
    else:
        first = datetime(year,month,1,1,0,0)
        if month == 12:
            last = datetime(year+1,1,1,0,0,0)
            if skip_year_end:
                last -= timedelta(days=1)
        else:
            last = datetime(year,month+1,1,0,0,0)

    nrec = (last.date()-first.date()).days*24 
    hours = [first + timedelta(hours=i) for i in range(nrec)]
    return hours

def read_results(varnames,year,month=None,skip_year_end=False):
    fn = 'results_{0.run}_{0.GOTM_version}.nc'.format(medsea)
    
    hours = hour_range(year,month=month,skip_year_end=False)
    first, last = hours[0], hours[-1]
    nrec = len(hours)
    
    if isinstance(varnames,str):
        varname = varnames
        if varname == 'temp' or varname == 'salt':
            dims = (nrec,medsea.nlev,medsea.M,medsea.N)
        else:
            dims = (nrec,medsea.M,medsea.N)
        varnames = [varname]
    else:
        # Do not allow these two 4D variables to be lumped in one file.
        assert 'temp' not in varnames
        assert 'salt' not in varnames
        
        # Common dimensions for every other 3D variables.
        dims = (nrec,medsea.M,medsea.N)

    # Allocate memory for a temporary array.
    def tmp_array():
        tic = time()
        print('Initializing a masked array of shape {!s} for temporary RAM storage...'.format(dims))
        tmp = masked_all(dims)
        print('Elapsed: {!s}'.format(time()-tic))
        return tmp

    # Define how to read from each grid point.
    failed_list = list()
    def read_each(i):
        nonlocal err_count
        m,n = medsea.sea_m[i], medsea.sea_n[i]
        fullfn = join(medsea.get_local_folder(m,n),fn) 
        try:
            with Dataset(fullfn,'r') as ds:
                # Compute indices, assuming GOTM gives continuously hourly output begininng at 1, 2, 3 hours...
                first_i = int((date2num(first,ds['time'].units) - ds['time'][0])/3600)
                last_i = int((date2num(last,ds['time'].units) - ds['time'][0])/3600)+1
                if first_i < 0:
#                     print(first,first_i,ds['time'][0],ds['time'].units)
                    raise ValueError("Time out of range. First record found is at {!s}".format(num2date(ds['time'][0],ds['time'].units)))
                if last_i >= ds['time'].size:
#                     print(last,last_i,ds['time'][-1],ds['time'].units)
                    raise ValueError("Time out of range. Last record found is at {!s}".format(num2date(ds['time'][-1],ds['time'].units)))
                if len(dims) == 3:
                    data = ds[varname][first_i:last_i,0,0]
#                     print(first_i,last_i,data.shape)
                    return masked_invalid(data)
                elif len(dims) == 4:
                    # Reversing the z dimension to agree with "depth"
                    return masked_invalid(ds[varname][first_i:last_i,::-1,0,0])
                else:
                    raise ValueError("Problematic array dimensions: {!s}.".format(dims))
        except:
            raise
           
            return masked_all(dims[:-2]) # Returning a masked array without the last two dimensions (lat, lon).
        
    # Reading in from scattered netCDF3 files.
    grid_size = medsea.sea_m.size
    
    data_dict = dict()
    for varname in varnames:
        tic = time()
        err_count = 0
        data = tmp_array()

        print('Reading {:s} for {:s} data, {:d} grid points, serially...'.format(fn,varname,grid_size))
        for i in range(grid_size):
            m,n = medsea.sea_m[i], medsea.sea_n[i]
            if len(dims) == 3:
                data[:,m,n] = read_each(i)
            elif len(dims) == 4:
                data[:,:,m,n] = read_each(i)
            else:
                raise ValueError("Problematic array dimensions: {!s}.".format(dims))
                    
            elapsed = time()-tic
            togo = elapsed / (i+1) * (grid_size-i)
            progress_msg = 'Elapsed {:02d}:{:02d}:{:02d} /// {:d}/{:d} grid points processed. /// About {:02d}:{:02d}:{:02d} to go... ' 
            print(progress_msg.format(int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60),
                    i+1,grid_size, int(togo/3600),int(togo%3600/60),int(togo%60)),end='\r')
        data_dict[varname] = data.copy()

    elapsed = time()-tic
    print(' '*len(progress_msg),end='\r')
    print('Elapsed {:02d}:{:02d}:{:02d}'.format(int(elapsed/3600),int(elapsed%3600/60),int(elapsed%60)))
    
    if len(varnames) == 1:
        return data_dict[varnames[0]]
    else:
        return data_dict

def write_results(data,year,month=None,varname=None,fn=None):
    from numpy import loadtxt, cumsum
    epoch = medsea.epoch
    grid_fn = 'grid_75m.dat' # Temporary solution.
    
    if isinstance(data,dict):
        varnames = list(data.keys())
    elif varname is not None:
        varnames = [varname]
        data = {varname: data}
    else:
        raise ValueError("Check 'data' and 'varname'")

    dr = join(medsea.project_folder,'results','{0.run}_{0.grid}'.format(medsea))
    if not os.path.isdir(dr):
        os.mkdir(dr)
        
    if fn is None:
        #1. product name
        fn = 'medsea_GOTM'
        #2. variable name
        if varname is not None:
            fn += '_' + varname        
        #3. run code and grid code
        fn += '_' + medsea.run + '_' + medsea.grid
        #4. time range
        if month is None:
            fn += '_{:d}.nc'.format(year)
        else:
            fn += '_{:d}{:02d}.nc'.format(year,month)

    fullfn = join(dr,fn)
    
    tic = time()
    begin = tic
    print('Writing dimensions data to {:s}...'.format(fullfn))
    with Dataset(fullfn,'w',format='NETCDF3_CLASSIC') as ds:
        nctime, nclat, nclon = medsea.create_dimensions(ds,lat=medsea.grid_lats,lon=medsea.grid_lons)
        if varname == 'temp' or varname == 'salt':
            ds.createDimension('depth', size = medsea.nlev)
            depth = ds.createVariable('depth', 'f4', dimensions=('depth',))
            depth[:] = cumsum(loadtxt(join(medsea.project_folder,'nml',grid_fn),skiprows=1))

        # Setting the time units only
        if epoch is None:
            nctime.units = 'hours since {:d}-01-01 00:00:00'.format(year)
            #nctime[:] = hour_range(year,month) # A "None" could be passed to hour_range.
        else:
            nctime.units = 'hours since {!s}'.format(epoch)
            #nctime[:] = hour_range(year,month=month,start=num2date(0,'hours since {:s}'.format(epoch)))
        hours = hour_range(year,month=month)
        nctime[:] = date2num(hours,nctime.units)
    print('Elapsed {!s}s.'.format(time()-tic))
    
    for varname in varnames:
        tic = time()
        print('Writng data for {:s}...'.format(varname))
        if data[varname].shape[0] != len(hours):
            raise ValueError("Number of hours specified is NOT equal to the time dimension of {!s}!".format(varname))

        # Actual write
        with Dataset(fullfn,'a',format='NETCDF3_CLASSIC') as ds:
            if varname == 'temp' or varname == 'salt':            
                ncvar = medsea.create_variable(ds,varname,'f4',dimensions=('time','depth','lat','lon'))
            else:
                ncvar = medsea.create_variable(ds,varname,'f4',dimensions=('time','lat','lon'))
            ncvar[:] = data[varname][:]
        elapsed = time()-tic
        print('Elapsed {:d} min {:d} sec'.format(int(elapsed/60),int(elapsed%60)))    
        
    total_elapsed = time()-begin
    print('\n Total time: {:d} min {:d} sec'.format(int(total_elapsed/60),int(total_elapsed%60)))
    
if __name__=='__main__':
    
    grid,run,ver,year,month = sys.argv[1:]
    year = int(year)
    month = int(month)
    from pyGOTM import medsea
    medsea.set_grid(grid)
    medsea.run = run
    medsea.ver = ver
    varnames = ['sst','skint','x-taus','y-taus','swr','heat','total','lwr','sens','latent','albedo']
    data = read_results(varnames,year,month=month)
    write_results(data,year,month=month)
    del data
    
    temp = read_results('temp',year,month=month)
    write_results(temp,year,month=month,varname='temp')
    del temp
    
    salt = read_results('salt',year,month=month)
    write_results(salt,year,month=month,varname='salt')
    del salt
    
