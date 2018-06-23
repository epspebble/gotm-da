#!/home/simontse/miniconda3/bin/python
#SBATCH --requeue
#SBATCH --time=0-00:05
#SBATCH --mem=8G
#SBATCH --mail-user=simon.tse@twu.ca
#SBATCH --mail-type=END
#SBATCH --profile=all
##SBATCH --mail-type=BEGIN,END
##SBATCH -e %x.e%A-%a
##SBATCH -o %x.o%A-%a

from time import time
import os, sys
from os.path import join
from netCDF4 import Dataset, MFDataset, num2date, date2num
from numpy import loadtxt
from numpy.ma import masked_all, masked_invalid, masked_equal
from datetime import date, datetime, timedelta

def print_usage():
    print('Usage: ' + sys.argv[0] + ' <grid> <run> <GOTM_version> <year> [month] [--biannual={True|False}]')
    sys.exit(0)

def get_args():
    if len(sys.argv) < 5:
        print_usage()
    else:
        grid,run,ver,year = sys.argv[1:5]
        year = int(year)
        
    if len(sys.argv) >= 6:
        month = int(sys.argv[5])
        unparsed_kwargs = sys.argv[6:]
    else:
        month = None
        unparsed_kwargs = None

    # Maybe don't reinvent the wheel... see: https://docs.python.org/3/library/argparse.html
    kwargs_syntax = dict(biannual=bool)
    kwargs = dict()
    if unparsed_kwargs:
        for each in unparsed_kwargs:
            sep = each.find('=') # returns -1 if not found
            if each[:2] != '--' or sep == -1: 
                print('Invalid argument: ' + each)
                print_usage()
            else:
                key = each[:sep]
                val = each[sep+1:]
                try:
                    val_type = kwargs_syntax(key)
                    kwargs[key] = val_type(val)
                except Exception as e:
                    print(e)
                    print_usage()
    return grid, run, ver, year, month, kwargs

### Getting a configured medsea module instance.
# This was a portion of original main program brought forward.
if __name__=='__main__':
    grid, run, ver, year, month, kwargs = get_args()
    # from pyGOTM import config
    # config.GOTM_version = ver

    #20180622 After a short while, I've forgotten why I had "Old method" vs "New method" and
    # There're some vestiges in pyGOTM.medsea for which the purpose is not clear. Before a
    # refactoring happens, maybe it's just easier to pass all arguments every time rather than
    # to set some "global settings". It might be more useful in an interactive environment, but
    # as a script to do a fixed purpose, this is not really helping.
    from pyGOTM import medsea
    medsea.set_grid(grid) # Still needed to set medsea.grid_lats, .grid_lons, .M, .N, .sea_mn, etc convenience grid related values.

    # Old method:
    #medsea.set_grid(grid)
    #medsea.run = run

    # New method:
    #assert run[:3] == 'ASM'
    #medsea.set_grid(grid,run[3:]) # the string skipping to leading characters 'ASM'

### End setting the medsea module instance into this global namespace.

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

def read_results(varnames,year,
#                 medsea, # temporary hack to pass on configs
                 month=None,skip_year_end=False): 
    fn = 'results_{:s}_{:s}.nc'.format(run,ver)
    
    hours = hour_range(year,month=month,skip_year_end=skip_year_end)
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
                    raise ValueError("Time out of range. First record found is at {!s}".format(num2date(ds['time'][0],ds['time'].units)))
                if last_i > ds['time'].size:
                    raise ValueError("Time out of range. Last record found is at {!s}".format(num2date(ds['time'][-1],ds['time'].units)))
                if len(dims) == 3:
                    data = ds[varname][first_i:last_i,0,0]
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

def write_results(data,year,month=None,varname=None,fn=None,skip_year_end=False):
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

    dr = join(medsea.results_folder,run + '_' + grid)

    if not os.path.isdir(dr):
        os.mkdir(dr)
        
    if fn is None:
        #1. product name
        fn = 'medsea_GOTM'
        #2. variable name
        if varname is not None:
            fn += '_' + varname        
        #3. codenames: run profile, grid, GOTM_version
        fn += '_{:s}_{:s}_{:s}'.format(run,grid,ver)
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
        hours = hour_range(year,month=month,skip_year_end=skip_year_end)
        nctime[:] = date2num(hours,nctime.units)
    print('Elapsed {!s}s.'.format(time()-tic))
    
    for varname in varnames:
        tic = time()
        print('Writing data for {:s}...'.format(varname))
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

def combine_stat(*args,biannual=False):
    """
    combine_stat(year) or
    combine_stat(year,month) or
    combine_stat(start,stop) where start, stop are datetime.date objects

    Version as of 2017-11-28.
    """
    # Checking arguments.
    if len(args) == 2:
        if isinstance(args[0],date) and isinstance(args[1],date):
            start = args[0]
            stop = args[1]
            assert stop > start
            year = None
            month = None
        else:
            year, month = args
            assert isinstance(year,int) and year > 2010
            if month is None:
                # Same as with len(args) == 1, i.e. only year is supplied.
                start = date(year,1,1)
                stop = date(year+1,1,1)
            else:
                assert isinstance(month,int)
                start = date(year,month,1)
                stop = date(year+1,1,1) if month == 12 else date(year,month+1,1)
    elif len(args) == 1:
        year = args[0]
        month = None
        assert isinstance(year,int) and year > 2010
        start = date(year,1,1)
        stop = date(year+1,1,1)
        
    dat_fn = 'daily_stat_{:s}_{:s}.txt'.format(run,ver)
               
    # Done defining start, stop.
    print('Working on {!s} to {!s}...\n\n'.format(start, stop))  
    print('Daily statistics file: ' + dat_fn)
    
    outdir = os.path.join(medsea.results_folder,run + '_' + grid)
    if not(os.path.isdir(outdir)):
          os.mkdir(outdir)

    basename = 'daily_stat_{:s}_{:s}_{:s}_'.format(run,grid,ver)
    if year is None: # A start-stop based time range:
        basename += '{0.year:d}{0.month:02d}{0.day:02d}'.format(start) + '_' + \
                    '{0.year:d}{0.month:02d}{0.day:02d}'.format(stop) + '.nc'
    else:
        if month is None:
            basename += '{:d}.nc'.format(year)
        else:
            basename += '{:d}{:02d}.nc'.format(year,month)
            
    outfn = os.path.join(outdir,basename)

    print('Writing dimensions and metadata...')
    elapsed = 0
    medsea.tic()

    with Dataset(outfn,'w',format='NETCDF3_CLASSIC') as ds:

        # Dimensions
        ds.createDimension('day')
        ds.createDimension('lat', size = len(medsea.grid_lats))
        ds.createDimension('lon', size = len(medsea.grid_lons))
        
        # Dimension variables.
        if year is None: 
            day_var = ds.createVariable('day','i4',dimensions=('day',))
            day_var.units = "days since " + str(start-timedelta(days=1))
        else:
            day_var = ds.createVariable('day','i4',dimensions=('day',))
            day_var.units = "days since " + str(date(year-1,12,31))
            
        lat_var = ds.createVariable('lat','f4',dimensions=('lat',))
        lat_var.units = 'degrees north'
        lon_var = ds.createVariable('lon','f4',dimensions=('lon',))
        lon_var.units = 'degrees east'
        lat_var[:] = medsea.grid_lats
        lon_var[:] = medsea.grid_lons
        print('Done initializing dimensions.')
        
        data = [None for i in range(7)]
        data[0]=ds.createVariable('assim_time','i4',dimensions=('day','lat','lon'))
        data[0].units = 'seconds since midnight in local time'
        data[4]=ds.createVariable('SST_max','f4',dimensions=('day','lat','lon'))
        data[4].units = 'degree Celsius'
        data[2]=ds.createVariable('SST_min_day','f4',dimensions=('day','lat','lon'))
        data[2].units = 'degree Celsius'
        data[6]=ds.createVariable('SST_min_night','f4',dimensions=('day','lat','lon'))
        data[6].units = 'degree Celsius'
        data[3]=ds.createVariable('SST_max_time','i4',dimensions=('day','lat','lon'))
        data[3].units = 'seconds since midnight in local time'
        data[1]=ds.createVariable('SST_min_day_time','i4',dimensions=('day','lat','lon'))
        data[1].units = 'seconds since midnight in local time'
        data[5]=ds.createVariable('SST_min_night_time','i4',dimensions=('day','lat','lon'))
        data[5].units = 'seconds since midnight in local time'
        print('Done creating variables.')

        yrdays = 366 if year%4==0 else 365

        # temporary fix for binannual 2013-2014 run.
        ndays = (stop-start).days
        
        start_daynum = (start - date(start.year-1,12,31)).days
        # The 12-th month stops at the start of next year
        stop_daynum = yrdays+1 if month == 12 else (stop - date(stop.year-1,12,31)).days
        # print(start_daynum,stop_daynum,ndays)
        #assert ndays == stop_daynum - start_daynum, "Might have crossed year boundary, don't know what to do..."
        date_range = num2date(range(ndays), 'days since ' + str(start))            
        for i, each in enumerate(date_range):
            day_var[i] = date2num(each,day_var.units)
        
        # redefine ndays and day_var
        temp = [masked_all((ndays,medsea.M,medsea.N)) for i in range(7)]
        
        medsea.tic()
        ndfv = [0,0,99.,0,-99.,0,99.] # fill values
        err_count = 0
        for i in range(medsea.sea_m.size):
            m = medsea.sea_m[i]
            n = medsea.sea_n[i]
            latlong = medsea.print_lat_lon(*medsea.get_lat_lon(m,n))
                
            infn = os.path.join(medsea.get_local_folder(m,n),dat_fn)
            if not(os.path.isfile(infn)):
                err_count += 1
                print('\nFile not found: ' + infn)
                continue
            try:
                tmp = loadtxt(infn)
            except ValueError as ve :
                print('\n Error reading {:s}: ValueError: {!s}'.format(infn,ve))
                err_count += 1
                continue
                
            if tmp.size == 0: # empty file
                print('\nEmpty file: ' + infn)
                err_count +=1
                continue

            day_of_year = tmp[:,0]                
            ind = (day_of_year >= start_daynum) & (day_of_year < stop_daynum)
            for k in range(7):
                #if dat_fn[:25] == 'daily_stat_20130012014365':
                
                ## Special code to combine a biannual run for 2013-2014.
                if biannual:
                    if tmp.shape[0] != 729:
                        err_count += 1
                        print('\nUnexpected line count: {:d} for the location {!s}'.format(tmp.shape[0],latlong))
                        break
                    if year == 2013:
                        temp[k][:364,m,n] = tmp[:364,k+1]
                    elif year == 2014:
                        temp[k][:364,m,n] = tmp[365:,k+1]
                else:
                    try:
                        temp[k][:,m,n] = tmp[ind,k+1]
                    except Exception as e:
                        err_count +=1
                        print(e)
                        raise
            # Beware, since the progress message overwrites itself, error messages need to do a newline first.
            print('Reading {0:d}/{1:d} grid points, {2:d}/{1:d} bypassed.'.format(i+1,medsea.sea_m.size,err_count),end='\r')
        print('\nDone reading in all data.')
        elapsed += medsea.toc()
        #save('daily_stats',temp)

        for i in range(7):
            medsea.tic()
            data[i][:] = masked_equal(temp[i],ndfv[i])
            #print("Done writing a variable.")
            #elapsed += medsea.toc()
            
        print('Done writing out data to {:s}.'.format(outfn))
        elapsed += medsea.toc()

# If called as a script.        
if __name__=='__main__':
    # Main program separated into two parts, just to initialize pyGOTM.medsea
    # earlier on.

    # Temporary code.
    if year in [2013, 2014]:
        biannual = True
        skip_year_end = False if year == 2013 else True
    else:
        biannual = False
        skip_year_end = False
    
    varnames = ['sst','skint','x-taus','y-taus','swr','heat','total','lwr','sens','latent','albedo','coszen']
    data = read_results(varnames,year,month=month,skip_year_end=skip_year_end)
    write_results(data,year,month=month,skip_year_end=skip_year_end)
    del data
    
    temp = read_results('temp',year,month=month,skip_year_end=skip_year_end)
    write_results(temp,year,month=month,varname='temp',skip_year_end=skip_year_end)
    del temp
    
    salt = read_results('salt',year,month=month,skip_year_end=skip_year_end)
    write_results(salt,year,month=month,varname='salt',skip_year_end=skip_year_end)
    del salt

    combine_stat(year,month,biannual=biannual)
    
