def run_setup(grid, ver, wd):
    # Configuring modules.
    from pyGOTM import config
    config.setup(GOTM_version=ver)

    from pyGOTM import medsea
#    assert run[:3] == 'ASM'
#    assert run[-4:] == '-75m'
#    medsea.set_grid(grid,run[3:-4])
    medsea.set_grid(grid)
    #medsea.run = run
    assert medsea.GOTM_version == ver
    medsea.set_folders(new_grid_folder=wd)
    medsea.verbose = False
    
    # grid_folder, latlong = os.path.split(fd)
    # assert grid_folder == medsea.grid_folder
    # # Use separate cache folders to avoid crashing.
    # medsea.cache_folder = os.path.join('/dev/shm','medsea_GOTM_' + latlong)
    # if not(os.path.isdir(medsea.cache_folder)):
    #     os.mkdir(medsea.cache_folder)
    
    return medsea

def write_dat(indices,run_instance):
    """    
    Usage:
    from pyGOTM import medsea # import a predefined run instance
    
    indices = medsea.sea_mn[0:100] # the predefined indices was stored as an attribute of the run instance, take only first 100 indices

    write_dat(indices,medsea) # passing the indices and the run instance to the write_date method, which should be implemented as a method
    of the run instance instead of a method of the pipeline... maybe

    Current design issue (2018-08-22), in realistic usage, the input data files may not be static and could come from different sources.

    For e.g. the chlo.dat file could come from either linear interpolations of 8-day average, 
    or a daily average with some sort of interpolation or reconstruction.

    Moreover, the same file might be used for more than one purpose: chlo.dat is used both as an input to the extinction (light absorption) 
    parametrization, and to the albedo parametrization. Opening the same file twice is acceptable if compiled by ifort, but not by gfortran
    as per Fortran standard.


    """
    from netCDF4 import Dataset, MFDataset, num2date
    from os.path import join
    from numpy.ma import masked_outside, is_masked
    import sys
    from os import getenv
    from time import sleep
    from time import time as toc
    from sys import argv

    # Create folders if necessary
    
    
    print('Generating input data files...')
    write_dat_begin = toc()
    #print('Creating folders if necessary...')
    for m,n in indices:
        #lat,lon = run_instance.get_lat_lon(m,n)
        dst_folder = run_instance.get_local_folder(m,n,create=True)
    #print('Elapsed {:.2f}s'.format(toc()-tic))
    
    # Read data for heat.dat
    tic = toc()
    print('Reading in data for heat.dat...')
    with run_instance.get_data_sources('heat') as ds:
        time = num2date(ds['time'][:],ds['time'].units)
        lwrd = ds['lwrd'][:]
        swrd = ds['swrd'][:]
        swrd_cs = ds['swrd_cs'][:]
        cloud_factor = masked_outside(ds['cloud_factor'][:],0,2)
    print('Elapsed {:.2f}s'.format(toc()-tic))
    
    # Write data to heat.dat
    tic = toc()
    print('Writing heat.dat to the grid folders...')
    for m,n in indices:
        lat,lon = run_instance.get_lat_lon(m,n)
        dst_folder = run_instance.get_local_folder(m,n,)
        dst_fn = 'heat.dat'
        with open(join(dst_folder,dst_fn),'w') as f:
            col = [None for i in range(4)]
            for i, t in enumerate(time):
                col[0] = str(t)
                col[1] = swrd[i,m,n]
                col[2] = 0 if is_masked(cloud_factor[i,m,n]) else cloud_factor[i,m,n] # Masked when nighttime.
                col[3] = lwrd[i,m,n]
                line = ('{:s}'+' {:10.5g}'*3 + '\n').format(*col)
                f.write(line)    
    print('Elapsed {:.2f}s'.format(toc()-tic))
    del lwrd, swrd, swrd_cs, cloud_factor

    # Read data for met.dat
    tic = toc()
    print('Reading in data for met.dat...')
    with run_instance.get_data_sources('met') as ds:
        time = num2date(ds['time'][:],ds['time'].units)
        u10m = ds['u10m'][:]
        v10m = ds['v10m'][:]
        sp = ds['sp'][:]
        t2m = ds['t2m'][:]
        q2m = ds['q2m'][:]
        precip = ds['precip'][:]    
    print('Elapsed {:.2f}s'.format(toc()-tic))

    # Write data for met.dat
    tic = toc()
    print('Writing met.dat to the grid folders...')
    for m,n in indices:
        #lat,lon = run_instance.get_lat_lon(m,n)
        dst_folder = run_instance.get_local_folder(m,n)
        dst_fn = 'met.dat'
        with open(join(dst_folder,dst_fn),'w') as f:
            col = [None for i in range(9)]
            for i, t in enumerate(time):
                col[0] = str(t)
                col[1] = u10m[i,m,n]
                col[2] = v10m[i,m,n]
                col[3] = sp[i,m,n]
                col[4] = t2m[i,m,n]
                col[5] = q2m[i,m,n]
                col[6] = 0 # cloud!?
                col[7] = precip[i,m,n]
                col[8] = 0 # snow!?
                line = ('{:s}'+' {:10.5g}'*8 + '\n').format(*col)
                f.write(line)    
    print('Elapsed {:.2f}s'.format(toc()-tic))
    del time, u10m, v10m, sp, t2m, q2m, precip

    # Read data for tprof.dat
    tic = toc()
    print('Reading in data for tprof.dat...')    
    with run_instance.get_data_sources('tprof') as ds: 
        depth = ds['depth'][:]
        dates = num2date(ds['time'][:],ds['time'].units)
        temp = masked_outside(ds['temp'][:],4,35)   
    print('Elapsed {:.2f}s'.format(toc()-tic))

    # Write data for tprof.dat
    tic = toc()
    print('Writing tprof.dat to the grid folders...')    
    for m,n in indices:
        #lat,lon = run_instance.get_lat_lon(m,n)
        dst_folder = run_instance.get_local_folder(m,n)
        dst_fn = 'tprof.dat'
        with open(join(dst_folder,dst_fn),'w') as f:
            col = [None for i in range(9)]
            for i in range(len(dates)):
                f.write(str(dates[i]) + ' {0:d} 2\n'.format(len(depth))) # Always two columns.
                for j in range(len(depth)):
                    line = ('{0:g} {1:g}\n').format(-depth[j],temp[i,j,m,n])
                    f.write(line)       
    print('Elapsed {:.2f}s'.format(toc()-tic))
    del depth, dates, temp

    # Read data for sprof.dat
    tic = toc()
    print('Reading in data for sprof.dat...')    
    with run_instance.get_data_sources('sprof') as ds: 
        depth = ds['depth'][:]
        dates = num2date(ds['time'][:],ds['time'].units)
        salt = masked_outside(ds['salt'][:],0,1000)    
    print('Elapsed {:.2f}s'.format(toc()-tic))

    # Write data for sprof.dat
    tic = toc()
    print('Writing sprof.dat to the grid folders...')       
    for m,n in indices:
        #lat,lon = run_instance.get_lat_lon(m,n)
        dst_folder = run_instance.get_local_folder(m,n)
        dst_fn = 'sprof.dat'
        with open(join(dst_folder,dst_fn),'w') as f:
            col = [None for i in range(9)]
            for each in dates:
                f.write(str(each) + ' {0:d} 2\n'.format(len(depth))) # Always two columns.
                for j in range(len(depth)):
                    line = ('{0:g} {1:g}\n').format(-depth[j],salt[i,j,m,n])
                    f.write(line)
    print('Elapsed {:.2f}s'.format(toc()-tic))
    del depth, dates, salt

    # Read data for chlo.dat
    tic = toc()
    print('Reading in data for chlo.dat...')    
    with run_instance.get_data_sources('chlo') as ds:
        time = num2date(ds['time'][:],ds['time'].units)
        chlo = masked_outside(ds['chlor_a'][:],0,10)
    print('Elapsed {:.2f}s'.format(toc()-tic))

    # Write data for chlo.dat    
    tic = toc()
    print('Writing chlo.dat to the grid folders...')       
    for m,n in indices:
        dst_folder = run_instance.get_local_folder(m,n)
        dst_fn = 'chlo.dat'
        with open(join(dst_folder,dst_fn),'w') as f:
            col = [None for i in range(2)]
            count = 0
            for i, t in enumerate(time):
                if is_masked(chlo[i,m,n]):
                    count +=1
                    continue
                else:
                    count = 0
                    line = '{:s} {:10.5g}\n'.format(str(t),chlo[i,m,n])
                    f.write(line)
                if count > 4:
                    print('WARNING: {:d} consecutive nan values while writing {:s}.'.format(count,join(dst_folder,dst_fn)))

    del chlo
    print('Elapsed {:.2f}s'.format(toc()-tic))
    
    # Read data for iop.dat
    tic = toc()
    print('Reading in data for iop.dat...')
    with run_instance.get_data_sources('iop') as ds:
        time = num2date(ds['time'][:],ds['time'].units)
        a = ds['a_488_giop'][:]
        bb = ds['bb_488_giop'][:]
    print('Elapsed {:.2f}s'.format(toc()-tic))

    # Write data for iop.dat    
    tic = toc()
    print('Writing iop.dat to the grid folders...')
    for m,n in indices:
        dst_folder = run_instance.get_local_folder(m,n)
        dst_fn = 'iop.dat'
        with open(join(dst_folder,dst_fn),'w') as f:
            col = [None for i in range(2)]
            count = 0
            for i, t in enumerate(time):
                if is_masked(a[i,m,n]) or is_masked(bb[i,m,n]):
                    count +=1
                    continue
                else:
                    count = 0
                    line = '{:s} {:10.5g} {:10.5g}\n'.format(str(t),a[i,m,n],bb[i,m,n])
                    f.write(line)
            if count > 4:
                print('WARNING: {:d} consecutive nan values while writing {:s}.'.format(count,join(dst_folder,dst_fn)))
    #print('Elapsed {:.2f}s'.format(toc()-tic))
    del a, bb
    print('Elapsed {:.2f}s'.format(toc()-tic))
    print('Total write_dat time: {:.2f}s'.format(toc()-write_dat_begin))
    
def create_header(fn,i0,i1,var3d,var4d=None,epoch='1981-01-01 00:00:00'):
    from netCDF4 import Dataset
    with Dataset(fn,'w',format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('time')
        ds.createDimension('index', size = i1-i0)
        if var4d is not None:
            ds.createDimension('depth', size = 122)
            ncdepth = ds.createVariable('depth','f4',dimensions=('depth',))
        
        nctime = ds.createVariable('time','i4',dimensions=('time',))
        nctime.units = 'hours since ' + epoch
        ncindex = ds.createVariable('index','i2',dimensions=('index',))
        ncindex = range(i0,i1)
        nclat = ds.createVariable('lat','f4',dimensions=('index',))
        nclat.units = 'degrees north'
        nclon = ds.createVariable('lon','f4',dimensions=('index',))
        nclon.units = 'degree east'

        for varname in var3d:
            ds.createVariable(varname,'f4',dimensions=('time','index'))
        for varname in var4d:
            ds.createVariable(varname,'f4',dimensions=('time','index','depth'))
            
def batch_run(grid, run, ver, start, stop, wd, i0=None,i1=None, c=None, k=None, prepared_input=True, cached_write=False):
    """
    grid = '144x'
    run = 'ASM3-75m'
    ver = 'v3g'
    start = '2013-01-01'
    stop = '2014-12-31'
    c = 100 % chunk size
    k = 0 % which chunk of indices? range(c*k,c*(k+1))
    wd = /dev/shm %working directory
    """
    from os.path import join, isfile, basename
    from time import time
    from math import log10
    import shutil
    from netCDF4 import Dataset
    
    job_start = time()
    
    # Setup the run instance.
    init_start = time()
    print('Setting up batch GOTM-{:s} run for {:s} grid in {:s}'.format(ver,grid,wd))
    medsea = run_setup(grid, ver, wd)

    assert medsea.grid == grid
    assert medsea.GOTM_version == ver
    
    # Calculate the subset of grid points to work on.
    I = len(medsea.sea_mn)
    if c is not None:
        assert k is not None, 'Neither nor both c and k should be provided.'
        assert i0 is None and i1 is None, 'Neither i0 nor i1 should be provided when c and k are provided.'
        i0 = c*k
        i1 = c*(k+1)
        if i0 >= I:
            raise RuntimeError('The specified index range is outside of the total number of grid points: {:d}.'.format(I))
        else:
            # during run_setup(), the i1 > i0 has been checked,
            i1 = min(i1,I) # Tolerate upper index that is too big
    subset_mn = medsea.sea_mn[i0:i1]

    # Prepare the individual grid point folders.
    index_fmt = '{:0' + str(int(log10(I))+1) + 'd}'
    fn_tar = join(medsea.data_folder, 'input_c{:d}'.format(c),('input_' + index_fmt + '_' + index_fmt + '.tgz').format(i0,i1))
    if prepared_input and isfile(fn_tar):
        print('Prepared input data for {:s} found at: '.format(grid) + fn_tar)
        import tarfile
        
        gunzip_start = time()
        print('Extracting content to {:s}...'.format(wd))
        with tarfile.open(fn_tar,'r:gz') as tar:
            tar.extractall(path=wd)
        print('Elapsed {:.2f}s.'.format(time()-gunzip_start))
    else:
        print('Creating folder structure for {:d} grid points...'.format(i1-i0))
        write_dat_start = time()
        write_dat(subset_mn,medsea)
        print('Elapsed {:.2f}s.'.format(time()-write_dat_start))

    # Handle memory caching of output file
    out_fn = join(medsea.results_folder,
        ('medsea_GOTM'+'_{!s}'*3+('_'+index_fmt)*2+'.nc').format(run,grid,ver,i0,i1))
    if cached_write == True:
        fn_batch = join('/dev/shm',basename(out_fn))
    else:
        fn_batch = out_fn
        
    # Populate the batch results netCDF file header.    
    var3d = ['sst','skint','x-taus','y-taus','albedo','coszen','heat','total','swr','lwr','sens','latent']
    var4d = ['temp']
    
    print('Creating nc file header to store results containing: \n{!s}'.format(var3d))
    create_header(fn_batch,i0,i1,var3d,var4d,str(medsea.epoch))

    print('Initialization finished. Elapsed {:.2f}s.'.format(time()-init_start))
    
    # Run serially and append results.
    count = 0
    print('Writing results to {:s}...'.format(fn_batch))
    run_start = time()
    for i in range(i0,min(i1,I)):
        m, n = medsea.sea_mn[i]
        fn_pt = join(medsea.get_local_folder(m,n),
                             'results_' + run + '_' + ver + '.nc')
        # Actual run.
        try:
            medsea.local_run(m,n,start,stop,run=run,cached=False)
            # Append data immediately.
            print('Appending data point #{:d}, index ({:d},{:d})...'.format(i,m,n))
            append_start = time()
            with Dataset(fn_batch,'a',disk_less=True) as ds_out:
                with Dataset(fn_pt,'r') as ds_in:
                    ds_out['lat'][i-i0] = ds_in['lat'][0]
                    ds_out['lon'][i-i0] = ds_in['lon'][0]                
                    for varname in var3d:
                        ds_out[varname][:,i-i0] = ds_in[varname][:,0,0]
                    for varname in var4d:
                        ds_out[varname][:,i-i0,:] = ds_in[varname][:,:,0,0]
            print('Elapsed {:.2f}s.'.format(time()-append_start))
        except Exception as e:
            print(e)
                    
        count += 1
        run_elapsed = time()-run_start

        print('Average runtime per grid point: {:.1f} seconds.'.format(run_elapsed/count))
    print('\nGOTM runs finished. Elasped {:02d}h:{:02d}m:{:02d}s.'.format(int(run_elapsed/3600),
                                                                          int((run_elapsed%3600)/60),
                                                                          int(run_elapsed%60)))
    print('\nResults written to {:s}'.format(fn_batch))

    # Move the prepared_input file to disk
    if cached_write == True:
        shutil.move(fn_batch,out_fn)
    total_elapsed = time()-job_start
    print('\nTotal walltime: {:02d}h:{:02d}m:{:02d}s.'.format(int(total_elapsed/3600),
                                                              int((total_elapsed%3600)/60),
                                                              int(total_elapsed%60)))
