from .config import *
from .gotmks import *

### Medsea run functions

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
        raise IOError('The base folder: ' + new_base_folder + 'is either not accessible or created.')
    base_folder = new_base_folder
    run_folder=os.path.join(base_folder,'run')
    if not(os.path.isdir(run_folder)):
        os.mkdir(run_folder)
       
def print_lat_lon(lat,lon,fmt_str='.2f'):
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
    "Return the grid index (m,n) given latlong."
    return int((lat-medsea_lats[0])/0.75), int((lon-medsea_lons[0])/0.75)

def get_lat_lon(m,n):
    "Return the latlong given the grid index (m,n)."
    return medsea_lats[m],medsea_lons[n]

def create_dimensions(nc, lat=medsea_lats, lon=medsea_lons):
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

def get_local_folder(lat,lon,run):
    "Return the corresponding local folder for given given grid point for the given run code."
    # Temporary hack, be forgiving if the provided lat, lon are actually indices of our medsea grid.
    if isinstance(lat,int):
        lat = medsea_lats[lat]
    if isinstance(lon,int):
        lon = medsea_lons[lon]
        
    latlong = print_lat_lon(lat,lon)
    local_folder = os.path.join(base_folder,run,latlong)
    if not(os.path.isdir(local_folder)):
        raise IOError("The local folder {:s} is not found. Have you run local_dat()?".format(local_folder))
    #     os.mkdir(local_folder)
    return local_folder 

def get_core_folder(year,month,lat,lon):
    """ Create folder structure initially or mid-way (if not yet done). 
        The innermost subfolder, which is called 'core_folder' is 
        where a GOTM run is executed for one grid point per time period.  
        It contains settings (*.inp), input data (*.dat) and output data (*.nc) """
    # Temporary hack, be forgiving if the provided lat, lon are actually indices of our medsea grid.
    if isinstance(lat,int):
        lat = medsea_lats[lat]
    if isinstance(lon,int):
        lon = medsea_lons[lon]
        
    monthly_folder = os.path.join(run_folder,'{:d}{:02d}'.format(year,month))
    if not(os.path.isdir(monthly_folder)):
        os.mkdir(monthly_folder)
    latlong = print_lat_lon(lat,lon)
    core_folder = os.path.join(monthly_folder,latlong)
    if not(os.path.isdir(core_folder)):
        os.mkdir(core_folder)
    return core_folder 

def write_dat(m,n,dat_fn,nc,outdir):
    " Write dat files for each lat/lon in medsea_lats/medsea_lons from a given netCDF Dataset or MFDataset."

    from netCDF4 import Dataset, MFDataset
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
        print(fn) # Print the filename for debug.
        if dat_fn == 'tprof':
            for i in range(len(time)):
                ndepth = nc['depth']
                f.write(timestr(time,i) + ' {0:d} 2\n'.format(len(ndepth))) # Always two columns.
                for j in range(len(ndepth)):
                    line = ('{0:g} {1:g}\n').format(-nc['depth'][j],nc['votemper'][i,j,m,n])
                    f.write(line)
        elif dat_fn == 'sprof':
            for i in range(len(time)):
                ndepth = nc['depth']
                f.write(timestr(time,i) + ' {0:d} 2\n'.format(len(ndepth))) # Always two columns.
                for j in range(len(ndepth)):
                    line = ('{0:g} {1:g}\n').format(-nc['depth'][j],nc['vosaline'][i,j,m,n])
                    f.write(line)
        elif dat_fn == 'heat':
            col = [None for i in range(4)]
            # Temporary hack #1: repeat the first record if it starts at 03:00:00 so that the 
            # simulation can start at midnight instead. Otherwise, GOTM will generate a
            # plethora of nan values.
            if timestr(time,0)[-8:] == '03:00:00':
                col[0] = timestr(time,0)[:-8] + '00:00:00'
                col[1] = nc['swrd'][0,m,n]
                col[2] = 1 if 'cloud_factor' not in nc.variables else nc['cloud_factor'][0,m,n]
                col[3] = nc['lwrd'][0,m,n]
                line = ('{:s}' + ' {:g}'*3 + '\n').format(*col)
                f.write(line)

            # The following loop is hacked temporarily to accomodate GOTM Fortran code assumption,
            # reinterpreting the timsteamp to mean the beginning of 3-hourly periods.
            for i in range(len(time)-1): #  Last record is not used.
                #col[0] = timestr(time,i)
                col[0] = timestr(time,i)
                col[1] = nc['swrd'][i+1,m,n]
                col[2] = 1 if 'cloud_factor' not in nc.variables else nc['cloud_factor'][i+1,m,n]
                col[3] = nc['lwrd'][i+1,m,n]
                line = ('{:s}' + ' {:g}'*3 + '\n').format(*col)
                f.write(line)
            
            # Temporary hack #3: include the last line being first day next month midnight, whose
            # value should not be used because of the new interpretation of timing. However, for some 
            # reason GOTM halted at the last hour of time. So let's just repeat the value 3 hours earlier
            # at 21:00:00 last day of month.
            assert timestr(time,-1)[-8:] == '00:00:00' # The last record is at midnight.
            col[0] = timestr(time,-1)
            col[1] = nc['swrd'][-1,m,n]
            col[2] = 1 if 'cloud_factor' not in nc.variables else nc['cloud_factor'][-1,m,n]
            col[3] = nc['lwrd'][-1,m,n]
            line = ('{:s}' + ' {:g}'*3 + '\n').format(*col)
            f.write(line)
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
        elif dat_fn == 'sst':
            for i in range(len(time)):
                line = '{:s} {:g}\n'.format(timestr(time,i),nc['analysed_sst'][i,m,n])
                f.write(line)
        elif dat_fn == 'chlo':
            for i in range(len(time)):
                line = '{:s} {:g}\n'.format(timestr(time,i),nc['chlor_a'][i,m,n])
                f.write(line)
        else:
            raise Exception("Requested {}.dat has no recipes defined in core_dat()".format(dat_fn))

    print('Done writing {}.\n'.format(fn))

def local_dat(m,n,run='default',dat=['heat','met','tprof','sprof']):
    """
    Generate *.dat files from all available data. See core_dat() for other explanations.
    """

    from netCDF4 import MFDataset
    from tempfile import mkdtemp
    from shutil import copyfile
    from glob import glob
    
    if isinstance(dat,str):
        dat = [dat]

    lat = medsea_lats[m]
    lon = medsea_lons[n]

    run_folder = os.path.join(base_folder,run)
    if not(os.path.isdir(run_folder)):
        os.mkdir(run_folder)

    local_folder = os.path.join(run_folder,print_lat_lon(lat,lon))
    if not(os.path.isdir(local_folder)):
        os.mkdir(local_folder)

    # temp_folder = mkdtemp(prefix=temp_base_folder)
    # ERA_files = glob(os.path.join(data_folder,'medsea_ERA-INTERIM','*.nc')) 
    # rea_files = glob(os.path.join(data_folder,'medsea_rea','*.nc'))
    # ERA_tempfiles = [copyfile(fn,os.path.join(temp_folder,os.path.basename(fn))) for fn in ERA_files]
    # rea_tempfiles = [copyfile(fn,os.path.join(temp_folder,os.path.basename(fn))) for fn in rea_files]

    print("Using nc files in " + data_folder + "...")
    nc_dict = data_sources(dat=dat)
    if not(isinstance(nc_dict, dict)):
        nc_dict = {dat[0]: nc_dict} 

    for dat_fn, nc in nc_dict.items():
        write_dat(m,n,dat_fn,nc,local_folder)
        
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
    lat = medsea_lats[m]
    lon = medsea_lons[n]
    # latlong = print_lat_lon(lat,lon)
    core_folder = get_core_folder(year,month,lat,lon) 
    
    for dat_fn, nc in nc_dict.items():
        # print(dat_fn)
        write_dat(m,n,dat_fn,nc,core_folder)
    return

def prepare_run(start,stop,out_dir,out_fn='results',m=None,n=None,lat=None,lon=None, **gotm_user_args):
    "Transfer config files and GOTM executable to the folder in which GOTM will be run."
    import shutil
    # Determine the grid point location.
    if (m is None) and (n is None):
        (m,n) = get_m_n(lat,lon)
    if (lat is None) and (lon is None):
        (lat,lon) = get_lat_lon(m,n)
    latlong = print_lat_lon(lat,lon)        
    run_name = 'medsea_GOTM, #(m,n)=({1:d},{2:d})'.format(latlong,m,n)    

    # Set up GOTM arguments.
    gotm_args = dict(name = run_name,
                     start = str(start), stop = str(stop), 
                     latitude = float(lat), longitude = float(lon), 
                     out_dir = str(out_dir), out_fn = out_fn)
    gotm_args.update(**gotm_user_args) 

    # Copy over GOTM executable and config files, the latter to be updated when gotm() is called.
    if not(os.path.exists(GOTM_executable)):
        os.symlink(GOTM_executable,os.path.join(folder,'gotm'))    
    for each in GOTM_nml_list:
        # All config files are overwritten every time GOTM is run.
        shutil.copyfile(os.path.join(GOTM_nml_path,each),os.path.join(out_dir,each))       
    
    updatecfg(path=out_dir, **gotm_args)
    return gotm_args

def local_run(year,month,m,n,run,verbose=False,**gotm_user_args):
    """ 
    
    Generate GOTM results for the (m,n)-th grid point at the specified month. Only *.dat files are expected to be 
    in the local folder. The config file, GOTM run time will be generated or copied over. 
    
    """
    from datetime import datetime
    from netCDF4 import Dataset, num2date
    lat, lon = get_lat_lon(m,n)
    local_folder = get_local_folder(lat,lon,run)
    start = datetime(year,month,1);
    stop = datetime(year,month+1,1) if month < 12 else datetime(year+1,1,1)
    gotm_args = prepare_run(start,stop,local_folder,lat=lat,lon=lon,
                            out_fn='results-{:d}{:02d}'.format(year,month),
                            assim_event_fn='assim_event-{:d}{:02d}.dat'.format(year,month),
                            sst_event_fn='sst_event-{:d}{:02d}.dat'.format(year,month),
                            **gotm_user_args)
    os.chdir(local_folder)
    stat = dict()
    try:
        tic()
        print('GOTM running at ' + local_folder + '...')
        logfn = 'GOTM_' + print_ctime(sep='_') + '.log'
        gotm(verbose=verbose, logfn=logfn, run_folder = local_folder, varsout = {})
        stat['elapsed'] = toc()
        statfn = 'stat_{:d}{:02d}.dat'.format(year,month)
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
    """ Generate GOTM results for the (m,n)-th grid point for a specified month in the 21x57 (lat,long) medsea grid. 
    
        All necessary files (*.inp, *.dat) are assumed to be present in `core_folder` (see get_core_folder(year,month,m,n)), 
        except the GOTM executable. The program changes directory into the core folder to run and the log and results are 
        both saved in core_folder. """
    from datetime import datetime
    import os, shutil

    ## Setup GOTM arguments for this run.
    start = datetime(year,month,1,0,0)
    stop = datetime(year,month+1,1,0,0) if month < 12 else datetime(year+1,1,1,0,0)
    if not m is None:
        lat = medsea_lats[m]
    if not n is None:
        lon = medsea_lons[n]
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

def combine_run(year, month, run,
                var3dnames = ['sst','skint'],
                var4dnames = ['temp'], 
                format = 'NETCDF3_CLASSIC',
                cleanup = False):
    " Combine GOTM results nc files from each grid point into a single monthly nc file. "
    from netCDF4 import Dataset
    from numpy.ma import masked_array, zeros
    
    print('Combining GOTM results for {0:d}{1:02d} ...'.format(year,month))
    outfn = os.path.join(base_folder,run,'medsea_GOTM_{0:d}{1:02d}.nc'.format(year,month))

    print('Writing dimensions and metadata...')
    elapsed = 0
    tic()

    with Dataset(outfn,'w',format=format) as nc:
        # Default dimensions for medsea.
        nctime, nclat, nclon = create_dimensions(nc)

        # Also create the depth dimension.
        nz = 150 # Could use a global setting when we refactor next time.
        nc.createDimension('depth', size = nz)
        ncdepth = nc.createVariable('depth','f4',dimensions=('depth',))
        
        # Create nc variables for each GOTM output variable specified.
        ncvar3d = {name: create_variable(nc,name,'f8', dimensions=('time','lat','lon')) for name in var3dnames}
        ncvar4d = {name: create_variable(nc,name,'f8',dimensions=('time','depth','lat','lon')) for name in var4dnames}
        
        # We actually know all our sea locations, and (30.75N,18.75E) is the first point in our medsea grid.
        fn = os.path.join(base_folder,run,print_lat_lon(*sea_locations[0]),'results-{0:d}{1:02d}.nc'.format(year,month))

        # Transfer units and dimensions.
        with Dataset(fn,'r') as first:
            nclat.units = first['lat'].units
            nclon.units = first['lon'].units
            nctime.units = first['time'].units
            ncdepth.units = first['z'].units
            for name in var3dnames:
                ncvar3d[name].units = first[name].units
            for name in var4dnames:
                ncvar4d[name].units = first[name].units
            nctime[:] = first['time'][:]
            ncdepth[:] = -first['z'][::-1] # reverse order and sign
        
        elapsed += toc()
        
        print('Begin reading data into memory...')
        tic()

        # Initialize temp arrays to store data.
        var3d_tmp = dict()
        var4d_tmp = dict()
        if month == 12:
            num_hr = (datetime(year+1,1,1)-datetime(year,month,1)).days*24;
        else:
            num_hr = (datetime(year,month+1,1)-datetime(year,month,1)).days*24;
        for name in var3dnames:
            var3d_tmp[name] = masked_array(zeros((num_hr,21,57)),mask=True)
        for name in var4dnames:
            var4d_tmp[name] = masked_array(zeros((num_hr,nz,21,57)),mask=True)
            
        # Now proceed to read from each GOTM result nc file.
        for m, n in sea_mn:
            fn = os.path.join(base_folder,run,print_lat_lon(*get_lat_lon(m,n)),'results-{0:d}{1:02d}.nc'.format(year,month))
            with Dataset(fn,'r') as each:
                for name in var3dnames:
                    var3d_tmp[name][:,m,n] = each[name][:,0,0]
                for name in var4dnames:
                    # Make sure the depth axis is reversed.
                    var4d_tmp[name][:,::-1,m,n] = each[name][:,:,0,0] 
            if cleanup:
                os.remove(fn)        
        elapsed += toc()

        print('Begin writing to {}'.format(outfn))
        tic()
        for name in var3dnames:
            ncvar3d[name][:] = var3d_tmp[name]
        for name in var4dnames:
            ncvar4d[name][:] = var4d_tmp[name]
        elapsed += toc()

        print('Finished combining GOTM results after {0:.0f}s'.format(elapsed))

### Medsea results visualization toolbox

# Import matplotlib colormap for assigning default values to functions.
from matplotlib import cm

def medsea_heatmap(data, # The only necessary argument: a netCDF Dataset to be opened by user.
                   varnames=['sst'],fun = lambda varargin: mean(varargin[0][:],axis=0),
                   ax=None, draw_colorbar=True, vlim = None,cmap=cm.coolwarm):
    
    # Expect a netCDF4 Data
    from netCDF4 import Dataset
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
    sst_max = ones((nday,21,57))*NaN
    sst_min = ones((nday,21,57))*NaN
    sst_range = ones((nday,21,57))*NaN
    for day in range(nday):
        for k in range(21*57):
            m = int(k/57)
            n = mod(k,57)
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
    sst_max = ones((nday,21,57))*NaN
    sst_min = ones((nday,21,57))*NaN
    sst_range = ones((nday,21,57))*NaN
    for day in range(nday):
        for k in range(21*57):
            m = int(k/57)
            n = mod(k,57)
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
