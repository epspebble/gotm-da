### Medsea run functions

## Helper functions
def timestr(nctime,i):
    " Return a formatted time string from a nc time variable at index i."
    from netCDF4 import datetime, num2date
    try:
        ts = datetime.strftime(num2date(nctime[i],nctime.units),'%Y-%m-%d %H:%M:%S')
    except:
        print(i,nctime[i],len(nctime))
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

def print_ctime(dt=datetime.now(),sep=' '):
    return dt.strftime('%Y%m%d' + sep + '%H%M%S')

def get_m_n(lat,lon):
    "Return the grid index (m,n) given latlong."
    return int((lat-medsea_lats[0])/0.75), int((lon-medsea_lons[0])/0.75)

def get_lat_lon(m,n):
    "Return the latlong given the grid index (m,n)."
    return medsea_lats[m],medsea_lons[n]

def get_ERA_yearly_data(folder,year,name,alias,lat_indices,lon_indices):
    from netCDF4 import Dataset
    fn = 'MEDSEA_ERA-INT_' + name + '_y' + str(year) + '.nc'
    with Dataset(os.path.join(folder,fn),'r') as nc:
        # First, confirm every time that the indices for lat and lon are correct.
        assert all(nc['lat'][ERA_lat_ind] == medsea_lats)
        assert all(nc['lon'][ERA_lon_ind] == medsea_lons)
        # Then return the data unpacked from the netCDF object.
        return nc[alias][:,lat_indices,lon_indices]

def create_dimensions(nc, epoch, lat=medsea_lats, lon=medsea_lons):
    " Declaring dimensions and creating the coordinate variables for each dimension."
    # Dimensions
    nc.createDimension('time') # unlimited
    nc.createDimension('lat', size = len(lat))
    nc.createDimension('lon', size = len(lon))
    
    # Dimension variables.
    nctime = nc.createVariable('time','i4',dimensions=('time',))
    nctime.units = 'hours since ' + str(epoch)
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
                f.write(timestr(time,i) + ' 18 2\n') # Always 18 readings and two columns.
                for j in range(18):
                    line = ('{:g} {:g}\n').format(-nc['depth'][j],nc['votemper'][i,j,m,n])
                    f.write(line)
        elif dat_fn == 'sprof':
            for i in range(len(time)):
                f.write(timestr(time,i) + ' 18 2\n') # Always 18 readings and two columns.
                for j in range(18):
                    line = ('{:g} {:g}\n').format(-nc['depth'][j],nc['vosaline'][i,j,m,n])
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

def local_dat(m,n,run='default'):
    """
    Generate *.dat files from all available data. See core_dat() for other explanations.
    """

    from netCDF4 import MFDataset
    from tempfile import mkdtemp
    from shutil import copyfile
    from glob import glob

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

    # nc_dict = dict(heat = MFDataset(os.path.join(temp_folder,'medsea_ERA_*.nc')), 
    #                met = MFDataset(os.path.join(temp_folder,'medsea_ERA_*.nc')), 
    #                tprof = MFDataset(os.path.join(temp_folder,'medsea_rea_votemper_*.nc')), 
    #                sprof = MFDataset(os.path.join(temp_folder,'medsea_rea_vosaline_*.nc')))

    print("Using nc files in " + data_folder + "...")
    nc_dict = dict(heat = MFDataset(os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_*.nc')), 
                   met = MFDataset(os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_*.nc')), 
                   tprof = MFDataset(os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_*.nc')), 
                   sprof = MFDataset(os.path.join(data_folder,'medsea_rea','medsea_rea_vosaline_*.nc')))

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
    gotm_args = prepare_run(start,stop,local_folder,lat=lat,lon=lon,out_fn='results-{:d}{:02d}'.format(year,month),**gotm_user_args)
    os.chdir(local_folder)
    stat = dict()
    try:
        tic()
        print('GOTM running at ' + local_folder + '...')
        logfn = 'GOTM_' + print_ctime(sep='_') + '.log'
        gotm(verbose=verbose, logfn=logfn, run_folder = local_folder, varsout = {})
        stat['elapsed'] = toc()
        statfn = 'stat_{:d}{:02d}.dat'.format(year,month)
        with open(statfn,'a') as f:
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
