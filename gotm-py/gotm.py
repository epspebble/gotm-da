# This code is written in Python 3, but for Python 2.6 or above, we need...
from __future__ import print_function
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

### Global settings
import os
from datetime import datetime
from numpy import pi, cos, sin

## For general GOTM setup

# Set the default GOTM executable and namelist locations.
userhome = os.getenv('HOME')
project_folder = os.path.join(userhome,'medsea_GOTM')
GOTM_executable = os.path.join(project_folder,'bin','gotm')
if not(os.path.isfile(GOTM_executable)):
    raise FileNotFoundError("The GOTM executable not found at " + GOTM_executable)
GOTM_nml_path = os.path.join(project_folder,'config')
GOTM_nml_list = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp']
for nml in GOTM_nml_list:
    GOTM_nml_template = os.path.join(GOTM_nml_path,nml)
    if not(os.path.isfile(GOTM_nml_template)):
           raise FileNotFoundError("The GOTM config namelist " + GOTM_nml_template + " is invalid.")

# GOTM namelist values
timestep = 30
            
## For medsea simulations

# Top-level project folders
data_folder = os.path.join('/global/scratch',os.getenv('USER'))
base_folder = os.path.join(data_folder,'medsea_GOTM')

if not(os.path.isdir(base_folder)):
    raise IOError('The base folder: ' + base_folder + ' is either not accessible or created.')
run_folder = os.path.join(base_folder,'run')
if not(os.path.isdir(run_folder)):
    os.mkdir(run_folder)

# Ocean and Satellite products datasets source folders.
p_sossta_folder = os.path.join(data_folder,'p_sossta')
ERA_folder = os.path.join(p_sossta_folder,'medsea_ERA-INTERIM','3-hourly')
rea_folder = os.path.join(p_sossta_folder,'medsea_rea')

# GOTM dat files' netCDF reformatted dataset sources.
def data_sources(year, month, mode = 'r', dat = ['heat','met','tprof','sprof']):
    from netCDF4 import Dataset
    import os

    if isinstance(dat,str):
        dat = [dat] # So that list comprehension still works.

    fn_dict = {'heat' : os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_{:d}{:02d}.nc'.format(year,month)),
               'met'  : os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_{:d}{:02d}.nc'.format(year,month)),
               'tprof': os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_{:d}{:02d}.nc'.format(year,month)),
               'sprof': os.path.join(data_folder,'medsea_rea','medsea_rea_vosaline_{:d}{:02d}.nc'.format(year,month)),
               'sst'  : os.path.join(data_folder,'medsea_OSTIA','medsea_OSTIA_sst_{:d}{:02d}.nc'.format(year,month)),
               'chlo' : os.path.join(data_folder,'medsea_MODIS','medsea_MODIS_chlor_a_{:d}{:02d}.nc'.format(year,month))}
    assert all([each in fn_dict.keys() for each in dat])

    for each in fn_dict.keys():
        try:
            ds_dict = {each : Dataset(fn_dict[each],mode) for each in dat}
        except OSError:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            raise
        except:
            raise


    if len(ds_dict.keys()) == 1:
        # If only one dataset requested, return the netcdf dataset unwrapped from the dict. 
        return ds_dict[dat[0]] 
    else:
        # Else return a dictionary of datasets.
        return ds_dict

# Global setting for the core_dat() routines (and possibly the ERA routines as well)
overwrite=True

# Our grid points. Maybe we can reduce dependence on numpy by just using a Python array.
medsea_lats = tuple(30.75+0.75*i for i in range(21))
medsea_lons = tuple(-6.0+0.75*i for i in range(57))

# The corresponding index ranges in medsea_ERA-INTERIM datasets.
ERA_lat_ind = slice(-8,4,-1)
ERA_lon_ind = slice(12,-18,1)

# The corresponding index ranges in medsea_rea datasets.
rea_lat_ind = slice(9,250,12)
rea_lon_ind = slice(0,673,12)
rea_depth_ind = slice(0,18)

# Enumerate the grid points 
import itertools
mm, nn = zip(*itertools.product(range(21),range(57)))
# Make use of a medsea_rea dataset to infer sea, shallow and land locations.
with data_sources(2014,1,dat='tprof') as rea_ds:
    votemper = rea_ds['votemper']
    # Preallocate
    is_sea = list(None for i in range(21*57))
    is_shallow = list(None for i in range(21*57))
    is_land = list(None for i in range(21*57))
    for i in range(21*57):
        # Since fill value is 1e20, and sea water should not be boiling...
        is_sea[i] = (votemper[0,-1,mm[i],nn[i]]<100) # deepest location in our data should be about 100m.
        is_land[i] = (votemper[0,0,mm[i],nn[i]]>100) # shallowest data
        is_shallow[i] = \
            (votemper[0,0,mm[i],nn[i]]<100) and \
            (votemper[0,-1,mm[i],nn[i]]>100)

# Check that there are no logical loopholes.
assert sum(is_sea) + sum(is_land) + sum(is_shallow) == 21*57

# Return the counters i for lat/lon index arrays mm and nn.
sea_i = tuple(itertools.compress(range(21*57),is_sea))
land_i = tuple(itertools.compress(range(21*57),is_land))
shallow_i = tuple(itertools.compress(range(21*57),is_shallow))

# Return the actual lat/lon index pairs (m,n)
sea_mn = tuple((mm[i],nn[i]) for i in sea_i)
shallow_mn = tuple((mm[i],nn[i]) for i in shallow_i)
land_mn = tuple((mm[i],nn[i]) for i in land_i)

# Return the actual (lat,lon) coorindates as well
sea_locations = tuple((medsea_lats[m],medsea_lons[n]) for (m,n) in sea_mn)
shallow_locations = tuple((medsea_lats[m],medsea_lons[n]) for (m,n) in shallow_mn)
land_locations = tuple((medsea_lats[m],medsea_lons[n]) for (m,n) in land_mn)
    
### General GOTM wrappers

# Running GOTM console through a subprocess as if we were in a linux terminal.
def run_command(cmd, output='PIPE'): 
    " Execute a command as a subprocess and directing output to a pipe or a log file. "
    from subprocess import Popen, PIPE, STDOUT
    import shlex,os,_io

    # For backward-compatibility.
    if isinstance(output,bool) and output == True:
        output = 'PIPE'
    
    # New code block for Python 2
    ## Open a subprocess.
    #if output == 'PIPE':
    #    p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True)
    #    for line in p.stdout:
    #        print(line,end='')
    #elif isinstance(output,_io.TextIOWrapper) or isinstance(output,file): # Python 2/3 discrepancy here.
    #    p = Popen(shlex.split(cmd), stdout=output, stderr=STDOUT, bufsize=1, universal_newlines=True)
    #    #p.wait()
    #    #output.flush()
    #else:
    #    p = Popen(shlex.split(cmd), stdout=None, stderr=None)
    #
    #exit_code = p.poll()
    #if not(exit_code==0) and not(exit_code is None): #debug
    #    print('exit_code: ', exit_code, type(exit_code)) 
    #    raise RuntimeError("Command: " + cmd + " failed.")

    # The following does not work with Python <= 3.2
    # To a console
    if output == 'PIPE':
        with Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
            exit_code = p.poll()
    # To a file.
    elif isinstance(output,_io.TextIOWrapper) or isinstance(output,file): # Python 2/3 discrepancy here.
        with Popen(shlex.split(cmd), stdout=output, stderr=STDOUT) as p:
            exit_code = p.poll()    
    # To nothing, just run it.
    else:
        with Popen(shlex.split(cmd), stdout=None, stderr=None) as p:
            exit_code = p.poll()     

    return exit_code

# This function should be extended to output a simulation object, encapsulating the run options and the results in 
# one class. This will in turn allow us to run continuation calls easily.
def gotm(varsout = {}, run_folder = '.', verbose = False, logfn = 'gotm.log',
         GOTM_executable = GOTM_executable, GOTM_nml_templates_path = None, inp_backup = False, **gotm_args):
    """ Runs GOTM in with extra functions. """
    import os, shutil, time
    
    # Use default templates if not user-specified.
    if GOTM_nml_templates_path is None:
        GOTM_nml_templates_path = GOTM_nml_path
    
    # Remember the current working folder then walk into the run folder.
    home = os.getcwd() 
    if os.path.isdir(run_folder):
        os.chdir(run_folder)   
    else:
        raise IOError("The folder path: " + run_folder + " is invalid.")
    
    # Check for GOTM config namelists in the local run_folder.
    #GOTM_nml_list = ['gotmrun.inp','gotmmean.inp','gotmturb.inp','airsea.inp','obs.inp'] # Moved to the top, as global var.
    for each in GOTM_nml_list:
        if not(os.path.isfile(each)):
            shutil.copyfile(os.path.join(GOTM_nml_templates_path,each),os.path.join(run_folder,each))
    
    # Update the config as well if user specified extra options.
    if gotm_args:
        new_cfg = updatecfg(path = run_folder, inp_backup = inp_backup, verbose = verbose, **gotm_args)
    
    # Now run the GOTM executable in this folder.
    if verbose:       
        run_command(GOTM_executable,output='PIPE')
    else:
        with open(logfn,'w') as logfile:
            run_command(GOTM_executable,output=logfile)
    # Return to the original working directory.
    os.chdir(home)
    
    # Check whether a result file is created. Should be replaced by a better exception handling structure.
    from netCDF4 import Dataset
    dr = getval(loadcfg(verbose=False),'out_dir')
    fn = getval(loadcfg(verbose=False),'out_fn') + '.nc'
    with Dataset(os.path.join(dr,fn),'r') as results:
        #print('len of time in results = ',len(results['time'])>0)
        nrec = len(results['time'])
        if nrec == 0:
            raise Exception("Invalid GOTM results! Time dimension is empty. GOTM failed?")
        
        #dims = set()
        #varsout = set(varsout)
        #for var in varsout:
        #    dims = dims.union(results[var].dimensions)
        # also save the dimension variables to output.
        #varsout = varsout.union(dims)

        #print(dims)
        #print(varsout)

        #dim_dict = dict()
        #for dim in dims:
        # print('Retrieving {}... with shape {}'.format(dim,shape(results[dim])))
        #    dim_dict[dim] = results[dim][:]
        #    
        #var_dict = {var: {'values': results[var][:], 
        #                  'dimensions': results[var].dimensions,
        #                  'units': results[var].units} for var in list(varsout)}
        # 
        #return (dim_dict, var_dict) 
    return
      
## Treating f90 namelists used by GOTM 3.0.0
    
def loadcfg(path='.', verbose = True):
    from f90nml import read
    from os.path import join
    config = dict()
    for eachnml in GOTM_nml_list:
        fn = join(path,eachnml)
        if verbose:
            print('Reading {} ...'.format(fn))
        ### Warning! Keep failing to an infinite loop if the file is empty for some reason!
        config[eachnml] = read(fn)
    return config

def writecfg(gotm_cfg, path='.', inp_backup = False, verbose = True):
    from f90nml import write
    from os.path import exists, join
    from os import rename
    from datetime import datetime
    timestr = print_ctime(sep='_')
    for eachnml in GOTM_nml_list:
        fullfile = join(path,eachnml)
        if exists(fullfile) and inp_backup:
            # Append a suffix with the current timestamp, almost ISO8601-like, sans '-', ':' and timezone.
            rename(fullfile,fullfile[:-4] + '_' + timestr + '.inp')    
        write(gotm_cfg[eachnml],fullfile,force=True)
    if verbose:
        if inp_backup:
            print('A backup set of namelists saved at ' + timestr)
        print('GOTM config written.')
        
        
def updatecfg(path='.', inp_backup = False, verbose = True, **kwargs):
    # NOTE: Currently, this method can update multiple key/value pairs if the key names repeat among the files, i.e.
    # We assume, in spite of the hierarchy of namelists, that the name of the keys are unique across hierarchies and
    # namelist files.
    #
    # For example: 
    #
    # updatecfg(gotm_cfg, start='2014-01-01 00:00:00') will update the 'start' value in `gotmrun.inp` but will also 
    # update the 'start' value in in, say, `obs.inp` as well if it exists (luckily. this is not true).
    # Though very unlikely, still it is better to perform a test flattening the nested namelist structure to 
    # confirm the key/value pairs at leaf node level do not repeat in the name of the keys, even across several files. 
    
    #print(kwargs, ' in updatecfg()')
    import f90nml
    from os import rename, remove
    from os.path import join
    from datetime import datetime
    list_of_keys_to_update = list(kwargs.keys())
    #print(kwargs) #debug
    def recursively_update(nml, **kwargs):
        for k,v in nml.items():
            has_key = list_of_keys_to_update.count(k)
            if isinstance(v,f90nml.namelist.Namelist):
            #if isinstance(v,f90nml.NmlDict):
                new_cfg = recursively_update(v, **kwargs)
            elif has_key == 0:
                continue
            else:
                assert has_key == 1
                nml[k] = kwargs[k]
        return nml
    newcfg = recursively_update(loadcfg(path=path, verbose=False), **kwargs)
    for eachnml in GOTM_nml_list:
        inp = join(path,eachnml)
        timestr = print_ctime(sep='_')
        inpbkp = inp[:-4] + '_' + timestr + '.inp'
        rename(inp,inpbkp) 
        f90nml.patch(inpbkp,newcfg[eachnml],inp)
#        f90nml.write(newcfg[eachnml],inp)
        if not(inp_backup):
            remove(inpbkp)
        if verbose:
            if inp_backup:
                print('A backup set of namelists saved at ' + timestr)
    return 

def getval(gotm_cfg, key):
    "Return the value and walking through the hierarchy of the namelists."
    import f90nml
    # This is quite Fortran style. Maybe should use return values instead.
    result = []; # List is mutable and the nonlocal keyword allow the recursive calls to bind to THIS variable.
    def recursively_find_in(nml):
        for k,v in nml.items():
            if isinstance(v,f90nml.namelist.Namelist):
                recursively_find_in(v)
            elif k == key:
                result.append(v)
    recursively_find_in(gotm_cfg)
    if len(result) == 0:
        return None
    else: 
        assert(len(result) == 1) # Expect unique key names
        return result[0] 

### Medsea run functions

## Helper functions
def change_base(new_base_folder):
    global base_folder, run_folder
    if not(os.path.isdir(new_base_folder)):
        raise IOError('The base folder: ' + new_base_folder + 'is either not accessible or created.')
    base_folder = new_base_folder
    run_folder=os.path.join(base_folder,'run')
    if not(os.path.isdir(run_folder)):
        os.mkdir(run_folder)

def timestr(nctime,i):
    " Return a formatted time string from a nc time variable at index i."
    from netCDF4 import datetime, num2date
    try:
        ts = datetime.strftime(num2date(nctime[i],nctime.units),'%Y-%m-%d %H:%M:%S')
    except:
        print(i,nctime[i],len(nctime))
    return ts
        
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


## Medsea parallel run toolbox
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

def medsea_dat(year=2014, month=1):
    " Generate vaious GOTM dat files for all medsea grid points, using load-balanced ipyparallel engines. "
    import ipyparallel as ipp
    import itertools as itt
    
    rc = ipp.Client()
    dv = rc[:]
    lv = rc.load_balanced_view()
    with dv.sync_imports():
        import gotm
        
    mm,nn = zip(*itt.product(range(21),range(57)))
    print('Generating dat files for {0:d}{1:02d} ...'.format(year,month))
    #print('Using netCDF files from ' + base_folder)

    dv.push(dict(year=year,month=month))
    dv.execute('ds = gotm.data_sources({},{})'.format(year,month))
    run = lambda m,n: gotm.core_dat(year,month,m,n,**ds)
    results = lv.map(run,mm,nn)

    #results.wait_interactive()
    return results

def core_results(year,month,m,n,outfn='results'):
    " Convenience function for getting a netCDF dataset handle for the current core folder result."
    import os
    from netCDF4 import Dataset
    return Dataset(os.path.join(get_core_folder(year,month,m,n),outfn+'.nc'),'r')

def medsea_data(year=2014,month=1,results_folder='ASM0'):
    " Convenience function for getting a netCDF dataset handle for a monthly output."
    import os
    from netCDF4 import Dataset
    fn = os.path.join(base_folder,results_folder, 'medsea_GOTM_{:d}{:02d}.nc'.format(year,month))
    return Dataset(nc_fn,'r')

### For SWRD Calcluations

# Temporary global variable
year = 2014

## Helper functions
def tic():
    import time
    global lap_time
    lap_time = time.time()

def toc():
    import time
    try:
        elapsed = time.time() - lap_time
        print("Elapsed: {:.4g} seconds.".format(elapsed))
        return elapsed 
    except NameError:
        print("Have you run tic()?")
        
def tz(lon):
    if lon > 0:
        return int((lon+7.5)/15)
    else:
        return int((lon-7.5)/15)

def yrdays(year):
    return 366 if year % 4 == 0 else 365

def date2num(dt):
    from datetime import datetime
    ndays = (dt - datetime(dt.year,1,1)).days
    nsecs = (dt - datetime(dt.year,dt.month,dt.day,0,0,0)).seconds
    return ndays, nsecs

def UTC_to_local_nv(ndays,nsecs,lon):
    local_nsecs = nsecs + tz(lon)*3600
    local_ndays = ndays
    if local_nsecs > 86400:
        local_nsecs -= 86400
        local_ndays +=1
        
    elif local_nsecs < 0:
        local_nsecs += 86400
        local_ndays -= 1
    #if local_ndays < 0:
    #    print("WARNING: local date is in the previous year!")
    #if local_ndays > yrdays(year):
    #    print("WANRING: local date is in the next year!")
    
    return local_ndays, local_nsecs
    
def UTC_to_local(ndays,nsecs,lon):
    from numpy import vectorize
    vfunc = vectorize(UTC_to_local_nv)
    return vfunc(ndays,nsecs,lon)

def local_to_UTC_nv(ndays,nsecs,lon):
    UTC_nsecs = nsecs - tz(lon)*3600
    UTC_ndays = ndays
    if UTC_nsecs > 86400:
        UTC_nsecs -= 86400
        UTC_ndays +=1
    elif UTC_nsecs < 0:
        UTC_nsecs += 86400
        UTC_ndays -= 1
        #if UTC_ndays < 0:
        #   print("WARNING: UTC date is in the previous year!")
        #if UTC_ndays > yrdays(year):
        #   print("WARNING: UTC date is in the next year!")
    return UTC_ndays, UTC_nsecs 
        
def local_to_UTC(ndays,nsecs,lon):
    from numpy import vectorize
    vfunc = vectorize(local_to_UTC_nv)
    return vfunc(ndays,nsecs,lon)

## Low-accuracy General Solar Position Calculations by NOAA Global Monitoring Division [https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF. The value we need is coszen.

def gamma(ndays,nsecs):
    "Fractional solar year (in radians), which begins on Jan 1st noon. Assumes input in UTC."
    return 2*pi/yrdays(year)*(ndays+(nsecs/3600-12)/24)

def sundec(t):
    "Sun declination (in radians) as a function of fractional solar year (in radians)"
    return 0.006918 - 0.399912*cos(t) + 0.070257*sin(t) \
                    - 0.006758*cos(2*t) + 0.000907*sin(2*t) \
                    - 0.002697*cos(3*t) + 0.00148*sin(3*t)

def eqtime(t):
    "Equation of time (in minutes) as a function of fractional solar year (in radians)."
    # The NOAA version has 0.000075 which, according to Spencer via ... is incorrect.
    return 229.18*(0.0000075 + 0.001868*cos(t)   - 0.032077*sin(t) \
                             - 0.014615*cos(2*t) - 0.040849*sin(2*t))

def coszen_nv(ndays,nsecs,lat,lon):
    from numpy import pi
    alat = lat/180*pi
    alon = lon/180*pi
    decl = sundec(gamma(*local_to_UTC(ndays,nsecs,lon)))
    #decl = sundec(gamma(ndays,nsecs))
    #lndays,lnsecs = UTC_to_local(ndays,nsecs,lon) # Which day it is should not matter.
    time_offset = eqtime(gamma(*local_to_UTC(ndays,nsecs,lon)))+4*lon-60*tz(lon)
    #time_offset = eqtime(gamma(ndays,nsecs))+4*lon-60*tz(lon)
    #tst = lnsecs/60.0 + time_offset
    #print(tz(lon),alon,time_offset,eqtime(gamma(*local_to_UTC(ndays,nsecs,lon))),eqtime(gamma(ndays,nsecs)))
    tst = nsecs/60.0 + time_offset
    ha = tst/4-180
    thsun = ha/180*pi
    return sin(alat)*sin(decl)+cos(alat)*cos(decl)*cos(thsun)        
        
def coszen(ndays,nsecs,lat,lon):
    "Cosine of the solar zenith angle. Inputs date and time should be local."
    from numpy import vectorize
    vfunc = vectorize(coszen_nv)
    return vfunc(ndays,nsecs,lat,lon)

## Requires the solar_utils package from PyPI
#def sp_coszen(ndays,nsecs,lat,lon):
#    dt = datetime(year,1,1) + timedelta(days=ndays) + timedelta(seconds=nsecs)
#    (angles, airmass) = solpos(location=[lat, lon, tz(lon)],
#                               datetime=[dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second],
#                               weather=[1015.62055, 40.0]) # Pressure and dry-bulb temp
#    zenith, azimuth = angles # in degrees
#    return cos(zenith/180*pi)
#
#def coszen_argmax_error_seconds(ndays,lat,lon):
#    ss = linspace(0,86400)
#    sp_cc = array([sp_coszen(ndays,s,lat,lon) for s in ss])
#    cc = coszen(ndays,ss,lat,lon)
#    return ss[argmax(sp_cc)]-ss[argmax(cc)]

# Compute our swr using Rosati's formulas.
def swr_nv(ndays,nsecs,lat,lon):
    tau = 0.7  # atmospheric transmission coefficient (0.7 used in Rossatti)
    aozone = 0.09 # water vapour plus ozone absorption (0.09 used in Rossatti)
    solar = 1370  # solar constant

    # Beware if called with a nsecs > 86400 or < 0 because of adjusting the timezone.
    if nsecs > 86400:
        ndays += 1
        nsecs -= 86400
        
    cz = coszen_nv(ndays,nsecs,lat,lon)
    if cz <= 0:
        cz = 0.0
        qatten = 0.0
    else:
        qatten = tau**(1/cz)
    
    qzer  = cz*solar                       
    
    # Rosati (88) eq. (3.5-3.7)
    qdir  = qzer*qatten
    qdiff = ((1.-aozone)*qzer-qdir)*0.5
    qtot  = qdir+qdiff   
    return qtot    

def swr(ndays,nsecs,lat,lon):
    "Calculate short wave radiation for the given moment, time assumed to be local."
    from numpy import vectorize
#     assert (min(array(nsecs)) >= 0) and (max(array(nsecs))<=86400)
    vfunc = vectorize(swr_nv)
    return vfunc(ndays,nsecs,lat,lon)

def swr_3hourly_mean(ndays,nsecs,lat,lon,timestep,method='quadrature'):
    """ Calculate short wave radiation, averaged for the subsequent 3-hourly period.
    Time period begins at the given moment (ndays, nsecs), assumed to be local.
    Average taken over equally spaced samples with duration of 'timestep' in seconds. 
    """
    from scipy.integrate import quadrature
    
    # Number of samples
    ns = 3*3600/timestep
    assert ns-int(ns) == 0.0
    if ns-int(ns) == 0.0:
        ns = int(ns) 
    else: 
        raise Exception('The timestep given does not divide 3*3600 seconds.')

    ## This method is too slow in Python.
    # Find cumulative sum, then divide by # of samples to get mean.
    cumsum = 0
    if method == 'cumsum':
        for i in range(ns):
            add_secs = timestep*i
            cumsum += swr_nv(ndays+int(nsecs/86400),
                             (nsecs+add_secs)%86400,
                             lat,lon)
            swr_mean = cumsum/ns
    elif method == 'quadrature':
        swr_mean =  quadrature(lambda nsecs: swr(ndays,nsecs,lat,lon),
                               nsecs,nsecs+10800,tol=1e-1,rtol=1e-3)[0]/10800
    else:
        raise NotImplementedError('The requested method of integration: {:s}'.format(method), 'is not implemented.')

    return swr_mean
    

def swr_3hourly_mean_monthly(year,month,m,n,method='quadrature'):
    """ 3-hourly mean computed for a month (UTC day 1 of month midnight to UTC day of of next month, midnight) """
    from time import time
    from datetime import datetime
    from numpy import ones

    # The following datetimes are to be interpreted as UTC, wich is also the timezone assumed for all date values used in GOTM 
    # medsea simulations.
    start = datetime(year,month,1)
    stop = datetime(year+1,1,1) if month == 12 else datetime(year,month+1,1)
    nrec = (stop-start).days*8
    ndays_start = (start-datetime(start.year,1,1)).days
    ndays_stop = (stop-datetime(start.year,1,1)).days

    #tic = time()
    swr_mean = ones((ndays_stop-ndays_start)*8)
    
    for ndays in range(ndays_start,ndays_stop):
        for i in range(8):
            nsecs = i*3*3600
            ndays_of_month = ndays-ndays_start            
            j = ndays_of_month*8+i

            # TAKE CARE BELOW.
            # Local time not necssarily at 0, 3, 6, 9 ... etc hours. 
            # Also, swr_3hourly_mean gives the mean over the SUBSEQUENT 3 hours. 
            I_0_calc = swr_3hourly_mean(*UTC_to_local_nv(ndays,nsecs,medsea_lons[n]),
                                        medsea_lats[m],medsea_lons[n],timestep,
                                        method=method) 
            swr_mean[j] = I_0_calc
    #toc = time()
    # How long does it take for one grid point and one 3-hourly period?
    #print('time elapsed for (m,n) = ({},{}):'.format(m,n), toc-tic)
    return swr_mean

def cloud_factor_calc(year,month,m,n,swr_ERA,method='quadrature'):
    from datetime import datetime
    from time import time
    from netCDF4 import Dataset
    import os
    from gotm import medsea_lats, medsea_lons
    from numpy import ones

    # The following datetimes are to be interpreted as UTC, wich is also the timezone assumed for all date values used in GOTM 
    # medsea simulations.
    start = datetime(year,month,1)
    stop = datetime(year+1,1,1) if month == 12 else datetime(year,month+1,1)
    nrec = (stop-start).days*8
    ndays_start = (start-datetime(start.year,1,1)).days
    ndays_stop = (stop-datetime(start.year,1,1)).days
    
    #tic = time()
    cloud_factor = ones((ndays_stop-ndays_start)*8) # Defaults to clear sky value.
    swr_mean = swr_3hourly_mean_monthly(year,month,m,n,method=method)
    for ndays in range(ndays_start,ndays_stop):
        for i in range(8):
            nsecs = i*3*3600
            ndays_of_month = ndays-ndays_start
            j = ndays_of_month*8+i

            # Instead of checking I_0_ERA, which can be zero for various reasons, check the clear sky value. 
            # If it's non-zero, compute the factor, otherwise, sun below horizon, so we just keep it as the defaul
            # value of 1 as initialized. Well, it should not matter.
            I_0_calc = swr_mean[j]
            I_0_obs = swr_ERA[j,m,n]
            if I_0_calc > 1: 
                # If I_0_calc is really small but positive, we get into some trouble...
                cloud_factor[j] = I_0_obs / I_0_calc
            if I_0_obs < 0:
                # Use cloud_factor to zero out negative values of data.
                cloud_factor[j] = 0
            if I_0_obs < 1:
                # Maybe we zero out small negligible values as well?
                cloud_factor[j] = 0
    #toc = time()
    # How long does it take for one grid point and one 3-hourly period?
    #print('time elapsed for (m,n) = ({},{}):'.format(m,n), toc-tic)
    
    return cloud_factor, swr_mean

def cloud_factor_calc_monthly(year,month,use_ipp=True,append_to_ERA_dataset_now=False):
    """ Calculate cloud factor from ERA swrd and internal algorithm, then append a cloud_factor to the ERA dataset. """
    
    import itertools as itt
    ## Get the cloud_factor value for the month
    swr_ERA = data_sources(year,month)['heat']['swrd'][:]

    mm,nn = zip(*itt.product(range(21),range(57)))

    if use_ipp:
        ## Start an ipcluster to speed up.
        rc, lv = prepare_engine()
        dv = rc[:]
        run = lambda m,n: cloud_factor_calc(year,month,m,n,swr_ERA)
        
        # Push these common values to global scope of each engine.
        dv.push(dict(year=year,month=month,swr_ERA=swr_ERA))
        tic()
        results = lv.map(run,mm,nn)
        results.wait()
        toc()
    else:
        results = [cloud_factor_calc(year,month,mm[i],nn[i],swr_ERA) for i in range(21*57)]
    
    ## Append to the ERA_file.
    if not(append_to_ERA_dataset_now):
        return results

    assert append_to_ERA_dataset_now
    with data_sources(year,month,mode='a')['heat'] as ds:
        if 'cloud_factor' in ds.varibles():
            cf = ds['cloud_factor']
        else:
            cf = ds.createVariable('cloud_factor','f8',dimensions=('time','lat','lon'),zlib=True,fill_value=1e+20)
        if 'swrd_clear_sky' in ds.variables():
            swr_cs = ds['swrd_clear_sky']
        else:
            swr_cs = ds.createVariable('swrd_clear_sky','f8',dimensions=('time','lat','lon'),zlib=True,fill_value=1e+20)
        for i,(cloud_factor, swr_mean) in enumerate(results):
            m = mm[i]
            n = nn[i]
            cf[:,m,n] = cloud_factor
            swr_cs[:,m,n] = swr_mean
