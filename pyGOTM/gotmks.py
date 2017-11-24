# Explcitly import what we need from pyGOTM.config
from pyGOTM.config import \
    project_name, project_folder, \
    GOTM_executable, GOTM_executable_path, \
    GOTM_nml, GOTM_nml_path, \
    GOTM_version, epoch
from datetime import datetime
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
def gotm(run_folder = '.', verbose = False, logfn = 'gotm.log',
         inp_backup = False, **gotm_args):
    """ Runs GOTM in with extra functions. """
    import time
    from os.path import isfile,isdir, join
    from os import chdir, getcwd
    from shutil import copyfile

    # Remember the current working folder then walk into the run folder.
    home = getcwd() 
    if isdir(run_folder):
        chdir(run_folder)   
    else:
        raise IOError("The folder path: " + run_folder + " is invalid.")

    # Store the full filenames of all generated files.
    run_files = dict()
    
    # Transfer GOTM config namelists in the local run_folder.
    for each in GOTM_nml:
        dst = join(run_folder,each)
        if not(isfile(dst)):
            src = join(GOTM_nml_path,each)
            copyfile(src,dst)
        run_files[each] = dst
    
    # Update the config as well if user specified extra options.
    if gotm_args:
        new_cfg = updatecfg(path = run_folder, inp_backup = inp_backup, verbose = verbose, **gotm_args)
    
    # Now run the GOTM executable in this folder.
    exe = join(GOTM_executable_path,GOTM_executable)
    if verbose:       
        run_command(exe,output='PIPE')
    else:
        with open(logfn,'w') as logfile:
            run_command(exe,output=logfile)
        run_files['log'] = logfn
        
    # Return to the original working directory.
    chdir(home)
    
    # Check whether a result file is created. Should be replaced by a better exception handling structure.
    from netCDF4 import Dataset
    dr = getval(loadcfg(verbose=False),'out_dir')
    fn = getval(loadcfg(verbose=False),'out_fn') + '.nc'
    with Dataset(join(dr,fn),'r') as results:
        #print('len of time in results = ',len(results['time'])>0)
        nrec = len(results['time'])
        if nrec == 0:
            raise Exception("Invalid GOTM results! Time dimension is empty. GOTM failed?")
        else:
            #TODO: Maybe not. Sometimes GOTM terminates midway. Need to catch that exit code.
            print("GOTM run completed successfully at {!s}.".format(datetime.now()))
            run_files['results'] = join(dr,fn)        
    return run_files
      
## Treating f90 namelists used by GOTM 3.0.0
def loadcfg(path='.', verbose = True):
    from f90nml import read
    from os.path import join
    config = dict()
    for eachnml in GOTM_nml:
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
    for eachnml in GOTM_nml:
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
    for eachnml in GOTM_nml:
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

# Miscellanenous helper functions.
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

def print_ctime(dt=datetime.now(),sep=' '):
    return dt.strftime('%Y%m%d' + sep + '%H%M%S')

def print_cdate(dt=datetime.now()):
    return dt.strftime('%Y%m%d')

def tz(lon):
    if lon > 0:
        return int((lon+7.5)/15)
    else:
        return int((lon-7.5)/15)

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

def print_lat_lon(lat,lon,fmt_str='g'):
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

def subgrid_slice(src,dst):
    """
    Return a slice object 'sl' so that src[sl] == dst.

        For e.g.
        
        1. 'src' (e.g. 1.5, 2., 2.5, ... 10.) is taken from
        an arithmetic sequence (e.g. 0., 0.5, 1., ...)
        
        2. 'dst' (e.g. 2., 4., 6., 8., 10.) is taken from an 
        arithmetic subsequence of the above (e.g. 0., 1., 2., ...)
        
        3. 'src' contains all the values of 'dst'
        
        It also works if src or dst were reversed in order of elements.
    """
    
    from numpy import where, array_equal, sign, array

    def get_index_of(val):
        "Check and return the index of a unique element 'val' in 'src'."
        found = where(array(src) == val)[0]
        if len(found) == 0:
            raise ValueError("Cannot find the value {!s} in the vector: \n{!s}.".format(val,src))
        elif len(found) != 1:
            raise ValueError("Repeated elements of {!s} found in the vector: \n{!s}".format(val,src))
        else:
            return found[0]
        
    first = get_index_of(dst[0])
    last = get_index_of(dst[-1])
    step_ratio = (dst[1]-dst[0])/(src[1]-src[0])
    if step_ratio%1. != 0.:
        raise ValueError("Step size of 'dst': {:g}, is not an integer multiple of that of 'src': {:g}.".format(dst[1]-dst[0],src[1]-src[0]))
    step = int(step_ratio)
    start, stop = first, last + sign(step)
    if not array_equal(src[start:stop:step],dst):
        raise ValueError("Assumption failed. Are 'dst' and 'src' equally sized and 'src' contains all values of 'dst'?")
    else:
        return slice(start, stop, step)
