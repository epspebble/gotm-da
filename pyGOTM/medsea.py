#Use tabulate to beatify output if possible
try:
    from tabulate import tabulate
except ImportError:
    pass # silently ignore if it is not installed.

## TODO: Make an explicit list of functions
from pyGOTM.gotmks import *

# Explcitly import what we need from pyGOTM.config
from pyGOTM.config import \
    project_name, project_folder, \
    GOTM_executable, GOTM_executable_path, \
    GOTM_nml, GOTM_nml_path, \
    GOTM_version, epoch, verbose

# Necessary library
import os, sys

## Global settings and initializations for medsea runs.

# If a change to the following is desired, do so right after loading the module and then call set_grid() and set_folders() in this order
## TODO The region and horizontal grid level (1x) and vertical grid level (122) should be include in the run tag?
max_depth = 75 # TO BE DEPRECATED. Propose to use 'grid' to denote the vertical grid instead.
nlev = 122 # Number of levels in the truncated grid for up to 75m.
grid = '1x' # TO BE DEPRECATED. Propose to use 'res' to denote the spatial resolution instead.
region = 'medsea' # practically just a prefix of filenames for now, the code is strongly tied to this assumption

# Code names of products, which are used by the pyGOTM.reformat scripts to generate folder structure using these names.
atm_product = 'ERA-INTERIM'
ocean_product = 'MFC_midnight'
remote_sensing_product = 'MODIS'

def set_sources(new_atm_product=None,
                new_ocean_product=None,
                new_remote_sensing_product=None):
    if new_atm_product is not None:
        global atm_product
        atm_product = new_atm_product
    if new_ocean_product is not None:
        global ocean_product
        ocean_product = new_ocean_product
    if new_remote_sensing_product is not None:
        global remote_sensing_product
        remote_sensing_product = new_remote_sensing_product
    

# A switch.
overwrite = True # True means running at the same grid point will overwrite files if already present (notably the *.inp etc...)

# This replaces the ambiguous base_folder or run_folder used in the past, each run folder is a
# print_lat_lon() named subfolder of the following, which should be symlinked to fast filesystem outside of the home folder,
# e.g. /scratch/[name] on clusters, or /dev/shm/[name] on a system with sufficient RAM.
external_folder = os.path.join(project_folder,'external')
#results_folder = os.path.join(project_folder,'results','ASM{!s}-{!s}m_{:s}'.format(ASM,max_depth,grid)
results_folder = os.path.join(project_folder,'results')
cache_folder = '/dev/shm'
plots_folder = os.path.join(project_folder,'plots')

# Folders that depend on the grid choice. Should be updated when set_grid() is run.
data_folder = os.path.join(project_folder,'data',grid)
grid_folder = os.path.join(project_folder,'grid',grid) 

def set_folders(new_external_folder=None,
                new_results_folder=None,
                new_cache_folder=None,
                new_plots_folder=None,
                new_data_folder=None,
                new_grid_folder=None):
    if new_external_folder is not None:
        global external_folder
        external_folder = new_external_folder
    if new_results_folder is not None:
        global results_folder
        results_folder = new_results_folder
    if new_cache_folder is not None:
        global cache_folder
        cache_folder = new_cache_folder
    if new_plots_folder is not None:
        global plots_folder
        plots_folder = new_plots_folder
    if new_data_folder is not None:
        global data_folder
        data_folder = new_data_folder
    if new_grid_folder is not None:
        global grid_folder
        grid_folder = new_grid_folder
    return

## Global config for medsea simulations 
# The following are not meant to be changed interactively, as functions in other modules in thisos.listdir(ms.base_folder)
# package will not see the changes. 
#
# Configs that are meant to be changed in runtime should go into medsea.py.
#

# TODO: We should specify the corresponding GOTM code version in the profiles as well.

# Use a for-loop to regenerate profiles:

ASM_level = {
    # Assimilate at local midnight.
    'ASM0': dict(assimilation_type=0, extinct_method=9), # Paulson-Simpson 9-band.
    'ASM1': dict(assimilation_type=0, extinct_method=12, extinct_file='chlo.dat'), # Ohlmann-Siegel (2000) chlorophyll-a based
    # Assimlate at I_0 > 1 after local midnight.
    # 'ASM2': dict(assimilation_type=2, assim_window=1, extinct_method=9), # Paulson-Simpson 9-band.
    # 'ASM3': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat'), # Ohlmann-Siegel (2000) chlorophyll-a based
    # 'ASM4': dict(assimilation_type=2, assim_window=1, extinct_method=15), # Paulson-Simpson 9-band with Jerlov type I modification due to Verevochkin (2005)
    # 'ASM5': dict(assimilation_type=2, assim_window=1, extinct_method=13, extinct_file='iop.dat'), # Lee et al. (2003) IOP-based
    # 'ASM6': dict(assimilation_type=2, assim_window=1, extinct_method=16), # Paulson-Simpson 9-band with Jerlov type I modification due to Soloviev et al. (2005)

    # switch to assimilation at "sunrise" by coszen \approx 0.
    'ASM2': dict(assimilation_type=2, assim_window=3, extinct_method=9), # Paulson-Simpson 9-band.
    'ASM3': dict(assimilation_type=2, assim_window=3, extinct_method=12, extinct_file='chlo.dat'), # Ohlmann-Siegel (2000) chlorophyll-a based
    'ASM4': dict(assimilation_type=2, assim_window=3, extinct_method=15), # Paulson-Simpson 9-band with Jerlov type I modification due to Verevochkin (2005)
    'ASM5': dict(assimilation_type=2, assim_window=3, extinct_method=13, extinct_file='iop.dat'), # Lee et al. (2003) IOP-based
    'ASM6': dict(assimilation_type=2, assim_window=3, extinct_method=16), # Paulson-Simpson 9-band with Jerlov type I modification due to Soloviev et al. (1996)

}

ASM_level_original = ASM_level
ASM_level_expanded = {k+'S': dict(**v,t_prof_file='tprof_SST_adjusted.dat') for k,v in ASM_level.items()} # Repeat the ASM levels with an extra argument.
ASM_level = {**ASM_level_original,**ASM_level_expanded}

albedo_method = {
    # 0 : dict(albedo_method=0), # DEFAULT. Payne (1976)
    '.1': dict(albedo_method=1, albedo_file='chlo.dat'), # Ohlmann-Siegel
    '.2': dict(albedo_method=2),
}

coolskin_method = {
    # 0 : dict(coolskin_method=0), # DEFAULT. Fairall et al. (1996a)
    'a': dict(coolskin_method=1), # Artale (2002)
}

vgrid_choice = {
    '75m': dict(depth = 74.539324233308434, nlev = 122, grid_method = 2, grid_file = 'grid_75m.dat'),
    'HYB-75m': dict(depth = 74.389762997627258, nlev = 41, grid_method = 2, grid_file = 'grid_HYB_75m.dat'),
    'MFC-75m': dict(depth = 74.389762997627258, nlev = 15, grid_method = 2, grid_file = 'grid_MFC_75m.dat'),
    '100m': dict(depth = 99.282236525788903, nlev = 132, grid_method = 2, grid_file = 'grid_100m.dat'),
    '75m_4t': dict(depth = 74.539324233308434, nlev = 122, grid_method = 2, grid_file = 'grid_75m.dat', nsave=30),
    'EXT-75m': dict(depth = 74.539324233308434, nlev = 125, grid_method = 2, grid_file = 'grid_EXT_75m.dat'),    
    'EXT-75m_4t': dict(depth = 74.539324233308434, nlev = 125, grid_method = 2, grid_file = 'grid_EXT_75m.dat', nsave=30),
}

run_profiles = dict()
# 5 grids x 7 x 3  = 105 combinations.
for asm in ASM_level.keys():
    for grd in vgrid_choice.keys():
        # Base profiles using default albedo_method and default coolskin_method.
        run_key = asm + '-' + grd
        run_profiles[run_key] = dict(**ASM_level[asm],**vgrid_choice[grd])

        # Add albedo sub-profiles:
        for alb in albedo_method.keys():
            run_key = asm + alb + '-' + grd
            run_profiles[run_key] = dict(**ASM_level[asm],**vgrid_choice[grd],**albedo_method[alb])

            # Add coolskin sub-sub-profiles.
            # NOTE: the combinations like ASM3a-75m is not included in such a loop.
            for csk in coolskin_method.keys():
                run_key = asm + alb + csk + '-' + grd
                run_profiles[run_key] = dict(**ASM_level[asm],**vgrid_choice[grd],**albedo_method[alb], **coolskin_method[csk])
            #    run_key = asm + alb + csk

### FIX ME: Temporarily put some module-wide default global settings here.
asm = 'ASM3' # This selects the GOTM extra parameters profile
alb = '.1'
csk = ''
grd = '75m'
run_key = asm + alb + csk + '-' + grd

def set_run(new_run_key=run_key,asm=None,alb=None,csk=None,grd=None):
    global run_keu
    if asm is None:
        assert alb is None
        assert csk is None
        assert grd is None
        assert new_run_key is not None
        run_key = new_run_key
    else:
        assert alb is not None
        assert csk is not None
        assert grd is not None
        assert new_run_key is None
        new_run_key = asm + alb + csk + '-' + grd
        if new_run_key in run_profiles.keys():
            run_config = run_profiles[new_run_key]
            if verbose:
                print('Preset run profile: ' + new_run_key + ' selected:')
                # Print run profile details.
                if 'tabulate' not in sys.modules.keys():
                    for key in run_config:
                        print('\t',key,':\t', run_config[key])
                else:
                    print(tabulate([[name,val] for name, val in run_config.items()]))
            else:
                print('WARNING: unknown run profile: ' + new_run_key)

        run_key = new_run_key
    return run_key

# List of run_profiles done in the past.

# run_profiles = {
#     # The V2 runs. They were run at 150m deep with 150 levels.
# #    'ASM0': dict(assimilation_type=0, extinct_method=9),
# #    'ASM1': dict(assimilation_type=0, extinct_method=12),
# #    'ASM2': dict(assimilation_type=2, assim_window=1, extinct_method=9),
# #    'ASM3': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat',),

#     # Alternative vertical grid runs.
#     # Maybe we should dynamically calculate 'nlev' and 'depth' from number of lines in the 'grid.dat' files.
#     'ASM3-100m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat',
#                       depth = 99.282236525788903, nlev = 132,
#                       grid_method = 2, grid_file = 'grid_100m.dat'), 
#     'ASM3-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
#                      depth = 74.539324233308434, nlev = 122,
#                      grid_method = 2, grid_file = 'grid_75m.dat'),
#     'ASM3-MFC-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
#                          depth = 74.389762997627258, nlev = 15,
#                          grid_method = 2, grid_file = 'grid_MFC_75m.dat'),
    
#     # Added for completion and comparison to ASM3-75m
#     'ASM0-75m': dict(assimilation_type=0, extinct_method=9, 
#                      depth = 74.539324233308434, nlev = 122,
#                      grid_method = 2, grid_file = 'grid_75m.dat'),
#     'ASM1-75m': dict(assimilation_type=0, extinct_method=12,
#                      depth = 74.539324233308434, nlev = 122,
#                      grid_method = 2, grid_file = 'grid_75m.dat'),
#     'ASM2-75m': dict(assimilation_type=2, assim_window=1, extinct_method=9, 
#                      depth = 74.539324233308434, nlev = 122,
#                      grid_method = 2, grid_file = 'grid_75m.dat'),
#     # 20171031 Previously ASM4-75m uses extinct_method=13, but the GOTM code
#     # has not incorporated this case and it fell through to the default case
#     # which is extinct_method=1, Jerlev (1976) water type 1.
#     'ASM4-75m': dict(assimilation_type=2, assim_window=1,
#                      extinct_method=15,  
#                      depth = 74.539324233308434, nlev = 122,
#                      grid_method = 2, grid_file = 'grid_75m.dat'),
#     # 20170912 Hybrid grid
#     'ASM3-HYB-75m': dict(assimilation_type=2, assim_window=1, extinct_method=12, extinct_file='chlo.dat', 
#                          depth = 74.389762997627258, nlev = 41,
#                          grid_method = 2, grid_file = 'grid_HYB_75m.dat'),
#     # 20171031 Two levels as per email with Sam Pimentel, 2017-10-27.
#     # This is identical to the previous ASM4-75m, now should properly use iop.dat
#     'ASM5-75m': dict(assimilation_type=2, assim_window=1,
#                      extinct_method=13, extinct_file='iop.dat', 
#                      depth = 74.539324233308434, nlev = 122,
#                      grid_method = 2, grid_file = 'grid_75m.dat'),

#     # This is identical to ASM3-75m, the difference lies in the GOTM code. The fixed GOTM
#     # code has the albedo issue fixed, so that Payne's albedo is not applied to Ohlmann-Siegel (2000)
#     # formulas which implicitly include albedo.
#     'ASM3.1-75m': dict(albedo_method=1, albedo_file='chlo.dat', # Now albedo method defaults to Payne's, so we need to change it explicitly.
#                        assimilation_type=2, assim_window=1, 
#                        extinct_method=12, extinct_file='chlo.dat', 
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),

#     # 20171106, same as the previous, except with a higher resolution output.
#     'ASM5-75m_4t': dict(assimilation_type=2, assim_window=1,
#                         extinct_method=13, extinct_file='iop.dat', 
#                         depth = 74.539324233308434, nlev = 122,
#                         grid_method = 2, grid_file = 'grid_75m.dat',
#                         nsave=30),
#     'ASM3.1-75m_4t': dict(assimilation_type=2, assim_window=1,
#                           extinct_method=12, extinct_file='chlo.dat', 
#                           depth = 74.539324233308434, nlev = 122,
#                           grid_method = 2, grid_file = 'grid_75m.dat',
#                           nsave=30),
    
#     # # assim_window=2 means use the tprof profile time as assimilation time.
#     # 'ASM5-75m_no_salt': dict(assimilation_type=2, assim_window=2,
#     #                            s_prof_method=0, # skip s_prof 
#     #                            extinct_method=13, extinct_file='iop.dat', 
#     #                            depth = 74.539324233308434, nlev = 122,
#     #                            grid_method = 2, grid_file = 'grid_75m.dat'),
#     # 'ASM3.1-75m_no_salt': dict(assimilation_type=2, assim_window=2,
#     #                              s_prof_method=0, # skip s_prof 
#     #                              extinct_method=12, extinct_file='chlo.dat', 
#     #                              depth = 74.539324233308434, nlev = 122,
#     #                              grid_method = 2, grid_file = 'grid_75m.dat'),
#     # 'ASM2-75m_no_salt': dict(assimilation_type=2, assim_window=2, extinct_method=9, 
#     #                          s_prof_method=0, # skip s_prof
#     #                          depth = 74.539324233308434, nlev = 122,
#     #                          grid_method = 2, grid_file = 'grid_75m.dat'),

#     # As per emails conversations with Sam on 2017-12-20,
#     # New GOTM switches:
#     # * 'albedo_method', default is the built-in Payne (1972)
#     # * 'coolskin_method', default is Fairall (1996a)
#     # Recall also that 'extinct_method' default is Jerlov water type I (1976)
#     'ASM3.1a-75m': dict(albedo_method = 1, albedo_file = 'chlo.dat', # Ohlmann-Siegel (2000)'s approximation
#                         extinct_method = 12, extinct_file = 'chlo.dat',
#                         coolskin_method = 1, # Artale (2002)
#                         assimilation_type=2, assim_window=1,
#                         depth = 74.539324233308434, nlev = 122,
#                         grid_method = 2, grid_file = 'grid_75m.dat'),
#     'ASM4.1-75m': dict(albedo_method = 1, albedo_file = 'chlo.dat', # Ohlmann-Siegel (2000)'s approximation
#                        extinct_method = 15, # Paulson-Simption 9-band with Jerlov type I modification due to Verevochkin & Startsev, J. Fluid Mech. (2005)
#                        coolskin_method = 0, # Fairall (1996a)
#                        # No extinct_file needed
#                        assimilation_type=2, assim_window=1,
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),
#     'ASM5.1-75m': dict(albedo_method = 1, albedo_file = 'chlo.dat', # Ohlmann-Siegel (2000)'s approximatio
#                        extinct_method = 13, extinct_file='iop.dat',
#                        coolskin_method = 0, # Fairall (1996a)
#                        assimilation_type=2, assim_window=1,
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),
#     'ASM5.2-75m': dict(albedo_method = 2, # Jin et al. (2011)
#                        extinct_method = 13, extinct_file='iop.dat',
#                        coolskin_method = 0, # Fairall (1996a)
#                        assimilation_type=2, assim_window=1,
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),
    

#     # 20180227 Accidentally ran "ASM2.1", might as well include it.
#     'ASM2.1-75m': dict(albedo_method = 1, albedo_file = 'chlo.dat', # Ohlmann-Siegel (2000)'s approximatio
#                        assimilation_type=2, assim_window=1, extinct_method=9, 
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),

#     # Completing the run matrix for Jin albedo.
#     'ASM2.2-75m': dict(albedo_method = 2, # Jin et al. (2011)
#                        assimilation_type=2, assim_window=1, extinct_method=9, 
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),
    
#     'ASM3.1-75m': dict(albedo_method=1, albedo_file='chlo.dat', # Now albedo method defaults to Payne's, so we need to change it explicitly.
#                        assimilation_type=2, assim_window=1, 
#                        extinct_method=12, extinct_file='chlo.dat', 
#                        depth = 74.539324233308434, nlev = 122,
#                        grid_method = 2, grid_file = 'grid_75m.dat'),
    
# }

# Alias
def setup(*args,**kwargs):
    return set_grid(*args,**kwargs)


# Routines to set global values in this module. Can be used in interactive session to change config.
def set_grid(new_grid=grid, #new_ASM=ASM,
             new_max_depth=max_depth, # These names just need to be different... Because we cannot declare an input name global below...
             subindices=None,
             return_grid_data = False,
             verbose = verbose, stat = True, plot = False
             ):
    """ Obtain the 1/16 degree grid used in medsea_rea, and classify each grid points according to 'max_depth'.
        A sub-grid is set to the global variables in the module, and also returned by specifying 'subindices':
            ** 9x test grid (the only one with both lat and lon being multiples of 0.25):
                    subindices=(slice(1,None,4), slice(0,None,4))
            ** 9x grids (16 of them, including the 9x test grid above):
                    subindices=(slice(i,None,4), slice(j,None,4))
               for (i,j) in {0,1,2,3} x {0,1,2,3}
            ** 1x grid that is co-locational with ERA data grid:
                    subindices=(slice(9,None,12), slice(None,None,12))
            ** a mini grid with only 23 sea locations with depth >= 75 that can be used for testing:
                    subindices=(slice(1,None,48), slice(None,None,48))
               NOTE: In medsea_ERA dataset, the latitudes are arranged, exceptionally, in descending order.

       Returns three tuples:
            subgrid = (grid_lats, grid_lons, medsea_flags, medsea_bdy, max_depth)
            rea_indices = (medsea_rea_lat_ind, medsea_rea_lon_ind, medsea_rea_ndepth)
            grid_indices = (M, N, sea_mn, sea_m, sea_n)
    """

    # Global folder names to be used in this function.
    global grid_folder, data_folder, results_folder

    # Declaring global is necessary for modifying them interactively after importing this module.
    global run, grid, max_depth#, ASM
    global grid_lats, grid_lons, medsea_flags, max_depth
    global medsea_rea_lat_ind, medsea_rea_lon_ind, medsea_rea_ndepth
    global M, N, sea_mn, sea_m, sea_n

    # For plotting heatmaps more conveniently
    global extent
    
    # Skip the whole function if a cache file is found.
    import numpy as np
    grid = new_grid
    new_grid_folder = os.path.join(project_folder,'grid',new_grid)
    new_data_folder = os.path.join(project_folder,'data',new_grid)
    set_folders(new_grid_folder = new_grid_folder, new_data_folder = new_data_folder)
    grid_cache_fn = os.path.join(new_grid_folder,'grid_data.npy')
    if os.path.isfile(grid_cache_fn):
        subgrid, rea_indices, grid_indices = np.load(grid_cache_fn)
        (grid_lats, grid_lons, medsea_flags, max_depth) = subgrid
        (medsea_rea_lat_ind, medsea_rea_lon_ind, medsea_rea_ndepth, loc_type) = rea_indices
        (M, N, sea_mn, sea_m, sea_n) = grid_indices
        if return_grid_data:        
            return subgrid, rea_indices, grid_indices
        else:
            #print('Using grid cache at {:s}...'.format(grid_cache_fn))
            return

    print('Grid cache not found, generating grid data...')
    from netCDF4 import Dataset
    from matplotlib import cm

    # Override the global variables using values passed from function call.
    grid = new_grid
    #ASM = new_ASM
    grid_folder = os.path.join(project_folder,'grid',grid)
    if not os.path.isdir(grid_folder):
        print('Creating new grid folder: {!s}'.format(grid_folder))
        try:
            os.mkdir(grid_folder)
        except Exception as e:
            print(e)
    data_folder = os.path.join(project_folder,'data',grid)
    if not os.path.isdir(data_folder):
        print('Creating new data folder: {!s}'.format(data_folder))
        try:
            os.mkdir(data_folder)
        except Exception as e:
            print(e)

    max_depth = new_max_depth

    # # Update the run
    # run = 'ASM{!s}-{!s}m'.format(ASM,max_depth)
    # if run in run_profiles.keys():
    #     run_config = run_profiles[run]
    #     if verbose:
    #         print('Preset run profile: ' + run + ' selected:')
    #         # Print run profile details.
    #         if 'tabulate' not in sys.modules.keys():
    #             for key in run_config:
    #                 print('\t',key,':\t', run_config[key])
    #             else:
    #                 print(tabulate([[name,val] for name, val in run_config.items()]))
    # else:
    #     print('WARNING: unknown run profile: ' + run)

    # Update the results folder
    #results_folder = os.path.join(project_folder,'results','ASM{!s}-{!s}m_{:s}'.format(ASM,max_depth,grid))
    
    # Set up the slices for the subgrid using the name.
    if new_grid in ['9x_{:d}'.format(i) for i in range(16)]:
        k = int(new_grid[3:])
        #print('k=',k)
        if subindices is not None:
            raise Exception("Mistakes? 'subindices' need not be given if a known 'grid' is provided.")
        else:
            i = int(k/4)
            j = k%4
            subindices = (slice(i,None,4), slice(j,None,4))

    if new_grid == 'aeg_144x':
        subindices = (slice(77,186,None), slice(456,553,None))

    if new_grid == 'aeg_36x':
        subindices = (slice(77,186,2), slice(456,553,2))
        
    if new_grid == 'aeg_9x':
        subindices = (slice(77,186,4), slice(456,553,4))
    
    if new_grid == '144x':
        subindices = (slice(None,None,None),slice(None,None,None))

    # 2017-10-24
    if new_grid == '36x':
        subindices = (slice(1,None,2), slice(None,None,2))

    if new_grid == '9x':
        subindices = (slice(1,None,4), slice(None,None,4))
        
    if new_grid == '1x':
        subindices = (slice(9,None,12), slice(None,None,12))

    if new_grid == 'mini':
        subindices=(slice(1,None,48), slice(None,None,48))

    medsea_rea_lat_ind = subindices[0]
    medsea_rea_lon_ind = subindices[1]
    #print(medsea_rea_lat_ind)
    #print(medsea_rea_lon_ind)

    if verbose:
        print('Initializing grid...')
    # Load the rea grid, and a sample set of data for its masks.
    with Dataset(os.path.join(external_folder, 'medsea_rea/2013/20130101_TEMP_re-fv6.nc'),'r') as ds:
        lat_rea = ds['lat'][:]
        lon_rea = ds['lon'][:]
        temp_rea = ds['votemper'][:]
        depth_rea = ds['depth'][:]

    # Without the 'safety' of adding one more level, we can get more grid points...
    #medsea_rea_ndepth = sum(depth_rea<max_depth)

    medsea_rea_ndepth = sum(depth_rea<max_depth)+1

    # 2 means deeper than max_depth, 1 means less than max_depth, 0 means land
    loc_type = 0 + ~temp_rea[0,0,:].mask + ~temp_rea[0,medsea_rea_ndepth,:].mask

    # Setting values to global names.
    grid_lats = lat_rea[subindices[0]]
    grid_lons = lon_rea[subindices[1]]
    medsea_flags = loc_type[subindices[0],subindices[1]]
    #print(grid_lats,grid_lons)
    assert grid_lats.shape, grid_lons.shape == medsea_flags.shape

    M, N = medsea_flags.shape
    assert M == grid_lats.size and N == grid_lons.size

    sea_m, sea_n = np.where(medsea_flags==2)
    sea_mn = [(sea_m[i],sea_n[i]) for i in range(sea_m.size)]
    assert len(sea_mn) == sea_m.size

    # Compute water body boundary used in this subgrid level.
    from numpy import zeros
    from itertools import product
    medsea_bdy = zeros(medsea_flags.shape)

    for m,n in product(range(M), range(N)):
        if medsea_flags[m,n] > 0:
            continue
        m_nbr = list()
        n_nbr = list()
        if m < M-1:
            m_nbr.append(m+1)
        if m > 0:
            m_nbr.append(m-1)
        if n < N-1:
            n_nbr.append(n+1)
        if n > 0:
            n_nbr.append(n-1)
        for mb, nb in product(m_nbr,n_nbr):
            if medsea_flags[mb,nb] > 0:
                medsea_bdy[m,n] = 1
    
    #print(grid_lats.size,grid_lons.size,grid_lats.min(),grid_lats.max(),grid_lons.min(),grid_lons.max())
    extent = (grid_lats.min(),grid_lats.max(),grid_lons.min(),grid_lons.max())
    if verbose:
        print('Finished setting up a subgrid of shape {!s} x {!s} '.format(M,N) +
              'with {!s} <= latitude <= {!s}, {!s} <= longitude <= {!s}.'.format(*extent))
    if verbose and stat:
        # The following are for the current subgrid.
        def print_stat(bl_array):
            drei = bl_array.size
            zwei = (bl_array == 2).sum()
            ein = (bl_array == 1).sum()
            null = (bl_array == 0).sum()
            assert null+ein+zwei == drei
            print('\t' + 'sea (>{!s}m) locations: {!s}'.format(max_depth,zwei) + ', {:4.2f}% out of {!s}'.format(zwei*100/drei,drei))
            print('\t' + 'sea (<{!s}m) locations: {!s}'.format(max_depth,ein) + ', {:4.2f}% out of {!s}'.format(ein*100/drei,drei))
            print('\t' + 'land locations: {!s}'.format(null) + ', {:4.2f}% out of {!s}'.format(null*100/drei,drei))

        print('In the current grid:')
        print_stat(medsea_flags)
        print('In the full medsea_rea grid:')
        print_stat(loc_type)

    if plot:
        import matplotlib.pyplot as plt
        # fig, axes = plt.subplots(1,2,figsize=(14,14))
        # ax1, ax2 = axes
        # ax1.imshow(loc_type,origin='lower',extent=[lon_rea.min(),lon_rea.max(),lat_rea.min(),lat_rea.max()],cmap=cm.Blues)
        # ax1.set_title('full medsea_rea grid')
        # ax2.imshow(medsea_flags,origin='lower',extent=[grid_lons.min(),grid_lons.max(),grid_lats.min(),grid_lats.max()],cmap=cm.Blues)
        # ax2.set_title('current subgrid')
        fig, ax1 = plt.subplots()
        ax1.imshow(loc_type,origin='lower',extent=[lon_rea.min(),lon_rea.max(),lat_rea.min(),lat_rea.max()],cmap=cm.Blues)
        ax1.set_title('full medsea_rea grid')
        fig, ax2 = plt.subplots()
        ax2.imshow(medsea_flags,origin='lower',extent=[grid_lons.min(),grid_lons.max(),grid_lats.min(),grid_lats.max()],cmap=cm.Blues)
        ax2.set_title('current subgrid')
        for ax in (ax1,ax2):
            ax.set_xlabel('longitude')
            ax.set_ylabel('latitude')

    # These values have been written directly to global variables as well:
    subgrid = (grid_lats, grid_lons, medsea_flags, max_depth)
    rea_indices = (medsea_rea_lat_ind, medsea_rea_lon_ind, medsea_rea_ndepth, loc_type)
    grid_indices = (M, N, sea_mn, sea_m, sea_n)

    # If the function has not returned in the begining... Pack up and generate a cache file.
    assert not os.path.isfile(grid_cache_fn)
    np.save(grid_cache_fn, (subgrid, rea_indices, grid_indices))    
    
    if return_grid_data:
        return subgrid, rea_indices, grid_indices

# Set default grid
set_grid('1x')
    
# Print some settings upon loading of this 
data_sources = dict(
    heat = 'medsea_ERA_INTERIM',
    met = 'medsea_ERA_INTERIM',
    tprof = 'medsea_MFC_midnight',
    sprof = 'medsea_MFC_sunrise',
    chlo = 'medsea_MODIS',
    iop = 'medsea_MODIS')

# Return a plot of the water mass:
def plot_water(ax=None):
    from netCDF4 import Dataset
    from matplotlib import cm
    with Dataset(os.path.join(external_folder, 'medsea_rea/2013/20130101_TEMP_re-fv6.nc'),'r') as ds:
        temp_rea = ds['votemper'][:]
    is_land  = temp_rea[0,0,:].mask
    
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots()
    ax.imshow(is_land,origin='lower',cmap=cm.binary)
    fig = ax.get_figure()
    return fig, ax

# GOTM dat files' netCDF reformatted dataset sources.
def get_data_sources(dat, year=None, month=None, mode='r'):
#                    dat=['heat','met','tprof','sprof','chlo','iop']):
    """
    Return the netCDF4 Dataset (or MFDataset) handles for the data source.
    Omitting month returns yearly data, omitting both year and month returns
    all available data.

    Requires the module global dictionary 'data_source' to be hard coded, and 'dat' 
    to be in the keys.

    The dataset returns in read-only mode.

    """
    from netCDF4 import Dataset, MFDataset
    import os

    if (year is None) or (month is None):
        MF = True
        # Need to use MFDataset
        NCDataset = MFDataset
    else:
        # When both year and month is given, we can specifically return handle to our reformatted datasets
        # organized by months.
        NCDataset = Dataset

    if isinstance(dat,str):
        # if only one type is requested, still make it into a list so that list comprehension still works.
        dat = [dat]

    # nc_dict = dict(heat = MFDataset(os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_*.nc')),
    #                met = MFDataset(os.path.join(data_folder,'medsea_ERA-INTERIM','medsea_ERA_*.nc')),
    #                tprof = MFDataset(os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_*.nc')),
    #                sprof = MFDataset(os.path.join(data_folder,'medsea_rea','medsea_rea_vosaline_*.nc')))

    def src(name):
        "Data source name for each 'name'.dat file."
        assert isinstance(name,str)

        if name == 'heat' or name == 'met':
            return atm_product
        elif name == 'tprof' or name == 'sprof':
            return ocean_product
        elif name == 'chlo' or name == 'iop':
            return remote_sensing_product
        else:
            raise NotImplementedError('Data source unknown for ' + name)
        return None

    def fullfn(name, subfolder=None):
        '''
        Return full absolute path to the nc file. Specify a subfolder for
        flexility of referring to different versions of the datasets.
        '''
        folder = os.path.join(data_folder, region + '_' + src(name))

        # Building the suffix that decreases the time range of coverage, use
        # wildcards if necessary.
        if year is not None:
            if month is None:
                # Both *_2013.nc and *_201301.nc through *201312.nc will be
                # caught be this pattern. Be careful in the folder content.
                suffix = '_{:d}*.nc'.format(year)
            else:
                suffix = '_{:d}{:02d}.nc'.format(year,month)
        else: # So, 'year' is None
            assert month is None
            # Hopefully, the time units of each nc file in the folder has
            # the same unit and epoch, so that MFDataset opens them properly.
            suffix = '_*.nc'

        fn = region + '_' + src(name) + '_' + name + suffix

        if subfolder:
            fn = os.path.join(subfolder,fn)

        return os.path.join(folder,fn)

    # Building up the dictionary to map the appropriate filenames.
    fn_dict = {name: fullfn(name) for name in dat}

    # Now we atta
    for each in fn_dict.keys():
        try:
            ds_dict = {each : NCDataset(fn_dict[each],mode) for each in dat}
        except OSError:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            raise
        except:
            print('Error accessing {:s}.'.format(fn_dict[each]))
            print('Requested list of dat files: {!s}.'.format(dat))
            print('Available built-in data sources: {!s}.'.format(fn_dict))
            raise

    if len(ds_dict.keys()) == 1:
        # If only one dataset requested, return the netcdf dataset unwrapped.
        return ds_dict[dat[0]]
    else:
        # Else return a dictionary of datasets.
        return ds_dict

## Helper functions

## TODO.
# def chunk(c,k,fun,*args,**kwargs):
#     """
#     calls fun() for the c-th chunk with k grid points (c-1)*10, (c-1)*10+1, ... (c-1)*10+9.
#     This partitions the 390 grid points on sea into 30 chunks.
#     """

def get_results(region=None,run_key=None,grid=None,ver=None,year=2014):
    "If any of region, run_key, grid, ver are unspecified, they are taken from global module state." 
    import sys
    current_state = sys.modules[__name__]
    if region is None:
        region = current_state.region
    if run_key is None:
        run_key = current_state.run_key
    if grid is None:
        grid = current_state.grid       
    if ver is None:
        ver = current_state.GOTM_version    
        
    from netCDF4 import Dataset
    fn = region + '_GOTM_' + run_key + '_' + grid + '_' + ver + '_' + str(year) + '.nc'
    
    return Dataset(os.path.join(project_folder,'results',run_key + '_' + grid,fn),'r')
    
def get_daily_stats(region=None,run_key=None,grid=None,ver=None,year=2014):
    "If any of region, run_key, grid, ver are unspecified, they are taken from global module state." 
    import sys
    current_state = sys.modules[__name__]
    if region is None:
        region = current_state.region
    if run_key is None:
        run_key = current_state.run_key
    if grid is None:
        grid = current_state.grid       
    if ver is None:
        ver = current_state.GOTM_version    
        
    from netCDF4 import Dataset
    fn = 'daily_stat_' + run_key + '_' + grid + '_' + ver + '_' + str(year) + '.nc'
    
    return Dataset(os.path.join(project_folder,'results',run_key + '_' + grid,fn),'r')

def get_m_n(*args,silent=False):
    """
    get_m_n(lat,lon) returns the closest grid point index (m,n).

    get_m,n(i) returns the i-th precomputed sea location in the grid.

    get_m,n(m,n) returns (m,n) with whether (m,n) correspond to a grid point on the sea checked.
    
    We assume a uniform and equal spacing in both directions of the grid.
    """
    if len(args) == 1:
        i = args[0]
        return sea_m[i],sea_n[i]
    elif len(args) == 2:
        # Check if it is a valid pair of (m,n)
        m, n = args
        #print(sea_m.size,m in sea_m, isinstance(m,int),sea_n.size,n in sea_n,isinstance(n,int))
        if (m%1.==0. and int(m) in sea_m) and (n%1.==0. and n in sea_n):
            return int(m), int(n)
        else:
            #print('Maybe {!s} are (lat, lon)?'.format(args))
            lat, lon = args
            spacing = grid_lats[1]-grid_lats[0]
            # Calculate m, n
            m = (lat-grid_lats[0])/spacing
            n = (lon-grid_lons[0])/spacing
            if m%1. != 0. or n%1. != 0.:
                if not silent:
                    print('Warning: given coordinate: ({!s},{!s}) does not correspond to a point in the {:s} grid.'.format(lat,lon,grid))
                    print('Attempting to return the closest grid point if it is on the sea.')
                    m = round(m)
                    n = round(n)
                
            # Convert to integer
            m = int(m)
            n = int(n)
            if m in sea_m and n in sea_n:
                return m, n
            else:
                raise RuntimeError('(lat,lon) = ({!s},{!s}) is neither on the {:s} grid nor is its closest grid point a sea location.'.format(lat,lon,grid))
    else:
        print(args)
        raise TypeError('Wrong number of positional arguments given.')
            
def get_lat_lon(*args,silent=False):
    """
    get_lat_lon(lat,lon) returns the (lat,lon) of the closest grid point.

    get_lat_lon(m,n) returns the (lat,lon) of the grid at index (m,n)
  
    get_lat_lon(i) returns the (lat,lon) of the i-th precomputed sea location in the grid.

    We assume a uniform and equal spacing in both directions of the grid.
    """
    try:
        m, n = get_m_n(*args,silent=silent)
    except TypeError:
        raise
    
    return grid_lats[m],grid_lons[n]

## Medsea serial / parallel run toolbox
def get_local_folder(*args, create=False):
    """
    Returning the path to a run folder.

    get_local_folder(m,n) return the run folder given the grid point indices (m,n)

    get_local_folder(lat,lon) return the specified or the closest grid point's run folder.

    get_local_folder(i) return the run folder cooresponding (i+1)-th precomputed sea location (i.e. m, n = sea_m[i], sea_n[i]).

    get_local_folder() return a random precomputed sea location.

    Specify 'create=True' to create the folder if it is not yet existent.
    """
    ## Setting default argument values.
    if len(args) == 0:
        # print('No grid indices given, defaulting to the favourite grid point (lat, lon) = (36.00, 25.50) near buoy 61277...')
        # m, n = get_m_n(36.,25.5)
        from random import randint
        i = randint(0,sea_m.size-1)
        m, n = sea_m[i], sea_n[i]
        print('Randomly picking the grid point with indices: m =',m,'n =',n)
    elif len(args) == 1 or len(args) == 2:
        m, n = get_m_n(*args)
    else:
        print('Positional arguments given: ', args)
        raise Exception('Wrong number of positional arguments given.')

    ## Set up local variables
    lat  = grid_lats[m]
    lon = grid_lons[n]
    latlong = print_lat_lon(lat,lon)
    local_folder = os.path.join(grid_folder,latlong)

    if not(os.path.isdir(local_folder)):
        if create:
            try:
              #print('The folder {:s} is not found, creating it.'.format(local_folder))
              os.mkdir(local_folder)
            except:
              raise
        else:
            raise IOError("The local folder {:s} is not found. Have you run local_dat()?".format(local_folder))
    return local_folder

def prepare_run(*args, # Necessary GOTM run arguments: location i, or (m,n), or (lat,lon); and start/stop,
                run = set_run(),
                year = None, month = None, # Provided to compute start, stop.
                out_dir = '.', # Specifiy where to write GOTM ouputs if not the current folder... for caching purpose.
                **gotm_user_args): # Extra GOTM arguments, e.g. from run_profiles):
            
    """
    Prepare the run folder (which itself should not be generated if not existent, because the .dat files should be generated before calling this): 
    1. transfer config files, and update them
    2. transfer GOTM executable with the correct version
    3. generate GOTM input data (heat.dat, met.dat, tprof.dat etc.) if not found. 
    4. to the folder in which GOTM will be run.

    Positional arguments:
     - prepare_run(m,n,start,stop) 
     - prepare_run(i,start,stop) 
     - prepare_run(lat,lon,start,stop)"
    
    chooses the grid point by (m,n) as specified or as returned by get_m_n(i) / get_m_n(lat,lon), and takes
    the string representation of the start and stop argument to be passed to GOTM.

    The 'start' and 'stop' arguments can be omitted if 'year' and 'month' are specified by keyword.

    All remaining keyword arguments are passed to update GOTM input namelists. 
    """

    ## Preparations.

    from os import symlink
    from os.path import join, isfile, getsize
    from shutil import copyfile
    
    # If 'year' and / or 'month' is specified, extend the positional arguments to match the standard call signature.
    if year is not None:
        from datetime import datetime
        start = datetime(year,1,1,0,0,0) if month is None else datetime(year,month,1,0,0,0)
        stop = datetime(year+1,1,1,0,0,0) if month == 12 or month is None else datetime(year,month+1,1,1,0,0,0)
        args += (start,stop)

    # Now we treat the standard caller with start, stop as the last two arguments.
    if len(args) == 3 or len(args) == 4:
        # Bad arugments will be caught in the next step.
        lat, lon = get_lat_lon(*args[:-2])
        # Passing except the last two arguments to get_local_folder(), which creates the folder if necessary
        # and do the error check on the locations as well.
        local_folder = get_local_folder(*args[:-2], create=False) # This calls get_m_n() too.
        # TODO: Error check for start/stop, which I have a code somewhere in scripts...
        start, stop = args[-2:]
    else:
        print(args)
        raise TypeError('Invalid GOTM run location or start/stop times.')

    # Same latlong could correspond to different (m,n) indices with different subgrids, not so helpful?
    #run_name = 'medsea_GOTM, {:s} grid, #(m,n)=({:d},{:d})'.format(grid,m,n) # debug info

    # Put a tag to generated file names.
    tag = run + '_' + GOTM_version
    out_fn = 'results_' + tag
    daily_stat_fn='daily_stat_'+tag+'.txt'
    
    # Set up GOTM arguments.
    gotm_args = dict(name = tag, # Just use the run name / profile instead of run_name
                     start = str(start), stop = str(stop),
                     latitude = float(lat), longitude = float(lon),
                     out_dir = out_dir, out_fn = out_fn,
                     daily_stat_fn = daily_stat_fn)
    gotm_args.update(**gotm_user_args)

    ## Step 1.
    import shutil
    for each in GOTM_nml:
        # All config files are overwritten every time GOTM is run.
        src = join(GOTM_nml_path,each)
        dst = join(local_folder,each)
        # print('Copying {:s} to {:s}'.format(src,dst))
        copyfile(src,dst)

    # The grid data file too if necessary.
    if 'grid_method' in gotm_user_args.keys() and gotm_user_args['grid_method'] == 2:
        assert 'grid_file' in gotm_user_args.keys()
        # Get the externally specified sea level widths if method is 2.
        grid_file = gotm_user_args['grid_file']
        src = join(GOTM_nml_path,grid_file)
        dst = join(local_folder,grid_file)
        # print('Copying {:s} to {:s}'.format(src,dst))
        copyfile(src,dst)

    ## Step 2.
    # GOTM_executable, default to be at [project_folder]/bin/[GOTM_executable], should have been checked for existence when pyGOTM.config was loaded
    if not(isfile(dst)):
        symlink(join(GOTM_executable_path,GOTM_executable),
                join(local_folder,GOTM_executable))

    ## Step 3.
    # Create a list of dat files we expect to see.
    dat_list = ['tprof','sprof','heat','met'] # Usual set.
    if ('t_prof_method' in gotm_user_args):
        tm = gotm_user_args['t_prof_method']
        if tm != 2: # 2 indicates read from file.
            dat_list.remove('tprof')
    if ('s_prof_method' in gotm_user_args):
        sm = gotm_user_args['s_prof_method']
        if sm != 2: # 2 indicates read from file.
            dat_list.remove('sprof')
    if ('extinct_method' in gotm_user_args):
        em = gotm_user_args['extinct_method']
        if em <= 11 or em == 15 or em == 16:
            pass # Either Jerlov or variants of Paulson-Simpson. Data file not needed.
        elif em == 12:
            dat_list.append('chlo')
        elif em == 13:
            dat_list.append('iop')
        else:
            raise NotImplementedError('Not implemented for extinct_method={:d} yet.'.format(em))

    # Now whether each dat file exists and has nonzero size (not just 'touched' by an erroneous read by GOTM)
    for each in dat_list:
        datfn = join(local_folder,each+'.dat')
        try:
            assert getsize(datfn) > 0
        except:
            raise OSError(each+'.dat not found or has zero length. Check overall grid data integrity. ')

    updatecfg(path=local_folder, **gotm_args)
    return local_folder, out_fn, gotm_args

def buoy_run(code,start,stop,run,cached=True):

    from datetime import datetime
    from netCDF4 import Dataset, num2date
    from os.path import join, isfile, isdir, getsize, basename
    import shutil
    from shutil import copyfile, rmtree
    from glob import glob
    from tempfile import mkdtemp

    if not isinstance(code,str):
        code = str(code)
            
    ### Configure the folder where GOTM is run
    buoy_dir = join(project_folder,'buoys',code)
    if cached:
        dat_src = glob(join(buoy_dir,'*.dat'))
        out_dir = mkdtemp(prefix='medsea_buoy_' + code + '_', dir=cache_folder)
        # Copy the .dat files over
        for src in dat_src:
            fn = basename(src)
            dst = join(out_dir,fn)
            copyfile(src,dst)
    else:
        out_dir = buoy_dir

    # Put a tag to generated file names.
    tag = run + '_' + GOTM_version
    out_fn = 'results_' + tag
    daily_stat_fn='daily_stat_'+tag+'.txt'
        
  
    # Set up GOTM arguments.
    lat = buoys[code].lat
    lon = buoys[code].lon
   
    gotm_args = dict(name = code + '_' + tag, # Just use the run name / profile instead of run_name
                     start = str(start), stop = str(stop),
                     latitude = float(lat), longitude = float(lon),
                     out_dir = out_dir, out_fn = out_fn,
                     daily_stat_fn = daily_stat_fn)

    if run in run_profiles.keys():
        gotm_args.update(run_profiles[run])
    else:
        print('Run profile name given: {:s}'.format(run))
        raise NotImplementedError('A profile must be specified and hard-coded.')

    ### Prepare the run_dir.
    
    ## Step 1.
    for each in GOTM_nml:
        # All config files are overwritten every time GOTM is run.
        src = join(GOTM_nml_path,each)
        dst = join(out_dir,each)
        # print('Copying {:s} to {:s}'.format(src,dst))
        copyfile(src,dst)

    # The grid data file too if necessary.
    if 'grid_method' in gotm_args.keys() and gotm_args['grid_method'] == 2:
        assert 'grid_file' in gotm_args.keys()
        # Get the externally specified sea level widths if method is 2.
        grid_file = gotm_args['grid_file']
        src = join(GOTM_nml_path,grid_file)
        dst = join(out_dir,grid_file)
        # print('Copying {:s} to {:s}'.format(src,dst))
        copyfile(src,dst)

    ## Step 2.
    # GOTM_executable, default to be at [project_folder]/bin/[GOTM_executable], should have been checked for existence when pyGOTM.config was loaded
    if not(isfile(dst)):
        symlink(join(GOTM_executable_path,GOTM_executable),
                join(out_dir,GOTM_executable))

    ## Step 3.
    # Create a list of dat files we expect to see.
    dat_list = ['tprof','sprof','heat','met'] # Usual set.
    if ('t_prof_method' in gotm_args):
        tm = gotm_args['t_prof_method']
        if tm != 2: # 2 indicates read from file.
            dat_list.remove('tprof')
    if ('s_prof_method' in gotm_args):
        sm = gotm_args['s_prof_method']
        if sm != 2: # 2 indicates read from file.
            dat_list.remove('sprof')
    if ('extinct_method' in gotm_args):
        em = gotm_args['extinct_method']
        if em <= 11 or em == 15 or em == 16:
            pass # Either Jerlov or variants of Paulson-Simpson. Data file not needed.
        elif em == 12:
            dat_list.append('chlo')
        elif em == 13:
            dat_list.append('iop')
        else:
            raise NotImplementedError('Not implemented for extinct_method={:d} yet.'.format(em))

    # Now whether each dat file exists and has nonzero size (not just 'touched' by an erroneous read by GOTM)
    for each in dat_list:
        datfn = join(out_dir,each+'.dat')
        try:
            assert getsize(datfn) > 0
        except:
            print(datfn)
            raise OSError(each+'.dat not found or has zero length. Check overall grid data integrity. ')

    # Actually update the namelists: after this, out_dir is really configured to run GOTM.
    updatecfg(path=out_dir, **gotm_args)

    ### Actual run.
    curdir = os.getcwd()
    os.chdir(out_dir)
    stat = dict()
    try:
        tic()

        # 1. Run GOTM keeping the output to log.
        logfn = 'GOTM_' + run + '_' + print_ctime(sep='_') + '.log'
        print('GOTM running... ')
        print('  Working from: {:s}...'.format(out_dir))
        gotm(verbose=verbose,logfn=logfn)
        print('  Results written to {:s}.nc...'.format(join(out_dir,out_fn)))

        # 2. Transfer the results nc and daily_stats txt files immediately.
        if cached:
            def move(src_dir,dst_dir,fn):
                src = join(src_dir,fn)
                dst = join(dst_dir,fn)
                try:
                    print('Moving {:s} to {:s}...'.format(src,dst))
                    shutil.move(src,dst)
                except:
                    raise
            move(out_dir,buoy_dir,out_fn+'.nc')
            move(out_dir,buoy_dir,daily_stat_fn)
        else:
            assert out_dir == buoy_dir
            print('GOTM has been run at {:s} uncached, with I/O at every nsave interval. '.format(buoy_dir))
            

        # 3. Print out GOTM execution statistics.
        stat['elapsed'] = toc()
        # statfn = 'stat_{:d}{:02d}.dat'.format(year,month)
        #statfn = 'STAT_' + run + '_' + print_ctime(sep='_') + '.log'
        statfn = 'STAT_' + run + '.log'
        with open(join(buoy_dir,statfn),'a+') as f:
            print('{:s}: Writing diagnostic statistics to {:s}...\n'.format(print_ctime(sep=' '),statfn))
            f.write('--------------------------------------------------------------\n')
            f.write('Run parameters supplied to override defaults in medsea_GOTM/config/*.inp :\n')
            for key, val in gotm_args.items():
                f.write('    {:s} = {!s}\n'.format(key,val))
            f.write('--------------------------------------------------------------\n')
            f.write('Run statistics:\n')
            f.write('Elapsed: {:.2f} seconds.\n'.format(stat['elapsed']))
            f.write('--------------------------------------------------------------\n')
            f.write('Data statistics:\n')
            with Dataset(join(buoy_dir,gotm_args['out_fn']+'.nc'),'r') as ds:
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
        os.chdir(curdir)
        raise
    
    os.chdir(curdir)
    if cached:
        print("Removing the cached files at {:s}, containing: ".format(out_dir))
        for each in glob(join(out_dir,'*')):
            print('\t {:s}'.format(basename(each)))
        rmtree(out_dir)
    return stat
        
def local_run(*args, # Necessary GOTM run arguments: location i, or (m,n), or (lat,lon); and start/stop,
              run=set_run(), # A name for this run, if it matches a key in run_profile, the settings will be loaded. Should specify when loading the module.
              year = None, month = None, # Provided to compute start, stop.
#              create=False, # Should not run if the folder is not prepared.
              cached = False,
              verbose = verbose, # Print GOTM output to screen? 
              plotvars = None, # Plot results immediately after run? Pass a list of GOTM output variable names.
              **gotm_user_args): # Extra GOTM arguments, e.g. from run_profiles):
    """

    Generate GOTM results for the (m,n)-th grid point at the specified month. Only *.dat files are expected to be
    in the local folder. The config file, GOTM run time will be generated or copied over by calling prepare_run().

    See prepare_run() for argument options.

    """
    from datetime import datetime
    from netCDF4 import Dataset, num2date
    from os.path import join
    import shutil

    if cached:
        out_dir = cache_folder
    else:
        out_dir = '.'
    if run in run_profiles.keys():
        # Should subclass an Exception to tell people what happened.
        if gotm_user_args != {}:
            print(('{:s} = {!s}' * len(gotm_user_args)).format(*(gotm_user_args.items())))
            raise RuntimeError("A recorded run profile {:s} is specified, rejecting all user arguments for GOTM.\n")
        print('Using pre-defined profile: {:s}...'.format(run))
        local_folder, out_fn, gotm_args = prepare_run(*args,run=run,
                                                      year=year, month=month,
                                                      out_dir=out_dir,
                                                      **run_profiles[run])
    else:
        print('Running without a pre-defined profile and tagging results with "{:s}"...'.format(run))
        local_folder, out_fn, gotm_args = prepare_run(*args, run=run,
                                                    year=year, month=month,
                                                      out_dir=out_dir,
                                                      **gotm_user_args)
    os.chdir(local_folder)
    stat = dict()
    try:
        tic()

        # 1. Run GOTM keeping the output to log.
        logfn = 'GOTM_' + run + '_' + print_ctime(sep='_') + '.log'
        print('GOTM running... ')
        print('  Working from: {:s}...'.format(local_folder))
        gotm(verbose=verbose,logfn=logfn)
        print('  Results written to {:s}.nc...'.format(join(out_dir,out_fn)))

        # 2. Transfer the results nc immediately.
        if cached:
            src = join(out_dir,out_fn+'.nc')
            dst = join(local_folder,out_fn+'.nc')
            try:
                print('Moving {:s} to {:s}...'.format(src,dst))
                shutil.move(src,dst)
            except:
                raise

        # 3. Print out GOTM execution statistics.
        stat['elapsed'] = toc()
        # statfn = 'stat_{:d}{:02d}.dat'.format(year,month)
        statfn = 'STAT_' + run + '_' + print_ctime(sep='_') + '.log'
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
            with Dataset(join(local_folder,gotm_args['out_fn']+'.nc'),'r') as ds:
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

    if plotvars:
        from matplotlib.pyplot import subplots
        assert isinstance(plotvars,list)
        ds = Dataset(out_fn+'.nc','r')
        time = num2date(ds['time'][:],ds['time'].units)

        # Read in data for plotting.
        var = dict()
        for varname in plotvars:
            # Some checking
            assert varname in ds.variables, 'Invalid variable specified: {!s}'.format(varname)
            assert len(time) == ds[varname].shape[0], 'Time dimension of variable: {!s} does not match with the time vector.'.format(varname)
            if len(ds[varname].shape) != 3:
                print('Skipping the variable: {!s} that depends on depth (for now...)'.format(varaname))
                continue

            # Now assuming we only plot against time.
            var[varname] = ds[varname][:,0,0]

        fig, axes = subplots(len(var),1,sharex=True)
        for i,(name, val) in enumerate(var.items()):
            ax = axes[i]
            ax.plot(time,val,label=name)
            ax.set_title(name)
            ax.grid('on')
        fig.autofmt_xdate()
        return stat, var # Return the plotted data as well.
    else:
        return stat


#Still used by combine.py and combine_multiple.py, to be deprecated
def create_dimensions(nc, lat=grid_lats, lon=grid_lons):
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

# Taking from scratch notebook. Good to saving buoy information in this module.
class buoy:
    def __init__(self,code,lat,lon,grid='1x'):
        self.code = code
        self.lat = lat
        self.lon = lon
        self.m, self.n = get_m_n(lat,lon,silent=True)
        self.grid_lat, self.grid_lon = get_lat_lon(self.m,self.n,silent=True)
        self.run_GOTM = lambda start, stop: buoy_run(self.code,start,stop,run)
        
# buoy_UTC0 = buoy('61430',39.56,2.1,'1x')
# buoy_UTC1 = buoy('ADN-E2M3A',41.5277,18.0824,'1x')
# buoy_UTC2 = buoy('61277',35.723,25.462,'1x')
# buoys = [buoy_UTC0,buoy_UTC1,buoy_UTC2]

buoys = { '61277' : buoy('61277',35.723,25.462),
          '61280' : buoy('61280',40.6875,1.472),
          '61281' : buoy('61281',39.5215,0.2075),
          '61430' : buoy('61430',39.555,2.105),
          '68422' : buoy('68422',36.829,21.608),
          'SARON' : buoy('SARON',37.61,23.569)
}

## Post-run helper functions.

# def get_daily_stats_by_year(year=2014):
#     from os.path import join
#     from netCDF4 import Dataset
#     return Dataset(join(results_folder,run_key,
#                         'daily_stat_{!s}_{!s}_{!s}_{!s}.nc'.format(run_key,grid,GOTM_version,year)))
# def get_hourly_results_by_year(year=2014):
#     from os.path import join
#     from netCDF4 import Dataset
#     return Dataset(join(results_folder,
#                         'medsea_GOTM_{!s}_{!s}_{!s}_{!s}.nc'.format(run_key,grid,GOTM_version,year)))   

## Rewrite and put in another file, not here.
# def local_dat(mm,nn,dat=['heat','met','tprof','sprof','chlo','iop']):
#     """
#     Generate *.dat files from all available data. See core_dat() for other explanations.
#     mm, nn can be a sequence of m,n's
#     """

#     from netCDF4 import MFDataset
#     from tempfile import mkdtemp
#     from shutil import copyfile
#     from glob import glob

#     if isinstance(dat,str):
#         dat = [dat]

# # 20170621, we no longer create separate folders for different runs, but share the same set of subfolders named
# # by print_lat_lon(), to avoid creating too many files and draining disk quota too fast.
# #    run_folder = os.path.join(base_folder,run)
#     run_folder = base_folder

# #    if not(os.path.isdir(run_folder)):
# #        os.mkdir(run_folder)

#     print("Looking for data sources from " + data_folder + "...")
#     nc_dict = data_sources(dat=dat)
#     if not(isinstance(nc_dict, dict)):
#         nc_dict = {dat[0]: nc_dict}
#     print(nc_dict)

#     for dat_fn, nc in nc_dict.items():
#         ## Assume m, n are iterable and of the same length:
#         if isinstance(mm,int) and isinstance(nn,int):
#             mm = [mm]
#             nn = [nn]
#         else:
#             assert len(m) == len(n)
#         try:
#             for m, n in zip(mm, nn):
#                 lat = grid_lats[m]
#                 lon = grid_lons[n]

#                 # Create the local folder if necessary.
#                 local_folder = os.path.join(run_folder,print_lat_lon(lat,lon))
#                 if not(os.path.isdir(local_folder)):
#                     os.mkdir(local_folder)

#                 # Actually write.
#                 write_dat(m,n,dat_fn,nc,local_folder)
#         except Exception:
#                 print('mm',mm)
#                 print('nn',nn)

#         ## When m,n are integers, exception occurs at zip() instead, and it not a TypeError.
#         # except TypeError as te:
#         #     if isinstance(mm,int) and isinstance(nn,int):
#         #         m = mm
#         #         n = nn
#         #         # So the provided sequences are just a single pair of indices. Just write.
#         #         write_dat(m,n,dat_fn,nc,local_folder)
#         #    else:

# #        write_dat(m,n,dat_fn,nc,local_folder)

#     for nc in nc_dict.values():
#         nc.close()
