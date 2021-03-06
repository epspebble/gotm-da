### Common settings
import sys, os
from .config import *
from .gotmks import *
from .medsea import *

## These global settings should now in pygotm.config, imported above.

# ### Global settings
# data_folder = os.path.join('/global/scratch',os.getenv('USER'))
# medsea_lats = tuple(30.75+0.75*i for i in range(21))
# medsea_lons = tuple(-6.0+0.75*i for i in range(57))
# print("Output folder: ", data_folder)
# print("Latitudes: ", medsea_lats)
# print("Longitudes: ",medsea_lons)
#
# ## ERA-INTERIM specific settings
# ERA_folder = os.path.join(data_folder,'p_sossta','medsea_ERA-INTERIM','3-hourly')

# The corresponding index ranges in medsea_ERA-INTERIM datasets.
medsea_ERA_lat_ind = slice(-8,4,-1)
medsea_ERA_lon_ind = slice(12,-18,1)

# The corresponding index ranges in medsea_rea datasets.
medsea_rea_lat_ind = slice(9,250,12)
medsea_rea_lon_ind = slice(0,673,12)

# In medsea_rea ocean reanalysis data:
# depth[24] = 176.82929993, we want up to 150m only.
# We used depth[18] = 101.78025055, when our max_depth was 100m.
#ndepth = 18
ndepth = 24 # Use this for any region? Double check!
rea_depth_ind = slice(0,ndepth)

met_names = ['u10m','v10m','sp','t2m','q2m','precip','snow']
met_alias = ['u10m','v10m','sp','t2m','q2m','var144','var228']
heat_names = ['lwrd','swrd']
heat_alias = ['var175','var169']

# ERA_names = met_names + heat_names
# ERA_alias = met_alias + heat_alias

def timings(year,month):
    from datetime import datetime,timedelta
    ### MAIN WORK: calculating various date and time offsets.
    
    ### INTERPRETATION of ERA-INTERIM data source.
    # In ERA-Interim (3-hourly), u10m, v10m, sp, t2m, q2m are instantaneous values at the end
    # of each period, while lwrd, swrd are mean values over each period.
    # NOTE: To change less Fortran, we opt to move the time stamps of the lwrd and swrd values
    # 3 hours early, so the the values "persist" for the whole 3-hourly period.
    ###

    ### GOTM timing of variables.

    ## Find no. of 3-hourly periods.
    # Midnight of first day of this month
    this_month = datetime(year,month,1,0,0,0)
    # Midnight of first day of next month
    next_month = datetime(year,month+1,1,0,0,0) if month < 12 else datetime(year+1,1,1,0,0,0)
    # No. of 3-hourly periods, must be integer, converting to int explicitly.
    nper = int((next_month - this_month).total_seconds() / (3*3600))
        
    ## Calculate the time variables start and end of GOTM records to write.
    #
    # 'start_sec' time value of our first GOTM record.
        
    ## WT 20170327 This extra padding of 00 hour (i.e. the instananeous data for met variables, but 21:00:00-00:00:00 
    # the previous day for heat variables is causing debugging trouble. Let's just skip it for now, and let GOTM 
    # handle missing data at start?
    #if month == 1:
    #    # First record at 3 hours on Jan 1st.
    #    start_hour = int((this_month + timedelta(hours=3) - epoch).total_seconds()/3600)
    #else:
    #    # First records at 00 hour on 1st of Feb, Mar, Apr ..., Dec.
    #    start_hour = int((this_month - epoch).total_seconds()/3600) 
        
    # A faithful copy of original ERA data, do not make changes. Just recalculate the time values of 3, 6, 9 hrs etc...
    start_hour = int((this_month + timedelta(hours=3) - epoch).total_seconds()/3600)
            
    # Last GOTM record time.
    end_hour = int((next_month - epoch).total_seconds()/3600)

    ## Find the number of days since the first day of year to infer indices of the records we want.        
    first_day_of_year = datetime(this_month.year,1,1,0,0,0)
    start_day = (this_month-first_day_of_year).days
    end_day = (next_month-first_day_of_year).days
        
    ## Calculate indices of time to read in ERA data.
    #
    ## WT 2017037 Do not include previous record. Debugging hazard.
    ## Include one previous record at 00:00:00 unless at beginning of year.
    #start_ind = start_day*8 if month == 1 else start_day*8-1
        
    # A faithful copy of original ERA data.
    start_ind = start_day*8
        
    # Total number of records should be (end_day-start_day)*8
    end_ind = end_day*8
        
    ## WT 2017037 Skip this hassle!
    ## Make doubly sure there's no time-shift by 3-hour periods!
    #if month == 1:
    #    assert end_ind - start_ind == nper
    #else:
    #    assert end_ind - start_ind == nper + 1
        
    assert end_ind - start_ind == nper

    return start_hour, end_hour, start_ind, end_ind


def ERA_reformat(year,region='medsea'):
    " Reformat the ERA data by combining variables needed for met.dat into monthly files. "
    import os
    from time import time
    from netCDF4 import Dataset
    from datetime import datetime, timedelta
    
    ## 'epoch' moved to global config, load it by import .config
    # The ERA yearly dataset has epoch at the start of that year... 
    # We unify them to 1981-01-01 00:00:00, which seems to be the convention used in UKMO and
    # Copernicus analysed SST products.
    # epoch = datetime(1981,1,1,0,0,0)
    # data = {name: get_ERA_yearly_data(ERA_folder,year,name,alias) \
    #        for name,alias in zip(ERA_names,ERA_alias)}

    # Create the monthly files first with the basic dimensions.
    #for month in range(1,13):
    #    # Output filename and full path.
    #    outfn = region+'_ERA_{0:d}{1:02d}.nc'.format(year,month)
    #    fullfile = os.path.join(data_folder,region+'_ERA-INTERIM',outfn)

    def get_ERA_yearly_data(name,alias,region='medsea'):
        if region == 'medsea':
            fn = 'MEDSEA_ERA-INT_' + name + '_y' + str(year) + '.nc'
            with Dataset(os.path.join(ERA_folder,fn),'r') as nc:
                # First, confirm every time that the indices for lat and lon are correct.
                assert all(nc['lat'][medsea_ERA_lat_ind] == medsea_lats)
                assert all(nc['lon'][medsea_ERA_lon_ind] == medsea_lons)
                # Then return the data unpacked from the netCDF object.
                data = nc[alias][:,medsea_ERA_lat_ind,medsea_ERA_lon_ind]
        else:
            raise NotImplementedError('The specified region: ' + region + ' is not implemented for yet.')
        return data
            
    # Iterate through the variables and append the data.
    def write_ERA(ERA_names, ERA_alias, subtype):
        for i,(name,alias) in enumerate(zip(ERA_names,ERA_alias)):
            # Fetch the yearly data in one go.
            tic()
            print("Reading in {:d}'s ERA-INTERIM data for {:s}...".format(year,name))
            data = get_ERA_yearly_data(name,alias)
            toc()
        
            # Output monthly and yearly files.
            elapsed = 0
            for month in range(1,13):
                tic()
                print("Writing {1:s} data for the month #{0:d}...".format(month,name))
                # Output filename and full path.
                outfn = region+'_ERA_{2:s}_{0:d}{1:02d}.nc'.format(year,month,subtype)
                fullfile = os.path.join(data_folder,region+'_ERA-INTERIM',outfn)
                # Interpretation of ERA data timings CRITICAL here to get these indices correct.
                start_hour, end_hour, start_ind, end_ind = timings(year,month)

                # Once: open nc file and create dimensions.
                if i == 0: # At the first variable only.
                    # Start the nc file with the basic dimensions, overwriting existing file.
                    with Dataset(fullfile,'w',format="NETCDF3_CLASSIC") as nc:
                        # Routine stuff delegated to a helper funciton.
                        nctime, nclat, nclon = create_dimensions(nc)
                        # Write the time values to the nc file.
                        nctime[:] = [hour for hour in range(start_hour, end_hour+3,3)]
                        #print(fullfile + ' with time({}), lat({}), lon({}) set up.'.format(len(nctime),len(nclat),len(nclon)))
                
                # Create each variable and append values.
                with Dataset(fullfile,"a") as nc:
                    ncvar = create_variable(nc,name,'f8')
                    ncvar[:] = data[start_ind:end_ind,:,:]
                elapsed += toc()
            print("Total time for writing {:s}: {:2f}s".format(name,elapsed))

    write_ERA(heat_names, heat_alias,'heat')
    write_ERA(met_names, met_alias,'met')

def medsea_rea_reformat(year,month,varname,fn_keyword):
    from netCDF4 import Dataset, MFDataset, num2date, date2num
    from numpy import linspace
    from os.path import isfile

    output_folder = os.path.join(data_folder,'medsea_rea')
    outfn = os.path.join(output_folder,'medsea_rea_{2}_{0:d}{1:02d}.nc'.format(year,month,varname,output_folder))
                       
    with Dataset(outfn,"w",format='NETCDF3_CLASSIC') as nc:
        nctime, nclat, nclon = create_dimensions(nc)
        nc.createDimension('depth', size = ndepth) 
        ncdepth = nc.createVariable('depth','f4',dimensions=('depth',))
        ncvar = create_variable(nc,varname,'f8',dimensions=('time','depth','lat','lon'))
        print('Done initializing dimensions and variables.')

        # Create the variable in question.
        infn = "{0:d}{1:02d}??_{2}_re-fv6.nc".format(year,month,fn_keyword)
        with MFDataset(os.path.join(rea_folder,str(year),infn),'r') as rea_data:
            # Copy over the depths up to the truncation level.
            ncdepth[:] = rea_data['depth'][0:ndepth]
            ncdepth.units = rea_data['depth'].units
            # Write in the time values, convert to our epoch and units.
            nctime[:] = date2num(num2date(rea_data['time'][:],rea_data['time'].units),'hours since '+str(epoch))
            # Copy over variable data.
            temp = rea_data[varname][:,rea_depth_ind,medsea_rea_lat_ind,medsea_rea_lon_ind]
            ncvar[:] = temp # Now random-access in RAM to copy over.
            ncvar.units = rea_data[varname].units
            print('Done copying over {} values of one month.'.format(varname))

    # Also copy over the data for the first day of the next month.
    if month == 12:
        year_new = year+1
        month_new = 1
    else:
        year_new = year
        month_new = month+1
    infn = "{0:d}{1:02d}01_{2}_re-fv6.nc".format(year_new,month_new,fn_keyword)
    if os.path.isfile(infn): # Check availability of data, e.g. cannot do 2014-12 without data for 2015-01
        with Dataset(outfn,"a") as nc:
            ncvar = nc[varname]
            with Dataset(os.path.join(rea_folder,str(year_new),infn),'r') as rea_data:
                print('Before appending first day value of the next month, {0}, has dimensions {1},'.format(varname,ncvar[:].shape))
                print('whereas time has length {}'.format(len(nc['time'][:])))
                last = len(nc['time'])
                nc['time'][last] = date2num(num2date(rea_data['time'][0],rea_data['time'].units),'hours since '+str(epoch))
                ncvar[last,:,:,:] = rea_data[varname][0,rea_depth_ind,medsea_rea_lat_ind,medsea_rea_lon_ind]
                print('After appending, the dimensions of {0} become {1},'.format(varname,ncvar[:].shape))
                print('whereas time has length {}'.format(len(nc['time'][:])))
            
    print("Done creating {}".format(outfn))

# New functions
def get_REA_grid():
    with Dataset('/global/scratch/simontse/p_sossta/medsea_rea/2013/20130101_TEMP_re-fv6.nc','r') as nc:
        lat_rea = nc['lat'][:]
        lon_rea = nc['lon'][:]
        temp_rea = nc['votemper'][:]
        is_sea = ~temp_rea[0,0,:].mask
    return lat_rea, lon_rea, is_sea

def get_ERA_yearly_data(name,alias,year=2013,region='medsea'):
    if region == 'medsea':
        fn = 'MEDSEA_ERA-INT_' + name + '_y' + str(year) + '.nc'
        with Dataset(os.path.join(ERA_folder,fn),'r') as nc:
            # These indices are selected to cover the rea grid just right.
            lat = nc['lat'][-7:3:-1]
            lon = nc['lon'][12:-17]
            data = nc[alias][:,-7:3:-1,12:-17] # Now we want to cover the rea grid
    else:
        raise NotImplementedError('The specified region: ' + region + ' is not implemented for yet.')
    return lat, lon, data

def write_ERA_on_rea_grid(ERA_names, ERA_alias, subtype, 
                          year=2013, months=range(1,13),
                          format="NETCDF4"):
    from scipy.interpolate import RectBivariateSpline
    from numpy.ma import array
    region = 'medsea'
    lat, lon, is_sea = get_REA_grid() 

    for i,(name,alias) in enumerate(zip(ERA_names,ERA_alias)):
        # Fetch the yearly data in one go.
        tic()
        print("Reading in {:d}'s ERA-INTERIM data for {:s}...".format(year,name))
        lat_ERA, lon_ERA, data = get_ERA_yearly_data(name,alias,year=year)
        toc()

        # Output monthly and yearly files.
        elapsed = 0
        for month in months:
            tic()
            print("Writing 144x (full rea grid) interpolated {1:s} data for the month #{0:d}...".format(month,name))
            # Output filename and full path.
            outfn = region+'_ERA_{2:s}_{0:d}{1:02d}.nc'.format(year,month,subtype)
            fullfile = os.path.join(data_folder,region+'_ERA-INTERIM_144x',outfn)
            # Interpretation of ERA data timings CRITICAL here to get these indices correct.
            start_hour, end_hour, start_ind, end_ind = timings(year,month)

            # Once: open nc file and create dimensions.
            if i == 0: # At the first variable only.
                # Start the nc file with the basic dimensions, overwriting existing file.
                with Dataset(fullfile,'w',format=format) as nc:
                    # Routine stuff delegated to a helper funciton.
                    nctime, nclat, nclon = create_dimensions(nc,lat=lat,lon=lon)
                    # Write the time values to the nc file.
                    nctime[:] = [hour for hour in range(start_hour, end_hour+3,3)]
                    #print(fullfile + ' with time({}), lat({}), lon({}) set up.'.format(len(nctime),len(nclat),len(nclon)))

            # Create each variable and append values.
            with Dataset(fullfile,"a") as nc:
                ncvar = create_variable(nc,name,'f8',zlib=True)
                data_slice = data[start_ind:end_ind,:]
                data_intp = array([RectBivariateSpline(lat_ERA,lon_ERA,data_slice[i,:])(lat,lon) for i in range(data_slice.shape[0])])
                ncvar[:] = data_intp[:]
            elapsed += toc()
        print("Total time for writing {:s}: {:2f}s".format(name,elapsed))

def write_SWR_CS_on_rea_grid(year,month):
    from pygotm.medsea import tic, toc
    from netCDF4 import Dataset
    from numpy.ma import array
    from numpy import load
    from scipy.interpolate import RectBivariateSpline
    swr_cs = load('swr_cs.npy')

    print('Interpolating the swrd_clear_sky values for {:d}-{:02d}'.format(year,month))
    tic()
    def swr_cs_interp(year,month):
        with Dataset('/home/simontse/scratch/medsea_ERA-INTERIM_20170331_withCF/medsea_ERA_{:d}{:02d}.nc'.format(year,month),'r') as ds:
            lat_ERA = ds['lat'][:]
            lon_ERA = ds['lon'][:]
            num = (year-2013)*12 + month-1
            vals = array([RectBivariateSpline(lat_new,lon_new,swr_cs[num][i,:])(lat_rea,lon_rea) for i in range(swr_cs[num].shape[0])]) 
        return vals

    with Dataset('/home/simontse/scratch/medsea_ERA-INTERIM_144x/medsea_ERA_heat_{:d}{:02d}.nc'.format(year,month),'a') as ds:
        new_data = swr_cs_interp(year,month)
        if 'swrd_clear_sky' not in ds.variables.keys():
            new_var = ds.createVariable('swrd_clear_sky','f8',dimensions=('time','lat','lon'),zlib=True,fill_value=1e20)
        else:
            new_var = ds['swrd_clear_sky']
        new_var[:] = new_data
    toc()

def write_CF_on_rea_grid(year,month):
    from pygotm.medsea import tic, toc
    from numpy.ma import divide, masked_outside
    print('Calcuating cloud_factor values for {:d}-{:02d}'.format(year,month))
    tic()
    fn = '/home/simontse/scratch/medsea_ERA-INTERIM_144x/medsea_ERA_heat_{:d}{:02d}.nc'.format(year,month)
    with Dataset(fn,'a') as ds:
        if 'swrd_clear_sky' not in ds.variables.keys():
            raise Exception('The variable swrd_clear_sky not found in ' + fn)
        else:
            swr_cs = ds['swrd_clear_sky'][:]
            swr_era = ds['swrd'][:]
        if 'cloud_factor' not in ds.variables.keys():
            cf = ds.createVariable('cloud_factor','f8',dimensions=('time','lat','lon'),zlib=True,fill_value=1e20)
        else:
            cf = ds['cloud_factor']
        cf[:] = masked_outside(divide(swr_era,swr_cs),0,1)
    toc()

        
# if __name__ == '__main__':
#     if len(sys.argv) == 1: # No arguments provided, assume 2013, 2014.
#         years = [2013, 2014]
#     else:
#         years = [int(sys.argv[i+1]) for i in range(len(argv)-1)]
#     for year in years:
#         ERA_reformat(year)
#         for i in range(24):
#             rea_reformat(2013+int(i/12),1+i%12,'votemper','TEMP')
#             rea_reformat(2013+int(i/12),1+i%12,'vosaline','PSAL')
