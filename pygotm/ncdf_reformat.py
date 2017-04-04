### Common settings
import sys, os
from .config import *
from .gotmks import *

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
ERA_lat_ind = slice(-8,4,-1)
ERA_lon_ind = slice(12,-18,1)

# The corresponding index ranges in medsea_rea datasets.
rea_lat_ind = slice(9,250,12)
rea_lon_ind = slice(0,673,12)
rea_depth_ind = slice(0,18)


met_names = ['u10m','v10m','sp','t2m','q2m','precip','snow']
met_alias = ['u10m','v10m','sp','t2m','q2m','var144','var228']
heat_names = ['lwrd','swrd']
heat_alias = ['var175','var169']

ERA_names = met_names + heat_names
ERA_alias = met_alias + heat_alias


def get_ERA_yearly_data(folder,year,name,alias,lat_indices,lon_indices):
    from netCDF4 import Dataset
    fn = 'MEDSEA_ERA-INT_' + name + '_y' + str(year) + '.nc'
    with Dataset(os.path.join(folder,fn),'r') as nc:
        # First, confirm every time that the indices for lat and lon are correct.
        assert all(nc['lat'][ERA_lat_ind] == medsea_lats)
        assert all(nc['lon'][ERA_lon_ind] == medsea_lons)
        # Then return the data unpacked from the netCDF object.
        return nc[alias][:,lat_indices,lon_indices]

def timings(month):
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


def ERA_reformat(year):
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
    #    outfn = 'medsea_ERA_{0:d}{1:02d}.nc'.format(year,month)
    #    fullfile = os.path.join(data_folder,'medsea_ERA-INTERIM',outfn)
            
    # Iterate through the variables and append the data.
    for i,(name,alias) in enumerate(zip(ERA_names,ERA_alias)):
        # Fetch the yearly data in one go.
        tic()
        print("Reading in {:d}'s ERA-INTERIM data for {:s}...".format(year,name))
        data = get_ERA_yearly_data(ERA_folder,year,name,alias,ERA_lat_ind,ERA_lon_ind)
        toc()
        
        # Output monthly and yearly files.
        for month in range(1,13):
            tic()
            print("Writing data for the month #{:d}...".format(month))
            # Output filename and full path.
            outfn = 'medsea_ERA_{0:d}{1:02d}.nc'.format(year,month)
            fullfile = os.path.join(data_folder,'medsea_ERA-INTERIM',outfn)
            # Interpretation of ERA data timings CRITICAL here to get these indices correct.
            start_hour, end_hour, start_ind, end_ind = timings(month,epoch)

            if i == 0: # Creating the first variable.
                # Start the nc file with the basic dimensions, overwriting existing file.
                with Dataset(fullfile,'w',format="NETCDF3_CLASSIC") as nc:
                    # Routine stuff delegated to a helper funciton.
                    nctime, nclat, nclon = create_dimensions(nc,epoch)
                    # Write the time values to the nc file.
                    nctime[:] = [hour for hour in range(start_hour, end_hour+3,3)]
                    #print(fullfile + ' with time({}), lat({}), lon({}) set up.'.format(len(nctime),len(nclat),len(nclon)))
                
            # Append the new variable.
            with Dataset(fullfile,"a") as nc:
                ncvar = create_variable(nc,name,'f8')
                ncvar[:] = data[start_ind:end_ind,:,:]
            toc()

def rea_reformat(year,month,varname,fn_keyword):
    from netCDF4 import Dataset, MFDataset
    from numpy import linspace
    from os.path import isfile
    output_folder = os.path.join(data_folder,medsea_rea)
 
    outfn = os.path.join(output_folder,'medsea_rea_{2}_{0:d}{1:02d}.nc'.format(year,month,varname,output_folder))
                         
    with Dataset(outfn,"w",data_model='NETCDF3_CLASSIC') as nc:
        nctime, nclat, nclon = create_dimensions(nc,epoch)
        nc.createDimension('depth', size = 18) # 18 was found manually to get just over 100m.
        ncdepth = new.createVariable('depth','f4',dimensions=('depth',))
    print('Done creating dimensions.')
                
    with MFDataset(os.path.join(rea_folder,"{0:d}{1:02d}??_{2}_re-fv6.nc".format(year,month,fn_keyword)),'r') as REA_data:
        ncdepth.units = REA_data['depth'].units
        rea_lat_ind = slice(9,250,12)
        rea_lon_ind = slice(0,673,12)
        rea_depth_ind = slice(0,18)
    
        temp = REA_data[varname][:,0:18,:,:] # Read continguously first.
        ncvar[:] = temp[:,:,rea_lat_ind,rea_lon_ind] # Now random-access in RAM to copy over.
        ncvar.units = REA_data[varname].units
        print('Done copying over {} values of one month.'.format(varname))
            
### STOPPED HERE.
        # Also copy over the data for the first day of the next month.
        if month == 12:
            year_new = year+1
            month_new = 1
        else:
            year_new = year
            month_new = month+1
        
        # Cannot do 2014-12 without data for 2015-01
        if (year == 2014) and (month == 12):
            return
        
        with Dataset("{3}/{0:d}{1:02d}01_{2}_re-fv6.nc".format(year_new,month_new,fn_keyword,REA_folder),'r') as REA_data:
            print('Before appending first day value of the next month, {0} has dimensions {1}'.format(varname,shape(nc_var)))
            print('Meanwhile, time has length {}'.format(len(new['time'][:])))
            last = len(new['time'])
            new['time'][last] = REA_data['time'][:]
            temp = REA_data[varname][:,0:18,:,:]
            nc_var[last,:,:,:] = temp[:,:,9:250:12,0:673:12]
            print('After appending, the dimensions of {0} become {1}'.format(varname,shape(nc_var)))
            print('Meanwhile, time has length {}'.format(len(new['time'][:])))
            
    print("Done creating {}".format(newfn))

        
if __name__ == '__main__':
    if len(sys.argv) == 1: # No arguments provided, assume 2013, 2014.
        years = [2013, 2014]
    else:
        years = [int(sys.argv[i+1]) for i in range(len(argv)-1)]
    for year in years:
        ERA_reformat(year)
