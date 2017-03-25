### Common settings
import sys, os
home = os.getenv('HOME')
project_folder = os.path.join(home,'gotm-dst')

### Global settings
output_folder = os.path.join('/scratch','simontse') 
GOTM_lat = tuple(30.75+0.75*i for i in range(21))
GOTM_lon = tuple(-6.0+0.75*i for i in range(57))
print("Output folder: ", output_folder)
print("Latitudes: ", GOTM_lat)
print("Longitudes: ",GOTM_lon)

## ERA-INTERIM specific settings
ERA_folder = os.path.join(project_folder,'p_sossta','medsea_ERA-INTERIM','3-hourly')
ERA_lat_ind = slice(-8,4,-1)
ERA_lon_ind = slice(12,-18,1)

met_names = ['u10m','v10m','sp','t2m','q2m','precip','snow']
met_alias = ['u10m','v10m','sp','t2m','q2m','var144','var228']
heat_names = ['lwrd','swrd']
heat_alias = ['var175','var169']

ERA_names = met_names + heat_names
ERA_alias = met_alias + heat_alias

### Helper functions
def get_ERA_yearly_data(folder,year,name,alias,lat_indices,lon_indices):
    from netCDF4 import Dataset
    fn = 'MEDSEA_ERA-INT_' + name + '_y' + str(year) + '.nc'
    with Dataset(os.path.join(folder,fn),'r') as nc:
        # First, confirm every time that the indices for lat and lon are correct.
        assert all(nc['lat'][ERA_lat_ind] == GOTM_lat)
        assert all(nc['lon'][ERA_lon_ind] == GOTM_lon)
        # Then return the data unpacked from the netCDF object.
        return nc[alias][:,lat_indices,lon_indices]

def create_dimensions(nc, epoch, lat=GOTM_lat, lon=GOTM_lon):
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

def ERA_reformat(year):
    " Reformat the ERA data by combining variables needed for met.dat into monthly files. "
    import os
    from netCDF4 import Dataset
    from datetime import datetime, timedelta
    
    # The ERA yearly dataset has epoch at the start of that year... 
    # We unify them to 1980-01-01 00:00:00.
    epoch = datetime(1980,1,1,0,0,0)
    #data = {name: get_ERA_yearly_data(ERA_folder,year,name,alias) \
    #        for name,alias in zip(ERA_names,ERA_alias)}

    def timings(month):
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
        if month == 1:
            # First record at 03 hour on Jan 1st.
            start_hour = int((this_month + timedelta(hours=3) - epoch).total_seconds()/3600)
        else:
            # First records at 00 hour on 1st of Feb, Mar, Apr ..., Dec.
            start_hour = int((this_month - epoch).total_seconds()/3600) 
        # Last GOTM record time.
        end_hour = int((next_month - epoch).total_seconds()/3600)

        ## Find the number of days since the first day of year to infer indices of the records we want.        
        first_day_of_year = datetime(this_month.year,1,1,0,0,0)
        start_day = (this_month-first_day_of_year).days
        end_day = (next_month-first_day_of_year).days
        
        ## Calculate indices of time to read in ERA data.
        #
        # Include one previous record at 00:00:00 unless at beginning of year.
        start_ind = start_day*8 if month == 1 else start_day*8-1
        # Total number of records should be (end_day-start_day)*8
        end_ind = end_day*8
        # Make doubly sure there's no time-shift by 3-hour periods!
        if month == 1:
            assert end_ind - start_ind == nper
        else:
            assert end_ind - start_ind == nper + 1

        return start_hour, end_hour, start_ind, end_ind

    # Create the monthly files first with the basic dimensions.
    for month in range(1,13):
        # Output filename and full path.
        outfn = 'medsea_ERA_{0:d}{1:02d}.nc'.format(year,month)
        fullfile = os.path.join(output_folder,'medsea_ERA-INTERIM',outfn)
        # Interpretation of ERA data timings CRITICAL here to get these indices correct.
        start_hour, end_hour, start_ind, end_ind = timings(month)
        with Dataset(fullfile,'w') as nc:
            # Routine stuff delegated to a helper funciton.
            nctime, nclat, nclon = create_dimensions(nc,epoch)
            # Write the time values to the nc file.
            nctime[:] = [hour for hour in range(start_hour, end_hour+3,3)]
        print(fullfile + ' created with dimensions set up.')
            
    # Iterate through the variables and append the data.
    for name,alias in zip(ERA_names,ERA_alias):
        # Fetch the yearly data in one go.
        data = get_ERA_yearly_data(ERA_folder,year,name,alias,ERA_lat_ind,ERA_lon_ind)

        for month in range(1,13):
            # Output filename and full path.
            outfn = 'medsea_ERA_{0:d}{1:02d}.nc'.format(year,month)
            fullfile = os.path.join(output_folder,'medsea_ERA-INTERIM',outfn)
            # Append the new variable.
            with Dataset(fullfile,"a") as nc:
                ncvar = create_variable(nc,name,'f8')
                ncvar[:] = data[start_ind:end_ind,:,:]
        print('Done copying over {} values'.format(name))
        
