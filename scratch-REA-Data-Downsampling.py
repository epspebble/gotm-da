
# coding: utf-8

# In[1]:

get_ipython().magic('pylab inline')


# In[2]:

from netCDF4 import Dataset, MFDataset


# In[3]:

REA_folder = '/home/simon/gotm-dst/medsea_rea'


# In[4]:

get_ipython().magic('ls ../medsea_ERA-INTERIM')


# In[5]:

# The following has been checked to be equal to the data in ERA.
ERA_lat = linspace(49.5,25.5,33) # Evenly spaced vertical grid lines.
ERA_lon = linspace(-15.,49.5,87) # Evenly spaced horizontal grid lines. 
print(ERA_lat)
print(ERA_lon)

# Observe that spacing is 0.75 in both directions.


# In[13]:

REA_fn = "{}/20140201_TEMP_re-fv6.nc".format(REA_folder)
with Dataset(REA_fn,'r') as REA_data:
    REA_lat = REA_data['lat'][:]
    REA_lon = REA_data['lon'][:]
print('REA grid bounding box, {0} <= lat <= {1}, {2} <= lon <= {3}'.format(
        min(REA_lat),max(REA_lat),min(REA_lon),max(REA_lon)))

lat_overlap = set(REA_lat).intersection(set(ERA_lat))
lon_overlap = set(REA_lon).intersection(set(ERA_lon))

lat_min = min(lat_overlap)
lat_max = max(lat_overlap)
lon_min = min(lon_overlap)
lon_max = max(lon_overlap)
lat_n = int((lat_max-lat_min)/0.75+1)
lon_n = int((lon_max-lon_min)/0.75+1)

print('ERA `lat` grid lines subset that overlaps with that of REA:',lat_min,lat_max)
print('ERA `lon` grid lines subset that overlaps with that of REA:',lon_min,lon_max)
print('Reconstructing a sub-grid of ERA with {0:.0f} x {1:.0f} (lat,lon) grid points.'.format(lat_n,lon_n))

GOTM_lat = linspace(lat_min,lat_max,lat_n)
GOTM_lon = linspace(lon_min,lon_max,lon_n)

print('GOTM_lat:',GOTM_lat)
print('GOTM_lon:',GOTM_lon)


# In[7]:

# Index offsets
ind = dict()
ind['lat_min'] = (min(REA_lat)-lat_min)/0.0625
ind['lat_max'] = (max(REA_lat)-lat_max)/0.0625
ind['lon_min'] = (min(REA_lon)-lon_min)/0.0625
ind['lon_max'] = (max(REA_lon)-lon_max)/0.0625

# Step size
ind['step'] = 0.75/0.0625

# Confirm that the indices are integers and convert datatype.
def chk_int(num):
    assert(mod(num,1) == 0.)
    return int(num)
for key in ind.keys():
    ind[key] = chk_int(ind[key])
    
# Generate subindices for later use when downsampling.
REA_lat_subindices = range(-ind['lat_min'],len(REA_lat)-ind['lat_max'],ind['step'])
REA_lon_subindices = range(-ind['lon_min'],len(REA_lon)-ind['lon_max'],ind['step'])

# Double check that we will be getting the same grid.
assert(all(GOTM_lat == REA_lat[REA_lat_subindices]))
assert(all(GOTM_lon == REA_lon[REA_lon_subindices]))

print("Subindices for REA dataset for the grid subset computed.")

print('REA_lat_subindices: ', REA_lat_subindices)
print('REA_lon_subindices: ', REA_lon_subindices)


# In[8]:

with Dataset("{}/20140202_TEMP_re-fv6.nc".format(REA_folder),'r') as REA_data:
    print(REA_data['depth'][0:18])


# In[26]:

def downsampling_monthly(inputs,override = False):
    from netCDF4 import Dataset, MFDataset
    from numpy import linspace
    from os.path import isfile
    year,month,varname,fn_keyword = inputs
    print("Inputs given: ", inputs)
    REA_folder = '/home/simon/gotm-dst/medsea_rea'
    output_folder = "/home/simon/gotm-dst/medsea_GOTM/profiles"
    newfn = "{3}/medsea_rea_{2}_{0:d}{1:02d}.nc".format(year,month,varname,output_folder);
    if isfile(newfn) and not(override):
        print('File {} is found. Skipping...'.format(newfn))
        return
    with Dataset(newfn,"w",data_model='NETCDF4_CLASSIC') as new:
        new.createDimension('time') # unlimited
        new.createDimension('depth', size = 18) # 18 was found manually to get just over 100m.
        new.createDimension('lat', size = 21) # use the same as in ERA dataset
        new.createDimension('lon', size = 57) # use the same as in ERA dataset
        print('Done creating dimensions.')

        nc_time = new.createVariable('time','f8',dimensions=('time',))
        nc_depth = new.createVariable('depth','f8',dimensions=('depth',))
        
        nc_lat = new.createVariable('lat','f8',dimensions=('lat',))
        nc_lon = new.createVariable('lon','f8',dimensions=('lon',))
        nc_var = new.createVariable(varname,'f8',dimensions=('time','depth','lat','lon'),
                           zlib=True, fill_value=1e+20) # Fill-value the same as in REA dataset.       
        print('Done creating variables.')
        
        # A subset of ERA grid that is within REA's grid bbox.
        nc_lat[:] = linspace(30.75,45.75,21)
        nc_lon[:] = linspace(-6.0,36.0,57)
        print('Done creating lat/lon grid.')
        
        with MFDataset("{3}/{0:d}{1:02d}??_{2}_re-fv6.nc".format(year,month,fn_keyword,REA_folder),'r') as REA_data:
            nc_lat.units = REA_data['lat'].units
            nc_lon.units = REA_data['lon'].units
            nc_time.units = REA_data['time'].units
            nc_depth.units = REA_data['depth'].units
            new['time'][:] = REA_data['time'][:]
            new['depth'][:] = REA_data['depth'][0:18]
            print('Done copying time and depth variable values.')
            
            temp = REA_data[varname][:,0:18,:,:] # Read continguously first.
            nc_var[:] = temp[:,:,9:250:12,0:673:12] # Now random-access in RAM to copy over.
            nc_var.units = REA_data[varname].units
            print('Done copying over {} values of one month.'.format(varname))
            
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


# In[27]:

for i in range(1,13):
    downsampling_monthly((2013,i,'votemper','TEMP'))
    downsampling_monthly((2013,i,'vosaline','PSAL'))
    downsampling_monthly((2014,i,'votemper','TEMP'))
    downsampling_monthly((2014,i,'vosaline','PSAL'))


# In[22]:

# from ipyparallel import Client
# rc = Client()
# dv = rc[:]
# lv = rc.load_balanced_view()
# inputs_1 = [(2013,month,'votemper','TEMP') for month in range(1,13)] + \
# [(2014,month,'votemper','TEMP') for month in range(1,13)]
# outputs_1 = lv.map(downsampling_monthly,inputs_1)
# inputs_2 = [(2013,month,'vosaline','PSAL') for month in range(1,13)] + \
# [(2014,month,'vosaline','PSAL') for month in range(1,13)]
# outputs_2 = lv.map(downsampling_monthly,inputs_2)
# msg_ids = list()
# msg_ids.extend(outputs_1.msg_ids+outputs_2.msg_ids)
# rc.get_result(msg_ids, owner=False).wait_interactive()

