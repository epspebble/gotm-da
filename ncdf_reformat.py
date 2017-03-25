### Common settings
import sys, os
home = os.getenv('HOME')
project_folder = os.path.join(home,'gotm-dst')

### Global settings
ERA_folder = os.path.join(project_folder,'p_sossta','medsea_ERA-INTERIM','3-hourly')
output_folder = os.path.join(project_folder,'medsea_GOTM')
GOTM_lat = tuple(30.75+0.75*i for i in range(21))
GOTM_lon = tuple(-6.0+0.75*i for i in range(57))

ERA_lat_ind = slice(-8,4,-1)
ERA_lon_ind = slice(12,-18,1)

with Dataset("{}/MEDSEA_ERA-INT_lwrd_y2013.nc".format(ERA_folder),'r') as ERA_data:
    print("Confirming whether the manually found lat,lon indices for ERA datasets really gives our GOTM grid...")
    print("lat: ", all(ERA_data['lat'][ERA_lat_ind] == GOTM_lat))
    print("lon: ", all(ERA_data['lon'][ERA_lon_ind] == GOTM_lon))
    
print("Latitudes: ", GOTM_lat)
print("Longitude: ",GOTM_lon)
