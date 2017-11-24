#!/home/simontse/miniconda3/bin/python3
import os, sys
from pygotm import medsea as ms
from netCDF4 import *
from numpy.ma import *
from numpy import loadtxt, nan, save

# Copied from pygotm/medsea.py on 2017-05-24
def combine_stat(year,format='NETCDF3_CLASSIC'):
    import os
    from numpy import loadtxt
    from numpy import ma

    dat_fn = dsfn.format(year)
    
    print('Combining GOTM daily statistics in {:s}, for {:d}...'.format(dat_fn,year))
    outdir = os.path.join(ms.data_folder,run+'_'+grid)
    if not(os.path.isdir(outdir)):
        os.mkdir(outdir)
    outfn = os.path.join(outdir,'daily_stat_{0:s}_{1:s}_{2:d}.nc'.format(run,grid,year))

    print('Writing dimensions and metadata...')
    elapsed = 0
    ms.tic()

    with Dataset(outfn,'w',format=format) as ds:

        # Dimensions
        ds.createDimension('day_of_year')
        ds.createDimension('lat', size = len(ms.grid_lats))
        ds.createDimension('lon', size = len(ms.grid_lons))
        
        # Dimension variables.
        day_var = ds.createVariable('day_of_year','i4',dimensions=('day_of_year',))
        day_var.units = "number of days since last new year's eve"
        lat_var = ds.createVariable('lat','f4',dimensions=('lat',))
        lat_var.units = 'degrees north'
        lon_var = ds.createVariable('lon','f4',dimensions=('lon',))
        lon_var.units = 'degrees east'
        lat_var[:] = ms.grid_lats
        lon_var[:] = ms.grid_lons
        print('Done initializing dimensions.')
        
        data = [None for i in range(7)]
        data[0]=ds.createVariable('assim_time','i4',dimensions=('day_of_year','lat','lon'))
        data[0].units = 'number of seconds since midnight in local time'
        data[4]=ds.createVariable('SST_max','f4',dimensions=('day_of_year','lat','lon'))
        data[4].units = 'degree Celsius'
        data[2]=ds.createVariable('SST_min_day','f4',dimensions=('day_of_year','lat','lon'))
        data[2].units = 'degree Celsius'
        data[6]=ds.createVariable('SST_min_night','f4',dimensions=('day_of_year','lat','lon'))
        data[6].units = 'degree Celsius'
        data[3]=ds.createVariable('SST_max_time','i4',dimensions=('day_of_year','lat','lon'))
        data[3].units = 'number of seconds since midnight in local time'
        data[1]=ds.createVariable('SST_min_day_time','i4',dimensions=('day_of_year','lat','lon'))
        data[1].units = 'number of seconds since midnight in local time'
        data[5]=ds.createVariable('SST_min_night_time','i4',dimensions=('day_of_year','lat','lon'))
        data[5].units = 'number of seconds since midnight in local time'
        print('Done creating variables.')

        yrdays = 366 if year%4==0 else 365
        day_var[:] = range(1,yrdays+1)
        ndfv = [0,0,99.,0,-99.,0,99.]
        temp = [nan*ones((365,ms.M,ms.N)) for i in range(7)]
        ms.tic()

        err_count = 0
        for i in range(ms.sea_m.size):
            m = ms.sea_m[i]
            n = ms.sea_n[i]
            latlong = ms.print_lat_lon(*ms.get_lat_lon(m,n))
            # This run spans two years.
            stat_fn = os.path.join(ms.get_local_folder(m,n),dat_fn)
            if not(os.path.isfile(stat_fn)):
                err_count += 1
                print('\nFile not found: ' + latlong)
                continue
            try:
                tmp = loadtxt(stat_fn)
            except ValueError as ve :
                print('\nValue error: {!s}'.format(ve) + latlong)
                err_count += 1
                continue
            if tmp.size == 0: # empty file
                print('\nEmpty file: ' + latlong)
                err_count +=1
                continue
            
            for k in range(7):
                if dat_fn[:25] == 'daily_stat_20130012014365':
                    if tmp.shape[0] != 729:
                        err_count += 1
                        print('\nUnexpected line count: {:d} for the location {!s}'.format(tmp.shape[0],latlong))
                        break
                    if year == 2013:
                        temp[k][:364,m,n] = tmp[:364,k+1]
                    elif year == 2014:
                        temp[k][:364,m,n] = tmp[365:,k+1]
                else:
                    temp[k][:tmp.shape[0],m,n] = tmp[:,k+1]
            # Beware, since the progress message overwrites itself, error messages need to do a newline first.
            print('Reading {0:d}/{1:d} grid points, {2:d}/{1:d} bypassed.'.format(i+1,ms.sea_m.size,err_count),end='\r')
        print('\nDone reading in all data.')
        elapsed += ms.toc()
        #save('daily_stats',temp)

        for i in range(7):
            ms.tic()
            data[i][:]=ma.masked_equal(temp[i],ndfv[i])
            #print("Done writing a variable.")
            #elapsed += ms.toc()
            
        print('Done writing out data to {:s}.'.format(outfn))
        elapsed += ms.toc()
        
if __name__ == '__main__':
    print('Number of command line arguments received: {:d}.'.format(len(sys.argv)))    
    if len(sys.argv) == 1:
        grid = '144x'        
        run = 'ASM3-75m'
        dsfn = 'daily_stat_20130012014365_{:s}.dat'.format(run)
        years = [2013,2014]
    elif len(sys.argv) == 3:
        grid = sys.argv[1]
        run = sys.argv[2]        
        dsfn = 'daily_stat_' + run + '_{:d}.dat'
        years = [2013,2014]
    elif len(sys.argv) == 4:
        grid, run, dsfn = sys.argv[1:]
        years = [2013,2014]
    else: # The trailing arguments are all years.
        grid, run, dsfn = sys.argv[1:4]
        years = [int(yr) for yr in sys.argv[4:]]
    ms.run = run
    ms.set_grid(grid)

    for year in years:
        combine_stat(year)
