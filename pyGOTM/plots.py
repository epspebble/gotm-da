## Medsea results visualization toolbox

## The following were the old plotters since v2 (before 2017-07-13)

# Import matplotlib colormap for assigning default values to functions.
from matplotlib import cm
from numpy.ma import mean

def medsea_heatmap(data, # The only necessary argument: a netCDF Dataset to be opened by user.
                   varnames=['sst'],fun = lambda varargin: mean(varargin[0][:],axis=0),
                   ax=None, draw_colorbar=True, vlim = None,cmap=cm.coolwarm):
    # Expect a netCDF4 Data
    from netCDF4 import Dataset
    from matplotlib.pyplot import subplots
    if not(isinstance(data,Dataset)):
        raise Exception('First argument must be a netCDF4 Dataset.')

    # Load essential data
    lat = data['lat'][:]
    lon = data['lon'][:]
    #epoch = data['time'].units
    #time = data['time'][:]
    #z = data['z'][:]

    # Evaluate data to be plotted.
    varargin = [data[each] for each in varnames]
    val = fun(varargin)

    # Create figures / axis if not provided.
    if ax is None:
        fig, ax = subplots(figsize=(11,8.5)) # Full letter landscape size.
    fig = ax.get_figure()

    # Generate heatmap
    if vlim is None:
        im = ax.imshow(val, extent=(lon.min(), lon.max(), lat.min(), lat.max()),
                        interpolation='nearest', origin='lower', cmap=cmap)
    else:
        im = ax.imshow(val, extent=(lon.min(), lon.max(), lat.min(), lat.max()),
                        vmin = vlim[0], vmax = vlim[1],
                        interpolation='nearest', origin='lower', cmap=cmap)

    if draw_colorbar:
        cbar = fig.colorbar(im)

    ax.grid('on')
    ax.axis('tight')

    return fig, ax, im # Return them for further annotations.

def medsea_plot_mean(varname='sst', ax=None, cmap = cm.coolwarm,
                     fig_fn_prefix = 'medsea_ASM0_avg_sst', show_fig = True,
                     year=2014, month=1,
                     nc_folder='/home/simon/gotm-dst/medsea_GOTM/results',nc_fn_prefix='no_assim/p_sossta/medsea_GOTM'):

    from numpy.ma import mean
    data = medsea_data(year=year,month=month,nc_folder=nc_folder,nc_fn_prefix=nc_fn_prefix)
    fig, ax = medsea_heatmap(data, ax=ax, cmap=cmap,
                             varnames=[varname],fun=lambda varargin: mean(varargin[0][:],axis=0))
    # Annotations
    ax.set_title('mean {} for {:d}-{:02d}'.format(varname,year,month))

    # Save the figure
    fig_fn = '{}_{}_mean_{:d}{:02d}.png'.format(fig_fn_prefix,varname,year,month)
    print('Saving {}...'.format(fig_fn))
    fig.savefig(fig_fn)

    if not(show_fig):
        close(fig)

    # Remember to close the file.
    data.close()

def sst_range_ASM0(time,lon,sst):
    nday = int(time[-1]/86400)
    sst_max = ones((nday,M,N))*NaN
    sst_min = ones((nday,M,N))*NaN
    sst_range = ones((nday,M,N))*NaN
    for day in range(nday):
        for k in range(M*N):
            m = int(k/N)
            n = mod(k,N)
            offset = int(ceil(lon[n]/15))

            h1 = 11 - offset + day*24
            h2 = 21 - offset + day*24
            sst_max[day,m,n] = max(sst[h1:h2,m,n])


            h3 = 3 - offset + day*24 # We increased the starting hour from 1 (HCMR ppt) to 3 to avoid going to the previous day.
            h4 = 12 - offset + day*24
            #sst_min[day,m,n] = min(sst[h1:h2,m,n]) # OH NO! TYPO! The min was taken from the same range as in taking max!!! That means we're measuring the amount of cooling, but it's not all the way to the coolest point yet...
            sst_min[day,m,n] = min(sst[h3:h4,m,n])
    sst_range = sst_max - sst_min
    return sst_range

def sst_range_ASM2(time,lon,sst,swr):
    import numpy
    nday = int(time[-1]/86400)
    sst_max = ones((nday,M,N))*NaN
    sst_min = ones((nday,M,N))*NaN
    sst_range = ones((nday,M,N))*NaN
    for day in range(nday):
        for k in range(M*N):
            m = int(k/N)
            n = mod(k,N)
            offset = int(ceil(lon[n]/15))

            h1 = 11 - offset + day*24
            h2 = 21 - offset + day*24
            sst_max[day,m,n] = max(sst[h1:h2,m,n])

            swr_hourly_today = swr[day*24:(day+1)*24,m,n] # This could be a masked array! So find could return empty array.
            if size(find(swr_hourly_today)) == 0:
                swr_first_nonzero = -1000 # So that max(3-offset,swr_first_nonzero)=3-offset
            else:
                swr_first_nonzero = find(swr_hourly_today)[0]
            # Without taking the max of swr_first_nonzero, it could cover a jagged portion of sunrise assim.
            h3 = max([3-offset,swr_first_nonzero]) + day*24
            h4 = 12 - offset + day*24
            sst_min[day,m,n] = min(sst[h3:h4,m,n])
    sst_range = sst_max - sst_min
    return sst_range

def medsea_monthly_mean_heatmaps(assim = 0, show_fig = False, save_fig = True,
                                 months = [(2013,each) for each in range(1,13)] + [(2014,each) for each in range(1,13)],
                                 fig_folder = 'fig/monthly_means'):
    """ assim = 0 (none), 1 (local midnight), 2 (swr onset after midnight, i.e. sunrise) """

    # Choose the right subfolder of nc files by the assim argument.
    if assim == 0:
        nc_folder = 'results/p_sossta/no_assim'
    elif assim == 1:
        nc_folder = 'results/p_sossta/assim_midnight'
    elif assim == 2:
        nc_folder = 'results/p_sossta/assim_sunrise'
    else:
        raise Exception('assim = 0/1/2')

    for year,month in months:
        data = medsea_data(year=year,month=month,nc_folder=nc_folder)

        # A template plotter for monthly means.
        def do_heatmap(fig_title_prefix, fig_filename_prefix, **kwargs):
            fig, ax = medsea_heatmap(data,**kwargs)
            ax.set_title(fig_title_prefix + ' for {:d}-{:02d}'.format(year,month))
            if save_fig:
                fig_fn = os.path.join(base_folder,fig_folder,fig_filename_prefix+'_ASM{:d}_{:d}{:02d}.png'.format(assim,year,month))
                print('Saving ' + format(fig_fn) + '...')
                fig.savefig(fig_fn)
            if not(show_fig):
                close(fig)

        # Mean values of "sea surface variables", i.e. those with dimensions time, lat, lon but not depth.
        for varname in data.variables:
            if data[varname].dimensions == ('time','lat','lon'):
                do_heatmap('mean ' + varname,
                           'medsea_mean_' + varname,
                           varnames=[varname],
                           fun=lambda varargin: mean(varargin[0][:],axis=0),
                           cmap=cm.coolwarm)
        # Mean value of cooling effect
        do_heatmap('mean skin cooling effect',
                   'medsea_mean_cooling',
                   varnames=['sst','skint'],
                   fun=lambda args: mean(args[0][:]-args[1][:],axis=0),
                   cmap=cm.cool)

        # Mean diurnal warming amount
        if assim == 0:
            varnames = ['time','lon','sst']
            fun = lambda args: mean(sst_range_ASM0(args[0][:],args[1][:],args[2][:]),axis=0)
        elif assim == 2:
            varnames = ['time','lon','sst','swr']
            fun = lambda args: mean(sst_range_ASM2(args[0][:],args[1][:],args[2][:],args[3][:]),axis=0)
        else:
            raise Exception('Assim = 1 not implemented.')

        do_heatmap('mean diural warming amount in SST',
                   'medsea_mean_SST_range',
                   varnames=varnames, fun=fun,
                   cmap=cm.hot)
        data.close()

def buoy_comparisons(ax=None,year=2014,month=4,station='61277',plot_GOTM_SST=False):
    import os
    home = os.getcwd()
    os.chdir(base_folder)

    if ax is None:
        fig, ax = subplots()

    # Buoy
    if station == '61280' or station == '61281' or station == '61430':
        buoy_fn = 'buoys/{2}/IR_{0:d}{1:02d}_TS_MO_{2}.nc'.format(year,month,station)
    elif station == '61277' or station == '68422' or station == 'SARON':
        buoy_fn = 'buoys/{2}/MO_{0:d}{1:02d}_TS_MO_{2}.nc'.format(year,month,station)
    if station == '61277' or station == '68422':
        ind = 1
    else:
        ind = 0

    with MFDataset('buoys/{0:}/*{0:}.nc'.format(station),'r') as buoy_data:
        buoy_lat = (float(buoy_data.geospatial_lat_max) + float(buoy_data.geospatial_lat_min))/2
        buoy_lon = (float(buoy_data.geospatial_lon_max) + float(buoy_data.geospatial_lon_min))/2
        depth_buoy = buoy_data['DEPH'][0,ind]

    data = list()
    if os.path.isfile(buoy_fn):
        with Dataset(buoy_fn,'r') as buoy_data:
            print('Buoy {0} at {1} found for {2:d}-{3:02d}'.format(station,print_lat_lon(buoy_lat,buoy_lon),year,month))
            temp_buoy = buoy_data['TEMP'][:,ind]
            temp_buoy[temp_buoy==5.0] = NaN # Spurious 5.0 degree readings.
            time_buoy = buoy_data['TIME']
            date_buoy = num2date(time_buoy[:],time_buoy.units)
            ax.plot(date_buoy,temp_buoy,'red',
                    label='{} at {} ({:d}m)'.format(station,print_lat_lon(buoy_lat,buoy_lon),int(depth_buoy)))
            data.append((date_buoy,temp_buoy))

    # GOTM - ASM0
    with Dataset(os.path.join('results/p_sossta/no_assim','medsea_GOTM_{:d}{:02d}.nc'.format(year,month)),'r') as medsea_data:
        medsea_lat = medsea_data['lat'][:]
        medsea_lon = medsea_data['lon'][:]
        m = argmin(abs(medsea_lat-buoy_lat),axis=0)
        n = argmin(abs(medsea_lon-buoy_lon),axis=0)
        print('Closest medsea_GOTM grid point at (m,n) = ({},{})'.format(m,n), print_lat_lon(medsea_lat[m],medsea_lon[n]))
        #medsea_sst = medsea_data['sst']
        time_GOTM = medsea_data['time']
        date_GOTM = num2date(time_GOTM[:],time_GOTM.units)

        # Plot GOTM SST
        if plot_GOTM_SST:
            sst_GOTM = medsea_data['temp'][:,0,m,n]
            ax.plot(date_GOTM, sst_GOTM,'gray',linestyle='dashed',label='GOTM ASM0 SST at ' + \
                    print_lat_lon(medsea_lat[m],medsea_lon[n]))
            data.append((date_GOTM,sst_GOTM))

        # 2016-12-19, temporary solution for plotting GOTM TEMP at buoy depth.
        from scipy.interpolate import RectBivariateSpline as interp
        depth_GOTM = medsea_data['depth'][:]
        temp_GOTM = medsea_data['temp'][:,:,m,n]
        f_GOTM = interp(time_GOTM,depth_GOTM,temp_GOTM,kx=3,ky=3, # bicubic interpolation
                        bbox=[0,time_GOTM[-1],0,max(depth_GOTM)])
        temp_GOTM_at_buoy_depth = f_GOTM(time_GOTM,depth_buoy)
        ax.plot(date_GOTM, temp_GOTM_at_buoy_depth,'gray', linestyle='solid',
                label='GOTM ASM0 TEMP ({0:d}m) '.format(int(depth_buoy)) + \
                print_lat_lon(medsea_lat[m],medsea_lon[n]))
        data.append((date_GOTM,temp_GOTM_at_buoy_depth))

    # GOTM - ASM2
    with Dataset(os.path.join('results/p_sossta/assim_sunrise','medsea_GOTM_{:d}{:02d}.nc'.format(year,month)),'r') as medsea_data:
        medsea_lat = medsea_data['lat'][:]
        medsea_lon = medsea_data['lon'][:]
        m = argmin(abs(medsea_lat-buoy_lat),axis=0)
        n = argmin(abs(medsea_lon-buoy_lon),axis=0)
        print('Closest medsea_GOTM grid point at (m,n) = ({},{})'.format(m,n), print_lat_lon(medsea_lat[m],medsea_lon[n]))
        #medsea_sst = medsea_data['sst']
        time_GOTM = medsea_data['time']
        date_GOTM = num2date(time_GOTM[:],time_GOTM.units)

        # Plot GOTM SST
        if plot_GOTM_SST:
            sst_GOTM = medsea_data['temp'][:,0,m,n]
            ax.plot(date_GOTM, sst_GOTM,'blue',linestyle='dashed',label='GOTM ASM2 SST at ' + \
                    print_lat_lon(medsea_lat[m],medsea_lon[n]))
            data.append((date_GOTM,sst_GOTM))

        # 2016-12-19,  temporary solution for plotting GOTM TEMP at buoy depth.
        from scipy.interpolate import RectBivariateSpline as interp
        depth_GOTM = medsea_data['depth'][:]
        temp_GOTM = medsea_data['temp'][:,:,m,n]
        f_GOTM = interp(time_GOTM,depth_GOTM,temp_GOTM,kx=3,ky=3, # bicubic interpolation
                        bbox=[0,time_GOTM[-1],0,max(depth_GOTM)])
        temp_GOTM_at_buoy_depth = f_GOTM(time_GOTM,depth_buoy)
        ax.plot(date_GOTM, temp_GOTM_at_buoy_depth,'blue', linestyle='solid',
                label='GOTM ASM2 TEMP ({0:d}m) '.format(int(depth_buoy)) + \
                print_lat_lon(medsea_lat[m],medsea_lon[n]))
        data.append((date_GOTM,temp_GOTM_at_buoy_depth))

    # REA
    with Dataset(os.path.join('profiles','medsea_rea_votemper_{:d}{:02d}.nc'.format(year,month)),'r') as rea_data:
        votemper = rea_data['votemper']
        time_rea = rea_data['time']
        date_rea = num2date(time_rea[:],time_rea.units)
        temp_REA = votemper[:,0,m,n]
        if votemper[0,0,m,n] > 1e10: # the first (1st index) temperature record at shallowest depth (2nd index)
            print('Possible land location!')
        elif votemper[0,-1,m,n] > 1e10: # the first (1st index) temperature record at deepest depth (2nd index)
            print('Possible shallow water location!')

        depth_REA = rea_data['depth'][0]
        ax.plot(date_rea,temp_REA,'green', label='REA ({:.2f}m) at '.format(depth_REA) + \
                print_lat_lon(rea_data['lat'][m],rea_data['lon'][n]))
        data.append((date_rea,temp_REA))

    os.chdir(home)
    return ax.get_figure(), ax, data

def buoy_comparisons_full_year(year=2014,station='61277',showfig=True,fig_folder='fig/buoys'):
    fig, axes = subplots(ncols=1,nrows=12,figsize=(17,22))
    for i,ax in enumerate(axes):
        buoy_comparisons(ax=ax,year=year,month=i+1,station=station)

    # 1.4721m is the shallowest depth in REA data
    suptitle('GOTM-ASM0 SST (gray), GOTM-ASM2 SST (blue), Buoy data (red) vs REA at 1.4721m (green) for near {}'.format(station),fontsize=16,y=0.92)
    fig_fn = os.path.join(base_folder,fig_folder,'SST_GOTM_vs_REA_near_{}_{:d}.png'.format(station,year))
    print('Saving ' + fig_fn + '...')
    fig.savefig(fig_fn)
    if not(showfig):
        close(fig)

def SST_monthly_comparisons(year=2013,month=2,m=7,n=42,runs=['ASM0','ASM2'],colors=None,
                            ax=None,pretty=True,output_data=True,fig_fn=None):
    from netCDF4 import Dataset, num2date
    import matplotlib.cm as cm
    from matplotlib.pyplot import subplots
    from numpy import linspace

    if ax is None:
        fig, ax = subplots(figsize=(16,4))

    if colors is None:
        cfun = cm.rainbow
        ax.set_prop_cycle('color',[cfun(x) for x in linspace(0,1,len(runs))]) # Extra "run" for rea data
        colors = [cfun(i/len(runs)) for i in range(len(runs))]

    data = list()
    # GOTM runs
    for i, run in enumerate(runs):
        with Dataset(os.path.join(base_folder,run,'medsea_GOTM_{:d}{:02d}.nc'.format(year,month)),'r') as medsea_data:
            time_GOTM = medsea_data['time']
            date_GOTM = num2date(time_GOTM[:],time_GOTM.units)
            sst_GOTM = medsea_data['sst'][:,m,n]
            ax.plot(date_GOTM, sst_GOTM, color=colors[i], label=run)
            if output_data:
                data.append({run: (date_GOTM,sst_GOTM)})

    # REA
    with Dataset(os.path.join(data_folder,'medsea_rea','medsea_rea_votemper_{:d}{:02d}.nc'.format(year,month)),'r') as rea_data:
        votemper = rea_data['votemper']
        time_rea = rea_data['time']
        date_rea = num2date(time_rea[:],time_rea.units)
        temp_REA = votemper[:,0,m,n]

        if votemper[0,0,m,n] > 1e10: # the first (1st index) temperature record at shallowest depth (2nd index)
            print('Possible land location!')
        elif votemper[0,-1,m,n] > 1e10: # the first (1st index) temperature record at deepest depth (2nd index)
            print('Possible shallow water location!')

        depth_REA = rea_data['depth'][0]
        ax.plot(date_rea,temp_REA,color='black',linewidth=2,label='REA'.format(depth_REA))
        ax.set_xlim(left=date_rea[0],right=date_rea[-1])
        ax.set_xticks(date_rea,minor=False)
        ax.set_xticklabels([date.day for date in date_rea])
        if output_data:
            data.append({'rea': (date_rea,temp_REA)})

    fig = ax.get_figure()
    if pretty:
        # fig.autofmt_xdate()
        ax.set_title('SST comparisons at {} for {:d}-{:02d}'.format(print_lat_lon(*get_lat_lon(m,n)),year,month))
        ax.grid('on')
        ax.legend()
    if output_data:
        return fig, ax, data
    else:
        return fig, ax

def SST_yearly_comparisons(year,m,n,runs=['ASM0','ASM2'],showfig=True,
                           fig_subfolder = 'fig/spotchk'):
    from matplotlib.pyplot import subplots, suptitle, legend, close

    fig, axes = subplots(ncols=1,nrows=12,figsize=(17,22),sharex=False)
    for i,ax in enumerate(axes):
        SST_monthly_comparisons(year=year,month=i+1,m=m,n=n,runs=['ASM0','ASM2'],ax=ax,pretty=False,output_data=False)
        ax.grid('on')
        ax.set_xticklabels([])
        if i == 0:
            ax.legend()
            lines = ax.get_lines()

    latlon = print_lat_lon(*get_lat_lon(m,n))
    # 1.4721m is the shallowest depth in REA data
    suptitle('SST comparisons for the year {:d} at {:s}'.format(year, latlon),
             fontsize=16,y=0.92)
    labels = (*runs, 'REA at 1.4721m')
    # legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0,-0.1,1,1),
    #        bbox_transform = fig.transFigure )
    fig_fn = os.path.join(base_folder,fig_subfolder,'SST_{}_GOTM_ASM0_vs_ASM2_vs_REA_{:d}.png'.format(latlon, year))
    print('Saving ' + fig_fn + '...')
    fig.savefig(fig_fn)
    if not(showfig):
        close(fig)
