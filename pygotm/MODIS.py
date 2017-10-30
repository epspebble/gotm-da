from pylab import *
from netCDF4 import Dataset
from numpy.ma import masked_invalid, masked_outside
import os, sys
from scipy.interpolate import *
from pygotm import medsea

target_grid = '144x'
medsea.set_grid(target_grid)
medsea.set_folders()

data_folder = medsea.data_folder
p_sossta_folder = medsea.p_sossta_folder
grid_lats = medsea.grid_lats
grid_lons = medsea.grid_lons

coverage_length = '8D'
MODIS_folder = os.path.join(p_sossta_folder,'glo_MODIS',coverage_length)
# usually treating one variable at a time, declare as module global variable to avoid passing the same arguments
# over and over again
obs_type = 'IOP'
varname = 'a_488_giop'

def fn(year, num):
    """ 
    The naming convention of MODIS data files seem to be the concatenation of the following pieces:
    1. [A] -> MODIS-AQUA
    2. [YYYYDDD] -> coverage start year & day of year (for 8day files)
    3. [YYYYDDD] -> coverage end year & day of year 
    4. .L3M (processing level)
    5. _8D -> coverage time length
    6. _IOP -> type of data (3 capital letters) (e.g. IOP, CHL)
    7. _a_488_giop -> variable name inside the netCDF dataset (so that ds[varname] works if ds is a netCDF4.Dataset)
    8. _9km -> approx. spatial resolution (1km / 4km / 9km)
    9. .nc -> file extension.

    """
    "num is the number of 8day periods elapsed, i.e. from 0 to 45 inclusive."

    assert isinstance(num,int) and num <= 45 and num >= 0
    start = 1+num*8
    yrdays = 366 if mod(year,4) == 0 else 365 
    end = yrdays if num == 45 else start+7
#    print('year,start,end',year,start,end)
    fn = os.path.join(MODIS_folder,'A{0:d}{1:03d}{0:d}{2:03d}.L3m_8D_{3:s}_{4:s}_9km.nc'.format(year,start,end,obs_type,varname))
    if not(os.path.isfile(fn)):
        print('year,start,end',year,start,end)
        raise OSError(fn + 'not found.')
    return fn

def get_ncattr(year=2013,num=13,ncattrs=['time_coverage_start','time_coverage_end']):
    with Dataset(fn(year,num),'r') as nc:
        print([nc.getncattr(name) for name in ncattrs])

def coverage_length(year=2013,num=4):
    from datetime import datetime
    def conv(date_string):
        return datetime.strptime(date_string,'%Y-%m-%dT%H:%M:%S.000Z')
    with Dataset(fn(year,num),'r') as nc:
        start = conv(nc.getncattr('time_coverage_start'))
        end = conv(nc.getncattr('time_coverage_end'))
        delta = end - start
    return delta, start, end
    
def coverage_midpoint(year=2014,num=4):
    delta, start, end = coverage_length(year,num)
    return start + delta/2        
        
def data(year,num):
    with Dataset(fn(year,num),'r') as nc:
        modis_lats = nc['lat'][711:529:-1]
        modis_lons = nc['lon'][2087:2593]
        modis_data = masked_outside(nc[varname][711:529:-1,2087:2593],0,3)
    return modis_lats, modis_lons, modis_data

def plot(year,num,ax=None):
    modis_lats, modis_lons, modis_data = data(year,num)
    if ax is None:
        fig, ax = subplots(figsize=(10,3))
    fig = ax.get_figure()
    im = ax.imshow(modis_data, extent=(modis_lons.min(), modis_lons.max(), modis_lats.min(), modis_lats.max()),
                   interpolation='nearest', origin='lower', cmap=cm.coolwarm)
    cbar = fig.colorbar(im)
#     ax.grid('on')
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    fig.tight_layout()
    return im, cbar

def animate(year):
    from matplotlib import animation
    from IPython.display import HTML

    # Set up figure, axis and plot elements.
    fig, ax = subplots(figsize=(10,3))
    im, cbar = plot(year,1,ax=ax)
    # ax.set_title('{:d}, 8-Day beginning on Day #{:003d}'.format(year,1))

    #Set up plot backgroudn
    def init():
        im.set_alpha(0)
        return im,
    def animate(i):
        im.set_alpha(1)
        im.set_array(data(year,i)[2])
        ax.set_title('{:d}, 8days {:s} mean beginning on day #{:003d}'.format(year,varname,(i-1)*8+1))
        return im,

    anim = animation.FuncAnimation(fig, animate, #init_func=init,
                                   frames=range(1,47), interval=28, blit=True)
    
    gif_fn = '{:d}_medsea_MODIS_{:s}.gif'.format(year,varname)
    anim.save(gif_fn, writer='imagemagick', fps=3)
    print('Animation saved to {:s}!'.format(gif_fn))
    return anim
# HTML(anim.to_html5_video())

def data_stat(years=[2013,2014,2015,2016]):
    cumavail = 0
    minavail = [1,2000,0]
    maxavail = [0,2000,0]

    minval = 1.e20
    maxval = -1.e20
    cumsum = 0
    
    if isinstance(years,int):
        years = [years]

    for year in years:
        for i in range(46):
            [lat,lon,modis_data] = data(year+int(i/45),i%45)

            # Find the min, max and mean value.
            new_min = modis_data.min()
            if new_min < minval:
                minval = new_min
            new_max = modis_data.max()
            if new_max > maxval:
                maxval = new_max
            cumsum += modis_data.mean()

            # Finding the min and max coverage period.
            avail = sum(~modis_data.mask)/(len(lat)*len(lon))
            cumavail += avail
            if avail < minavail[0]:
                minavail = [avail,year,i]
            if avail > maxavail[0]:
                maxavail = [avail,year,i]
                
    avgavail = cumavail/(len(years)*46)
    avgval = cumsum/(len(years)*46)

    print('Statistics for MODIS {:s} over {!s}...'.format(varname,years))
    print('Value:')
    print('  Min: {!s}'.format(minval))
    print('  Max: {!s}'.format(maxval))
    print('  Mean: {!s}'.format(avgval))
    print('Availability:')
    print('  Worst period: {:.4f} at {:d} 8day-period #{:d}'.format(*minavail))
    print('  Best period: {:.4f} at {:d} 8day-period #{:d}'.format(*maxavail))
    print('  Average: {:.2f}%'.format(avgavail*100))
    return {'stat': (minval,maxval,avgval), 'avail': (minavail,maxavail,avgavail)}

# Reference ranges of data for plotting over the medsea.
def set_vmin_vmax():
    global vmin, vmax
    if varname == 'chlor_a':
        # Max and min of Ohlmann et al.'s formulas for case 12.
        vmin = 0
        vmax = 3
    if varname == 'a_488_giop':
        # 2 sig fig from MODIS medsea data
        vmin = 0.017
        vmax = 0.10
    if varname == 'bb_488_giop':
        # 2 sig fig from MODIS medsea data
        vmin = 0.0014
        vmax = 0.017

def mapcomp(z1,z2):
    set_vmin_vmax()
    
    fig, axes = subplots(1,2,figsize=(14,6))
    ax1, ax2 = axes
    
    im1 = ax1.imshow(z1, extent=(-6,36.25,30.75,45.75), origin='lower',
                     vmin=vmin,vmax=vmax, cmap=cm.coolwarm)
    # cbar1 = fig.colorbar(im1)
    im2 = ax2.imshow(z2, extent=(-6,36.25,30.75,45.75), origin='lower',
                     vmin=vmin,vmax=vmax, cmap=cm.coolwarm)
    # cbar2 = fig.colorbar(im2)
    #     ax.grid('on')
    for ax in axes:
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
    fig.tight_layout()

def interp(year,num,method='linear'):
    from numpy.ma import masked_all

    set_vmin_vmax()
    # Data
    y,x,z = data(year,num)
    xx, yy = np.meshgrid(x,y)
    loc = array([d for d in zip(xx.reshape(-1),yy.reshape(-1))])
    val = z.reshape(-1)
    if method == 'linear':
        interpolant = LinearNDInterpolator(loc[~val.mask],val[~val.mask])
    else: # fallback to nearest:
        interpolant = NearestNDInterpolator(loc[~val.mask],val[~val.mask])
    
    # New grid
    xx_new, yy_new = meshgrid(grid_lons,grid_lats)
    zz_new = masked_all(xx_new.shape)
    for m,n in medsea.sea_mn:
        zz_new[m,n] = interpolant(xx_new[m,n],yy_new[m,n])
    return zz_new

def plot_interp(year,num, ax=None):

    set_vmin_vmax()
    from netCDF4 import Dataset
    with Dataset(os.path.join(data_folder,'medsea_MODIS','8days','medsea_MODIS_{:s}_8D_{:d}.nc'.format(obs_type,year))) as nc:
        modis_data = nc[varname][num,:]
    if ax is None:
        fig, ax = subplots(figsize=(10,3))
    fig = ax.get_figure()
    im = ax.imshow(modis_data, extent=(-6., 36.25, 30.75, 45.75),
                   vmin=vmin,vmax=vmax,
                   origin='lower',cmap=cm.coolwarm)
    cbar = fig.colorbar(im)
#     ax.grid('on')
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    fig.tight_layout()
    return im, cbar, modis_data

def animate_interp(year,num=46):
    from matplotlib import animation
    from IPython.display import HTML

    # Set up figure, axis and plot elements.
    fig, ax = subplots(figsize=(10,3))
    im, cbar, modis_data = plot_interp(year,0,ax=ax)
    # ax.set_title('{:d}, 8-Day beginning on Day #{:003d}'.format(year,1))

    #Set up plot backgroudn
    def init():
        im.set_alpha(0)
        return im,
    def animate(i):
        im.set_alpha(1)
        with Dataset(os.path.join(data_folder,'medsea_MODIS','8days','medsea_MODIS_{:s}_8D_{}.nc'.format(obs_type,year))) as nc:
            modis_data = nc[varname][i,:]
        im.set_array(modis_data)
        ax.set_title('{:d}, 8days mean beginning on day #{:003d}'.format(year,i*8+1))
        return im,

    anim = animation.FuncAnimation(fig, animate, #init_func=init,
                                   frames=range(num), interval=28, blit=True)

    fn = '{:d}_medsea_MODIS_{:s}.gif'.format(year,varname)
    anim.save(fn, writer='imagemagick', fps=3)
    print('Animation saved to {:s}!'.format(fn))
    return anim
# HTML(anim.to_html5_video())

# def format_modis_8days(year,outdir=os.path.join(data_folder,'medsea_MODIS','8days'),
#                        lat=medsea.grid_lats,lon=medsea.grid_lons):
#     from netCDF4 import Dataset, date2num
#     from numpy.ma import masked_invalid
#     global varname

#     epoch = 'seconds since 1981-01-01'

#     # Writing to an nc file.
#     outfn = os.path.join(outdir,'medsea_MODIS_{:s}_8D_{:d}.nc'.format(obs_type,year))
#     print('Writing to {:s}...'.format(outfn))
#     with Dataset(outfn, 'w', format='NETCDF3_CLASSIC') as new: 

#         new.createDimension('time') # unlimited
#         new.createDimension('lat', size = len(lat)) 
#         new.createDimension('lon', size = len(lon)) 
#         print('Done creating dimensions.')
        
#         nctime = new.createVariable('time','f8',dimensions=('time',))
#         nctime.units = epoch
#         nclat = new.createVariable('lat','f8',dimensions=('lat',))
#         nclat.units = 'degrees north'
#         nclon = new.createVariable('lon','f8',dimensions=('lon',))
#         nclon.units = 'degrees east'
#         nclat[:] = lat
#         nclon[:] = lon
#         print('Done creating lat/lon grid.')

#         for i in range(46):
#             nctime[i] = date2num(coverage_midpoint(year,i),epoch)
#         print("Done writing 'time' using coverage midpoints of the nc files for {:s}".format(varname))
            
#         if obs_type == 'IOP':
#             # ignore the current global varname, these two are what we want
#             for varname in ['a_488_giop','bb_488_giop']:
#                 ncvar = new.createVariable(varname,'f8',dimensions=('time','lat','lon'),
#                                            zlib=True, fill_value=nan) # Fill-value the same as in REA dataset.       
#                 for i in range(46):
#                     ncvar[i,:] = masked_invalid(interp(year,i,method='linear'))
#                 print("Done writing '{:s}'.".format(varname))

#         if obs_type == 'CHL':
#             assert varname == 'chlor_a' # only one variable in this obs_type
#             ncvar = new.createVariable('chlor_a','f8',dimensions=('time','lat','lon'),
#                                        zlib=True, fill_value=nan) # Fill-value the same as in REA dataset.
#             for i in range(46):
#                 ncvar[i,:] = masked_invalid(interp(year,i,method='linear'))
#             print("Done writing '{:s}'.".format(varname))

def format_modis_8days(year,num=46,outdir=os.path.join(data_folder,'medsea_MODIS','8days'),
                       lat=MODIS.grid_lats,lon=MODIS.grid_lons):
    from netCDF4 import Dataset, date2num
    from numpy.ma import masked_invalid
    
    epoch = 'seconds since 1981-01-01'
    
    # Writing to an nc file.
    outfn = os.path.join(outdir,'medsea_MODIS_{:s}_8D_{:d}.nc'.format(MODIS.obs_type,year))
    print('Writing to {:s}...'.format(outfn))
    mode = 'a' if os.path.isfile(outfn) else 'w'
    with Dataset(outfn, mode, format='NETCDF3_CLASSIC') as new:
        if mode == 'w':
            new.createDimension('time') # unlimited
            new.createDimension('lat', size = len(lat))
            new.createDimension('lon', size = len(lon))
            print('Done creating dimensions.')
            
            nctime = new.createVariable('time','f8',dimensions=('time',))
            nctime.units = epoch
            nclat = new.createVariable('lat','f8',dimensions=('lat',))
            nclat.units = 'degrees north'
            nclon = new.createVariable('lon','f8',dimensions=('lon',))
            nclon.units = 'degrees east'
            
            nclat[:] = lat
            nclon[:] = lon
            print('Done creating lat/lon grid.')

        # Common to mode == 'w' and mode == 'a'    
        ncvar = new.createVariable(MODIS.varname,'f8',dimensions=('time','lat','lon'),
                                   zlib=True, fill_value=nan) # Fill-value the same as in REA dataset.
        print('Done creating variables.')
        
        for i in range(num):
            if mode == 'w':
                nctime[i] = date2num(MODIS.coverage_midpoint(year,i),epoch)
                
            ncvar[i,:] = masked_invalid(MODIS.interp(year,i,method='linear'))
                
        print("Done writing 'time' and '{:s}'.".format(MODIS.varname))
