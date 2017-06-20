from pylab import *
from netCDF4 import Dataset
from numpy.ma import masked_invalid, masked_outside
import os, sys
from scipy.interpolate import *

MODIS_folder = '/global/scratch/simontse/p_sossta/glo_MODIS/8days'
data_folder = '/global/scratch/simontse'

def fn(year,num):
    "num is the number of 8day periods elapsed, i.e. from 0 to 45 inclusive."
    assert isinstance(num,int) and num <= 45 and num >= 0
    start = 1+num*8
    yrdays = 366 if mod(year,4) == 0 else 365 
    end = yrdays if num == 45 else start+7
#    print('year,start,end',year,start,end)
    fn = os.path.join(MODIS_folder,'A{0:d}{1:03d}{0:d}{2:03d}.L3m_8D_CHL_chlor_a_9km.nc'.format(year,start,end))
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
        modis_lat = nc['lat'][711:529:-1]
        modis_lon = nc['lon'][2087:2593]
        modis_chlo = masked_outside(nc['chlor_a'][711:529:-1,2087:2593],0,3)
    return modis_lat, modis_lon, modis_chlo

def plot(year,num, ax=None):
    modis_lat, modis_lon, modis_chlo = data(year,num)
    if ax is None:
        fig, ax = subplots(figsize=(10,3))
    fig = ax.get_figure()
    im = ax.imshow(modis_chlo, extent=(modis_lon.min(), modis_lon.max(), modis_lat.min(), modis_lat.max()),
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
    im, cbar = plot(year,1, ax=ax)
    # ax.set_title('{:d}, 8-Day beginning on Day #{:003d}'.format(year,1))

    #Set up plot backgroudn
    def init():
        im.set_alpha(0)
        return im,
    def animate(i):
        im.set_alpha(1)
        im.set_array(data(year,i)[2])
        ax.set_title('{:d}, 8days mean beginning on day #{:003d}'.format(year,(i-1)*8+1))
        return im,

    anim = animation.FuncAnimation(fig, animate, #init_func=init,
                                   frames=range(1,47), interval=28, blit=True)

    anim.save('{:d}_medsea_MODIS.gif'.format(year), writer='imagemagick', fps=3)
    print('Animation saved to {:d}_medsea_MODIS.gif!'.format(year))
    return anim
# HTML(anim.to_html5_video())

def data_stat(years=[2013,2014,2015]):
    cumsum = 0
    minval = [1,2000,0]
    maxval = [0,2000,0]
    if isinstance(years,int):
        years = [years]
    for year in years:
        for i in range(46):
            [lat,lon,chlo] = data(year+int(i/45),int(1+mod(i,45)))
            val = sum(~chlo.mask)/(len(lat)*len(lon))
            cumsum += val
            if val < minval[0]:
                minval = [val,year,i]
            if val > maxval[0]:
                maxval = [val,year,i]
    avgval = cumsum/(len(years)*46)
    print('Minimum coverage: {!s} at {:d} 8day-period #{:d}'.format(*minval))
    print('Maximum coverage: {!s} at {:d} 8day-period #{:d}'.format(*maxval))
    print('Average coverage: {!s}'.format(avgval))
    return minval,maxval,avgval

def mapcomp(z1,z2):
    fig, axes = subplots(1,2,figsize=(14,6))
    ax1, ax2 = axes
    im1 = ax1.imshow(z1, extent=(-6,36.25,30.75,45.75), origin='lower', vmin=0,vmax=3, cmap=cm.coolwarm)
    # cbar1 = fig.colorbar(im1)
    im2 = ax2.imshow(z2, extent=(-6,36.25,30.75,45.75), origin='lower', vmin=0,vmax=3, cmap=cm.coolwarm)
    # cbar2 = fig.colorbar(im2)
    #     ax.grid('on')
    for ax in axes:
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
    fig.tight_layout()

def get_REA_grid():
    with Dataset('/global/scratch/simontse/p_sossta/medsea_rea/2013/20130101_TEMP_re-fv6.nc','r') as nc:
        lat_rea = nc['lat'][:]
        lon_rea = nc['lon'][:]
        temp_rea = nc['votemper'][:]
        is_sea = ~temp_rea[0,0,:].mask
    return lat_rea, lon_rea, is_sea

def interp(year,num,method='linear'):
    # Data
    y,x,z = data(year,num)
    xx, yy = np.meshgrid(x,y)
    loc = array([d for d in zip(xx.reshape(-1),yy.reshape(-1))])
    val = z.reshape(-1)
    if method == 'linear':
        interpolant = LinearNDInterpolator(loc[~val.mask],val[~val.mask])
    else: # defaults to nearest:
        interpolant = NearestNDInterpolator(loc[~val.mask],val[~val.mask])
    
    # New grid
    lat_rea, lon_rea, is_sea = get_REA_grid()
    xx_rea, yy_rea = meshgrid(lon_rea,lat_rea)
    zz_rea = ones(xx_rea.shape)
    for m in range(len(lat_rea)):
        for n in range(len(lon_rea)):
            if is_sea[m,n]:
                zz_rea[m,n] = interpolant(xx_rea[m,n],yy_rea[m,n])
            else:
                zz_rea[m,n] = nan     
    return zz_rea

def plot_interp(year,num, ax=None):
    from netCDF4 import Dataset
    with Dataset(os.path.join(data_folder,'medsea_MODIS','medsea_MODIS_chlor_a_8D_{}.nc'.format(year))) as nc:
        chlo = nc['chlor_a'][num,:]
    if ax is None:
        fig, ax = subplots(figsize=(10,3))
    fig = ax.get_figure()
    im = ax.imshow(chlo, extent=(-6., 36.25, 30.75, 45.75),vmin=0,vmax=3,
                   origin='lower',cmap=cm.coolwarm)
    cbar = fig.colorbar(im)
#     ax.grid('on')
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    fig.tight_layout()
    return im, cbar, chlo

def animate_interp(year):
    from matplotlib import animation
    from IPython.display import HTML

    # Set up figure, axis and plot elements.
    fig, ax = subplots(figsize=(10,3))
    im, cbar, chlo = plot_interp(year,0,ax=ax)
    # ax.set_title('{:d}, 8-Day beginning on Day #{:003d}'.format(year,1))

    #Set up plot backgroudn
    def init():
        im.set_alpha(0)
        return im,
    def animate(i):
        im.set_alpha(1)
        with Dataset(os.path.join(data_folder,'medsea_MODIS','medsea_MODIS_chlor_a_8D_{}.nc'.format(year))) as nc:
            chlo = nc['chlor_a'][i,:]
        im.set_array(chlo)
        ax.set_title('{:d}, 8days mean beginning on day #{:003d}'.format(year,i*8+1))
        return im,

    anim = animation.FuncAnimation(fig, animate, #init_func=init,
                                   frames=range(46), interval=28, blit=True)

    anim.save('{:d}_medsea_MODIS_rea.gif'.format(year), writer='imagemagick', fps=3)
    print('Animation saved to {:d}_medsea_MODIS_rea.gif!'.format(year))
    return anim
# HTML(anim.to_html5_video())
