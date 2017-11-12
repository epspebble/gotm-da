from .config import *

from math import sin, cos, asin, acos, tan, atan

### For SWRD Calcluations

# Reference: https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF

# Temporary global variable
year = 2014

## Helper functions        
def tz(lon):
    if lon > 0:
        return int((lon+7.5)/15)
    else:
        return int((lon-7.5)/15)

def yrdays(year):
    return 366 if year % 4 == 0 else 365

def dt2daysecs(dt):
    from datetime import datetime
    ndays = (dt - datetime(dt.year,1,1)).days
    nsecs = (dt - datetime(dt.year,dt.month,dt.day,0,0,0)).seconds
    return ndays, nsecs

def UTC_to_local_nv(ndays,nsecs,lon):
    lsecs = nsecs + tz(lon)*3600
    ldays = ndays
    if lsecs > 86400:
        lsecs -= 86400
        ldays +=1
        
    elif lsecs < 0:
        lsecs += 86400
        ldays -= 1
    #if ldays < 0:
    #    print("WARNING: local date is in the previous year!")
    #if ldays > yrdays(year):
    #    print("WANRING: local date is in the next year!")
    
    return ldays, lsecs
    
def UTC_to_local(ndays,nsecs,lon):
    from numpy import vectorize
    vfunc = vectorize(UTC_to_local_nv)
    return vfunc(ndays,nsecs,lon)

def local_to_UTC_nv(ldays,lsecs,lon):
    nsecs = nsecs - tz(lon)*3600
    ndays = ndays
    if nsecs > 86400:
        nsecs -= 86400
        ndays +=1
    elif nsecs < 0:
        nsecs += 86400
        ndays -= 1
        #if ndays < 0:
        #   print("WARNING: UTC date is in the previous year!")
        #if ndays > yrdays(year):
        #   print("WARNING: UTC date is in the next year!")
    return ndays, nsecs 
        
def local_to_UTC(ldays,lsecs,lon):
    from numpy import vectorize
    vfunc = vectorize(local_to_UTC_nv)
    return vfunc(ldays,lsecs,lon)

## Low-accuracy General Solar Position Calculations by NOAA Global Monitoring Division [https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF. The value we need is coszen.

def gamma(ndays,nsecs):
    "Fractional solar year (in radians), which begins on Jan 1st noon. Assumes input in UTC."
    return 2*pi/yrdays(year)*(ndays+(nsecs/3600-12)/24)

def sundec(t):
    "Sun declination (in radians) as a function of fractional solar year (in radians)."
    return 0.006918 - 0.399912*cos(t) + 0.070257*sin(t) \
                    - 0.006758*cos(2*t) + 0.000907*sin(2*t) \
                    - 0.002697*cos(3*t) + 0.00148*sin(3*t)

def eqtime(t):
    "Equation of time (in minutes) as a function of fractional solar year (in radians)."
    # The NOAA version has 0.000075 which, according to Spencer via ... is incorrect.
    return 229.18*(0.0000075 + 0.001868*cos(t)   - 0.032077*sin(t) \
                             - 0.014615*cos(2*t) - 0.040849*sin(2*t))

def coszen_nv(ndays,nsecs,lat,lon):
    "Cosine of the solar zenith angle. Non-vectorized version. Input date and time should be in UTC."
    from numpy import pi
    alat = lat/180*pi
    alon = lon/180*pi
    #decl = sundec(gamma(*local_to_UTC(ndays,nsecs,lon)))
    decl = sundec(gamma(ndays,nsecs))

    #time_offset = eqtime(gamma(*local_to_UTC(ndays,nsecs,lon)))+4*lon-60*tz(lon)
    time_offset = eqtime(gamma(ndays,nsecs))+4*lon-60*tz(lon)
    #tst = lnsecs/60.0 + time_offset
    #print(tz(lon),alon,time_offset,eqtime(gamma(*local_to_UTC(ndays,nsecs,lon))),eqtime(gamma(ndays,nsecs)))
    ldays,lsecs = UTC_to_local(ndays,nsecs,lon) # Which day it is should not matter.
    tst = lsecs/60.0 + time_offset
    ha = (tst/4-180)/180*pi # convert tst from minutes to degrees, then recenter, then convert to radians.
    #thsun = ha/180*pi
    #return sin(alat)*sin(decl)+cos(alat)*cos(decl)*cos(thsun)
    return sin(alat)*sin(decl)+cos(alat)*cos(decl)*cos(ha)
        
def coszen(ndays,nsecs,lat,lon):
    "Cosine of the solar zenith angle. Vectorized version. Inputs date and time should be in UTC."
    from numpy import vectorize
    vfunc = vectorize(coszen_nv)
    return vfunc(ndays,nsecs,lat,lon)

# Requires the solar_utils package from PyPI
def sp_coszen(ldays,lsecs,lat,lon):
    """
    SOLPOS assumes Standard Time input, i.e. local time zone without
     daylight savings. Our function assumes input is in theoretical local timezone.
    """
    from solar_utils import solposAM as solpos
    from datetime import datetime, timedelta
    dt = datetime(year,1,1) + timedelta(days=ldays) + timedelta(seconds=lsecs)
    (angles, airmass) = solpos(location=[lat, lon, tz(lon)],
                               datetime=[dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second],
                               weather=[1013.0, 10.0]) # Standard pressue and temperature, but we could have used 'sp' and 't2m' data.
    zenith, azimuth = angles # in degrees
    return cos(zenith/180*pi)

def coszen_argmax_error_seconds(ldays,lat,lon):
    from numpy import linspace, array, argmax
    ss = linspace(0,86400)
    sp_cc = array([sp_coszen(ldays,s,lat,lon) for s in ss])
    cc = coszen(*UTC_to_local(ldays,ss,lon),lat,lon) # Need to convert to UTC before calling coszen.
    return ss[argmax(sp_cc)]-ss[argmax(cc)]

# Compute our swr using Rosati's formulas.
def swr_nv(ndays,nsecs,lat,lon):
    "Downward solar short wave radiation for the given location and time, which is assumed to be UTC. Non-vectorized version."
    tau = 0.7  # atmospheric transmission coefficient (0.7 for the latitudes of medsea used in Rosati)
    aozone = 0.09 # water vapour plus ozone absorption (0.09 used in Rosati)
    solar = 1370  # solar constant

    # Beware if called with a nsecs > 86400 or < 0 because of adjusting the timezone.
    if nsecs > 86400:
        ndays += 1
        nsecs -= 86400
        
    cz = coszen_nv(ndays,nsecs,lat,lon)
    if cz <= 0:
        cz = 0.0
        qatten = 0.0
    else:
        qatten = tau**(1/cz)
    
    qzer  = cz*solar                       
    
    # Rosati (88) eq. (3.5-3.7)
    qdir  = qzer*qatten
    qdiff = ((1.-aozone)*qzer-qdir)*0.5
    qtot  = qdir+qdiff   
    return qtot    

def swr(ndays,nsecs,lat,lon):
    "Downward solar short wave radiation for the given location and time, which is assumed to be UTC. Vectorized version."
    from numpy import vectorize
#     assert (min(array(nsecs)) >= 0) and (max(array(nsecs))<=86400)
    vfunc = vectorize(swr_nv)
    return vfunc(ndays,nsecs,lat,lon)

def swr_3hourly_mean(ndays,nsecs,lat,lon,timestep=15,method='quadrature'):
    """ Calculate short wave radiation, averaged for the subsequent 3-hourly period.
    
    3-hourly period begins at the given moment (ndays, nsecs), assumed to be given in UTC.
    
    Average taken either over equally spaced samples with duration of 'timestep' in seconds, or estimated by a quadrature.
    """
    from scipy.integrate import quadrature
    
    # Number of samples
    ns = 3*3600/timestep
    assert ns-int(ns) == 0.0
    if ns-int(ns) == 0.0:
        ns = int(ns) 
    else: 
        raise Exception('The timestep given does not divide 3*3600 seconds.')

    ## This method is too slow in Python.
    # Find cumulative sum, then divide by # of samples to get mean.
    cumsum = 0
    if method == 'cumsum':
        for i in range(ns):
            add_secs = timestep*i
            cumsum += swr_nv(ndays+int(nsecs/86400),
                             (nsecs+add_secs)%86400,
                             lat,lon)
            swr_mean = cumsum/ns
    elif method == 'quadrature':
        # CHECK whether swr() assumes LT or UT.
        swr_mean =  quadrature(lambda nsecs: swr(ndays,nsecs,lat,lon),
                               nsecs,nsecs+10800,tol=1e-1,rtol=1e-3)[0]/10800
    else:
        raise NotImplementedError('The requested method of integration: {:s}'.format(method), 'is not implemented.')

    return swr_mean
    

def swr_3hourly_mean_monthly(year,month,lat,lon,method='quadrature'):
    """ 3-hourly mean computed for a month (UTC day 1 of month midnight to UTC day of of next month, midnight) """
    from time import time
    from datetime import datetime
    from numpy import ones

    # The following datetimes are to be interpreted as UTC, wich is also the timezone assumed for all date values used in GOTM 
    # medsea simulations.
    start = datetime(year,month,1)
    stop = datetime(year+1,1,1) if month == 12 else datetime(year,month+1,1)
    nrec = (stop-start).days*8
    ndays_start = (start-datetime(start.year,1,1)).days
    ndays_stop = (stop-datetime(start.year,1,1)).days

    swr_mean = ones((ndays_stop-ndays_start)*8)
    
    for ndays in range(ndays_start,ndays_stop):
        for i in range(8):
            nsecs = i*3*3600
            ndays_of_month = ndays-ndays_start            
            j = ndays_of_month*8+i

            # TAKE CARE BELOW. Check the timezone and interval assumptions of swr_3hourly_mean()
            swr_mean[j] = swr_3hourly_mean(ndays,nsecs,lat,lon,timestep,method=method) 
    return swr_mean

def cloud_factor_local(year,month,lat,lon,swr_obs,method='quadrature'):

    "Compute cloud factor per grid point: lat/lon required and I_0 and I_0_calc expected to be given for the month with same number of records."
    from datetime import datetime
    from time import time
    from netCDF4 import Dataset
    import os
    from numpy import ones

    # The following datetimes are to be interpreted as UTC, wich is also the timezone assumed for all date values used in GOTM 
    # medsea simulations.
    start = datetime(year,month,1)
    stop = datetime(year+1,1,1) if month == 12 else datetime(year,month+1,1)
    nrec = (stop-start).days*8
    ndays_start = (start-datetime(start.year,1,1)).days
    ndays_stop = (stop-datetime(start.year,1,1)).days
    
    #tic = time()
    cloud_factor = ones((ndays_stop-ndays_start)*8) # Defaults to clear sky value.
    swr_mean = swr_3hourly_mean_monthly(year,month,lat,lon,method=method)
    for ndays in range(ndays_start,ndays_stop):
        for i in range(8):
            nsecs = i*3*3600
            ndays_of_month = ndays-ndays_start
            j = ndays_of_month*8+i

            # Instead of checking I_0_ERA, which can be zero for various reasons, check the clear sky value. 
            # If it's non-zero, compute the factor, otherwise, sun below horizon, so we just keep it as the defaul
            # value of 1 as initialized. Well, it should not matter.
            I_0_calc = swr_mean[j]
            I_0_obs = swr_obs[j] # To be taken from swr_ERA[:,m,n], where lat = lat[m], lon = lon[n].
            if I_0_calc > 1: 
                # If I_0_calc is really small but positive, we get into some trouble...
                cloud_factor[j] = I_0_obs / I_0_calc
            if I_0_obs < 0:
                # Use cloud_factor to zero out negative values of data.
                cloud_factor[j] = 0
            if I_0_obs < 1:
                # Maybe we zero out small negligible values as well?
                cloud_factor[j] = 0
    #toc = time()
    # How long does it take for one grid point and one 3-hourly period?
    #print('time elapsed for (m,n) = ({},{}):'.format(m,n), toc-tic)
    
    return cloud_factor, swr_mean

def cloud_factor_all(year,month,grid='1x',use_ipp=True,append_to_ERA_dataset_now=False):
    """ Calculate cloud factor from ERA swrd and internal algorithm, then append a cloud_factor to the ERA dataset. """
    
    import itertools as itt
    swr_ERA = data_sources(year,month,grid=grid)['heat']['swrd'][:]
    lat = data_sources(year,month,grid=grid)['heat']['lat'][:]
    lon = data_sources(year,month,grid=grid)['heat']['lon'][:]

    M = len(lat)
    N = len(lon)
    mm,nn = zip(*itt.product(range(M),range(N)))

    run = lambda m,n: cloud_factor_calc(year,month,m,n,swr_ERA[:,m,n])
    if use_ipp:
        ## Start an ipcluster to speed up.
        rc, lv = prepare_engine()
        dv = rc[:]
        
        # Push these common values to global scope of each engine.
        dv.push(dict(year=year,month=month,swr_ERA=swr_ERA))
        tic()
        results = lv.map(run,mm,nn)
        results.wait()
        toc()
    else:
        results = [run(mm[i],nn[i]) for i in range(M*N)]
    
    ## Append to the ERA_file.
    if not(append_to_ERA_dataset_now):
        return results

    assert append_to_ERA_dataset_now
    with data_sources(year,month,mode='a')['heat'] as ds:
        if 'cloud_factor' in ds.varibles():
            cf = ds['cloud_factor']
        else:
            cf = ds.createVariable('cloud_factor','f8',dimensions=('time','lat','lon'),zlib=True,fill_value=1e+20)
        if 'swrd_clear_sky' in ds.variables():
            swr_cs = ds['swrd_clear_sky']
        else:
            swr_cs = ds.createVariable('swrd_clear_sky','f8',dimensions=('time','lat','lon'),zlib=True,fill_value=1e+20)
        for i,(cloud_factor, swr_mean) in enumerate(results):
            m = mm[i]
            n = nn[i]
            cf[:,m,n] = cloud_factor
            swr_cs[:,m,n] = swr_mean

# New stuff on 2017-11-08.
def solar_times(day,lat,lon,events=['sunrise','snoon','sunset']):
    """Calculate UTC time of sunrise / solar noon / sunset at a given day."""
    from datetime import date, datetime
    if isinstance(day,date) or isinstance(day,datetime):
        year = day.year
        month = day.month
        day = day.day
        ndays= (date(year,month,day)-date(year-1,12,31)).days
    elif isinstance(day,int): # Assume it is directly ndays, not day_of_year = ndays + 1
        ndays = day
    else:
        raise Exception('Unhandled argument type: {!s}'.format(type(day)))
            
    t = gamma(ndays,3600*12) # Use the UTC noon time sun declination, this was not specified in reference.
    decl = sundec(t)
    alat = lat/180*pi
    zen = 90.833/180*pi # The approx correction for atmospheric refraction at sunrise and sunset, and size of solar disk.
    cos_ha = cos(zen)/cos(alat)/cos(decl) - tan(alat)*tan(decl)
    ha = abs(acos(cos_ha))/pi*180 # Convert to degrees
    eqt = eqtime(t)

    UTC_times = dict(
        sunrise = 720 - 4*(lon+ha) - eqt,
        snoon = 720 - 4*lon - eqt,
        sunset = 720 - 4*(lon-ha) - eqt,
        )
    if isinstance(events,list):
        return [UTC_times[event] for event in events]
    else:
        event = events
        assert event in UTC_times.keys()
        return UTC_times[event]
