from .config import *

### For SWRD Calcluations

# Temporary global variable
year = 2014

## Helper functions        
def yrdays(year):
    return 366 if year % 4 == 0 else 365

def dt2daysecs(dt):
    from datetime import datetime
    ndays = (dt - datetime(dt.year,1,1)).days
    nsecs = (dt - datetime(dt.year,dt.month,dt.day,0,0,0)).seconds
    return ndays, nsecs

def UTC_to_local_nv(ndays,nsecs,lon):
    local_nsecs = nsecs + tz(lon)*3600
    local_ndays = ndays
    if local_nsecs > 86400:
        local_nsecs -= 86400
        local_ndays +=1
        
    elif local_nsecs < 0:
        local_nsecs += 86400
        local_ndays -= 1
    #if local_ndays < 0:
    #    print("WARNING: local date is in the previous year!")
    #if local_ndays > yrdays(year):
    #    print("WANRING: local date is in the next year!")
    
    return local_ndays, local_nsecs
    
def UTC_to_local(ndays,nsecs,lon):
    from numpy import vectorize
    vfunc = vectorize(UTC_to_local_nv)
    return vfunc(ndays,nsecs,lon)

def local_to_UTC_nv(ndays,nsecs,lon):
    UTC_nsecs = nsecs - tz(lon)*3600
    UTC_ndays = ndays
    if UTC_nsecs > 86400:
        UTC_nsecs -= 86400
        UTC_ndays +=1
    elif UTC_nsecs < 0:
        UTC_nsecs += 86400
        UTC_ndays -= 1
        #if UTC_ndays < 0:
        #   print("WARNING: UTC date is in the previous year!")
        #if UTC_ndays > yrdays(year):
        #   print("WARNING: UTC date is in the next year!")
    return UTC_ndays, UTC_nsecs 
        
def local_to_UTC(ndays,nsecs,lon):
    from numpy import vectorize
    vfunc = vectorize(local_to_UTC_nv)
    return vfunc(ndays,nsecs,lon)

## Low-accuracy General Solar Position Calculations by NOAA Global Monitoring Division [https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF. The value we need is coszen.

def gamma(ndays,nsecs):
    "Fractional solar year (in radians), which begins on Jan 1st noon. Assumes input in UTC."
    return 2*pi/yrdays(year)*(ndays+(nsecs/3600-12)/24)

def sundec(t):
    "Sun declination (in radians) as a function of fractional solar year (in radians)"
    return 0.006918 - 0.399912*cos(t) + 0.070257*sin(t) \
                    - 0.006758*cos(2*t) + 0.000907*sin(2*t) \
                    - 0.002697*cos(3*t) + 0.00148*sin(3*t)

def eqtime(t):
    "Equation of time (in minutes) as a function of fractional solar year (in radians)."
    # The NOAA version has 0.000075 which, according to Spencer via ... is incorrect.
    return 229.18*(0.0000075 + 0.001868*cos(t)   - 0.032077*sin(t) \
                             - 0.014615*cos(2*t) - 0.040849*sin(2*t))

def coszen_nv(ndays,nsecs,lat,lon):
    from numpy import pi
    alat = lat/180*pi
    alon = lon/180*pi
    decl = sundec(gamma(*local_to_UTC(ndays,nsecs,lon)))
    #decl = sundec(gamma(ndays,nsecs))
    #lndays,lnsecs = UTC_to_local(ndays,nsecs,lon) # Which day it is should not matter.
    time_offset = eqtime(gamma(*local_to_UTC(ndays,nsecs,lon)))+4*lon-60*tz(lon)
    #time_offset = eqtime(gamma(ndays,nsecs))+4*lon-60*tz(lon)
    #tst = lnsecs/60.0 + time_offset
    #print(tz(lon),alon,time_offset,eqtime(gamma(*local_to_UTC(ndays,nsecs,lon))),eqtime(gamma(ndays,nsecs)))
    tst = nsecs/60.0 + time_offset
    ha = tst/4-180
    thsun = ha/180*pi
    return sin(alat)*sin(decl)+cos(alat)*cos(decl)*cos(thsun)        
        
def coszen(ndays,nsecs,lat,lon):
    "Cosine of the solar zenith angle. Inputs date and time should be local."
    from numpy import vectorize
    vfunc = vectorize(coszen_nv)
    return vfunc(ndays,nsecs,lat,lon)

## Requires the solar_utils package from PyPI
#def sp_coszen(ndays,nsecs,lat,lon):
#    dt = datetime(year,1,1) + timedelta(days=ndays) + timedelta(seconds=nsecs)
#    (angles, airmass) = solpos(location=[lat, lon, tz(lon)],
#                               datetime=[dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second],
#                               weather=[1015.62055, 40.0]) # Pressure and dry-bulb temp
#    zenith, azimuth = angles # in degrees
#    return cos(zenith/180*pi)
#
#def coszen_argmax_error_seconds(ndays,lat,lon):
#    ss = linspace(0,86400)
#    sp_cc = array([sp_coszen(ndays,s,lat,lon) for s in ss])
#    cc = coszen(ndays,ss,lat,lon)
#    return ss[argmax(sp_cc)]-ss[argmax(cc)]

# Compute our swr using Rosati's formulas.
def swr_nv(ndays,nsecs,lat,lon):
    tau = 0.7  # atmospheric transmission coefficient (0.7 used in Rossatti)
    aozone = 0.09 # water vapour plus ozone absorption (0.09 used in Rossatti)
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
    "Calculate short wave radiation for the given moment, time assumed to be local."
    from numpy import vectorize
#     assert (min(array(nsecs)) >= 0) and (max(array(nsecs))<=86400)
    vfunc = vectorize(swr_nv)
    return vfunc(ndays,nsecs,lat,lon)

def swr_3hourly_mean(ndays,nsecs,lat,lon,timestep,method='quadrature'):
    """ Calculate short wave radiation, averaged for the subsequent 3-hourly period.
    Time period begins at the given moment (ndays, nsecs), assumed to be local.
    Average taken over equally spaced samples with duration of 'timestep' in seconds. 
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
        swr_mean =  quadrature(lambda nsecs: swr(ndays,nsecs,lat,lon),
                               nsecs,nsecs+10800,tol=1e-1,rtol=1e-3)[0]/10800
    else:
        raise NotImplementedError('The requested method of integration: {:s}'.format(method), 'is not implemented.')

    return swr_mean
    

def swr_3hourly_mean_monthly(year,month,m,n,method='quadrature'):
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

    #tic = time()
    swr_mean = ones((ndays_stop-ndays_start)*8)
    
    for ndays in range(ndays_start,ndays_stop):
        for i in range(8):
            nsecs = i*3*3600
            ndays_of_month = ndays-ndays_start            
            j = ndays_of_month*8+i

            # TAKE CARE BELOW.
            # Local time not necssarily at 0, 3, 6, 9 ... etc hours. 
            # Also, swr_3hourly_mean gives the mean over the SUBSEQUENT 3 hours. 
            I_0_calc = swr_3hourly_mean(*UTC_to_local_nv(ndays,nsecs,medsea_lons[n]),
                                        medsea_lats[m],medsea_lons[n],timestep,
                                        method=method) 
            swr_mean[j] = I_0_calc
    #toc = time()
    # How long does it take for one grid point and one 3-hourly period?
    #print('time elapsed for (m,n) = ({},{}):'.format(m,n), toc-tic)
    return swr_mean

def cloud_factor_calc(year,month,m,n,swr_ERA,method='quadrature'):
    from datetime import datetime
    from time import time
    from netCDF4 import Dataset
    import os
    from gotm import medsea_lats, medsea_lons
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
    swr_mean = swr_3hourly_mean_monthly(year,month,m,n,method=method)
    for ndays in range(ndays_start,ndays_stop):
        for i in range(8):
            nsecs = i*3*3600
            ndays_of_month = ndays-ndays_start
            j = ndays_of_month*8+i

            # Instead of checking I_0_ERA, which can be zero for various reasons, check the clear sky value. 
            # If it's non-zero, compute the factor, otherwise, sun below horizon, so we just keep it as the defaul
            # value of 1 as initialized. Well, it should not matter.
            I_0_calc = swr_mean[j]
            I_0_obs = swr_ERA[j,m,n]
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

def cloud_factor_calc_monthly(year,month,use_ipp=True,append_to_ERA_dataset_now=False):
    """ Calculate cloud factor from ERA swrd and internal algorithm, then append a cloud_factor to the ERA dataset. """
    
    import itertools as itt
    ## Get the cloud_factor value for the month
    swr_ERA = data_sources(year,month)['heat']['swrd'][:]

    mm,nn = zip(*itt.product(range(21),range(57)))

    if use_ipp:
        ## Start an ipcluster to speed up.
        rc, lv = prepare_engine()
        dv = rc[:]
        run = lambda m,n: cloud_factor_calc(year,month,m,n,swr_ERA)
        
        # Push these common values to global scope of each engine.
        dv.push(dict(year=year,month=month,swr_ERA=swr_ERA))
        tic()
        results = lv.map(run,mm,nn)
        results.wait()
        toc()
    else:
        results = [cloud_factor_calc(year,month,mm[i],nn[i],swr_ERA) for i in range(21*57)]
    
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
