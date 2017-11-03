!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
subroutine short_wave_radiation(jul,secs,alat,alon,swr)
  !
  ! !DESCRIPTION:
  !  This subroutine calculates the short--wave net radiation based on 
  !  latitude, longitude, time, fractional cloud cover and albedo.
  !  The albedo monthly values from \cite{Payne72} are given here
  !  as means of the values between 
  !  at 30$^{\circ}$ N and 40$^{\circ}$ N for the Atlantic Ocean 
  !  (hence the same latitudinal band of the Mediterranean Sea).
  !
  ! !USES:
  use time, only: calendar_date, julian_day, UTC_to_local, tz
  use airsea, only: albedo, solar

  ! WT 20171102, the following are the intended output variables after calling this subroutine.
  use airsea, only: qtot, cosr, coszen, sunbet, qdir_frac, qdiff_frac

  IMPLICIT NONE
  !
  ! !INPUT PARAMETERS:
  integer, intent(in)                 :: jul,secs
  double precision, intent(in)                :: alat,alon !WT assume radian
  !
  ! !OUTPUT PARAMETERS:
  double precision, optional, intent(out)     :: swr
  !   
  ! !REVISION HISTORY:
  !  Original author(s): Karsten Bolding
  !
  ! SH - 2003 - retain fractions of direct and diffuse radiation
  ! SH 08/2003 - keep a calculated 'clear sky' value of I_0 
  !
  !  See log for airsea module
  !
  !EOP
  !
  ! !LOCAL VARIABLES:
  !
  !   double precision                  :: solar=1350.

  ! WT local version of jul and secs, also lon in degrees.   
  !integer                           :: ljul,lsecs
  double precision                  :: lon,lat

  ! WT 20171102 Convenience parameters copied from airsea.f90
  double precision, parameter                :: pi=3.14159265358979323846
  double precision, parameter                :: deg2rad=pi/180.
  double precision, parameter                :: rad2deg=180./pi

  
  ! WT 20171102 Public module variables copied from airsea.f90
  !double precision, public	:: airt
  !double precision, public	:: solar
  !double precision, public	:: I_0_calc,coszen,cosr,qdir_frac,qdiff_frac
  !double precision, public      :: albedo,qtot,sunbet

  double precision                  :: eclips=23.439*deg2rad
  ! tau = 0.7 gives a better fit compared to 0.66 in medsea.
  double precision                  :: tau=0.7  !Arab=0.74,COARE=0.63,Sub=0.7
  double precision                  :: aozone=0.09

  double precision                  :: th0,th02,th03,thsun,solar_time,sundec !, coszen
  double precision                  :: gamma,tst,tst_offset,eqtime,ha,decl!, coszen, NOAA
  double precision                  :: zen,dzen  !,sunbet
  double precision                  :: qatten,qzer,qdir,qdiff,qshort
  double precision                  :: altitude !, qtot
  integer                   :: jab,count1,count2,k
  ! WT 20170315 Modifying new code by SP
  integer                   :: yyyy,mm,dd
  integer                           :: ljul, lsecs !local time
  integer                           :: jul0, jul1 !temp vars
  double precision                  :: yrdays,days,hours !WT renamed hour to hours

  double precision                  :: tjul
  double precision           ::alpha(1:480)

  !SP-cummulative days at each month
  !integer                   :: yday(12) = &
  !     (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)

  !SP-values of albedo taken from Table 1 (Payne,72), with atmospheric transmittance of T=.70
  double precision                  :: alb1(46) = &
       (/.719,.656,.603,.480,.385,.300,.250,.193,.164,.145, &
       .131,.116,.103,.092,.084,.076,.071,.065,.061,.057, &
       .054,.051,.049,.039,.037,.035,.035,.035,.034,.033, &
       .033,.032,.032,.032,.029,.029,.029,.029,.028,.028, &
       .028,.028,.027,.027,.028,.028/)

  !SP-corresponding values of sun altitude
  double precision                  :: alt(46) = &
       (/0.0,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28., &
       30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,52.,54.,56., &
       58.,60.,62.,64.,66.,68.,70.,72.,74.,76.,78.,80.,82.,84.,86., &
       88.,90./)
  !
  !
  !-----------------------------------------------------------------------
  !BOC

  !Convert back to degrees
  lon = alon / deg2rad
  lat = alat / deg2rad

!!! WT 20170316 The sundec formula below is defined from the UTC time zone.
  !   The following effort to to find local values of date and times in hours,
  !   number of seconds since midnight etc... are no longer necessary.
  !
  !
  ! 1. Find the local calendar date: yyyy-mm-dd.
  !print*, "UTM ", jul,secs
  !call UTC_to_local(jul,secs,lon,ljul,lsecs)
  !print *, "jul,secs,lon,ljul,lsecs", jul,secs,lon,ljul,lsecs
  !print*, "local ", ljul,lsecs
  !call calendar_date(ljul,yyyy,mm,dd)
  !call calendar_date(jul,yyyy,mm,dd)
  !
  ! 2. Find the julian day of the final day of last year, save to ljul0.
  !call julian_day(yyyy-1,12,31,ljul0)
  !
  ! 3. Now get the day number of the local day of year.
  !days=float(yday(mm))+float(dd) ! SP's version.
  !days = float(ljul-ljul0) ! int to float, for later formulas
  !print *, "yyyy,mm,dd,days,ljul,ljul0", yyyy,mm,dd,days,ljul,ljul0
  !print *, "days, ljul-ljul0", days, ljul-ljul0 ! Should be the same.

  !kbk   if (mod(yy,4) .eq. 0 ! leap year I forgot
  !yrdays=365. ! GOTM's version.

  ! 4. Find the julian day number of the final day of this year, save to ljul1.
  !call julian_day(yyyy,12,31,ljul1)
  !
  ! 5. Find the number of days in this year.
  !yrdays = float(ljul1-ljul0) ! int to float, for later formulas
  !print *, "yrdays",yrdays
  !
  ! !!! WT End of the previous calculations using local timezone.

!!! Calculation of true solar time to find clear sky value of solar SWR, I_0_calc
  !
  ! Sources:
  ! https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
  ! https://arxiv.org/pdf/1102.3825.pdf

  ! 1. Find calendar date from julian day number in UTC: yyyy-mm-dd, and yrdays of that year
  call calendar_date(jul,yyyy,mm,dd)
  if (mod(yyyy,4) .eq. 0) then
     yrdays = 366
  else
     yrdays = 365
  end if

  ! 2. Find number of days since yyyy-01-01, and the fractional hour since midnight.
  call julian_day(yyyy,1,1,jul0) ! jul0 = julian day number of the first day of year.
  days = jul - jul0 ! = day_of_year-1 = 0 for the first day yyyy:01-01
  hours = 1.*secs/3600. ! hours should be UTC decimal time

  ! 3. The fractional solar year (\gamma) in solareqns.pdf, which begins noon on civil New Year Day. 
  gamma = (2.*pi/yrdays)*(days+((hours-12)/24))  

  ! 4. Finding sun declination, the angle between the equator and sun ray.
  ! The Spencer formula (Spencer, 1971):
  decl = 0.006918 - 0.399912*cos(gamma)    + 0.070257*sin(gamma)        &
       - 0.006758*cos(2.*gamma) + 0.000907*sin(2.*gamma)     &
       - 0.002697*cos(3.*gamma) + 0.001480*sin(3.*gamma)  ! in radians

  !th0 = (2.*pi/yrdays)*(days-1+((hour-12)/24))  ! hour should be UTC decimal time
  !! An alternative (crude)
  !!th0 = 2.*pi*days/yrdays
  !
  !th02 = 2.*th0
  !th03 = 3.*th0
  !sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)         &
  !     - 0.006758*cos(th02) + 0.000907*sin(th02)                 &
  !     - 0.002697*cos(th03) + 0.001480*sin(th03)  ! in radians
  !An alternative
  !sundec = 23.45*pi / 180.*sin(2.*pi*(284. + days) / 365.)

  ! 5. Equation of time, i.e. difference between apparent solar time (by sundials),
  ! and mean solar time (by civil calendar & clock).
  ! * the coefficient 0.0000075 is correct: http://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html
  ! * the erroneous coefficient 0.000075 is reproduced in multiple documents, including the NOAA pdf quoted above.
  !eqtime = 229.18*(0.0000075+0.001868*cos(th0)-0.032077*sin(th0)-0.014615*cos(th02)-0.040849*sin(th02))  ! in minutes
  eqtime = 229.18*(0.0000075 + 0.001868*cos(gamma)    - 0.032077*sin(gamma) &
       - 0.014615*cos(2.*gamma) - 0.040849*sin(2.*gamma))  ! in minutes
  ! An alternative
  !IF ((days.GE.1).AND.(days.LE.106)) THEN
  !    eqtime = -14.2*sin(pi*(days+7.)/111.)
  !ELSE IF ((days.GE.107).AND.(days.LE.166)) THEN
  !    eqtime = 4.0*sin(pi*(days-106.)/59.)
  !ELSE IF ((days.GE.167).AND.(days.LE.246)) THEN
  !    eqtime = -6.5*sin(pi*(days-166.)/80.)
  !ELSE IF ((days.GE.247).AND.(days.LE.365)) THEN
  !    eqtime = 16.4*sin(pi*(days-247.)/113.)
  !END IF


  ! WT Be careful, does hour mean the hour of day, or the fractional number of hours since midnight?
  !solar_time = hour + (eqtime/60) + ((30.-lon)/15) ! in hours
  !solar_time = hour + (eqtime/60) - lon/15 ! in hours (UTC time), lon=degrees
  !solar_time = hour*60 +eqtime + 4*lon - (lsecs-secs)/60 ! The value should not exceed 1440.

  ! 5. Find true solar time in local time at a position with longitude = `lon`.
  ! Along prime meridian, with theoretical timezone tz = 0, the only offset is due to eqtime. Elsewhere, it's the difference between the longitudinal minute and the timezone offset, which are also zeros at multiples of 15 degress if tz is the theoretical timezone only depending on longitude.
  tst_offset = eqtime + 4.0*lon - 60.0*tz(lon)
  call UTC_to_local(jul,secs,lon,ljul,lsecs) ! Get the local time no. of seconds since midnight.
  tst = lsecs/60.0 + tst_offset ! True solar time, in minutes.

  !WT No need for getting the true tst, ha differ by 360 degrees
  ! after adjusting tst by 1440 minutes. It will not change value of cos(ha)
  ! if (tst > 24*60) then
  !    ! ignore the change in calendar date which does not affect the following
  !    tst = tst - 24*60 
  ! else if (tst < 0 ) then
  !    tst = tst + 24*60
  ! end if

  !print *, "lsecs,tst_offset,tst",lsecs,tst_offset,tst
  !print *,"tst,solar_time",tst,solar_time

  !  sun hour angle :
  !   thsun = (hour-12.)*15.*deg2rad + alon
  !thsun = (hour-12.)*15.*deg2rad
  !thsun = (12.-hour)*pi/12.
  !thsun = (12.-solar_time)*pi/12.
  !thsun = pi*((solar_time/12)-1)  ! radians
  !thsun = pi*((solar_time/(4*180))-1)  ! radians

  ! 6. Find the solar hour angle.
  ha = (tst/4-180)*deg2rad ! radians
  !print *,"thsun,ha,thsun-ha",thsun,ha,thsun-ha
  !PRINT*, thsun

  ! 7. Cosine of the solar zenith angle
  !(Rosati(88) eq. 3.4 :
  !coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)
  coszen =sin(alat)*sin(decl)+cos(alat)*cos(decl)*cos(ha) ! ha, decl not needed anymore from now on
  !WT 20170316 Consider moving the above to time.f90 or a new solar.f90 for later reuse (e.g. assimilation
  ! at sunrise or sunset, compare the SST turnaround times with solar sunrise and sunset...)
  !We could also use NREL's solpos C code also if we know how to recompile everything in separate dynamic linkable
  !library files. 


  if (coszen .le. 0.0) then
     coszen = 0.0           !SP-this is when sun is below horizon
     qatten = 0.0           !i.e between dusk and dawn
  else
     qatten = tau**(1./coszen)
  end if

  qzer  = coszen * solar                       
  qdir  = qzer * qatten                        !Rosati(88) eq. (3.5)
  qdiff = ((1.-aozone)*qzer - qdir) * 0.5      !Rosati(88) eq. (3.6)
  qtot  =  qdir + qdiff                        !Rosati(88) eq. (3.7)

!!!!! SH - 2003 - retain fractions of direct and diffuse radiation
  ! ensure fractions remain finite
  if (qtot .gt. 0) then 
     qdir_frac = qdir/qtot
     qdiff_frac = qdiff/qtot
  else
     qdir_frac = 0.
     qdiff_frac = 0.
  end if

!!!!!

  tjul = (days-81.)/yrdays*2.*pi

  !  sin of the solar noon altitude in radians (Rosati, 88) eq. 3.9:
  sunbet=sin(alat)*sin(eclips*sin(tjul))+cos(alat)*cos(eclips*sin(tjul))
  !  solar noon altitude in degrees :
  sunbet = asin(sunbet)*rad2deg

  !  calculates the albedo as a function of sun altitude :
  !  (after Payne jas 1972)
  !  solar zenith angle in degrees :
  zen=(180./pi)*acos(coszen)
  !  sun altitude :
  altitude=90.-zen

  jab=0.5*altitude + 1.

  !linear interpolation
  albedo=alb1(jab)+.5*(alb1(jab+1)-alb1(jab))*(altitude-alt(jab))      

  !SH  calculate cosine of angle of direct refracted entrant radiation 
  cosr = cos(asin((3./4.)*sin(acos(coszen))))

  !  radiation as from Reed(1977), Simpson and Paulson(1979)
  !  calculates SHORT WAVE FLUX ( watt/m*m )
  !  Rosati,Miyakoda 1988 ; eq. 3.8
  !  clouds from COADS perpetual data set


  !170301 Reed formula used here.
  !   if(cloud .lt.0.3) then
  !   if(cloud.eq.0.0) then
  !      qshort  = qtot*(1.-albedo)        !SP albedo factor needed here
  !   else !170301 consider removing this case... over the top?
  !      qshort  = qtot*(1.-.62*cloud + .0019*sunbet)*(1.-albedo)
  !   end if
  ! print *, 'qtot, qshort, albedo', qtot, qshort, albedo
  qshort=qtot*(1.-albedo)

  if (present(swr)) then
     swr = qshort
  else
     !     I_0_calc = qshort
     !     I_0_cs   = qtot*(1.-albedo)
     write(0,*) 'WARNING: short_wave_radiation() called without specifying an output variable name'
  end if

  return
end subroutine short_wave_radiation
     !EOC
