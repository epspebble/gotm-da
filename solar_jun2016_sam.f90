subroutine short_wave_radiation(jul,secs,lon,lat,swr)
  !
  ! !DESCRIPTION:
  !  This subroutine calculates the short--wave net radiation based on 
  !  latitude, longitude, time, fractional cloud cover and albedo.
  !  The albedo monthly values from \cite{Payne72} are given  here
  !  as means of the values between 
  !  at 30$^{\circ}$ N and 40$^{\circ}$ N for the Atlantic Ocean 
  !  (hence the same latitudinal band of the Mediterranean Sea).
  !
  ! !USES:
  IMPLICIT NONE
  !
  ! !INPUT PARAMETERS:
  integer, intent(in)                 :: jul,secs
  double precision, intent(in)                :: lon,lat
  !
  ! !OUTPUT PARAMETERS:
  double precision, optional, intent(out)     :: swr
  !
  ! !REVISION HISTORY:
  !  Original author(s): Karsten Bolding
  !
!!!!! SH - 2003 - retain fractions of direct and diffuse radiation
!!!!! SH 08/2003 - keep a calculated 'clear sky' value of I_0 
  !
  !  See log for airsea module
  !
  !EOP
  !
  ! !LOCAL VARIABLES:
  !
  !   double precision                  :: solar=1350.
  double precision                  :: eclips=23.439*deg2rad
  double precision                  :: tau=0.63  !Arab=0.74,COARE=0.63,Sub=0.7
  double precision                  :: aozone=0.09


  double precision                  :: th0,th02,th03,sundec
  double precision                  :: thsun,zen,dzen  !,sunbet
  double precision                  :: qatten,qzer,qdir,qdiff,qshort
  double precision                  :: altitude !, qtot
  integer                   :: jab,count1,count2,k
  integer                   :: yy,mm,dd
  double precision                  :: yrdays,days,hour,tjul
  double precision           ::alpha(1:480)

  !SP-cummulative days at each month
  integer                   :: yday(12) = &
       (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)

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

  !  number of days in a year :
  call calendar_date(jul,yy,mm,dd)
  days=float(yday(mm))+float(dd)
  hour=1.0*secs/3600.
  !kbk   if (mod(yy,4) .eq. 0 ! leap year I forgot
  yrdays=365.

  th0 = 2.*pi*days/yrdays
  th02 = 2.*th0
  th03 = 3.*th0
  !  sun declination :
  sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)         &
       - 0.006758*cos(th02) + 0.000907*sin(th02)                 &
       - 0.002697*cos(th03) + 0.001480*sin(th03)

  !  sun declination :
  !SP from http://solardat.uoregon.edu/SolarRadiationBasics.html
  !   sundec = 1.00011 + 0.034221 * cos(th0) + 0.001280 * sin(th0)   &           !                + 0.000719 * cos(th02) + 0.000077 * sin(th02)

  !  sun hour angle :
  thsun = (hour-12.)*15.*deg2rad + alon

  !  cosine of the solar zenith angle (Rosati(88) eq. 3.4 :
  coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)
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

  ! !SP-this represents the possible errors in the derived cloud value
  ! !   if((cloud.LT.0.1).AND.(ABS(swr_error+0.1).LT.0.0001)) then
  ! !       swr_error=0.0
  ! !   else if((cloud.LT.0.1).AND.(ABS(swr_error-0.1).LT.0.0001)) then
  ! !       swr_error=0.2
  ! !   else if ((cloud.GT.0.9).AND.(ABS(swr_error-0.1).LT.0.0001)) then
  ! !       swr_error=0.0
  ! !   else if((cloud.GT.0.9).AND.(ABS(swr_error+0.1).LT.0.0001)) then
  ! !       swr_error=-0.2
  ! !   end if

  ! if(swr_error.ne.0.0) then
  !    if(ABS(cloud-0.001).LT.0.0001) then  !SP IR obs therefore no cloud
  !       cloud=0.0
  !       border=1
  !    else if(ABS(swr_error-9.0).LT.0.0001) then
  !       swr_error=0.0
  !       if(ABS(border).LT.0.0001) then
  !          border=2
  !       end if
  !    else 
  !       cloud=cloud+swr_error
  !       if (cloud.GT.1.0) then
  !          cloud=1.0
  !          border=3
  !       else if (cloud.LT.0.0) then
  !          cloud=0.0
  !          border=3
  !       else
  !          if((ABS(border).LT.0.0001).OR.(ABS(border-2.0).LT.0.0001)) then
  !             border=1
  !          end if
  !       end if
  !    end if
  ! end if
  ! !SP end swr_error section

  if(cloud .lt.0.3) then
     qshort  = qtot*(1.-albedo)        !SP albedo factor needed here
  else
     qshort  = qtot*(1.-.72*cloud + .0019*sunbet)*(1.-albedo) 
  end if

  if (present(swr)) then
     swr = qshort
  else
     I_0_calc = qshort
     I_0_cs   = qtot*(1.-albedo)
  end if

  return
end subroutine short_wave_radiation
