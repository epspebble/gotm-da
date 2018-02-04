!SP
!This module has the exchange_coefficients subroutine changed to calculate
!fluxes according to the COARE Bulk Flux Algorithm version 3.0a modified 
!slightly by Simon Josey and then adapted by Sam Pimentel in order to fit into
!GOTM. 
! 
!In addition to airsea_solar1.f90 this has the full albedo table from Payne, 
!and I feel this is slightly neater.  Also if the Arabian sea observations for
!incoming solar radiation are used, then the albedo factor needs to be included
!to give a net solar radiation, so this is rectified.
!
!this includes improvements from airsea_solar3.f90, I have added the ability to
!read in the incoming longwave radiation measurements available for the arabian
!sea data.
!09/04
!Improvements have been made to the calculation of the skin temp, it has 
!changed to use previous step values of swr and lwr to calculate skin temp.
!Also the skin temp is now used in the lwr calculation.
!
!changes to the calculation of q, the specific humidity of the air.  Now 
!calculated using RH observations.
!
!The different background radiation methods available are now read in from the !namelists airsea.inp
!11/04
!Albedo factor include for parameterised SWR with cloud less than 0.3
!Minor changes also for when using LWR obs. and when using clouds in LWR-Clarke!
!10/05 : for assimilation at local midnight, run must start local midnight 
!-----------------------------------------------------------------------
!BOP
!
! !GOTM 3.0
!
! !MODULE: airsea --- atmopheric fluxes \label{sec:airsea}
!
! !INTERFACE:
module airsea
  !
  ! !DESCRIPTION:
  !  This module provides various ways to calculate the heat, momentum 
  !  and freshwater fluxes. They are either prescribed as constant values,
  !  see {\tt airsea.inp}, read in from files or calculated by means of
  !  bulk formulae, using observed or modelled meteorological parameters.
  !
  ! !USES:
  use time, only: julian_day, time_diff, calendar_date, write_time_string 
  use time, only: UTC_to_local, tz !WT new function
  use observations, only: read_obs

  !WT Added the following for the copy-and-pasted extinct_method=9, 15, 16 codes
  use observations, only: fsen,zdeta
  
  !
  IMPLICIT NONE
  !  default: all is private.
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public                              :: init_air_sea
  public                              :: air_sea_interaction
  public                              :: set_sst
  public                              :: integrated_fluxes
  public                              :: do_calc_fluxes

  ! !PUBLIC DATA MEMBERS:
  logical,  public                    :: calc_fluxes=.false.
  !  tx and ty are the surface stress components:
  double precision, public                    :: tx,ty
  !  I_0 and heat are short-wave radiation and heat flux at the surface:
  double precision, public                    :: I_0,heat,qb,qh,qe,cloud
  !  sst and sss are sea surface temperature and sea surface salinity:
  double precision, public                    :: sst,sss,skint
  !  integrated short-wave radiation, heat flux and total heat flux:
  double precision, public                    :: int_sw=0.,int_hf=0.
  double precision, public                    :: int_total=0.,int_cs=0.
  !
!!!!! SH - 14/08/2003 - additional parameters for thermal skin etc
  !
  double precision, public	:: I_0_calc,coszen,cosr,qdir_frac,qdiff_frac
  double precision, public	:: deltaTo=0., sssto=0., deltaT=0.
  double precision, public	:: sssT=0., delta=0.,I_0_cs
  double precision, public	:: I_cool = 0.,albedo,qtot,sunbet
  !   double precision, public	:: qb,qe,qh
  double precision, public	:: dlr,w_10

  double precision, public       :: seviri_diff=0.,amsre_diff=0.,tmi_diff=0.
  double precision, public       :: ostia_diff=0.,ostia_seviri_diff=0.
  double precision, public       :: ostia_amsre_diff=0.,ostia_tmi_diff=0.
  double precision, public       :: seviri_sq_diff=0.,amsre_sq_diff=0.
  double precision, public       :: tmi_sq_diff=0.,ostia_sq_diff=0.
  double precision, public       :: ostia_seviri_sq_diff=0.
  double precision, public       :: ostia_amsre_sq_diff=0.
  double precision, public       :: ostia_tmi_sq_diff=0.,seviri_obs=0.
  double precision, public       :: amsre_obs=0.,tmi_obs=0.
  !
  ! !DEFINED PARAMETERS:
  integer, parameter                  :: meteo_unit=20
  integer, parameter                  :: heat_unit=21
  integer, parameter                  :: momentum_unit=22
  integer, parameter                  :: p_e_unit=23
  integer, parameter                  :: sst_unit=24
  integer, parameter                  :: sst_unit2=27
  integer, parameter                  :: sss_unit=25
  integer, parameter                  :: airt_unit=26
  integer, parameter                  :: albedo_unit=28

  double precision, parameter                :: cpa=1004.67 !J/kg/K specific heat of dry air (Businger 1982)
  !   double precision, parameter                :: cp=3995   !3985.
  double precision, parameter                :: emiss=0.98
  double precision, parameter                :: bolz=5.67e-8
  double precision, parameter                :: Kelvin=273.16
  double precision, parameter                :: const06=0.62198
  double precision, parameter                :: pi=3.14159265358979323846
  double precision, parameter                :: deg2rad=pi/180.
  double precision, parameter                :: rad2deg=180./pi

  integer, parameter                  :: CONSTVAL=1
  integer, parameter                  :: FROMFILE=2

!!!!! SH
  integer, parameter                  :: SET_TO_SST=1
  !
  !
  ! !REVISION HISTORY:
  !  Original author(s): Karsten Bolding, Hans Burchard
  !
  !  $Log: airsea_tested.f90,v $
  !Revision 1.1  2003/10/24  10:24:54  helen
  !Initial revision
  !
  !  Revision 1.6  2003/03/28 09:20:34  kbk
  !  added new copyright to files
  !
  !  Revision 1.5  2003/03/28 08:13:47  kbk
  !  removed tabs
  !
  !  Revision 1.4  2003/03/10 08:37:56  gotm
  !  HB fixed the Kondo calculations
  !
  !  Revision 1.3  2001/11/18 11:43:48  gotm
  !  Cleaned
  !
  !  Revision 1.2  2001/06/13 07:40:39  gotm
  !  Lon, lat was hardcoded in meteo.F90 - now passed via init_meteo()
  !
  !  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
  !  initial import into CVS
!!!!! SH - added solar and non solar heat methods
  ! retained original heat method for backwards compatibility
  !
  !EOP
  !
  ! private data members
  integer, public                   :: flux_method
  integer, public                   :: swr_method
  integer, public             	     :: lwr_method
  integer, public                   :: longwave_method
  integer, public                   :: momentum_method
  integer, public                   :: p_e_method
  integer, public                   :: sst_method
  integer, public                   :: sst_method2
  integer, public                   :: sss_method
  integer, public                   :: airt_method
  !WT New stuff as per email with SP 2017-12-22
  integer, public                   :: albedo_method=-999, coolskin_method=-999 ! Implies case default.

  !!HK - make meteo_file public)
  character(len=255), public   :: meteo_file
  character(len=255), public   :: heatflux_file
  character(len=255), public   :: momentumflux_file
  character(len=255), public   :: p_e_flux_file
  character(len=255), public   :: sss_file
  character(len=255), public   :: sst_file
  character(len=255), public   :: sst_file2
  character(len=255), public   :: airt_file
  character(len=255), public   :: albedo_file ! WT When using Ohlmann-Siegel (2000), needs chlo.dat

  double precision                  :: wx,wy
  double precision                  :: wx_obs,wy_obs
  !HK added :
  double precision, public :: u10,v10
  !end of added

  double precision                  :: w
  double precision                  :: airp
!!!!! SH Make public airt - leave twet private to this module
  !   double precision                  :: airt,twet
  double precision                  :: twet
  !   double precision                  :: cloud
  double precision                  :: rh
  double precision                  :: spec_hum,dew_pt
  double precision, public                  :: rho_air
  double precision, public                  :: const_tx,const_ty
  double precision, public                  :: const_qin,const_qout
  double precision, public                  :: swr_error,wind_error  !SP
  double precision, public                  :: wind_h,rh_h,airt_h   !SP
  double precision, public                  :: border,ostia

  double precision                  :: es,ea,e,qs,q,mr,xlv,rnl
  double precision                  :: cdd,chd,ced,Du,zt,dqer

  double precision                  :: alat,alon

!!!!! SH 14/08/2003
  double precision, public	:: airt
  double precision, public	:: solar
  double precision, public	:: net_ir=1353 !1350.

  !WT 2016-09-24 Debug
  character(len=19) :: assim_timestr
  
  !WT 20171026
  integer                           ::iter_max=60 ! was hard-coded to be 30
  real                              ::err_max=0.0001 ! was hard-coded to be 0.001


  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: Initialise the air--sea interaction module
  !
  ! !INTERFACE:
  subroutine init_air_sea(namlst,lat,lon,julday,secs)
    !
    ! !DESCRIPTION:
    !  This routine initialises the air--sea module by reading various variables
    !  from a namelist and open relevant files.
    !
    ! !USES:
    IMPLICIT NONE
    !
    ! !INPUT PARAMETERS:
    integer, intent(in)                 :: namlst
    integer, intent(in)                 :: julday,secs
    double precision, intent(in)                :: lat,lon
    !
    ! !REVISION HISTORY:
    !  Original author(s): Karsten Bolding
    !
    !  See log for airsea module
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    namelist /airsea/ calc_fluxes, &
         meteo_file, &
         flux_method,	  &
         swr_method,		  &
         lwr_method,             &
         longwave_method,             &     !SP
         wind_h,                      &     !SP
         rh_h,                        &     !SP
         airt_h,                      &     !SP
         solar,                       &
         net_ir,  &
         swr_error, &
         wind_error, &
         const_qin,&
         const_qout, &
         heatflux_file, &
         momentum_method, &
         const_tx,&
         const_ty, &
         momentumflux_file, &
         p_e_method,p_e_flux_file, &
         sst_method, sst_file, &
         sst_method2, sst_file2, &
         sss_method, sss_file, &
         airt_method, airt_file, &
         albedo_method, albedo_file, & !WT
         coolskin_method !WT
    
    namelist /coolskin/ iter_max,err_max !WT
    !
    !-----------------------------------------------------------------------
    !BOC
    write(0,*) '   ', 'init_air_sea'

!WT 20171026
    open(namlst,file='coolskin.inp',action='read',status='old',err=901)
    read(namlst,nml=coolskin,err=902)
    close(namlst)
    
999 write(0,*) 'Finished temporary hack for coolskin control.'

    open(namlst,file='airsea.inp',action='read',status='old',err=90)
    read(namlst,nml=airsea,err=91)
    close(namlst)


    !WT 20180120 Edited the behavior when exception occurs: STOP and don't
    ! silently run it!
    !------------------------------------------------------------------------
    select case (albedo_method)
    case (0)
       write(0,*) '       ', 'Using default albedo calculation by table due to Payne (1976).'
       write(0,*) '       ', ' No external data needed.'
    case (1)
       open(albedo_unit,file=albedo_file,action='read',status='old',err=100)
       write(0,*) '       ', 'Using the implicit albedio due to Ohlmann-Siegel (2000).'
       write(0,*) '       ', 'Chlorophyll-a data requied for computing Ohlmann-Siegel (2000) albedo. Reading from:'
       write(0,*) '           ', trim(albedo_file)
    case (2)
       write(0,*) '       ', 'Using four-component albedo by Jin et al. (2011).'
       write(0,*) '       ', '   No external data needed.'
    case default
       print *, 'Unexpected albedo_method = ', albedo_method
       stop 'NotImplementedError'
    end select

    select case (coolskin_method)
    case (0) ! The default Fairall method
    case (1) ! Altare's modification of thickness
    case default
       print *, 'Unexpected coolskin_method = ', coolskin_method
       stop 'NotImplementedError'
    end select
    !------------------------------------------------------------------------

    if (calc_fluxes) then
       open(meteo_unit,file=meteo_file,action='read',status='old',err=92)
       write(0,*) '       ', 'Reading meteo data from:'
       write(0,*) '           ', trim(meteo_file)
    end if

    !   else

    if (flux_method .eq. FROMFILE) then
       open(heat_unit,file=heatflux_file,action='read',status='old',err=93)
       write(0,*) '       ', 'Reading heat fluxes from:'
       write(0,*) '           ', trim(heatflux_file)
    end if

    !     The momentum fluxes
    select case (momentum_method)
    case (FROMFILE)
       open(momentum_unit,file=momentumflux_file,action='read', &
            status='old',err=94)
       write(0,*) '       ', 'Reading momentum fluxes from:'
       write(0,*) '           ', trim(momentumflux_file)
    case default
    end select

    !     The fresh water fluxes
    select case (p_e_method)
    case (FROMFILE)
       open(p_e_unit,file=p_e_flux_file,action='read', &
            status='old',err=95)
       write(0,*) '       ', 'Reading precipitation / evaporation data from:'
       write(0,*) '           ', trim(p_e_flux_file)
    case default
    end select

    !     The sea surface temperature
    select case (sst_method)
    case (FROMFILE)
       open(sst_unit,file=sst_file,action='read',status='old',err=96)
       write(0,*) '       ', 'Reading sea surface temperature from:'
       write(0,*) '           ', trim(sst_file)
       !  call read_sst(julday,secs,sst)
    case default
    end select
    select case (sst_method2)
    case (FROMFILE)
       open(sst_unit2,file=sst_file2,action='read',status='old',err=99)
       write(0,*) '       ', 'Reading sea surface temperature from:'
       write(0,*) '           ', trim(sst_file2)
    case default
    end select

    !     The sea surface salinity
    select case (sss_method)
    case (FROMFILE)
       open(sss_unit,file=sss_file,action='read',status='old',err=97)
       write(0,*) '       ', 'Reading sea surface salinity from:'
       write(0,*) '           ', trim(sss_file)
    case default
    end select

    !     The air temperature
    select case (airt_method)
    case (FROMFILE)
       open(airt_unit,file=airt_file,action='read',status='old',err=98)
       write(0,*) '       ', 'Reading air temperatur from from:'
       write(0,*) '           ', trim(airt_file)
    case default
    end select

    !   end if

    twet=0.
    rh=0.
    spec_hum=0.
    dew_pt=0.
    cloud=0.
    sss=0.
    airt=0.

    alon = deg2rad*lon
    alat = deg2rad*lat

    return

90  write(0,*) 'FATAL ERROR: ', 'I could not open airsea.inp'
    stop 'init_airsea'
91  write(0,*) 'FATAL ERROR: ', 'I could not read airsea namelist'
    stop 'init_airsea'
92  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(meteo_file)
    stop 'init_airsea'
93  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(heatflux_file)
    stop 'init_airsea'
94  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(momentumflux_file)
    stop 'init_airsea'
95  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(p_e_flux_file)
    stop 'init_airsea'
96  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(sst_file)
    stop 'init_airsea'
97  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(sss_file)
    stop 'init_airsea'
98  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(airt_file)
    stop 'init_airsea'
99  write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(sst_file2)
    stop 'init_airsea'
100 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(albedo_file)
    stop 'init_airsea'

!WT 20171026
901 write(0,*) 'WARNING: ', 'I could not open coolskin.inp.\n Resorting to default values.'
    iter_max = 30
    err_max = 0.001
    goto 999
902 write(0,*) 'WARNING: ', 'I could not read coolskin namelist.\n Resorting to default values.'
    iter_max = 30
    err_max = 0.001
    goto 999

    
  end subroutine init_air_sea
  !EOC

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: Obtain the air--sea fluxes
  !
  ! !INTERFACE:
  subroutine air_sea_interaction(jul,secs)
    !
    ! !DESCRIPTION:
    !
    !  Depending on the value of the boolean variable {\tt calc\_fluxes},
    !  the calculation of the fluxes and the short wave radiation are
    !  called or the fluxes are directly read in from the namelist
    !  {\tt airsea.inp} as constants or read in from files. With the present
    !  version of GOTM, the 
    !  surface momentum flux and the surface heat flux can be caluclated.
    !  The surface salinity flux is not coded yet. 
    !
    !  On the long run this will be the routine to call, to calculate and
    !  to obtain air--sea related variables. 
    !
    ! !USES:
    IMPLICIT NONE
    !
    ! !INPUT PARAMETERS:
    integer, intent(in)                 :: jul,secs
    !
    ! !REVISION HISTORY:
    !  Original author(s): Karsten Bolding
    !
    !  See log for airsea module

    !  See airsea module
    !

    !
    ! !LOCAL VARIABLES:
    double precision  :: dummy,adjustment,qb_down,lwrcloud,top,bottom,cloud_factor
    integer           :: i,ios,count=1,count2=1,count3=0
    !
    !
    !EOP
    !-----------------------------------------------------------------------
    !BOC
    !SP moved SST up the order (was after wind_error correction) 28/04/06
    !     The sea surface temperature
    select case (sst_method)
    case (FROMFILE)
       call read_sst(jul,secs,sst)
    case default
       select case (sst_method2)
       case (FROMFILE)
          call read_sst(jul,secs,sst)
       case default
       end select
    end select

    if (calc_fluxes) then 
       call flux_from_meteo(jul,secs)     
    end if

    select case (flux_method)
    case (CONSTVAL)
       heat=const_qout
    case (FROMFILE)
       !WT 20170402 With the precomputation of cloud_factor, the "adjustment" column, which is
       ! the ERA swrd values are not used at all, and they have been used when cloud_factor was
       ! calculated.
       call read_heat_flux(jul,secs,adjustment,cloud_factor,qb_down)
    case default
    end select
    
    select case (swr_method)
    case (CONSTVAL)
       I_0=const_qin
    case (FROMFILE)
       if (flux_method .ne. FROMFILE ) then
          !PRINT*,'ERROR'
          !READ*
          !WT Better error output and should stop the program without holding it up for input.
          print *, "swr_method =", swr_method
          print *, "BUT flux_method =", flux_method
          print *, "The flux_method must equal ", FROMFILE, " as well."
          stop "RuntimeError"
       else
          
          ! WT Calculate "clear-sky downward SWR according to Rosati (88)"
          ! Output saved as I_0_calc. It should not include albedo calculation here. Revision needed.
          call short_wave_radiation(jul,secs,alat,alon,I_0_calc)

          ! WT Use an adjustment factor to find the actuall downward SWR
          !!! WT Possible misnomer below!!! 
          I_0=cloud_factor*I_0_calc  ! cloud_factor comes from heat.dat 2nd column
          cloud = min(1.,max(0.,cloud_factor))  ! fixed fraction cloud values between 0 and 1
          !!! WT Possible misnomer above.
          
          ! WT Some other previous comments not by me.
          !Net shortwave radiation
          !  I_0=I_0_calc
          !     The heat fluxes
          
          !              I_0=adjustment*(1.-albedo) ! adjustment comes from heat.dat 3rd column

          ! SP June 2016 - determine cloud fraction
          !		if(I_0_calc .ne. 0) then
          !			cloud = ((1 - (I_0/I_0_calc) + 0.0019*sunbet)/0.62)
          !			if(cloud .LT. 0.0) then
          !				cloud = 0.0
          !			end if
          !			if(cloud .GT. 1.0) then
          !				cloud = 1.0
          !			end if
          !		else
          !			cloud = 0.0;
          !		end if
          ! SP June 2016 - recalculate swr (now with cloud values)

       end if
    case default
       print *,"swr_method=",swr_method
       stop "NotImplementedError"
    end select

    select case (lwr_method)
    case (CONSTVAL)
       print*,'ERROR'
    case (FROMFILE)
       if (flux_method .ne. FROMFILE ) then
          PRINT*,'ERROR'
          READ*
       else
          qb=(emiss*bolz*(sst+kelvin)**4)-(0.955*qb_down)
       end if
    case default
    end select

    !     The momentum fluxes
    select case (momentum_method)
    case (CONSTVAL)
       tx=const_tx
       ty=const_ty
    case (FROMFILE)
       call read_momentum_flux(jul,secs,tx,ty)
    case default
    end select


    !     The sea surface salinity
    select case (sss_method)
    case (FROMFILE)
    case default
    end select

    !     The air temperature
    select case (airt_method)
    case (FROMFILE)
    case (SET_TO_SST)
       airt=sst
    case default
    end select

    if(calc_fluxes) then
    else
       !   Calculate cool skin effect (Wick, 96) SP: 13/03/06
       call skin_temp(sst,skint) ! WT 2017-10-31 Obsolete parametrization.
    end if

    return
  end subroutine air_sea_interaction
  !EOC
  !-----------------------------------------------------------------------
  subroutine skin_temp(sst,skint)

    ! USES:

    use meanflow,     only:  gravity,rho_0,cp

    IMPLICIT NONE

    double precision, intent(in)   :: sst
    double precision, intent(out)  :: skint

    !LOCAL VARIABLES:
    double precision               :: al,visw,tcw,charn
    double precision               :: u_star,delta_t,zo
    double precision               :: visa,airt
    double precision               :: c_shear,c_conv,Ri_crit,Ri

    !-----------------------------------------------------------------------   

    visw=1.e-6                   !m2/s kinematic viscosity water
    tcw=0.6                      !W/m/K   Thermal conductivity water
    charn=0.011
    al=2.1e-5*(sst+3.2)**0.79      !water thermal expansion coefft.

    airt=sst-1 !airt not known, so guess!
    ! Kinematic viscosity of dry air, in m2/s - Andreas (1989) CRREL Rep. 89-11
    visa=1.326e-5*(1.+6.542e-3*airt+8.301e-6*airt**2-4.84e-9*airt**3)


    !mean of table 4 values in Wick, 1996 JPO
    c_shear=226.5
    c_conv=2.71
    Ri_crit=-1.6e-4

    u_star=sqrt(sqrt(tx*tx+ty*ty)/rho_0)
    zo=charn*u_star*u_star/gravity + 0.11*visa/u_star    !after Smith 1988

    !surface Richardson number as defined by Soloviev and Schlussel, 1994
    Ri=-al*gravity*(-heat-I_0*.1)*visw/(rho_0*cp*(u_star**4))

    if(-heat-I_0*.1.GT.0) then
       delta_t=((-heat-I_0*.1)/(rho_0*cp*sqrt(tcw/(rho_0*cp))))* &
            sqrt(c_shear*sqrt(visw*zo/(u_star**3)) + &
            (c_conv*sqrt(visw*rho_0*cp/(al*gravity*(-heat-I_0*.1))) - &
            c_shear*sqrt(visw*zo/(u_star**3)))*exp(-Ri_crit/Ri))
    else
       delta_t=0.
    end if

    skint=sst-delta_t

  end subroutine skin_temp
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: Finish the air--sea interactions
  !
  ! !INTERFACE:
  subroutine finish_air_sea_interaction
    !
    ! !DESCRIPTION:
    !  Various files are closed in this routine.
    !
    ! !USES:
    IMPLICIT NONE
    !
    ! !REVISION HISTORY:
    !  Original author(s): Karsten Bolding
    !
    !  See log for airsea module
    !
    !EOP
    !-----------------------------------------------------------------------
    !BOC
    if (calc_fluxes) then
       close(meteo_unit)
    else
       if (flux_method .eq. FROMFILE) close(heat_unit)
       if (momentum_method .eq. FROMFILE) close(momentum_unit)
       if (p_e_method      .eq. FROMFILE) close(p_e_unit)
       if (sst_method      .eq. FROMFILE) close(sst_unit)
       if (sst_method2      .eq. FROMFILE) close(sst_unit2)
       if (sss_method      .eq. FROMFILE) close(sss_unit)
       if (airt_method     .eq. FROMFILE) close(airt_unit)
    end if
    return
  end subroutine finish_air_sea_interaction
  !EOC

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: Compute the exchange coefficients \label{sec:calcCoeff}
  !
  ! !INTERFACE:
  subroutine exchange_coefficients()
    !
    ! !DESCRIPTION:
    !  
    !
    ! !USES:

    use meanflow,     only:  gravity,rho_0,cp
    use observations, only: extinct_method,chlo,abp_coe,bb

    IMPLICIT NONE
    !
    ! !DEFINED PARAMETERS:

    !
    ! !REVISION HISTORY:
    !  Original author(s): Sam Pimentel, using code taken from COARE Bulk Flux 
    !                      Algorithm version 3.0a.  See Fairall et al.,2003
    !                      J. Climate. Code then slightly adapted by Simon Josey.
    !
    !  See log for the airsear module
    !
    !WT 20171027 The COARE Bulk-Flux Algorithm is now found at
    !   ftp://ftp.etl.noaa.gov/BLO/Air-Sea/bulkalg/
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    double precision                  ::zu,zq,sst_depth,zs,zi 
    double precision                  ::esatu,qsatu,qa 
    double precision                  ::Rgas,Cpv,al,be,cpw,rhow,von  !,grav 
    double precision                  ::Rns,charn 
    double precision                  ::beta,visw,tcw,bigc,wetc,visa,ta 
    double precision                  ::dter,zo,Wg,Bf,tkt,Dt,Dq
    double precision                  ::u10,usr,zo10,Cd10,Ch10,Ct10,zot10
    double precision                  ::Ct,CC,Ribcu,Ribu,zetu,L10,qsr,tsr
    double precision                  ::rr,zoq,zot,zL,L,psu,pst
    double precision                  ::usrold,tsrold,qsrold,conusr,contsr,conqsr
    double precision                  ::hsb,hlb,qout,dels,qcol,alq,xlamx
    double precision                  ::f_c=0.
! Formula coefficients for extinct_method = 12
double precision                  :: C1(16) = &
(/.026,-.009,-.015,-.003,.063,.278,3.91,16.64, &
.033,-.01,-.019,-.006,.066,.396,7.68,51.27/)
double precision                  :: C2(16) = &
(/.112,.034,-.006,-.131,-.015,-.562,-12.91,-478.28, &
.0,.0,.0,.0,.0,.0,.0,.0/)
double precision                  :: C3(16) = &
(/.0,.0,.0,.0,.0,.0,.0,.0, &
-.025,-.007,-.003,-.004,.006,-.027,-2.49,13.14/)
double precision                  :: C4(16) = &
(/.366,.207,.188,.169,.082,1.02,16.62,736.56, &
.419,.231,.195,.154,.066,.886,17.81,665.19/)

! Parameters for extinct_method 12
double precision		:: para_A,para_K
! Parameters for extinct_method 13
double precision              :: K_ir,K_vis,K1,K2

    integer                           ::iter

    integer                           ::i,j,ios,count=1,count2=1,count3=0

    !
    !
    !-----------------------------------------------------------------------
    !BOC

    !SP-this represents the possible errors in obs
    IF(border.NE.0) THEN
       wx=(1+wind_error)*wx
       wy=(1+wind_error)*wy
    ELSE
       if((wind_error.NE.0.0).AND.(ABS(wind_error-0.25).GT.0.0001).AND.(ABS(wind_error+0.25).GT.0.0001)) then
          if(ABS(wind_error-9.0).LT.0.0001) then
             wind_error=0.0
             border=2
             !         else if(wind_error.LT.-0.999) then
             !            wind_error=-0.999
             !            border=3
             !         else if(wind_error.LT.-0.95) then
             !            wind_error=-0.95
             !            border=3
             !         else if(wind_error.GE.3) then
             !            wind_error=3.0
             !            border=3
          else
             border=1
          end if
       end if
       wx=(1+wind_error)*wx
       wy=(1+wind_error)*wy
    END IF
    !SP end wind errors section

    w = sqrt(wx*wx+wy*wy) ! wind speed

    !  set up fixed instrument levels
    zu=wind_h                 !height of wind measurement Arabian Sea
    zt=airt_h                 !height of air temp Arabian Sea
    zq=rh_h                   !height of humidity Arabian Sea
    sst_depth=0.005           !top layer in GOTM's grid

    !  default values for standard height (normally 10m) pressure and mixed layer height
    zs=10.
    zi=600. 
    airp=airp*.01

    call humidity(airt,airp,ea)         !Teten's returns sat. vapour pressure, at air temp., ea, in mb

    qa=.62197*(ea/(airp-0.378*ea))   !convert vapour pressure to mixing ratio (saturated at air temp.)

    !The below 2 formula are used in relative humidity is observed (rh)
    !      e=ea*rh*0.01                     !relative humidity = vapour pressure/satuated vapour pressure
    !      q=.62197*(e/(airp-0.378*e))     !mixing ratio   (kg/kg)

    !W-H. Tse 160910 debug.      
    !      print *, "e,ea,rh,airp,q"
    !      print *,e,ea,rh,airp,q    


    !W-H. Tse 160910 The information in the link is good, but the following formula was wrongly input by the previous guy,
    !where q and spec_hum is reversed in relationship.

    ! PREVIOUS 
    !The below formula is used if specific humidity is observed (spec_hum)
    !from http://amsglossary.allenpress.com/glossary/search?id=specific-humidity1 using specific humidity (spec_hum) in Kg/Kg
    !      q=spec_hum/(1+spec_hum)           !mixing ratio  (kg/kg)
    ! END PREVIOUS

    ! NEW
    ! The formula by solving for q in spec_hum = q/(1+q), which gives q = spec_hum/(1-spec_hum).

    q = spec_hum/(1-spec_hum)

    ! END NEW      
    rh=q*100/qa                       !relative humidity

    !W-H. Tse 160910 debug.            
    !      print *, "spec_hum,q,qa,rh"
    !      print *, spec_hum,q,qa,rh

    !The below is used if dew point temperature is observed (dew_pt)
    !from http://www.srh.noaa.gov/elp/wxcalc/formulas/vaporPressure.html using dew point temperature (dew_pt) in Celsius
    !      e=6.11*10**(7.5*dew_pt/(237.7+dew_pt))
    !      q=.62197*(e/(airp-0.378*e))     !mixing ratio   (kg/kg)
    !      rh=e*100/ea

    !OPEN(UNIT=77,FILE="OBS/q_with_rh_new2.asc")
    !WRITE(UNIT=77,FMT='(F9.8)') q

    call humidity(sst,airp,es)        !returns saturated vapour pressure at SST, es, in mb
    es=es*0.98                     !reduced for salinity Kraus 1972 p. 46
    qs=.62197*(es/(airp-0.378*es)) !convert from mb to mixing ratio  kg/kg

    ! Constants and coefficients (Stull 1988 p640). 
    Rgas=287.1                    !J/kg/K     gas const. dry air
    Cpv=Cpa*(1.+0.84*q)         !Moist air - currently not used (Businger 1982)
    rho_air=airp*100./(Rgas*(airt+Kelvin)*(1.+0.61*q)) !kg/m3  Moist air density
    !      grav=9.78401528445819         !value of gravity at latitude 15.5 
    al=2.1e-5*(sst+3.2)**0.79      !water thermal expansion coefft.
    be=0.026                      !salinity expansion coefft.
    cpw=cp                        !J/kg/K specific heat water
    rhow=rho_0                    !kg/m3  density water
    von=0.4                       !von Karman's "constant

    ! Factors
    Beta=1.2     !Given as 1.25 in Fairall et al.(1996)

    ! Additional constants needed for cool skin
    visw=1.e-6                   !m2/s kinematic viscosity water
    tcw=0.6                      !W/m/K   Thermal conductivity water
    bigc=16.*gravity*cpw*(rhow*visw)**3/(tcw*tcw*rho_air*rho_air)
    xlv=(2.501-0.00237*sst)*1e+6 !J/kg latent heat of vaporization at sst (3C warming=0.3%)
    wetc=0.622*xlv*qs/(rgas*(sst+Kelvin)**2) !Clausius-Clapeyron
    ! Kinematic viscosity of dry air, in m2/s - Andreas (1989) CRREL Rep. 89-11
    visa=1.326e-5*(1.+6.542e-3*airt+8.301e-6*airt**2-4.84e-9*airt**3) 
    ta=airt+Kelvin      !air temperature K

    ! Initial guesses
    dter=0.3                    !cool skin Dt
    dqer=wetc*dter              !cool skin Dq
    zo=0.0001
    Wg=0.5                      !Gustiness factor initial guess
    tkt= 0.001                  !Cool skin thickness first guess

    ! Air-sea differences - includes warm layer in Dt and Dq
    Du=(w**2.+Wg**2.)**.5       !include gustiness in wind spd. difference
    Dt=sst-airt-0.0098*zt       !potential temperature difference.
    Dq=qs-q                     !mixing ratio difference

    ! **************** neutral coefficients ******************

    u10=Du*dlog(10/zo)/dlog(zu/zo)
    usr=0.035*u10
    zo10=0.011*usr*usr/gravity+0.11*visa/usr
    Cd10=(von/dlog(10/zo10))**2
    Ch10=0.00115
    Ct10=Ch10/sqrt(Cd10)
    zot10=10./dexp(von/Ct10)
    cdd=(von/dlog(zu/zo10))**2

    ! ************* Grachev and Fairall (JAM, 1997) **********

    Ct=von/dlog(zt/zot10)         ! Temperature transfer coefficient
    CC=von*Ct/cdd                  ! z/L vs Rib linear coefficient
    Ribcu=-zu/(zi*0.004*Beta**3)  ! Saturation or plateau Rib 
    Ribu=-gravity*zu*((Dt-dter)+0.61*ta*Dq)/(ta*Du**2)
    if (Ribu.lt.0.) then
       zetu=CC*Ribu/(1.+Ribu/Ribcu)   ! Unstable G and F
    else
       zetu=CC*Ribu*(1.+27./9.*Ribu/CC) ! Stable
    end if
    L10=zu/zetu                       ! MO length

    ! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate zo and z/L

    usr= Du*von/(dlog(zu/zo10)-psiu(zu/L10))
    tsr=-(Dt-dter)*von/(dlog(zt/zot10)-psit(zt/L10))
    qsr=-(Dq-dqer)*von/(dlog(zq/zot10)-psit(zq/L10))

    charn=0.011     !then modify Charnock for high wind speeds Chris' data
    if(Du.gt.10) charn=0.011+(0.018-0.011)*(Du-10)/(18-10)
    if(Du.gt.18) charn=0.018

    ! **** Iterate across u*(t*,q*),zo(zot,zoq) and z/L including cool skin ****

    !(sxj) original loop to nits commented out and convergence test introduced
    do 10 iter=1,iter_max

       zo=charn*usr*usr/gravity + 0.11*visa/usr    !after Smith 1988
       rr=zo*usr/visa

       ! *** zoq and zot fitted to results from several ETL cruises ************

       zoq=min(1.15e-4,5.5e-5/rr**0.6)
       zot=zoq

       zL=von*gravity*zu*(tsr*(1.+0.61*q)+0.61*ta*qsr)/((airt+Kelvin)*usr*usr*(1.+0.61*q))
       L=zu/zL
       psu=psiu(zu/L)
       pst=psit(zt/L)
       dqer=wetc*dter

       !(sxj) store previous usr,tsr,qsr values and calculate convergence stats
       usrold=usr
       tsrold=tsr
       qsrold=qsr

       usr=Du*von/(dlog(zu/zo)-psiu(zu/L))
       tsr=-(Dt-dter)*von/(dlog(zt/zot)-psit(zt/L))
       qsr=-(Dq-dqer)*von/(dlog(zq/zoq)-psit(zq/L))

       !(sxj) calculate convergence stats
       conusr=(usr-usrold)/usrold
       contsr=(tsr-tsrold)/tsrold
       conqsr=(qsr-qsrold)/qsrold

       Bf=-gravity/ta*usr*(tsr+0.61*ta*qsr)
       if (Bf.gt.0) then
          Wg=Beta*(Bf*zi)**.333
       else
          Wg=0.2
       end if
       Du=sqrt(w**2.+Wg**2.)        !include gustiness in wind spd.
       !      
       !use net swr and net lwr calculated from previous step
       rns=I_0
       rnl=qb   
       !   Cool skin
       hsb=-rho_air*cpa*usr*tsr
       hlb=-rho_air*xlv*usr*qsr
       qout=rnl+hsb+hlb

       ! calculate the fraction of solar radiation absorbed in the skin layer
       ! this depends on extinction method
        select case (extinct_method)

        case(9)
            !Eq. (16) from Fairall et al., 1996b, using the full Paulson & Simpson (1981) 9-band scheme
            do j=1,9
                f_c = f_c + fsen(j)*( 1 - (zdeta(j)/tkt)*(1-exp(-tkt/zdeta(j))) )
            end do
            dels=rns*f_c
            f_c=0.

        case(12)
            !Eq. (7) in Ohlmann & Siegel (2000)
            if(coszen.lt. 0.2588) then
                coszen=0.2588
            end if
            if(cloud.gt.0.1) then
                do j=1,4
                    para_A=C1(j)*chlo+C2(j)*cloud+C4(j)
                    para_K=C1(j+4)*chlo+C2(j+4)*cloud+C4(j+4)
                    f_c = f_c + ( para_A*tkt - (para_A/para_K)*(1-exp(-para_K*tkt)) )/tkt
                end do
            else
                do j=9,12
                    para_A=C1(j)*chlo+(C3(j)/coszen)+C4(j)
                    para_K=C1(j+4)*chlo+(C3(j+4)/coszen)+C4(j+4)
                    f_c = f_c + ( para_A*tkt - (para_A/para_K)*(1-exp(-para_K*tkt)) )/tkt
                end do
            end if
            dels=rns*f_c
            f_c=0.

        case(13)
            !uses Eq. (11) in Lee et al., 2005
            K1 = (-0.057+0.482*sqrt(abp_coe)+4.221*bb)*(1+0.09*sin(acos(coszen)))
            K2 = (0.183+0.702*abp_coe-2.567*bb)*(1.465-0.667*coszen)

            K_vis = K1+K2/sqrt(1+(tkt/2))
            K_ir = (0.560+2.304/(0.001+(tkt/2))**0.65)*(1+0.002*acos(coszen)*180/(3.1415926))
            
            f_c = 1 - ( 0.576*exp(-K_ir*(tkt/2)) + 0.424*exp(-K_vis*(tkt/2)) ) ! crude numerical approx. of integral, i.e. 1-T(z=tkt/2)
            dels=rns*f_c
            f_c=0.

        case(15)
            !Eq. (16) from Fairall et al., 1996b, using the full Paulson & Simpson (1981) 9-band scheme
            zdeta(1)=24.1;  ! Jerlov water type I
            zdeta(2)=0.673; ! Jerlov water type I
            !zdeta(1)=11;    ! Jerlov water type II
            !zdeta(2)=0.664; ! Jerlov water type II
            !zdeta(1)=6.54;  ! Jerlov water type III
            !zdeta(2)=0.654; ! Jerlov water type III
            !zdeta(1)=3.36;  ! Jerlov water type 1
            !zdeta(2)=0.653; ! Jerlov water type 1
            !zdeta(1)=2.27;  ! Jerlov water type 3
            !zdeta(2)=0.645; ! Jerlov water type 3
            !zdeta(1)=1.56;  ! Jerlov water type 5
            !zdeta(2)=0.627; ! Jerlov water type 5
            !zdeta(1)=1.13;  ! Jerlov water type 7
            !zdeta(2)=0.606; ! Jerlov water type 7
            !zdeta(1)=0.736;  ! Jerlov water type 9
            !zdeta(2)=0.58; ! Jerlov water type 9
            do j=1,9
                f_c = f_c + fsen(j)*( 1 - (zdeta(j)/tkt)*(1-exp(-tkt/zdeta(j))) )
            end do
            dels=rns*f_c
            f_c=0.

        case(16)
            !Eq. (16) from Fairall et al., 1996b, using the full Paulson & Simpson (1981) 9-band scheme
            zdeta(1)=1/0.066;  ! Jerlov water type I
            !zdeta(1)=1/0.076;  ! Jerlov water type IA
            !zdeta(1)=1/0.088;  ! Jerlov water type IB
            !zdeta(1)=1/0.132;  ! Jerlov water type II
            !zdeta(1)=1/0.382;  ! Jerlov water type III
            !zdeta(1)=1/0.49;   ! Jerlov water type 1
            !zdeta(1)=1/0.7;    ! Jerlov water type 3
            !zdeta(1)=1/1.0;    ! Jerlov water type 5
            !zdeta(1)=1/1.09;   ! Jerlov water type 7
            !zdeta(1)=1/1.6;    ! Jerlov water type 9
            do j=1,9
                f_c = f_c + fsen(j)*( 1 - (zdeta(j)/tkt)*(1-exp(-tkt/zdeta(j))) )
            end do
            dels=rns*f_c
            f_c=0.0

        case default

            !Eq. (17) from Fairall et al., 1996b (an approximation of Paulson & Simpson (1981) 9-band scheme)
            !f_c = .137+11.*tkt-6.6e-5/tkt*(1.-dexp(-tkt/8.0e-4))

            !Eq. (17) from Fairall et al., 1996b, first term is changed based on suggestion by Wick et al., 2005 (who used the Ohlmann & Siegel (2000) scheme)
            f_c = .067+11.*tkt-6.6e-5/tkt*(1.-dexp(-tkt/8.0e-4))

            dels=rns*f_c
            f_c=0.

        end select

!------ Coolskin parametrization options  ------------------------------------        
        qcol=qout-dels
        
        if (qcol.le.0) then
           !no cool-skin, we do not model a possible warm-skin
           dter = 0
        else
           select case (coolskin_method) 
           case (1)
              !---------------------------------------- Artale et al (2002)    HX
              !WT Split casing w.r.t to wind speed.
              if (w.le.7.5) then
                 xlamx = (sqrt(rho_air/rhow)*usr)*von*86400/((0.2*w+0.5)*rhow*10*cpw*visw)
              end if
              if (w.gt.7.5.and.w.lt.10) then
                 xlamx = (sqrt(rho_air/rhow)*usr)*von*86400/((1.6*w-10)*rhow*10*cpw*visw)
              end if
              if (w.ge.10) then
                 xlamx = (sqrt(rho_air/rhow)*usr)*von*86400/(6*rhow*10*cpw*visw)
              end if
              tkt=xlamx*visw/(sqrt(rho_air/rhow)*usr)
              
              !WT Seem to be common to Fairall et al (1996b) and Artale et al (2002).
              dter=qcol*tkt/tcw                                 ! Cool skin (Eq. 13 in Fairall et al. (1996b))
              !------------------------end Artale et al (2002)------------

!---------------------------------------------------PS81  !HX
!        xlamx = 6.5
!        tkt=xlamx*visw/(sqrt(rho_air/rhow)*usr)
!-----------------------end PS81-------------------------

!---------------------------------------------------W85   !HX
!        if (u10.le.7) then
!           xlamx = 2+(5/7)*u10
!        else
!           xlamx = 7.
!        end if
!        tkt=xlamx*visw/(sqrt(rho_air/rhow)*usr)
!----------------------end W85--------------------------

           case default
             !---------------------------------------- Fairall et al (1996b)
             alq=Al*qcol+be*hlb*cpw/xlv                      !Eq. 8 in Fairall et al. (1996b)
             if(alq.gt.0.) then                              !originally (qcol.gt.0)
                xlamx=6./(1.+(bigc*alq/usr**4)**.75)**.333      !Eq 14 Saunders coeff.
                tkt=xlamx*visw/(sqrt(rho_air/rhow)*usr)          !Eq.12 Sublayer thickness
             else
                xlamx=6.                                      !prevent excessive warm skins
                tkt=min(.01,xlamx*visw/(sqrt(rho_air/rhow)*usr)) !Limit tkt
             end if
             dter=qcol*tkt/tcw                                 ! Cool skin (Eq. 13 in Fairall et al. (1996b))
             !------------------------end Fairall et al (1996b)------------
            
         end select

         end if
         
       dqer=wetc*dter


       !WT Following should be immediately after computations of conusr, contsr, conqsr
       ! for readability up to the second 'continue' statement (jump out point)
       
       !(sxj) check for convergence and leave loop if met
       IF((iter==iter_max).AND.(max(abs(conusr),abs(contsr),abs(conqsr)).gt.err_max)) THEN
          PRINT*,'convergence error: consur, constr, conqsr = ', conusr, contsr, conqsr
       END IF
       if (max(abs(conusr),abs(contsr),abs(conqsr)).lt.err_max) then
          goto 912
       end if
10     continue                                           ! end iterations

       !(sxj) jump out point
912    continue

       ! compute surface fluxes and other parameters
       skint=sst-dter                    !final skin temperature this timestep

       ! compute transfer coefficients
       cdd=(USR/Du)**2
       chd=USR*TSR/(Du*(airt-skint+.0098*zt)) 
       ced=USR*QSR/(Du*(Q-QS+dqer))
       return 

     end subroutine exchange_coefficients
     !EOC
     !------------------------------------------------------------------
     subroutine humidity(T,P,esat)                                 

       ! Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 

       double precision :: T,P,esat

       esat = (1.0007+3.46e-6*P)*6.1121*dexp(17.502*T/(240.97+T)) !mb
       return
     end subroutine humidity

     !------------------------------------------------------------------
     function psiu(zL)

       ! psiu and psit evaluate stability function for wind speed and scalars
       ! matching Kansas and free convection forms with weighting f
       ! convective form follows Fairall et al (1996) with profile constants
       ! from Grachev et al (2000) BLM
       ! stable form from Beljaars and Holtslag (1991)

       double precision :: zL,x,y,psik,psic,f,psiu,c
       if(zL.lt.0) then
          x=(1-15.*zL)**.25                        !Kansas unstable
          psik=2.*dlog((1.+x)/2.)+dlog((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.)
          y=(1.-10.15*zL)**.3333                   !Convective
          psic=1.5*dlog((1.+y+y*y)/3.)-sqrt(3.)*atan((1.+2.*y)/sqrt(3.)) &
               +4.*atan(1.)/sqrt(3.)
          f=zL*zL/(1.+zL*zL)
          psiu=(1.-f)*psik+f*psic
       else
          c=min(50.,0.35*zL)                       !Stable
          psiu=-((1.+1.*zL)**1.+.6667*(zL-14.28)/dexp(c)+8.525)
       end if
       return
     end function psiu

     !--------------------------------------------------------------  
     function psit(zL)
       double precision :: zL,x,y,psik,psic,f,psit,c
       if(zL.lt.0) then
          x=(1-15.*zL)**.5                          !Kansas unstable
          psik=2.*dlog((1.+x)/2.)
          y=(1.-34.15*zL)**.3333                    !Convective
          psic=1.5*dlog((1.+y+y*y)/3.)-sqrt(3.)*atan((1.+2.*y)/sqrt(3.))   &
               +4.*atan(1.)/sqrt(3.)
          f=zL*zL/(1.+zL*zL)
          psit=(1.-f)*psik+f*psic
       else
          c=min(50.,0.35*zL)                        !Stable
          psit=-((1.+2.*zL/3.)**1.5+.6667*(zL-14.28)/dexp(c)+8.525)
       end if
       return
     end function psit
     !-----------------------------------------------------------------------

     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Calculate the heat fluxes \label{sec:calcFluxes}
     !
     ! !INTERFACE:
     !   subroutine do_calc_fluxes(heatf,taux,tauy) 
     subroutine do_calc_fluxes(qb,qh,qe,taux,tauy)
       !
       ! !DESCRIPTION:
       !  The latent and the sensible heat flux, the long-wave back
       !  radiation (and thus the total net surface heat flux) and
       !  the surface momentum flux are calclated here. For the
       !  long--wave back radiation, the formulae of \cite{Clarketal74} 
       !  and \cite{HastenrathLamb78} may be used as alternatives. 
       !
       ! !USES:
       IMPLICIT NONE
       !
       ! !OUTPUT PARAMETERS:
       double precision, optional, intent(out)     :: qb,qh,qe,taux,tauy

       !
       ! !DEFINED PARAMETERS:
       !   integer, parameter        :: clark=1  ! Clark et. al, 1974
       !   integer, parameter        :: hastenrath=2  ! Hastenrath and Lamb, 1978
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module
       ! SH Nov 2002
       !! Added modifications to include cloud physics (unfinished)
       !  See airsea module
       !
!!!!! SH - added parameterisations for long wave radiation, 3,4,5,6
       ! 6 - based on parameterisation in skin_2_bulk.pro
       !
       !
       ! !LOCAL VARIABLES:
       double precision                  :: tmp
       real                              :: clark_lambda

!!!!! SH
       double precision		:: flwclr ! (clearsky component Swinbank,1963)
       double precision		:: airtk ! air temperature in Kelvin
       double precision		:: cl_tk ! cloud temperature in Kelvin
       double precision		:: alr, cl_ht, n_eff  ! 
       double precision		:: lw_in, lw_out ! 
       double precision		:: q10,dtemp

       !   integer, parameter	:: clark=1	! Clark et. al, 1974
       !   integer, parameter	:: hastenrath=2	! Hastenrath and Lamb, 1978
       !   integer, parameter	:: zillman=3    ! Zillman, 1972 
       !   integer, parameter	:: Brunt_type=4    ! Coud physics
       !   integer, parameter   :: binami=5      !Binami et al, 1995
       !   integer, parameter   :: josey=6       !Josey and Pascal, 2003

       !
       !EOP

       !BOC
       qh=cpa*rho_air*chd*Du*(skint-airt-.0098*zt)     !sensible W/m2
       qe=xlv*rho_air*ced*Du*(qs-q-dqer)               !latent W/m2 

       !OPEN(UNIT=77,FILE="inputs_detail.asc")
       !WRITE(UNIT=77,FMT='(F18.16,2X,F5.2,2X,F6.3,2X,F18.13)') w,airt,sst,airp

       tmp=skint+Kelvin
       select case(longwave_method)    ! back radiation

       case(1)
          !values from Table 3 in Josey et al, 1997, JGR
          if ((alat/deg2rad.GT.-2.5).AND.(alat/deg2rad.LT.2.5)) then
             clark_lambda=.51
          else if ((alat/deg2rad.GT.-7.5).AND.(alat/deg2rad.LT.7.5)) then
             clark_lambda=.53
          else if ((alat/deg2rad.GT.-15).AND.(alat/deg2rad.LT.15)) then
             clark_lambda=.56
          else if ((alat/deg2rad.GT.-25).AND.(alat/deg2rad.LT.25)) then
             clark_lambda=.6
          else if ((alat/deg2rad.GT.-35).AND.(alat/deg2rad.LT.35)) then
             clark_lambda=.64
          else if ((alat/deg2rad.GT.-45).AND.(alat/deg2rad.LT.45)) then
             clark_lambda=.69
          else if ((alat/deg2rad.GT.-55).AND.(alat/deg2rad.LT.55)) then
             clark_lambda=.73
          else if ((alat/deg2rad.GT.-65).AND.(alat/deg2rad.LT.65)) then
             clark_lambda=.77
          else if ((alat/deg2rad.GT.-75).AND.(alat/deg2rad.LT.75)) then
             clark_lambda=.81
          else
             clark_lambda=.85
          end if

          qb=(1.0-clark_lambda*cloud*cloud)                                     &
               *emiss*bolz*(tmp**4)*(0.39-0.05*sqrt(es*rh/100.))          &
               +4.0*emiss*bolz*(tmp**3)*(sst-airt)


       case(2)                    ! qa in g(water)/kg(wet air)
          qb=(1.0-.8*cloud*cloud)                                     &
               *emiss*bolz*(tmp**4)*(0.39-0.056*sqrt(1000*q))          &
               +4.0*emiss*bolz*(tmp**3)*(sst-airt)

!!!!! SH - added options
       case(3)
          dlr= ((9.2e-6)*(airtk**2)*bolz*(airtk**4) &
               -0.96*n_eff*bolz*(airtk**4)*(1-(9.2e-6)*(airtk**2))) 
          airtk=airt+Kelvin
          qb=emiss*bolz*(tmp**4) - dlr        ! Zillman, 1972
          ! From IDL code
       case(4)      !(QJR met soc. (58), pp 389-420, year: 1932)
          airtk=airt+Kelvin

          !        q10 = RH * 610.8*exp(19.85*(1.0 - 273.16 / airtk)) / (554.0 * airtk)
          ! gas constant for moist air
          !        Rgma = 8.31436 / (( 1.0 - q10 ) * 0.028966 + q10 * 0.018016) 
          q10 = 0.00001 * rh * &
               (6.1094*exp(17.625*(airtk-273.16)/(243.04+airtk-273.16))*1.00071*exp(0.0000045*(airp)))* &
               (100000.0/(8.31451*airtk/0.018016)) / rho_air
          dlr=(0.24+4.33*sqrt(q10))*bolz*(airtk**4)

          qb=emiss*bolz*(tmp**4) - dlr

       case(5)         !JGR, 1995
          airtk=airt+Kelvin
          dlr=bolz*(0.684 +0.0056*(es*rh/100.))*(1.+0.1762*(cloud**2))*(airtk**4)
          qb=emiss*bolz*(tmp**4) - (1.-0.045)*dlr
       case(6)      !Josey and Pascal JGR 2003
          airtk=airt+Kelvin
          dtemp=34.07 +4157/LOG(2.1718e8/(es*rh/100.))   !dew point temp 
          !(Henderson-Sellers,1984)
          dlr=(1.-0.045)*bolz* &
               (airtk+10.77*(cloud**2)+2.34*cloud-18.44+0.84*(dtemp-airtk+4.01))**4 
          qb=emiss*bolz*(tmp**4) - dlr

       case default

       end select

       !   if(present(heatf)) then
       !     heatf = -(qe+qh+qb)
       !   else
       !     heat = -(qe+qh+qb)
       !   end if

       tmp=rho_air*cdd*Du                            !stress N/m2


       !OPEN(UNIT=78,FILE="data_feb_hourly_mean.asc")
       !WRITE(UNIT=78,FMT='(F5.2,2X,F17.14,2X,2(F6.2,2X),F18.13,2X,F5.2,2X,F17.12,2x,F17.13)') &
       !     airt,sst,wx,wy,airp,rh,I_0,rnl
       !OPEN(UNIT=79,FILE="fluxes_feb_hourly_mean_data.asc")
       !WRITE(UNIT=79,FMT='(F17.13,2X,F17.13,2X,F17.13,2X,F17.14)') & 
       !     qe,qh,qb,tmp


       if(present(taux)) then
          taux  = tmp*wx
       else
          tx = tmp*wx
       end if
       if(present(tauy)) then
          tauy  = tmp*wy
       else
          ty = tmp*wy
       end if

       return
     end subroutine do_calc_fluxes

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
       use observations, only: chlo
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
       !-----------------Albedo Options-------------------------HX
       
       double precision             ::fmiusigma,alphasdif,alphasdir,rpara,rperp,Rtotal,sigma,n_a,n_o
       ! WT EI = 'extraterrestrial radiation', k_t 'clearness index' for ERBS model code by SP, see https://plantpredict.com/algorithm/irradiance-radiation/#erbs-model for reference.
       double precision             :: f_dir, f_dif, k_t, EI ! WT Note f_dir is related to qdir / qtot from Rosati's formulas. 
       double precision             :: sinzen, sinzent, coszent !WT Extra auxiliaries for clarity
        double precision         ::sigmasquare,wind,fwc,albedoe,para_A,Trans_1
        integer                   ::j
        double precision                  :: C1(16) = &
            (/.026,-.009,-.015,-.003,.063,.278,3.91,16.64, &
            .033,-.01,-.019,-.006,.066,.396,7.68,51.27/)
        double precision                  :: C2(16) = &
            (/.112,.034,-.006,-.131,-.015,-.562,-12.91,-478.28, &
            .0,.0,.0,.0,.0,.0,.0,.0/)
        double precision                  :: C3(16) = &
            (/.0,.0,.0,.0,.0,.0,.0,.0, &
            -.025,-.007,-.003,-.004,.006,-.027,-2.49,13.14/)
        double precision                  :: C4(16) = &
            (/.366,.207,.188,.169,.082,1.02,16.62,736.56, &
            .419,.231,.195,.154,.066,.886,17.81,665.19/)
!-------------------------------------------HX
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

       !SH  calculate cosine of angle of direct refracted entrant radiation 
       cosr = cos(asin((3./4.)*sin(acos(coszen)))) !WT An output as a public variable for variants of 9-band Paulson et al. light extinction.
       
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
      
! --------------------------------- Albedo parametrization options ---------------
       select case (albedo_method)
          
       case (0)
          !------------------------------------- Payne (1972)
          !  calculates the albedo as a function of sun altitude :
          !  (after Payne jas 1972)
          !  solar zenith angle in degrees :
          zen=(180./pi)*acos(coszen)
          !  sun altitude :
          altitude=90.-zen
          
          jab=0.5*altitude + 1.
          
          !linear interpolation
          albedo=alb1(jab)+.5*(alb1(jab+1)-alb1(jab))*(altitude-alt(jab))      
          !-------------------------------------- end Payne (1972)

       case (1)
          !---------------------------------------- Ohlmann & Siegel (2000)
          !WT The O-S. formula implicitly includes surface albedo, and when evaluated at z=0 gives surface albedo.

          ! code by HX
          if(cloud.gt.0.1) then
             do j=1,4
                para_A=C1(j)*chlo+C2(j)*cloud+C4(j)
                Trans_1 = Trans_1 + para_A
             end do
          else
             do j=9,12
                if(coszen.lt. 0.2588) then
                   coszen=0.2588
                end if
                para_A=C1(j)*chlo+(C3(j)/coszen)+C4(j)
                Trans_1 = Trans_1+para_A
             end do
          end if
          albedo = 1 - Trans_1

          ! WT The above should be equivalent to the following, if the function is properly imported from
          ! the temperature subroutine under meanflow module. But we better wait until a better organization
          ! of code.
          !albedo = trans(0,extinct_method)
          !---------------------------------------- end Ohlmann & Siegel (2000)

       case (2)
          !---------------------------------------- Jin et al. (2011)
          !new albedo calculation for broadband from Jin(2011) (HX 30/05/2017)
          n_a = 1.               !refractive index of air
          n_o = 1.34             !refractive index of seawater

          !WT According to last rows of Fig. 8, Fig. 9 in Jin et al.  (2011),
          ! when sky is clear, f_dir is close to 1, and
          ! when overcast, f_dir is could be close to 0.
          ! Compare last rows of Fig.8, Fig.9.
          !f_dir = 0.7            !coefficient for direct radiance
          !f_dif = 0.3            !coefficient for diffuse radiance
          ! WT An alternative is to use Rosati (88) formula and put f_dir = qdir / qtot, f_dif = 1 - f_dir
          
          !SP compute fraction of diffuse radiation
          ! Erbs, et al, Estimation of the diffuse radiation fraction for hourly, daily and monthly-average global radiation. Solar Energy, 28:293-302, 1982.
          ! extraterrestrial irradiance:
          EI = solar*(1.00011 + 0.034221*cos(gamma) + 0.00128*sin(gamma) + 0.000719*cos(2.*gamma) + 0.000077*sin(2.*gamma))
          ! clearness index:
          !WT I_0 yet to be defined, which will have a cloud factor applied, but NOT the albedo (circular reference).
          !k_t = I_0/(EI*coszen)

          !WT Replacing I_0 by qtot*cloud to be the (Global Horizontal Irradiance (i.e. total downward swr in our case))
          !WT Note that cloud is restriction of cloud_factor in [0,1] but it is 1 when there is clear-sky, i.e. no cloud, misnomer!!!
          k_t = qtot*cloud/(EI*coszen) !WT Suggested change that accounts for the overall reduction of radiation.
          
          if (k_t.le.0.22) then
             f_dif = 1. - 0.09*k_t
          else if (k_t.gt.0.22 .and. k_t.le.0.8) then
             f_dif = 0.9511-0.1604*k_t + 0.4388*k_t**2. - 16.638*k_t**3. + 12.338*k_t**4.
          else
             f_dif = 0.165
          end if
          f_dir = 1. - f_dif
          ! end compute fraction of diffuse radiation

          if (qtot.gt.0) then
             if (coszen==0.0) then
                coszent=0.0
                albedo=0.0
             else
                ! WT This seems overly complicated, e.g. acos is nonnegative, why abs?.
                ! I think abs(sin(abs(acos(coszen)))) = sqrt(1-coszen**2)
                !coszent = abs(cos(abs(asin(abs(sin(abs(acos(coszen))))*n_a/n_o))))
                ! WT 20180120 Make the formula more transparent...
                sinzen = sqrt(1-coszen**2)
                sinzent = sinzen*n_a/n_o
                coszent = sqrt(1-sinzent**2)
                
                rpara = (n_a*coszen-n_o*coszent)/(n_a*coszen+n_o*coszent)     !Fresnel's equations for reflection
                rperp = (n_o*coszen-n_a*coszent)/(n_o*coszen+n_a*coszent)
                
                Rtotal = 0.5*(rpara**2. + rperp**2.)    !Unpolarised light
                wind = sqrt(wx_obs**2.+wy_obs**2.)
                sigma = sqrt(0.003+0.00512*wind)                  !eq.2
                !WT The following case splitting seem unnecessary and also inconsistent with the IDL code from Jin (2001)'s paper.
                ! if (wind.eq.0) then
                !    sigmasquare = 0
                !    sigma = 0
                ! else
                !    sigmasquare = 0.003+0.00512*wind                  !eq.2
                !    sigma = sqrt(sigmasquare)
                !end if
                
                fmiusigma = (0.0152-1.7873*coszen+6.8972*coszen**2.0  &
                     - 8.5778*coszen**3.0+4.071*sigma-7.6446*coszen*sigma) &
                     *exp(0.1643-7.8409*coszen-3.5639*coszen**2.0-2.3588*sigma  &
                     +10.0538*coszen*sigma)                                         !eq.4
                
                alphasdir = Rtotal-fmiusigma
                if (cloud.gt.0.8) then !WT Revision on when to split may be needed.
                   ! WT For "clear-sky" according to Jin et al.
                   alphasdif = -0.1482-0.012*sigma+0.1608*n_o-0.0244*n_o*sigma    !eq.5a !WT Modified last coefficient from -0.0193 to -0.0244 (the equation is found on p.5)
                else ! WT For "overcast" condition according to Jin et al.
                   alphasdif = -0.1479+0.1502*n_o-0.0176*n_o*sigma                !eq.5b
                end if
                
                albedo = f_dir*alphasdir+f_dif*alphasdif+0.006                   !eq.15
                
                ! foam corrected alternative
                !fwc = 2.95e-6*wind**3.52                                          !eq.16
                !albedoe = 0.55*fwc+albedo*(1-fwc)    !foam corrected albedo  (Koepkw,1984)    !eq.17
                !       albedo = qdir_frac*alphasdir+qdiff_frac*alphasdif+0.006  !?
                
             end if
             
          else
             albedo = 0.
          end if
          !---------------------------------------- end Jin et al. (2011)
          
       case default
          stop 'NotImplementedError'
          
       end select
       ! --------------------------------- End albedo parametrization options ------------
          
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
          I_0_calc = qshort
          I_0_cs   = qtot*(1.-albedo)
       end if

       return
     end subroutine short_wave_radiation
     !EOC

     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Read meteo data, interpolate in time
     !
     ! !INTERFACE:
     subroutine flux_from_meteo(jul,secs)
       !
       ! !DESCRIPTION:
       !  This routine reads meteo data from a file and calculates the 
       !  fluxes of heat and momentum, and the
       !  short--wave radiation, from these data as described in 
       !  \sect{sec:calcCoeff}, \sect{sec:calcFluxes}, and \sect{sec:swr}.
       !  Then, the results are interpolated in time.

       !
       ! !USES:
       IMPLICIT NONE
       !
       ! !INPUT PARAMETERS:
       integer, intent(in)                 :: jul,secs
       !
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module
       !
!!!!! SH - 11/08/2003 - removed calculated short wave from this routine - not needed
       ! model calculated I_0 is better that linear interpolation
       ! Probably was here for completeness but interferes with full range of flux permutations 
       !

       !EOP
       !
       ! !LOCAL VARIABLES:
       integer                   :: yy,mm,dd,hh,min,ss
       double precision                  :: t
       double precision, SAVE            :: dt
       integer, save             :: meteo_jul1,meteo_secs1
       integer, save             :: meteo_jul2=0,meteo_secs2=0
       !double precision, save        :: obs(6)
       !WT added for SP's new code below:
       double precision, save            :: obs1(6), obs2(6)

       !HK: changed alpha(4) to alpha(6)
       double precision, save            :: alpha(9)

       !HK added :
       double precision         :: wx1,wx2,wy1,wy2
       !SP added :
       double precision,save         :: qb1,qb2=0.,qh1,qh2=0.,qe1,qe2=0.
       integer, save                 :: count,count2
       integer                       :: count3,ios
       double precision              :: sst_obs

       double precision, save            :: I1,h1,tx1,ty1
       double precision, save            :: I2=0.,h2=0.,tx2=0.,ty2=0.
       logical, save             :: first=.true.
       integer                   :: rc
       integer             ::loop_counter=0

       !


       !-----------------------------------------------------------------------
       !BOC

       !WT SP's new code as of 20171115, replaces the whole subroutine.
       !  This part initialises and reads in new values if necessary
       if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .lt. 0) then
          do
             meteo_jul1 = meteo_jul2
             meteo_secs1 = meteo_secs2
             obs1 = obs2
             call read_obs(meteo_unit,yy,mm,dd,hh,min,ss,6,obs2,rc)
             call julian_day(yy,mm,dd,meteo_jul2)
             meteo_secs2 = hh*3600 + min*60 + ss
             if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .gt. 0) EXIT
          end do
          dt = time_diff(meteo_jul2,meteo_secs2,meteo_jul1,meteo_secs1)
       end if
       !  Do the time interpolation
       t  = time_diff(jul,secs,meteo_jul1,meteo_secs1)

       alpha(1) = (obs2(1)-obs1(1))/dt
       wx    = obs1(1) + t*alpha(1)
       wx_obs = obs1(1) + t*alpha(1)
       alpha(2) = (obs2(2)-obs1(2))/dt
       wy    = obs1(2) + t*alpha(2)
       wy_obs = obs1(2) + t*alpha(2)
       alpha(3) = (obs2(3)-obs1(3))/dt
       airp  = (obs1(3) + t*alpha(3))*100. !kbk mbar/hPa --> Pa
       alpha(4) = (obs2(4)-obs1(4))/dt
       airt  = obs1(4) + t*alpha(4)
       alpha(5) = (obs2(5)-obs1(5))/dt
       spec_hum = obs1(5) + t*alpha(5)
       !      dew_pt = obs1(5) + t*alpha(5)
       !      rh    = obs1(5) + t*alpha(5)
       alpha(6) = (obs2(6)-obs1(6))/dt
       cloud = obs1(6) + t*alpha(6)

       call exchange_coefficients()

       call do_calc_fluxes(qb,qh,qe,taux=tx,tauy=ty)

     end subroutine flux_from_meteo
     !EOC

     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Read heat flux data, interpolate in time
     !
     ! !INTERFACE:
     subroutine read_heat_flux(jul,secs,I_0,heat,qb)
       !
       ! !DESCRIPTION:
       !  This routine reads heat fluxes from a file and interpolates them in
       !  time.
       !
       ! !USES:
       IMPLICIT NONE
       !
       ! !INPUT PARAMETERS:
       integer, intent(in)                 :: jul,secs
       !
       ! !OUTPUT PARAMETERS:
       double precision, intent(out)               :: I_0,heat,qb
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module

       !
       !EOP
       !
       ! !LOCAL VARIABLES:
       integer                   :: yy,mm,dd,hh,min,ss
       double precision                  :: t,alpha
       double precision, SAVE            :: dt
       integer, save             :: heat_jul1,heat_secs1
       integer, save             :: heat_jul2=0,heat_secs2=0
       double precision, save            :: obs1(3),obs2(3)=0.
       integer                   :: rc
       !
       !-----------------------------------------------------------------------
       !BOC
       !  This part initialise and read in new values if necessary.
       if(time_diff(heat_jul2,heat_secs2,jul,secs) .lt. 0) then
          do
             heat_jul1 = heat_jul2
             heat_secs1 = heat_secs2
             obs1 = obs2
             call read_obs(heat_unit,yy,mm,dd,hh,min,ss,3,obs2,rc)
             call julian_day(yy,mm,dd,heat_jul2)
             heat_secs2 = hh*3600 + min*60 + ss
             if(time_diff(heat_jul2,heat_secs2,jul,secs) .gt. 0) EXIT
          end do
          dt = time_diff(heat_jul2,heat_secs2,heat_jul1,heat_secs1)
       end if

       !  Do the time interpolation
       t  = time_diff(jul,secs,heat_jul1,heat_secs1)

       !**I'VE COMMENTED OUT INTERPOLATION**
       alpha = (obs2(1)-obs1(1))/dt
       I_0 = obs1(1) !+ t*alpha
       alpha = (obs2(2)-obs1(2))/dt
       heat = obs1(2) !+ t*alpha
       alpha = (obs2(3)-obs1(3))/dt
       qb = obs1(3) !+ t*alpha

     end subroutine read_heat_flux
     !EOC

     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Read momentum flux data, interpolate in time
     !
     ! !INTERFACE:
     subroutine read_momentum_flux(jul,secs,tx,ty)
       !
       ! !DESCRIPTION:
       !  This routine reads momentum fluxes from a file and interpolates them in
       !  time.
       !
       ! !USES:
       IMPLICIT NONE
       !
       ! !INPUT PARAMETERS:
       integer,intent(in)                  :: jul,secs
       !
       ! !OUTPUT PARAMETERS:
       double precision,intent(out)                :: tx,ty
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module
       !
       ! !LOCAL VARIABLES:
       integer                   :: yy,mm,dd,hh,min,ss
       double precision                  :: t,alpha
       double precision, save            :: dt
       integer, save             :: mom_jul1,mom_secs1
       integer, save             :: mom_jul2=0,mom_secs2=0
       double precision, save            :: obs1(2),obs2(2)=0.
       integer                   :: rc
       !
       !EOP
       !-----------------------------------------------------------------------
       !BOC

       !  This part initialise and read in new values if necessary.
       if(time_diff(mom_jul2,mom_secs2,jul,secs) .lt. 0) then
          do
             mom_jul1 = mom_jul2
             mom_secs1 = mom_secs2
             obs1 = obs2
             call read_obs(momentum_unit,yy,mm,dd,hh,min,ss,2,obs2,rc)
             call julian_day(yy,mm,dd,mom_jul2)
             mom_secs2 = hh*3600 + min*60 + ss
             if(time_diff(mom_jul2,mom_secs2,jul,secs) .gt. 0) EXIT
          end do
          dt = time_diff(mom_jul2,mom_secs2,mom_jul1,mom_secs1)
       end if

       !**I'VE COMMENTED OUT INTERPOLATION**
       !  Do the time interpolation
       t  = time_diff(jul,secs,mom_jul1,mom_secs1)
       alpha = (obs2(1)-obs1(1))/dt
       tx = obs1(1) !+ t*alpha
       alpha = (obs2(2)-obs1(2))/dt
       ty = obs1(2) !+ t*alpha

       return
     end subroutine read_momentum_flux
     !EOC

     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Read SST, interpolate in time
     !
     ! !INTERFACE:
     subroutine read_sst(jul,secs,sst)
       !
       ! !DESCRIPTION:
       !  This routine reads sea surface temperature (SST) from a file 
       !  and interpolates in time.
       !
       ! !USES:
       use meanflow, only                : T,S,h

       IMPLICIT NONE

       ! !INPUT PARAMETERS:
       integer, intent(in)                 :: jul,secs
       !
       ! !OUTPUT PARAMETERS:
       double precision,intent(out)                :: sst
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module
       !
       !EOP
       !
       ! !LOCAL VARIABLES:
       integer                   :: yy,mm,dd,hh,min,ss
       integer                   :: yy2,mm2,dd2,hh2,min2,ss2
       double precision                  :: time,alpha
       double precision, save            :: dt
       integer, save             :: sst_jul1,sst_secs1
       integer, save             :: sst_jul2=0,sst_secs2=0
       integer, save             :: sst2_jul1,sst2_secs1
       integer, save             :: sst2_jul2=0,sst2_secs2=0
       double precision, save            :: obs1(1),obs2(1)=0.,obs3(4),obs4(4)=0.
       logical, save             :: first=.true.,first2=.true.
       integer                   :: rc
       !SP
       double precision          :: delta
       !
       !-----------------------------------------------------------------------
       !BOC
       if(sst_method==2) then
          !use sst observation to assimilate into mixed layer
          if(time_diff(sst_jul2,sst_secs2,jul,secs) .lt. 0) then
             do
                sst_jul1 = sst_jul2
                sst_secs1 = sst_secs2
                obs1 = obs2
                call read_obs(sst_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
                call julian_day(yy,mm,dd,sst_jul2)
                sst_secs2 = hh*3600 + min*60 + ss
                if(time_diff(sst_jul2,sst_secs2,jul,secs) .gt. 0) EXIT
             end do
             dt = time_diff(sst_jul2,sst_secs2,sst_jul1,sst_secs1)

             if(first) then
                first=.false.
                !SP 21/03/06
                !Assimilate SST observations at observation time
                sst=obs1(1)     !sets observation to sst value
                ostia=obs1(1)
                !PRINT *,obs1(1),'OSTIA'
                ostia_diff=ostia_diff+(T(150)-obs1(1))
                ostia_sq_diff=ostia_sq_diff+(T(150)-obs1(1))**2
                call assimilate_satellite_obs(sst,T(1:150),S(1:150),h(1:150))
                !END SP
             else
             end if

             call write_time_string(jul,secs,assim_timestr)
             print *,'SST values assimilated at ', assim_timestr
             print *,'Value:',sst
          end if
       end if
       if(sst_method2==2) then 
          !use sst observation to calculate a cost function
          if(time_diff(sst2_jul2,sst2_secs2,jul,secs) .lt. 0) then
             do
                sst2_jul1 = sst2_jul2
                sst2_secs1 = sst2_secs2
                obs3 = obs4
                call read_obs(sst_unit2,yy2,mm2,dd2,hh2,min2,ss2,4,obs4,rc)
                call julian_day(yy2,mm2,dd2,sst2_jul2)
                sst2_secs2 = hh2*3600 + min2*60 + ss2
                if(time_diff(sst2_jul2,sst2_secs2,jul,secs) .gt. 0) EXIT
             end do
             dt = time_diff(sst2_jul2,sst2_secs2,sst2_jul1,sst2_secs1)

             !SP 28/03/06
             !use obs to calculate cost function
             if(first2) then
                seviri_diff=0.
                amsre_diff=0.
                tmi_diff=0.
                ostia_seviri_diff=0.
                ostia_amsre_diff=0.
                ostia_tmi_diff=0.
                seviri_sq_diff=0.
                amsre_sq_diff=0.
                tmi_sq_diff=0.
                ostia_seviri_sq_diff=0.
                ostia_amsre_sq_diff=0.
                ostia_tmi_sq_diff=0.
                seviri_obs=0.
                amsre_obs=0.
                tmi_obs=0.
                first2=.false.
             else
                if(ABS(obs3(4)-2.0).LT.0.0001) then
                   PRINT*,obs3(1),'SEVIRI OBS'
                   seviri_diff=seviri_diff+(skint-obs3(1))
                   ostia_seviri_diff=ostia_seviri_diff+(ostia-obs3(1))
                   seviri_sq_diff=seviri_sq_diff+(skint-obs3(1))**2
                   ostia_seviri_sq_diff=ostia_seviri_sq_diff+(ostia-obs3(1))**2
                   seviri_obs=seviri_obs+1
                else if(ABS(obs3(4)-1.0).LT.0.001) then
                   PRINT*,obs3(1),'AMSRE OBS'
                   amsre_diff=amsre_diff+(T(150)-obs3(1))
                   ostia_amsre_diff=ostia_amsre_diff+(ostia-obs3(1))
                   amsre_sq_diff=amsre_sq_diff+(T(150)-obs3(1))**2
                   ostia_amsre_sq_diff=ostia_amsre_sq_diff+(ostia-obs3(1))**2
                   amsre_obs=amsre_obs+1
                else
                   PRINT*,obs3(1),'TMI OBS'
                   tmi_diff=tmi_diff+(T(150)-obs3(1))
                   ostia_tmi_diff=ostia_tmi_diff+(ostia-obs3(1))
                   tmi_sq_diff=tmi_sq_diff+(T(150)-obs3(1))**2
                   ostia_tmi_sq_diff=ostia_tmi_sq_diff+(ostia-obs3(1))**2
                   tmi_obs=tmi_obs+1
                end if
             end if

          end if
       end if

       !  Do the time interpolation
       !   time  = time_diff(jul,secs,sst_jul1,sst_secs1)
       !   alpha = (obs2(1)-obs1(1))/dt
       !   sst = obs1(1) + time*alpha

       return
     end subroutine read_sst
     !EOC
     !-----------------------------------------------------------------------
     subroutine assimilate_satellite_obs(sst,T,S,h)

       use eqstate, only: unesco

       IMPLICIT NONE

       double precision, intent(in)    :: sst,S(1:150),h(1:150)
       double precision, intent(inout) :: T(1:150)

       double precision  :: delta,assim_depth,z,density,change
       integer           :: i,count,tendepth,ml_level

       !-----------------------------------------------------------------------
       !method used from june 2006
       !find mixed layer level using definition from Kara, 2000, JGR
       tendepth=83
       z=0.0
       do i=1,tendepth
          z=z+h(i)
       end do
       density=unesco(S(tendepth),T(tendepth),z/10.,.true.)
       change=unesco(S(tendepth),T(tendepth)+0.8,z/10.,.true.)-density

       do i=tendepth-1,1,-1
          z=z-h(i)
          if(unesco(S(i),T(i),z/10.,.true.).LE.(density+change)) then
             ml_level=i
             exit
          end if
          if(unesco(S(i),T(i),z/10.,.true.).GE.(density-change)) then
             ml_level=i
             exit
          end if
          if(i==1) then
             ml_level=i
          end if
       end do
       !adjust mixed layer
       PRINT*,'ASSIMILATE SST'

       delta=sst-T(150)
       do i=ml_level,150
          T(i)=T(i)+delta
       end do


       !old method 
       !Assimilate SST into mixed layer
       !         delta=sst-T(150)
       !         assim_depth=mld(T(1:150))
       !         if(sst.LT.T(nint(assim_depth))) then
       !            do i=1,150
       !               assim_depth=assim_depth-1
       !               if((T(nint(assim_depth)).LE.sst).OR.(assim_depth.LT.1.1)) then
       !                  exit
       !               end if
       !            end do
       !         end if
       !         PRINT*,'assimilate'
       !         DO i=nint(assim_depth),150
       !            T(i)=T(i)+delta
       !         END DO

     end subroutine assimilate_satellite_obs
     !-----------------------------------------------------------------------
     !-----------------------------------------------------------------------
     function mld(T)
       double precision, intent(in):: T(1:150)
       integer                     :: twenty_level, max_T_level
       integer                     :: i,count,mld

       max_T_level=150
       twenty_level=66

       !Find the grid level of the maximum temperature in the top 20m
       DO i=150-1,twenty_level,-1
          if (T(max_T_level).lt.T(i)) then
             max_T_level=i
          end if
       END DO
       !end fine the grid level of the max temp in top 20m

       !find the grid level of the mixed layer depth            
       Do count=max_T_level-1,1,-1
          if ((T(max_T_level)-T(count)).GE.0.1) then
             exit
          end if
       end do
       mld=count
       if(mld==0) then
          mld=1
       end if
       !end find the grid level of the mixed layer depth

     end function mld
     !-----------------------------------------------------------------------
     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Integrate short--wave and sea surface fluxes
     !
     ! !INTERFACE:
     subroutine integrated_fluxes(dt,int_cs)
       !
       ! !DESCRIPTION:
       !  This utility routine integrates the short--wave radiation 
       !  and heat--fluxes over time.
       !
       ! !USES:
       IMPLICIT NONE
       !
       ! !INPUT PARAMETERS:
       double precision, intent(in)                :: dt
       double precision, intent(inout)             :: int_cs
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module
       !
       !EOP
       !-----------------------------------------------------------------------
       !BOC
       int_sw = int_sw + dt*I_0
       int_cs = int_cs + dt*I_0_cs      !SP, integrated clear sky swr
       int_hf = int_hf + dt*heat
       int_total = int_sw + int_hf
       return
     end subroutine integrated_fluxes
     !EOC

     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Set the SST to be used from model.
     !
     ! !INTERFACE:
     subroutine set_sst(temp)
       !
       ! !DESCRIPTION:
       !  This routine sets the sea surface temperature (SST) to be used for 
       !  the surface flux calculations.
       !
       ! !USES:
       IMPLICIT NONE
       !
       ! !INPUT PARAMETERS:
       double precision, intent(in)                :: temp
       !
       ! !REVISION HISTORY:
       !  Original author(s): Karsten Bolding
       !
       !  See log for airsea module
       !
       !EOP
       !-----------------------------------------------------------------------
       !BOC
       sst = temp
       return
     end subroutine set_sst
     !EOC

     !-----------------------------------------------------------------------

   end module airsea

   !-----------------------------------------------------------------------
   ! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
   !----------------------------------------------------------------------- 
