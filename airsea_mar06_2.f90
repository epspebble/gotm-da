!SP
!This module has the exchange_coefficients subroutine changed to calculate
!fluxes according to the COARE Bulk Flux Algorithm version 3.0a modified 
!slightly by Simon Josey and then adapted by Sam Pimentel in order to fit into
!GOTM.  
!In addition to airsea_solar1.f90 this has the full albedo table from Payne, 
!and I feel this is slightly neater.  Also if the Arabian sea observations for
!incoming solar radiation are used, then the albedo factor needs to be included
!to give a net solar radiation, so this is rectified.
!
!this includes improvements from airsea_solar3.f90, I have added the ability to!read in the incoming longwave radiation measurements available for the arabian!sea data.
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
   use time, only: julian_day, time_diff, calendar_date
   use observations, only: read_obs
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

!
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
 

!
! !DEFINED PARAMETERS:
   integer, parameter                  :: meteo_unit=20
   integer, parameter                  :: heat_unit=21
   integer, parameter                  :: momentum_unit=22
   integer, parameter                  :: p_e_unit=23
   integer, parameter                  :: sst_unit=24
   integer, parameter                  :: sss_unit=25
   integer, parameter                  :: airt_unit=26

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
   integer, public                   :: sss_method
   integer, public                   :: airt_method


!!HK - make meteo_file public)
   character(len=255), public   :: meteo_file
   character(len=255), public   :: heatflux_file
   character(len=255), public   :: momentumflux_file
   character(len=255), public   :: p_e_flux_file
   character(len=255), public   :: sss_file
   character(len=255), public   :: sst_file
   character(len=255), public   :: airt_file

   double precision                  :: wx,wy
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
   double precision, public                  :: rho_air
   double precision, public                  :: const_tx,const_ty
   double precision, public                  :: const_qin,const_qout
   double precision, public                  :: clark_lamda,swr_error,wind_error  !SP
   double precision, public                  :: wind_h,rh_h,airt_h   !SP

   double precision                  :: es,ea,e,qs,q,mr,xlv,rnl
   double precision                  :: cdd,chd,ced,Du,zt,dqer

   double precision                  :: alon,alat

!!!!! SH 14/08/2003
   double precision, public	:: airt
   double precision, public	:: solar
   double precision, public	:: net_ir=1353 !1350.

!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the air--sea interaction module
!
! !INTERFACE:
   subroutine init_air_sea(namlst,lat,lon)
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
                     clark_lamda,                 &     !SP
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
                     sss_method, sss_file, &
                     airt_method, airt_file
!
!-----------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'init_air_sea'

   open(namlst,file='airsea.inp',action='read',status='old',err=90)
   read(namlst,nml=airsea,err=91)
   close(namlst)

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
            write(0,*) '       ', 'Reading precipitatio/evaporation data from:'
            write(0,*) '           ', trim(p_e_flux_file)
         case default
      end select

!     The sea surface temperature
      select case (sst_method)
         case (FROMFILE)
            open(sst_unit,file=sst_file,action='read',status='old',err=96)
            write(0,*) '       ', 'Reading sea surface temperature from:'
            write(0,*) '           ', trim(sst_file)
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
!   cloud=0.
   sss=0.
   airt=0.

   alon = deg2rad*lon
   alat = deg2rad*lat

   return

90 write(0,*) 'FATAL ERROR: ', 'I could not open airsea.inp'
   stop 'init_airsea'
91 write(0,*) 'FATAL ERROR: ', 'I could not read airsea namelist'
   stop 'init_airsea'
92 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(meteo_file)
   stop 'init_airsea'
93 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(heatflux_file)
   stop 'init_airsea'
94 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(momentumflux_file)
   stop 'init_airsea'
95 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(p_e_flux_file)
   stop 'init_airsea'
96 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(sst_file)
   stop 'init_airsea'
97 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(sss_file)
   stop 'init_airsea'
98 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(airt_file)
   stop 'init_airsea'

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
  double precision  :: dummy,adjustment,qb_down,lwrcloud,top,bottom
  integer           :: i,ios,count=1,count2=1,count3=0
!
!
!EOP
!-----------------------------------------------------------------------
!BOC

  if (calc_fluxes) call flux_from_meteo(jul,secs)     

   call short_wave_radiation(jul,secs,alon,alat)  !SP placed comment here

   !Net shortwave radiation
   I_0=I_0_calc
!     The heat fluxes

      select case (flux_method)
         case (CONSTVAL)
            heat=const_qout
         case (FROMFILE)
             call read_heat_flux(jul,secs,adjustment,heat,qb_down)
         case default
      end select

      select case (swr_method)
         case (CONSTVAL)
            I_0=const_qin
         case (FROMFILE)
            if (flux_method .ne. FROMFILE ) then
               PRINT*,'ERROR'
               READ*
            else
               I_0=I_0_calc*adjustment
               cloud=1-adjustment  !cloud index
!               I_0=adjustment*(1.-albedo)
!               cloud=1-I_0/I_0_calc
            end if 
         case default
      end select

!SP-this represents the possible errors in obs
              I_0=(1.+swr_error)*I_0
!SP end swr_error section

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

!SP-this represents the possible errors in obs
              tx=(1.+wind_error)*tx
              ty=(1.+wind_error)*ty
!SP end wind_error section

!     The sea surface temperature
      select case (sst_method)
         case (FROMFILE)
            call read_sst(jul,secs,sst)
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

!   Calculate cool skin effect (Wick, 96) SP: 13/03/06
      call skin_temp(sst,skint)

!   end if

   return
   end subroutine air_sea_interaction
!EOC
!-----------------------------------------------------------------------
   subroutine skin_temp(sst,skint)

! !USES:

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
   double precision                  ::hsb,hlb,qout,dels,qcol,alq,xlamx,tkt
   integer                           ::iter

   integer                           ::i,ios,count=1,count2=1,count3=0

!
!
!-----------------------------------------------------------------------
!BOC
   w = sqrt(wx*wx+wy*wy)

!SP - A file to read in wind errors so that wind obs can be adjusted 29/07/05              

!SP-this represents the possible errors in obs
!               w=(1.+wind_error)*w
!SP end wind errors section

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
      e=ea*rh*0.01                     !relative humidity = vapour pressure/satuated vapour pressure
       q=.62197*(e/(airp-0.378*e))     !mixing ratio   (kg/kg)

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
      endif
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
      do 10 iter=1,30
         
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
       endif
         Du=sqrt(w**2.+Wg**2.)        !include gustiness in wind spd.
!      
!use net swr and net lwr calculated from previous step
      rns=I_0
      rnl=qb   
!   Cool skin

           hsb=-rho_air*cpa*usr*tsr
           hlb=-rho_air*xlv*usr*qsr
           qout=rnl+hsb+hlb
           dels=rns*(.065+11.*tkt-6.6e-5/tkt*(1.-dexp(-tkt/8.0e-4))) !Eq.16 Ohlmann 
           qcol=qout-dels
         alq=Al*qcol+be*hlb*cpw/xlv                      !Eq. 7 Buoy flux water
         if(alq.gt.0.) then                              !originally (qcol.gt.0)
           xlamx=6./(1.+(bigc*alq/usr**4)**.75)**.333      !Eq 13 Saunders coeff.
           tkt=xlamx*visw/(sqrt(rho_air/rhow)*usr)          !Eq.11 Sublayer thickness
         else
           xlamx=6.                                      !prevent excessive warm skins
           tkt=min(.01,xlamx*visw/(sqrt(rho_air/rhow)*usr)) !Limit tkt
         endif
       dter=qcol*tkt/tcw                                 ! Eq.12 Cool skin
       dqer=wetc*dter

!(sxj) check for convergence and leave loop if met
       IF((iter==30).AND.(max(abs(conusr),abs(contsr),abs(conqsr)).gt.0.001)) THEN
            PRINT*,'convergence error'
         END IF
       if (max(abs(conusr),abs(contsr),abs(conqsr)).lt.0.001) then
          goto 912
       endif
   10 continue                                           ! end iterations

!(sxj) jump out point
912   continue

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
      endif
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
      endif
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
         qb=(1.0-clark_lamda*cloud*cloud)                                     &
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

   tmp=-rho_air*cdd*Du                            !stress N/m2


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
   double precision                  :: tau=0.7
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
   double precision, save            :: obs(6)

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

!  This part initialises and reads in new values if necessary
   if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .lt. 0) then
      do
         meteo_jul1 = meteo_jul2
         meteo_secs1 = meteo_secs2
         call read_obs(meteo_unit,yy,mm,dd,hh,min,ss,6,obs,rc)
         call julian_day(yy,mm,dd,meteo_jul2)
         meteo_secs2 = hh*3600 + min*60 + ss
         if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .gt. 0) EXIT
      end do
   end if           !SP I added this so that the fluxes are calculated each time step with a new sst, if using must comment out end if below.
      wx    = obs(1)
      wy    = obs(2)

!SP- this represents the possible 10% error in wind obs
!wx=1.1*wx
!wy=1.1*wy

      airp  = obs(3)*100. !kbk mbar/hPa --> Pa
      airt  = obs(4)
      rh    = obs(5)
!      cloud = obs(6)

      call exchange_coefficients()
!SP FEB05: cloud data assimilation scheme using skin sst obs
!count=count+1
!count2=count2+1 
!if(count==96) then
!   OPEN(UNIT=9,FILE="OBS/SUB/sst.dat",IOSTAT=ios)
!		IF(ios/=0) THEN
!			PRINT*,"Error during opening input file, stopping"; STOP
!		END IF
!	DO count3=1,count2+85               	
!		READ (UNIT=9,FMT='(3X,E13.12)')  sst_obs
!	END DO
!   CLOSE (UNIT=9)
   
!   cloud=cloud+(skint-sst_obs)/0.1
!   PRINT*,cloud
!   count=0
!   if (cloud.lt.0.0) then
!            cloud=0.0
!   else if (cloud.gt.1.0) then
!            cloud=1.0
!   end if
!end if     

!HK add in code to get time interpolated u10 and v10 out
!so they can be used in co2transfer.f90

      if (first) then
!         call do_calc_fluxes(heatf=h1,taux=tx1,tauy=ty1)
         call do_calc_fluxes(qb1,qh1,qe1,taux=tx1,tauy=ty1)
!         call short_wave_radiation(jul,secs,alon,alat,swr=I1)
!         I2  = I1                
!         h2  = h1
         tx2 = tx1
         ty2 = ty1
         qb2 = qb1
         qh2 = qh1
         qe2 = qe1
!HK added:
       wx1=wx
       wy1=wy        
       wx2 = wx1
       wy2 = wy1
!end of added

         first = .false.
      else
!         I1  = I2               
!         h1  = h2
         tx1 = tx2
         ty1 = ty2
         qb1 = qb2
         qh1 = qh2
         qe1 = qe2
!         call do_calc_fluxes(heatf=h2,taux=tx2,tauy=ty2)
         call do_calc_fluxes(qb2,qh2,qe2,taux=tx2,tauy=ty2)
	 
!         call short_wave_radiation(jul,secs,alon,alat,swr=I2) 
!HK added:
       wx1=wx2
       wy1=wy2
       wx2=wx
       wy2=wy
!end of added

      end if

      dt = time_diff(meteo_jul2,meteo_secs2,meteo_jul1,meteo_secs1)

!      alpha(1) = (I2-I1)/dt
!      alpha(2) = (h2-h1)/dt
      alpha(3) = (tx2-tx1)/dt
      alpha(4) = (ty2-ty1)/dt
      alpha(7) = (qb2-qb1)/dt
      alpha(8) = (qh2-qh1)/dt
      alpha(9) = (qe2-qe1)/dt
!HK added:
      alpha(5) =  (wx2-wx1)/dt
      alpha(6) =  (wy2-wy1)/dt
!end of added

!      end if           ! SP - commented out as it placed earlier

!  Do the time interpolation
   t  = time_diff(jul,secs,meteo_jul1,meteo_secs1)
!   I_0  = I1  + t*alpha(1)    
!   heat = h1  + t*alpha(2)
   tx   = tx1 + t*alpha(3)
   ty   = ty1 + t*alpha(4)
   qb   = qb1 + t*alpha(7)
   qh   = qh1 + t*alpha(8)
   qe   = qe1 + t*alpha(9)
   
!HK added:
   u10   = wx1 + t*alpha(5)
   v10   = wy1 + t*alpha(6)
!end of added

   return
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

   return
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

   IMPLICIT NONE
!
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
   double precision                  :: time,alpha
   double precision, save            :: dt
   integer, save             :: sst_jul1,sst_secs1
   integer, save             :: sst_jul2=0,sst_secs2=0
   double precision, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialise and read in new values if necessary.
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
   end if

!  Do the time interpolation
   time  = time_diff(jul,secs,sst_jul1,sst_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   sst = obs1(1) + time*alpha

   return
   end subroutine read_sst
!EOC

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
