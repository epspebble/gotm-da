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
   integer, parameter                  :: sst_unit2=27
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
   integer, public                   :: sst_method2
   integer, public                   :: sss_method
   integer, public                   :: airt_method


!!HK - make meteo_file public)
   character(len=255), public   :: meteo_file
   character(len=255), public   :: heatflux_file
   character(len=255), public   :: momentumflux_file
   character(len=255), public   :: p_e_flux_file
   character(len=255), public   :: sss_file
   character(len=255), public   :: sst_file
   character(len=255), public   :: sst_file2
   character(len=255), public   :: airt_file

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
   double precision, public                  :: border

   double precision                  :: es,ea,e,qa,L,qs,q,mr,xlv,rnl
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
99 write(0,*) 'FATAL ERROR: ', 'I could not open ',trim(sst_file2)
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
!               I_0=I_0_calc*adjustment
!               cloud=1-adjustment  !cloud index
               I_0=adjustment*(1.-albedo)
!               cloud=1-I_0/I_0_calc
            end if 
         case default
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
         call skin_temp(sst,skint)
      end if


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
!  Based on the wind vector at 10 m, the sea surface temperature (either 
!  from the model or from a data file), the 
!  relative humidity, the air temperature and the air pressure at 2 m,
!  this function computes the surface momentum flux and the latent and 
!  sensible heat flux coefficients according to the \cite{Kondo75} 
!  bulk formulae.
!
! !USES:
   IMPLICIT NONE
!
! !DEFINED PARAMETERS:
   double precision, parameter                 :: a1=6.107799961
   double precision, parameter                 :: a2=4.436518521e-1
   double precision, parameter                 :: a3=1.428945805e-2
   double precision, parameter                 :: a4=2.650648471e-4
   double precision, parameter                 :: a5=3.031240396e-6
   double precision, parameter                 :: a6=2.034080948e-8
   double precision, parameter                 :: a7=6.136820929e-11
   double precision, parameter                 :: eps=1.0e-12
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the airsear module
!
!EOP
!
! !LOCAL VARIABLES:
   double precision                  :: tvirt,s,s0
   double precision                  :: ae_d,be_d,pe_d
   double precision                  :: ae_h,be_h,ce_h,pe_h
   double precision                  :: ae_e,be_e,ce_e,pe_e
   double precision                  :: x
   double precision                  :: w_ten,cdd_ten,chd_ten,ced_ten,z_zero
   integer                           ::iterate
!
!
!-----------------------------------------------------------------------
!BOC
   w = sqrt(wx*wx+wy*wy)
L=2.5e6
!saturated vapour pressure
es=EXP( -L/(461.7*(sst+273.15)) + LOG(611.0) + L/(273.15*461.7))

!vapour pressure
ea=0.01*rh*es

!specific humidity
qa=(const06*ea/(airp-ea))/(1.+(const06*ea/(airp-ea)))

!specific humidity saturated at the sst
qs=(const06*es/(airp-es))/(1.+(const06*es/(airp-es)))           

!air density
rho_air=airp/(287.05*(airt+273.15)*(1.+.61*qa))

!virtual temperature
tvirt=(1.+0.61*qa)*(airt+273.15)

!  Transfer coefficient for heat and momentum

   if (w .lt. 2.2) then
      ae_d=0.0;   be_d=1.08;                  pe_d=-0.15;
      ae_h=0.0;   be_h=1.185;  ce_h=0.0;      pe_h=-0.157;
      ae_e=0.0;   be_e=1.23;   ce_e=0.0;      pe_e=-0.16;
   else if (w .lt. 5.0) then
      ae_d=0.771; be_d=0.0858;                pe_d=1.0;
      ae_h=0.927; be_h=0.0546; ce_h=0.0;      pe_h=1.0;
      ae_e=0.969; be_e=0.0521; ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 8.0) then
      ae_d=0.867; be_d=0.0667;                pe_d=1.0;
      ae_h=1.15;  be_h=0.01;   ce_h=0.0;      pe_h=1.0;
      ae_e=1.18;  be_e=0.01;   ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 25.0) then
      ae_d=1.2;   be_d=0.025;                 pe_d=1.0;
      ae_h=1.17;  be_h=0.0075; ce_h=-0.00045; pe_h=1.0;
      ae_e=1.196; be_e=0.008;  ce_e=-0.0004;  pe_e=1.0
   else
      ae_d=0.0;   be_d=0.073;                 pe_d=1.0;
      ae_h=1.652; be_h=-0.017; ce_h=0.0;      pe_h=1.0;
      ae_e=1.68;  be_e=-0.016; ce_e=0;        pe_e=1.0;
   end if

!transform C(10m) into C(z)
w_ten=w

if(abs(w_ten).lt.0.0001) then
      cdd=0.0
      chd=0.0
      ced=0.0
else

DO iterate=1,3
      cdd_ten=(ae_d+be_d*w_ten**pe_d)*1.0e-3
      chd_ten=(ae_h+be_h*w_ten**pe_h+ce_h*(w_ten-8.0)**2)*1.0e-3
      ced_ten=(ae_e+be_e*w_ten**pe_e+ce_e*(w_ten-8.0)**2)*1.0e-3 

      z_zero=EXP(LOG(10.)-.4/SQRT(cdd_ten))
      w_ten=w*LOG(10./z_zero)/LOG(3.63/z_zero)
!READ*
END DO 

cdd=(.4**2)*(.4*(cdd_ten**(-.5))-LOG(10./3.63))**(-2)
chd=.4*SQRT(cdd)/(.4*SQRT(cdd_ten)/chd_ten +LOG(3.63/10.))
ced=.4*SQRT(cdd)/(.4*SQRT(cdd_ten)/ced_ten +LOG(3.63/10.))

!  Stability
s0=((sst-airt)+0.61*(qs-qa)*(airt+0.01*3.63))/(w*1.+LOG10(10/3.63))**2
s=s0*abs(s0)/(abs(s0)+0.01)

!PRINT*,sst-airt
!PRINT*,(sst-airt)+0.61*(qs-qa)*(airt+0.01*3.63)

   if((sst-airt)+0.61*(qs-qa)*(airt+0.01*3.63) .lt. 0.) then
      if (s .gt. -3.3) then
         x = 0.1+0.03*s+0.9*exp(4.8*s)
      else
         x = 0.0
      end if
      cdd=x*cdd
      chd=x*chd
      ced=x*ced
   else
      cdd=cdd*(1.0+0.47*sqrt(s))
      chd=chd*(1.0+0.63*sqrt(s))
      ced=ced*(1.0+0.63*sqrt(s))
   end if
!PRINT*,cdd
!PRINT*,chd
!PRINT*,ced
!PRINT*,integer2
end if

   return
   end subroutine exchange_coefficients
!EOC
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
       qe=ced*L*rho_air*w*(qs-qa)            ! latent - original
   qh=chd*cpa*rho_air*w*(sst-airt)       ! sensible - original

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

   tmp = -cdd*rho_air*w


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

!SP-this represents the possible errors in the derived cloud value
!   if((cloud.LT.0.1).AND.(ABS(swr_error+0.1).LT.0.0001)) then
!       swr_error=0.0
!   else if((cloud.LT.0.1).AND.(ABS(swr_error-0.1).LT.0.0001)) then
!       swr_error=0.2
!   else if ((cloud.GT.0.9).AND.(ABS(swr_error-0.1).LT.0.0001)) then
!       swr_error=0.0
!   else if((cloud.GT.0.9).AND.(ABS(swr_error+0.1).LT.0.0001)) then
!       swr_error=-0.2
!   end if

         if(swr_error.ne.0.0) then
            if(ABS(cloud-0.001).LT.0.0001) then  !SP IR obs therefore no cloud
               swr_error=0.0
               cloud=0.0
               border=1
            else if(ABS(swr_error-9.0).LT.0.0001) then
               swr_error=0.0
               if(ABS(border).LT.0.0001) then
                  border=2
               end if
            else 
               cloud=cloud+swr_error
               if (cloud.GT.1.0) then
                  cloud=1.0
                  border=3
               else if (cloud.LT.0.0) then
                  cloud=0.0
                  border=3
               else
                  if((ABS(border).LT.0.0001).OR.(ABS(border-2.0).LT.0.0001)) then
                     border=1
                  end if
               end if
            end if
         end if
!SP end swr_error section

   if(cloud .lt.0.3) then
      qshort  = qtot*(1.-albedo)        !SP albedo factor needed here
   else
      qshort  = qtot*(1.-.62*cloud + .0019*sunbet)*(1.-albedo) 
   end if

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
      wx_obs=obs(1)
      wy_obs=obs(2)
      airp  = obs(3)*100. !kbk mbar/hPa --> Pa
      airt  = obs(4)
!      spec_hum = obs(5)
!      dew_pt = obs(5)
      rh    = obs(5)
      cloud = obs(6)

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
   use meanflow, only                : T,S,h
   use meanflow, only                : diff,abs_diff,sq_diff,number_obs

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
      
      !SP 21/03/06
      !Assimilate SST observations at observation time
           sst=obs1(1)     !sets observation to sst value
           PRINT*,obs1(1),'OSTIA'
           call assimilate_satellite_obs(sst,T(1:150),S(1:150),h(1:150))
      !END SP

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
         diff=0.
         abs_diff=0.
         sq_diff=0.
         number_obs=0.
         first2=.false.
      else
         if(ABS(obs3(4)-2.0).LT.0.0001) then
            PRINT*,obs3(1),'IR OBS'
            diff=diff+(skint-obs3(1))
            abs_diff=abs_diff+ABS(T(150)-obs3(1))
            sq_diff=sq_diff+(T(150)-obs3(1))**2
            number_obs=number_obs+1
         else
            PRINT*,obs3(1),'MW OBS'
            diff=diff+(T(150)-obs3(1))
            abs_diff=abs_diff+ABS(T(150)-obs3(1))
            sq_diff=sq_diff+(T(150)-obs3(1))**2
            PRINT*,(T(150)-obs3(1))
            number_obs=number_obs+1
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
!ajust mixed layer
PRINT*,'ASSIMILATE'

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
