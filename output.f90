!-----------------------------------------------------------------------
!BOP
!
! !MODULE: output -- saving the results
!
! !INTERFACE:
   module output
!
! !DESCRIPTION:
!  This module acts as an interface between GOTM and modules/routines
!  doing the actual output. In order to add a new output format it is only
!  necessary to add hooks in this module and write the actual output
!  routines. It is not necessary to change anything in GOTM itself.
!
! !USES:
   use time, ONLY: write_time_string,julianday,secondsofday,timestep
   use asciiout

   use ncdfout, ONLY:  init_ncdf,do_ncdf_out,close_ncdf


   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   logical                             :: write_results
   character(len=19)                   :: ts
   integer                             :: out_fmt=1
   character(len=255)             :: out_dir='.'
   character(len=255)             :: out_fn='gotm'
   character(len=255)             :: daily_stat_fn='daily_stat.dat', sst_event_fn='sst_event.dat'
   integer                             :: nsave=1
   logical                             :: variances=.false.
   logical                             :: diagnostics=.false.
   integer                             :: mld_method=1
   double precision                            :: diff_k=1.e-5
   double precision                            :: Ri_crit=0.5
   logical                             :: rad_corr=.true.

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!
!  $Log: output.F90,v $
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:24:19  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:53:05  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/11/18 11:51:52  gotm
!  Fixed a typo
!
!  Revision 1.1.1.1  2001/02/12 15:55:59  gotm
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
   integer, private, parameter         :: ascii_unit=50
   integer, private, parameter         :: grads_unit=51
!  Used for diagnostic output
   integer, private, parameter         :: temp_unit=54
   integer, private, parameter         :: mld_unit=55
   integer, private, parameter         :: sf_unit=56
   integer, private, parameter         :: fric_unit=57
   integer, private, parameter         :: heat_unit=58
   integer, private, parameter         :: energy_unit=59
   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the output module
!
! !INTERFACE:
   subroutine init_output(title,nlev,latitude,longitude)
!
! !DESCRIPTION:
!  Calls the initialization routine based on output format selected by 
!  the user.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*), intent(in)        :: title
   integer, intent(in)                 :: nlev
   double precision, intent(in)                :: latitude,longitude
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See output module
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=255)   :: ext,fname,tmp_fn
!
!-------------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'init_output'
   call write_time_string(julianday,secondsofday,ts)

   select case (out_fmt)
      case (1)
         EXT = 'out'
         fname = TRIM(out_dir) //'/'// TRIM(out_fn) // '.' // ext
         write(0,*) '       ', 'Output in ASCII:'
         write(0,*) '       ', TRIM(fname)
         call init_ascii(fname,title,ascii_unit)
      case (2)

         EXT = 'nc'
         fname = TRIM(out_dir) //'/'// TRIM(out_fn) // '.' // ext
         write(0,*) '       ', 'Output in NetCDF (time unit is set to seconds):'
         write(0,*) '       ', TRIM(fname)
         call init_ncdf(fname,title,latitude,longitude,nlev,ts,0)
      case (3)
         EXT = 'nc'
         fname = TRIM(out_dir) //'/'// TRIM(out_fn) // '.' // ext
         write(0,*) '       ', 'Output in GrADS specific NetCDF (time unit is hours):'
         call init_ncdf(fname,title,latitude,longitude,nlev,ts,2)

      case default
        write(0,*) '   ', 'Fatal error: A non valid output format has been chosen'
        stop 'init_output'
   end select

   if (diagnostics) then
      tmp_fn = TRIM(out_dir) //'/'// TRIM(out_fn) // '.temp'
      open(temp_unit,file=trim(tmp_fn),err=100,status='unknown')
      tmp_fn = TRIM(out_dir) //'/'// TRIM(out_fn) // '.mld'
      open(mld_unit,file=trim(tmp_fn),err=100,status='unknown')
      tmp_fn = TRIM(out_dir) //'/'// TRIM(out_fn) // '.sf'
      open(sf_unit,file=trim(tmp_fn),err=100,status='unknown')
      tmp_fn = TRIM(out_dir) //'/'// TRIM(out_fn) // '.fric'
      open(fric_unit,file=trim(tmp_fn),err=100,status='unknown')
      tmp_fn = TRIM(out_dir) //'/'// TRIM(out_fn) // '.heat'
      open(heat_unit,file=trim(tmp_fn),err=100,status='unknown')
      tmp_fn = TRIM(out_dir) //'/'// TRIM(out_fn) // '.energy'
      open(energy_unit,file=trim(tmp_fn),err=100,status='unknown')
   end if
   return
100 write(0,*) 'FATAL ERROR: ', 'unable to open ',trim(tmp_fn),' for output.'
   stop

   end subroutine init_output
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set some variables related to output
!
! !INTERFACE:
   subroutine prepare_output(n)
!
! !DESCRIPTION:
!  This routine check whether output should  be written at
!  the current time step. If this is the case, the model
!  time is written to a string for display on the screen.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See output module
!
!EOP
!-------------------------------------------------------------------------
!BOC
   write_results = mod(n,nsave).eq.0
   call write_time_string(julianday,secondsofday,ts)

   return
   end subroutine prepare_output
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save the model results in file
!
! !INTERFACE:
   subroutine do_output(n,nlev)
!
! !DESCRIPTION:
!  Calls the routine, which will do the actual storing of results, depending
!  on the output format.
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n,nlev
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: secs
!
!-------------------------------------------------------------------------
!BOC
   if (write_results) then

      if(variances) then
         call do_variances(nlev)
      end if

!kbk     call write_time_string(julianday,secondsofday,ts)
!kbk     call ss_nn_obs 
      write(0,137) '       ', 'Saving......',ts
137   format(TL38,A7,A12,A19,$) ! 38=7+12+19
      secs = nint(n*timestep)
      select case (out_fmt)
         case (1)
            call do_ascii_out(nlev,ts,variances,ascii_unit)

         case (2, 3)
            call do_ncdf_out(nlev,secs)

         case default
           write(0,*) '   ', 'Fatal error: A non valid output format has been chosen'
           stop 'do_output'
      end select
   end if

!  Diagnostic output
!   if(diagnostics) then
!      STDERR 'Diagnostics'
!      call do_diagnostics(n,nlev,BuoyMeth,dt,u_taus,u_taub)
!   end if

   return
   end subroutine do_output
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Close files used for saving model results
!
! !INTERFACE:
   subroutine close_output()
!
! !DESCRIPTION:
!  Call routines for closing any open output files.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See output module
!
!EOP
!-------------------------------------------------------------------------
!BOC
   select case (out_fmt)
      case (1)
         call close_ascii(ascii_unit)

      case (2, 3)
         call close_ncdf()

      case default
        write(0,*) '   ', 'Fatal error: A non valid output format has been chosen'
        stop 'do_output'
   end select

   if(diagnostics) then
      close(temp_unit); close(mld_unit); close(sf_unit);
      close(fric_unit); close(heat_unit); close(energy_unit)
   end if

   return
   end subroutine close_output
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate variances of turbulent quantities 
!
! !INTERFACE:
   subroutine do_variances(nlev)
!
! !DESCRIPTION:
!  This subroutine diagnostically computes variances of $u'$, 
!  $v'$, $w'$, $\theta'$ for output. Presently, the Algebraic
!  Stress Model of \cite{Canutoetal2001a} is used.
!
! !USES:
   use meanflow,   only: depth0,h,u,v,t,B
   use turbulence, only: tke,eps,num,nuh,uu,vv,ww,tt,chi
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard
!
!  See output module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   double precision                  :: dzu2,dzv2,dzt2
   double precision, parameter       :: L1=0.1070
   double precision, parameter       :: L2=0.0032
   double precision, parameter       :: L3=0.0864
   double precision, parameter       :: L4=0.1200
   double precision, parameter       :: L5=11.9000
   double precision, parameter       :: L6=0.4000
   double precision, parameter       :: L7=0.0000
   double precision, parameter       :: L8=0.4800
   double precision                  :: zz(0:nlev)
!
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      dzu2=( (u(i+1)-u(i)) / (0.5*(h(i+1)+h(i))) )**2
      dzv2=( (v(i+1)-v(i)) / (0.5*(h(i+1)+h(i))) )**2
      dzt2=( (t(i+1)-t(i)) / (0.5*(h(i+1)+h(i))) )**2
      uu(i)=2./3.*tke(i)+2./3.*tke(i)/eps(i)*(                      &
            num(i)*((L2+3*L3)*dzu2-2*L2*dzv2)-2*L4*B(i))
      vv(i)=2./3.*tke(i)+2./3.*tke(i)/eps(i)*(                      &
            num(i)*((L2+3*L3)*dzv2-2*L2*dzu2)-2*L4*B(i))
      ww(i)=2./3.*tke(i)+2./3.*tke(i)/eps(i)*(                      &
            num(i)*((L2-3*L3)*(dzu2+dzv2))+4*L4*B(i))
      tt(i)=1.44*nuh(i)*tke(i)/eps(i)*dzt2
      chi(i)=tt(i)*eps(i)/tke(i)*2./1.44
      if (abs(uu(i)+vv(i)+ww(i)-2*tke(i)) .gt. 1.e-8) then
         write(0,*) 'Velocity variances do not add up to 2k.'
         write(0,*) 'Program aborted.'
         write(0,*) i,uu(i),vv(i),ww(i)
         write(0,*) tke(i),abs(uu(i)+vv(i)+ww(i)-2*tke(i))
         stop
      end if
   end do

   return
   end subroutine do_variances
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute various diagnostic/integrated variables
!
! !INTERFACE:
   subroutine do_diagnostics(n,nlev,BuoyMeth,dt,u_taus,u_taub,I_0,heat)
!
! !DESCRIPTION:
!  This subroutine calculates the following diagnostic/integrated variables.
!
! !USES:
   use airsea,       only: sst
   use meanflow,     only: gravity,rho_0,cp
   use meanflow,     only: h,u,v,s,t,NN,SS,buoy,rad
   use turbulence,   only: kappa
   use turbulence,   only: tke
   use observations, only: tprof,b_obs_sbf
   use eqstate,      only: eqstate1
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n,nlev,BuoyMeth
   double precision, intent(in)                :: dt
   double precision, intent(in)                :: u_taus,u_taub
   double precision, intent(in)                :: I_0,heat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See output module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   double precision                  :: mld_surf,mld_bott
   double precision                  :: ekin,epot,eturb
   double precision, save            :: epot0
   double precision                  :: heat_sim,heat_obs
   double precision, save            :: heat_sim0,heat_obs0,heat_flux
   double precision                  :: z,dtt,dtb,x
   double precision                  :: wstar,tstar
   double precision                  :: sbf,stf,MOL
   double precision                  :: Ri(0:nlev)
   logical, save             :: first=.true.
!
!-----------------------------------------------------------------------
!BOC
   select case(mld_method)
      case(1)          ! MLD according to TKE criterium
         mld_surf    = 0.0
         i=nlev
100      i=i-1
         mld_surf=mld_surf+h(i+1)
         if ((tke(i) .gt. diff_k) .and. (i .gt. 0)) goto 100
         mld_bott    = 0.0
         i=0
101      i=i+1
         mld_bott=mld_bott+h(i)
         if ((tke(i) .gt. diff_k) .and. (i .lt. nlev)) goto 101
      case(2)          ! MLD according to critical Ri number
         do i=1,nlev-1
            Ri(i)=NN(i)/(SS(i)+1.e-10)
         end do
         mld_surf    = 0.0
         i=nlev
200      i=i-1
         mld_surf=mld_surf+h(i+1)
         if ((Ri(i) .lt. Ri_crit) .and. (i .gt. 0)) goto 200
         mld_bott    = 0.0
         i=0
201      i=i+1
         mld_bott=mld_bott+h(i)
         if ((Ri(i) .lt. Ri_crit) .and. (i .lt. nlev)) goto 201
      case default
   end select

!  Here, the surface buoyancy flux (sbf) and the surface temperature
!  flux (stf) are calculated.

   if (BuoyMeth .eq. 1) then
      dtt=0.01
      x=0.5*h(nlev)/10.
      dTb=(eqstate1(S(nlev),T(nlev)+0.5*dtt,x,gravity,rho_0)         &
          -eqstate1(S(nlev),T(nlev)-0.5*dtt,x,gravity,rho_0))/dtt
      sbf=-heat/cp/rho_0*dTb
!     Correction of surface buoyancy and temperature flux
!     for solar radiation penetration
      if (rad_corr) then
         z=0.
         i=nlev+1
444      i=i-1
         z=z+h(i)
         if (z.lt.mld_surf) goto 444
         sbf=sbf-dTb*(rad(nlev)-rad(i))
!Using dTb such as here is not fully
!correct, but we assume that the
!thermal expansion coefficient does not
!vary significantly over the mixed layer.
      end if
      stf=sbf/dTb
   else
      sbf=b_obs_sbf
      stf=0.
   end if
   if (sbf.ge.0) then
      wstar=(sbf*mld_surf)**(1./3.)
   else
      wstar=-(-sbf*mld_surf)**(1./3.)
   end if
   if (wstar.eq.0) then
      TStar=1.e15
   else
      tstar=stf/wstar
   end if
!  Calculation of Monin-Obukhov length MOL:
   if (abs(sbf).lt.1.e-10) then
      if (sbf.ge.0) MOL=-1.e10
      if (sbf.lt.0) MOL= 1.e10
   else
      MOL=-u_taus**3/kappa/sbf
   end if

   heat_sim=0
   heat_obs=0
   do i=1,nlev
      heat_sim=heat_sim+T(i)*h(i)*rho_0*cp
      heat_obs=heat_obs+tprof(i)*h(i)*rho_0*cp
   end do
   if (first) then
      heat_sim0=heat_sim
      heat_obs0=heat_obs
      heat_flux=0.
   end if
   heat_sim=heat_sim-heat_sim0
   heat_obs=heat_obs-heat_obs0

   heat_flux=heat_flux+dt*(I_0+heat)

   ekin=0.
   epot=0.
   eturb=0.
   z=0.
   do i=1,nlev
      z=z-0.5*h(i)
      ekin=ekin+0.5*h(i)*(u(i)**2+v(i)**2)
      eturb=eturb+h(i)*(tke(i)+tke(i-1))
      epot=epot+h(i)*buoy(i)*z
      z=z-0.5*h(i)
   end do
   if (first) then
      epot0=epot
   end if
   epot=epot-epot0
   ekin=ekin*rho_0
   epot=epot*rho_0
   eturb=eturb*rho_0

!  The output parameters are:
!  mld_surf:  Surface mixed layer depth
!  mld_bott:  Bottom mixed layer depth
!  sbf    :  Surface buoyancy flux
!  stf    :  Surface temperature flux
!  MOL    :  Monin-Obukhov length
!  wstar  :  Deardorff convective velocity scale
!  tstar  :  Deardorff convective temperature scale
!  u_taus :  Surface friction velocity
!  u_taub :  Bottom friction velocity
!  heat_flux: Accumulated surface heat flux J/m^2
!  heat_sim: Relative heat content from simulated T-profiles J/m^2
!  heat_obs: Relative heat content from observed T-profiles J/m^2
!  ekin   : Kinetic energy of the water column J/m^2
!  epot   : Potential energy of the water column J/m^2
!  eturb  : Turbulent energy of the water column J/m^2

   x = N*dt/(86400.)
   write(temp_unit,111)   x,sst,t(nlev)
   write(mld_unit,111)    x,mld_surf,mld_bott
   write(sf_unit,111)     x,sbf,stf,MOL
   write(fric_unit,111)   x,wstar,tstar,u_taus,u_taub
   write(heat_unit,111)   x,heat_flux,heat_sim,heat_obs
   write(energy_unit,111) x,ekin,epot,eturb

   first=.false.
111 format(F10.5,1x,4(E12.5,1x))

   return
   end subroutine do_diagnostics
!EOC

!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
