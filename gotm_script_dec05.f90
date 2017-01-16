!SP
!25/11/04 : cloud=cloud+(sst(model)-sst(obs))/.21 use this formulae once a day
!02/05 : use sst(obs) to change mixed layer depth by adjusting profiles.
!05/05 : subroutine called assimilation coded up which include 4 options
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gotm --- the general framework \label{sec:gotm}
!! !INTERFACE:
   module gotm
!
! !DESCRIPTION:
! This is 'where it all happens'. This module provides the internal
! routines {\tt init\_gotm()} to initialise the whole model and 
! {\tt time\_loop()} to manage the time-stepping of all fields. These 
! two routines in turn call more specialised routines e.g.\ of the 
! {\tt meanflow} and {\tt turbulence} modules to delegate the job.
!
!  Here is also the place for a few words on Fortran `units' we used. 
!  The method of Fotran units is quite rigid and also a bit dangerous,
!  but lacking a better alternative we adopted it here. This requires
!  the definition of ranges of units for different purposes. In GOTM
!  we strongly suggest to use units according to the following
!  conventions.
!  \begin{itemize}
!     \item unit=10 is reserved for reading namelists.
!     \item units 20-29 are reserved for the {\tt airsea} module.
!     \item units 30-39 are reserved for the {\tt meanflow} module.
!     \item units 40-49 are reserved for the {\tt turbulence} module.
!     \item units 50-59 are reserved for the {\tt output} module.
!     \item units 60-69 are reserved for the {\tt extra} modules 
!           like those dealing with sediments or sea--grass.
!     \item units 70- are \emph{not} reserved and can be used as you 
!           wish.
!  \end{itemize}
!
! !USES:
   use airsea, only: init_air_sea,air_sea_interaction
   use airsea, only: set_sst,integrated_fluxes
   use airsea, only: calc_fluxes,swr_error,wind_error
   use airsea, only: tx,ty,I_0,heat,qh,qb,qe,cloud,int_cs
!HK 
   use airsea, only: qdir_frac, qdiff_frac,cosr

   use meanflow
   use turbulence
   use observations
   use output
   use time
   use mtridiagonal
   use eqstate

   use sediment


   use seagrass

!
   IMPLICIT NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_gotm, time_loop, clean_up
!
! !DEFINED PARAMETERS:
   integer, parameter                  :: namlst=10

   integer, parameter                  :: unit_sediment=61


   integer, parameter                  :: unit_seagrass=62
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: gotm.F90,v $
!  Revision 1.7  2003/03/28 09:20:34  kbk
!  added new copyright to files
!
!  Revision 1.6  2003/03/28 09:11:30  kbk
!  removed tabs
!
!  Revision 1.5  2003/03/10 09:20:27  gotm
!  Added new Generic Turbulence Model + improved documentation and cleaned up code
!
!  Revision 1.3  2001/11/18 15:58:02  gotm
!  Vertical grid can now be read from file
!
!  Revision 1.2  2001/06/13 07:40:39  gotm
!  Lon, lat was hardcoded in meteo.F90 - now passed via init_meteo()
!
!  Revision 1.1.1.1  2001/02/12 15:55:59  gotm
!  initial import into CVS
!
!EOP
!
!  private data members initialised via namelists
   character(len=80)         :: title
   integer                   :: nlev
   double precision                  :: dt
   double precision                  :: cnpar
   integer                   :: buoy_method
!  station description
   character(len=80)         :: name
   double precision                  :: latitude,longitude
   double precision            :: cloud_gradient
   integer                     :: sst_obs,profile_obs,obs_level,assimilation_type,assim_window
   double precision, public    :: advect(1:150)
!SP 19/09/05 test   
   double precision            :: j_one,j_one_b,j_two,j_three,j_four,j_five
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model \label{initGOTM}
!
! !INTERFACE:
   subroutine init_gotm()
!
! !DESCRIPTION:
!  This internal routine triggers the initialization of the model.
!  The first section reads the namelists of {\tt gotmrun.inp} with
!  the user specifications. Then, one by one each of the modules are
!  initialised with help of more specialised routines like 
!  {\tt init\_meanflow()} or {\tt init\_turbulence()} defined inside 
!  their modules, respectively.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!
! !LOCAL VARIABLES:
   namelist /model_setup/ title,nlev,dt,cnpar,buoy_method
   namelist /station/     name,latitude,longitude,depth
   namelist /time/        timefmt,MaxN,start,stop
   namelist /output/      out_fmt,out_dir,out_fn,nsave,variances, &
                          diagnostics,mld_method,diff_k,Ri_crit,rad_corr, &
                          assimilation_type,cloud_gradient,sst_obs,profile_obs,obs_level,assim_window    !SP
!   integer                ::count,ios
!
!-----------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'init_gotm'
   write(0,*) "------------------------------------------------------------------------"

!  open the namelist file.
   write(0,*) '       ', 'reading model setup namelists..'
   open(namlst,file='gotmrun.inp',status='old',action='read',err=90)
   read(namlst,nml=model_setup,err=91)
   timestep = dt ! timestep comes from the time module.
   read(namlst,nml=station,err=92)
   read(namlst,nml=time,err=93)
   read(namlst,nml=output,err=94)
   depth0=depth
   write(0,*) '       ', 'done.'

   write(0,*) '       ', trim(title)
   write(0,*) '       ', 'Using ',nlev,' layers to resolve a depth of',depth
   write(0,*) '       ', 'The station ',trim(name),' is situated at (lat,long) ', &
           latitude,longitude
   write(0,*) '       ', trim(name)

   write(0,*) '       ', 'initializing modules....'
   call init_time(MinN,MaxN)
   call init_eqstate(namlst)
   close (namlst)

!  From here - each init_? is responsible for opening and closing the
!  namlst - unit.
   call init_meanflow(namlst,'gotmmean.inp',nlev,latitude)
   call init_tridiagonal(nlev) 
   call updategrid(nlev,dt,zeta)
   call init_turbulence(namlst,'gotmturb.inp',nlev)
   call init_observations(namlst,'obs.inp',julianday,secondsofday, &
                          depth,nlev,z,h)
   
   s = sprof
   t = tprof
   u = uprof
   v = vprof
   call init_output(title,nlev,latitude,longitude)
   call init_air_sea(namlst,latitude,longitude)

!  Initialise each of the extra features/modules

   call init_sediment(namlst,'sediment.inp',unit_sediment,nlev, &
                      gravity,rho_0)


   call init_seagrass(namlst,'seagrass.inp',unit_seagrass,nlev,h)

   write(0,*) '       ', 'done.'
   write(0,*) "------------------------------------------------------------------------"

   return

90 write(0,*) 'FATAL ERROR: ', 'I could not open gotmrun.inp for reading'
   stop 'init_gotm'
91 write(0,*) 'FATAL ERROR: ', 'I could not read the "model_setup" namelist'
   stop 'init_gotm'
92 write(0,*) 'FATAL ERROR: ', 'I could not read the "station" namelist'
   stop 'init_gotm'
93 write(0,*) 'FATAL ERROR: ', 'I could not read the "time" namelist'
   stop 'init_gotm'
94 write(0,*) 'FATAL ERROR: ', 'I could not read the "output" namelist'
   stop 'init_gotm'
95 write(0,*) 'FATAL ERROR: ', 'I could not read the "eqstate" namelist'
   stop 'init_gotm'
96 write(0,*) 'FATAL ERROR: ', 'I could not read the "turbulence" namelist'
   stop 'init_gotm'
   end subroutine init_gotm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Manage global time--stepping \label{timeLoop}
!
! !INTERFACE:
   subroutine time_loop()
!
! !DESCRIPTION:
! This internal routine is the heart of the code. It contains
! the main time--loop inside of which all routines required 
! during the time step are called. The following main processes are 
! successively triggered.
! \begin{enumerate}
!  \item The model time is updated and the output is prepared.
!  \item Air--sea interactions (flux, SST) are computed.
!  \item The time step is performed on the mean-flow equations
!        (momentum, temperature).
!  \item Some quantities related to shear and stratification are updated 
!        (shear--number, buoyancy frequency, etc).
!  \item Turbulence is updated depending on what turbulence closure
!        model has been specified by the user.
!  \item The results are written to the output files.
! \end{enumerate}
!
! Depending on macros set for the Fortran pre--processor, extra features
! like the effects of sea--grass or sediments are considered in this routine 
! (see \sect{sec:extra}).
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   ::n,count,mark,ios,day_cycle
!
!-----------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'time_loop'


   day_cycle=0  !This assumes that gotmrun.inp starts at local midnight
   cloud=0.5
   
   do n=MinN,MaxN
512   continue     

!SP 14/10/05 read in state vector 
      if (n==MinN) then
        OPEN(UNIT=73,FILE="current_optimal_state.dat",IOSTAT=ios)
	IF(ios/=0) THEN
		PRINT*,"Error during opening input file, stopping"; STOP
       	END IF
	DO count=0,150                 		
		READ (UNIT=73,*) T(count)
	END DO
        DO count=0,150                 		
		READ (UNIT=73,*) u(count)
	END DO
        DO count=0,150                 		
		READ (UNIT=73,*) v(count)
	END DO
        DO count=0,150                 		
		READ (UNIT=73,*) num(count)
	END DO
        DO count=0,150                 		
		READ (UNIT=73,*) nuh(count)
	END DO
        DO count=0,150                 		
		READ (UNIT=73,*) tke(count)
	END DO
        DO count=0,150
                READ (UNIT=73,*) eps(count)
        END DO
        DO count=0,150
                READ (UNIT=73,*) L(count)
        END DO
                READ (UNIT=73,*) I_0
                READ (UNIT=73,*) qb
	CLOSE (UNIT=73)
     end if
!SP 14/10/05 end read in state vector

      call update_time(n)
      call prepare_output(n)

!     all observations/data
      call get_all_obs(julianday,secondsofday,nlev,z)    

!SP 24/10/05 assimilate local midnight sst to initialise run
   if (n==MinN) then   
      if (assimilation_type.ne.0) then
            call assimilation(T(1:nlev),tprof(1:nlev),cloud,advect,int_cs)
      end if
   end if
!SP 24/10/05 end assimilate local midnight sst to initialise run

!     external forcing
      if( calc_fluxes ) then
         call set_sst(T(nlev))
         !call set_sst(tprof(nlev)) !use observed temperatures for flux calculations
      end if
    
      call air_sea_interaction(julianday,secondsofday)

      heat = -(qb+qe+qh)    !comment this line if using heat obs
      tx = tx/rho_0
      ty = ty/rho_0

!     meanflow integration starts
      call updategrid(nlev,dt,zeta)
      call coriolis(nlev,dt)
      SS = 0.
      call uequation(nlev,dt,cnpar,tx,num,PressMethod)
      call vequation(nlev,dt,cnpar,ty,num,PressMethod)
      call extpressure(PressMethod,nlev)
      if (int_press_method .ne. 0) call intpressure(nlev)
      call friction(kappa,avmolu,tx,ty)

      call calc_seagrass(nlev,dt)

!SP - below if statement added to allow TKE to build up before heat is added
!      if(n.lt.MinN+480) then
!         PRINT*,'in the loop'
!         heat=0.0
!      end if    

      if (s_prof_method .ne. 0.) call salinity(nlev,dt,cnpar,nuh)
!start SP : use if want all available observed salt profiles
!      s=sprof
!end SP
      if (t_prof_method .ne. 0.) then

         call temperature(nlev,dt,cnpar,I_0,heat,advect,qdir_frac,qdiff_frac,cosr,nuh,rad)
end if
!start SP : use if want all available observed temp profiles
!      T=tprof
!end SP

call integrated_fluxes(dt,int_cs)       !SP previously below do_output 18/05/05


     call stratification(nlev,buoy_method,dt,cnpar,gravity,rho_0,nuh)

      call calc_sediment(nlev,dt)

      select case (turb_method)
         case (0)
            call convectiveadjustment(nlev,num,nuh,const_num,const_nuh,&
                                      buoy_method,gravity,rho_0)
         case (1)
            write(0,*) '... turb_method=1 is not coded yet.'
            write(0,*) 'Choose  turb_method=0 or turb_method=2.'
            write(0,*) 'Program execution stopped ...'
            stop 'time_loop'
         case (2)
            call production(nlev,alpha,num,nuh)
            call stabilityfunctions(nlev,NN,SS,1)
            call do_tke(nlev,dt,u_taus,u_taub,z0s,z0b,h,SS,NN,P,B)
            call lengthscale(nlev,dt,z0b,z0s,u_taus,u_taub,depth,h,NN,P,B)
            call turbulence_adv(nlev,dt,h)
            call kolpran(nlev)
         case default
      end select

      call internal_wave(nlev,NN,SS)
      tx=tx*rho_0
      ty=ty*rho_0


!SP 17/10/05 evaluate cost function
      day_cycle=day_cycle+1
      if (day_cycle==1440) then
         cost_fctn=T(obs_level)-tprof(obs_level)
         OPEN(UNIT=75,FILE="costfunction_new.dat",IOSTAT=ios)
         IF(ios/=0) THEN
		PRINT*,"Error during opening input file, stopping"; STOP
       	 END IF
         WRITE (UNIT=75,*) T(obs_level),tprof(obs_level),swr_error,wind_error
         CLOSE (UNIT=75)
      end if
      if (day_cycle==2880) then
         cost_fctn=cost_fctn**2+(T(obs_level)-tprof(obs_level))**2
         day_cycle=0
      end if
!SP 17/10/05 end evaluate cost function


      call do_output(n,nlev)


!kbk - next version      call heat_content()

!     write out variances
      if(variances) then
         call do_variances(nlev)
      end if

!     Diagnostic output
      if(diagnostics) then
         call do_diagnostics(n,nlev,buoy_method,dt, &
                             u_taus,u_taub,I_0,heat)
      end if


!SP 14/10/05 print out state vector
      if (n==1440) then
        OPEN(UNIT=73,FILE="end_state.dat",IOSTAT=ios)
	IF(ios/=0) THEN
		PRINT*,"Error during opening input file, stopping"; STOP
       	END IF
	DO count=0,150                 		
		WRITE (UNIT=73,*) T(count)
	END DO
        DO count=0,150                 		
		WRITE (UNIT=73,*) u(count)
	END DO
        DO count=0,150                 		
		WRITE (UNIT=73,*) v(count)
	END DO
        DO count=0,150                 		
		WRITE (UNIT=73,*) num(count)
	END DO
        DO count=0,150                 		
		WRITE (UNIT=73,*) nuh(count)
	END DO
        DO count=0,150                 		
		WRITE (UNIT=73,*) tke(count)
	END DO
        DO count=0,150
                WRITE (UNIT=73,*) eps(count)
        END DO
        DO count=0,150
                WRITE (UNIT=73,*) L(count)
        END DO
                WRITE (UNIT=73,*) I_0
                WRITE (UNIT=73,*) qb
	CLOSE (UNIT=73)
      end if
!SP 14/10/05 print out state vector

   end do
   write(0,*) "------------------------------------------------------------------------"

   return
   end subroutine time_loop
!EOC
!-----------------------------------------------------------------------
   subroutine assimilation(T,tprof,cloud,advect,int_cs)

     IMPLICIT NONE

     double precision, intent(in)        :: tprof(1:nlev)
     double precision, intent(inout)     :: T(1:nlev),cloud,int_cs
     double precision, intent(out)       :: advect(1:nlev)
 
     double precision        :: delta,heat_adv,cloud_new,sum_h_sq=0.
     double precision        :: pen,h_content=0.,sum_advect=0.
     integer                 :: k,i,n,assim_depth,z_depth

!-----------------------------------------------------------------------   
select case (assimilation_type)

case(1)
   call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)

case(2)
   call assimilate_profile(tprof(1:nlev),T(1:nlev))

case(3)
   call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)
   call assimilate_profile(tprof(1:nlev),T(1:nlev))

case(4)
!for cloud adjustment
         cloud=cloud+(T(obs_level)-tprof(obs_level))/cloud_gradient   
         if (cloud.lt.0.0) then
            cloud=0.0
         else if (cloud.gt.1.0) then
            cloud=1.0
         end if
         !PRINT*,cloud
!end for cloud adjustment

      call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)
      call assimilate_profile(tprof(1:nlev),T(1:nlev))

case(5)

   call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)
   
!find heat supplied in data assimlation
   DO n=assim_depth,nlev
      h_content=h_content+delta*h(n)
   END DO
   h_content=cp*rho_0*h_content
!end heat supplied in data assimilation

!calculate cloud value
!   pen=(1.-A*exp(-depth0/g1)-(1-A)*exp(-depth0/g2))
   pen=(1.-A*exp(-assim_depth/g1)-(1-A)*exp(-assim_depth/g2)) !mixed layer
   cloud_new=cloud+(h_content/(-.62*int_cs*pen))
PRINT*,cloud_new
   if(cloud_new.gt.1) then
      cloud=1.
   else if(cloud_new.lt.0) then
      cloud=0.
   else
      cloud=cloud_new
   end if
!end calculate cloud value

advect=0.
h_content=0.
int_cs=0.

case(6)

   call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)

!find heat supplied in data assimlation
   DO n=assim_depth,nlev
      h_content=h_content+delta*h(n)
   END DO
   h_content=cp*rho_0*h_content
!end heat supplied in data assimilation

!calculate advection term 
   DO n=assim_depth,nlev              !1,nlev
      sum_h_sq=sum_h_sq+(h(n)**2)
   END DO

   DO n=assim_depth,nlev              !1,nlev
      advect(n)=advect(n)+h(n)*h_content/(cp*rho_0*86400*sum_h_sq)
   END DO
!only for mixed layer
   DO n=1,assim_depth-1        
      advect(n)=0.
   END DO
!end only for mixed layer
!end calculate advection term 

sum_h_sq=0
h_content=0.
int_cs=0.

case(7)

   call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)

!find heat supplied in data assimlation
   DO n=assim_depth,nlev
      h_content=h_content+delta*h(n)
   END DO
   h_content=cp*rho_0*h_content
!end heat supplied in data assimilation

   pen=(1.-A*exp(-depth0/g1)-(1-A)*exp(-depth0/g2)) 
   DO n=assim_depth,nlev               !1,nlev
      sum_advect=sum_advect+(advect(n)*h(n))
   END DO
   cloud_new=cloud+((h_content+sum_advect*cp*rho_0*86400)/(-.62*int_cs*pen))
   
   DO n=assim_depth,nlev               !1,nlev
      sum_h_sq=sum_h_sq+(h(n)**2)
   END DO

!calculate cloud & advection term
   if(cloud_new.gt.1) then
      cloud_new=1.
      heat_adv=h_content+.62*int_cs*pen*(cloud_new-cloud)
      DO n=assim_depth,nlev            !1,nlev
         advect(n)=advect(n)+h(n)*heat_adv/(cp*rho_0*86400*sum_h_sq)
      END DO
!only for mixed layer
   DO n=1,assim_depth-1        
      advect(n)=0.
   END DO
!end only for mixed layer
   else if(cloud_new.lt.0) then
      cloud_new=0.
      heat_adv=h_content+.62*int_cs*pen*(cloud_new-cloud)
      DO n=assim_depth,nlev            !1,nlev
         advect(n)=advect(n)+h(n)*heat_adv/(cp*rho_0*86400*sum_h_sq)
      END DO
!only for mixed layer
   DO n=1,assim_depth-1        
      advect(n)=0.
   END DO
!end only for mixed layer
   else
      advect=0.

   end if
   cloud=cloud_new
!end calculate cloud & advection term

sum_h_sq=0.
sum_advect=0.   
h_content=0.
heat_adv=0.
int_cs=0.

end select

     end subroutine assimilation
!-----------------------------------------------------------------------
     subroutine assimilate_sst(tprof,T,delta,assim_depth)

       IMPLICIT NONE

       double precision, intent(in)    :: tprof(1:nlev)
       double precision, intent(inout) :: T(1:nlev)
       double precision, intent(out)   :: delta
       integer, intent(out)            :: assim_depth

!       double precision  :: delta,assim_depth
       integer           :: i,count
     
!-----------------------------------------------------------------------
count=count+1
!Assimilate SST into mixed layer
      if(count==sst_obs) then
         delta=tprof(obs_level)-T(obs_level)
         if(mld(T(1:nlev)).lt.obs_level) then
            assim_depth=mld(T(1:nlev))
         else
            assim_depth=obs_level
         end if
         
         DO i=assim_depth,nlev
            T(i)=T(i)+delta
         END DO
!end assimilate SST into mixed layer
         count=0
      end if


     end subroutine assimilate_sst
!-----------------------------------------------------------------------
     subroutine assimilate_profile(tprof,T)

       IMPLICIT NONE

       double precision, intent(in)    :: tprof(1:nlev)
       double precision, intent(inout) :: T(1:nlev)

       double precision  :: tprof_assim(1:nlev)
       integer           :: i,count
     
!-----------------------------------------------------------------------
count=count+1

!Re-initialise temperature profile 
!      if(count.gt.profile_obs-7) then
!            do i=1,nlev
!               !sum profile at each time step for a week
!               tprof_assim(i)=tprof_assim(i)+tprof(i)
!            end do
!         if(count==profile_obs) then
!            count=0
!            do i=1,nlev
!               !previous weeks mean is used to reinitialise
!               tprof_assim(i)=tprof_assim(i)/7
!               T(i)=tprof_assim(i)          
!            end do
!            tprof_assim=0.
!         end if
!      end if
!end reinitialise temp profile

!without averaging
      if(count==profile_obs) then
            count=0
            do i=1,nlev
               T(i)=tprof(i)          
            end do
         end if
!end without averaging

     end subroutine assimilate_profile
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
     function mld(T)
       double precision, intent(in):: T(1:nlev)
       integer                     :: twenty_level, max_T_level
       integer                     :: i,count,mld

       max_T_level=nlev
       twenty_level=66

!Find the grid level of the maximum temperature in the top 20m
         DO i=nlev-1,twenty_level,-1
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
!BOP
!
! !IROUTINE: The run is over --- now clean up.
!
! !INTERFACE:
   subroutine clean_up()
!
! !DESCRIPTION:
! This function is just a wrapper for the external routine 
! {\tt close\_output()} discussed in \sect{sec:output}. All open files
! will be closed after this call.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'clean_up'
  
   call close_output()

   return
   end subroutine clean_up
!EOC

!-----------------------------------------------------------------------

   end module gotm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
