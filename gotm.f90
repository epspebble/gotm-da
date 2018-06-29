!SP
!25/11/04 : cloud=cloud+(sst(model)-sst(obs))/.21 use this formulae once a day
!02/05 : use sst(obs) to change mixed layer depth by adjusting profiles.
!05/05 : subroutine called assimilation coded up which include 4 options
!10/05 : for assimilation at local midnight, run must start local midnight
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gotm --- the general framework \label{sec:gotm}
!
! !INTERFACE:
module gotm

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
  use airsea, only: calc_fluxes
  use airsea, only: tx,ty,heat,qh,qb,qe,cloud,int_cs,I_0, coszen !WT I_0, coszen used in demarcating assimilation cycle boundaries.
  use airsea, only: qdir_frac, qdiff_frac,cosr !HK 

  use meanflow
  use turbulence
  use observations !WT 20171121 v3c, This now includes next_tprof_secondsofday
  use output
  use time
  use mtridiagonal
  use eqstate
  !WT Uncomment this after moving the subroutines over, checking what used to be USE-ed, and what subroutines to export.
!  use assimilation

  use sediment
  use seagrass

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

  ! WT 20170331
  integer, parameter                  :: unit_daily_stat=71
  integer, parameter                  :: unit_sst_event=72
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
  integer                     :: sst_obs,profile_obs,obs_level
  double precision, public    :: advect(1:150)

  !WT 20171112 Stuff related to assimilation should be moved to a new module. Now for expediency.
  integer, public :: assimilation_type, assim_window

  !SP 19/09/05 test   
  double precision            :: j_one,j_one_b,j_two,j_three,j_four,j_five

  !WT 20170331 temporary variables
  character(len=19) :: tmp_str
  integer :: jul0,jul1,yyyy,mm,dd 

  !WT 20170406 daily stat output variables
!  integer :: day_number
  double precision :: daily_SST_max = -99, daily_SST_min_day = 99, daily_SST_min_night = 99, daily_SST_mean = 0.
  integer :: lsecs_assim_time, lsecs_solar_noon, lsecs_SST_max, lsecs_SST_min_day, lsecs_SST_min_night
  
  logical :: first
  double precision            :: sst_save=0., d_sst=0.
  

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
         assimilation_type,cloud_gradient,sst_obs,profile_obs,obs_level,assim_window, &    !SP
         daily_stat_fn, sst_event_fn !WT 2017046
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
    write(0,*) '       ', 'The station ',trim(name),' is located at (lat,long) ', &
         latitude,longitude
    write(0,*) '       ', trim(name)

    ! WT 20170331
    ! Open the extra files to store daily assimilation and SST events.
    !print *, 'Writing events to ',trim(daily_stat_fn), ' and ', trim(sst_event_fn), '... '
    open(unit=unit_daily_stat,file=daily_stat_fn,status='replace')
!    open(unit=unit_sst_event,file=sst_event_fn,status='replace')

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
         depth,nlev,z,h,latitude,longitude)

    ! Copy the first profiles to be the values of s, t, u, v before time loop begins.
    s = sprof
    t = tprof
    u = uprof
    v = vprof
    call init_output(title,nlev,latitude,longitude)
    call init_air_sea(namlst,latitude,longitude,julianday,secondsofday)

    ! WT 20171112 Place holder for moving stuff into a new assimilation module. 
    call init_assimilation() 
    
    !  Initialise each of the extra features/modules

    call init_sediment(namlst,'sediment.inp',unit_sediment,nlev, &
         gravity,rho_0)

    call init_seagrass(namlst,'seagrass.inp',unit_seagrass,nlev,h)

    write(0,*) '       ', 'done.'
    write(0,*) "------------------------------------------------------------------------"

    return

90  write(0,*) 'FATAL ERROR: ', 'I could not open gotmrun.inp for reading'
    stop 'init_gotm'
91  write(0,*) 'FATAL ERROR: ', 'I could not read the "model_setup" namelist'
    stop 'init_gotm'
92  write(0,*) 'FATAL ERROR: ', 'I could not read the "station" namelist'
    stop 'init_gotm'
93  write(0,*) 'FATAL ERROR: ', 'I could not read the "time" namelist'
    stop 'init_gotm'
94  write(0,*) 'FATAL ERROR: ', 'I could not read the "output" namelist'
    stop 'init_gotm'
95  write(0,*) 'FATAL ERROR: ', 'I could not read the "eqstate" namelist'
    stop 'init_gotm'
96  write(0,*) 'FATAL ERROR: ', 'I could not read the "turbulence" namelist'
    stop 'init_gotm'
    
  end subroutine init_gotm
  !EOC

  subroutine init_assimilation()
    !WT Place holder.
    ! Uses:
    IMPLICIT None
    print *, 'init_assimilation() called. Doing nothing as of now.'
  end subroutine init_assimilation
  
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
    integer                   ::n,ios,day_cycle
    character(len=45)          ::fname
    character(len=24)     ::c1='OBS/ECMWF/initial_state_'
    character(len=2)       ::c2='N_'
    character(len=12)       ::c3='E_010106.dat'
    character(len=2)      ::c4='S_'
    integer                   ::lat,lon
    logical :: mark=.false., first=.true.
    !
    !-----------------------------------------------------------------------
    !BOC

    ! Initialize by theoretical timezone.
    day_cycle = tz(longitude)*3600/timestep 
    if (day_cycle < 0) then
       day_cycle = day_cycle + 86400/timestep !86400/timestep=2880 if timestep=30 
    endif


    write(0,*) '   ', 'time_loop'
    do n=MinN,MaxN
512    continue     

       !!WT 20180623 Not necessary.  
       ! if (first .and. mark) then
       !    print *, 'Line 328 is called.'
       !    call write_daily_stats(unit_daily_stat,0)
       !    first = .false.
       ! endif
       
       ! This is the final step for something below... 
       if (n .eq. MaxN) then
          ! Print the last set of stat.
          call write_daily_stats(unit_daily_stat,-1)
       endif
       
       call update_time(n)

       !WT Attempt a daily assimilation.
       call assimilate_daily(mark,first,day_cycle)
       
       !SP 17/10/05 evaluate cost function

       !SP 17/10/05 end evaluate cost function

       call prepare_output(n)

       !     all observations/data
       call get_all_obs(julianday,secondsofday,nlev,z)    

       !     external forcing
       !      if( calc_fluxes ) then
       !SP 16/03/06 needed for the cool skin parameterisation
       call set_sst(T(nlev))
       ! call set_sst(tprof(nlev)) !use observed temperatures for flux calculations?
       !      end if

       call air_sea_interaction(julianday,secondsofday)

       heat = -(qb+qe+qh)     !comment this line if heat from file SP      
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
          !SP comment below if want persistence only (i.e. null hypothesis)
          call temperature(nlev,dt,cnpar,I_0,heat,advect,qdir_frac,qdiff_frac,cosr,nuh,rad)
       end if

       !start SP : use if want all available observed temp profiles
       !      T=tprof
       !end SP

       !call sst_event(T(nlev))

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

    end do
    write(0,*) "------------------------------------------------------------------------"

    return
  end subroutine time_loop
  !EOC
  !-----------------------------------------------------------------------

  ! BEGIN Prospective assimilation module section.
  subroutine assimilate_daily(mark,first,day_cycle)
    logical, intent(inout) :: mark, first
    integer, intent(inout) :: day_cycle

    ! LOCAL VARIABLES
    integer :: ljul, lsecs


    ! Find date and time in local timezone.
    call UTC_to_local(julianday,secondsofday,longitude,ljul,lsecs)
    
    if (assimilation_type.ne.0) then

       ! Advance the counter whenever time_loop is called.
       day_cycle = day_cycle + 1

       ! find daily SST max / mins

       ! search SST min in local time 1:00 - 11:00, same as HCMR criterion when checking SEVIRI data
       if ((day_cycle .ge. 1*3600/timestep) .and. (day_cycle .le. 11*3600/timestep)) then
          ! Will only be run if assimilation does not occur at midnight, (assim_window .gt. 0), 
          ! because mark is reset to 0 after momentarily after being setting to 1.
          ! Also note the order of statements below.
          if (mark .eqv. .true.) then ! After assimilation event (i.e. sunrise if assimilating at sunrise)
             if (T(nlev) < daily_SST_min_day) then
                daily_SST_min_day = T(nlev)
                lsecs_SST_min_day = lsecs
             end if
          else if (mark .eqv. .false.) then ! After local midnight, before sunrise if assimilating at sunrise
             if (T(nlev) < daily_SST_min_night) then
                daily_SST_min_night = T(nlev)
                lsecs_SST_min_night = lsecs
             end if
          end if
       end if

       ! search SST max local time 11:00 - 21:00, same as HCMR criterion when checking SEVIRI data
       if ((day_cycle .ge. 11*3600/timestep) .and. (day_cycle .le. 21*3600/timestep)) then
          if (T(nlev) > daily_SST_max) then
             daily_SST_max = T(nlev)
             lsecs_SST_max = lsecs
          end if
       end if

       ! Aggregate SST over an assimilation / day cycle. 
       ! if (count .eq. 0) then
       !    print *, daily_SST_mean
       ! end if
       daily_SST_mean = daily_SST_mean + T(nlev)
       !count = count + 1
       ! if (T(nlev) .lt. 20) then
       !    print *, T(nlev)
       ! end if

       ! 'mark' == 0 means assimilation not yet done in this cycle.
       if (.not. mark) then
          !WT We assume simulation begins at 00:00:00.

          select case (assim_window)

          case (0) ! assimilate at "midnight"
             !WT "midnight" means day_cycle = 0, i.e. local, not UTC, but it may not
             ! correspond to 00:00:00 in tprof.dat, which more likely in UTC than in local time.
             ! Check needed.
             if (day_cycle .eq. 86400/timestep) then
                !!! The daily_stats is not outputted for this assim_window. Not implemented.
                call assimilation(T(1:nlev),S(1:nlev),tprof(1:nlev),sprof(1:nlev),cloud,advect,int_cs)
                lsecs_assim_time = lsecs
                mark = .true. 
             endif

          case (1) ! assimiluate at "sunrise"
             !WT 2016/09/13 "sunrise" means the first moment after
             ! 00:00:00 in tprof.dat that I_0 is found to be nonzero.
             if (I_0.gt.1) then
                if (.not. first) then
                   ! Theses result are for the assim cycle that is just completed, i.e. the previous day.
                   call write_daily_stats(unit_daily_stat)
                end if

                ! Also reset aggregator variables
                !count = 0
                daily_SST_mean = 0 
                daily_SST_max = -99
                daily_SST_min_day = 99
                daily_SST_min_night = 99

                call assimilation(T(1:nlev),S(1:nlev),tprof(1:nlev),sprof(1:nlev),cloud,advect,int_cs)
                ! Here the new assimilation cycle really begins. Record the time.
                lsecs_assim_time = lsecs
                mark = .true.
             endif

          case (2) ! assimilate at the tprof timestamp
             ! Basic logic:
             !      If the secondsofday greater than tprof_secondsofday, assimlate.
             if (time_diff(julianday,secondsofday,next_tprof_julianday,next_tprof_secondsofday) .ge. 0) then
                !   print *, julianday, next_tprof_secondsofday
                !   print *, secondsofday, next_tprof_secondsofday
                call assimilation(T(1:nlev),S(1:nlev),tprof(1:nlev),sprof(1:nlev),cloud,advect,int_cs)
                mark = .true.
             endif

          case (3) ! assimilate at "sunrise" version two, when coszen changes sign from - to +.
             if (coszen.gt.0) then
                !! WT 20180319 The following code is a copy-and-paste of case(1).

                ! The previous assimilation cycle has just completed. Write down the daily stats now.
                if (first) then
                   ! This result is for the part of the day before the first assimilation event.
                   call write_daily_stats(unit_daily_stat,0)
                else
                   ! Theses result are for the assim cycle that is just completed, i.e. the previous day.
                   call write_daily_stats(unit_daily_stat)
                end if

                ! Also reset aggregator variables
                !count = 0
                daily_SST_mean = 0 
                daily_SST_max = -99
                daily_SST_min_day = 99
                daily_SST_min_night = 99

                call assimilation(T(1:nlev),S(1:nlev),tprof(1:nlev),sprof(1:nlev),cloud,advect,int_cs)
                ! WT Here the new assimilation cycle really begins. Record the time.
                ! This value may be negative (more likely) or greater than 60*60*24 = 86400 (less likely) which depend on timezone,
                ! when this happens, the true lsecs_assim_time in that day is the positive remainder after division by 86400, but the
                ! sign or value the quotient has significance being that this occurs in a day before or later compared to the date in
                ! UTC timezone.
                lsecs_assim_time = lsecs
                mark = .true.
             end if
          end select
       endif

       ! restart day_cycle
       if (day_cycle .eq. 86400/timestep) then

          ! All these variables should have been updated by now: 
          !       solar_noon, daily_SST_max, lsecs_SST_max, daily_SST_min, lsecs_SST_min

          ! the summing of SST values have been done from day_cycle=0, ... day_cycle=86400/timestep-1, so we find average now.

          if (first) then
             daily_SST_mean = daily_SST_mean/(86400-tz(longitude)*3600)*timestep 
          else
             daily_SST_mean = daily_SST_mean/86400*timestep
          end if

          ! Reset counters at local midnight every day.
          if (first) then ! Do not yell error if it's the first day.
             first = .false.
             day_cycle = 0
          else
             ! Check whether the assimilation was done or not every day for
             ! assim_window that depend on local day cycle, currently:
             ! 0 (local midnight) and
             ! 1 (local sunrise), AND
             ! 2 (data timestamp)
             if ((.not.mark) .and. (assim_window.le.3)) then
                stop "Assimilation not performed in one full cycle. STOP"
             endif

             ! Resetting both the day cycle number and the marker for daily assimilation.
             day_cycle = 0
             mark = .false.
          endif
       endif
       
    endif

    ! Obsolete comments
        !SP 30/01/06 calculate local time based on longitude
    !This assumes that gotmrun.inp starts at midnight
    !if (longitude.gt.352) then
    !   day_cycle=0
    !else if ((longitude.gt.337).AND.(longitude.le.352)) then
    !   day_cycle=120
    !else if ((longitude.gt.322).AND.(longitude.le.337)) then
    !   day_cycle=240
    !else if ((longitude.gt.307).AND.(longitude.le.322)) then
    !   day_cycle=360
    !else if ((longitude.gt.292).AND.(longitude.le.307)) then
    !   day_cycle=480
    !else if ((longitude.gt.277).AND.(longitude.le.292)) then
    !   day_cycle=600
    !else if ((longitude.gt.262).AND.(longitude.le.277)) then
    !   day_cycle=720
    !else if ((longitude.gt.247).AND.(longitude.le.262)) then
    !   day_cycle=840
    !else if ((longitude.gt.232).AND.(longitude.le.247)) then
    !   day_cycle=960
    !else if ((longitude.gt.217).AND.(longitude.le.232)) then
    !   day_cycle=1080
    !else if ((longitude.gt.202).AND.(longitude.le.217)) then
    !   day_cycle=1200
    !else if ((longitude.gt.187).AND.(longitude.le.202)) then
    !   day_cycle=1320
    !else if ((longitude.gt.172).AND.(longitude.le.187)) then
    !   day_cycle=1440
    !else if ((longitude.gt.157).AND.(longitude.le.172)) then
    !   day_cycle=1560
    !else if ((longitude.gt.142).AND.(longitude.le.157)) then
    !   day_cycle=1680
    !else if ((longitude.gt.127).AND.(longitude.le.142)) then
    !   day_cycle=1800
    !else if ((longitude.gt.112).AND.(longitude.le.127)) then
    !  day_cycle=1920
    !else if ((longitude.gt.97).AND.(longitude.le.112)) then
    !   day_cycle=2040
    !else if ((longitude.gt.82).AND.(longitude.le.97)) then
    !   day_cycle=2160
    !else if ((longitude.gt.67).AND.(longitude.le.82)) then
    !   day_cycle=2280
    !else if ((longitude.gt.52).AND.(longitude.le.67)) then
    !   day_cycle=2400
    !else if ((longitude.gt.37).AND.(longitude.le.52)) then
    !   day_cycle=2520
    !else if ((longitude.gt.22).AND.(longitude.le.37)) then
    !   day_cycle=2640
    !else if ((longitude.gt.7).AND.(longitude.le.22)) then
    !   day_cycle=2760
    !else

    !day_cycle=0
    !mark=0

    !end if
    !SP 30/01/06 end calculate local time


    !SP 03/04/06 print out initial state
    !lat=nint(latitude)
    !lon=nint(longitude)
    !if(lat.gt.0) then
    !   write(fname,'(a24,i3.3,a2,i3.3,a12)') c1,lat,c2,lon,c3
    !else
    !   lat=-lat
    !   write(fname,'(a24,i3.3,a2,i3.3,a12)') c1,lat,c4,lon,c3
    !end if

    !        OPEN(UNIT=73,FILE=fname,IOSTAT=ios)
    !	IF(ios/=0) THEN
    !		PRINT*,"Error during opening input file, stopping"; STOP
    !       END IF
    !                WRITE (UNIT=73,*) tprof(nlev)-T(nlev)
    !	DO count=0,150                 		
    !		WRITE (UNIT=73,*) T(count)
    !	END DO
    !        	DO count=0,150                 		
    !		WRITE (UNIT=73,*) S(count)
    !	END DO
    !        DO count=0,150                 		
    !		WRITE (UNIT=73,*) u(count)
    !	END DO
    !        DO count=0,150                 		
    !		WRITE (UNIT=73,*) v(count)
    !	END DO
    !        DO count=0,150                 		
    !		WRITE (UNIT=73,*) num(count)
    !	END DO
    !        DO count=0,150                 		
    !		WRITE (UNIT=73,*) nuh(count)
    !	END DO
    !        DO count=0,150                 		
    !		WRITE (UNIT=73,*) tke(count)
    !	END DO
    !        DO count=0,150
    !                WRITE (UNIT=73,*) eps(count)
    !        END DO
    !        DO count=0,150
    !                WRITE (UNIT=73,*) L(count)
    !        END DO
    !                WRITE (UNIT=73,*) I_0
    !                WRITE (UNIT=73,*) heat
    !                WRITE (UNIT=73,*) tx
    !                WRITE (UNIT=73,*) ty
    !	CLOSE (UNIT=73)
    !SP 03/04/06 end print out initial state


    !SP 24/10/05 assimilate midnight
    ! if (assimilation_type.ne.0) then
    !    if((mark.eq.0).AND.(day_cycle>2160).AND.(I_0>1)) then
    !       call assimilation(T(1:nlev),S(1:nlev),tprof(1:nlev),sprof(1:nlev),cloud,advect,int_cs)
    !    elseif(day_cycle==0) then
    !       call assimilation(T(1:nlev),S(1:nlev),tprof(1:nlev),sprof(1:nlev),cloud,advect,int_cs)
    !    end if
    !    day_cycle=day_cycle+1
    !    if((mark.eq.0).AND.(day_cycle>2160).AND.(I_0>1)) then
    !       day_cycle=1
    !       mark=1
    !    elseif(day_cycle==2880) then
    !       day_cycle=0
    !    end if
    ! end if
    !SP 24/10/05 end assimilate midnight

    !WT 2016-09-24 Assimilate sunrise/midnight. Version as of 2016-09-11,
    ! but with case (2): SST turnaround removed, and the `event` signal
    ! variable removed, and re-using `mark`.
    !
    ! Key Ideas
    !   1. We want to assimilate exactly once in every diurnal
    !      cycle, which is approximately 24 hours (it could be
    !      slightly less than 24 hours if the sun rises later in a
    !      day after summer solstice.)
    !   2. Both the timezone (determined by longitude) and the
    !      distance from equator (latitude) affects how we interpret
    !      the data. For the sake of eventually running GOTM on a
    !      grid over the whole world, our code will assume the data
    !      is given in GMT.
    !   3. As a result of assuming GMT, at every location not in
    !      this very timezone, 00:00:00 will not correspond to local
    !      midnight.
    !
    ! Potential problems
    !   1. The `day_cycle = 2880`, i.e. number of half-minutes per day, which depends on the assumption
    !      that `dt = 30` in the model_setup namelist in gotmrun.inp. It will be wrong if `dt` is not 30.
    !      A safer approach is to either find out what dt is or use the actual time variable value
    !      in the time loop, which is a top todo in the future.
    !   2. Is using I_0 > 1 (why not I_0 > 0? maybe for some good reason, but let's keep it as it is) a
    !      robust way to indicate sunrise?
    !   3. We need to confirm whether using the condition (day_cycle.eq.2880) to reset both the counter
    !      and marker really correspond to midnight, and whether I_0 would be nonzero after this point
    !      but not due to sunrise.
    !   4. This code will fail to run the historical buoy site, say, at ARABIAN.
    
    
  end subroutine assimilate_daily
  
  subroutine assimilation(T,S,tprof,sprof,cloud,advect,int_cs)

    IMPLICIT NONE

    double precision, intent(in)        :: tprof(1:nlev),sprof(1:nlev)
    double precision, intent(inout)     :: T(1:nlev),S(1:nlev),cloud,int_cs
    double precision, intent(out)       :: advect(1:nlev)

    double precision        :: delta,heat_adv,cloud_new,sum_h_sq=0.
    double precision        :: pen,h_content=0.,sum_advect=0.
    integer                 :: k,i,n,assim_depth,z_depth

    !-----------------------------------------------------------------------   

    select case (assimilation_type)

    case(1)
       call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)

    case(2)
       call assimilate_profile(tprof(1:nlev),T(1:nlev),sprof(1:nlev),S(1:nlev))

    case(3)
       call assimilate_profile(tprof(1:nlev),T(1:nlev),sprof(1:nlev),S(1:nlev))
       call assimilate_sst(tprof(1:nlev),T(1:nlev),delta,assim_depth)       

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
       call assimilate_profile(tprof(1:nlev),T(1:nlev),sprof(1:nlev),S(1:nlev))

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
       !delta=22.5750-T(obs_level)
       if(mld(T(1:nlev)).lt.obs_level) then
          assim_depth=mld(T(1:nlev))
       else
          assim_depth=obs_level
       end if
       print *,'Assimilate SST into mixed layer from grid level ', assim_depth, ' till ', nlev
       DO i=assim_depth,nlev
          T(i)=T(i)+delta
       END DO
       !end assimilate SST into mixed layer
       count=0
    end if


  end subroutine assimilate_sst
  !-----------------------------------------------------------------------
  subroutine assimilate_profile(tprof,T,sprof,S)

    IMPLICIT NONE

    double precision, intent(in)    :: tprof(1:nlev),sprof(1:nlev)
    double precision, intent(inout) :: T(1:nlev),S(1:nlev)

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
          S(i)=sprof(i)
       end do

       !! Organize under day_cycle loop, and print several information.
       ! write(unit_daily_stat,*) tmp_str, day_number, secondsofday, T(nlev)

       !WT 2016-09-24
       call write_time_string(julianday,secondsofday,tmp_str)
       write(0,*) 'Temperature and salinity profiles assimilated at ', tmp_str, dayofyear,secondsofday
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
  ! END Prospective assimilation module portion.
  
  subroutine write_daily_stats(unit, special)
    integer :: unit, yrdays
    integer, optional :: special

    if (present(special)) then
       select case (special)
       case (0) ! First record
          ! This line of stat is for the current incomplete assimilation cycle.
          ! Should also expect daily_SST_min_day and daily_SST_max to still be the fill in values of 99 and -99 repectively.
          write(unit_daily_stat,714) dayofyear, lsecs_assim_time, & 
               lsecs_SST_min_day, daily_SST_min_day, lsecs_SST_max, daily_SST_max, &
               lsecs_SST_min_night, daily_SST_min_night
       case (-1) ! Last record
          !This line of stats is for the current incomplete assim cycle
          write(unit_daily_stat,714) dayofyear, lsecs_assim_time, & 
               lsecs_SST_min_day, daily_SST_min_day, lsecs_SST_max, daily_SST_max, lsecs_SST_min_night, daily_SST_min_night
       case default
          stop 'You should not see this. Something is wrong.'
       end select
    else if (dayofyear .eq. 1) then
       ! Would have written dayofyear-1 (which is 0) by next block. Instead, use the number of days of of last year.
       call calendar_date(julianday,yyyy,mm,dd)
       if (mod(yyyy-1,4) .eq. 0) then
          yrdays = 366
       else
          yrdays = 365
       end if
       write(unit_daily_stat,714) yrdays, lsecs_assim_time, & 
            lsecs_SST_min_day, daily_SST_min_day, lsecs_SST_max, daily_SST_max, &
            lsecs_SST_min_night, daily_SST_min_night
    else
       ! This give the main lines of stats, which are for the assimilation cycle that is just completed, so dayofyear-1
       write(unit_daily_stat,714) dayofyear-1, lsecs_assim_time, & 
            lsecs_SST_min_day, daily_SST_min_day, lsecs_SST_max, daily_SST_max, &
            lsecs_SST_min_night, daily_SST_min_night
    endif
    
!713 format(I3,4(1x,I5,1x,F10.6))
714 format(I3,1x,I5,3(1x,I5,1x,F10.6))
    
  end subroutine write_daily_stats
  
  ! BEGIN Prospective event handling module.
  subroutine sst_event(sst)
    ! !DESCRIPTION:
    ! This function detects sst turnaround, and write time and temperature to a file.
    double precision :: sst
    if (sign(1.0d0,sst-sst_save) * sign(1.0d0,d_sst) .eq. -1.0d0) then ! sign change in two immedaite time steps
       call write_time_string(julianday,secondsofday,tmp_str)
       ! write(0,*) "SST turnaround occurs at ", tmp_str, sst
       write(unit_sst_event,*) tmp_str, dayofyear, secondsofday, sst
    endif
    
    d_sst = sst - sst_save    
    sst_save = sst
    
  end subroutine sst_event
    
  ! END Prospective event handling module.

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
