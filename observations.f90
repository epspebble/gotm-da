!-----------------------------------------------------------------------
!BOP
!
! !MODULE: observations --- the 'real' world \label{sec:observations}
!
! !INTERFACE:
   module observations
!
! !DESCRIPTION: 
!  This module provides the necessary subroutines for communicating
!  `observations' to GOTM. 
!  The module operates according to the general philosophy used in GOTM, 
!  i.e.\ it provides {\tt init\_observ\-ations()} to be called in the overall 
!  initialisation routine and {\tt get\_all\_obs()} to be called in the time 
!  loop to actually obtain the `observations'.
!  In addition to these subroutines the module also provides two routines 
!  for reading scalar-type observations and profile-type observations.
!  Each observation has a date stamp with the format {\tt yyyy-mm-dd hh:dd:mm}. 
!  The module uses the {\tt time} module (see \sect{sec:time}) 
!  to convert the time string to the 
!  internal time representation of GOTM.
!  Profiles are interpolated to the actual GOTM model grid.
!  Free format is used for reading-in the actual data.
!
! !USES:
   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_observations,get_all_obs,read_obs,read_profiles
!
! !PUBLIC DATA MEMBERS:
!
!  'observed' salinity profile
   double precision, public, dimension(:), allocatable   :: sprof

!  'observed' temperature profile
   double precision, public, dimension(:), allocatable   :: tprof

!  'observed' horizontal salinity  gradients
   double precision, public, dimension(:), allocatable   :: dsdx,dsdy

!  'observed' horizontal temperarure  gradients
   double precision, public, dimension(:), allocatable   :: dtdx,dtdy

!  internal horizontal pressure gradients
   double precision, public, dimension(:), allocatable   :: idpdx,idpdy

!  horizontal velocity profiles
   double precision, public, dimension(:), allocatable   :: uprof,vprof

!  observed profile of turbulent dissipation rates
   double precision, public, dimension(:), allocatable   :: epsprof

!  ralaxation times for salinity and temperature
   double precision, public, dimension(:), allocatable   :: SRelaxTau,TRelaxTau

!  sea surface elevation, sea surface gradients and height of velocity obs.
   double precision, public          :: zeta=0.,dpdx=0.,dpdy=0.,h_press=0

!  vertical advection velocity 
   double precision, public          :: w_adv=0.,w_height

!  Parameters for water classification - default Jerlov type I
   double precision, public          :: A=0.58,g1=0.35,g2=23.0
   double precision, public          :: chlo      !SP 20/02/06
   double precision, public          :: abp_coe,bb      !HX 11/05/2017 
   double precision, public          :: Trans_1
!HK/SH
!  Parameters for absorption of solar energy - Paulson and Simpson 1981
!  Clear water - Defant,A. 1961. - SH 07/11/02!
   double precision, public, dimension(9) :: &
     fsen = (/0.237,0.360,0.179,0.087,0.08 &
             ,0.0246,0.025,0.007,0.0004/)
   double precision, public, dimension(9) :: &
     zdeta = (/3.4849e1,2.2661,3.1486e-2,5.4831e-3,8.317e-4 &
              ,1.2612e-4,3.1326e-4,7.8186e-5,1.4427e-5/)

!------------------------------------------------------------------------------
!
! the following data are not all public,
! but have been included for clarity here
!
!------------------------------------------------------------------------------

!  Salinity profile(s)
   integer, public           :: s_prof_method=0
   character(LEN=255)   :: s_prof_file='sprof.dat'
   character(LEN=255)   :: init_s_prof='init_sprof.dat'  !SP 15/07/05
   double precision                  :: z_s1,s_1,z_s2,s_2
   double precision                  :: SRelaxTauM=0.
   double precision                  :: SRelaxTauS=0.
   double precision                  :: SRelaxTauB=0.
   double precision                  :: SRelaxSurf=0.
   double precision                  :: SRelaxBott=0.

!  Temperature profile(s)
   integer, public           :: t_prof_method=0
   character(LEN=255)   :: t_prof_file='tprof.dat'
   character(LEN=255)   :: init_t_prof='init_tprof.dat'  !SP 19/05/05
   double precision                  :: z_t1,t_1,z_t2,t_2
   double precision                :: TRelaxTauM=0.
   double precision                  :: TRelaxTauS=0.
   double precision                  :: TRelaxTauB=0.
   double precision                  :: TRelaxSurf=0.
   double precision                  :: TRelaxBott=0.

!  External pressure - 'press' namelist
   integer, public           :: ext_press_method=0,PressMethod=0
   character(LEN=255)   :: ext_press_file=''
   double precision, public          :: PressConstU=0.
   double precision, public          :: PressConstV=0.
   double precision, public          :: PressHeight=0.
   double precision, public          :: PeriodM=44714.
   double precision, public          :: AmpMu=0.
   double precision, public          :: AmpMv=0.
   double precision, public          :: PhaseMu=0.
   double precision, public          :: PhaseMv=0.
   double precision, public          :: PeriodS=43200.
   double precision, public          :: AmpSu=0.
   double precision, public          :: AmpSv=0.
   double precision, public          :: PhaseSu=0.
   double precision, public          :: PhaseSv=0.

!  Internal pressure - 'internal_pressure' namelist
   integer, public           :: int_press_method=0
   character(LEN=255)   :: int_press_file=''
   double precision, public          :: const_dsdx=0.
   double precision, public          :: const_dsdy=0.
   double precision, public          :: const_dtdx=0.
   double precision, public          :: const_dtdy=0.
   double precision                  :: const_idpdx=0.
   double precision                  :: const_idpdy=0.
   logical, public           :: s_adv=.false.
   logical, public           :: t_adv=.false.

!  Light extinction - the 'extinct' namelist
   integer, public                   :: extinct_method=1
   character(LEN=255)   :: extinct_file='extinction.dat'

!  Vertical advection velocity - 'w_advspec' namelist
   integer, public           :: w_adv_method=0
   double precision, public          :: w_adv0=0.
   character(LEN=255)   :: w_adv_file='w_adv.dat'
   integer, public           :: w_adv_discr=1

!  Sea surface elevations - 'zetaspec' namelist
   integer,public            :: zeta_method=0
   character(LEN=255)   :: zeta_file='zeta.dat'
   double precision, public          :: zeta_0=0.
   double precision, public          :: period_1=44714.
   double precision, public          :: amp_1=0.
   double precision, public          :: phase_1=0.
   double precision, public          :: period_2=43200.
   double precision, public          :: amp_2=0.
   double precision, public          :: phase_2=0.

!  Observed velocity profile profiles - typically from ADCP
   integer                   :: vel_prof_method=0
   CHARACTER(LEN=255)   :: vel_prof_file='velprof.dat'
   double precision, public          :: vel_relax_tau=3600.
   double precision, public          :: vel_relax_ramp=86400.

!  Observed dissipation profiles
   integer                   :: e_prof_method=0
   double precision                  :: e_obs_const=1.e-12
   CHARACTER(LEN=255)   :: e_prof_file='eprof.dat'

!  Buoyancy - 'bprofile' namelist
   double precision, public          :: b_obs_surf=0.,b_obs_NN=0.
   double precision, public          :: b_obs_sbf=0.

   double precision,public, parameter:: pi=3.141592654

! !DEFINED PARAMETERS:

!  Unit numbers for reading observations/data.
   integer, parameter        :: s_prof_unit=30
   integer, parameter        :: t_prof_unit=31
   integer, parameter        :: ext_press_unit=32
   integer, parameter        :: int_press_unit=33
   integer, parameter        :: extinct_unit=34
   integer, parameter        :: w_adv_unit=35
   integer, parameter        :: zeta_unit=36
   integer, parameter        :: vel_prof_unit=37
   integer, parameter        :: e_prof_unit=38




!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: observations.F90,v $
!  Revision 1.7  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.6  2003/03/28 08:08:21  kbk
!  removed tabs
!
!  Revision 1.5  2003/03/10 13:51:08  lars
!  changed intent(out) to intent(inout) for lines in read_profiles()
!
!  Revision 1.4  2003/03/10 08:51:58  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.3  2001/11/27 15:35:55  gotm
!  zeta_method now public - used by updategrid()
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer, parameter        :: READ_SUCCESS=1
   integer, parameter        :: END_OF_FILE=-1
   integer, parameter        :: READ_ERROR=-2
   integer, parameter        :: NOTHING=0
   integer, parameter        :: ANALYTICAL=1
   integer, parameter        :: CONSTANT=1
   integer, parameter        :: FROMFILE=2
   character(len=82)         :: cbuf
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the observation module
!
! !INTERFACE:
   subroutine init_observations(namlst,fn,julday,secs,depth,nlev,z,h,latitude,longitude)
!
! !DESCRIPTION:
!  The {\tt init\_observations()} subroutine basically reads the {\tt obs.inp}
!  file with a number of different namelists and takes actions according
!  to the specifications in the different namelists.
!  In this routine also memory is allocated to hold the 'observations'.
!  Finally, all variables are initialised to sane values, either by
!  reading from files, by prescribing constant values, or by using analytical 
!  expressions.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: julday,secs
   double precision, intent(in)                :: depth
   integer, intent(in)                 :: nlev
   double precision, intent(in)                :: z(0:nlev),h(0:nlev)
   double precision, intent(in)        :: latitude,longitude
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!
! !LOCAL VARIABLES:
   namelist /sprofile/ s_prof_method,z_s1,s_1,z_s2,s_2,s_prof_file,init_s_prof, &
                       SRelaxTauM,SRelaxTauB,SRelaxTauS,SRelaxBott,SRelaxSurf
   namelist /tprofile/ t_prof_method,z_t1,t_1,z_t2,t_2,t_prof_file,init_t_prof, &
                       TRelaxTauM,TRelaxTauB,TRelaxTauS,TRelaxBott,TRelaxSurf
   namelist /ext_pressure/                                      &
            ext_press_method,PressMethod,ext_press_file,        &
            PressConstU,PressConstV,PressHeight,                &
            PeriodM,AmpMu,AmpMv,PhaseMu,PhaseMv,                &
            PeriodS,AmpSu,AmpSv,PhaseSu,PhaseSv
   namelist /int_pressure/                                      &
            int_press_method,int_press_file,                    &
            const_dsdx,const_dsdy,const_dtdx,const_dtdy,        &
            const_idpdx,const_idpdy,s_adv,t_adv
   namelist /extinct/ extinct_method,extinct_file
   namelist /w_advspec/                                         &
            w_adv_method,w_adv_file,w_adv0,w_adv_discr                      
   namelist /zetaspec/                                          &
            zeta_method,zeta_file,zeta_0,                       &
            period_1,amp_1,phase_1,period_2,amp_2,phase_2
   namelist /velprofile/ vel_prof_method,vel_prof_file,         &
            vel_relax_tau,vel_relax_ramp
   namelist /eprofile/ e_prof_method,e_obs_const,e_prof_file
   namelist /bprofile/ b_obs_surf,b_obs_NN,b_obs_sbf

   integer                   :: rc,i
   double precision                  :: ds,db
!   SP
   integer                     ::ios, count
!
!-----------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'init_observations'

   open(namlst,file=fn,status='old',action='read',err=80)
   read(namlst,nml=sprofile,err=81)
   read(namlst,nml=tprofile,err=82)
   read(namlst,nml=ext_pressure,err=83)
   read(namlst,nml=int_pressure,err=84)
   read(namlst,nml=extinct,err=85)
   read(namlst,nml=w_advspec,err=86)
   read(namlst,nml=zetaspec,err=87)
   read(namlst,nml=velprofile,err=88)
   read(namlst,nml=eprofile,err=89)
   read(namlst,nml=bprofile,err=90)
   close(namlst)

   allocate(sprof(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (sprof)'
   sprof = 0.

   allocate(tprof(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (tprof)'
   tprof = 0.

   allocate(dsdx(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dsdx)'
   dsdx = 0.

   allocate(dsdy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dsdy)'
   dsdy = 0.

   allocate(dtdx(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dtdx)'
   dtdx = 0.

   allocate(dtdy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dtdy)'
   dsdy = 0.

   allocate(idpdx(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (idpdx)'
   idpdx = 0.

   allocate(idpdy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (idpdy)'
   idpdy = 0.

   allocate(SRelaxTau(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (SRelaxTau)'
   SRelaxTau = 0.

   allocate(TRelaxTau(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (TRelaxTau)'
   TRelaxTau = 0.

   allocate(uprof(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_observations: Error allocating (uprof)'
   uprof = 0.

   allocate(vprof(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_observations: Error allocating (vprof)'
   vprof = 0.

   allocate(epsprof(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_observations: Error allocating (epsprof)'
   epsprof = 0.

!  Setting relaxation time profiles TRelaxTau and SRelaxTau for T and S
!  SRelaxBott : height of bottom relaxation layer for S in meters
!  SRelaxSurf : height of surface relaxation layer for S in meters
!  SRelaxTauB : relaxation time in bottom relaxation layer for S in seconds
!  SRelaxTauS : relaxation time in surface relaxation layer for S in seconds
!  SRelaxTauM : relaxation time outside relaxation layers for S in seconds
!  TRelaxBott : height of bottom relaxation layer for T in meters
!  TRelaxTurf : height of surface relaxation layer for T in meters
!  TRelaxTauB : relaxation time in bottom relaxation layer for T in seconds
!  TRelaxTauS : relaxation time in surface relaxation layer for T in seconds
!  TRelaxTauM : relaxation time outside relaxation layers for T in seconds
   db=0.
   ds=depth
   do i=1,nlev
      TRelaxTau(i)=TRelaxTauM
      SRelaxTau(i)=SRelaxTauM
      db=db+0.5*h(i)
      ds=ds-0.5*h(i)
      if (db.le.SRelaxBott) SRelaxTau(i)=SRelaxTauB
      if (ds.le.SRelaxSurf) SRelaxTau(i)=SRelaxTauS
      if (db.le.TRelaxBott) TRelaxTau(i)=TRelaxTauB
      if (ds.le.TRelaxSurf) TRelaxTau(i)=TRelaxTauS
      db=db+0.5*h(i)
      ds=ds-0.5*h(i)
      if ((s_prof_method.ne.0).and.(SRelaxTau(i).le.0.)) then
         write(0,*) '       ', ''
         write(0,*) '       ', '***************************************************'
         write(0,*) '       ', 'SRelaxTau at i=',i,' is not a positive value.'
         write(0,*) '       ', 'Please correct namelist.inp and rerun.'
         write(0,*) '       ', 'Program aborted.'
         write(0,*) '       ', '***************************************************'
         stop 'init_observations'
      end if
      if ((t_prof_method.ne.0).and.(TRelaxTau(i).le.0.)) then
         write(0,*) '       ', ''
         write(0,*) '       ', '***************************************************'
         write(0,*) '       ', 'TRelaxTau at i=',i,' is not a positive value.'
         write(0,*) '       ', 'Please correct namelist.inp and rerun.'
         write(0,*) '       ', 'Program aborted.'
         write(0,*) '       ', '***************************************************'
         stop 'init_observations'
      end if
   end do

!  The salinity profile
   select case (s_prof_method)
      case (NOTHING)
         sprof = 0.
      case (ANALYTICAL)
         call analytical_profile(nlev,z,z_s1,s_1,z_s2,s_2,sprof)
      case (FROMFILE)
         open(s_prof_unit,file=s_prof_file,status='unknown',err=101)
         write(0,*) '       ', 'Reading salinity profiles from:'
         write(0,*) '           ', trim(s_prof_file)
!SP-this is were the initialization of the salt profile happens
         call get_s_profile(s_prof_unit,julday,secs,nlev,z)

!SP -this is used for a mean salt profile initialization, you need to also comment out the get_s_profile line above
!   OPEN(UNIT=72,FILE=init_s_prof,IOSTAT=ios)
!	IF(ios/=0) THEN
!		PRINT*,"Error during opening input file, stopping"; STOP
!       END IF
!	DO count=1,150                 		
!		READ (UNIT=72,FMT='(3X,E13.11)') sprof(count)
!	END DO             
!	CLOSE (UNIT=72)
!end SP

      case default
   end select

!  The temperature profile
   select case (t_prof_method)
      case (NOTHING)
         tprof = 0.
      case (ANALYTICAL)
         call analytical_profile(nlev,z,z_t1,t_1,z_t2,t_2,tprof)
      case (FROMFILE)
         open(t_prof_unit,file=t_prof_file,status='unknown',err=102)
         write(0,*) '       ', 'Reading temperature profiles from:'
         write(0,*) '           ', trim(t_prof_file)
!SP-this is were the initialization of the temp profile happens
         call get_t_profile(t_prof_unit,julday,secs,nlev,z)

!SP -this is used for a mean temperature profile initialization, you need to also comment out the get_t_profile line above
!   OPEN(UNIT=71,FILE=init_t_prof,IOSTAT=ios)
!	IF(ios/=0) THEN
!		PRINT*,"Error during opening input file, stopping"; STOP
!       	END IF
!	DO count=1,150                 		
!		READ (UNIT=71,FMT='(3X,E13.11)') tprof(count)
!	END DO             
!	CLOSE (UNIT=71)
!end SP
      
     case default
   end select

!  The external pressure
   select case (ext_press_method)
      case (FROMFILE)
         open(ext_press_unit,file=ext_press_file,status='unknown',err=103)
         write(0,*) '       ', 'Reading external pressure from:'
         write(0,*) '           ', trim(ext_press_file)
      case default
   end select
   call get_ext_pressure(ext_press_method,ext_press_unit,julday,secs)

!  The internal pressure
   select case (int_press_method)
      case (CONSTANT)
         dsdx=const_dsdx
         dsdy=const_dsdy
         dtdx=const_dtdx
         dtdy=const_dtdy
         idpdx=const_idpdx
         idpdy=const_idpdy
      case (FROMFILE)
         open(int_press_unit,file=int_press_file,status='unknown',err=104)
         write(0,*) '       ', 'Reading internal pressure from:'
         write(0,*) '           ', trim(int_press_file)
      case default
   end select
   call get_int_pressure(int_press_method,int_press_unit,julday,secs,nlev,z)

!  The light extinction profiles
   select case (extinct_method)
      case (1)
         A=0.58;g1=0.35;g2=23.0
      case (2)
         A=0.68;g1=1.20;g2=28.0
      case (3)
         A=0.62;g1=0.60;g2=20.0
      case (4)
         A=0.67;g1=1.00;g2=17.0
      case (5)
         A=0.77;g1=1.50;g2=14.0
      case (6)
         A=0.78;g1=1.40;g2=7.9
      case (7)
         A=0.7;g1=0.40;g2=8.0 ! Adolf Stips - Lago Maggiore

!!!!!!WT Case (0), (12), (13), (14) temporarily moved to case default.
      ! case (12)
      !    open(extinct_unit,file=extinct_file,status='unknown',err=105)
      !    write(0,*) '       ', 'Reading chlorophyll-a data from:'
      !    write(0,*) '           ', trim(extinct_file)
      !    call read_chlo(extinct_unit,julday,secs)
      !    !SP 22/02/06  The parameterisation holds for values of chlorophyll between 0.03 and 3.0
      !    IF(chlo.gt.3.0) THEN
      !       chlo=3.0
      !    END IF
      !    IF (chlo.lt.0.03) THEN
      !       chlo=0.03
      !    END IF

      ! case(13)  !HX 12/May/2017
      !    open(extinct_unit,file=extinct_file,status='unknown',err=105)
      !    write(1,*) '       ', 'Reading absorption-c data from:'
      !    write(1,*) '           ', trim(extinct_file)
      !    call read_IOP(extinct_unit,julday,secs)  
      !    open(extinct_unit,file=extinct_file,status='unknown',err=105)
      !    write(2,*) '       ', 'Reading backscattering-c data from:'
      !    write(2,*) '           ', trim(extinct_file)
      !    call read_IOP(extinct_unit,julday,secs)  
        
      !  case (14) !HX 16/June/2017
      !    open(extinct_unit,file=extinct_file,status='unknown',err=105)
      !    write(0,*) '       ', 'Reading chlorophyll-a data from:'
      !    write(0,*) '           ', trim(extinct_file)
      !    call read_chlo(extinct_unit,julday,secs)
         
       case default
         open(extinct_unit,file=extinct_file,status='unknown',err=105)
         write(0,*) '       ', 'Reading extinction data from:'
         write(0,*) '           ', trim(extinct_file)
         call read_extinction(extinct_unit,julday,secs,extinct_method)

         !WT 20170726, postprocessing code (e.g. imposes 0.03 .lt. chlo .lt. 3.0)
         ! reloacted to read_extinct.f90
         
   end select
   

!  The vertical advection velocity
   select case (w_adv_method)
      case (FROMFILE)
         open(w_adv_unit,file=w_adv_file,status='unknown',err=106)
         write(0,*) '       ', 'Reading vertical velocity observations from:'
         write(0,*) '           ', trim(w_adv_file)
      case default
   end select
   call get_w_adv(w_adv_method,w_adv_unit,julday,secs)

!  The sea surface elevation
   select case (zeta_method)
      case (FROMFILE)
         open(zeta_unit,file=zeta_file,status='unknown',err=107)
         write(0,*) '       ', 'Reading sea surface elevations from:'
         write(0,*) '           ', trim(zeta_file)
      case default
   end select
   call get_zeta(zeta_method,zeta_unit,julday,secs)

!  The observed velocity profile
   select case (vel_prof_method)
      case (0)
         uprof = 0.
         vprof = 0.
      case (2)
         open(vel_prof_unit,file=vel_prof_file,status='UNKNOWN',err=108)
         write(0,*) '       ', 'Reading velocity profiles from:'
         write(0,*) '           ', trim(vel_prof_file)
         call get_vel_profile(vel_prof_unit,julday,secs,nlev,z)
      case default
   end select

!  The observed dissipation profile
   select case (e_prof_method)
      case (0)
         epsprof = 0.
      case (2)
         open(e_prof_unit,file=e_prof_file,status='UNKNOWN',err=109)
         write(0,*) '       ', 'Reading dissipation profiles from:'
         write(0,*) '           ', trim(e_prof_file)
         call get_eps_profile(e_prof_unit,julday,secs,nlev,z)
      case default
   end select

   return

80 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(fn),'" for reading'
   stop 'init_observations'
81 write(0,*) 'FATAL ERROR: ', 'I could not read "sprofile" namelist'
   stop 'init_observations'
82 write(0,*) 'FATAL ERROR: ', 'I could not read "tprofile" namelist'
   stop 'init_observations'
83 write(0,*) 'FATAL ERROR: ', 'I could not read "ext_pressure" namelist'
   stop 'init_observations'
84 write(0,*) 'FATAL ERROR: ', 'I could not read "int_pressure" namelist'
   stop 'init_observations'
85 write(0,*) 'FATAL ERROR: ', 'I could not read "extinct" namelist'
   stop 'init_observations'
86 write(0,*) 'FATAL ERROR: ', 'I could not read "w_advspec" namelist'
   stop 'init_observations'
87 write(0,*) 'FATAL ERROR: ', 'I could not read "zetaspec" namelist'
   stop 'init_observations'
88 write(0,*) 'FATAL ERROR: ', 'I could not read "velprofile" namelist'
   stop 'init_observations'
89 write(0,*) 'FATAL ERROR: ', 'I could not read "eprofile" namelist'
   stop 'init_observations'
90 write(0,*) 'FATAL ERROR: ', 'I could not read "bprofile" namelist'
   stop 'init_observations'

101 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(s_prof_file),'" for reading'
   stop 'init_observations'
102 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(t_prof_file),'" for reading'
   stop 'init_observations'
103 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(ext_press_file),'" for reading'
   stop 'init_observations'
104 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(int_press_file),'" for reading'
   stop 'init_observations'
105 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(extinct_file),'" for reading'
   stop 'init_observations'
106 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(w_adv_file),'" for reading'
   stop 'init_observations'
107 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(zeta_file),'" for reading'
   stop 'init_observations'
108 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(vel_prof_file),'" for reading'
   stop 'init_observations'
109 write(0,*) 'FATAL ERROR: ', 'Unable to open "',trim(e_prof_file),'" for reading'
   stop 'init_observations'

   return
   end subroutine init_observations 
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_all_obs
!
! !INTERFACE:
   subroutine get_all_obs(julday,secs,nlev,z)
!
! !DESCRIPTION:
!  During the time integration this subroutine is called each time step
!  to update all 'observation'. The routine is basically a wrapper
!  routine which calls the variable specific routines.
!  The only input to this routine is the time (in internal GOTM 
!  representation) and the vertical grid. It is up to each of the individual
!  routines to use this information and to provide updated 'observations'.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: julday,secs
   integer, intent(in)                 :: nlev
   double precision, intent(in)                :: z(:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if(s_prof_method .eq. 2) then
      call get_s_profile(s_prof_unit,julday,secs,nlev,z)
   end if

   if(t_prof_method .eq. 2) then
      call get_t_profile(t_prof_unit,julday,secs,nlev,z)
   end if

   call get_ext_pressure(ext_press_method,ext_press_unit,julday,secs)

   call get_int_pressure(int_press_method,int_press_unit,julday,secs,nlev,z)

   if(extinct_method .eq. 0) then
      call read_extinction(extinct_unit,julday,secs)
   end if
    if(extinct_method .eq. 12) then
        call read_chlo(extinct_unit,julday,secs)
    end if
    if(extinct_method .eq. 13) then
        call read_IOP(extinct_unit,julday,secs)
    end if
    if(extinct_method .eq. 14) then
        call read_chlo(extinct_unit,julday,secs)
    end if
   call get_w_adv(w_adv_method,w_adv_unit,julday,secs)

   call get_zeta(zeta_method,zeta_unit,julday,secs)

   if(vel_prof_method .eq. 2) then
      call get_vel_profile(vel_prof_unit,julday,secs,nlev,z)
   end if

   if(e_prof_method .eq. 2) then
      call get_eps_profile(e_prof_unit,julday,secs,nlev,z)
   end if

   return
   end subroutine get_all_obs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_obs
!
! !INTERFACE:
   subroutine read_obs(unit,yy,mm,dd,hh,min,ss,N,obs,ierr)
!
! !DESCRIPTION:
!  This routine will read all non-profile observations.
!  The routine allows for reading more than one scalar variable at a time.
!  The number of data to be read is specified by {\tt N}.
!  Data read-in are returned
!  in the 'obs' array. It is up to the calling routine to assign
!  meaning full variables to the individual elements in {\tt obs}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: N
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: yy,mm,dd,hh,min,ss
   double precision,intent(out)                :: obs(:)
   integer, intent(out)                :: ierr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ierr=0
   read(unit,'(A82)',ERR=100,END=110) cbuf
   read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   read(cbuf(20:),*,ERR=100,END=110) (obs(i),i=1,N)
   return
100 ierr=READ_ERROR 
   return
110 ierr=END_OF_FILE 
   return
900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)
   end subroutine read_obs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_profiles
!
! !INTERFACE:
   subroutine read_profiles(unit,nlev,cols,yy,mm,dd,hh,min,ss,z, &
                            profiles,lines,ierr)
!
! !DESCRIPTION:
!  Similar to {\tt read\_obs()} but used for reading profiles instead of 
!  scalar data.
!  The data will be interpolated on the grid specified by nlev and z.
!  The data can be read 'from the top' or 'from the bottom' depending on
!  a flag in the actual file. 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: nlev,cols
   double precision, intent(in)                :: z(:)
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)              :: lines
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: yy,mm,dd,hh,min,ss
   double precision, intent(out)               :: profiles(:,:)
   integer, intent(out)                :: ierr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: i,j,rc
   integer                   :: N,up_down
   double precision,allocatable,dimension(:)   :: tmp_depth
   double precision,allocatable,dimension(:,:) :: tmp_profs
!
!-----------------------------------------------------------------------
!BOC
   ierr=0
   read(unit,'(A72)',ERR=100,END=110) cbuf
   read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   read(cbuf(20:),*,ERR=100,END=110) N,up_down

   lines = lines+1

   allocate(tmp_depth(0:N),stat=rc)
   if (rc /= 0) stop 'read_profiles: Error allocating memory (tmp_depth)'
   allocate(tmp_profs(0:N,cols),stat=rc)
   if (rc /= 0) stop 'read_profiles: Error allocating memory (tmp_profs)'

   if(up_down .eq. 1) then
      do i=1,N
         lines = lines+1
         read(unit,*,ERR=100,END=110) tmp_depth(i),(tmp_profs(i,j),j=1,cols)
      end do
   else
      do i=N,1,-1 
         lines = lines+1
         read(unit,*,ERR=100,END=110) tmp_depth(i),(tmp_profs(i,j),j=1,cols)
      end do
   end if
!   PRINT*,tmp_depth
!   PRINT*,tmp_profs
!   READ*
   call gridinterpol(N,cols,tmp_depth,tmp_profs,nlev,z,profiles)
!   PRINT*,tmp_profs
!   PRINT*,profiles
   deallocate(tmp_depth)
   deallocate(tmp_profs)

   return
100 ierr=READ_ERROR 
   return
110 ierr=END_OF_FILE 
   return
900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)

   end subroutine read_profiles
!EOC

!-----------------------------------------------------------------------

   end module observations

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
