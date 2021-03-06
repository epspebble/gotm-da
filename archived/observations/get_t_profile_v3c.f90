!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_t_profile
!
! !INTERFACE:
   subroutine get_t_profile(unit,jul,secs,nlev,z)
!
! !DESCRIPTION:
!  This routine is responsible for providing sane values to an `observed'
!  temperature profile.
!  The subroutine is called in the {\tt get\_all\_obs()} subroutine
!  as part of the main integration loop.
!  In case of observations from file the temporal interpolation is
!  done in this routine.
!
! !USES:
   use time
   use observations, only: read_profiles, tprof, next_tprof_julianday, next_tprof_secondsofday

   ! WT 20171112 Temporary. Make it work first. Should update this after assimilation related code moved to new module.
   ! Probably not a good idea to ask get_t_profile to be 'aware' of assimilation events. Maybe just pass an extra flag
   ! about interpolation preferences somewhere.
   use gotm, only: assim_window ! WT 

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: jul,secs
   integer, intent(in)                 :: nlev

   double precision, intent(in)                :: z(0:nlev)

   ! !OUTPUT PARAMETERS:
   
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_t_profile.F90,v $
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 08:51:57  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: yy,mm,dd,hh,min,ss
   double precision                  :: t,dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   integer, parameter        :: cols=1
   integer, save             :: lines=0
   integer, save             :: nprofiles=0
   logical, save             :: one_profile=.false.
   double precision, save, dimension(:,:), allocatable :: prof1,prof2,alpha

   ! WT 20171112 For instructing the subroutine NOT to interpolate in time.
   ! This is used in the case when daily assimilation occurs with assim_window .eq. 2,
   ! where the timestamps are the intended time for assimilation, and so the profile
   ! should be read in exactly without interpolation.
   logical                   :: interpolate=.true.
   if (assim_window .eq. 2) then
      interpolate = .false.
   end if

!
!-----------------------------------------------------------------------
!BOC
!SP - the first 3 if statements allocate memory during the initialisation
   if ( .not. allocated(prof1)) then
      allocate(prof1(0:nlev,cols),stat=rc)
      if (rc /= 0) stop 'get_t_profile: Error allocating memory (prof1)'
      prof1 = 0.
   end if
   if ( .not. allocated(prof2)) then
      allocate(prof2(0:nlev,cols),stat=rc)
      if (rc /= 0) stop 'get_t_profile: Error allocating memory (prof2)'
      prof2 = 0.
   end if
   if ( .not. allocated(alpha)) then
      allocate(alpha(0:nlev,cols),stat=rc)
      if (rc /= 0) stop 'get_t_profile: Error allocating memory (alpha)'
   end if

!  This part initialise and read in new values if necessary
   if(.not. one_profile .and. time_diff(jul2,secs2,jul,secs) .lt. 0) then 
      do
         jul1 = jul2
         secs1 = secs2
         prof1 = prof2
         call read_profiles(unit,nlev,cols,yy,mm,dd,hh,min,ss,z,prof2,lines,rc)
         if(rc .ne. 0) then
            if(nprofiles .eq. 1) then
               write(0,*) '           ', 'Only one temperature profile present.'
               one_profile = .true.
               tprof = prof1(:,1)
            else
               write(0,*) 'FATAL ERROR: ', 'Error reading temperature profile around line #',lines
            end if    
            EXIT
         else
            !WT When reading multiple records from file, obtain next record time and find time difference.
            nprofiles = nprofiles + 1  
            call julian_day(yy,mm,dd,jul2)
            secs2 = hh*3600 + min*60 + ss

            !WT 20171112 Communicate the time of the next record for the assimilation routines to be aware of.
            !print *, 'get_t_profile() called'
            !print *, 'prev:',jul1,secs1,'next:',jul2,secs2
            next_tprof_julianday = jul2
            next_tprof_secondsofday = secs2
     
            !WT If now (jul,secs) is later than the next record time, no more new information, exit DO.
            if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT            
         end if
      end do
      if( .not. one_profile) then
         if (interpolate) then
            dt = time_diff(jul2,secs2,jul1,secs1)
            alpha = (prof2-prof1)/dt
         else
            !WT 20171112 When assim_window .eq. 2, this is when the tprof is updated and ready for assimilation.
            tprof = prof2(:,1)
         end if
      end if
   end if

!  Do the time interpolation
   if (( .not. one_profile) .and. interpolate) then
      t  = time_diff(jul,secs,jul1,secs1)
      tprof = prof1(:,1) + t*alpha(:,1)     
   end if

   return
   end subroutine get_t_profile
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
