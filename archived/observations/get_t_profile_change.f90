!SP 01/12/04 changes made
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
   use observations, only: read_profiles,tprof
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: jul,secs
   integer, intent(in)                 :: nlev
   double precision, intent(in)                :: z(0:nlev)
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

!if loop added by SP 01-12-04
if(jul==2449643.and.secs==30)then
   jul2=2449643
   secs2=900
   dd=17
   lines=0
   close(unit=31)
   open(unit=31,file="OBS/ARABIAN/arabian_t.dat")
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
            nprofiles = nprofiles + 1
            call julian_day(yy,mm,dd,jul2)
            secs2 = hh*3600 + min*60 + ss
            if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
         end if
      end do
      if( .not. one_profile) then
         dt = time_diff(jul2,secs2,jul1,secs1)
         alpha = (prof2-prof1)/dt
      end if
   end if

!  Do the time interpolation
   if( .not. one_profile) then
      t  = time_diff(jul,secs,jul1,secs1)
      tprof = prof1(:,1) + t*alpha(:,1)
     
   end if


   return
   end subroutine get_t_profile
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
