!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_extinction
!
! !INTERFACE:
   subroutine read_extinction(unit,jul,secs,extinct_method)
!
! !DESCRIPTION:
!  This routine will provide the light extinction coefficients. It
!  is only called if no Jerlov class has been specified in {\tt obs.inp}.
!
! !USES:
   use time
   use observations, only : read_obs
   use observations, only : A,g1,g2

   !WT 20170727 for case(CHLOR_A)
   use observations, only : chlo
   !WT 20170727 for case(IOP_A_BB)
   use observations, only : abp_coe, bb
   
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit,jul,secs
   integer, optional, intent(in)       :: extinct_method
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: read_extinction.F90,v $
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 09:02:09  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:51:58  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   double precision                  :: t
   double precision, save            :: dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   ! WT 20170726 Nothing changed. Just note that the default max number of
   ! columns to read in is 3. Need to change the following if more needed.
   double precision, save            :: alpha(3)
   double precision, save            :: obs(3),obs1(3),obs2(3)=0.

   ! WT 20170726 For more flexibility and readability.
   integer                   :: n ! # of obs entries (columns) that we'll really use.
   integer, parameter        :: CHLOR_A=12
   integer, parameter        :: IOP_A_BB=13
   integer, parameter        :: WIP=14 ! WIP = 'Work In Progress'
   
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialise and read in new values if necessary.

   !WT Prior call signature does not have 'extinct_method' as input,
   ! Setting a default value that matches that in observations.f90

   !if (.not. present(extinct_method)) then
   !   extinct_method = 0
   !end if
   
   !WT Determine real size of obs, i.e. the no. of columns to read from file.
   select case (extinct_method)
   case (CHLOR_A)
      n = 1
   case (IOP_A_BB)
      n = 2
   case (WIP)
      STOP 'An incomplete work in progress in read_extinction.f90'
   case default !WT The original value in Revision 1.4 by kbk is 0.
      n = 3
   end select

   !WT Read from file only if the simulation time (jul, secs) is later than the
   ! previous observation time (jul2, secs2), and find linear interpolation
   ! coefficient(s), alpha.
   if (time_diff(jul2,secs2,jul,secs) .lt. 0) then 
      do
         jul1 = jul2
         secs1 = secs2
         obs1 = obs2
         call read_obs(unit,yy,mm,dd,hh,min,ss,n,obs2,rc)
         call julian_day(yy,mm,dd,jul2)
         secs2 = hh*3600 + min*60 + ss
         if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(jul2,secs2,jul1,secs1)

      !WT Be reminded this is vector arithmetic.
      alpha = (obs2-obs1)/dt
   end if

   !WT Find the length of time from previous and do the time interpolation
   t  = time_diff(jul,secs,jul1,secs1)
   obs = obs1 + t*alpha


   !WT Now write to the appropriate, public variables, depending on
   ! extinct_method.

   select case (extinct_method)
      
   case (CHLOR_A) !WT New case by SP, 20170106
      chlo = obs(1)

      !WT 20170726 "postprocessing code" due to SP 20160222
      if (chlo .gt. 3.0) then
         chlo = 3.0
      else if (chlo .lt. 0.03) then
         chlo = 0.03
      end if
            
   case (IOP_A_BB) !WT New case by HX, 20170512
      abp_coe = obs(1)
      bb = obs(2)
      
   case (WIP) !WT 201707026 A new work in progress by HX, 20170616
      ! Expect something new here. Update case (WIP) above as well.
      
   case default !WT Original code in Rev1.4 (kbk) follows, which is case 0.
      A = obs(1)
      g1 = obs(2)
      g2 = obs(3)

   end select
     
   return
   end subroutine read_extinction

   !EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
