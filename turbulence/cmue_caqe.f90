!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  \cite{Canutoetal2001a} quasi-eq.\ stability func.\ ver.\ A
! 
! !INTERFACE:
   subroutine cmue_caqe(nlev)
!
! !DESCRIPTION:
!  This subroutine computes quasi-equilibrium version of the stability functions
!  according to \cite{Canutoetal2001a}, see section \ref{cmue_ca}.
!  The version `A' is implemented. 
!  Here, the
!  quasi-equilibrium assumption (\ref{quasieq}) is applied in addition.
!
! !USES:
   use turbulence, only: an
   use turbulence, only: cmue1,cmue2,cm0
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY: 
!  Original author(s): Hans Burchard & Karsten BoldingQ
!
!  $Log: cmue_caqe.F90,v $
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 09:02:03  gotm
!  Added new Generic Turbulence Model + improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
! 
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   double precision, parameter       :: L1=0.1070
   double precision, parameter       :: L2=0.0032
   double precision, parameter       :: L3=0.0864
   double precision, parameter       :: L4=0.1200
   double precision, parameter       :: L5=11.9000
   double precision, parameter       :: L6=0.4000
   double precision, parameter       :: L7=0.0000
   double precision, parameter       :: L8=0.4800
   double precision, parameter       :: tnmin=-12.27
   double precision                  :: a2_cm03
   double precision                  :: s0,s1,s2,s4,s5,s6
   double precision                  :: d,d0,d1,d2,d3,d4,d5
   double precision                  :: tsmax
   double precision                  :: tn,ts,sm,sh
   double precision                  :: PP,QQ
!
!-----------------------------------------------------------------------
!BOC
   a2_cm03=2./cm0**3
   s0 = 1.5*L1*L5*L5
   s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
   s2 = -0.375*L1*(L6*L6-L7*L7)
   s4 = 2.*L5
   s5 = 2.*L4
   s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7)

   d0 = 3.*L5*L5
   d1 = L5*(7.*L4+3.*L8)
   d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
   d3 = L4*(4.*L4+3.*L8)
   d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
   d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

   do i=1,nlev-1
      tn = 4./cm0**6 * an(i)
      if (tn .lt. tnmin) tn = tnmin

      PP=(s0+(s1-s6)*tn-2.*(d2+d4*tn))/(s2-2.*d5)
      QQ=-(2.*(d0+d1*tn+d3*tn*tn)+(s4+s5*tn)*tn)/(s2-2.*d5)

      ts=-0.5*PP-sqrt(PP**2/4.-QQ)

      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1(i) = a2_cm03 * sm
      cmue2(i) = a2_cm03 * sh
   end do

   return
   end subroutine cmue_caqe
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
