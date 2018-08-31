!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The light absorbtion
!
! !INTERFACE:
   subroutine light_absorbtion(nlev,I_0) 
!
! !DESCRIPTION:
!  The irradiance profile is used in the temperature equation as an
!  inner source term (however the sum of latent, sensible, and longwave
!  radiation is treated as a boundary condition) and used in the biological
!  equation to compute the depth-profile of PAR (Photosyntehic Active
!  Radiation). Absorbtion of radiation at interface levels is calculated
!  following various exponential law:
!
!  \begin{equation}
!  \begin{array}{l}
!   Rad(z)=Qlw+Qsw \\ 
!   Qlw=I_0*A*e^{-K1.zint} \\
!   Qsw=I_0*(1-A)*e^{-(K2+K3*Phy_av).zint}, \mbox{where } Phy_av=Phyt/zint
!   \end{array}
!  \end{equation}
!
!  where "A"=the weighting function for spectral range
!  and the extinction coefficients read:
!
!        \begin{tabular}{ll}
!        "K1"=for long-wave radiation  "Qlw" &   (red) [/m] \\
!        "K2"=for short-wave radiation "Qsw" &   (visible blue-green) [/m] \\
!        "K3"=for biotic self-shading substance & (Cholrophyll a) [m2/mmolN]
!        \end{tabular}
!
! !USES:
   use meanflow, only: h
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	                :: nlev
   double precision, intent(in)	                :: I_0
!
! !REVISION HISTORY:
!  Original author(s): Pierre-Phillipe Mathieu 
!
!  $Log: lightabsorbtion.F90,v $
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 08:50:06  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
! !LOCAL VARIABLES:
   logical, save             :: kbk_dummy = .true.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (kbk_dummy ) then
!kbk      STDERR 'light_absorption is not finished - we wait for PP'
      kbk_dummy = .false. 
   end if



   return
   end subroutine light_absorbtion 
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
