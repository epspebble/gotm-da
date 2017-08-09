!-----------------------------------------------------------------------
!BOP
!
!
! !ROUTINE: The temperature equation \label{sec:temperature}
!
! !INTERFACE:
   
  
!------------------------------------------------------------
    
!------------------------------------------------------------
    subroutine temperature(nlev,dt,cnpar,I_0,heat,advect,qdir_frac,qdiff_frac,cosr,nuh,rad)
!
!
! !DESCRIPTION:
! This subroutine computes the balance of heat in the form
!  \begin{equation}
!   \label{TEq}
!    \dot{\theta}
!    = {\cal D}_\theta
!    - \frac{1}{\tau_R(\theta)}(\theta-\theta_{obs})
!    + \frac{1}{C_p \rho_0} \partder{I}{z}
!    \comma
!  \end{equation}
!  where $\dot{\theta}$ denotes the material derivative of the potential 
!  temperature $\theta$, and
!  ${\cal D}_\theta$ is the sum of the turbulent and viscous transport
!  terms modelled according to
!  \begin{equation}
!   \label{DT}
!    {\cal D}_\theta 
!    = \frstder{z} 
!     \left( 
!        \left( \nu'_t + \nu^\theta \right) \partder{\theta}{z}
!      \right) 
!    \point
!  \end{equation}
!  In this equation, $\nu'_t$ and $\nu^\theta$ are the turbulent and 
!  molecular diffusivities of heat, respectively. The computation
!  of $\nu'_t$ is discussed in \sect{sec:turbulenceIntro}.
!
!  Horizontal advection is optionally
!  included  (see {\tt obs.inp}) by means of prescribed
!  horizontal gradients $\partial_x\theta$ and $\partial_y\theta$ and 
!  calculated horizontal velocities $u$ and $v$.
!  Relaxation with the time scale $\tau_R (\theta)$ 
!  towards a precribed (changing in time)
!  profile $\theta_{obs}$ is possible. 
!
!  The sum of latent, sensible, and longwave radiation is treated
!  as a boundary condition. Solar radiation is treated as an inner 
!  source, $I(z)$. It is computed according the
!  exponential law (see \cite{PaulsonSimpson77})
!  \begin{equation}
!    \label{Iz}
!    I(z) = I_0 \bigg(Ae^{-\eta_1z}+(1-A)e^{-\eta_2z}\bigg).
!  \end{equation}
!  The absorbtion coefficients $\eta_1$ and $\eta_2$ depend on the water type
!  and have to be prescribed either by means of choosing a \cite{Jerlov68} class
!  (see \cite{PaulsonSimpson77}) or by reading in a file through the namelist
!  {\tt extinct} in {\tt obs.inp}. 

!  Diffusion is numerically treated implicitly, see equations (\ref{sigmafirst})-
!  (\ref{sigmalast}).
!  The tri--diagonal matrix is solved then by a simplified Gauss elimination.
!  Vertical advection is included for accounting for adaptive grids,
!  see {\tt adaptivegrid.F90}.
!
! !USES:
   use meanflow, only: avmolt,rho_0,cp
   use meanflow, only: h,ho,u,v,T,avh,w,grid_method,w_grid
   use observations, only: dtdx,dtdy,t_adv,w_adv,w_adv_discr
   use observations, only: tprof,TRelaxTau,w_adv_method
   use airsea, only:cloud,coszen
! Sh 03/12/2003
   use observations, only: A,g1,g2,chlo
   use observations, only: extinct_method
   use observations, only: fsen,zdeta
! enable this routine to know local time - for relaxation purposes
!   use time, only: local_secondsofday

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   double precision, intent(in)        :: dt,cnpar
   double precision, intent(in)        :: I_0,heat,cosr,advect(1:nlev)
!HK cosr is the angle of incidence of the solar radiation
   double precision, intent(in)        :: qdir_frac,qdiff_frac
   double precision, intent(in)        :: nuh(0:nlev)
!
! !OUTPUT PARAMETERS:
   double precision       :: rad(0:nlev)
!HK/SH added Qsour here for output to ncdf for output 22/10/2003 @ 15:00
   double precision       :: Q_source(0:nlev)
!SH - 03/12/2003 - added tprofdum to assimilate tprof and smootprof below TRelaxSmoo
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: temperature_tested.f90,v $
!Revision 1.2  2003/12/03  17:09:02  stephenh
!Relax at night and/or below some level with variable RelaxTau
!
!Revision 1.1  2003/10/24  10:24:54  helen
!Initial revision
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:50:07  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/11/18 11:50:37  gotm
!  Cleaned
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
! Added mods to distribute insolation through the model depth in 9 wavebands
! Includes diffuse and direct components
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,j,Bcup,Bcdw,flag
! SH removed Qsource from local since now public for output 22/10/2003 @ 15:00
!   double precision          :: Qsour(0:nlev) 
   double precision          :: Tup,Tdw,z
   logical                   :: surf_flux,bott_flux
   double precision		:: absfac_dir=0.,absfac_diff=0.
   double precision		:: trans,para_A,para_K
   double precision                  :: C1(16) = &
                 (/.026,-.009,-.015,-.003,.063,.278,3.91,16.64, &
                   .033,-.01,-.019,-.006,.066,.396,7.68,51.27/)
   double precision                  :: C2(16) = &
                 (/.112,.034,-.006,-.131,-.015,-.562,-12.91,-478.28, &
                   .0,.0,.0,.0,.0,.0,.0,.0/)
   double precision                  :: C3(16) = &
                 (/.0,.0,.0,.0,.0,.0,.0,.0, &
                   -.025,-.007,-.003,-.004,.006,-.027,-2.49,13.14/)
   double precision                  :: C4(16) = &
                 (/.366,.207,.188,.169,.082,1.02,16.62,736.56, &
                   .419,.231,.195,.154,.066,.886,17.81,665.19/)
!
!-----------------------------------------------------------------------
!BOC

!  hard coding of parameters, to be included into namelist for gotm2.0 
   Bcup=1 !BC Neumann
   Tup=-heat/(rho_0*cp)!Heat flux (positive upward)
   Bcdw=1 !BC Neumann
   Tdw=0.!No flux
   surf_flux=.false.                     
   bott_flux=.false.

! could argue here that adjustment must be made for delta??

   rad(nlev)=I_0/(rho_0*cp)
   z=0.
   
! add delta to z??
   do i=nlev-1,0,-1 
      
      z=z+h(i+1)
      
      select case (extinct_method)
         

         case (8)              
! usually use this case
! 9 stream solar radiation distribution accounting for direct and 
! diffuse (5/3 approx - Andrews(book)) components and solar zenith angle 
! and refraction

            if (cosr > 0.0) then

               do j=1,9
                  absfac_dir=absfac_dir+fsen(j)*exp(-z/(zdeta(j)*cosr))
! Is this right? how do we decouple diffuse from direct?
                  Absfac_diff=absfac_diff+fsen(j)*exp(-z*(5.0/3.0)/(zdeta(j)))
               end do

            else
            absfac_dir = 0.
            absfac_diff = 0.
            end if

          rad(i)=I_0/(rho_0*cp)*(absfac_dir*qdir_frac + absfac_diff*qdiff_frac)
          absfac_dir=0.
          absfac_diff=0.
         
        case (9)   
! change fewer things - i.e. leave out cosr correction.

            do j=1,9
               absfac_dir=absfac_dir+fsen(j)*exp(-z/(zdeta(j)))
!              absfac_diff=absfac_diff+fsen(j)*exp(-z/(zdeta(j)*(5./3.)))
            end do
          rad(i)=I_0*absfac_dir/(rho_0*cp)
          absfac_dir=0.
      
        case (10)  
!just missing cosr

               do j=1,9
                  absfac_dir=absfac_dir+fsen(j)*exp(-z/(zdeta(j)))
! Is this right? how do we decouple diffuse from direct?
                  absfac_diff=absfac_diff+fsen(j)*exp(-z*(5.0/3.0)/(zdeta(j)))
               end do
              
         rad(i)=I_0/(rho_0*cp)*(absfac_dir*qdir_frac + absfac_diff*qdiff_frac)
         absfac_dir=0.
         absfac_diff=0.

        case (11)   

! case 11 may be more realistic vers. of case 8. May make little diff however?
! more realistic: get some diffuse when cosr<0
               do j=1,9
  
               if (cosr > 0.0) then
                  absfac_dir=absfac_dir+fsen(j)*exp(-z/(zdeta(j)*cosr))
               else
                  absfac_dir = 0.
               end if

! Is this right? how do we decouple diffuse from direct?
               absfac_diff=absfac_diff+fsen(j)*exp(-z*(5.0/3.0)/(zdeta(j)))
               end do

          rad(i)=I_0/(rho_0*cp)*(absfac_dir*qdir_frac + absfac_diff*qdiff_frac)
          absfac_dir=0.
          absfac_diff=0.


          case(12)

! case 12 is the solar radiation transmission through the upper ocean based on the parameterisation developed from Ohlmann and Siegel, J. Phys. Oceangr. 2000
! it uses satellite derived chlorophyll cencentrations, cloud amount and solar zenith angle.
     IF (I_0.le.0) then 
        rad(i)=0
     ElSE      
             trans=0.0
             IF(cloud.gt.0.1) then
             !revised SP 20/06/17 setting clear sky limit as 0.1 following Table 1 in Ohlmann&Siegel(2000)
                DO j=1,4
                   para_A=C1(j)*chlo+C2(j)*cloud+C4(j)
                   para_K=C1(j+4)*chlo+C2(j+4)*cloud+C4(j+4)
                   trans=trans+para_A*exp(-para_K*z)
                END DO
             ELSE
                DO j=9,12
                   IF(coszen.lt.0.2588) then
            !for very low value of coszen the parameterisation breaks down!  SP 22/02/06 revised 28/02/17
            !revised again SP 20/06/17 setting the limit as theta=75degrees (cos(theta)=0.2588) the largest angle in Ohlmann&Siegel(2000)
                      para_A=C1(j)*chlo+(C3(j)/0.2588)+C4(j)
                      para_K=C1(j+4)*chlo+(C3(j+4)/0.2588)+C4(j+4)
                      trans=trans+para_A*exp(-para_K*z)
                   ELSE                    
                      para_A=C1(j)*chlo+(C3(j)/coszen)+C4(j)
                      para_K=C1(j+4)*chlo+(C3(j+4)/coszen)+C4(j+4)
                      trans=trans+para_A*exp(-para_K*z)		
                   END IF
                END DO
             END IF
             ! 'trans' includes albedo implicitly due to Ohlmann, check whether I_0 includes albedo as well
             rad(i)=I_0*trans/(rho_0*cp)
      END IF    
       case default 
! extinct_methods 1-7 - these are defined in observations.f90
! leave out diffuse / direct discrimination
 
       rad(i)=I_0*(A*exp(-z/g1)+(1-A)*exp(-z/g2))/(rho_0*cp)
!PRINT*,A*exp(-z/g1)+(1-A)*exp(-z/g2)
   end select
   avh(i)=nuh(i)+avmolT 
   end do

   do i=1,nlev
      Q_source(i)=(rad(i)-rad(i-1))/h(i)           !SP 16/05/05
      !include advection source
      Q_source(i)=Q_source(i)+advect(i)                    !SP 08/05/05
      if (t_adv) Q_source(i)=Q_source(i)-u(i)*dtdx(i)-v(i)*dtdy(i) 
   end do

   flag=1  ! divergence correction for vertical advection

!Null hypothesis
!   avh=0.0
!   Q_source=0.0
!   Tup=0.0
!end Null hypothesis  

      call Yevol(nlev,Bcup,Bcdw,dt,cnpar,Tup,Tdw,TRelaxTau,h,ho,avh,w,        &
              Q_source,tprof,w_adv_method,w_adv_discr,T,surf_flux,bott_flux,  &
              grid_method,w_grid,flag)      

   return
   end subroutine temperature
!EOC
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
