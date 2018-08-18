!-----------------------------------------------------------------------
!BOP
!
!
! !ROUTINE: The temperature equation \label{sec:temperature}
!
! !INTERFACE:


!------------------------------------------------------------


!WT Part of the 'meanflow' module.

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
  use airsea, only:cloud,coszen,albedo
  ! Sh 03/12/2003
  use observations, only: A,g1,g2,chlo,abp_coe,bb
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
  !WT Q_source not an output of this subroutine, but is used when calling yevol()
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
  double precision          :: trans ! From function
  ! SH removed Qsource from local since now public for output 22/10/2003 @ 15:00
  !   double precision          :: Qsour(0:nlev) 
  double precision          :: Tup,Tdw,z
  logical                   :: surf_flux,bott_flux
  double precision		:: absfac_dir=0.,absfac_diff=0.
  
  !
  !-----------------------------------------------------------------------
  !BOC

  !  hard coding of parameters, to be included into namelist for gotm2.0 
  Bcup=1 !BC Neumann
  Tup=-heat/(rho_0*cp)!Heat flux (positive upward)
  !Bcdw=1 !BC Neumann
  Bcdw=2 !BC Dirichlet 
  Tdw=0.!No flux
  surf_flux=.false.                     
  bott_flux=.false.

  !-----------------------------------------------------------------------
  !WT Compute rad(nlev) for the top level, taking into account the reduction due to albedo.
  
  !WT What does 'delta' mean in the previous comments below?
  ! Is it the depth of top level? Is it z(nlev) (i.e. h(nlev))?  
  !
  !PREVIOUS COMMENTS:
  ! could argue here that adjustment must be made for delta??
  ! add delta to z??
  z=0. ! WT surmise this is depth of the top layer, sea surface.
  select case (extinct_method)
  case (12)
     !WT Surface radiation in Ohlmann-Siegel (2000)'s formulation is recovered by
     ! setting z=0 in the trasmission coefficient, and it 'agrees well with Payne
     ! (1972)'s value' (p.1859, last paragraph in section 5.)
     
     !WT Use below if I_0 has albedo factored out. 
     rad(nlev)= I_0/(1-albedo)/(rho_0*cp) * trans(z,extinct_method) !WT Last factor equals Trans1 in HX's version.
     !WT Use below if I_0 does not take albedo into account 
     ! rad(nlev)= I_0/(rho_0*cp)*OS_trans(0,extinct_method)

  !WT Does case (13), (14) use the Payne albedo? Should we modify these?
     
  case default
     !WT 2017-10-19 In the original code, default value of 'rad(nlev)'
     ! is assumed to be the following, REGARDLESS of 'extinct_method'     
     ! The term I_0 already has albedo effect considered (using Payne (1972)'s values).
     rad(nlev) = I_0/(rho_0*cp)
  end select

  !WT Now compute rad(i) for every level below the top.
  do i=nlev-1,0,-1
     ! WT surmise that z is the depth of the top of the layer in concern,
     ! and that both h and z are positive here, unlike the z in the output.
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
           ! absfac_diff=absfac_diff+fsen(j)*exp(-z/(zdeta(j)*(5./3.)))
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
        ! case 12 is the solar radiation transmission through the upper ocean
        ! based on the parameterisation developed from
        ! Ohlmann and Siegel, J. Phys. Oceangr. 2000
        ! it uses satellite derived chlorophyll cencentrations, cloud amount
        ! and solar zenith angle.

        IF (I_0.le.0) then 
           rad(i)=0
        ElSE
           !WT Read comments for the top layer.
           rad(i) = I_0/(1-albedo)*trans(z,extinct_method)/(rho_0*cp)
        END IF
        
     case (13) !HX 12/05/2017
        if (I_0 .le. 0) then
           rad(i) = 0
        else
           rad(i) = I_0*trans(z,extinct_method)/(rho_0*cp)                    
        end if

     case (14) !HX 16/06/2017
        if (I_0 .le. 0) then
           rad(i) = 0
        else
           rad(i) = I_0*trans(z,extinct_method)/(rho_0*cp)                    
        end if

     case (15)
        ! same as case (9) [Paulson & Simpson (1981)], except change the first two absorption coefficients based on data from Verevochkin & Startsev, J. Fluid Mech. (2005).
        zdeta(1)=24.1;  ! Jerlov water type I
        zdeta(2)=0.673; ! Jerlov water type I
        !zdeta(1)=11;    ! Jerlov water type II
        !zdeta(2)=0.664; ! Jerlov water type II
        !zdeta(1)=6.54;  ! Jerlov water type III
        !zdeta(2)=0.654; ! Jerlov water type III
        !zdeta(1)=3.36;  ! Jerlov water type 1
        !zdeta(2)=0.653; ! Jerlov water type 1
        !zdeta(1)=2.27;  ! Jerlov water type 3
        !zdeta(2)=0.645; ! Jerlov water type 3
        !zdeta(1)=1.56;  ! Jerlov water type 5
        !zdeta(2)=0.627; ! Jerlov water type 5
        !zdeta(1)=1.13;  ! Jerlov water type 7
        !zdeta(2)=0.606; ! Jerlov water type 7
        !zdeta(1)=0.736;  ! Jerlov water type 9
        !zdeta(2)=0.58; ! Jerlov water type 9
        do j=1,9
            absfac_dir=absfac_dir+fsen(j)*exp(-z/(zdeta(j)))
        end do
        rad(i)=I_0*absfac_dir/(rho_0*cp)
        absfac_dir=0.

    case (16)
        ! same as case (9) [Paulson & Simpson (1981)], except change the first absorption coefficient based on data from Soloviev & Schlussel, Boundary-Layer Meteorology (1996).
        zdeta(1)=1/0.066;  ! Jerlov water type I
        !zdeta(1)=1/0.076;  ! Jerlov water type IA
        !zdeta(1)=1/0.088;  ! Jerlov water type IB
        !zdeta(1)=1/0.132;  ! Jerlov water type II
        !zdeta(1)=1/0.382;  ! Jerlov water type III
        !zdeta(1)=1/0.49;   ! Jerlov water type 1
        !zdeta(1)=1/0.7;    ! Jerlov water type 3
        !zdeta(1)=1/1.0;    ! Jerlov water type 5
        !zdeta(1)=1/1.09;   ! Jerlov water type 7
        !zdeta(1)=1/1.6;    ! Jerlov water type 9
        do j=1,9
            absfac_dir=absfac_dir+fsen(j)*exp(-z/(zdeta(j)))
        end do
        rad(i)=I_0*absfac_dir/(rho_0*cp)
        absfac_dir=0.

     case default
        ! extinct_methods 1-7 - these are defined in observations.f90
        ! leave out diffuse / direct discrimination

        !WT From SP, the following is a 2-band parameterization that uses
        ! Jerlov water type (an obsolete index of ocean turbidity)
        ! N.G. Jerlov. Marine Optics. Elsevier, 1976.
        
        ! WT Maybe we should give a warning when this obsolete choice is being used!
        rad(i)=I_0*(A*exp(-z/g1)+(1-A)*exp(-z/g2))/(rho_0*cp)
        !PRINT*,A*exp(-z/g1)+(1-A)*exp(-z/g2)
     end select
     avh(i)=nuh(i)+avmolT 
  end do
  !------------------------------------------------------------------------
  

  !------------------------------------------------------------------------
  ! Compute Q_source
  
  do i=1,nlev
     Q_source(i)=(rad(i)-rad(i-1))/h(i)           !SP 16/05/05
     !include advection source
     Q_source(i)=Q_source(i)+advect(i)                    !SP 08/05/05
     if (t_adv) Q_source(i)=Q_source(i)-u(i)*dtdx(i)-v(i)*dtdy(i) 
  end do
  !------------------------------------------------------------------------

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

double precision function trans(z,extinct_method)

  ! ! Description
  ! Solar transmission rate for extinct_method = 12, 13, 14
  !
  ! For extinct_method 12, see eq(5a,b) and table 5 in
  !   Ohlmann and Siegel, J. Phys. Oceangr. 2000
  ! For extinct_method 13, see eq(5) in
  !   Lee et al, J. Geophysical Research (110), 2005
  ! For exintct_method 14, see eq(5), eq(6a,b,c,d), eq(7) in
  !   Ohlmann, J. Climate (16), 2003
  !
  ! WT 2017-10-19

  integer, intent(in) :: extinct_method
  double precision, intent(in) :: z

  ! Formula coefficients for extinct_method = 12
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

  ! Parameters for extinct_method 12
  double precision		:: para_A,para_K

  ! Parameters for extinct_method 13
  double precision              :: K_ir,K_vis,K1,K2

  ! Parameters for extinct_method 14
  double precision              :: A1,A2,B1,B2,Tr,F_theta,G_CI


  select case (extinct_method)

  case(12)
     !WT The below keep 'coszen' away from zero before taking reciprocal
     ! in O-S.'s formula below (valid up to 75 degress).

     !IF (coszen.lt. 2.49/(7.68*chlo+17.81)) then 
     !   coszen =2.49/(7.68*chlo+17.81)

     ! Begin summing for the tranmission coefficient.
     trans=0.0
     IF(cloud.gt.0.1) then
        !SP revised 20/06/17 setting clear sky limit as 0.1 following Table 1 in Ohlmann&Siegel(2000)
        DO j=1,4
           para_A=C1(j)*chlo+C2(j)*cloud+C4(j)
           para_K=C1(j+4)*chlo+C2(j+4)*cloud+C4(j+4)
           trans=trans+para_A*exp(-para_K*z)
           !WT For z=0, i.e. i=nlev, above is equivalent to:
           !para_A=C1(j)*chlo+C2(j)*cloud+C4(j)               
           !trans = trans+para_A
        END DO
     ELSE
        DO j=9,12
           if (coszen.lt.0.2588) THEN
              !for very low value of coszen the parameterisation breaks down!  SP 22/02/06 revised 28/02/17
              !revised again SP 20/06/17 setting the limit as theta=75degrees (cos(theta)=0.2588) the largest angle in Ohlmann&Siegel(2000)
              
              para_A=C1(j)*chlo+(C3(j)/0.2588)+C4(j)
              para_K=C1(j+4)*chlo+(C3(j+4)/0.2588)+C4(j+4)
           else
              para_A=C1(j)*chlo+(C3(j)/coszen)+C4(j)
              para_K=C1(j+4)*chlo+(C3(j+4)/coszen)+C4(j+4)
           endif
           trans=trans+para_A*exp(-para_K*z)
           !WT For z=0, i.e. i=nlev, above is equivalent to:
           !para_A=C1(j)*chlo+(C3(j)/coszen)+C4(j)
           !trans = trans+para_A
        END DO
     END IF
     
  case (13)
     K1 = (-0.057+0.482*sqrt(abp_coe)+4.221*bb)*(1+0.09*sin(acos(coszen)))
     K2 = (0.183+0.702*abp_coe-2.567*bb)*(1.465-0.667*coszen) 
     K_vis = K1+K2/sqrt(1+z)
     K_ir = (0.560+2.304/(0.001+z)**0.65)*(1+0.002*acos(coszen)*180/(3.1415926))
     !IF (z.gt. 3.0) then  ! E_IR is absorbed within the layer of top 3m.
     !   trans = 0.424*exp(-K_vis*z)
     !ELSE
     trans = 0.576*exp(-K_ir*z)+0.424*exp(-K_vis*z)
     !end if
     
  ! case (14)
  !    !WT 20180815 Warning. Untested case implemented by Joyce.
  !    ! Originally commenting B1, B1 in version 1, but commenting A1, A1 in version 2.
        
  !    ! A1 = 0.0268*log(chlo) + 0.5581                 !version 1: including albedo
  !    ! A2 = -0.017*log(chlo) + 0.2246
  !    ! B1 = -0.0184*chlo**2 + 0.1141*chlo+0.042
  !    ! B2 = 0.5641*chlo**0.13

  !    A1 = 0.571+0.025*log(0.149*chlo)               !version 2: not including albedo
  !    A2 = 0.223+0.010*log(2.329*chlo)
  !    B1 = 0.015+0.176*sqrt(0.462*chlo)
  !    B2 = 0.688+0.060*log(0.125*chlo)
     
  !    Tr = A1*exp(-B1*z) + A2*exp(-B2*z)
  !    F_theta = 0.42*coszen-0.34
  !    IF (cloud.gt. 0.1) then
  !       G_CI = 0
  !    ELSE 
  !       G_CI = 1  
  !    END IF
  !    trans = Tr*(1+F_theta*G_CI)
     
  end select

end function trans

!EOC
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
