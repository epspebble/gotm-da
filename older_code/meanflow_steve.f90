!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Mean Flow
!
! !INTERFACE:
   module meanflow
!
! !DESCRIPTION:
!  This module provides all variables necessary for the meanflow
!  calculation and also makes the proper initialisations.
!
! !USES:
   IMPLICIT NONE
!  Default all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_meanflow
!
! !PUBLIC DATA MEMBERS:

!  coordinate z, layer thicknesses
   double precision, public, dimension(:), allocatable  :: z,h,ho

!  the velocity components
   double precision, public, dimension(:), allocatable  :: u,v,w

!  potential temperature, salinity
   double precision, public, dimension(:), allocatable  :: T,S

!  buoyancy frequency and shear-frequency
   double precision, public, dimension(:), allocatable  :: NN,SS

!  shear and buoyancy production of tke
   double precision, public, dimension(:), allocatable  :: P,B

!  buoyancy, short-wave radiation, 
!  extra production of tke by see-grass etc
   double precision, public, dimension(:), allocatable  :: buoy,rad,xP

!  a temporal diffusivity
   double precision, public, dimension(:), allocatable  :: avh

!  grid-related vertical velocity
   double precision, public, dimension(:), allocatable  :: w_grid

!  extra friction terms due to e.g. seagrass
   double precision, public, dimension(:), allocatable  :: fric,drag

!  Internal heat source
   double precision, public, dimension(:), allocatable	:: Qsour

!  the 'meanflow' namelist
   double precision, public                    :: h0b=0.05
   double precision, public                    :: z0s_min=0.02
   logical,  public                    :: charnok=.false.
   double precision, public                    :: charnok_val=1400.
   double precision, public                    :: ddu=0.
   double precision, public                    :: ddl=0.
   integer,  public                    :: grid_method=1
   double precision, public                    :: c1ad=0.8
   double precision, public                    :: c2ad=0.0
   double precision, public                    :: c3ad=0.1
   double precision, public                    :: c4ad=0.1
   double precision, public                    :: Tgrid=3600.
   double precision, public                    :: NNnorm=0.2
   double precision, public                    :: SSnorm=0.2
   double precision, public                    :: dsurf=10.0
   double precision, public                    :: dtgrid=5.
   character(LEN=255), public     :: grid_file='grid.dat'
   double precision, public                    :: gravity=9.81
   double precision, public                    :: rho_0=1027.
   double precision, public                    :: cp=3985.
   double precision, public                    :: avmolu=1.3e-6
   double precision, public                    :: avmolT=1.4e-7
   double precision, public                    :: avmolS=1.1e-9
   integer,  public                    :: MaxItz0b=10
   logical,  public                    :: no_shear=.false.

   double precision, public	  :: u_taus_cr !HK added for SH skin correction


!  the surface roughness length
   double precision, public                    :: z0b,z0s

!  the coriolis parameter
   double precision, public                    :: cori

!  the friction velocities
   double precision, public                    :: u_taub,u_taus

!  other stuff
   integer,  public                    :: eq_state_method
   double precision, public                    :: depth0=0.
   double precision, public                    :: depth
   double precision, public                    :: obs_heat_content=0.
   double precision, public                    :: calc_heat_content=0.
!
! !DEFINED PARAMETERS:
   double precision, public, parameter         :: pi=3.141592654
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: meanflow.F90,v $
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:15:01  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:50:06  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/11/18 15:58:02  gotm
!  Vertical grid can now be read from file
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!
!  private date members
   double precision, parameter       :: omega=2*pi/86400.
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialisation of the mean flow variables
!
! !INTERFACE:
   subroutine init_meanflow(namlst,fn,nlev,latitude)
!
! !DESCRIPTION:
!  Allocates memory and initialises everything related 
!  to the `meanflow' component of GOTM.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                  :: namlst
   character(len=*), intent(in)         :: fn
   integer, intent(in)                  :: nlev
   double precision, intent(in)                 :: latitude
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the meanflow module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc

   namelist /meanflow/  h0b,z0s_min,charnok,charnok_val,ddu,ddl,       &
                        grid_method,c1ad,c2ad,c3ad,c4ad,Tgrid,NNnorm,  &
                        SSnorm,dsurf,dtgrid,grid_file,gravity,rho_0,cp,&
                        avmolu,avmolT,avmolS,MaxItz0b,no_shear,u_taus_cr
!
!-----------------------------------------------------------------------
!BOC
   write(0,*) '   ', 'init_meanflow'

   open(namlst,file=fn,status='old',action='read',err=80)
   write(0,*) '       ', 'reading meanflow namelists..'
   read(namlst,nml=meanflow,err=81)
   close (namlst)
   write(0,*) '       ', 'done.'

   write(0,*) '       ', 'allocation meanflow memory..'
   allocate(z(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (z)'
   z = 0.

   allocate(h(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (h)'
   h = 0.

   allocate(ho(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (ho)'
   ho = 0.

   allocate(u(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (u)'
   u = 0.

   allocate(v(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (v)'
   v = 0.

   allocate(w(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (w)'
   w = 0.

   allocate(fric(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (fric)'
   fric = 0.

   allocate(drag(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (drag)'
   drag = 0.

   allocate(T(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (T)'
   T = 0.

   allocate(S(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (S)'
   S = 0.

   allocate(NN(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (NN)'
   NN = 0.

   allocate(SS(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (SS)'
   SS = 0.

   allocate(P(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_meanflow: Error allocating (P)'
   P = 0.

   allocate(B(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_meanflow: Error allocating (B)'
   B = 0.

   allocate(xP(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (xP)'
   xP = 0.

   allocate(buoy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (buoy)'
   buoy = 0.

   allocate(rad(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (rad)'
   rad = 0.

   allocate(avh(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (avh)'
   avh = 0.

   allocate(w_grid(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_meanflow: Error allocating (w_grid)'
   w_grid = 0.


   write(0,*) '       ', 'done.'

   depth0=depth
   z0b=0.03*h0b

   z0s=z0s_min    ! lu (otherwise z0s is not initialised

   cori=2*omega * sin(2*pi*latitude/360.)

   return
80 write(0,*) 'FATAL ERROR: ', 'I could not open: ',trim(fn)
   stop 'init_meanflow'
81 write(0,*) 'FATAL ERROR: ', 'I could not read "meanflow" namelist'
   stop 'init_meanflow'

   end subroutine init_meanflow
!EOC

!-----------------------------------------------------------------------

   end module meanflow

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
