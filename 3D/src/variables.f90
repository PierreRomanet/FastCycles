module variables
implicit none
!
!
!   PARAMETERS
integer, parameter                                 :: dp=kind(1.d0)
real(kind=dp), parameter                           :: pi=3.141592653589793238_dp
!
!   OUTPUT 
character(len=50)                                  :: simulation_name              ! Read in config.in
character(len=100)                                 :: config_file = ""             !"config.in"
character(len=200)                                 :: path_problem = ""             !"config.in"
integer                                            :: freq_writing_file            ! Read in config.in
!
!   SIMULATION VARIABLES
integer				                               :: omp_threads
real(kind=dp)                                      :: time
real(kind=dp)                                      :: time_step_max 
integer                                            :: max_it
integer                                            :: stride_time                  ! Save time each stride_time
character(len=20)                                  :: friction_law                     ! Read in config.in
character(len=20)                                  :: state_evolution               ! Read in config.in
integer											   :: quasi_dynamic                      ! Read in config.in                
!
!   MESHING VARIABLES
integer           	                               :: nb_element_x, nb_element_y                 ! Calculated in interpolater
integer											   :: nb_frequency_x,nb_frequency_y 
real(kind=dp)									   :: L_fault_x
real(kind=dp)									   :: L_fault_y
real(kind=dp), dimension(:,:), allocatable        :: height                   ! height of the fault at each x and y coordinates
real(kind=dp), dimension(:,:), allocatable           :: L1                      ! coefficients of the second fundamental form.
real(kind=dp), dimension(:,:), allocatable           :: M1                      ! coefficients of the second fundamental form.
real(kind=dp), dimension(:,:), allocatable           :: N1                      ! coefficients of the second fundamental form.
real(kind=dp), dimension(:,:), allocatable           :: curvature             ! Allocated and calculated in initialisation
integer      , dimension(:,:), allocatable           :: mask                    ! Allocated and calculated in initialisation
real(kind=dp)                                      :: dx                           
real(kind=dp)                                      :: dy   
!      
!  FOURIER TRANSFORM
integer*8											:: plan_fwd,plan_bwd
real(kind=dp), dimension(:,:), allocatable           :: Kx,Ky,Ktot               ! Array for wavenumbers (initialise in initialisation)

!
!   PROBLEM VARIABLES
real(kind=dp), dimension(:,:), allocatable           :: Vx                    ! Allocated in initialisation
real(kind=dp), dimension(:,:), allocatable           :: Vy                    ! Allocated in initialisation
complex(kind=dp), dimension(:,:), allocatable        :: Vx_f                  ! Allocated in initialisation
complex(kind=dp), dimension(:,:), allocatable        :: Vy_f                  ! Allocated in initialisation
real(kind=dp), dimension(:,:), allocatable        :: curvV                 ! Allocated in initialisation
complex(kind=dp), dimension(:,:), allocatable        :: curvV_f               ! Allocated in initialisation

real(kind=dp), dimension(:,:), allocatable           :: slipX                 ! Allocated in initialisation
real(kind=dp), dimension(:,:), allocatable           :: slipY                 ! Allocated in initialisation

complex(kind=dp), dimension(:,:), allocatable        :: slipX_f               ! Allocated in initialisation 
complex(kind=dp), dimension(:,:), allocatable        :: slipY_f               ! Allocated in initialisation 

real(kind=dp), dimension(:,:), allocatable           :: curvSlip               ! Allocated in initialisation 
complex(kind=dp), dimension(:,:), allocatable        :: curvSlip_f               ! Allocated in initialisation 

real(kind=dp)                                      :: mu                   ! Shear modulus
real(kind=dp)                                      :: cs
real(kind=dp)                                      :: cp

real(kind=dp), dimension(:,:), allocatable           :: shearX_loading_dot        ! Background Loading (calculated in initialisation)
real(kind=dp), dimension(:,:), allocatable           :: shearY_loading_dot        ! Background Loading (calculated in initialisation)
real(kind=dp), dimension(:,:), allocatable           :: normal_loading_dot       ! Background Loading (calculated in initialisation)
real(kind=dp)									     :: Vx_mask
real(kind=dp)									     :: Vy_mask
real(kind=dp), dimension(:,:), allocatable           :: SigmaN0                  ! Total Normal traction on the faults   
real(kind=dp), dimension(:,:), allocatable           :: SigmaNel               ! Total Normal traction on the faults   
real(kind=dp), dimension(:,:), allocatable           :: SigmaN_dot           ! Derivative of the normal traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaXel               ! Tangential traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaX_dot           ! Derivative of the tangential traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaYel               ! Tangential traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaY_dot           ! Derivative of the tangential traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaT0              ! Tangential traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaX0              ! Tangential traction on the faults
real(kind=dp), dimension(:,:), allocatable           :: SigmaY0              ! Tangential traction on the faults

!    
! Friction law
! procedure(friction_laws), pointer                  :: fric_law => null()   ! Allocated in initialisation
! procedure(state_laws), pointer                     :: state_evol   => null()   ! Allocated in initialisation
! abstract interface
!     subroutine friction_laws(f,f_dot,V_i,theta_i,nb_el)
!     	integer, parameter                           :: dp=kind(1.d0)
!     	integer           	                         :: nb_el            
!     	real(kind=dp), dimension(nb_el)               :: f,f_dot
!     	real(kind=dp), dimension(nb_el)               :: V_i, theta_i
!     end subroutine
!     function state_laws(V_i,theta_i,nb_el)
!    		integer, parameter                           :: dp=kind(1.d0)
!     	integer           	                         :: nb_el                  
!     	real(kind=dp), dimension(nb_el)              :: state_laws
!     	real(kind=dp), dimension(nb_el)              :: V_i, theta_i
!     end function
! end interface
                     
!   RATE AND STATE VARIABLES
real(kind=dp), dimension(:,:), allocatable           :: theta
real(kind=dp), dimension(:,:), allocatable           :: theta_dot
real(kind=dp), dimension(:,:), allocatable           :: V0, f0
real(kind=dp), dimension(:,:), allocatable           :: a, b
real(kind=dp), dimension(:,:), allocatable           :: Dc
!
! TYPE
type mode
    integer                                       :: nb_time
    integer                                       :: current_time_id                                                 ! Calculated in interpolater
    complex(kind=dp), dimension(:), allocatable   :: Vx_f
    complex(kind=dp), dimension(:), allocatable   :: Vy_f
    complex(kind=dp), dimension(:), allocatable   :: CurvV_f
    real(kind=dp), dimension(:), allocatable      :: integrand_shear_modeII
    real(kind=dp), dimension(:), allocatable      :: integrand_shear_modeIII
	real(kind=dp), dimension(:), allocatable      :: integrand_normal
end type
!   INTERNAL VARIABLES: type(fault) faults
type(mode), dimension(:,:), allocatable              :: modes                 ! Allocated in input
type(mode), dimension(:,:), allocatable              :: modes_temp            ! Allocated in input
!
! SPECIFIC TO LAPUSTA 2000
real(kind=dp)                                      :: time_min
real(kind=dp), dimension(:,:), allocatable         :: Tw          ! Allocated in initialisation, this is the mode dependent maximum time of integration, see Lapusta 2009
real(kind=dp)                                      :: beta_min
real(kind=dp)                                      :: gamma
real(kind=dp)                                      :: rho_c      ! Allocated in initialisation, cutting frequency entering the definition of Tw
real(kind=dp)                                      :: eta_w      ! Allocated in initialisation, cutting frequency entering the definition of Tw
real(kind=dp)                                      :: nr_accuracy   ! Convergence of Newton raphson
real(kind=dp)                                      :: tol_solver   ! Convergence of Newton raphson

!
!


end module variables
