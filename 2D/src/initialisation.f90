!
!
!
!
!
module initialisation
use variables
use special_functions, only: W, phi
use friction
implicit none

include 'fftw3.f'

!
!
!
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Main subroutine 
!
subroutine initialise()
!
!   Check if output directory exist OK
call mkdir_simulation_name()
!
!   copy input files in the output directory OK
call cp_inputFile()
!
!   Allocate all the arrays OK
call allocate_ini()
!
! Preconfigure fft calls
call configure_fft()
!
!   Initiate absolute value of wave number use in the static term 
call ini_wavenumber()
!
! Calculate coefficient of the second fundamental form
call curvature_calc()
!
!   Create Initial values for each faults   
call ini_value()
!
! Pre-calculation of integrand
! Pre-calculation of integrand
if (quasi_dynamic == 0) then
	call allocate_fully_dynamic()
	call ini_fully_dynamic()
	call precalculation()
endif
!
! Issue warning if hyperparameter or discretisation is wrong
call issue_warning()
!


end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Make a directory which name is the simulation name
!
subroutine mkdir_simulation_name()
character(len=1000)              :: command_line
logical                         :: dir_exist
!
command_line = 'mkdir ./problems/'//trim(simulation_name)
inquire(FILE='./problems/'//trim(simulation_name),exist=dir_exist) 
if(.NOT. dir_exist) then
    call system(command_line)
end if
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    Copy config.in, data_mesh.in, config_fault.in in the output_dir
!
subroutine cp_inputFile()
character(len=1000)   :: command_line
!
!call execute_command_line ('rm -f ./problems/'//trim(simulation_name)//'/*')
!call execute_command_line ('clear')
print*,' //////////////////////////////////////////////////'
print*,' //                                              //'
print*,' //                 FASTCYCLES2D                 //'
print*,' //  Pierre Romanet, So Ozawa, and Harsha Bhat   //'
print*,' //                                              //'
print*,' //////////////////////////////////////////////////'
print*,''
call OMP_SET_NUM_THREADS(omp_threads);
print*,'	Running on ', omp_threads, 'threads'
	

command_line = 'cp -r '//trim(path_problem)//'/* '//'./problems/'//trim(simulation_name)
call system(trim(command_line))

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   allocate arrays and pointer needed by fault-type and solver 
!   
subroutine allocate_ini()
!
!   ALLOCATE COMMON ARRAYS
!
! Number of frequency
nb_frequency = nb_element/2+1
!
! Allocate wave number 
allocate(Kwave(nb_frequency))
! 
! Allocate slip
allocate(slip(nb_element))   
allocate(slip_f(nb_frequency))  
allocate(curvSlip(nb_element))   
allocate(curvSlip_f(nb_frequency))  

! Allocate V
allocate(V_f(nb_frequency))
allocate(curvV(nb_element))
allocate(curvV_f(nb_frequency))
!
! Allocate elastic stresses on the fault
allocate(SigmaNel(nb_element))              ! Derivative of the normal traction on the faults
allocate(SigmaN_dot(nb_element))              ! Derivative of the normal traction on the faults
allocate(SigmaTel(nb_element))
allocate(SigmaT0(nb_element))
allocate(SigmaT_dot(nb_element))              ! Derivative of the tangential traction on the faults
SigmaNel = 0._dp
SigmaN_dot = 0._dp
SigmaT0 = 0._dp
SigmaTel = 0._dp
SigmaT_dot = 0._dp
!       
! Allocate geometry
allocate(curvature(nb_element))
         

! Pointer to friction law
print*,'friction_law',friction_law
if(trim(friction_law) == 'RateState_R') then
    fric_law => rate_state_R
elseif(friction_law == 'RateState_fh') then
	fric_law => rate_state_fh
elseif(friction_law == 'RateState_fh') then
	fric_law => rate_state
else
    print*,'Problem in config.in: wrong friction law, please choose between "RateState_R", "RateState_fh" or "RateState"'
    print*,'Maybe the friction with the associated mode of fracture you wanted is not available yet'
    stop
end if
!
! Pointer to state evolution law
if(state_evolution == 'aging') then
    state_evol => aging_law
elseif(state_evolution == 'slip') then
    state_evol => slip_law  
else
    print*,'Problem in config.in: wrong friction law, please choose between "aging" or "slip"'
    print*,'Maybe the state evolution law you wanted is not available yet'
    stop
end if
!
end subroutine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine ini_value()
integer	                                           :: ier
real(kind=dp), dimension(nb_element) 			   :: coeff
integer 										   :: i, status
!
! Time
time = 0._dp
!
! Calculate ds 
ds = L_fault/nb_element
!
! Calculate time_min
time_min = beta_min*ds/cp
! 
! Calculate V with mask
V = V*mask+(1-mask)*V_mask
! Calculate slip_f and V_f
!Status = DftiComputeForward( Desc_fwd, slip(:),slip_f(:) )
call dfftw_execute_dft_r2c(plan_fwd,slip(:),slip_f(:) )
curvSlip = curvature*slip
!Status = DftiComputeForward( Desc_fwd, curvSlip(:),curvSlip_f(:) )
!Status = DftiComputeForward( Desc_fwd, V(:),V_f(:) )
call dfftw_execute_dft_r2c(plan_fwd, curvSlip(:),curvSlip_f(:) )
call dfftw_execute_dft_r2c(plan_fwd,V(:),V_f(:)  )
curvV = curvature*V
!Status = DftiComputeForward( Desc_fwd, curvV(:),curvV_f(:) )
call dfftw_execute_dft_r2c(plan_fwd,curvV(:),curvV_f(:))

!
! Calculate theta dot
theta_dot = V0/Dc*exp(-theta) - (V)/Dc
!
! Initialisation sigmaT
coeff = exp(f0/a+(b/a)*theta)
sigmaT0 = -sigmaN0*a*asinh(V/(2*V0)*coeff) + mu/(2._dp*cs)*V
!
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Allocate array needed for fully-dynamic
!   
subroutine allocate_fully_dynamic()
     
! type mode
allocate(modes(nb_frequency))
allocate(modes_temp(nb_frequency))
!
! Allocate variable specific to Lapusta 2000
allocate(Tw(nb_frequency))
!
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!	Initialize window size for fully dynamic
!
subroutine ini_fully_dynamic()
integer(kind=dp) :: mode_cut
integer :: i
!
!
! Calculate Tw (Tw is smaller for higher modes)
Tw(1) = eta_w*L_fault/cs
mode_cut = ceiling(rho_c/(2*pi*eta_w))
if (mode_cut>nb_frequency) then
	Tw(1:nb_frequency) = Tw(1)
else
	Tw(2:mode_cut) = Tw(1)
	! Calculate cutting time 
	do i=mode_cut+1,nb_frequency
		Tw(i) = Tw(1)*mode_cut/i
	end do
end if
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Initiate wavenumber
!
subroutine ini_wavenumber()
integer(kind=dp)                           :: k1
!
! Array for derivatives
do k1=1,nb_frequency
		Kwave(k1) = (2*pi*(k1-1)/L_fault)
end do
!
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Precalculate the integrand
!
subroutine precalculation()
integer :: k,i
integer :: nb_time



! For each mode
do k=1,nb_frequency
	
	!
	! Find the number of point to compute
	nb_time = ceiling(Tw(k) / time_min)
	!
	! Allocate
	modes(k)%nb_time = nb_time
	allocate(modes(k)%integrand_shear(nb_time))
	allocate(modes(k)%integrand_normal(nb_time))
	allocate(modes(k)%V_f(nb_time))
	allocate(modes(k)%curvV_f(nb_time))
	
	modes_temp(k)%nb_time = nb_time
	allocate(modes_temp(k)%integrand_shear(nb_time))
	allocate(modes_temp(k)%integrand_normal(nb_time))
	allocate(modes_temp(k)%V_f(nb_time))
	allocate(modes_temp(k)%curvV_f(nb_time))
	!
	! Current time id 
	modes(k)%current_time_id = 1
	!
	! Precalculate integrand for shear (midpoint integration)
	do i =1,nb_time
		modes(k)%integrand_shear(i) = time_min*calculate_integrand_shear_modeII(2*pi*(k-1)/L_fault*(time_min/2+(i-1)*time_min))
		modes(k)%V_f(i) = 0._dp
		modes(k)%curvV_f(i) = 0._dp
	end do
	!
	! Precalculate integrand for normal (midpoint integration)
	do i =1,nb_time
		modes(k)%integrand_normal(i) = time_min*calculate_integrand_normal(2*pi*(k-1)/L_fault*(time_min/2+(i-1)*time_min))
	end do
	!
	! Initiate first value of V_f and curvV_f
	modes(k)%V_f(modes(k)%current_time_id) = V_f(k)
	modes(k)%curvV_f(modes(k)%current_time_id) = curvV_f(k)
	!
	! Initiate modes_temp
	modes_temp(k)%V_f = modes(k)%V_f
	modes_temp(k)%curvV_f = modes(k)%curvV_f
	modes_temp(k)%integrand_shear = modes(k)%integrand_shear
	modes_temp(k)%integrand_normal = modes(k)%integrand_normal
	modes_temp(k)%current_time_id = 1
	
end do

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Function that calculate the integrand for shear traction
!
function calculate_integrand_shear_modeII(x) result(y)
real(kind=dp) :: x
real(kind=dp) :: y

! TODO: double check this part
y = (cs**2*x**2+cs**2/cp**2)*W(cp*x)-(1+cs**2*x**2)*W(cs*x)      &
    -cs**2/cp**2*bessel_J1(cp*x)+(0.5)*bessel_J1(cs*x)        &
    -(cs**2/cp*x)*bessel_J0(cp*x)+(cs*x)*bessel_J0(cs*x)

 end function calculate_integrand_shear_modeII
 
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Function that calculate the integrand for shear traction
!
function calculate_integrand_shear_modeIII(x) result(y)
real(kind=dp) :: x
real(kind=dp) :: y

! TODO: double check this part
y = 1/2._dp*W(cs*x)    

 end function calculate_integrand_shear_modeIII
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Function that calculate the integrand for normal traction
!
function calculate_integrand_normal(x) result(y)
real(kind=dp) :: x
real(kind=dp) :: y

! TODO
y = 0._dp


end function calculate_integrand_normal
 
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Calculate curvature
!
subroutine curvature_calc()
integer :: k1, Status
complex(kind=dp), dimension(nb_frequency) :: height_dot_f,height_2dot_f
real(kind=dp), dimension(nb_element) :: height_dot,height_2dot
complex(kind=dp), dimension(nb_frequency)  :: height_f                   ! height of the fault at each x and y coordinates
!
! Calculate Fourier transform of the surface
!Status = DftiComputeForward( Desc_fwd, height(:),height_f(:) )
call dfftw_execute_dft_r2c(plan_fwd,height(:),height_f(:))

!
! Calculate derivatives (in Fourier domain)
height_dot_f = cmplx(0._dp,1._dp)*Kwave*height_f
height_2dot_f = -Kwave**2._dp*height_f
! 
! Compute backward Fourier transform
!Status = DftiComputeBackward( Desc_bwd, height_dot_f(:),height_dot(:) )
!Status = DftiComputeBackward( Desc_bwd, height_2dot_f(:),height_2dot(:) )
call dfftw_execute_dft_c2r(plan_bwd,height_dot_f(:)/nb_element,height_dot(:))
call dfftw_execute_dft_c2r(plan_bwd, height_2dot_f(:)/nb_element,height_2dot(:))


!
! Calculation of curvature
curvature = height_2dot/(1._dp+height_dot**2._dp)**(1.5_dp)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Configure FFT
!
subroutine configure_fft()
integer :: Status



CALL dfftw_plan_dft_r2c_1d(plan_fwd,nb_element,V,V_f,FFTW_ESTIMATE)
CALL dfftw_plan_dft_c2r_1d(plan_bwd,nb_element,V_f,V,FFTW_ESTIMATE)




! Describe the dimension of the Fourier transform to perform
! Status = DftiCreateDescriptor(Desc_fwd, DFTI_DOUBLE,DFTI_REAL, 1, (/nb_element/) )
! Status = DftiCreateDescriptor(Desc_bwd, DFTI_DOUBLE,DFTI_REAL, 1, (/nb_element/) )
! 
!                                                                   
! Prepare stide
! Status = DftiSetValue(Desc_fwd, DFTI_INPUT_STRIDES,(/0,1,nb_element/))
! Status = DftiSetValue(Desc_fwd,DFTI_OUTPUT_STRIDES, (/0,1,nb_frequency/))
! 
! Status = DftiSetValue(Desc_bwd, DFTI_INPUT_STRIDES, (/0,1,nb_frequency/))
! Status = DftiSetValue(Desc_bwd,DFTI_OUTPUT_STRIDES,(/0,1,nb_element/))
! 
! 
! The output array is consider as conjugate even (because the input array is real)
! Status = DftiSetValue(Desc_fwd,DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
! Status = DftiSetValue(Desc_bwd,DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
! 
! Compute the FFT in a different array as input
! Status = DftiSetValue( Desc_fwd, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
! Status = DftiSetValue( Desc_bwd, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
! 
! Set factor
! Status = DftiSetValue( Desc_bwd,DFTI_BACKWARD_SCALE,1._dp/(nb_element))
! 
! Commit descriptor
! Status = DftiCommitDescriptor( Desc_bwd )
! Status = DftiCommitDescriptor( Desc_fwd )


end subroutine
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Make some checking at the beginning of the calculation
!
subroutine issue_warning()
real(kind=dp), dimension(nb_element)  :: Lb,Lnuc
!
! Calculate minimum process zone size (Rice 1993)
Lb = 2*mu*Dc/(abs(sigmaN0)*b)*(1._dp-cs**2/cp**2)
!
! calculate nucleation lenghscale minimum
Lnuc = 2*Lb/(pi*(1._dp-a/b)**2)
!
! Check if discretisation is ok compared to process zone size
if (ds>minval(Lb)/4) then
	print*,'ds is usually required to be lower than Lb/5'
	print*,'ds:', ds
	print*,'Lb:', minval(Lb)
	print*,'Are you sure you want to continue?'
	!read(*,*)
end if
!
! Check if L_fault is ok compared to length of the fault
if (L_fault<maxval(Lnuc)) then
	print*,'L_fault is small than Lnuc, your fault may not rupture'
	print*,'L_fault:', L_fault
	print*,'Lnuc:', maxval(Lnuc)
	print*,'Are you sure you want to continue?'
	!read(*,*)
end if
!
! Check hyperparameter
! beta_min
if (beta_min.gt.0.5) then
	print*,'beta_min is usually require to be lower than 1/2.'
	print*,'beta_min :', beta_min
	print*,'beta_min controls the minimum time discretisation.'
	print*,'Are you sure you want to continue?'
! 	read(*,*)
end if 


end subroutine 
end module 