!
!
!
!
!
module initialisation
use variables
use special_functions, only: W, phi

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
integer :: Status
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
!   Create Initial values for each faults   
call ini_value()
!
! Pre-calculation of integrand
if (quasi_dynamic == 0) then
	call allocate_fully_dynamic()
	call ini_fully_dynamic()
	call precalculation()
endif
!
! Calculate coefficient of the second fundamental form
call second_fundamental_form_calc()
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
inquire(file='./problems/'//trim(simulation_name),exist=dir_exist) 
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

	

command_line = 'cp '//trim(config_file)//' '//'./problems/'//trim(simulation_name)
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
nb_frequency_x= nb_element_x/2+1
nb_frequency_y= nb_element_y
!
! Allocate wave numbers used in static term 
allocate(Kx(nb_frequency_x,nb_frequency_y))
allocate(Ky(nb_frequency_x,nb_frequency_y))
allocate(Ktot(nb_frequency_x,nb_frequency_y))

! Allocate slip
allocate(slipX(nb_element_x,nb_element_y))   
allocate(slipY(nb_element_x,nb_element_y)) 
allocate(slipX_f(nb_frequency_x,nb_frequency_y)) 
allocate(slipY_f(nb_frequency_x,nb_frequency_y))     
allocate(curvSlip(nb_element_x,nb_element_y)) 
allocate(curvSlip_f(nb_frequency_x,nb_frequency_y)) 

! Allocate V
allocate(Vx_f(nb_frequency_x,nb_frequency_y))
allocate(Vy_f(nb_frequency_x,nb_frequency_y))
allocate(curvV(nb_element_x,nb_element_y))
allocate(curvV_f(nb_frequency_x,nb_frequency_y))

!
! Allocate elastic stresses on the fault
allocate(SigmaNel(nb_element_x,nb_element_y))              ! Derivative of the normal traction on the faults
allocate(SigmaN_dot(nb_element_x,nb_element_y))              ! Derivative of the normal traction on the faults
allocate(SigmaXel(nb_element_x,nb_element_y))
allocate(SigmaYel(nb_element_x,nb_element_y))
allocate(SigmaT0(nb_element_x,nb_element_y))
allocate(SigmaX0(nb_element_x,nb_element_y))
allocate(SigmaY0(nb_element_x,nb_element_y))
allocate(SigmaX_dot(nb_element_x,nb_element_y))              ! Derivative of the tangential traction on the faults
allocate(SigmaY_dot(nb_element_x,nb_element_y))              ! Derivative of the tangential traction on the faults

SigmaNel = 0._dp
SigmaN_dot = 0._dp
SigmaXel = 0._dp
SigmaX_dot = 0._dp
SigmaT0 = 0._dp
SigmaX0 = 0._dp
SigmaY0 = 0._dp
SigmaYel = 0._dp
SigmaY_dot = 0._dp
!                         
!
! Allocate geometry
allocate(curvature(nb_element_x, nb_element_y))
allocate(L1(nb_element_x, nb_element_y))
allocate(M1(nb_element_x, nb_element_y))
allocate(N1(nb_element_x, nb_element_y))
L1 = 0._dp
M1 = 0._dp
N1 = 0._dp

end subroutine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine ini_value()
integer	                                           :: ier
real(kind=dp), dimension(nb_element_x,nb_element_y) :: coeff,V
real(kind=dp) :: k_cut
integer :: k1,k2, Status
!
! Time
time = 0._dp
!
! Calculate ds 
dx = L_fault_x/nb_element_x
dy = L_fault_y/nb_element_y
!
! Calculate time_min
time_min = beta_min*min(dx,dy)/cp
! 
! Calculate V with mask
Vx = Vx*mask+(1-mask)*Vx_mask
Vy = Vy*mask+(1-mask)*Vy_mask


! Calculate slip_f and V_f
call dfftw_execute_dft_r2c(plan_fwd,slipX,slipX_f)
call dfftw_execute_dft_r2c(plan_fwd,slipY,slipY_f)
call dfftw_execute_dft_r2c(plan_fwd,curvSlip,curvSlip_f)
call dfftw_execute_dft_r2c(plan_fwd,Vx,Vx_f)
call dfftw_execute_dft_r2c(plan_fwd,Vy,Vy_f )
call dfftw_execute_dft_r2c(plan_fwd,curvV,curvV_f)
 
!
! Calculate V
V = (Vx**2+Vy**2)**(0.5_dp)
!
! Calculate theta dot
theta_dot = V0/Dc*exp(-theta)- V/Dc
!
! Initialisation sigmaT
coeff = dexp(f0/a+(b/a)*theta)
!
! Initiate sigmaT0, sigmaX0, and sigmaY0 so that the momentum balance equation is respected
sigmaT0 = -sigmaN0*a*asinh(V/(2*V0)*coeff)+ mu/(2._dp*cs)*V
sigmaX0 = sigmaT0*Vx/V*mask
sigmaY0 = sigmaT0*Vy/V*mask
!
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    Allocate array needed for fully-dynamic
!   
subroutine allocate_fully_dynamic()
     
! type mode
allocate(modes(nb_frequency_x,nb_frequency_y))
allocate(modes_temp(nb_frequency_x,nb_frequency_y))

!
! Allocate variable specific to Lapusta 2000
allocate(Tw(nb_frequency_x,nb_frequency_y))
!
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!	Initialize window size for fully dynamic
!
subroutine ini_fully_dynamic()
real(kind=dp) :: k_cut
integer :: k1,k2
!
!
! Calculate Tw (Tw is smaller for higher modes, it correspond to the window over which the time integration is performed)
! TODO:: double check this part
Tw(1,1) = eta_w*max(L_fault_x,L_fault_y)/cs
k_cut = rho_c/(max(L_fault_x,L_fault_y)*eta_w)
!
! For each wavenumber
do k1 = 1,nb_frequency_x
	do k2 = 1,nb_frequency_y
		!
		if (Ktot(k1,k2)<k_cut) then
			Tw(k1,k2) = Tw(1,1)
		else 
			Tw(k1,k2) =  Tw(1,1)*k_cut/Ktot(k1,k2)
		end if
	end do
end do
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Initiate wavenumber
!
subroutine ini_wavenumber()
real(kind=dp), dimension(nb_frequency_x,1) :: ones_line,lines
real(kind=dp), dimension(1,nb_frequency_y) :: ones_column,columns
integer(kind=dp)                           :: k1,k2
!
! Array for derivatives
ones_line = 1._dp
ones_column = 1._dp
do k1=1,nb_element_x/2+1
		lines(k1,1) = (2*pi*(k1-1)/L_fault_x)
end do
do k2=1,nb_element_y/2
		columns(1,k2) = (2*pi*(k2-1)/L_fault_y)
		columns(1,nb_element_y-k2+1) = -(2*pi*(k2)/L_fault_y)
end do
!
! KX KY
Kx = matmul(lines,ones_column)
Ky = matmul(ones_line,columns)
Ktot = (Kx**2+Ky**2)**(0.5_dp)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Precalculate the integrand
!
subroutine precalculation()
integer :: k1,k2,i
integer :: nb_time


! For each mode
do k1=1,nb_frequency_x
	do k2=1,nb_frequency_y
		
		! Find the number of point to compute
		nb_time = ceiling(Tw(k1,k2) / time_min)
		!
		! Allocate
		modes(k1,k2)%nb_time = nb_time
		allocate(modes(k1,k2)%integrand_shear_modeII(nb_time))
		allocate(modes(k1,k2)%integrand_shear_modeIII(nb_time))
		allocate(modes(k1,k2)%integrand_normal(nb_time))
		allocate(modes(k1,k2)%Vx_f(nb_time))
		allocate(modes(k1,k2)%Vy_f(nb_time))
		allocate(modes(k1,k2)%curvV_f(nb_time))
	
		modes_temp(k1,k2)%nb_time = nb_time
		allocate(modes_temp(k1,k2)%integrand_shear_modeII(nb_time))
		allocate(modes_temp(k1,k2)%integrand_shear_modeIII(nb_time))
		allocate(modes_temp(k1,k2)%integrand_normal(nb_time))
		allocate(modes_temp(k1,k2)%Vx_f(nb_time))
		allocate(modes_temp(k1,k2)%Vy_f(nb_time))
		allocate(modes_temp(k1,k2)%curvV_f(nb_time))
		!
		! Current time id 
		modes(k1,k2)%current_time_id = 1
		!
		! Precalculate integrand for shear (midpoint integration)
		do i =1,nb_time
			modes(k1,k2)%integrand_shear_modeII(i) = time_min*calculate_integrand_shear_modeII(Ktot(k1,k2)*(time_min/2+(i-1)*time_min))
			modes(k1,k2)%integrand_shear_modeIII(i) = time_min*calculate_integrand_shear_modeIII(Ktot(k1,k2)*(time_min/2+(i-1)*time_min))
			! modes(k1,k2)%Vx_f(i) = 0._dp
! 			modes(k1,k2)%Vy_f(i) = 0._dp
			modes(k1,k2)%Vx_f(i) = Vx_f(k1,k2)
 			modes(k1,k2)%Vy_f(i) = Vx_f(k1,k2)


			modes(k1,k2)%curvV_f(i) = 0._dp
		end do
		!
		! Precalculate integrand for normal (midpoint integration)
		do i =1,nb_time
			modes(k1,k2)%integrand_normal(i) = time_min*calculate_integrand_normal(Ktot(k1,k2)*(time_min/2+(i-1)*time_min))
		end do
		!print*,'2*pi*k/L_fault*(time_min/2+(i-1)*time_min)',2*pi*k/L_fault*(time_min/2+(i-1)*time_min)
		!print*,'calculate_integrand_shear(2*pi*k/L_fault*(time_min/2+(i-1)*time_min))',calculate_integrand_shear(2*pi*k/L_fault*(time_min/2+(i-1)*time_min))
	
	
		!
		! Initiate first value of V_f and curvV_f
		! modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id) = Vx_f(k1,k2)
! 		modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id) = Vy_f(k1,k2)
		modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id) = curvV_f(k1,k2)
		!
		! Initiate modes_temp
		modes_temp(k1,k2)%Vx_f = modes(k1,k2)%Vx_f
		modes_temp(k1,k2)%Vy_f = modes(k1,k2)%Vy_f
		modes_temp(k1,k2)%curvV_f = modes(k1,k2)%curvV_f
		modes_temp(k1,k2)%integrand_shear_modeII = modes(k1,k2)%integrand_shear_modeII
		modes_temp(k1,k2)%integrand_shear_modeIII = modes(k1,k2)%integrand_shear_modeIII
		modes_temp(k1,k2)%integrand_normal = modes(k1,k2)%integrand_normal
		modes_temp(k1,k2)%current_time_id = 1
		
	end do 
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
y = 0.0


end function calculate_integrand_normal
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Calculate second fundamental form coefficients
!
subroutine second_fundamental_form_calc()
integer :: k1,k2, Status
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y) :: L1_f, M1_f, N1_f
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)  :: height_f                   ! height of the fault at each x and y coordinates
complex(kind=dp), dimension(nb_frequency_x,1) :: ones_line,lines
complex(kind=dp), dimension(1,nb_frequency_y) :: ones_column,columns

! Initialise arrays
! Array for second fundamental form
L1_f = 0._dp
M1_f = 0._dp
N1_f = 0._dp
!
! Calculate Fourier transform of the surface
call dfftw_execute_dft_r2c(plan_fwd,height,height_f)
!
! Calculate derivatives (in Fourier domain)
L1_f = -KX*KX*height_f
M1_f = -KX*KY*height_f
N1_f = -KY*KY*height_f
! 
! Compute backward Fourier transform
call dfftw_execute_dft_c2r(plan_bwd,L1_f/(nb_element_x*nb_element_y),L1)
call dfftw_execute_dft_c2r(plan_bwd,M1_f/(nb_element_x*nb_element_y),M1)
call dfftw_execute_dft_c2r(plan_bwd,N1_f/(nb_element_x*nb_element_y),N1)

end subroutine
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Configure FFT
!
subroutine configure_fft()
!
! Create plant for 2d fourier transform using FFTW
CALL dfftw_plan_dft_r2c_2d(plan_fwd,nb_element_x,nb_element_y,Vx,Vx_f,FFTW_ESTIMATE)
CALL dfftw_plan_dft_c2r_2d(plan_bwd,nb_element_y,nb_element_x,Vx_f,Vx,FFTW_ESTIMATE)


end subroutine
 
 
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       
!     Make some checking at the beginning of the calculation
!
subroutine issue_warning()
real(kind=dp), dimension(nb_element_x,nb_element_y)  :: Lb,Lnuc
!
! Calculate minimum process zone size (Rice 1993)
Lb = 2*mu*Dc/(abs(sigmaN0)*b)*(1._dp-cs**2/cp**2)
!
! calculate nucleation lenghscale minimum
Lnuc = 2*Lb/(pi*(1._dp-a/b)**2)
!
! Check if discretisation is ok compared to process zone size
if (min(dx,dy)>minval(Lb)/4) then
	print*,'dx,dy is usually required to be lower than Lb/5'
	print*,'dx,dy:', dx,dy
	print*,'Lb:', minval(Lb)
	print*,'Are you sure you want to continue?'
	!read(*,*)
end if
!
! Check if L_fault is ok compared to length of the fault
! if (max(L_fault_1,<maxval(Lnuc)) then
! 	print*,'L_fault is small than Lnuc, your fault may not rupture'
! 	print*,'L_fault:', L_fault
! 	print*,'Lnuc:', maxval(Lnuc)
! 	print*,'Are you sure you want to continue?'
! 	!read(*,*)
! end if
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