!
!       SOLVER FOR THE TIME STEP, ADAPTED FROM LAPUSTA ET AL. 2000
!
module solver
use variables
use hdf5

!
!   Variables
implicit none
private



! Saving
real(kind=dp), dimension(:,:), allocatable         :: V_save, theta_save, slip_save,tractionNel_save,tractionTel_save  ! Tmp arrays to save data
real(kind=dp), dimension(:), allocatable           :: time_save                      ! Tmp array for solver (scaling)
integer                                            :: it_count_save        ! Number of calls to bsstep


! Time variable
real(kind=dp)                                      :: dtime_try,dtime_next
integer(kind=dp)                                   :: n_try
integer											   :: time_id


! Variable for the solver
real(kind=dp), dimension(:), allocatable           :: slip_half1, slip_half2, slip_full   
real(kind=dp), dimension(:), allocatable           :: theta_half1,theta_half2, theta_full
real(kind=dp), dimension(:), allocatable           :: V_half1, V_half2, V_full
real(kind=dp), dimension(:), allocatable           :: shear_loading_half1, shear_loading_half2, shear_loading_full

real(kind=dp), dimension(:), allocatable           :: theta_dot_half1, theta_dot_half2,theta_dot_full

complex(kind=dp), dimension(:), allocatable        :: slip_half1_f, slip_half2_f, slip_full_f
complex(kind=dp), dimension(:), allocatable        :: V_half1_f, V_half2_f, V_full_f
complex(kind=dp), dimension(:), allocatable        :: curvV_half1_f, curvV_half2_f, curvV_full_f

!
!   Public 
public :: solve
! 
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Solver routine 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine solve()
!
! Allocate arrays
call allocate_arrays()
!
! Initialisation
it_count_save = 0
dtime_try = time_min
dtime_next = 0._dp
!
! Enter loop 
do time_id=1,max_it
	print*,'------------------------------------------------------------------------------'
	print*,'------------------------------------------------------------------------------'
	print*,'------------------------------------------------------------------------------'

	print*,'time id: ',time_id
	print*,'time: ',time,' s'
 	print*,'V',minval(V),maxval(V)
 	print*,'theta',minval(theta),maxval(theta)
 	print*,'theta_dot',minval(theta_dot),maxval(theta_dot)

 	print*,'slip',minval(slip),maxval(slip)
 		print*,'dtime_try: ',dtime_try
	print*,'time_min',2*time_min

	! save data
	call save_data()

	!
	! Adaptive time stepping
	call stepsolver(dtime_try,dtime_next)   
	print*,'dtime_try,dtime_next',dtime_try,dtime_next

	time = time + dtime_try
	dtime_try = dtime_next
	!
	! Update the new value of the fields slip, theta, theta_dot, V
	call update_values()
	!
	! Save variable
    if(it_count_save==freq_writing_file) then
        call write_data()
        it_count_save = 0
    end if     
    !
end do
! Deallocate arrays.
call deallocate_arrays()

end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   One step of the solver
!
subroutine stepsolver(dt_try,dt_next)   


! Time variable
real(kind=dp)                                      :: dt_try,dt_next

! Internal

real(kind=dp)                                      :: error
integer											   :: istop, status,k
  
  	
  	! Initilisation
  	istop = 0
  	status = 0
  
  	! Until a good time step is not found
  	Do while(istop/=1)
  	
  	
  	! If the simulation is dynamic, dt_try is a multiple of time_min
  	if (quasi_dynamic == 0) then
  		! Transform dt_try in a multiple of t_min
  		! Calculate the number of time time_ev_temp contains time_min)
		n_try = int(dt_try/(2._dp*time_min),dp)
		!
		! Calculate time evolution final (ensure that time_ev is a multiple of time min)
		dt_try = max(time_min,n_try*time_min)
		if (dt_try-time_min<1d-9) then
			n_try = 1
		end if
	! If the simulation is quasi-dynamic, dt_try is whatever it needs to be
	elseif (quasi_dynamic ==1) then
		dt_try = dt_try/2._dp
	else
		print*,'Please chose quasi dynamin or fully-dynamic simulation'
		stop
	endif

	! If the simulation is fully dynamic
  	if(quasi_dynamic ==0) then
		! Initiate history
		do k=1,nb_frequency
  			modes_temp(k)%current_time_id = modes(k)%current_time_id
  			modes_temp(k)%V_f = modes(k)%V_f
  			modes_temp(k)%curvV_f = modes(k)%curvV_f
  		end do
  	endif
  	
  	
  	! Initiate shear_loading
  	shear_loading_half1 = shear_loading_dot*(time+dt_try)
    shear_loading_half2 = shear_loading_dot*(time+2._dp*dt_try)
    shear_loading_full = shear_loading_dot*(time+2._dp*dt_try)
  	
  	!! The Order of prediction matters here, because the history (modes_temps) is updated when using two steps. 
  	! Make prediction with one single time step
  	
  	
  	
  	!$OMP PARALLEL SECTIONS
	!$OMP SECTION
  ! 	print*,''
!   	print*,'----------------------------------'
!   	print*,'Making full step calculation'
!   	print*,'----------------------------------'
	
	
! 	print*,'slip',minval(slip),maxval(slip)
! 	print*,'slip',minval(abs(slip_f)),maxval(abs(slip_f))
! 	print*,'theta',minval(theta),maxval(theta)
! 	print*,'V',minval(V),maxval(V)
! 	print*,'V_f',minval(abs(V_f)),maxval(abs(V_f))
! 	print*,'theta_dot',minval(theta_dot),maxval(theta_dot)
	


	
	call make_prediction(slip_full,slip,   &
						slip_full_f,slip_f, &
  						theta_full,theta, &
  						V_full, V, &
  						V_full_f, V_f, &
  						curvV_full_f, curvV_f, &
  						theta_dot_full,theta_dot, &
  						shear_loading_full, &
  						2._dp*dt_try,1,status)



  	! Make first prediction for half the time step
  	! Make prediction of state variable OK 
  	
  	!$OMP SECTION
  	! print*,''
!   	print*,'----------------------------------'
!   	print*,'Making first half step calculation'
! 	print*,'----------------------------------'
	

	
  	call make_prediction(slip_half1,slip,   &
						slip_half1_f,slip_f, &
  						theta_half1,theta, &
  						V_half1, V, &
  						V_half1_f, V_f, &
  						curvV_half1_f, curvV_f, &
  						theta_dot_half1,theta_dot, &
  						shear_loading_half1, &
  						dt_try,0,status)
  						

	
		
	
	! Update history contained in modes_temps for the second step
	if (quasi_dynamic == 0 ) then
		call update_values_half()
	endif
! 	print*,''
! 	print*,'----------------------------------'
! 	print*,'Making second half step calculation'
! 	print*,'----------------------------------'



	
	
	! Make second prediction for half the time step
	call make_prediction(slip_half2,slip_half1,   &
						slip_half2_f,slip_half1_f, &
  						theta_half2,theta_half1, &
  						V_half2, V_half1, &
  						V_half2_f, V_half1_f, &
  						curvV_half2_f, curvV_half1_f, &
  						theta_dot_half2,theta_dot_half1, &
  						shear_loading_half2, &
  						dt_try,0,status,SigmaNel,SigmaTel)
  						
  						
  						
 
	
	!$OMP END PARALLEL SECTIONS
	
	
	! Compare error between full and half time step and find the next time
	call compare_error_find_next_time_step(dt_try,dt_next,istop,status)
	

end do  
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   make prediction 
!
subroutine compare_error_find_next_time_step(dt_try,dt_next,istop,status)
! Time variable
real(kind=dp)                                      :: dt_try,dt_next

! Internal
real(kind=dp)                                      :: error
integer											   :: istop, status,k


	! Compare the two results, calculation of the error
	error = 0
	do k=1,nb_element
		if (mask(k) ==1) then
			error = max(error,abs((slip_full(k)-slip_half2(k))/slip_half2(k)))
			error = max(error,abs(exp(theta_full(k)-theta_half2(k))-1._dp))
		end if
	end do

	print*,'error:', error
	
	! What do we do next?
	if ((error<tol_solver) .and. (status == 0)) then ! Everything went fine
  		print*,'Time step was correctly resolved.'
  		dt_next = min(2.0*0.9* dt_try * (tol_solver/error)**(1/3._dp),90.0*86400.0)
  		print*,'dt_next',dt_next
  		dt_try = 2._dp*dt_try 
  		istop = 1
  		
  	else if ( (quasi_dynamic == 0) .and. (dt_try == time_min) .and. (status == 0)) then ! If the time is already at the minimum, do not account for error in dynamic case
  	  	print*,'Time step did not take into account the error (already at minimum time step).'
  		dt_next = min(2.0*0.9* dt_try * (tol_solver/error)**(1/3._dp),90.0*86400.0)
  		print*,'dt_next',dt_next
  		dt_try = 2._dp*dt_try
  		istop = 1
  		
  		
  	else if ((error>tol_solver) .and. (status == 0)) then ! Decrease the time step to meet the requirement on error
  	  	print*,'Time step was decreased to meet the error requirement'
  		dt_try = 2.0*0.9* dt_try * (tol_solver/error)**(1/3._dp)
  		print*,'dt_next',dt_try
  		
  
  		
    else if ((dt_try > time_min) .and. (status == 1)) then    ! In case one newton-raphson failed decrease the time step
      	print*,'Time step was decreased (Newton-Raphson failed)'
    	dt_try = dt_try/2
    	status = 0
    			
  	else if ((quasi_dynamic == 0) .and. (dt_try == time_min) .and. (status == 1)) then ! If the time is already at the minimum, and Newton Raphson fails, and we are in dynamic case stop the calculation
  	  	print*,'Calculation stoped (Newton-Raphson failed and already at minimum time step'
  	  	stop
  	end if 

end subroutine




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   make prediction 
!
subroutine make_prediction(slip_e, slip_i,   &
						 slip_f_e, slip_f_i, &
						 theta_e,theta_i, &
						 V_e, V_i, &
						 V_f_e,V_f_i, &
						 curvV_f_e,curvV_f_i, &
						 theta_dot_e, theta_dot_i, &
						 shear_loading_i, &
						 dt_try,i_full, status,sigmaN_e, sigmaT_e ) 
						 
	! IN 			 
	real(kind=dp), dimension(nb_element)        :: slip_e,slip_i
	complex(kind=dp), dimension(nb_frequency)     :: slip_f_e, slip_f_i
	complex(kind=dp), dimension(nb_frequency)     :: curvSlip_f_e
	real(kind=dp), dimension(nb_element)        :: theta_e, theta_i
	real(kind=dp), dimension(nb_element)        :: V_e,V_i
	complex(kind=dp), dimension(nb_frequency)     :: V_f_e,V_f_i
	complex(kind=dp), dimension(nb_frequency)     :: curvV_f_e,curvV_f_i
	real(kind=dp), dimension(nb_element)        :: theta_dot_e, theta_dot_i

	complex(kind=dp), dimension(nb_frequency)     :: integral_history_shear
	complex(kind=dp), dimension(nb_frequency)     :: integral_history_normal
	
	real(kind=dp), dimension(nb_element), optional :: sigmaN_e, sigmaT_e
	real(kind=dp), dimension(nb_element)        :: sigmaN_temp, sigmaT_temp
	real(kind=dp), dimension(nb_element)        :: shear_loading_i

	real(kind=dp)								:: dt_try		 

	integer										:: status, i_full
	
	! Internal
	complex(kind=dp), dimension(nb_frequency)     :: integral_shear		
		print*,''




	! Initialisation
	integral_history_shear = 0._dp
	integral_history_normal = 0._dp

	! Make prediction slip and theta
	call make_prediction_theta_slip(slip_e,theta_e,slip_f_e,curvSlip_f_e,slip_i,theta_i,slip_f_i,V_i,V_f_i,theta_dot_i,dt_try)
	!
	! Make prediction of normal and shear traction
	call make_prediction_traction(sigmaN_temp, sigmaT_temp,slip_f_e,curvSlip_f_e,V_f_i,curvV_f_i,integral_history_shear,integral_history_normal,i_full,1,status)
	!
	! Find predictive values of V and theta	 
	call make_prediction_V_thetaDot(sigmaN_temp, sigmaT_temp,theta_e,V_e,V_i,theta_dot_e,V_f_e,curvV_f_e,shear_loading_i, status)
	!
	! Calculate final prediction of slip and state variable
	V_e = (V_e+V_i)/2._dp
	V_f_e = (V_f_e+V_f_i)/2._dp ! CHECK IF V_f_e IS OK
	theta_dot_e = (theta_dot_e+theta_dot_i)/2._dp

	! Make prediction slip and theta
	call make_prediction_theta_slip(slip_e,theta_e,slip_f_e,curvSlip_f_e,slip_i,theta_i,slip_f_i,V_e,V_f_e,theta_dot_e,dt_try)
	!
	! Make final prediction of normal and shear traction	
	call make_prediction_traction(sigmaN_temp, sigmaT_temp,slip_f_e,curvSlip_f_e,V_f_e,curvV_f_e,integral_history_shear,integral_history_normal,i_full,0,status)
    !
	! Calculate V2 and theta2_dot	 
	call make_prediction_V_thetaDot(sigmaN_temp, sigmaT_temp,theta_e,V_e,V_i,theta_dot_e,V_f_e,curvV_f_e, shear_loading_i, status)
	!
	! Eventually copy shear and normal traction
	if (present(sigmaN_e).and.present(sigmaT_e)) then
		sigmaN_e = sigmaN_temp
		sigmaT_e = sigmaT_temp
	end if

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   make_prediction_theta
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine make_prediction_theta_slip(slip_e,theta_e,slip_f_e,curvSlip_f_e,slip_i,theta_i,slip_f_i,V_i,V_f_i,theta_dot_i,dt_try)



! IN
real(kind=dp), dimension(nb_element)        :: slip_i,V_i
real(kind=dp), dimension(nb_element)        :: theta_i,theta_dot_i

complex(kind=dp), dimension(nb_frequency)     :: slip_f_i,V_f_i
!
real(kind=dp)								:: dt_try

! OUT
real(kind=dp), dimension(nb_element)        :: slip_e, theta_e, curvSlip_e
complex(kind=dp), dimension(nb_frequency)     :: slip_f_e, curvSlip_f_e
!
!integer										:: status_fft

!
! Prediction theta
theta_e = theta_i + dt_try*theta_dot_i*mask


! Prediction slip
slip_e = slip_i + dt_try*V_i
slip_f_e = slip_f_i + dt_try*V_f_i
!
! Prediction curvSlip 
curvSlip_e = curvature*slip_e
! status_fft = DftiComputeForward( Desc_fwd, curvSlip_e(:),curvSlip_f_e(:) ) 
call dfftw_execute_dft_r2c(plan_fwd,curvSlip_e(:),curvSlip_f_e(:) )
!
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   make_prediction_traction
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine make_prediction_traction(sigmaN_temp, sigmaT_temp,slip_f_e,curvSlip_f_e,V_f_e,curvV_f_e,integral_history_shear,integral_history_normal,i_full,i_history,status)

! Internal
complex(kind=dp), dimension(nb_frequency)          :: integral_shear
complex(kind=dp), dimension(nb_frequency)          :: integral_normal

! Inout
real(kind=dp), dimension(nb_element)        :: sigmaN_temp, sigmaT_temp
complex(kind=dp), dimension(nb_frequency)     :: integral_history_shear
complex(kind=dp), dimension(nb_frequency)     :: integral_history_normal
integer(kind=dp)        :: n_ev
integer :: k, id_time,id_slip, i_history, i_full
complex(kind=dp), dimension(nb_frequency)     :: slip_f_e,curvSlip_f_e,V_f_e,curvV_f_e
integer										:: status

!

! Calculate integral 
integral_shear = 0._dp
integral_normal = 0._dp





! If the calculation is fully dynamic
if (quasi_dynamic == 0) then

	! Select full or half step
	if (i_full==1) then
		n_ev = 2 * n_try
	else if (i_full==0) then
		n_ev = n_try
	end if



	!For each mode 
	do k=1,nb_frequency

		!Calculate the integral on past history
		if(i_history==1)then
		
	
			if (n_ev+1.le.modes(k)%nb_time) then
				do id_time=n_ev+1,modes(k)%nb_time				
			
					! Calculate id_slip
					id_slip = modulo(modes(k)%current_time_id-id_time+n_ev,modes(k)%nb_time)+1
				
					! Shear
					integral_history_shear(k) = integral_history_shear(k) + &
					modes(k)%integrand_shear(id_time) * modes(k)%V_f(id_slip)
		
					! Normal 
					integral_history_normal(k) = integral_history_normal(k) + &
					modes(k)%integrand_normal(id_time) * modes(k)%curvV_f(id_slip)
		
				end do	
			end if
		end if
    
		! Calculate the integral on current history
		do id_time=1,min(modes(k)%nb_time,n_ev)
		 
			! Shear
			integral_shear(k) = integral_shear(k) + &
					modes(k)%integrand_shear(id_time) * V_f_e(k) 
					      
			! Normal 
			integral_normal(k) = integral_normal(k) + &
							modes(k)%integrand_normal(id_time) * curvV_f_e(k) 

		end do
	
	end do
endif




! Add static term and history term
! Shear 
integral_shear =-mu*Kwave*((1-cs**2/cp**2)*slip_f_e + integral_shear + integral_history_shear)
 
! Normal
integral_normal = mu*cmplx(0._dp,1._dp)*((1-cs**2/cp**2)*curvSlip_f_e+integral_normal + integral_history_normal) 
	          	        
!

! Make sure that the 0 frequency is 0
integral_normal(1) = 0._dp
integral_shear(1) = 0._dp

!
! Inverse fourier transform of the shear traction
! status_fft = DftiComputeBackward( Desc_bwd, integral_shear(:),SigmaT_temp(:) )
! status_fft = DftiComputeBackward( Desc_bwd, integral_normal(:),SigmaN_temp(:) )
call dfftw_execute_dft_c2r(plan_bwd,integral_shear(:)/nb_element,SigmaT_temp(:) )
call dfftw_execute_dft_c2r(plan_bwd, integral_normal(:)/nb_element,SigmaN_temp(:) )

!
! Check that normal traction  is not becoming positive
if (maxval(SigmaN0+SigmaNel).ge.0) then
	print*,'Error: traction became negative'
	if(n_try == 1) then
		stop
	else
		status = 1
	end if
end if


end subroutine




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   find_predicted_slip
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine make_prediction_V_thetaDot(sigmaN_temp, sigmaT_temp,theta_e,V_e,V_i,theta_dot_e,V_f_e,curvV_f_e,shear_loading,status)
real(kind=dp), dimension(nb_element)        ::sigmaN_temp, sigmaT_temp
real(kind=dp), dimension(nb_element)        :: theta_e,V_e,V_i
real(kind=dp), dimension(nb_element)        :: theta_dot_e
real(kind=dp), dimension(nb_element)        :: shear_loading

complex(kind=dp), dimension(nb_frequency)     :: V_f_e,curvV_f_e
! Newton-raphson
real(kind=dp), dimension(nb_element)        :: res, res_dot
real(kind=dp), dimension(nb_element)        :: curvV_e

real(kind=dp), dimension(nb_element)        :: coeff1


integer   :: i, status
!
! Initialise
res = 10.0*nr_accuracy/V_e
res_dot =  1._dp

V_e = V_i
! Calculate coeff1
coeff1 = exp(f0/a+(b/a)*theta)* mask
!
! Initialise Newton raphson
i=0
!
! Use Newton raphson 

do while (maxval(abs((res/res_dot)*mask/V_e))>nr_accuracy)

	! Calculate friction coefficient and V-derivative of friction coefficient 
	!call fric_law(f,f_dot,V,theta,nb_element)
! 
	! Calculate new function to minimize and its derivative
	
	res = (sigmaT0+ sigmaT_temp + shear_loading - mu/(2._dp*cs)*V_e & 
	   +(sigmaN0+sigmaN_temp)*a*asinh(V_e/(2._dp*V0)*coeff1))*mask
	res_dot = - mu/(2*cs) + (SigmaN0+SigmaN_temp)*a/sqrt(1._dp/coeff1**2+(V_e/(2*V0))**2)/(2._dp*V0)
	
	! Update V1
	V_e = V_e -  (res/res_dot)*mask

	   	   
	! Check convergence
	i=i+1
! 	print*,'it nr: ',i
	if (i > 20) then
		print*,'number of iteration (NF):',i
	    print*,'maxval(abs((f/f_dot)*mask))',maxval(abs((res/res_dot)*mask/V_e))
	    print*,'V_e',minval(V_e),maxval(V_e)
	end if 
	
	
	if (i > 40) then
		print*,'Newton-Raphson (1) did not converge'
		status = 1
		exit
	end if   
end do
!
! 


! Calculate theta_dot
theta_dot_e = V0/Dc*exp(-theta) - (V_e)/Dc
!
! Calculate V_f_e and curvV_f_e
curvV_e = curvature*V_e
! status_fft = DftiComputeForward( Desc_fwd, V_e(:),V_f_e(:) )
! status_fft = DftiComputeForward( Desc_fwd, curvV_e(:),curvV_f_e(:) )
call dfftw_execute_dft_r2c(plan_fwd,  V_e(:),V_f_e(:) )
call dfftw_execute_dft_r2c(plan_fwd, curvV_e(:),curvV_f_e(:) )

end subroutine




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   save data 
!
subroutine save_data()   
   if (modulo(time_id-1,stride_time) == 0) then 

        it_count_save = it_count_save + 1
        V_save(it_count_save,:) = V
        slip_save(it_count_save,:) = slip
        theta_save(it_count_save,:) = theta
        time_save(it_count_save) = time 
        tractionNel_save(it_count_save,:) = sigmaNel
        tractionTel_save(it_count_save,:) = sigmaTel-mu/(2._dp*cs)*V
    end if
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Write data in a file
!
subroutine write_data()       


character(len=30)   :: num2str
integer :: error
integer(HSIZE_T) :: data_dims_array(2), data_dims_time(1)
integer(HID_T)   :: file_id, dspace_id_V,dspace_id_time, dset_id_V,dset_id_time
integer(HID_T)   :: dspace_id_theta,dset_id_theta
integer(HID_T)   :: dspace_id_tractionTel,dset_id_tractionTel
integer(HID_T)   :: dspace_id_tractionNel, dset_id_tractionNel



if (it_count_save.gt.0) then

	! Initialize Fortran interface
	CALL h5open_f(error)   
	!
	! Create a new file
	write(num2str,'(I7.7)') time_id
	CALL h5fcreate_f('problems/'//trim(simulation_name)//'/output'//trim(num2str), H5F_ACC_TRUNC_F, file_id, error)
	!
	! Specify ranks and dims
	
	data_dims_array(1) = freq_writing_file
	data_dims_array(2) = nb_element
	data_dims_time = freq_writing_file
	!
	! Create dataspace (the dataset is next) "dspace_id" is returned
	CALL h5screate_simple_f(2, data_dims_array, dspace_id_V, error)
	CALL h5screate_simple_f(2, data_dims_array, dspace_id_theta, error)
	CALL h5screate_simple_f(2, data_dims_array, dspace_id_tractionTel, error)
	CALL h5screate_simple_f(2, data_dims_array, dspace_id_tractionNel, error)
	CALL h5screate_simple_f(1, data_dims_time, dspace_id_time, error)

	
	! Create dataset with default properties "dset_id" is returned
	CALL h5dcreate_f(file_id, 'V', H5T_NATIVE_DOUBLE, dspace_id_V, dset_id_V, error)
	CALL h5dcreate_f(file_id, 'theta', H5T_NATIVE_DOUBLE, dspace_id_theta, dset_id_theta, error)
	CALL h5dcreate_f(file_id, 'tractionTel', H5T_NATIVE_DOUBLE, dspace_id_tractionTel, dset_id_tractionTel, error)
	CALL h5dcreate_f(file_id, 'tractionNel', H5T_NATIVE_DOUBLE, dspace_id_tractionNel, dset_id_tractionNel, error)
	CALL h5dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, dspace_id_time, dset_id_time, error)

                     
	! Write dataset 
	CALL h5dwrite_f(dset_id_V, H5T_NATIVE_DOUBLE, V_save(1:it_count_save,1:nb_element), data_dims_array, error)
	CALL h5dwrite_f(dset_id_theta, H5T_NATIVE_DOUBLE, theta_save(1:it_count_save,1:nb_element), data_dims_array, error)
	CALL h5dwrite_f(dset_id_tractionTel, H5T_NATIVE_DOUBLE, tractionTel_save(1:it_count_save,1:nb_element), data_dims_array, error)
	CALL h5dwrite_f(dset_id_tractionNel, H5T_NATIVE_DOUBLE, tractionNel_save(1:it_count_save,1:nb_element), data_dims_array, error)
	CALL h5dwrite_f(dset_id_time, H5T_NATIVE_DOUBLE, time_save(1:it_count_save), data_dims_time, error)
                      
	! Close access to dataset 1
	CALL h5dclose_f(dset_id_V, error)
	CALL h5dclose_f(dset_id_theta, error)
	CALL h5dclose_f(dset_id_tractionTel, error)
	CALL h5dclose_f(dset_id_tractionNel, error)
	CALL h5dclose_f(dset_id_time, error)

	! Close access to data space 1
	CALL h5sclose_f(dspace_id_V, error)
	CALL h5sclose_f(dspace_id_theta, error)
	CALL h5sclose_f(dspace_id_tractionTel, error)
	CALL h5sclose_f(dspace_id_tractionNel, error)
	CALL h5sclose_f(dspace_id_time, error)
	
	! Close access to file 1
	CALL h5fclose_f(file_id, error)
	
	! Close fortran interface
	call h5close_f(error) 
    
end if


end subroutine


  



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   update_values
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_values_half()
integer :: end_id,k



! Overwrite history
do k=1,nb_frequency

	! No problem to overwrite
	if (modes_temp(k)%current_time_id+n_try .le. modes_temp(k)%nb_time) then 
	
	
		! Overwrite
		modes_temp(k)%V_f(modes_temp(k)%current_time_id+1:modes_temp(k)%current_time_id+n_try) =  (V_half1_f(k)+	modes_temp(k)%V_f(modes_temp(k)%current_time_id))/2._dp
		modes_temp(k)%curvV_f(modes_temp(k)%current_time_id+1:modes_temp(k)%current_time_id+n_try) =  (curvV_half1_f(k)+modes_temp(k)%curvV_f(modes_temp(k)%current_time_id))/2._dp
		
		
		
		
		! Set new current time
		modes_temp(k)%current_time_id = modes_temp(k)%current_time_id + n_try	


		!
	! Start exactly at the beginning of the file
	else if (modes_temp(k)%current_time_id+1 .gt. modes_temp(k)%nb_time) then
		!
		! Calculate new end at the begining of the array
		end_id = modulo(modes_temp(k)%current_time_id+n_try,modes_temp(k)%nb_time)
		! Overwrite data at the beginning
		modes_temp(k)%V_f(1:end_id) =  (V_half1_f(k)+modes_temp(k)%V_f(modes_temp(k)%current_time_id))/2._dp
		modes_temp(k)%curvV_f(1:end_id) =  (curvV_half1_f(k)+modes_temp(k)%curvV_f(modes_temp(k)%current_time_id))/2._dp
		
	
		! Set new current time
		modes_temp(k)%current_time_id = end_id
	! All other cases
	else 
		! Calculate new end at the beginning of the array
		end_id = modulo(modes_temp(k)%current_time_id+n_try,modes_temp(k)%nb_time)
		
		
		! Overwrite data at the end
		modes_temp(k)%V_f(1:end_id) =  (V_half1_f(k)+modes_temp(k)%V_f(modes_temp(k)%current_time_id))/2._dp
		modes_temp(k)%curvV_f(1:end_id) =  (curvV_half1_f(k)+modes_temp(k)%curvV_f(modes_temp(k)%current_time_id))/2._dp


		! Overwrite data at the beginning
		modes_temp(k)%V_f(modes_temp(k)%current_time_id+1:modes_temp(k)%nb_time) =  (V_half1_f(k)+modes_temp(k)%V_f(modes_temp(k)%current_time_id))/2._dp
		modes_temp(k)%curvV_f(modes_temp(k)%current_time_id+1:modes_temp(k)%nb_time) =  (curvV_half1_f(k)+modes_temp(k)%curvV_f(modes_temp(k)%current_time_id))/2._dp

		! Set new current time
		modes_temp(k)%current_time_id = end_id
		

	end if 
end do


end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   update_values
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_values()
integer :: end_id,k


slip  = slip_half2
slip_f = slip_half2_f
theta = theta_half2
theta_dot = theta_dot_half2

! Update V
V = V_half2
V_f = V_half2_f
curvV_f = curvV_half2_f

if (quasi_dynamic ==0) then
	! Update to half1
	do k=1,nb_frequency
  			modes(k)%current_time_id = modes_temp(k)%current_time_id
  			modes(k)%V_f = modes_temp(k)%V_f
  			modes(k)%curvV_f = modes_temp(k)%curvV_f
	end do


	! Update to half2
	do k=1,nb_frequency
		! No problem to overwrite
		if (modes(k)%current_time_id+n_try .le. modes(k)%nb_time) then 
			! Overwrite
			modes(k)%V_f(modes(k)%current_time_id+1:modes(k)%current_time_id+n_try) = (V_f(k)+modes(k)%V_f(modes(k)%current_time_id))/2._dp
			modes(k)%curvV_f(modes(k)%current_time_id+1:modes(k)%current_time_id+n_try) = (curvV_f(k)+modes(k)%curvV_f(modes(k)%current_time_id))/2._dp


			! Set new current time
			modes(k)%current_time_id = modes(k)%current_time_id + n_try	
			!
		! Start exactly at the beginning of the file
		else if (modes(k)%current_time_id+1 .gt. modes(k)%nb_time) then
			! Calculate new end at the begining of the array
			end_id = modulo(modes(k)%current_time_id+n_try,modes(k)%nb_time)
			! Overwrite data at the beginning
			modes(k)%V_f(1:end_id) = (V_f(k)+modes(k)%V_f(modes(k)%current_time_id))/2._dp
			modes(k)%curvV_f(1:end_id) = (curvV_f(k)+modes(k)%curvV_f(modes(k)%current_time_id))/2._dp
	
			! Set new current time
			modes(k)%current_time_id = end_id
		! All other cases
		else 
			! Calculate new end at the beginning of the array
			end_id = modulo(modes(k)%current_time_id+n_try,modes(k)%nb_time)
			! Overwrite data at the beginning
			modes(k)%V_f(1:end_id) = (V_f(k)+modes(k)%V_f(modes(k)%current_time_id))/2._dp
			modes(k)%curvV_f(1:end_id) = (curvV_f(k)+modes(k)%curvV_f(modes(k)%current_time_id))/2._dp

			! Overwrite data at the beginning
			modes(k)%V_f(modes(k)%current_time_id+1:modes(k)%nb_time) = (V_f(k)+modes(k)%V_f(modes(k)%current_time_id))/2._dp
			modes(k)%curvV_f(modes(k)%current_time_id+1:modes(k)%nb_time) = (curvV_f(k)+modes(k)%curvV_f(modes(k)%current_time_id))/2._dp


			! Set new current time
			modes(k)%current_time_id = end_id
		end if 
	end do
endif


end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Allocate arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine allocate_arrays()
!
! Allocation for array try
allocate(slip_half1(nb_element))
allocate(slip_half2(nb_element))
allocate(slip_full(nb_element))

allocate(slip_half1_f(nb_element))
allocate(slip_half2_f(nb_element))
allocate(slip_full_f(nb_element))



allocate(theta_half1(nb_element))
allocate(theta_half2(nb_element))
allocate(theta_full(nb_element))

allocate(V_half1(nb_element))
allocate(V_half2(nb_element))
allocate(V_full(nb_element))

allocate(V_half1_f(nb_element))
allocate(V_half2_f(nb_element))
allocate(V_full_f(nb_element))

allocate(curvV_half1_f(nb_element))
allocate(curvV_half2_f(nb_element))
allocate(curvV_full_f(nb_element))

allocate(theta_dot_half1(nb_element))
allocate(theta_dot_half2(nb_element))
allocate(theta_dot_full(nb_element))

allocate(shear_loading_full(nb_element))
allocate(shear_loading_half1(nb_element))
allocate(shear_loading_half2(nb_element))
!
! Save data
allocate(time_save(freq_writing_file))
allocate(V_save(freq_writing_file,nb_element))
allocate(slip_save(freq_writing_file,nb_element))
allocate(tractionNel_save(freq_writing_file,nb_element))
allocate(tractionTel_save(freq_writing_file,nb_element))


allocate(theta_save(freq_writing_file,nb_element))

!
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Deallocate arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine deallocate_arrays()
integer :: status_fft

deallocate(slip_half1)
deallocate(slip_half2)
deallocate(slip_full)
deallocate(theta_half1)
deallocate(theta_half2)
deallocate(theta_full)
deallocate(V_half1)
deallocate(V_half2)
deallocate(V_full)
deallocate(theta_dot_half1)
deallocate(theta_dot_half2)
deallocate(theta_dot_full)


! Free Descriptor
! status_fft = DftiFreeDescriptor( Desc_bwd )
! status_fft = DftiFreeDescriptor( Desc_fwd )
call dfftw_destroy_plan(plan_bwd)
call dfftw_destroy_plan(plan_fwd)
end subroutine


end module

!%%
