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
real(kind=dp), dimension(:,:,:), allocatable         :: Vx_save,Vy_save, theta_save, slipX_save,slipY_save,tractionNel_save,tractionXel_save,tractionYel_save  ! Tmp arrays to save data
real(kind=dp), dimension(:), allocatable           :: time_save                      ! Tmp array for solver (scaling)
integer                                            :: it_count_save        ! Number of calls to bsstep


! Time variable
real(kind=dp)                                      :: dtime_try,dtime_next
integer(kind=dp)                                   :: n_try
integer											   :: time_id


! Variable for the solver
real(kind=dp), dimension(:,:), allocatable           :: slipX_half1, slipX_half2, slipX_full   
real(kind=dp), dimension(:,:), allocatable           :: slipY_half1, slipY_half2, slipY_full  
real(kind=dp), dimension(:,:), allocatable           :: theta_half1,theta_half2, theta_full
real(kind=dp), dimension(:,:), allocatable           :: Vx_half1, Vx_half2, Vx_full
real(kind=dp), dimension(:,:), allocatable           :: Vy_half1, Vy_half2, Vy_full


real(kind=dp), dimension(:,:), allocatable           :: theta_dot_half1, theta_dot_half2,theta_dot_full

complex(kind=dp), dimension(:,:), allocatable        :: slipX_half1_f, slipX_half2_f, slipX_full_f
complex(kind=dp), dimension(:,:), allocatable        :: slipY_half1_f, slipY_half2_f, slipY_full_f
complex(kind=dp), dimension(:,:), allocatable        :: Vx_half1_f, Vx_half2_f, Vx_full_f
complex(kind=dp), dimension(:,:), allocatable        :: Vy_half1_f, Vy_half2_f, Vy_full_f

complex(kind=dp), dimension(:,:), allocatable        :: curvV_half1_f, curvV_half2_f, curvV_full_f
!
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
! SPECIFIC TO BP7SEAS
! real(kind=dp), dimension(:,:), allocatable       :: perturb_space
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

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
use omp_lib
real(kind=dp)           :: t_beg,t_end

!
! Allocate arrays
call allocate_arrays()
!



!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
! call calculate_perturb_space()
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

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
 	print*,'Vx',minval(Vx),maxval(Vx)
 	print*,'Vy',minval(Vy),maxval(Vy)

 	print*,'theta',minval(theta),maxval(theta)
 	print*,'theta_dot',minval(theta_dot),maxval(theta_dot)

 	print*,'slipX',minval(slipX),maxval(slipX)
 	print*,'slipY',minval(slipY),maxval(slipY)
 	print*,'SigmaXel',minval(SigmaXel),maxval(SigmaXel)
 	print*,'SigmaYel',minval(SigmaYel),maxval(SigmaYel)
 	
 	print*,'dtime_try: ',dtime_try
	print*,'time_min',2*time_min

	! save data
	call save_data()
	!
	! Adaptive time stepping
	t_beg = omp_get_wtime()
	call stepsolver(dtime_try,dtime_next)   
	t_end = omp_get_wtime()
	print*,'Total time one time step = ',t_end-t_beg
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

!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   PERTURBATION SPECIFIC TO BP7
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! subroutine calculate_perturb_space()
! integer(kind=dp)										   :: k1,k2
! real(kind=dp)										   :: Lx,Ly,dR
! 
! !$OMP PARALLEL DO
! do k2=1,nb_element_y
! 	do k1=1,nb_element_x
! 		Lx = k1*dx-nb_element_x/2*dx-dx/2
!         Ly = k2*dy-nb_element_y/2*dy-dy/2
!         dR = sqrt((Lx+50._dp)**2._dp+(Ly+50._dp)**2._dp)
!         
!         if (dR<=150._dp) then
! 			perturb_space(k1,k2) = exp(dR**2/(dR**2-150._dp**2))
! 		else
! 			perturb_space(k1,k2) = 0._dp
! 		endif
! 	enddo
! enddo
! !$OMP END PARALLEL DO
! end subroutine
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   One step of the solver
!
subroutine stepsolver(dt_try,dt_next)   


! Time variable
real(kind=dp)                                      :: dt_try,dt_next

! Internal
real(kind=dp)                                      :: error, slip_full, slip_half2
real(kind=dp)                                      :: error1, error2

integer											   :: istop, status,k1,k2
  
  	
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
		do k1=1,nb_frequency_x
			do k2= 1,nb_frequency_y
  				modes_temp(k1,k2)%current_time_id = modes(k1,k2)%current_time_id
  				modes_temp(k1,k2)%Vx_f = modes(k1,k2)%Vx_f
  				modes_temp(k1,k2)%Vy_f = modes(k1,k2)%Vy_f
  				modes_temp(k1,k2)%curvV_f = modes(k1,k2)%curvV_f
  			end do
  		end do
  	endif
  	
  	
  	!! The Order of prediction matters here, because the history (modes_temps) is updated when using two steps. 
  	! Make prediction with one single time step
  	
  	
  	
  	!$OMP PARALLEL SECTIONS
	!$OMP SECTION
! 	print*,''
!   	print*,'----------------------------------'
!   	print*,'Making full step calculation'
!   	print*,'----------------------------------'
	
	
	
	call make_prediction(slipX_full,slipX,   &
						slipX_full_f,slipX_f, &
						slipY_full,slipY,   &
						slipY_full_f,slipY_f, &
  						theta_full,theta, &
  						Vx_full, Vx, &
  						Vx_full_f, Vx_f, &
  						Vy_full, Vy, &
  						Vy_full_f, Vy_f, &
  						curvV_full_f, curvV_f, &
  						theta_dot_full,theta_dot, &
  						2._dp*dt_try,1,status)



  	! Make first prediction for half the time step
  	! Make prediction of state variable OK 
  	
  	!$OMP SECTION
  	! print*,''
!   	print*,'----------------------------------'
!   	print*,'Making first half step calculation'
! 	print*,'----------------------------------'
	

	
  	call make_prediction(slipX_half1,slipX,   &
						slipX_half1_f,slipX_f, &
						slipY_half1,slipY,   &
						slipY_half1_f,slipY_f, &
  						theta_half1,theta, &
  						Vx_half1, Vx, &
  						Vx_half1_f, Vx_f, &
  						Vy_half1, Vy, &
  						Vy_half1_f, Vy_f, &
  						curvV_half1_f, curvV_f, &
  						theta_dot_half1,theta_dot, &
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
	call make_prediction(slipX_half2,slipX_half1,   &
						slipX_half2_f,slipX_half1_f, &
						slipY_half2,slipY_half1,   &
						slipY_half2_f,slipY_half1_f, &
  						theta_half2,theta_half1, &
  						Vx_half2, Vx_half1, &
  						Vx_half2_f, Vx_half1_f, &
  						Vy_half2, Vy_half1, &
  						Vy_half2_f, Vy_half1_f, &
  						curvV_half2_f, curvV_half1_f, &
  						theta_dot_half2,theta_dot_half1, &
  						dt_try,0,status,SigmaNel,SigmaXel,SigmaYel)
  						
  						
					
 
	
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
real(kind=dp)                                      :: error, slip_full, slip_half2
real(kind=dp)                                      :: error1, error2
integer											   :: istop, status,k1,k2



	! Compare the two results
	error = 0._dp
	error1 = 0._dp
	error2 = 0._dp
	do k1=1,nb_element_x
		do  k2=1,nb_element_y
			if (mask(k1,k2) ==1) then

				! Calculate the norm of slip_full and slip_half2
				slip_full = (slipX_full(k1,k2)**2+slipY_full(k1,k2)**2)**(0.5_dp)
				slip_half2 = (slipX_half2(k1,k2)**2+slipY_half2(k1,k2)**2)**(0.5_dp)
				error = max(error,abs((slip_full-slip_half2)/slip_half2))				
				
				! TODO remove this
				error = max(error,abs(exp(theta_full(k1,k2)-theta_half2(k1,k2))-1._dp))				
				error1 = max(error1,abs((slip_full-slip_half2)/slip_half2))
				error2 = max(error2,abs(exp(theta_full(k1,k2)-theta_half2(k1,k2))-1._dp))
				
				
			end if
		end do
	end do
	

	print*,'error:', error
	print*,'error1:', error1
	print*,'error2:', error2
	
	
	
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
subroutine make_prediction(slipX_e, slipX_i,   &
						 slipX_f_e, slipX_f_i, &
						 slipY_e, slipY_i,   &
						 slipY_f_e, slipY_f_i, &
						 theta_e,theta_i, &
						 Vx_e, Vx_i, &
						 Vx_f_e,Vx_f_i, &
						 Vy_e, Vy_i, &
						 Vy_f_e,Vy_f_i, &
						 curvV_f_e,curvV_f_i, &
						 theta_dot_e, theta_dot_i, &
						 dt_try,i_full, status,sigmaN_e, sigmaX_e,sigmaY_e ) 
	use omp_lib
				 
	! IN 			 
	real(kind=dp), dimension(nb_element_x,nb_element_y)        :: slipX_e,slipX_i,slipY_e,slipY_i
	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: slipX_f_e, slipX_f_i,slipY_f_e, slipY_f_i
	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: curvSlip_f_e
	real(kind=dp), dimension(nb_element_x,nb_element_y)         :: theta_e, theta_i
	real(kind=dp), dimension(nb_element_x,nb_element_y)         :: Vx_e,Vx_i,Vy_e,Vy_i
	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: Vx_f_e,Vx_f_i,Vy_f_e,Vy_f_i
	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: curvV_f_e,curvV_f_i
	real(kind=dp), dimension(nb_element_x,nb_element_y)        :: theta_dot_e, theta_dot_i

	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)  :: integral_history_slipX_modeII,integral_history_slipY_modeII
	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)  :: integral_history_slipX_modeIII,integral_history_slipY_modeIII
	complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)      :: integral_history_normal
	
	real(kind=dp), dimension(nb_element_x,nb_element_y) , optional :: sigmaN_e, sigmaX_e,sigmaY_e
	real(kind=dp), dimension(nb_element_x,nb_element_y)         :: sigmaN_temp, sigmaX_temp,sigmaY_temp


	real(kind=dp)								:: dt_try		 
	integer										:: status, i_full
	
	! Internal
	real(kind=dp), dimension(nb_element_x,nb_element_y)			:: tX,tY	 ! tangent vector	
	real(kind=dp), dimension(nb_element_x,nb_element_y)			:: V	     ! norm of V
	real(kind=dp), dimension(nb_element_x,nb_element_y)			:: SigmaT_tot_temp
	integer														:: k1,k2

	
	
	
	real(kind=dp)           :: t_beg,t_end


	! Initialisation
	integral_history_slipX_modeII = 0._dp
	integral_history_slipY_modeII = 0._dp
	integral_history_slipX_modeIII = 0._dp
	integral_history_slipY_modeIII = 0._dp

	integral_history_normal = 0._dp
	
	
	! Make prediction slip and theta
	call make_prediction_theta_slip(slipX_e,slipY_e,theta_e,slipX_f_e,slipY_f_e,curvSlip_f_e,slipX_i,slipY_i,theta_i,slipX_f_i,slipY_f_i,Vx_i,Vy_i,Vx_f_i,Vy_f_i,theta_dot_i,dt_try)
	!
	! Make prediction of normal and shear traction
	call make_prediction_traction(sigmaN_temp, sigmaX_temp, sigmaY_temp,slipX_f_e,slipY_f_e,curvSlip_f_e,Vx_f_i,Vy_f_i,curvV_f_i,integral_history_slipX_modeII,integral_history_slipX_modeIII,integral_history_slipY_modeII,integral_history_slipY_modeIII,integral_history_normal,i_full,1,status)	
	!
	! Calculate total traction without radiation damping term
	call calculate_traction( sigmaX_temp,sigmaY_temp, SigmaT_tot_temp,dt_try)	
	!
	! Calculate tangent vector to the previous total traction
	call calculate_tangent( sigmaX_temp,sigmaY_temp,SigmaT_tot_temp, tX, tY,dt_try)
	!
	! Find predictive values of V and theta	 
	call make_prediction_V_thetaDot(sigmaN_temp, SigmaT_tot_temp,theta_e,Vx_e,Vy_e,Vx_i,Vy_i,theta_dot_e,Vx_f_e,Vy_f_e,curvV_f_e,tX,tY, status)
					
					
					

	!
	! Calculate final prediction of slip and state variable
	!$OMP PARALLEL DO
	do k2=1,nb_element_y
		do k1=1,nb_element_x
			Vx_e(k1,k2) = (Vx_e(k1,k2)+Vx_i(k1,k2))/2._dp
			Vy_e(k1,k2) = (Vy_e(k1,k2)+Vy_i(k1,k2))/2._dp
			theta_dot_e(k1,k2) = (theta_dot_e(k1,k2)+theta_dot_i(k1,k2))/2._dp
		enddo
	enddo
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO
	do k2=1,nb_frequency_y
		do k1=1,nb_frequency_x
			Vx_f_e(k1,k2) = (Vx_f_e(k1,k2)+Vx_f_i(k1,k2))/2._dp ! CHECK IF V_f_e IS OK
			Vy_f_e(k1,k2) = (Vy_f_e(k1,k2)+Vy_f_i(k1,k2))/2._dp ! CHECK IF V_f_e IS OK
		enddo
	enddo
	!$OMP END PARALLEL DO
	

	! Make prediction slip and theta
	call make_prediction_theta_slip(slipX_e,slipY_e,theta_e,slipX_f_e,slipY_f_e,curvSlip_f_e,slipX_i,slipY_i,theta_i,slipX_f_i,slipY_f_i,Vx_e,Vy_e,Vx_f_e,Vy_f_e,theta_dot_e,dt_try)
	!
	! Make final prediction of normal and shear traction	
	call make_prediction_traction(sigmaN_temp, sigmaX_temp,sigmaY_temp,slipX_f_e,slipY_f_e,curvSlip_f_e,Vx_f_e,Vy_f_e,curvV_f_e,integral_history_slipX_modeII,integral_history_slipX_modeIII,integral_history_slipY_modeII,integral_history_slipY_modeIII,integral_history_normal,i_full,0,status)
    !
    ! Calculate total traction without radiation damping term
	call calculate_traction( sigmaX_temp,sigmaY_temp, SigmaT_tot_temp,dt_try)	
	!
	! Calculate tangent vector to the previous total traction
	call calculate_tangent( sigmaX_temp,sigmaY_temp,SigmaT_tot_temp, tX, tY,dt_try)	
	!
	! Calculate V2 and theta2_dot	 
	call make_prediction_V_thetaDot(sigmaN_temp, SigmaT_tot_temp,theta_e,Vx_e,Vy_e,Vx_e,Vy_e,theta_dot_e,Vx_f_e,Vy_f_e,curvV_f_e,tX,tY, status)
	          
     
	                               
	!
	! Eventually copy shear and normal traction
	if (present(sigmaN_e).and.present(sigmaX_e).and.present(sigmaY_e)) then
		sigmaN_e = sigmaN_temp
		sigmaX_e = sigmaX_temp
		sigmaY_e = sigmaY_temp
	end if
	
	
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   make_prediction_theta
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine make_prediction_theta_slip(slipX_e,slipY_e,theta_e,slipX_f_e,slipY_f_e,curvSlip_f_e,slipX_i,slipY_i,theta_i,slipX_f_i,slipY_f_i,Vx_i,Vy_i,Vx_f_i,Vy_f_i,theta_dot_i,dt_try)



! IN
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: slipX_i,Vx_i,slipY_i,Vy_i
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: theta_i,theta_dot_i

complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: slipX_f_i,Vx_f_i,slipY_f_i,Vy_f_i
!
real(kind=dp)								:: dt_try

! OUT
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: slipX_e,slipY_e,curvSlip_e, theta_e
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: tX_slip,tY_slip

complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: slipX_f_e, slipY_f_e, curvSlip_f_e
!
integer										:: Status_fft,k1,k2
!


!$OMP PARALLEL DO
do k2=1,nb_element_y
	do k1=1,nb_element_x
		! Prediction theta
		theta_e(k1,k2) = theta_i(k1,k2) + dt_try*theta_dot_i(k1,k2)*mask(k1,k2)
		!
		! Prediction slip
		slipX_e(k1,k2) = slipX_i(k1,k2) + dt_try*Vx_i(k1,k2)
		slipY_e(k1,k2) = slipY_i(k1,k2) + dt_try*Vy_i(k1,k2)
		!
		! Prediction curvSlip 
		tX_slip(k1,k2) = slipX_e(k1,k2)/(slipX_e(k1,k2)**2+slipY_e(k1,k2)**2)**(0.5_dp)
		tY_slip(k1,k2) = slipY_e(k1,k2)/(slipX_e(k1,k2)**2+slipY_e(k1,k2)**2)**(0.5_dp)
		curvature(k1,k2) = L1(k1,k2)*tX_slip(k1,k2)**2+2*M1(k1,k2)*tX_slip(k1,k2)*tY_slip(k1,k2)+N1(k1,k2)*tY_slip(k1,k2)**2

		curvSlip_e(k1,k2) = curvature(k1,k2)*(slipX_e(k1,k2)**2+slipY_e(k1,k2)**2)**(0.5_dp)
	enddo
enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
do k2=1,nb_frequency_y
	do k1=1,nb_frequency_x	
		slipX_f_e(k1,k2) = slipX_f_i(k1,k2) + dt_try*Vx_f_i(k1,k2)
		slipY_f_e(k1,k2) = slipY_f_i(k1,k2) + dt_try*Vy_f_i(k1,k2)
	enddo
enddo
!$OMP END PARALLEL DO
			
			
call dfftw_execute_dft_r2c(plan_fwd,  curvSlip_e,curvSlip_f_e)


end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   make_prediction_traction
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine make_prediction_traction(sigmaN_temp, sigmaX_temp,sigmaY_temp,slipX_f_e,slipY_f_e,curvSlip_f_e,Vx_f_e,Vy_f_e,curvV_f_e,integral_history_slipX_modeII,integral_history_slipX_modeIII,integral_history_slipY_modeII,integral_history_slipY_modeIII,integral_history_normal,i_full,i_history,status)

! Internal
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)           :: integral_slipX_modeII,integral_slipX_modeIII
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)           :: integral_slipY_modeII,integral_slipY_modeIII
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)           :: integral_shearX,integral_shearY

complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)           :: integral_normal
! Inout
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: sigmaN_temp, sigmaX_temp,sigmaY_temp
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: integral_history_slipX_modeII, integral_history_slipY_modeII
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: integral_history_slipX_modeIII, integral_history_slipY_modeIII

complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: integral_history_normal
integer(kind=dp)        :: n_ev
integer :: k1, k2, id_time,id_slip, i_history, i_full
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)     :: slipX_f_e,slipY_f_e,curvSlip_f_e,Vx_f_e,Vy_f_e,curvV_f_e
integer										:: status, status_fft

!
!$OMP PARALLEL DO
do k2=1,nb_frequency_y
	do k1=1,nb_frequency_x
		! Calculate integral 
		integral_slipX_modeII(k1,k2) = 0._dp
		integral_slipX_modeIII(k1,k2) = 0._dp
		integral_slipY_modeII(k1,k2) = 0._dp
		integral_slipY_modeIII(k1,k2) = 0._dp
		integral_shearX(k1,k2) = 0._dp
		integral_shearY(k1,k2) = 0._dp

		integral_normal(k1,k2) = 0._dp
	enddo
enddo
!$OMP END PARALLEL DO




! If the calculation is fully dynamic
if (quasi_dynamic == 0) then

	! Select full or half step
	if (i_full==1) then
		n_ev = 2 * n_try
	else if (i_full==0) then
		n_ev = n_try
	end if




	!For each mode 
	!$OMP PARALLEL DO
	do k2=1,nb_frequency_y
		do k1=1,nb_frequency_x
	

			! Calculate the integral on past history
			if(i_history==1)then
! 		
! 	
 				if (n_ev+1.le.modes(k1,k2)%nb_time) then
 					do id_time=n_ev+1,modes(k1,k2)%nb_time	
 			
! 						! Calculate id_slip
 						id_slip = modulo(modes(k1,k2)%current_time_id-id_time+n_ev,modes(k1,k2)%nb_time)+1
 					
! 						! Shear
						integral_history_slipX_modeII(k1,k2) = integral_history_slipX_modeII(k1,k2) + &
						modes(k1,k2)%integrand_shear_modeII(id_time)  * modes(k1,k2)%Vx_f(id_slip)
						
						integral_history_slipX_modeIII(k1,k2) = integral_history_slipX_modeIII(k1,k2) + &
						modes(k1,k2)%integrand_shear_modeIII(id_time)  * modes(k1,k2)%Vx_f(id_slip)
					
						integral_history_slipY_modeII(k1,k2) = integral_history_slipY_modeII(k1,k2) + &
						modes(k1,k2)%integrand_shear_modeII(id_time) * modes(k1,k2)%Vy_f(id_slip)
					
						integral_history_slipY_modeIII(k1,k2) = integral_history_slipY_modeIII(k1,k2) + &
						modes(k1,k2)%integrand_shear_modeIII(id_time) * modes(k1,k2)%Vy_f(id_slip)
! 		
! 						! Normal 
! 						integral_history_normal(k1,k2) = integral_history_normal(k1,k2) + &
! 						modes(k1,k2)%integrand_normal(id_time) * modes(k1,k2)%curvV_f(id_slip)
! 		
! 	
 					end do	
 				end if
 			end if
    	
			! Calculate the integral on current history
 			do id_time=1,min(modes(k1,k2)%nb_time,n_ev)
! 		 
! 				! Shear
				integral_slipX_modeII(k1,k2) = integral_slipX_modeII(k1,k2) + &
							modes(k1,k2)%integrand_shear_modeII(id_time) * Vx_f_e(k1,k2) 
						
				integral_slipX_modeIII(k1,k2) = integral_slipX_modeIII(k1,k2) + &
							modes(k1,k2)%integrand_shear_modeIII(id_time) * Vx_f_e(k1,k2) 
						
				integral_slipY_modeII(k1,k2) = integral_slipY_modeII(k1,k2) + &
							modes(k1,k2)%integrand_shear_modeII(id_time) * Vy_f_e(k1,k2) 
						
				integral_slipY_modeIII(k1,k2) = integral_slipY_modeIII(k1,k2) + &
							modes(k1,k2)%integrand_shear_modeIII(id_time) * Vy_f_e(k1,k2) 
! ! 					      
! 			! Normal 
! 			integral_normal(k1,k2) = integral_normal(k1,k2) + &
! 								modes(k1,k2)%integrand_normal(id_time) * curvV_f_e(k1,k2) 
! 
 			end do


		end do
	end do
	!$OMP END PARALLEL DO
end if
		
! Add static term and history term
!$OMP PARALLEL DO
do k2=1,nb_frequency_y
	do k1=1,nb_frequency_x
		! Shear  
			integral_shearX(k1,k2) =   -mu/Ktot(k1,k2)*Kx(k1,k2)**2*((1-cs**2/cp**2)*slipX_f_e(k1,k2) + integral_slipX_modeII(k1,k2) + integral_history_slipX_modeII(k1,k2)) &
	  	            				   -mu/Ktot(k1,k2)*Ky(k1,k2)**2*(slipX_f_e(k1,k2)/2 - integral_slipX_modeIII(k1,k2) - integral_history_slipX_modeIII(k1,k2))              &
	  	            			       -mu/Ktot(k1,k2)*Kx(k1,k2)*Ky(k1,k2)*((1-cs**2/cp**2)*slipY_f_e(k1,k2)-slipY_f_e(k1,k2)/2+ integral_slipY_modeII(k1,k2) + integral_history_slipY_modeII(k1,k2)+ integral_slipY_modeIII(k1,k2) + integral_history_slipY_modeIII(k1,k2))
	  	            
	  	          
	  	
	  	                  
			integral_shearY(k1,k2) =   -mu/Ktot(k1,k2)*Kx(k1,k2)**2*(slipY_f_e(k1,k2)/2 - integral_slipY_modeIII(k1,k2) - integral_history_slipY_modeIII(k1,k2)) &
	  	          					   -mu/Ktot(k1,k2)*Ky(k1,k2)**2*((1-cs**2/cp**2)*slipY_f_e(k1,k2) + integral_slipY_modeII(k1,k2) + integral_history_slipY_modeII(k1,k2)) &
	  	           				      -mu/Ktot(k1,k2)*Kx(k1,k2)*Ky(k1,k2)*((1-cs**2/cp**2)*slipX_f_e(k1,k2)-slipX_f_e(k1,k2)/2 + integral_slipX_modeII(k1,k2) + integral_history_slipX_modeII(k1,k2) + integral_slipX_modeIII(k1,k2) + integral_history_slipX_modeIII(k1,k2))


			! Normal (static only)
			integral_normal(k1,k2) = mu*cmplx(0._dp,1._dp)*(1-cs**2/cp**2)*curvSlip_f_e(k1,k2)
	
	enddo
enddo
!$OMP END PARALLEL DO


! Make sure that the 0 frequency is 0
integral_shearX(1,1) = 0._dp
integral_shearY(1,1) = 0._dp
integral_normal(1,1) = 0._dp
!
! Inverse fourier transform of the shear traction
call dfftw_execute_dft_c2r(plan_bwd,integral_shearX/(nb_element_x*nb_element_y),SigmaX_temp )
call dfftw_execute_dft_c2r(plan_bwd, integral_shearY/(nb_element_x*nb_element_y),SigmaY_temp)
call dfftw_execute_dft_c2r(plan_bwd,integral_normal/(nb_element_x*nb_element_y),SigmaN_temp )
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
!   Calculate total traction norm (without radiation damping) on the fault
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine calculate_traction( sigmaX_temp,sigmaY_temp, SigmaT_tot_temp,dt_try)
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: sigmaT_tot_temp, sigmaX_temp,sigmaY_temp
real(kind=dp)                                              :: dt_try
integer(kind=dp)										   :: k1,k2


!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
! real(kind=dp)                                              :: perturb_time
! real(kind=dp), dimension(nb_element_x,nb_element_y)        :: perturb_tot


! if (time+dt_try <=1._dp) then
! 	perturb_time = exp((time+dt_try-1._dp)**2/((time+dt_try)*(time+dt_try-2._dp)))
! else
! 	perturb_time = 1._dp
! end if
! perturb_tot = 2000000._dp*perturb_time*perturb_space
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%




!$OMP PARALLEL DO
do k2=1,nb_element_y
	do k1=1,nb_element_x
		! Calculate the norm of the force acting on the fault
		SigmaT_tot_temp(k1,k2) = ((sigmaX0(k1,k2)+sigmaX_temp(k1,k2) + shearX_loading_dot(k1,k2)*(time+dt_try))**2+(sigmaY0(k1,k2)+sigmaY_temp(k1,k2) + shearY_loading_dot(k1,k2)*(time+dt_try))**2)**(0.5_dp)

! 		SigmaT_tot_temp(k1,k2) = ((sigmaX0(k1,k2)+sigmaX_temp(k1,k2) +perturb_tot(k1,k2) + shearX_loading_dot(k1,k2)*(time+dt_try))**2+(sigmaY0(k1,k2)+sigmaY_temp(k1,k2) + shearY_loading_dot(k1,k2)*(time+dt_try))**2)**(0.5_dp)
	enddo
enddo
!$OMP END PARALLEL DO


end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Calculate tangent direction to the total traction on the fault (without radiation damping)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine calculate_tangent( sigmaX_temp,sigmaY_temp,SigmaT_tot_temp, tX_temp, tY_temp,dt_try)
real(kind=dp), dimension(nb_element_x,nb_element_y)        :: sigmaT_tot_temp, sigmaX_temp,sigmaY_temp, tX_temp, tY_temp
real(kind=dp) 											   :: dt_try
integer(kind=dp)										   :: k1,k2



!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
! real(kind=dp)                                              :: perturb_time
! real(kind=dp), dimension(nb_element_x,nb_element_y)        :: perturb_tot
! 
! 
! if (time+dt_try <=1._dp) then
! 	perturb_time = exp((time+dt_try-1._dp)**2/((time+dt_try)*(time+dt_try-2._dp)))
! else
! 	perturb_time = 1._dp
! end if
! perturb_tot = 2000000._dp*perturb_time*perturb_space
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

!$OMP PARALLEL DO
do k2=1,nb_element_y
	do k1=1,nb_element_x
		! Calculate the norm of the force acting on the fault
		!tX_temp(k1,k2) = (perturb_tot(k1,k2)+sigmaX0(k1,k2)+ sigmaX_temp(k1,k2) + shearX_loading_dot(k1,k2)*(time+dt_try))/sigmaT_tot_temp(k1,k2)*mask(k1,k2)
		tX_temp(k1,k2) = (sigmaX0(k1,k2)+ sigmaX_temp(k1,k2) + shearX_loading_dot(k1,k2)*(time+dt_try))/sigmaT_tot_temp(k1,k2)*mask(k1,k2)
		tY_temp(k1,k2) = (sigmaY0(k1,k2)+ sigmaY_temp(k1,k2) + shearY_loading_dot(k1,k2)*(time+dt_try))/sigmaT_tot_temp(k1,k2)*mask(k1,k2)
	enddo
enddo
!$OMP END PARALLEL DO
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   find_predicted_slip
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine make_prediction_V_thetaDot(sigmaN_temp, sigmaT_temp,theta_e,Vx_e,Vy_e,Vx_i,Vy_i,theta_dot_e,Vx_f_e,Vy_f_e,curvV_f_e, tX,tY, status)
real(kind=dp), dimension(nb_element_x,nb_element_y)        ::sigmaN_temp, sigmaT_temp
real(kind=dp), dimension(nb_element_x,nb_element_y)       :: theta_e,Vx_e,Vx_i,Vy_e,Vy_i,curvV_e
real(kind=dp), dimension(nb_element_x,nb_element_y)      :: theta_dot_e, tX, tY
complex(kind=dp), dimension(nb_frequency_x,nb_frequency_y)    :: Vx_f_e,Vy_f_e,curvV_f_e

! Newton-raphson
real(kind=dp), dimension(nb_element_x,nb_element_y)       :: res, res_dot,res_dot2
real(kind=dp), dimension(nb_element_x,nb_element_y)       :: coeff1

! Internal
integer   :: i, status, status_fft,k1,k2
real(kind=dp), dimension(nb_element_x,nb_element_y)       :: V_e


! Initialise Newton raphson
i=0
!
! Initialise
!$OMP PARALLEL DO
do k2=1,nb_element_y
	do k1=1,nb_element_x
		Vx_e(k1,k2) = Vx_i(k1,k2)
		Vy_e(k1,k2) = Vy_i(k1,k2)
		V_e(k1,k2) = (Vx_i(k1,k2)**2+Vy_i(k1,k2)**2)**(0.5_dp)
		!
		! Calculate coeff1
		coeff1(k1,k2) = exp(f0(k1,k2)/a(k1,k2)+(b(k1,k2)/a(k1,k2))*theta_e(k1,k2))*mask(k1,k2)
		!
		! Calculate res and res_dot
		res(k1,k2) = 10.0*nr_accuracy/V_e(k1,k2)
		res_dot(k1,k2) =  1._dp
	enddo
enddo
!$OMP END PARALLEL DO	

!
! Use Newton raphson 
do while (maxval(abs((res/res_dot)*mask/V_e))>nr_accuracy)


	! Calculate new function to minimize and its derivative
	!$OMP PARALLEL DO
	do k2=1,nb_element_y
		do k1=1,nb_element_x
			res(k1,k2) = ( sigmaT_temp(k1,k2) - mu/(2._dp*cs)*V_e(k1,k2) & 
	   		+(sigmaN0(k1,k2)+sigmaN_temp(k1,k2))*a(k1,k2)*asinh(V_e(k1,k2)/(2._dp*V0(k1,k2))*coeff1(k1,k2)))*mask(k1,k2)
	
			res_dot(k1,k2) = - mu/(2._dp*cs) + (SigmaN0(k1,k2)+SigmaN_temp(k1,k2))*a(k1,k2)/sqrt(1._dp/coeff1(k1,k2)**2+(V_e(k1,k2)/(2*V0(k1,k2)))**2)/(2._dp*V0(k1,k2))
		!	
			! Update V1
			V_e(k1,k2) = V_e(k1,k2) -  (res(k1,k2)/res_dot(k1,k2))*mask(k1,k2)

		enddo
	enddo
	!$OMP END PARALLEL DO
	!
	! Check convergence
	i=i+1
	!
	! In case Newton-Raphson does not converge
	if (i > 20) then
		print*,'number of iteration (NF):',i
	    print*,'maxval(abs((f/f_dot)*mask))',maxval(abs((res/res_dot)*mask/V_e))
	    print*,'res',minval(res*mask),maxval(res*mask)
	    print*,'res_dot',minval(max(res_dot,1._dp-mask)),maxval(res*mask)
	    print*,'V_e',minval(V_e),maxval(V_e)
	end if 
	
	
	if (i > 40) then
		print*,'Newton-Raphson (1) did not converge'
		status = 1
		exit
	end if   
	
end do
! 
!$OMP PARALLEL DO
do k2=1,nb_element_y
	do k1=1,nb_element_x
		! Compute Vx and Vy (same direction as traction)
		Vx_e(k1,k2) = tX(k1,k2)*V_e(k1,k2)*mask(k1,k2) + (1-mask(k1,k2))*Vx_mask
		Vy_e(k1,k2) = tY(k1,k2)*V_e(k1,k2)*mask(k1,k2) + (1-mask(k1,k2))*Vy_mask
		!
		! Calculate theta_dot
		theta_dot_e(k1,k2) = V0(k1,k2)/Dc(k1,k2)*exp(-theta(k1,k2)) - (V_e(k1,k2))/Dc(k1,k2)
	enddo
enddo
!$OMP END PARALLEL DO
!
! Calculate V_f_e and curvV_f_e
call dfftw_execute_dft_r2c(plan_fwd,  Vx_e,Vx_f_e )
call dfftw_execute_dft_r2c(plan_fwd, Vy_e,Vy_f_e)
!curvV_e = curvature*(Vx_e**2+Vy_e**2)**(0.5_dp)
!status_fft = DftiComputeForward( Desc_fwd, curvV_e(1,:),curvV_f_e(1,:) )


end subroutine




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   save data 
!
subroutine save_data()   
   if (modulo(time_id-1,stride_time) == 0) then 

        it_count_save = it_count_save + 1
        Vx_save(it_count_save,:,:) = Vx
        Vy_save(it_count_save,:,:) = Vy
        slipX_save(it_count_save,:,:) = slipX
        slipY_save(it_count_save,:,:) = slipY
        theta_save(it_count_save,:,:) = theta
        time_save(it_count_save) = time 
        tractionNel_save(it_count_save,:,:) = sigmaNel
        tractionXel_save(it_count_save,:,:) = sigmaXel!-mu/(2._dp*cs)*Vx
        tractionYel_save(it_count_save,:,:) = sigmaYel!-mu/(2._dp*cs)*Vy

    end if
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Write data in a file
!
subroutine write_data()       
character(len=30)   :: num2str
integer :: error
integer(HSIZE_T) :: data_dims_array(3), data_dims_time(1)
integer(HID_T)   :: file_id, dspace_id_Vx,dspace_id_time, dset_id_Vx,dset_id_time
integer(HID_T)   :: dspace_id_Vy,dset_id_Vy,dspace_id_theta,dset_id_theta
integer(HID_T)   :: dspace_id_tractionXel,dset_id_tractionXel,dspace_id_tractionYel,dset_id_tractionYel
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
	data_dims_array(2) = nb_element_x
	data_dims_array(3) = nb_element_y
	data_dims_time = freq_writing_file
	!
	! Create dataspace (the dataset is next) "dspace_id" is returned
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_Vx, error)
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_Vy, error)
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_theta, error)
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_tractionXel, error)
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_tractionYel, error)
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_tractionNel, error)
	CALL h5screate_simple_f(1, data_dims_time, dspace_id_time, error)

	
	! Create dataset with default properties "dset_id" is returned
	CALL h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, dspace_id_Vx, dset_id_Vx, error)
	CALL h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, dspace_id_Vy, dset_id_Vy, error)
	CALL h5dcreate_f(file_id, 'theta', H5T_NATIVE_DOUBLE, dspace_id_theta, dset_id_theta, error)
	CALL h5dcreate_f(file_id, 'tractionXel', H5T_NATIVE_DOUBLE, dspace_id_tractionXel, dset_id_tractionXel, error)
	CALL h5dcreate_f(file_id, 'tractionYel', H5T_NATIVE_DOUBLE, dspace_id_tractionYel, dset_id_tractionYel, error)
	CALL h5dcreate_f(file_id, 'tractionNel', H5T_NATIVE_DOUBLE, dspace_id_tractionNel, dset_id_tractionNel, error)
	CALL h5dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, dspace_id_time, dset_id_time, error)

                     
	! Write dataset 
	CALL h5dwrite_f(dset_id_Vx, H5T_NATIVE_DOUBLE, Vx_save(1:it_count_save,1:nb_element_x,1:nb_element_y), data_dims_array, error)
	CALL h5dwrite_f(dset_id_Vy, H5T_NATIVE_DOUBLE, Vy_save(1:it_count_save,1:nb_element_x,1:nb_element_y), data_dims_array, error)
	CALL h5dwrite_f(dset_id_theta, H5T_NATIVE_DOUBLE, theta_save(1:it_count_save,1:nb_element_x,1:nb_element_y), data_dims_array, error)
	CALL h5dwrite_f(dset_id_tractionXel, H5T_NATIVE_DOUBLE, tractionXel_save(1:it_count_save,1:nb_element_x,1:nb_element_y), data_dims_array, error)
	CALL h5dwrite_f(dset_id_tractionYel, H5T_NATIVE_DOUBLE, tractionYel_save(1:it_count_save,1:nb_element_x,1:nb_element_y), data_dims_array, error)
	CALL h5dwrite_f(dset_id_tractionNel, H5T_NATIVE_DOUBLE, tractionNel_save(1:it_count_save,1:nb_element_x,1:nb_element_y), data_dims_array, error)
	CALL h5dwrite_f(dset_id_time, H5T_NATIVE_DOUBLE, time_save(1:it_count_save), data_dims_time, error)
                      
	! Close access to dataset 1
	CALL h5dclose_f(dset_id_Vx, error)
	CALL h5dclose_f(dset_id_Vy, error)
	CALL h5dclose_f(dset_id_theta, error)
	CALL h5dclose_f(dset_id_tractionXel, error)
	CALL h5dclose_f(dset_id_tractionYel, error)
	CALL h5dclose_f(dset_id_tractionNel, error)
	CALL h5dclose_f(dset_id_time, error)

	! Close access to data space 1
	CALL h5sclose_f(dspace_id_Vx, error)
	CALL h5sclose_f(dspace_id_Vy, error)
	CALL h5sclose_f(dspace_id_theta, error)
	CALL h5sclose_f(dspace_id_tractionXel, error)
	CALL h5sclose_f(dspace_id_tractionYel, error)
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
integer :: end_id,k1,k2



! Overwrite history
do k1=1,nb_frequency_x
	do k2=1,nb_frequency_y

		! No problem to overwrite
		if (modes_temp(k1,k2)%current_time_id+n_try .le. modes_temp(k1,k2)%nb_time) then 
	
	
			! Overwrite
			modes_temp(k1,k2)%Vx_f(modes_temp(k1,k2)%current_time_id+1:modes_temp(k1,k2)%current_time_id+n_try) =  (Vx_half1_f(k1,k2)+	modes_temp(k1,k2)%Vx_f(modes_temp(k1,k2)%current_time_id))/2._dp
			modes_temp(k1,k2)%Vy_f(modes_temp(k1,k2)%current_time_id+1:modes_temp(k1,k2)%current_time_id+n_try) =  (Vy_half1_f(k1,k2)+	modes_temp(k1,k2)%Vy_f(modes_temp(k1,k2)%current_time_id))/2._dp

			!modes_temp(k1,k2)%curvV_f(modes_temp(k1,k2)%current_time_id+1:modes_temp(k1,k2)%current_time_id+n_try) =  (curvV_half1_f(k1,k2)+modes_temp(k1,k2)%curvV_f(modes_temp(k1,k2)%current_time_id))/2._dp
		
		
		
		
			! Set new current time
			modes_temp(k1,k2)%current_time_id = modes_temp(k1,k2)%current_time_id + n_try	


			!
		! Start exactly at the beginning of the file
		else if (modes_temp(k1,k2)%current_time_id+1 .gt. modes_temp(k1,k2)%nb_time) then
			!
			! Calculate new end at the begining of the array
			end_id = modulo(modes_temp(k1,k2)%current_time_id+n_try,modes_temp(k1,k2)%nb_time)
			
			! Overwrite data at the beginning
			modes_temp(k1,k2)%Vx_f(1:end_id) =  (Vx_half1_f(k1,k2)+modes_temp(k1,k2)%Vx_f(modes_temp(k1,k2)%current_time_id))/2._dp
			modes_temp(k1,k2)%Vy_f(1:end_id) =  (Vy_half1_f(k1,k2)+modes_temp(k1,k2)%Vy_f(modes_temp(k1,k2)%current_time_id))/2._dp
			!modes_temp(k1,k2)%curvV_f(1:end_id) =  (curvV_half1_f(k1,k2)+modes_temp(k1,k2)%curvV_f(modes_temp(k1,k2)%current_time_id))/2._dp
		
	
			! Set new current time
			modes_temp(k1,k2)%current_time_id = end_id
		! All other cases
		else 
			! Calculate new end at the beginning of the array
			end_id = modulo(modes_temp(k1,k2)%current_time_id+n_try,modes_temp(k1,k2)%nb_time)
			
			! Overwrite data at the end
			modes_temp(k1,k2)%Vx_f(1:end_id) =  (Vx_half1_f(k1,k2)+modes_temp(k1,k2)%Vx_f(modes_temp(k1,k2)%current_time_id))/2._dp
			modes_temp(k1,k2)%Vy_f(1:end_id) =  (Vy_half1_f(k1,k2)+modes_temp(k1,k2)%Vy_f(modes_temp(k1,k2)%current_time_id))/2._dp

			!modes_temp(k1,k2)%curvV_f(1:end_id) =  (curvV_half1_f(k1,k2)+modes_temp(k1,k2)%curvV_f(modes_temp(k1,k2)%current_time_id))/2._dp

			! Overwrite data at the beginning
			modes_temp(k1,k2)%Vx_f(modes_temp(k1,k2)%current_time_id+1:modes_temp(k1,k2)%nb_time) =  (Vx_half1_f(k1,k2)+modes_temp(k1,k2)%Vx_f(modes_temp(k1,k2)%current_time_id))/2._dp
			modes_temp(k1,k2)%Vy_f(modes_temp(k1,k2)%current_time_id+1:modes_temp(k1,k2)%nb_time) =  (Vy_half1_f(k1,k2)+modes_temp(k1,k2)%Vy_f(modes_temp(k1,k2)%current_time_id))/2._dp
			!modes_temp(k1,k2)%curvV_f(modes_temp(k1,k2)%current_time_id+1:modes_temp(k1,k2)%nb_time) =  (curvV_half1_f(k1,k2)+modes_temp(k1,k2)%curvV_f(modes_temp(k1,k2)%current_time_id))/2._dp

			! Set new current time
			modes_temp(k1,k2)%current_time_id = end_id
		

		end if 
	end do
end do


end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   update_values
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_values()
integer :: end_id,k1,k2


slipX  = slipX_half2
slipX_f = slipX_half2_f
slipY  = slipY_half2
slipY_f = slipY_half2_f
theta = theta_half2
theta_dot = theta_dot_half2

! Update V
Vx = Vx_half2
Vx_f = Vx_half2_f
Vy = Vy_half2
Vy_f = Vy_half2_f
curvV_f = curvV_half2_f

if (quasi_dynamic ==0) then
	! Update to half1
	do k1=1,nb_frequency_x
		do k2=1,nb_frequency_y
  			modes(k1,k2)%current_time_id = modes_temp(k1,k2)%current_time_id
  			modes(k1,k2)%Vx_f = modes_temp(k1,k2)%Vx_f
  			modes(k1,k2)%Vy_f = modes_temp(k1,k2)%Vy_f
  			!modes(k1,k2)%curvV_f = modes_temp(k1,k2)%curvV_f
  		end do
	end do


	! Update to half2
	do k1=1,nb_frequency_x
		do k2=1,nb_frequency_y
			! No problem to overwrite
			if (modes(k1,k2)%current_time_id+n_try .le. modes(k1,k2)%nb_time) then 

				! Overwrite
				modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id+1:modes(k1,k2)%current_time_id+n_try) = (Vx_f(k1,k2)+modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id))/2._dp
				modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id+1:modes(k1,k2)%current_time_id+n_try) = (Vy_f(k1,k2)+modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id))/2._dp
				!modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id+1:modes(k1,k2)%current_time_id+n_try) = (curvV_f(k1,k2)+modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id))/2._dp


				! Set new current time
				modes(k1,k2)%current_time_id = modes(k1,k2)%current_time_id + n_try	
				!
				! Start exactly at the beginning of the file
			else if (modes(k1,k2)%current_time_id+1 .gt. modes(k1,k2)%nb_time) then
				!
				! Calculate new end at the begining of the array
				end_id = modulo(modes(k1,k2)%current_time_id+n_try,modes(k1,k2)%nb_time)
				!
				! Overwrite data at the beginning			
				modes(k1,k2)%Vx_f(1:end_id) = (Vx_f(k1,k2)+modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id))/2._dp
				modes(k1,k2)%Vy_f(1:end_id) = (Vy_f(k1,k2)+modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id))/2._dp
			

				!modes(k1,k2)%curvV_f(1:end_id) = (curvV_f(k1,k2)+modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id))/2._dp
				! Set new current time
				modes(k1,k2)%current_time_id = end_id
			! All other cases
			else 

				! Calculate new end at the beginning of the array
				end_id = modulo(modes(k1,k2)%current_time_id+n_try,modes(k1,k2)%nb_time)
				! Overwrite data at the beginning
				modes(k1,k2)%Vx_f(1:end_id) = (Vx_f(k1,k2)+modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id))/2._dp
				modes(k1,k2)%Vy_f(1:end_id) = (Vy_f(k1,k2)+modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id))/2._dp
				!modes(k1,k2)%curvV_f(1:end_id) = (curvV_f(k1,k2)+modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id))/2._dp

				! Overwrite data at the beginning
				modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id+1:modes(k1,k2)%nb_time) = (Vx_f(k1,k2)+modes(k1,k2)%Vx_f(modes(k1,k2)%current_time_id))/2._dp
				modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id+1:modes(k1,k2)%nb_time) = (Vy_f(k1,k2)+modes(k1,k2)%Vy_f(modes(k1,k2)%current_time_id))/2._dp
				!modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id+1:modes(k1,k2)%nb_time) = (curvV_f(k1,k2)+modes(k1,k2)%curvV_f(modes(k1,k2)%current_time_id))/2._dp


				! Set new current time
				modes(k1,k2)%current_time_id = end_id
			end if 
		end do
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
allocate(slipX_half1(nb_element_x,nb_element_y))
allocate(slipX_half2(nb_element_x,nb_element_y))
allocate(slipX_full(nb_element_x,nb_element_y))

allocate(slipX_half1_f(nb_frequency_x,nb_frequency_y))
allocate(slipX_half2_f(nb_frequency_x,nb_frequency_y))
allocate(slipX_full_f(nb_frequency_x,nb_frequency_y))

allocate(slipY_half1(nb_element_x,nb_element_y))
allocate(slipY_half2(nb_element_x,nb_element_y))
allocate(slipY_full(nb_element_x,nb_element_y))

allocate(slipY_half1_f(nb_frequency_x,nb_frequency_y))
allocate(slipY_half2_f(nb_frequency_x,nb_frequency_y))
allocate(slipY_full_f(nb_frequency_x,nb_frequency_y))


allocate(theta_half1(nb_element_x,nb_element_y))
allocate(theta_half2(nb_element_x,nb_element_y))
allocate(theta_full(nb_element_x,nb_element_y))

allocate(Vx_half1(nb_element_x,nb_element_y))
allocate(Vx_half2(nb_element_x,nb_element_y))
allocate(Vx_full(nb_element_x,nb_element_y))

allocate(Vx_half1_f(nb_frequency_x,nb_frequency_y))
allocate(Vx_half2_f(nb_frequency_x,nb_frequency_y))
allocate(Vx_full_f(nb_frequency_x,nb_frequency_y))

allocate(Vy_half1(nb_element_x,nb_element_y))
allocate(Vy_half2(nb_element_x,nb_element_y))
allocate(Vy_full(nb_element_x,nb_element_y))

allocate(Vy_half1_f(nb_frequency_x,nb_frequency_y))
allocate(Vy_half2_f(nb_frequency_x,nb_frequency_y))
allocate(Vy_full_f(nb_frequency_x,nb_frequency_y))

allocate(curvV_half1_f(nb_frequency_x,nb_frequency_y))
allocate(curvV_half2_f(nb_frequency_x,nb_frequency_y))
allocate(curvV_full_f(nb_frequency_x,nb_frequency_y))

allocate(theta_dot_half1(nb_element_x,nb_element_y))
allocate(theta_dot_half2(nb_element_x,nb_element_y))
allocate(theta_dot_full(nb_element_x,nb_element_y))

!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
! allocate(perturb_space(nb_element_x,nb_element_y))
!#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%

!
! Save data
allocate(time_save(freq_writing_file))

allocate(Vx_save(freq_writing_file,nb_element_x,nb_element_y))
allocate(Vy_save(freq_writing_file,nb_element_x,nb_element_y))

allocate(slipX_save(freq_writing_file,nb_element_x,nb_element_y))
allocate(slipY_save(freq_writing_file,nb_element_x,nb_element_y))

allocate(tractionNel_save(freq_writing_file,nb_element_x,nb_element_y))
allocate(tractionXel_save(freq_writing_file,nb_element_x,nb_element_y))
allocate(tractionYel_save(freq_writing_file,nb_element_x,nb_element_y))


allocate(theta_save(freq_writing_file,nb_element_x,nb_element_y))
!
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Deallocate arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine deallocate_arrays()
integer :: Status_fft

deallocate(slipX_half1)
deallocate(slipX_half2)
deallocate(slipX_full)
deallocate(slipY_half1)
deallocate(slipY_half2)
deallocate(slipY_full)
deallocate(theta_half1)
deallocate(theta_half2)
deallocate(theta_full)
deallocate(Vx_half1)
deallocate(Vx_half2)
deallocate(Vx_full)
deallocate(Vy_half1)
deallocate(Vy_half2)
deallocate(Vy_full)
deallocate(theta_dot_half1)
deallocate(theta_dot_half2)
deallocate(theta_dot_full)

! Free Descriptor
call dfftw_destroy_plan(plan_bwd)
call dfftw_destroy_plan(plan_fwd)

end subroutine

end module

!%%
