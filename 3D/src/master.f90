!
!
!
!
!
!
program earthquake_cycle
use input_output
use initialisation
use solver
use omp_lib
real(kind=dp)           :: t_beg,t_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%% Reading data %%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print*, "Reading configuration file" 
call config()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%% Initialisation %%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print*, "Initialising arrays"
call initialise()
!
!call show_variables()
!print*, "	COMPUTATION STARTED"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%% Running cycles %%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!!t_beg = omp_get_wtime()
call solve()
!t_end = omp_get_wtime()
!write(6,"(A,F16.4,A)"), " 	COMPUTATION FINISHED IN",t_end-t_beg,'(s)'
!print*,''
!print*,'	/////////////////////////////////////////////////'
!print*,''

end program