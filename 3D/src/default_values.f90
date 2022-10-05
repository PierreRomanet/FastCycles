!
!
!
module default_values
use variables
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine assign_default_values()
integer                 :: k
!
!   Default values
!   Main
simulation_name = 'test/'
!
! Friction law
friction_law = 'RateState'
state_evolution = 'aging'

!
!   material and loading
mu = 30e9_dp
cs = 3000.0_dp
cp = 5000.0_dp
!
!
! Simulation parameters
max_it = 10
!
! Output
stride_time = 1                       ! Save data each stride_time time step
freq_writing_file = 10000             ! Frequence at which we write the output (every ~time step)
!
! Specific variables to Lapusta & al 2000
rho_c = 100
eta_w = 1
beta_min = 0.5d0
gamma = 2.0
nr_accuracy = 1e-9
tol_solver = 1e-6
quasi_dynamic = 1
!
end subroutine
end module