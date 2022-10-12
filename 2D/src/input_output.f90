!
!
!
!
!
!
!
module input_output
use variables
use default_values, only:assign_default_values
!
implicit none
!

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Main routine that read config, data, and then allocate everything
!
subroutine config()
!
! Get path to config file and geometry file 
call get_location_config()
!
! Read config file
call read_config()
!
! Allocate list
call allocate_arrays()
!
! Read 
call read_ini()
!
! Show simulation parameter
call show_variables()

end subroutine 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read location of config.in and geometry.in from command line argument
!
subroutine get_location_config()
integer         :: narg,stat1
!
!   Get number of argument
narg = command_argument_count()
!
!   If there is one argument, get it
if(narg==1) then
        call get_command_argument(1,path_problem)
else
    print*,'Error: please specify a directory containing config.in, geometry.in and curvature.in'
    stop
end if
!
!   Build path to config_file
config_file = trim(path_problem)//'config.in'
!
! Try to open it
OPEN(UNIT=1, FILE=config_file,ACTION='read',IOSTAT=stat1,STATUS='OLD')
close(1)
if(stat1.ne.0)then
    print*,'Error: please check the path to config.in'
    stop
end if
end subroutine





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read the main configuration file
!
subroutine read_config()
    namelist /main/ simulation_name
    namelist /friction/ friction_law, state_evolution, flash_heating, Vw, fw, n_fh 
    namelist /material_and_loading/ mu, cp, cs,V_mask
    namelist /geometry/ L_fault, nb_element
    namelist /output/               stride_time,freq_writing_file
    namelist /hyperparameters/      rho_c, eta_w, beta_min, nr_accuracy, tol_solver,max_it, omp_threads, quasi_dynamic 
!
!   Read variables from the input file
    open(unit=1,file=trim(config_file))
        read(1,nml=main) 
        read(1,nml=friction)
        read(1,nml=material_and_loading)
        read(1,nml=geometry)
        read(1,nml=output)
        read(1,nml=hyperparameters)
    close(1)
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Allocate all the arrays
!
subroutine allocate_arrays()
! Friction
allocate(a(nb_element))
allocate(b(nb_element))
allocate(Dc(nb_element))
allocate(f0(nb_element))
allocate(V0(nb_element))

! Loading
allocate(normal_loading_dot(nb_element))
allocate(shear_loading_dot(nb_element))

! Initial parameters
allocate(V(nb_element))
allocate(theta(nb_element))
allocate(sigmaN0(nb_element))

! Geometry
allocate(height(nb_element))
allocate(mask(nb_element))
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Read all the given nodes
!
subroutine read_ini()
! 
! Read all the values of friction
call read_from_file1(trim(path_problem)//'friction/a.in',a)
call read_from_file1(trim(path_problem)//'friction/b.in',b)
call read_from_file1(trim(path_problem)//'friction/Dc.in',Dc)
call read_from_file1(trim(path_problem)//'friction/f0.in',f0)
call read_from_file1(trim(path_problem)//'friction/V0.in',V0)
!
! Read all the values of loading 
call read_from_file1(trim(path_problem)//'loading/normal_loading_dot.in',normal_loading_dot)
call read_from_file1(trim(path_problem)//'loading/shear_loading_dot.in',shear_loading_dot)
!
! Read all the values of initial_parameters
call read_from_file1(trim(path_problem)//'initial_parameters/V.in',V)
call read_from_file1(trim(path_problem)//'initial_parameters/theta.in',theta)
call read_from_file1(trim(path_problem)//'initial_parameters/sigmaN.in',sigmaN0)
!
! Real all the values of geometry
call read_from_file1(trim(path_problem)//'geometry/height.in',height)
call read_from_file2(trim(path_problem)//'geometry/mask.in',mask)


end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Read all the given nodes
!
subroutine read_from_file1(path,variable)
integer                               :: iflag,element_count
real(kind=dp), dimension(nb_element)  :: variable
character(len=*)                    :: path
print*,'path: ', trim(path)
!
! Open file
open(unit=1,file=trim(path))
! Initializing 
iflag = 0
! Loop over faults and nodes 
do element_count=1,nb_element
    read(unit=1,FMT=*,iostat=iflag) variable(element_count)
    if(iflag/=0) then
        print*,'Error: check your data format: ', path
        stop
    end if
end do

! Close file
close(1)
end subroutine
subroutine read_from_file2(path,variable)
integer                               :: iflag,element_count
integer, dimension(nb_element)        :: variable
character(len=*)                      :: path
!
! Open file
open(unit=1,file=trim(path))
! Initializing 
iflag = 0
! Loop over faults and nodes 
do element_count=1,nb_element
    read(unit=1,FMT=*,iostat=iflag) variable(element_count)
    if(iflag/=0) then
        print*,'Error: check your data format: ', path
        stop
    end if
end do

! Close file
close(1)
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       Print every parameter at screen
!
subroutine show_variables()
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%% Modeling %%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*
print*,"You are modeling one fault with a ",trim(friction_law)," friction law."
print*,"State evolution law: ",trim(state_evolution)
print*
print*,"%%%%%%%%%%%%%% Material properties %%%%%%%%%%%%%"
print*
print*,"mu = ", mu
print*,"Cs = ", cs
print*,"Cp = ", cp
print*
if (flash_heating == 1) then
print*,"%%%%%%%%%%%%%% Flash Heating %%%%%%%%%%%%%%%%%%%%"
print*
print*,"fw = ", fw
print*,"Vw = ", Vw
print*,"n = ", n_fh
print*
end if
print*,"%%%%%%%%%%%%%%% Hyperparameter %%%%%%%%%%%%%%%%%"
print*
print*,"eta_w = ", eta_w
print*,"rho_c = ", rho_c
print*,"beta_min = ", beta_min
print*,"nr_accuracy = ", nr_accuracy
print*,"tol_solver = ", tol_solver
print*,"quasi_dynamic = ", quasi_dynamic
print*
print*,"%%%%%%%%%%%%%%% Time to compute %%%%%%%%%%%%%%%%"
print*
print*
print*
print*,'Nb of elements = ', nb_element
print*
print*
print*,'Please press ENTER'


! 
! print*,'a'
! print*,a
! print*,'b'
! print*,b
! print*,'Dc'
! print*,Dc
! print*,'f0'
! print*,f0
! print*,'V0'
! print*,V0
! 
! print*,'normal_loading_dot'
! print*,normal_loading_dot
! print*,'shear_loading_dot'
! print*,shear_loading_dot
! 
! print*,'curvature'
! print*,curvature
! print*,'mask'
! print*,mask
! 
! print*,'V'
! print*,V
! print*,'theta'
! print*,theta
! print*,'sigmaN'
! print*,sigmaN





end subroutine

end module
