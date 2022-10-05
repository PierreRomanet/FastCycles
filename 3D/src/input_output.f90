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
use hdf5
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
print*,omp_threads
call OMP_SET_NUM_THREADS(omp_threads)
call omp_set_nested(1)

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
    namelist /friction/ friction_law, state_evolution
    namelist /material_and_loading/ mu, cp, cs,Vx_mask,Vy_mask
    namelist /geometry/ L_fault_x, nb_element_x, L_fault_y, nb_element_y
    namelist /output/               stride_time,freq_writing_file
    namelist /hyperparameters/      rho_c, eta_w, beta_min, gamma, nr_accuracy, tol_solver,max_it, omp_threads, quasi_dynamic
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
allocate(a(nb_element_x, nb_element_y))
allocate(b(nb_element_x, nb_element_y))
allocate(Dc(nb_element_x, nb_element_y))
allocate(f0(nb_element_x, nb_element_y))
allocate(V0(nb_element_x, nb_element_y))

! Loading
allocate(normal_loading_dot(nb_element_x, nb_element_y))
allocate(shearX_loading_dot(nb_element_x, nb_element_y))
allocate(shearY_loading_dot(nb_element_x, nb_element_y))


! Initial parameters
allocate(Vx(nb_element_x, nb_element_y))
allocate(Vy(nb_element_x, nb_element_y))

allocate(theta(nb_element_x, nb_element_y))
allocate(sigmaN0(nb_element_x, nb_element_y))

! Geometry
allocate(height(nb_element_x, nb_element_y))
allocate(mask(nb_element_x, nb_element_y))
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Read all the given nodes
!
subroutine read_ini()
integer(HID_T)  :: file, dset
integer :: hdferr
integer(HSIZE_T), DIMENSION(2) :: dims 
	
! Calculate dims
dims(1) = nb_element_x
dims(2) = nb_element_y
CALL h5open_f(hdferr) ! Initialise fortran interface
CALL h5fopen_f( trim(path_problem)//"parameters.h5", H5F_ACC_RDWR_F , file, hdferr) ! Open the file 

!



! Read all the values of friction
CALL h5dopen_f (file, "friction/a", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,a,dims, hdferr)  ! Read the data using the default properties.
CALL h5dopen_f (file, "friction/b", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,b,dims, hdferr)
CALL h5dopen_f (file, "friction/Dc", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Dc,dims, hdferr)
CALL h5dopen_f (file, "friction/f0", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,f0,dims, hdferr)
CALL h5dopen_f (file, "friction/V0", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,V0,dims, hdferr)

!
! Read all the values of loading 
CALL h5dopen_f (file, "loading/normal_loading_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,normal_loading_dot,dims, hdferr)
CALL h5dopen_f (file, "loading/shearX_loading_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,shearX_loading_dot,dims, hdferr)
CALL h5dopen_f (file, "loading/shearY_loading_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,shearY_loading_dot,dims, hdferr)
! 
! ! Read all the values of initial_parameters
CALL h5dopen_f (file, "/initial_parameters/Vx", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Vx,dims, hdferr)
CALL h5dopen_f (file, "/initial_parameters/Vy", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Vy,dims, hdferr)
CALL h5dopen_f (file, "/initial_parameters/theta", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,theta,dims, hdferr)
CALL h5dopen_f (file, "/initial_parameters/sigmaN", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,sigmaN0,dims, hdferr)
! 
! ! Real all the values of geometry
CALL h5dopen_f (file, "/geometry/height", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,height,dims, hdferr)
CALL h5dopen_f (file, "/geometry/mask", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,mask,dims, hdferr)

! Close and release resources.
CALL h5dclose_f(dset , hdferr)
CALL h5fclose_f(file , hdferr)

end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       Print every parameter at screen
!
subroutine show_variables()
use omp_lib
integer :: nb_thread_omp

!

print*,' //////////////////////////////////////////////////'
print*,' //                                              //'
print*,' //                  FASTCYCLE3D                 //'
print*,' //      Pierre Romanet, So Ozawa, Harsha Bhat,  //'
print*,' //      Tatsuhiko Saito, and Eiichi Fukuyama    //'
print*,' //                                              //'
print*,' //////////////////////////////////////////////////'
print*,''
nb_thread_omp = omp_get_max_threads()
print*,'Thread OMP = ', nb_thread_omp
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
print*,"%%%%%%%%%%%%%%% Hyperparameter %%%%%%%%%%%%%%%%%"
print*
print*,"eta_w = ", eta_w
print*,"rho_c = ", rho_c
print*,"beta_min = ", beta_min
print*,"gamma = ", gamma
print*,"nr_accuracy = ", nr_accuracy
print*,"tol_solver = ", tol_solver
print*,"quasi_dynamic = ", quasi_dynamic

print*
print*,"%%%%%%%%%%%%%%% Time to compute %%%%%%%%%%%%%%%%"
print*
print*
print*
print*,'Nb of elements = ', nb_element_x, nb_element_y
print*
print*
print*,'Please press ENTER'





end subroutine

end module
