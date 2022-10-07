# GNU compilation

INCPATH=-I/Users/pierre/Dropbox/OnGoingResearch/FastCycles_openSource/library/fftw-3.3.10/api
LIBS=-L/usr/local/lib -lfftw3 -Llibhdf5_fortran.a

F90=/Users/pierre/Dropbox/OnGoingResearch/FastCycles_openSource/library/HDF5/hdf5-1.12.2/hdf5/bin/h5fc

F90FLAGS=$(INCPATH) -ffree-line-length-none -O3 -w -fopenmp -c
CFLAGS=$(INCPATH)


SRC=src
DST=build


OBJ = variables.o default_values.o input_output.o special_functions.o initialisation.o solver.o

software: $(OBJ) master.f90
	$(F90) -ffree-line-length-none -O3 -w -fopenmp $(INCPATH) $(LIBS) $(OBJ) master.f90 -o ../$(DST)/fastcycles 

input_output.o: variables.o default_values.o input_output.f90
	$(F90) $(F90FLAGS) $(FLAG) input_output.f90

default_values.o: variables.o default_values.f90
	$(F90) $(F90FLAGS) default_values.f90

solver.o: variables.o solver.f90
	$(F90) $(F90FLAGS) solver.f90
	
friction.o: friction.f90
	$(F90) $(F90FLAGS) friction.f90

initialisation.o: variables.o special_functions.o initialisation.f90
	$(F90) $(F90FLAGS) initialisation.f90 

%.o: %.f90
	$(F90) $(F90FLAGS)   -c $< -o  $@

%.o: %.f
	$(F90) $(F90FLAGS)   -c $< -o  $@

variables.o: variables.f90
	$(F90) $(F90FLAGS) variables.f90

clean:
	rm *.o *.mod ../$(DST)/fastcycles
