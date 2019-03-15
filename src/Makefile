FC=gfortran
CFLAGS = -Wall -O2 -fcheck=bounds -Wunused-parameter #-fbounds-check

#***name of program and objects
PROGRAM = natsure
OBJS    = variables.o \
  		  parameters.o \
				init.o \
				analysis.o \
		  io.o \
			PBC.o\
          lattice.o \
		  potentials.o \
			bookKeeping.o \
		  thermostat.o \
  	      integrators.o \
          md.o

#.SUFFIXES: %.o %.f90
#*.o: *.f03
#	$(FC) -c $<

#***Making program
$(PROGRAM): $(OBJS)
	$(FC) -o $@ $^ $(CFLAGS)

variables.o: variables.f08
	$(FC) -c $<
parameters.o: parameters.f08
	$(FC) -c $<
init.o: init.f08
	$(FC) -c $<
analysis.o: analysis.f08
	$(FC) -c $<
PBC.o: PBC.f08
	$(FC) -c $<
io.o: io.f08
	$(FC) -c $<
lattice.o: lattice.f08
	$(FC) -c $<
potentials.o: potentials.f08
	$(FC) -c $<
bookKeeping.o: bookKeeping.f08
	$(FC) -c $<
thermostat.o: thermostat.f08
	$(FC) -c $<
integrators.o: integrators.f08
	$(FC) -c $<
md.o: md.f08
	$(FC) -c $<
#
.PHONY:clean
clean:
	rm $(OBJS) $(PROGRAM) *.mod
