FFLAGS = -lblas -llapack -O3 -fopenmp -finit-real=zero -g -fcheck=all -Wall -fbacktrace -fcheck=bounds
FC = gfortran
F77FLAGS = $(FFLAGS)
HPROG1 = iypoint32

 
DEP1  = mathmod.o module.o  crystalplasticity.o iypoint32.f95
DEP2  = mathmod.o module.f95
DEP3  = mathmod.o module.o  crystalplasticity.f95

 
 
PROGS = $(HPROG1) 
 
all          : $(PROGS)
 
$(HPROG1)    : $(DEP1)
	$(FC) $(DEP1) -o $@ $(FFLAGS)
 

 
module.o     : $(DEP2)
	$(FC) $(FFLAGS) -c module.f95
mathmod.o     : mathmod.f95
	$(FC) $(FFLAGS) -c mathmod.f95
crystalplasticity.o : $(DEP3)
	$(FC) $(FFLAGS) -c crystalplasticity.f95 
clean:
	rm module.o mathmod.o crystalplasticity.o
	rm Dp_*


#iypoint32: iypoint32.f95 module.f95
#	gfortran mathmod.f95 module.f95 iypoint32.f95 -o iypoint32 -lblas -llapack -O3 -fopenmp -fcheck=bounds	-finit-real=zero
#module
 