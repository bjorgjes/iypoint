FFLAGS = -lblas -llapack -O3 -fopenmp -finit-real=zero -g -fcheck=all -Wall -fbacktrace -fcheck=bounds
FC = gfortran
F77FLAGS = $(FFLAGS)
HPROG1 = iypoint32

 
DEP1  = mathmod.o module.o bjorns_mod.o crystalplasticity.o iypoint32.f95
DEP2  = mathmod.o bjorns_mod.o module.f95 
DEP3  = mathmod.o module.o bjorns_mod.o crystalplasticity.f95
DEP4  = bjorns_mod.f95
 
 
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
bjorns_mod.o : $(DEP4)
	$(FC) $(FFLAGS) -c bjorns_mod.f95

clean:
	rm module.o mathmod.o crystalplasticity.o
	rm Dp_*

git:
	git add module.f95
	git add crystalplasticity.f95
	git add mathmod.f95
	git add iypoint32.f95
	git add Makefile
	@read -p "Enter commit message:" messag; \
	git commit -m "$$messag"
	git push



#iypoint32: iypoint32.f95 module.f95
#	gfortran mathmod.f95 module.f95 iypoint32.f95 -o iypoint32 -lblas -llapack -O3 -fopenmp -fcheck=bounds	-finit-real=zero
#module
 