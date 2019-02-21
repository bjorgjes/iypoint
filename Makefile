iypoint32: iypoint32.f95 module.f95
	gfortran module.f95 iypoint32.f95 -o iypoint32 -lblas -llapack -O3 -fopenmp -fcheck=bounds	-finit-real=zero
