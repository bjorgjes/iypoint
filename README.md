# iypoint
Framework for computation of instantaneous yield surfaces 
README
The code is written in fortran 95 , and compiled using gfortran.
The compilation parameters are provided in the makefile, however some packages are needed:
--
  -- blas , linear algebra package
  -- lapack , Another more sophisticated linear algebra package
  -- openmp , multithread package. Compilation without this package is ok, 
      but there is considerable timesavings to be made by running mulitple cores. 
      
  Running Ubuntu, these packages, including the gfortran compiler may be installed using apt-get. 
  The code have been written using Visual Studio Code(VSC), which have worked well. 
  
