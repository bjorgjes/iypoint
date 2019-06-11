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
  
 Module.f95 contents and description
  
 line 21 : subroutine init()
           subroutine assigining values to all global parameters:
            - Elasticity
            - timestep
            - Hardening
            - Consitutive model parameters
            - Slip systems
            - Reads euler angles from file: "euleranglesin.txt"
            
  line 106: function voigt(A)
            Multiplies a strain tensor with the elasticity tensor using Voigt notation.
            
  line 141: subroutine eulang(R,phi1,Phi,phi2)
            Computes Bunge Euler angles from rotation matrix R
 
  line 159 : subroutine Acoeff(alpha,beta, Ctr, coeff)
             Calculate the coefficient A^(alpha beta) used for building the system of equations in the crystal plasticity model
             
  line 179: subroutine eqvstress(tag,sigma)
            Calculate equivalent stress of a Hoshford/Hersey yield surface from the stress tensor tag. 
            Following the procedure derived in Barlat 1991
   
  line 217: function kronecker(i,j)
            calculates the Kronecker delta
   
  line 229: subroutine rotmat(phi1,Phi, phi2,R)
            Calculates rotation matrix R from bunge euler angles
            
  line 247: subroutine slipsys(Schmid,m,n,snum)
            Calculate the Schmid matrix, m and n given a slip system number
            
  line 270: function gaveps(epsp)
            Calculates g(epsp) used in the constitutive model, where the fitting paramaeters input.
            
  line 286: function haveps(epsp)
            Calculates h(epsp) used in the constitutive model, where the fitting paramaeters input.
            
  line 303 : function alphacoeff(theta, epsp, modelnum, G)
             Calculates the parameters alpha used in the constitutive model
             
  line 332 : subroutine hoshfordnormal(Tag,grad)
  
  
  mathmod.f95 : Contents and description
  
  
  line 9 :   function contract2(T,S)
            Calculates double dot product of two second order tensors
            
  line 24:  subroutine deter(A,det) 
            calculates the determinant of a 3x3 matrix 
            
  line 42 : subroutine POLAR(F,R)
            Calculates polar decomposition 
  
  line 118: subroutine minv3(A,B)
            Calculates inverse of a 3x3 matrix
            
  line 142: subroutine matinv4(A,B)
            Calculates inverse of a 4x4 matrix
            
  line 179: function tens2vec(T)
            Converts a second order symmetric tensor to a six dimensional vector.
        
  line 192: function vec2tens(v)
            Converts a six dimensional vector to a second order symmetric tensor
            
  line 210: function macauley(h)
            Macauley brackets operator
            
            
            
  iypoint32.f95: Contents and description 
  
  line 1 : main program
  
  line 240: subroutine steepestdecent(solution, initial)
  
  
  
  
  constitutive_model.f95 Contents and description
  
  line 12: subroutine constexpr(l,part,bryter,bcond,strain,tag,epsp,propconst,fid)
            The subroutine used for the constitutive model, boundary condition and initial stress and strain is input. 
  
  line 404: subroutine elasticsolution(D,tag)
            Subroutine to perform elastic step using crystal plasticity model and calculate elastic parameters.
     
  line 470: subroutine yoshi2(tag,D,epsp,dt0,consistent,Dp) 
            Old version of subroutine which calculates a timestep from the constitutive model.
  
  
  
  
  
            
  
 
            
