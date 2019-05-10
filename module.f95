module global
    use mathmod
    use Bjorns_numgrad
    implicit none

    real(8) , dimension(6,12)               :: slip
    real(8) , dimension(6,6)                :: Cel
    real(8) , dimension(:,:), Allocatable   :: eul
    integer                                 :: nlines
    real(8), dimension(3,3)                 :: id
    real(8)                                 ::  dt, pi , dgamma
    real(8) , dimension(:,:,:), Allocatable :: R    
    logical                                 :: hardening  
    real(8)                                 :: c1, c2 , c3 ,theta0                            
    
    !!$OMP THREADPRIVATE(slip,Cel,eul,nlines,id,R, hardening, pi, dt)
    contains
   
   !!!!! SUBROUTINE FOR INITIALIZATION OF MATERIAL AND MODEL PARAMETERS
   
    subroutine init()
    integer :: i
    real(8) :: phi1, Phi, phi2
        pi = 4.D0*DATAN(1.D0)
        dt = 0.001
        dgamma = 0.00001
        hardening = .false.
     !!!! parameters for constitutive model
        c1 = 0.3 
        c2 = 0.5
        c3 = 0.4
        theta0 = pi/18

    !Create identity matrix
    id = 0
    forall(i = 1:3) id(i,i) = 1
    
    !Elasticity tensor
    Cel = 0
    forall (i = 1:3) Cel(i,i)= 170*1000
    forall (i = 4:6) Cel(i,i)= 75*1000
    forall (i = 1:2) Cel(i+1,i) = 124*1000
    forall (i = 1:2)  Cel(i,i+1) = 124*1000
    Cel(3,1) = 124*1000
    Cel(1,3) = 124*1000
    
    ! The crystal slip systems
    slip(1,1:12) = (/ 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0),-1/sqrt(3.0),-1/sqrt(3.0),-1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0),&
     1/sqrt(3.0),-1/sqrt(3.0), &
    -1/sqrt(3.0),-1/sqrt(3.0)/)
    slip(2,1:12) = (/ 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0),-1/sqrt(3.0),-1/sqrt(3.0),&
    -1/sqrt(3.0),-1/sqrt(3.0), &
    -1/sqrt(3.0),-1/sqrt(3.0)/)
    slip(3,1:12) = (/ 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0),&
     1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0), 1/sqrt(3.0)/)
    slip(4,1:12) = (/ 1/sqrt(2.0),-1/sqrt(2.0),0.0          , 1/sqrt(2.0),-1/sqrt(2.0), 0.0          , 0.0          ,-1/sqrt(2.0),&
     1/sqrt(2.0),-1/sqrt(2.0), 1/sqrt(2.0), 0.0          /)
    slip(5,1:12) = (/-1/sqrt(2.0), 0.0          , 1/sqrt(2.0),0.0          ,-1/sqrt(2.0), 1/sqrt(2.0),-1/sqrt(2.0),0.0          ,&
     1/sqrt(2.0), 1/sqrt(2.0),0.0          ,-1/sqrt(2.0)/)
    slip(6,1:12) = (/0.0          ,1/sqrt(2.0),-1/sqrt(2.0), 1/sqrt(2.0),0.0          ,-1/sqrt(2.0),-1/sqrt(2.0), 1/sqrt(2.0),&
     0.0          , 0.0          , 1/sqrt(2.0),-1/sqrt(2.0)/)
    
     call countlin('euleranglesin.txt',nlines)

     open(action='read',unit=15,file="euleranglesin.txt",status='old')
     !open(action='read',unit=15,file="oysurf001of.txt",status='old')
     allocate(eul(nlines,3))
     allocate(R(3,3,nlines))   
     do i = 1,nlines
         read(15,*) eul(i,1:3)
            phi1 = eul(i,1)
            phi = eul(i,2)
            phi2 = eul(i,3)
    
    !!!!!! Variable assignments end
    ! Calculate rotation matrix
    call rotmat(phi1,Phi,phi2,R(1:3,1:3,i))
    
        ! write(*,*) eul(i,1:3)
     end do
    end subroutine init
    

    Subroutine countlin(filename,nlines)
    implicit none
    character(len=*)    :: filename
    integer             :: nlines 
    integer             :: io
  
    open(10,file=filename, iostat=io, status='old')
    if (io/=0) stop 'Cannot open file! '
  
    nlines = 0
    do
      read(10,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
    end do
    close(10)
    return
end subroutine countlin



!Subrotuine which multiply tensor A with the 4th order elasticity tensor, following Voigt notation,  and returns the resulting tensor
function voigt(A)
        real(8) , dimension(6) :: vec
        real(8) , dimension(3,3) :: A
        real(8) , dimension(3,3) :: voigt
        
        vec = matmul(Cel,(/A(1,1), A(2,2),A(3,3), 2*A(2,3),2*A(1,3),2*A(1,2)/))
        
    
        voigt(1,1) = vec(1)
        voigt(2,2) = vec(2)
        voigt(3,3) = vec(3)
        voigt(1,2) = vec(6)
        voigt(2,1) = vec(6)
        voigt(1,3) = vec(5)
        voigt(3,1) = vec(5)
        voigt(3,2) = vec(4)
        voigt(2,3) = vec(4)
end function voigt


!!! Calculate the von Mises equivalent stress given a cauchy stress tensor T
subroutine eqvstr(T,sigmaeq)
    implicit none
    real(8) , intent(in), dimension(3,3) :: T
    real(8) , intent(out) :: sigmaeq


sigmaeq = 1.00/2.00*((T(1,1)-T(2,2))**2+(T(2,2)-T(3,3))**2+(T(3,3)-T(1,1))**2+6.00*(T(1,2)**2+T(2,3)**2+T(3,1)**2))

sigmaeq = sqrt(sigmaeq)
return
end subroutine eqvstr


!!! SUBROUTINE TO FIND EULERANGLES FROM ROTATION MATRIX
subroutine eulang(R,phi1,Phi,phi2)
    ! Declear variables
        real(8), dimension(3,3), intent(in) :: R
        real(8), intent(out) :: phi1, Phi, phi2
    ! Set limit for changes 
        if (R(3,3) > 0.999999 .or. R(3,3) < -0.999999 ) then
            Phi = 0 
            phi1 = 0
            phi2 = phi1 + atan2(R(1,2),R(1,1))
        else
          Phi = ACOS(R(3,3))
          phi2 = ATAN2(R(1,3),R(2,3))
          phi1 = ATAN2(-R(3,1),R(3,2))  
        end if
    
    end subroutine eulang

!!! Calculate the coefficient A^(alpha beta) used for building the system of equations in the crystal plasticity model
    subroutine Acoeff(alpha,beta, Ctr, coeff)

        !Declear variables
        integer, intent(in) :: alpha, beta
        real(8), dimension(3,3),intent(in) :: Ctr
        real(8), dimension(3,3) :: Salpha, Sbeta, rm, csv, tot
        real(8), dimension(3) :: m,n
        real(8), intent(out) :: coeff
    !Computation
        
        Call slipsys(Salpha,m,n,alpha)
        Call slipsys(Sbeta,m,n,beta)
        csv = matmul(Ctr,Sbeta)
        csv = (csv+transpose(csv))/2
        !call voigt(csv,rm)
    tot = matmul(transpose(Salpha),voigt(csv))
    coeff = tot(1,1)+tot(2,2)+tot(3,3)
    return
    end subroutine Acoeff

    subroutine eqvstress(tag,sigma)
        !! Calculate equivalent stress of a Hoshford/Hersey yield surface from the stress tensor tag. Following the procedure derived in Barlat 1991
        
        implicit none 
        real(8), dimension(3,3), intent(in) :: Tag
        real(8) :: A, B, C, F, G, H, I2, I3, theta, m = 9.d+0 ,sum, sFi, sigma
       if (norm2(tag) == 0) then
        sigma = 0
       else 

        A = tag(2,2)-tag(3,3)
        B = tag(3,3)-tag(1,1)
        C = tag(1,1)-tag(2,2)
        F = tag(2,3)
        G = tag(1,3)
        H = tag(1,2)
        I2 = (F**2+G**2+H**2)/3. + ((A-C)**2+(C-B)**2+(B-A)**2)/54.
        I3 = ((C-B)*(A-C)*(B-A))/54 + F*G*H - ((C-B)*F**2+(A-C)*G**2+(B-A)*H**2)/6
        if (I3/I2**(3./2.) > 1 ) then
            !write(*,*) I3/I2**(3./2.)
            theta = 0
        else if ( I3/I2**(3./2.) < -1 ) then
            !write(*,*) I3/I2**(3./2.)
            theta = pi
        else
        theta = acos(I3/I2**(3./2.))
        end if
        !write(*,*) theta

        sum = (2*cos((2*theta+pi)/6.))**(m)+((2*cos((2*theta -  3*pi)/6.)))**(m)  &
            + (-2*cos((2*theta +  5*pi)/6.))**(m)
            sFi = (3*I2)**(m/2.)*sum    
       
        sigma = (sFi/2.)**(1./m)
       end if 
    return
    end subroutine eqvstress

    function kronecker(i,j)
        implicit none
        integer :: i,j
        real(8) :: kronecker
        if (i == j ) then 
            kronecker = 1
        else 
            kronecker = 0
        end if 
        return
    end function kronecker

    subroutine rotmat(phi1,Phi, phi2,R)
        implicit none 
    real(8) :: phi1, Phi, phi2
    
    real(8) , dimension(3,3) :: R
    R(1,1) = cos(phi1)*COS(phi2)-sin(phi1)*sin(phi2)*COS(Phi)
    R(1,2)= sin(phi1)*COS(phi2)+cos(phi1)*sin(phi2)*COS(Phi)
    R(1,3)= sin(phi2)*sin(Phi)
    R(2,1)= -cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*COS(Phi)
    R(2,2)=-sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(Phi)
    R(2,3)=cos(phi2)*sin(Phi)
    R(3,1)=sin(phi1)*sin(Phi)
    R(3,2)=-cos(phi1)*sin(Phi)
    R(3,3)=cos(Phi)
    return
    end subroutine rotmat


    subroutine slipsys(Schmid,m,n,snum)

        !Declear all returning variables
        real(8), dimension(3,3),intent(out) :: Schmid
        
        integer :: snum, h, p
        real(8), dimension(3), intent(out) :: m , n 
        !Declear all parameters used inside the subroutine
        
        
        
        !Calculations
        n = slip(1:3,snum)
        m = slip(4:6,snum)
        
            do h = 1,3
                do p = 1, 3
                    Schmid (h,p) = m(h)*n(p)
                end do
            end do
        return
        end subroutine slipsys
        
function gaveps(epsp)
    implicit none
    real(8) :: sigma0, eps0, n, gaveps, epsp
    sigma0 = 49.103558072070918    
    !n = 0.645430868
    !eps0 =  0.016142686
    n = 0.0027724527
    eps0 = 0.002438999583505
    if (hardening) then
    gaveps = sigma0*(1+epsp/eps0)**n
    else if (.not. hardening) then 
        gaveps = sigma0
    end if
    return
end function gaveps

function haveps(epsp)
    implicit none
    real(8) :: sigma0, eps0, n, haveps,epsp
    sigma0 = 49.103558072070918    
    !n = 0.645430868
    !eps0 =  0.016142686
    n = 0.0027724527
    eps0 = 0.002438999583505
    
    if (hardening) then
    haveps = n*sigma0*(1+epsp/eps0)**(n-1)/eps0
    else if (.not. hardening) then
        haveps = 0
    end if 
    return
end function haveps

function alphacoeff(theta, epsp, modelnum, G)
implicit none
real(8) :: c1= 0.3, c2 = 0.5, c3 = 0.4 ,theta0
real(8) :: alphacoeff, theta, epsp, modelnum, G

if (modelnum == 1) then
    theta0 = pi/18.
else if (modelnum == 2) then
    theta0 = pi/18./2.
end if 


if (theta >= 0 .and. theta <= theta0 ) then
        alphacoeff = 1-c1*sqrt(gaveps(epsp)/G)-c2*macauley(haveps(epsp))/G
else if (theta > theta0 .and. theta < pi/2) then
    if (modelnum == 1) then 
        alphacoeff = (1-c1*sqrt(gaveps(epsp)/G)-c2*macauley(haveps(epsp))/G)*((pi/2-theta)/(pi/2-theta0))
    else if (modelnum == 2) then
        alphacoeff = (1-c1*sqrt(gaveps(epsp)/G)-c2*macauley(haveps(epsp))/G)*Tan(c3*(theta-theta0)+theta0)/tan(theta)
    end if
    !alpha = (pi/2-theta)/(pi/2-theta0)
else if ( theta >= pi/2 .and. theta < pi ) then
    alphacoeff = 0
    write(*,*) 'Warning unloading'
end if   
return
end function alphacoeff

!!! SUBROUTINE TO CALCULATE THE GRADIENT OF A GIVEN HOSHFORD YIELD SURFACE. 
subroutine hoshfordnormal(Tag,grad)
    !Subroutine for calculation of the GRADIENT of the Hoshford/Hersey yield surface
    ! The method is adapted from Barlat et al (1991) "A six-component yield function for anisotropic materials"
    ! Since the hoshford/Hersey yield surface is given in principle stresses, the eigenvalue problem has to be analytically solved in order to calculate the gradient directly
    ! The partial derivatives of A,B,C,F,G,H (Bishop-Hill notation) with respect to the six stress component are writen in the matrix partials.
    !   
    !               dA/ds11 dB/ds11    .....  dH/ds11  
    !               dA/ds22 dB/ds22             .
    ! partials =    dA/ds33    .    .           .
    !               dA/ds23    .      .         .    
    !               dA/ds13    .        .        .            
    !               dA/ds12    .        ..... dH/ds12
    !
  
    implicit none 
    real(8), dimension(3,3), intent(in) :: Tag
    real(8), dimension(3,3), intent(out) :: grad
    real(8) :: A, B, C, F, G, H, dI2,dI3, dtheta, I2, I3, theta, m = 9.0 , dsum,sum, sFi, fraction
    real(8) , dimension(6,6) :: partials
    integer :: i 
    real(8) , dimension(6) :: n
    
    grad = 0
    n = 0 
   
    partials(1,1:6) = (/  0.0 , -1.0,  1.0,  0.0,  0.0,  0.0/)
    partials(2,1:6) = (/  1.0 ,  0.0, -1.0,  0.0,  0.0,  0.0/)
    partials(3,1:6) = (/ -1.0 ,  1.0,  0.0,  0.0,  0.0,  0.0/)
    partials(4,1:6) = (/  0.0 ,  0.0,  0.0,  1.0,  0.0,  0.0/)
    partials(5,1:6) = (/  0.0 ,  0.0,  0.0,  0.0,  1.0,  0.0/)
    partials(6,1:6) = (/  0.0 ,  0.0,  0.0,  0.0,  0.0,  1.0/)

    A = tag(2,2)-tag(3,3)
    B = tag(3,3)-tag(1,1)
    C = tag(1,1)-tag(2,2)
    F = tag(2,3)
    G = tag(1,3)
    H = tag(1,2)
    I2 = (F**2+G**2+H**2)/3. + ((A-C)**2+(C-B)**2+(B-A)**2)/54.
    I3 = ((C-B)*(A-C)*(B-A))/54 + F*G*H - ((C-B)*F**2+(A-C)*G**2+(B-A)*H**2)/6
    if (I3/I2**(3./2.) > 1 ) then
        fraction = 1
    else if (I3/I2**(3./2.) < -1) then
        fraction = -1
    else
    fraction = (I3/I2**(3./2.))
    end if
    
    theta = acos(fraction)
    sum = ((2*cos((2*theta+pi)/6.))**(m))+((2*cos((2*theta -  3*pi)/6.))**(m)) + ((-2*cos((2*theta +  5*pi)/6.))**(m))
    sFi = (3*I2)**(m/2.)*sum
    
    !write(*,*) theta, sum, sFi  
    do i = 1,6
        dI2 = 2./3.*(F*partials(i,4)+G*partials(i,5)+H*partials(i,6)) &
        +2./54.*((A-C)*(partials(i,1)-partials(i,3))) &
        +2./54.*((C-B)*(partials(i,3)-partials(i,2))) &
        +2./54.*((B-A)*(partials(i,2)-partials(i,1)))
        !write(*,*) -1/sqrt(1-fraction**2)

        dI3 = 1./54.*((partials(i,3)-partials(i,2))*(A-C)*(B-A) & 
        + ((partials(i,1)-partials(i,3))*(B-A)+(A-C)*(partials(i,2)-partials(i,1)))*(C-B)) &
        + partials(i,4)*G*H + F*partials(i,5)*H + F*H*partials(i,6) &
        - (2*F/6*partials(i,4)*(C-B)+F**2/6*(partials(i,3)-partials(i,2))) &
        - (2*G/6*partials(i,5)*(A-C)+G**2/6*(partials(i,1)-partials(i,3))) &
        - (2*H/6*partials(i,6)*(B-A)+H**2/6*(partials(i,2)-partials(i,1))) 
        
        if (fraction == 1 .or. fraction == -1) then
            dtheta = 0
        else
        dtheta = -1/sqrt(1-fraction**2)*(dI3*I2**(3./2.)-3./2.*I2**(1./2.)*dI2*I3)/I2**3
        end if
        !write(*,*) (dI3*I2**(3./2.)-3./2.*I2**(1./2.)*dI2*I3)/I2**3
       
        !write(*,*) dtheta
        
        dsum = m*dtheta* &
        (  ( 2*cos((2*theta +   pi)/6.))**(m-1) * (-4./6.*sin((2*theta +   pi)/6)) &
         + ( 2*cos((2*theta - 3*pi)/6.))**(m-1) * (-4./6.*sin((2*theta - 3*pi)/6.)) &
         + (-2*cos((2*theta + 5*pi)/6.))**(m-1) * ( 4./6.*sin((2*theta + 5*pi)/6)))
    
        n(i) = m/2*3**(m/2.)*I2**(m/2.-1.)*dI2*sum+(3*I2)**(m/2.)*dsum

        !write(*,*) dI2, dI3, dtheta, dsum
    end do
       !n = n/norm2(n)
!write(*,*) n/norm2(n)
!write(8,*) tag(1,1), tag(2,2), n(1), n(2)
!call vec2tens(grad,n)
grad = vec2tens(n)


grad = grad*1/(2*m)*(sFi/2)**(1/m-1)
!write(*,*) grad

return
end subroutine hoshfordnormal


end module