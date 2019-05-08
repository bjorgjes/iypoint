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
    
    !!$OMP THREADPRIVATE(slip,Cel,eul,nlines,id,R, hardening, pi, dt)
    contains
   
   !!!!! SUBROUTINE FOR INITIALIZATION OF MATERIAL AND MODEL PARAMETERS
   
    subroutine init()
    integer :: i
    real(8) :: phi1, Phi, phi2
        pi = 4.D0*DATAN(1.D0)
        dt = 0.0001
        dgamma = 0.00001
        hardening = .false.
        
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
            write(*,*) I3/I2**(3./2.)
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


end module