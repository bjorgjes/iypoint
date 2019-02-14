module global
    implicit none

    real(8) , dimension(6,12)               :: slip
    real(8) , dimension(6,6)                :: Cel
    real(8) , dimension(:,:), Allocatable   :: eul
    integer                                 :: nlines
    real(8), dimension(3,3)                 :: id
    real(8)                                 ::  dt, pi 
    real(8) , dimension(:,:,:), Allocatable :: R    
    logical                                 :: hardening                              
    

    contains
   
   !!!!! SUBROUTINE FOR INITIALIZATION OF MATERIAL AND MODEL PARAMETERS
   
    subroutine init()
    integer :: i
    real(8) :: phi1, Phi, phi2
        pi = 4.D0*DATAN(1.D0)
        dt = 0.0001
        hardening = .true.
        
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


!! Subroutine which inverts 3x3 matrix directly
subroutine minv3(A,B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        real(8), dimension(3,3)  :: A   !! Matrix
        real(8), dimension(3,3)             :: B  !! Inverse matrix
        real(8)            :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
        return
end subroutine minv3

!! Subroutine which calculate the determinant of a 3x3 matrix
subroutine deter(A,det) 
        !Declear variables
            real(8), dimension(3,3),intent(in) :: A
            real(8), intent(out) :: det
        
            det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        return
end subroutine deter


!Subrotuine which multiply tensor A with the 4th order elasticity tensor, following Voigt notation,  and returns the resulting tensor
subroutine voigt(A,B)
        real(8) , dimension(6) :: vec
        real(8) , dimension(3,3), intent(in) :: A
        real(8) , dimension(3,3), intent(out) :: B
        
        !vec = (/A(1,1),A(2,2),A(3,3),A(2,3),A(1,3),A(1,2)/)
        !vec(4:6) = 2*vec(4:6)
        vec = matmul(Cel,(/A(1,1), A(2,2),A(3,3), 2*A(2,3),2*A(1,3),2*A(1,2)/))
        
        !B(1,1) = DOT_PRODUCT(Cel(1,1:6),vec)
        !B(2,2) = DOT_PRODUCT(Cel(2,1:6),vec)
        !B(3,3) = DOT_PRODUCT(Cel(3,1:6),vec)
        !B(1,2) = DOT_PRODUCT(Cel(6,1:6),vec)
        !B(2,1) = DOT_PRODUCT(Cel(6,1:6),vec)
        !B(1,3) = DOT_PRODUCT(Cel(5,1:6),vec)
        !B(3,1) = DOT_PRODUCT(Cel(5,1:6),vec)
        !B(3,2) = DOT_PRODUCT(Cel(4,1:6),vec)
        !B(2,3) = DOT_PRODUCT(Cel(4,1:6),vec)
    
        B(1,1) = vec(1)
        B(2,2) = vec(2)
        B(3,3) = vec(3)
        B(1,2) = vec(6)
        B(2,1) = vec(6)
        B(1,3) = vec(5)
        B(3,1) = vec(5)
        B(3,2) = vec(4)
        B(2,3) = vec(4)
end subroutine voigt


!C**********************************************************************
!C                         FUNCTION POLAR                              *
!C**********************************************************************
!C Polar decomposition of F = R.U using the Cayley-Hamilton theorem    *
!C see Nemat-Nasser's book p.55                                        *
!C Returns R                                                           *
!C**********************************************************************
subroutine POLAR(F,R)
    !C
    IMPLICIT NONE
    !C
    REAL(kind=8)    :: F(3,3),R(3,3),C(3,3),CS(3,3),U(3,3),UI(3,3),C1,C3,P, &
                       CD11,CD22,CD33,CD2,CD3,U1,U2,U3,A,B,PHI,L1,D,E,A3,B2
    INTEGER         :: I,J,K
   
    intent(in)      :: F

    !C Compute stretch tensor: C = F^T.F
    DO J=1,3
        DO I=1,3
            C(I,J)=0.D0
            DO K=1,3
                C(I,J)=C(I,J)+F(K,I)*F(K,J)
            END DO
        END DO
    END DO
    !C Compute C^2
    CS = MatMUL(C, C)
    !C Compute invariants
    C1=C(1,1)+C(2,2)+C(3,3)

    C3=C(1,1)*(C(2,2)*C(3,3)-C(2,3)*C(3,2))+ &
       C(1,2)*(C(2,3)*C(3,1)-C(2,1)*C(3,3))+ &
       C(1,3)*(C(2,1)*C(3,2)-C(2,2)*C(3,1))
   
    !C Invariants of the deviatoric part CD of tensor C
    P=(C(1,1)+C(2,2)+C(3,3))/3.D0
    
    CD11=C(1,1)-P
    CD22=C(2,2)-P
    CD33=C(3,3)-P
   
    CD2=CD11*CD22+CD11*CD33+CD22*CD33- &
        (C(1,2)*C(2,1)+C(1,3)*C(3,1)+C(2,3)*C(3,2))
    CD3=CD11*(CD22*CD33-C(2,3)*C(3,2))+ &
        C(1,2)*(C(2,3)*C(3,1)-C(2,1)*CD33)+ &
        C(1,3)*(C(2,1)*C(3,2)-CD22*C(3,1))
    !C Invariants of U
    U3=sqrt(C3)
    A=-CD2/3.D0
    B=CD3/2.D0
    A3=A**3.D0
    B2=B*B
   
    IF (ABS(A3-B2).GT.1.D-12) THEN
        PHI=ACOS(B/A**(3.D0/2.D0))
        L1=SQRT(C1/3.D0+2.D0*SQRT(A)*COS(PHI/3.D0))
    ELSE
        L1=SQRT(C1/3.D0)
    END IF
   
    U1=L1+SQRT(C1-L1*L1+2.D0*U3/L1)
    U2=0.5D0*(U1*U1-C1)
    !C Computes U
    D=U3-U1*U2
    E=U1**2.D0-U2
    DO I=1,3
        DO J=1,3
            U(I,J)=(CS(I,J)-E*C(I,J))/D
            IF (I.EQ.J) THEN
                U(I,J)=U(I,J)-U1*U3/D
            END IF
        END DO
    END DO
    call minv3(U,UI)
    R = matmul(F,UI)
return
end subroutine POLAR


!!! Calculate the von Mises equivalent stress givan a cauchy stress tensor T
subroutine eqvstr(T,sigmaeq)
    implicit none
    real(8) , intent(in), dimension(3,3) :: T
    real(8) , intent(out) :: sigmaeq


sigmaeq = 1.00/2.00*((T(1,1)-T(2,2))**2+(T(2,2)-T(3,3))**2+(T(3,3)-T(1,1))**2+6.00*(T(1,2)**2+T(2,3)**2+T(3,1)**2))

sigmaeq = sqrt(sigmaeq)
return
end subroutine eqvstr


!! Directly calculates the inverse of a 4x4 matrix
subroutine matinv4(A,B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(8), intent(in) :: A(4,4)   !! Matrix
    real(8)             :: B(4,4)   !! Inverse matrix
    real(8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    return 
end subroutine


!! PERFORMS DOBLE DOT PRODUCT ON TWO 3X3 TENSORS 
!!
!!
subroutine contract2(T,S,dprod)
    implicit none

    real(8), dimension(3,3), Intent(in) :: T, S
    real(8), Intent(out) :: dprod
    dprod = 0
    dprod = T(1,1)*S(1,1) + T(2,2)*S(2,2) + T(3,3)*S(3,3) + &
            T(1,2)*S(1,2) + T(2,1)*S(2,1) + T(1,3)*S(1,3) + &
            T(3,1)*S(3,1) + T(2,3)*S(2,3) + T(3,2)*S(3,2)
    
return
end subroutine contract2

!! CONVERTS A SYMETRIC TENSOR INTO A VECTOR
!!
subroutine tens2vec(T,v)
    implicit none 
    real(8), intent(in), dimension(3,3) :: T
    real(8), intent(out), dimension(6) :: v

    v(1:6) = (/ T(1,1),T(2,2), T(3,3),T(1,2), T(1,3),&
                T(2,3) /)

return
end subroutine tens2vec

!! CONVERTS A VECTOR INTO A SYMMETRIC TENSOR
!!
subroutine vec2tens(T,v)
    implicit none 
    real(8), intent(out), dimension(3,3) :: T
    real(8), intent(in), dimension(6) :: v

T(1,1) = v(1)
T(2,2) = v(2)
T(3,3) = v(3)
T(2,3) = v(4)
T(3,2) = v(4)
T(1,3) = v(5)
T(3,1) = v(5)
T(1,2) = v(6)
T(2,1) = v(6)

return
end subroutine vec2tens

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


    subroutine Acoeff(alpha,beta, Ctr, coeff)

        !Declear variables
        integer, intent(in) :: alpha, beta
        real(8), dimension(3,3),intent(in) :: Ctr
        real(8), dimension(3,3) :: Salpha, Sbeta, rm, csv, tot
        real(8), dimension(3) :: m,n
        real(8), intent(out) :: coeff
    !Computation
        
        Call slipsys(slip,Salpha,m,n,alpha)
        Call slipsys(slip,Sbeta,m,n,beta)
        csv = matmul(Ctr,Sbeta)
        csv = (csv+transpose(csv))/2
        call voigt(csv,rm)
        !CS = (/csv(1,1),csv(2,2),csv(3,3),csv(2,3),csv(1,3),csv(1,2)/)
        !rm(1,1) = DOT_PRODUCT(Cel(1,1:6),CS)
        !rm(2,2) = DOT_PRODUCT(Cel(2,1:6),CS)
        !rm(3,3) = DOT_PRODUCT(Cel(3,1:6),CS)
        !rm(1,2) = DOT_PRODUCT(Cel(6,1:6),CS)
        !rm(2,1) = DOT_PRODUCT(Cel(6,1:6),CS)
        !rm(1,3) = DOT_PRODUCT(Cel(5,1:6),CS)
        !rm(3,1) = DOT_PRODUCT(Cel(5,1:6),CS)
        !rm(3,2) = DOT_PRODUCT(Cel(4,1:6),CS)
        !rm(2,3) = DOT_PRODUCT(Cel(4,1:6),CS)
    tot = matmul(transpose(Salpha),rm)
    coeff = tot(1,1)+tot(2,2)+tot(3,3)
    return
    end subroutine Acoeff

end module