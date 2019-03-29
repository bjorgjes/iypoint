module mathmod


contains    

!! PERFORMS DOBLE DOT PRODUCT ON TWO 3X3 TENSORS 
!!
!!
    function contract2(T,S)
        implicit none
    
        real(8), dimension(3,3):: T, S
        real(8) ::  contract2
        contract2 = 0
        contract2 = T(1,1)*S(1,1) + T(2,2)*S(2,2) + T(3,3)*S(3,3) + &
                T(1,2)*S(1,2) + T(2,1)*S(2,1) + T(1,3)*S(1,3) + &
                T(3,1)*S(3,1) + T(2,3)*S(2,3) + T(3,2)*S(3,2)
        
    return
    end function contract2


!! Subroutine which calculate the determinant of a 3x3 matrix
    subroutine deter(A,det) 
        !Declear variables
            real(8), dimension(3,3),intent(in) :: A
            real(8), intent(out) :: det
        
            det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        return
end subroutine deter



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


!! CONVERTS A SYMETRIC TENSOR INTO A VECTOR
!!
function tens2vec(T)
    implicit none 
    real(8),  dimension(3,3) :: T
    real(8),  dimension(6) :: tens2vec

    tens2vec = (/ T(1,1),T(2,2), T(3,3),T(2,3), T(1,3),&
                T(1,2) /)

return
end function tens2vec

!! CONVERTS A VECTOR INTO A SYMMETRIC TENSOR
!!
function vec2tens(v)
    implicit none 
    real(8),  dimension(3,3) :: vec2tens
    real(8), dimension(6) :: v

vec2tens(1,1) = v(1)
vec2tens(2,2) = v(2)
vec2tens(3,3) = v(3)
vec2tens(2,3) = v(4)
vec2tens(3,2) = v(4)
vec2tens(1,3) = v(5)
vec2tens(3,1) = v(5)
vec2tens(1,2) = v(6)
vec2tens(2,1) = v(6)

return
end function vec2tens

end module mathmod