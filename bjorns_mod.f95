module Bjorns_numgrad

contains
      subroutine Hersey(sigeq,sig,n)
!     returns equivalent stress sigeq
!     sigeq: input, equivalent stress, unchanged
!     sig(3,3): deviatoric stress tensor, not changed.
!     s(3): major stresses, calculated
!     n: exponent of the steady yield surface
      
      implicit none
      real(8) :: n
      double precision sig(3,3), sig_Y
      double precision s(3)
      double precision sigeq

      call majors(sig,s)

!     Hersey
      
      sigeq = abs(s(1)-s(2))**n &
          +abs(s(2)-s(3))**n &
          +abs(s(3)-s(1))**n
      sigeq =(0.5*sigeq)**(1./n)

      
      return
      end
!
!--------------------------------------
!      
      subroutine grad_phi(grad,sig,n)
!      
!     returns grad(phi), where phi is the steady yield function.
!     sig(3,3): stress tensor, unchanged
!     n: exponent of the steady yield surface
!     
      implicit none
      integer i, j
      real(8) :: n
      real(8) , dimension(3,3) :: grad, sig, sig1, sig2
      double precision sigeq,sigeq1,sigeq2,dsig, tmp, tmp1

      tmp=0.
      do i=1,3
         do j=1,3
            tmp1= abs(sig(i,j))
            tmp=max(tmp,tmp1)
         enddo
      enddo
      dsig = tmp*1.E-6
      
      do i = 1,3
         do j = 1,3
            sig1(i,j)=sig(i,j)
            sig2(i,j)=sig(i,j)
         enddo
      enddo

      call Hersey(sigeq,sig,n)

      do i = 1,3
         do j = 1,i
            sig1(i,j) = sig1(i,j) - dsig
            call Hersey(sigeq1,sig1,n)
            sig1(i,j) = sig1(i,j) + dsig
            sig2(i,j) = sig2(i,j) + dsig
            call Hersey(sigeq2,sig2,n)
            sig2(i,j) = sig2(i,j) - dsig
            grad(i,j) = (sigeq2 - sigeq1)/(2.*dsig)
            if(i.ne.j) grad(j,i) = grad(i,j)
         enddo
      enddo

      return
      end
!      
!
!--------------------------------------
!
      subroutine majors(sig,s)
!     sig is the stress tensor (destroyed on output)
!     s is the eigenstresses in descending order
      implicit none
      integer i,j,nrot
      double precision sig(3,3), v(3,3), s(3), tmp,store(3,3)
      do i = 1,3
         do j=1,3
            store(i,j)=sig(i,j)
         enddo
      enddo

      call b_jacobi(sig,3,3,s,v,nrot)

      do i = 1,3
         do j=1,3
            sig(i,j)=store(i,j)
         enddo
      enddo
      
      do i =2,3
         if(s(i).gt.s(1))then
            tmp = s(1)
            s(1)=s(i)
            s(i)=tmp
         endif
      enddo

      if(s(3).gt.s(2))then
         tmp = s(2)
         s(2)=s(3)
         s(3)=tmp
      endif

      RETURN
      end
!
!--------------------------------------
!
      subroutine dotprod(A,B,dot)
      implicit none
      integer i,j
      double precision A(3,3),B(3,3),dot
      dot=0.
      do i =1,3
         do j = 1,3
            dot=dot+A(i,j)*B(i,j)
         enddo
      enddo
      RETURN
      end
!     
!
!
      subroutine subtract(A,B,C)
      implicit none
      integer i,j
      double precision A(3,3),B(3,3),C(3,3)
      do i =1,3
         do j = 1,3
            C(i,j)=A(i,j)-B(i,j)
         enddo
      enddo
      RETURN
      end
!
!--------------------------------------
!
      SUBROUTINE B_JACOBI(A,N,NP,D,V,NROT)

!  Purpose: Computes all eigenvalues and eigenvectors of a real
!     symmetric matrix A, which is of size N by N, stored in a
!     physical NP by NP array.  On output, elements of A above the
!     diagonal are destroyed.  D returns the eigenvalues of A in
!     its first N elements.  V is a matrix with the same logical and
!     physical dimensions as A whose columns contain, on output, the
!     normalized eigenvectors of A.  NROT returns the number of Jacobi
!     rotations which were required. 

!  Source: W. H. Press et al., "Numerical Recipes", 1989, p. 346.

!  Modifications:

!     1. Double precision version

!  Prepared by J. Applequist, 10/23/91

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

!     Initialize the identity matrix.

      DO 12 IP=1,N
      DO 11 IQ=1,N
      V(IP,IQ)=0.D0
 11   CONTINUE
      V(IP,IP)=1.D0
 12   CONTINUE

!     Initialize B and D to the diagonal of A.

      DO 13 IP=1,N
      B(IP)=A(IP,IP)
      D(IP)=B(IP)
      Z(IP)=0.D0
 13   CONTINUE
      NROT=0
      DO 24 I=1,50
      SM=0.D0

!     Sum off-diagonal elements.

      DO 15 IP=1,N-1
      DO 14 IQ=IP+1,N
      SM=SM+DABS(A(IP,IQ))
 14   CONTINUE
 15   CONTINUE
      IF (SM.EQ.0.D0) RETURN
      IF (I.LT.4) THEN
      TRESH=0.2D0*SM/N**2
      ELSE
      TRESH=0.D0
      ENDIF
      DO 22 IP=1,N-1
      DO 21 IQ=IP+1,N
      G=100.D0*DABS(A(IP,IQ))

!     After four sweeps, skip the rotation if the off-diagonal
!     element is small.

      IF ((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP))) &
       .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ)))) THEN
      A(IP,IQ)=0.D0
      ELSE IF (DABS(A(IP,IQ)).GT.TRESH) THEN
        H=D(IQ)-D(IP)
        IF (DABS(H)+G.EQ.DABS(H)) THEN
        T=A(IP,IQ)/H
        ELSE
        THETA=0.5D0*H/A(IP,IQ)
        T=1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
        IF (THETA.LT.0.D0) T=-T
        ENDIF
      C=1.D0/DSQRT(1.D0+T**2)
      S=T*C
      TAU=S/(1.D0+C)
      H=T*A(IP,IQ)
      Z(IP)=Z(IP)-H
      Z(IQ)=Z(IQ)+H
      D(IP)=D(IP)-H
      D(IQ)=D(IQ)+H
      A(IP,IQ)=0.D0
      DO 16 J=1,IP-1
      G=A(J,IP)
      H=A(J,IQ)
      A(J,IP)=G-S*(H+G*TAU)
      A(J,IQ)=H+S*(G-H*TAU)
 16   CONTINUE
      DO 17 J=IP+1,IQ-1
      G=A(IP,J)
      H=A(J,IQ)
      A(IP,J)=G-S*(H+G*TAU)
      A(J,IQ)=H+S*(G-H*TAU)
 17   CONTINUE
      DO 18 J=IQ+1,N
      G=A(IP,J)
      H=A(IQ,J)
      A(IP,J)=G-S*(H+G*TAU)
      A(IQ,J)=H+S*(G-H*TAU)
 18   CONTINUE
      DO 19 J=1,N
      G=V(J,IP)
      H=V(J,IQ)
      V(J,IP)=G-S*(H+G*TAU)
      V(J,IQ)=H+S*(G-H*TAU)
 19   CONTINUE
      NROT=NROT+1
      ENDIF
 21   CONTINUE
 22   CONTINUE
      DO 23 IP=1,N
      B(IP)=B(IP)+Z(IP)
      D(IP)=B(IP)
      Z(IP)=0.D0
 23   CONTINUE
 24   CONTINUE
      WRITE (6,600) 
 600  FORMAT(/'50 ITERATIONS OCCURRED IN SUBROUTINE JACOBI.')
      RETURN
      END
      
   
end module 