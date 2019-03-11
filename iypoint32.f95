program ypoint
use crystalplasticity
    implicit none

integer :: part,bryter, k, i,bcond
real(8) :: t1,t2,omp_get_wtime,pw1,pw2,epsp, dt0
real(8) , dimension(3,3) :: tag, D
real(8) , dimension(5) :: propconst
!real(8) , dimension(:,:), Allocatable ::eul
real(8) , dimension(:,:,:), Allocatable  :: F0,Fp0
real(8),  dimension(:,:), Allocatable :: S0
t1 = omp_get_wtime()
open(unit=11,file="result.txt",status='replace')
open(unit=3,file="Stress.txt",status='replace')
open(unit=8,file="eulerangles.txt",status='replace')
open(unit=13,file="Dp_cp.txt",status='replace')
open(unit=18,file="Dp_con.txt",status='replace')
open(unit=14,file="Dp2_cp.txt",status='replace')
open(unit=16,file="Grad_cp.txt",status='replace')
open(unit=4,file="hardeningrate.txt",status='replace')
call init()
write(*,*) nlines

Allocate(F0(3,3,nlines),Fp0(3,3,nlines))
Allocate(s0(nlines,12))
!!!!! This is the section where the experiment design is setup. The  variable "bryter" determines in which mode the code is run. 
!!!
!!!   bryter = 1 - Calculates point on the yield surface corresponding to the prescribed L11,L22 and plastic work and updates F0,F0p,S0
!!!                Does not write to file, Initialize an undeformed crystal
!!!
!!!   bryter = 2 - Calculates point on the yield surface given L11,L22 to a prescribed plastic work, given a prestrained and relaxed crystal is given(F0,Fp0,S0)
!!!                Writes out the calculated point in the file stress.txt 
!!!                Used when calculating instantaneous yield surfaces.    
!!!   
!!!   bryter = 3 - Calculates point on the yield surface corresponding to the prescribed L11,L22 and plastic work, RELAXES, and updates F0,F0p,S0
!!!                Used when calculating instantaneous yield surfaces, prior to setting bryter = 2.  
!!!
!!!   bryter = 4 - Used for strain path change, takes a prestrained crystal(F0,Fp0,S0) calculated using eg bryter = 1, and perform strain in a prescribed direction for a given plastic work. 
!!!     



    
k = 1
    
    propconst = (/0.0, 0.0, 0.0, 0.0, 0.1*k/)
    bryter = 5
    tag = 0
    epsp = 0


    propconst = (/0.0, 0.0, 0.0, 0.0, 0.1*k/)
Tag = 0
epsp = 0
pw1 = 0.001
bryter = 5
bcond = 2
call constexpr(k,2,bryter,bcond,pw1, tag, epsp,propconst)
!call newton(1,2,bryter,bcond,F0,Fp0,S0,pw1,propconst) 
!write(*,*) bryter

bryter = 6
pw1 = 0.001
k = 0
!call constexpr(k,2,bryter,pw1, tag, epsp)
!call newton(0,3,bryter,F0,Fp0,S0,pw1)   
!write(*,*) 'check1'
!F0 = F0i
!S0 = S0i
!Fp0 = Fp0i
!pw2 = 0.003
!F0 = 0
!Fp0 = 0
!s0 = 16
bcond = 2
part = 5

call OMP_SET_NUM_THREADS(7)
!$OMP PARALLEL PRIVATE(F0,S0,Fp0,bryter,propconst,bcond,pw1,k,slip,Cel,eul,nlines,id,R, hardening, pi, dt, tag,epsp)
!$OMP DO
do k = 0,part
    bryter = 5
    pw1 = 0.001
    bcond = 2
    write(*,*) 'start'
    propconst = (/0.0, 0.0, 0.0, 0.0, 0.1*k/)
    call newton(k,2,bryter,bcond,F0,Fp0,S0,pw1,propconst) 
    bryter = 5
    tag = 0
    epsp = 0
    call constexpr(k,2,bryter,bcond,pw1, tag, epsp,propconst)
    write(*,*) tag(1,1), tag(2,2) , k
    write(*,*) 'cycle'
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL


bcond = 2
part = 5

call OMP_SET_NUM_THREADS(7)
!$OMP PARALLEL PRIVATE(F0,S0,Fp0,bryter,propconst,bcond,pw1,k,slip,Cel,eul,nlines,id,R, hardening, pi, dt, tag,epsp)
!$OMP DO
do k = 0,part
    bryter = 5
    pw1 = 0.001
    bcond = 2
    write(*,*) 'start'
    propconst = (/0.0, 0.0, 0.0, 0.0, -0.1*k/)
    call newton(k,2,bryter,bcond,F0,Fp0,S0,pw1,propconst) 
    bryter = 5
    tag = 0
    epsp = 0
    call constexpr(k,2,bryter,bcond,pw1, tag, epsp,propconst)
    write(*,*) tag(1,1), tag(2,2) , k
    write(*,*) 'cycle'
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL


!bryter = 6
!call OMP_SET_NUM_THREADS(3)
!!$OMP PARALLEL PRIVATE(k,pw1,F0i,S0i,Fp0i,bryter)
!!$OMP DO 
!do k = -100,-30,10
!    pw1 = 10.0**(k*1.0/10)
!    F0i = F0
!    Fp0i = Fp0
!    S0i = S0
!    write(*,*) k
!    bryter = 6
!    call newton(1,4,nlines,eul,bryter,F0i,Fp0i,S0i,pw1)  
!end do
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL
!F0i = F0
!    Fp0i = Fp0
!    S0i = S0
!write(*,*) F0i
!write(*,*) F0i
!write(*,*) F0i

pw1 = 0.001
bryter = 4
!call newton(1,4,nlines,eul,bryter,F0i,Fp0i,S0i,pw1)   



pw2 = 0.0000001
part = 40


!!$OMP PARALLEL PRIVATE(k)
!!$OMP DO
!do k = 0,2*part
!    call newton(k,part,nlines,eul,2,F0i,Fp0i,S0i,pw2)
!end do
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

close(11)
close(3)
t2 = omp_get_wtime()

write(*,*) (t2-t1)/60
end program ypoint







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
    use global
    implicit none 
    real(8), dimension(3,3), intent(in) :: Tag
    real(8), dimension(3,3), intent(out) :: grad
    real(8) :: A, B, C, F, G, H, dI2,dI3, dtheta, I2, I3, theta, m = 8.8 , dsum,sum, sFi, fraction
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
    if (I3/I2**(3./2.) > 1 .or. I3/I2**(3./2.) < -1) then
        fraction = 1
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
call vec2tens(grad,n)



grad = grad*1/(2*m)*(sFi/2)**(1/m-1)
!write(*,*) grad

return
end subroutine hoshfordnormal



subroutine Yoshidamodel(Tag,D,Dp)
    use crystalplasticity
    implicit none

    real(8) , dimension(3,3) :: Dp,Dpt,Dpn,Tag,D,N, DdevT, tensprod,Ddev,Nnorm
    integer :: i,j,k,l
    real(8) :: lambdadot,lambdadottemp, theta,alpha, G, theta0, c1 = 0.3
    real(8) , dimension(6,6)              :: Chook
    real(8) , dimension(6)                :: vec
    
    Dp = 0
    Dpt = 0
    Dpn = 0
    DdevT = 0
    tensprod = 0
    Ddev = 0
    Nnorm = 0

    
    
    
    theta0 = 4.*pi/180.
    call hoshfordnormal(tag,N)
    call Elasticconstant(Chook,G)
    !write(*,*) N
    Nnorm = N/norm2(N)
    
    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3
            DdevT(i,j) = DdevT(i,j) + (1./2.*(kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k))-1./3.*(id(i,j)*id(k,l)) &
            -Nnorm(i,j)*Nnorm(k,l))*D(k,l)
                end do
            end do
        end do
    end do
    Ddev = D-1./3.*(D(1,1)+D(2,2)+D(3,3))*id
    call contract2(Ddev, N, theta)
    theta = acos(theta/norm2(Ddev)/norm2(N))
    
    if (theta > 0 .and. theta <= theta0 ) then
        alpha = 1-c1*sqrt(3.6*16/G)
    else if (theta > theta0 .and. theta <= pi/2) then
        alpha = (1-c1*sqrt(3.6*16/G))*((pi/2-theta)/(pi/2-theta0))
    else if ( theta > pi/2 .and. theta < pi ) then
        alpha = 0
    end if   
    
    Dpt = alpha*DdevT 

   vec = (/ D(1,1), D(2,2), D(3,3), 2.*D(2,3),2.*D(1,3),2.*D(1,2)/)
   vec = matmul(Chook,vec)
   call vec2tens(tensprod,vec)
   call contract2(tensprod,N,lambdadot)
   vec = (/ N(1,1), N(2,2), N(3,3), 2.*N(2,3),2.*N(1,3),2.*N(1,2)/)
   vec = matmul(Chook,vec)
   call vec2tens(tensprod,vec)
   call contract2(tensprod,N,lambdadottemp)
   
   lambdadot = lambdadot/lambdadottemp

Dpn = lambdadot*N

Dp = Dpn+Dpt


end subroutine Yoshidamodel






    
    

    subroutine yoshi(tag,D,epsp,dt0)
      use crystalplasticity
        implicit none

        real(8), dimension(3,3,3,3) :: Cep, dyadic, T, CT,Cel4
        real(8), dimension(3,3) :: N, D, tag, dtag,Nnorm, CN, NC, Ddev, tagcheck
        real(8) :: h , G,NCN, theta, c1= 0.3, theta0, alpha, sigma, sigma0, dt0, epsp, lamb, sigmacheck, dt1, consistency
        real(8) , dimension(6,6) :: Chook
        integer :: i,j,k,l,p,q
        logical :: consistent

        
        
        theta0 = pi/18.
        sigma0 = gaveps(epsp)
        h = haveps(epsp)
        dtag = 0


     
        call Elasticconstant(Chook,G)
  
        call eqvstress(tag,sigma)
        
        
    
   
   
    

   !! Convert 6x6 elastic matrix to 4th order tensor
    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                    Cel4(i,j,k,l) = Chook(1,2)*kronecker(i,j)*kronecker(k,l)+(Chook(1,1)-Chook(1,2))/2* & 
                                    (kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k))
                end do
            end do
        end do
    end do

    !write(*,*) epsp, abs(sigma0 - sigma)
    if (epsp == 0 .and. abs(sigma - sigma0) > 0.0000000000001) then
    !call eqvstress(tag,sigma)
    !write(*,*) sigma
    !write(*,*) 'Elastic'
        
        do i = 1,3
            do j = 1,3
                do k =1,3
                    do l = 1,3
                        dtag(i,j) = dtag(i,j) + Cel4(i,j,k,l)*D(l,k)
                    end do
                end do
            end do
        end do
        tagcheck = tag+dtag*dt0
    call eqvstress(tagcheck,sigmacheck)
    if (sigmacheck > sigma0) then
        consistency = abs(sigmacheck - sigma0)
        dt1 = dt0
        !write(*,*) consistency
        do while (consistency > 0.000000000001) 
        !do i = 1,20
  
            tagcheck = tag+dtag*dt1
            call eqvstress(tagcheck,sigmacheck)
            consistency = abs(sigmacheck - sigma0)
       !     write(*,*) dt1, consistency
            if (sigmacheck > sigma0 ) then 
            dt1 = (sigma0-sigma)/(sigmacheck-sigma)*dt1
            !write(*,*) dt1
            
            else
            tag = tagcheck
            sigma = sigmacheck
            end if
        
        end do
    else
        tag = tagcheck
    end if 
    !write(*,*) sigmacheck - sigma0
else
    !write(*,*) 'consistency'
    !Build T tensor ¨
    !Calculate Ce:N
    !Calculate N:Ce

    !write(*,*) tag
    call hoshfordnormal(tag,N)
    !write(*,*) N
    Nnorm = N/norm2(N)
        do i = 1,3
            do j = 1,3
                do k =1,3
                    do l = 1,3
                        T(i,j,k,l) = 1./2.*(kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k))-1./3.*(id(i,j)*id(k,l)) &
                                    -Nnorm(i,j)*Nnorm(k,l)  

                        CN(i,j) = CN(i,j)+Cel4(i,j,k,l)*N(l,k)  
                    
                        NC(i,j) = NC(i,j)+N(l,k)*Cel4(k,l,i,j)               
                    end do
                end do
            end do
        end do
       ! write(*,*) 'CN'
       ! write(*,*) 'N'
       ! 
       ! write(*,*) CN
        !write(*,*) N
        !Calculate N:Ce:N

        call contract2(N,CN,NCN)
        call contract2(D,NC,lamb)

        lamb = lamb/(NCN+sqrt(2./3.)*norm2(N)*h)
        
        epsp = epsp + sqrt(2./3.)*norm2(N)*lamb*dt0

        !Deviatoric D
        Ddev = D - id*(D(1,1)+D(2,2)+D(3,3))/3
        !Calculate theta
        call contract2(Ddev,N,theta)
        theta = acos(theta/norm2(Ddev)/Norm2(N))
    
      
        if (theta >= 0 .and. theta <= theta0 ) then
            alpha = 1-c1*sqrt(sigma0/G)
        else if (theta > theta0 .and. theta <= pi/2) then
            alpha = (1-c1*sqrt(sigma0/G))*((pi/2-theta)/(pi/2-theta0))
        else if ( theta > pi/2 .and. theta < pi ) then
            alpha = 0
        end if   
        !write(*,*) alpha
        !Calculate outer product of CN and NC

        do i = 1,3
            do j = 1,3
                do k =1,3
                    do l = 1,3
                        dyadic(i,j,k,l) = CN(i,j)*NC(k,l)
                    end do
                end do
            end do
        end do  

        !! Calculate CT

        do i = 1,3
            do j = 1,3
                do k =1,3
                    do l = 1,3
                        do p = 1,3
                            do q= 1,3
                                    CT(i,j,k,l) = CT(i,j,k,l) + alpha*Cel4(i,j,q,p)*T(p,q,k,l)
                            end do
                        end do
                    end do
                end do
            end do
        end do  
    
     Cep = Cel4 - dyadic/(NCN+sqrt(2./3.)*norm2(N)*h) - CT  
     if (theta >= pi/2) then
        !write(*,*) Cel4-Cep
     end if 

     do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                    dtag(i,j) = dtag(i,j) + Cep(i,j,k,l)*D(l,k)
                end do
            end do
        end do
    end do
    tag = tag+dtag*dt0
end if 
   ! call eqvstress(tag,sigma)
   ! write(*,*) sigma,gaveps(epsp), sigma-gaveps(epsp)
   !
   ! write(*,*) dtag
   ! write(*,*) 


    return
    end subroutine yoshi

    subroutine constexpr(l,part,bryter,bcond,strain,tag,epsp,propconst)
        use crystalplasticity
        implicit none

        real(8), dimension(3,3) :: tag, tagi
        real(8) :: epsp, epspi, strain
        real(8) , dimension(4)  :: sigma, offdl
        real(8), dimension(5) :: propconst, offdl2, sigma2, IPIV2
        integer, dimension(5) :: pos1, pos2
        real(8) :: dt0,gammaskrank, dl
        real(8) , dimension(3,3)  :: Lb, Tagb, La, Dp,tagc
        
        integer ::  switch , p,k,h,l, bryter,secit,part
        real(8) , dimension(4,4)  :: Jacob, Jinv
        
        integer :: LDA = 4,NRHS = 1, Info,  minl, maxl,nit,bcond
        integer , dimension(4) :: IPIV
        real(8) , dimension(5,5)  :: Jacob2, Jinv2
        real(8) :: pwpercision
        gammaskrank = 0
        pwpercision = 0.0000000001
        secit = 0
        pos1 = (/1, 1, 2, 3, 2/)
        pos2 =(/2, 3, 3, 3, 2/)
        dl = 0.001
    
    ! Sets the initial velocity gradient for the given boundary condition
    select case (bcond)
    case(1)
                    La = 0 
            La(1,1) = cos(pi*l/part)
            La(2,2) = sin(pi*l/part)
            La(3,3) =-0.3*(La(1,1)+La(2,2))
            La(1,2) = 0
            La(2,1) = 0
            La(1,3) = 0
            La(3,1) = 0
            La(2,3) = 0
            La(3,2) = 0
    case(2) 
        La(1,1) = 1/sqrt(1.0+propconst(5)**2)
        La(2,2) = propconst(5)**2/sqrt(1.0+propconst(5)**2)
        La(1,2) = 0
        La(2,1) = 0
        La(1,3) = 0
        La(3,1) = 0
        La(2,3) = 0
        La(3,2) = 0
        La(3,3) = -1.0/3.0*(La(1,1)+La(2,2))
    end select    
    
            dt0 = dt
iter: do while (switch < 100000)



Select case (bcond)
case (1)    
!!! Iteration to ensure boundary condition is satisfied at througout all timesteps. 
    nit = 0   
    boundarycond: do  while (nit < 100)  
    tagi = tag
    epspi = epsp
       call yoshi(Tagi,La,epspi,dt0)
 
       do h = 1,4
        sigma(h) = Tagi(pos1(h),pos2(h))
       end do
  
      ! write(*,*) sigma
         minl = minloc(sigma, DIM = 1)
         maxl = maxloc(sigma, DIM = 1)
         
         if (abs(sigma(minl)) < 0.00000000001 .and. abs(sigma(maxl)) < 0.00000000001) then
            !write(*,*) sigma(4) , Tagi(1,1), Tagi(2,2), epspi
            exit boundarycond
         end if 
!!!!!! Calculate Jacobian matrix in order to satisfy boundary condition
       do k = 1,4
            do p = 1,4
                tagb = tag
                epspi = epsp
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) + dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl 
                end if

                call yoshi(Tagb,Lb,epspi,dt0)
            jacob(k,p) = (Tagb(pos1(k),pos2(k))-Tagi(pos1(k),pos2(k)))/dl
            
            end do
        end do
        !write(*,*) tag
        !write(*,*) tagb
        Jinv = jacob
        offdl = -sigma

        call dgesv(LDA,NRHS,Jinv,LDA,IPIV,offdl,lda , Info)
        La(1,2) = La(1,2) + offdl(1)
        La(2,1) = La(2,1) + offdl(1)
        La(1,3) = La(1,3) + offdl(2)
        La(3,1) = La(3,1) + offdl(2)
        La(2,3) = La(2,3) + offdl(3)
        La(3,2) = La(3,2) + offdl(3)
        La(3,3) = La(3,3) + offdl(4)
        
        minl = minloc(sigma, DIM = 1)
        maxl = maxloc(sigma, DIM = 1)
       
       

   nit = nit+1
   !write(*,*) nit
    end do boundarycond
   
case (2)
!dl = 0.001  

    nit = 0
   boundary2: do while (nit < 100)
   
   tagi = tag
   epspi = epsp
   
   
   call yoshi(Tagi,La,epspi,dt0)
  ! write(*,*) tagi
    do h = 1,5
        sigma2(h) = Tagi(pos1(h),pos2(h))-propconst(h)*Tagi(1,1)
    end do
       !write(*,*) sigma2
         minl = minloc(sigma2, DIM = 1)
         maxl = maxloc(sigma2, DIM = 1)
         if (abs(sigma2(minl)) < 0.00000000001 .and. abs(sigma2(maxl)) < 0.00000000001) then
          ! write(*,*) sigma2 , epspi
          ! write(*,*) tagi
            exit boundary2
         end if 
     do k = 1,5
            do p = 1,5
                tagb = tag
                epspi = epsp
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) + dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl 
                end if

                call yoshi(Tagb,Lb,epspi,dt0)

                tagc = tag
                epspi = epsp
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) - dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) - dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) - dl 
                end if
                call yoshi(Tagc,Lb,epspi,dt0)

            jacob2(k,p) = ((Tagb(pos1(k),pos2(k))-propconst(k)*Tagb(1,1))-(Tagc(pos1(k),pos2(k))-propconst(k)*tagc(1,1)))/dl/2
            end do
        end do
        Jinv2 = jacob2
        offdl2 = -sigma2
        call dgesv(5,1,Jinv2,5,IPIV2,offdl2, 5 , Info)
      
        La(1,2) = La(1,2) + offdl2(1)
        La(2,1) = La(2,1) + offdl2(1)
        La(1,3) = La(1,3) + offdl2(2)
        La(3,1) = La(3,1) + offdl2(2)
        La(2,3) = La(2,3) + offdl2(3)
        La(3,2) = La(3,2) + offdl2(3)
        La(3,3) = La(3,3) + offdl2(4)
        La(2,2) = La(2,2) + offdl2(5)
      
        nit = nit+1

        !write(*,*) La
        !write(*,*) tagi
        !write(*,*) tagb
        !write(*,*) 
        !write(*,*) Jacob2(1,1:5)
        !write(*,*) Jacob2(2,1:5)
        !write(*,*) Jacob2(3,1:5)
        !write(*,*) Jacob2(4,1:5)
        !write(*,*) Jacob2(5,1:5)
        !write(*,*) 


        
        !call sleep(3)
   end do boundary2
end select
!call sleep(2)
if (bryter == 1 .or. bryter == 5 .or. bryter == 4) then
       
        
       
    if (epspi > strain .and. abs((epspi - strain)/strain) > pwpercision) then
        if (strain == 0) then
            dt0 = dt0/2
        else
        dt0 = (strain-epsp)*dt0/(epspi-epsp)*0.7
        end if
        secit = secit +1
        !write(*,*) dt0, epspi, epsp
        switch = switch + 1
        if (secit > 30) then 
            write(*,*) epspi
            write(*,*) 'early exit'
            exit iter
        end if 
        cycle iter
    end if 
    
    tag = tagi
    !write(*,*) epspi, nit, dt0, l
    epsp = epspi
    !write(*,*) epsp
    if (bryter == 5) then
        if (epsp > gammaskrank) then
            write(8,*) bryter, epspi
            
        write(11,*) Tag(1,3), Tag(1,2), Tag(2,3), Tag(3,3) ,epsp
        call Yoshidamodel(tag,La,Dp)
        
        write(18,*) Tag(1,1), Tag(2,2) , Dp(1,1), Dp(2,2) 
        gammaskrank = gammaskrank + 0.0000001
        
        end if
    end if 

    if (abs((epsp -strain)/strain) <= pwpercision) then
       exit iter
    end if

else if (bryter == 2 .or. bryter == 6) then
    
    if (strain /= 0) then
    if (epspi > strain .and. abs((epspi - strain)/strain) > pwpercision) then
        if (strain == 0) then
            dt0 = dt0/2
        else
        dt0 = (strain-epsp)*dt0/(epspi-epsp)*0.5
        end if
        write(*,*) dt0, sigma, epspi-strain
        switch = switch +1
        secit = secit +1
        if (secit > 15) then 
            exit iter
        end if 
        cycle iter
    end if 
    end if

    if (strain == 0) then
        if (epspi > strain .and. abs(epspi) > pwpercision) then
        dt0 = dt0/2
        !dt0 = (pw-epsp)*dt0/(epspi-epsp) 
        switch = switch +1
        secit = secit +1
        if (secit > 10) then 
            exit iter
        end if 
        cycle iter
        end if
    end if

    tag = tagi
    epsp = epspi 
    write(*,*) epsp

    if (bryter == 6) then
        if (epsp > gammaskrank) then
            write(8,*) bryter, epspi
            
        write(11,*) Tag(1,1), tag(2,2), Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3),  epsp
            call Yoshidamodel(tag,La,Dp)

        gammaskrank = gammaskrank + 0.00000005
        write(13,*) Tag(1,1), Tag(2,2) , Dp(1,1), Dp(2,2) 
        !write(8,*) Tag(1,1),Tag(2,2), Dp(1,1)/sqrt(Dp(1,1)**2+Dp(2,2)**2),Dp(2,2)/sqrt(Dp(1,1)**2+Dp(2,2)**2.)
        
        end if
    end if 
    
    if (strain /= 0 .and. abs((epsp -strain)/strain)<= pwpercision) then
        exit iter
    else if (strain == 0 .and. epsp <= pwpercision .and. epsp > 0) then
       ! write(*,*) 'check'
        exit iter

    end if
end if 
switch = switch +1 
    end do iter

end subroutine constexpr

subroutine elasticsolution(D,tag)
use crystalplasticity
    implicit none
real(8) sigma, sigma0, G, dt0, consistency, sigma1, dt1, dtt
real(8), dimension(6,6) :: Chook 
real(8), dimension(3,3) :: D, tag, dtag, tag1
real(8), dimension(3,3,3,3) :: Cel4
integer :: i,j,k,l

sigma= 0
sigma0 = gaveps(sigma)

    call elasticconstant(Chook,G)
     
     !! Convert 6x6 elastic matrix to 4th order tensor
    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                    Cel4(i,j,k,l) = Chook(1,2)*kronecker(i,j)*kronecker(k,l)+(Chook(1,1)-Chook(1,2))/2* & 
                                    (kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k))
                end do
            end do
        end do
    end do
    
    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                    dtag(i,j) = dtag(i,j) + Cel4(i,j,k,l)*D(l,k)
                end do
            end do
        end do
    end do
 tag1 = 0
 sigma = 0
 dt0 = 0
 dt1 = dt 
 consistency = abs(sigma-sigma0)
 do while (sigma < sigma0)
    tag1 = tag1+dtag*dt1
    call eqvstress(tag1,sigma)
 end do
 
 tag1 = tag1-dtag*dt1
 call eqvstress(tag1,sigma)
 

do while (consistency > 0.000001) 
  
    tag1 = tag+dtag*dt1
    call eqvstress(tag1,sigma1)
    consistency = abs(sigma1 - sigma0)
    if (sigma1 > sigma0 ) then 
    dt1 = (sigma0-sigma)/(sigma1-sigma)*dt1
    
    else
    tag = tag1
    end if

end do
 tag = tag1
 
end subroutine