program ypoint
use global
    implicit none

integer :: part,bryter, k
real(8) :: t1,t2,omp_get_wtime,pw1,pw2
!real(8) , dimension(:,:), Allocatable ::eul
real(8) , dimension(:,:,:), Allocatable  :: F0,Fp0
real(8),  dimension(:,:), Allocatable :: S0
real(8), dimension(6,6) :: Chook
t1 = omp_get_wtime()
open(unit=11,file="result.txt",status='replace')
open(unit=3,file="Stress.txt",status='replace')
open(unit=8,file="eulerangles.txt",status='replace')
open(unit=13,file="Dp.txt",status='replace')
open(unit=14,file="Dp2.txt",status='replace')
open(unit=16,file="Grad.txt",status='replace')
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






pw1 = 0.003
bryter = 5
call newton(1,2,bryter,F0,Fp0,S0,pw1) 
!write(*,*) bryter

bryter = 6
pw1 = 0.005
call newton(0,3,bryter,F0,Fp0,S0,pw1)   
!write(*,*) 'check1'
!F0 = F0i
!S0 = S0i
!Fp0 = Fp0i
!pw2 = 0.003
!part = 200
!call OMP_SET_NUM_THREADS(7)
!!$OMP PARALLEL PRIVATE(k,F0,S0,Fp0,bryter)
!!$OMP DO
!do k = 0,2*part
!    bryter = 7
!    call newton(k,part,bryter,F0,Fp0,S0,pw2)
!end do
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL



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



subroutine newton(k,part,bryter,F0i,Fp0i,S0i,pw)
    use global
    implicit none
real(8) , dimension(3,3)  :: La, Tag, Dp
real(8) :: pw
integer ::  k, part,teller,bry
integer  :: bryter
!real(8), dimension(nlines,3), intent(in) :: eul
real(8) , dimension(3,3,nlines)  :: Fp0i,F0i
real(8),  dimension(nlines,12) :: S0i

bry = 0
if (bryter == 4) then 
    bry = 1
else if (bryter == 3) then
    bryter = 1
    bry = 1
else if (bryter == 7) then
    bryter = 1
    bry = 2
else if (bryter == 8) then
    bryter = 4
    bry = 5
end if





teller = 0

!Initialize La
!strainrate
La = 0 
La(1,1) = cos(pi*k/part)
La(2,2) = sin(pi*k/part)
La(3,3) =-0.3*(La(1,1)+La(2,2))
La(1,2) = 0
La(2,1) = 0
La(1,3) = 0
La(3,1) = 0
La(2,3) = 0
La(3,2) = 0

    call taylor(La,Tag,bryter,F0i,Fp0i,S0i,pw,Dp)

    if (bry == 1) then ! If one want relaxed state. 
        !Performes relaxation
        call taylor(La,Tag,3,F0i,Fp0i,S0i,pw,Dp)
    end if
if (bryter == 2) then
write(3,*) Tag(1,1), Tag(2,2), k
!write(11,*) Tag(1,1), Tag(2,2)
write(*,*) Tag(1,1), Tag(2,2), k
end if 
if (bry == 2) then
    write(3,*) Tag(1,1), Tag(2,2), k
    !write(11,*) Tag(1,1), Tag(2,2)
    write(*,*) Tag(1,1), Tag(2,2), k
    bryter = 5
end if 
if (bry == 5) then
    bryter = 6
    write(11,*) Tag(1,1), Tag(2,2), Dp(1,1)/sqrt(dp(1,1)**2+dp(2,2)**2), dp(2,2)/sqrt(dp(1,1)**2+dp(2,2)**2)
    !write(11,*) Tag(1,1), Tag(2,2)
    write(*,*) Tag(1,1), Tag(2,2), Dp(1,1)/sqrt(dp(1,1)**2+dp(2,2)**2), dp(2,2)/sqrt(dp(1,1)**2+dp(2,2)**2)
    end if 
!end do
    if (bryter == 4) then 

    end if
end subroutine newton


subroutine taylor(La,Tag,bryter,F0i,Fp0i,S0i,pw,Dp)
    use global
        implicit none
    
     
    !Declear all variables
 
    real(8) :: phi1, Phi, phi2,dt0, sigmaeq, gammatot,pw,gammatoti,gammaskrank,dot, dl, nor1,nor2
    real(8) , dimension(3,3)  :: grad,T0 , Dp, Lb, Tagb, La0,Dp2
                        
    real(8) , dimension(3,3) :: La
    real(8) , dimension(3,3), intent(out)  :: Tag
    real(8) , dimension(3,3,nlines) ::  F0, Fp0,Fp0i,Fp0int,F0i,F0int, Tagc, Tagcint
    real(8),  dimension(nlines,12) :: s0,S0i,S0in
    integer :: i, switch , o,p,k,h, bryter,secit
    real(8) , dimension(4,4)  :: Jacob, Jinv
    real(8) , dimension(4)  :: sigma, offdl
    integer :: LDA = 4,NRHS = 1, Info,  minl, maxl,nit,bcond=1
    integer , dimension(4) :: IPIV
    real(8), dimension(5) :: propconst, offdl2, sigma2, IPIV2
    real(8) , dimension(5,5)  :: Jacob2, Jinv2
    integer, dimension(5) :: pos1, pos2
    real(8) :: pwpercision
    
    ! The percision of the plastic work given in relative fraction
    pwpercision = 0.0000000001
    !Timeincrement
    !dt0 = 0.0000001

    pos1 = (/1, 1, 2, 3, 2/)
    pos2 =(/2, 3, 3, 3, 2/)
    dl = 0.0000001
    La0 = La
    !Define velocity gradient
    !strainrate
        gammatot = 0
    !Copies of the input variables, in order not to update the initial condition when calculating instantaneous yield surface.
        S0 = s0i 
        Fp0 = Fp0i  
        F0 = F0i  

        propconst = (/0.0, 0.0, 0.0, 0.0, -10.0/)
    
    
 
   
    if (bryter == 1 .or. bryter == 5) then
        s0(:,1:12) = 16
        initalize: do o = 1,nlines ! Loop to calculate values at time 0
        
    !read(9,*) phi1, Phi, phi2
        phi1 = eul(o,1)
        phi = eul(o,2)
        phi2 = eul(o,3)
     !Deformation gradient at time 0 is I
        F0(1:3,1:3,o) = id
    !Cauchy Stress at time 0, Crystal orientation
        T0 = 0
    !Prior to all deformation Fp0 = I, Crystal orientation
        Fp0(1:3,1:3,o) = id
    
    !!!!!! Variable assignments end
    ! Calculate rotation matrix
    call rotmat(phi1,Phi,phi2,R(1:3,1:3,o))
    
   
    end do initalize
    
else if (bryter == 2 .or. bryter == 4 .or. bryter == 6) then
    initalize2: do o = 1,nlines ! Loop to calculate values at time 0
    phi1 = eul(o,1)
    phi = eul(o,2)
    phi2 = eul(o,3)
    
    !!!!!! Variable assignments end
    ! Calculate rotation matrix
    call rotmat(phi1,Phi,phi2,R(1:3,1:3,o))
    
    end do initalize2
    
    
else if (bryter == 3 ) then
    initalize3: do o = 1,nlines ! Loop to calculate values at time 0
    phi1 = eul(o,1)
    phi = eul(o,2)
    phi2 = eul(o,3)
    
    !!!!!! Variable assignments end
    ! Calculate rotation matrix
    call rotmat(phi1,Phi,phi2,R(1:3,1:3,o))
    
    
    end do initalize3
    write(*,*) bryter
    end if 



    gammaskrank = 0.00001
 
    secit = 0
    dt0 = dt
    switch = 1
    if (bryter == 3) then
        goto 13
    end if


    iter: do while (switch < 5000000)     
    
       
   
   
Select case (bcond)
case (1)    
!!! Iteration to ensure boundary condition is satisfied at througout all timesteps. 
    nit = 0   
    boundarycond: do  while (nit < 10)  
       ! write(*,*) La
    nor1 = norm2(La)
       call timestep(Tag, Dp, La, gammatot, gammatoti, Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
       write(*,*) sigma , Tag(1,1), Tag(2,2), gammatoti
       do h = 1,4
        sigma(h) = Tag(pos1(h),pos2(h))
       end do
      ! write(*,*) sigma
         minl = minloc(sigma, DIM = 1)
         maxl = maxloc(sigma, DIM = 1)
         if (abs(sigma(minl)) < 0.00000000001 .and. abs(sigma(maxl)) < 0.00000000001) then
            
         !   write(*,*) tag
            exit boundarycond
         end if 
!!!!!! Calculate Jacobian matrix in order to satisfy boundary condition
       do k = 1,4
            do p = 1,4
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) + dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl 
                end if

                call timestep(Tagb, Dp, Lb, gammatot, gammatoti, Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
            jacob(k,p) = (Tagb(pos1(k),pos2(k))-Tag(pos1(k),pos2(k)))/dl
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
       
       nor2 = norm2(La)

   nit = nit+1
   !write(*,*) nit
    end do boundarycond
   
case (2)
    nit = 0
   
    !La(1,1) = 1/sqrt(1.0+propconst(1)**2)
    !La(2,2) = propconst(1)**2/sqrt(1.0+propconst(1)**2)
    !La(1,2) = 0
    !La(2,1) = 0
    !La(1,3) = 0
    !La(3,1) = 0
    !La(2,3) = 0
    !La(3,2) = 0
    !La(3,3) = -1.0/3.0*(La(1,1)+La(2,2))
   boundary2: do while (nit < 20)
    call timestep(Tag, Dp, La, gammatot, gammatoti , Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
    do h = 1,5
        sigma2(h) = Tag(pos1(h),pos2(h))-propconst(h)*Tag(1,1)
    end do
      ! write(*,*) sigma2
    if (gammatoti > pw) then
      !  write(*,*) sigma2 , Tag(1,1), Tag(2,2), gammatoti
    end if
         minl = minloc(sigma2, DIM = 1)
         maxl = maxloc(sigma2, DIM = 1)
         if (abs(sigma2(minl)) < 0.00000000001 .and. abs(sigma2(maxl)) < 0.00000000001) then
           write(*,*) sigma2 , Tag(1,1), Tag(2,2), gammatoti
         !   write(*,*) tag
            exit boundary2
         end if 
     do k = 1,5
            do p = 1,5
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) + dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl 
                end if

                call timestep(Tagb, Dp, Lb, gammatot, gammatoti , Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
            jacob2(k,p) = ((Tagb(pos1(k),pos2(k))-propconst(k)*Tagb(1,1))-(Tag(pos1(k),pos2(k))-propconst(k)*tag(1,1)))/dl
            end do
        end do
        Jinv2 = jacob2
        offdl2 = -sigma2
        !write(*,*) offdl2
        call dgesv(5,1,Jinv2,5,IPIV2,offdl2, 5 , Info)
       ! write(*,*) info
       ! write(*,*) offdl2
        La(1,2) = La(1,2) + offdl2(1)
        La(2,1) = La(2,1) + offdl2(1)
        La(1,3) = La(1,3) + offdl2(2)
        La(3,1) = La(3,1) + offdl2(2)
        La(2,3) = La(2,3) + offdl2(3)
        La(3,2) = La(3,2) + offdl2(3)
        La(3,3) = La(3,3) + offdl2(4)
        La(2,2) = La(2,2) + offdl2(5)
      
        nit = nit+1
   end do boundary2
end select
  !if (gammatoti > pw) then 
  !  write(*,*) 'reached plastic work'
  !  exit iter
  ! end if
   
   
   

   
   
       !write(*,*) Tag
    if (mod(switch,5000-1)==0) then 
        !   write(8,*) phi1, Phi, phi2
       end if 
    call eqvstr(Tag,sigmaeq)
    !write(10,*) sigmaeq

    if (bryter == 1 .or. bryter == 5 .or. bryter == 4) then
       
        
       
        if (gammatoti > pw .and. abs((gammatoti - pw)/pw) > pwpercision) then
            if (pw == 0) then
                dt0 = dt0/2
            else
            dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot)*0.5
            end if
            write(*,*) dt0, sigma, gammatoti-pw, nor1, nor2
            switch = switch +1
            secit = secit +1
            if (secit > 30) then 
                write(*,*) gammatoti
                write(*,*) 'early exit'
                exit iter
            end if 
            cycle iter
        end if 
        
        F0 = F0int
        Fp0 = Fp0int
        Tagc = Tagcint
        s0 = s0in
        gammatot = gammatoti 
        
        if (bryter == 5) then
            if (gammatot > gammaskrank) then
                write(8,*) bryter, gammatoti
                call contract2(La,Dp,dot)
                dot = dot/norm2(La)/norm2(Dp)
            write(11,*) Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3) , dot , acos(dot), gammatot
            
            gammaskrank = gammaskrank + 0.00001
            !write(8,*) Tag(1,1),Tag(2,2), Dp(1,1)/sqrt(Dp(1,1)**2+Dp(2,2)**2),Dp(2,2)/sqrt(Dp(1,1)**2+Dp(2,2)**2.)
            write(8,*) 'La'
            write(8,*) La/norm2(La)
            call hoshfordnormal(tag,grad)
            write(8,*) 'Dp'
            write(8,*) Dp/norm2(Dp)
            write(8,*) 'Gradient'
            write(8,*) grad/norm2(grad)
            write(8,*) 'Dp-Yoshida'
            call Yoshidamodel(Tag,La,Dp2)
            write(8,*) Dp2/norm2(Dp2)
            write(8,*) 'Dp:La'
            write(8,*) dot
            call contract2(Dp,grad,dot)
            write(8,*) 'Dp:Grad'
            write(8,*) dot/norm2(Dp)/norm2(grad)
            write(8,*) 'La:Grad'
            call contract2(La,grad,dot)
            write(8,*) dot/norm2(La)/norm2(grad)
           
            
            call contract2(Dp,Dp2,dot)
            write(8,*) 'Dp:DpY'
            write(8,*) dot/norm2(Dp2)/norm2(Dp)
            write(8,*)
            write(8,*)

            write(13,*) Tag(1,1), Tag(2,2), Dp(1,1)  /sqrt(  Dp(1,1)**2+Dp(2,2)**2),Dp(2,2)/sqrt(Dp(1,1)**2+Dp(2,2)**2)
            write(14,*) Tag(1,1), Tag(2,2), Dp2(1,1) /sqrt( Dp2(1,1)**2+Dp2(2,2)**2), Dp2(2,2)/sqrt(Dp2(1,1)**2+Dp2(2,2)**2)
            write(16,*) Tag(1,1), Tag(2,2), Grad(1,1)/sqrt(Grad(1,1)**2+Grad(2,2)**2),Grad(2,2)/sqrt(Grad(1,1)**2+Grad(2,2)**2)
            end if
        end if 

        if (abs((gammatot -pw)/pw) <= pwpercision) then
            Fp0i = Fp0
            F0i  = F0
            S0i = S0   
            call Yoshidamodel(Tag,La,Dp2)
           ! call hoshfordnormal(tag,grad)
           ! call contract2(Dp,grad,dot)
           ! write(*,*) dot/norm2(Dp)/norm2(grad)
           ! write(*,*) grad
           !write(*,*) Dp
           exit iter
        end if
    
    else if (bryter == 2 .or. bryter == 6) then
        
        if (pw /= 0) then
        if (gammatoti > pw .and. abs((gammatoti - pw)/pw) > pwpercision) then
            if (pw == 0) then
                dt0 = dt0/2
            else
            dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot)*0.5
            end if
            write(*,*) dt0, sigma, gammatoti-pw
            switch = switch +1
            secit = secit +1
            if (secit > 15) then 
                exit iter
            end if 
            cycle iter
        end if 
        end if

        if (pw == 0) then
            if (gammatoti > pw .and. abs(gammatoti) > pwpercision) then
            dt0 = dt0/2
            !dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot) 
            switch = switch +1
            secit = secit +1
            if (secit > 10) then 
                exit iter
            end if 
            cycle iter
            end if
        end if

        F0 = F0int
        Fp0 = Fp0int
        Tagc = Tagcint
        s0 = s0in
        gammatot = gammatoti 

        if (bryter == 6) then
            if (gammatot > gammaskrank) then
                write(8,*) bryter, gammatoti
                call contract2(La,Dp,dot)
                dot = dot/norm2(La)/norm2(Dp)
            write(11,*) Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3) , dot , acos(dot), gammatot
          
            gammaskrank = gammaskrank + 0.000005
            !write(8,*) Tag(1,1),Tag(2,2), Dp(1,1)/sqrt(Dp(1,1)**2+Dp(2,2)**2),Dp(2,2)/sqrt(Dp(1,1)**2+Dp(2,2)**2.)
            write(8,*) 'La'
            write(8,*) La/norm2(La)
            call hoshfordnormal(tag,grad)
            write(8,*) 'Dp'
            write(8,*) Dp/norm2(Dp)
            write(8,*) 'Gradient'
            write(8,*) grad/norm2(grad)
            write(8,*) 'Dp-Yoshida'
            call Yoshidamodel(Tag,La,Dp2)
            write(8,*) Dp2/norm2(Dp2)
            write(8,*) 'Dp:La'
            write(8,*) dot
            call contract2(Dp,grad,dot)
            write(8,*) 'Dp:Grad'
            write(8,*) dot/norm2(Dp)/norm2(grad)
            write(8,*) 'La:Grad'
            call contract2(La,grad,dot)
            write(8,*) dot/norm2(La)/norm2(grad)
           
            
            call contract2(Dp,Dp2,dot)
            write(8,*) 'Dp:DpY'
            write(8,*) dot/norm2(Dp2)/norm2(Dp)
            write(8,*)
            write(8,*)

            write(13,*) Tag(1,1), Tag(2,2), Dp(1,1)  /sqrt(  Dp(1,1)**2+Dp(2,2)**2),Dp(2,2)/sqrt(Dp(1,1)**2+Dp(2,2)**2)
            write(14,*) Tag(1,1), Tag(2,2), Dp2(1,1) /sqrt( Dp2(1,1)**2+Dp2(2,2)**2), Dp2(2,2)/sqrt(Dp2(1,1)**2+Dp2(2,2)**2)
            write(16,*) Tag(1,1), Tag(2,2), Grad(1,1)/sqrt(Grad(1,1)**2+Grad(2,2)**2),Grad(2,2)/sqrt(Grad(1,1)**2+Grad(2,2)**2)
            end if
        end if 
        
        if (pw /= 0 .and. abs((gammatot -pw)/pw)<= pwpercision) then
            exit iter
        else if (pw == 0 .and. gammatot <= pwpercision .and. gammatot > 0) then
           ! write(*,*) 'check'
            exit iter

        end if
    end if 
        switch = switch +1
        end do iter
       
       
       !For unloading, iterate a given set of timesteps backwards in order to relax the crystal
        13 continue
        if (bryter == 3) then
            La = 0 
            La(1,1) = Tag(1,1)/sqrt(Tag(1,1)**2+Tag(2,2)**2)
            La(2,2) = Tag(2,2)/sqrt(Tag(1,1)**2+Tag(2,2)**2)
            La(3,3) =-0.3*(La(1,1)+La(2,2))
            La(1,2) =0
            La(2,1) = 0
            La(1,3) = 0
            La(3,1) = 0
            La(2,3) = 0
            La(3,2) = 0
       
        do i = 1,10
            
            call timestep(Tag, Dp, -La, gammatot, gammatoti , Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
            
            !if (sum(abs(x)) == 0) then 
            F0 = F0int
            Fp0 = Fp0int
            Tagc = Tagcint
            s0 = s0in
            gammatot = gammatoti 
            !else if (sum(abs(x)) /= 0) then
            !    La = -La
            !end if 
        end do
        
        
        Fp0i = Fp0
        F0i  = F0
        S0i = S0   
  
        end if



    return
    end subroutine taylor



subroutine timestep(Tag, Dp, La, gammatot, gammatoti , Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
    use global
    implicit none
    
    !Declear all variables
 
    real(8) :: phi1, Phi, phi2, det,dt0, cons, gammatot,gammatoti,h
    real(8) , dimension(3,3)  :: F1, Fp1, Fp1inv, Fe1, Fetr,Fp0inv, Ctr,T1,Ttr,Etr , & 
                         Schmid,TS,Tint, CS, Tst, Fpint ,Rn,  Rtest,Dpc,Dp
    real(8) , dimension(3,3), intent(in)  :: La
    real(8) , dimension(3,3), intent(out)  :: Tag
    real(8) , dimension(3,3,nlines) ::   Tagcint, Lc
    real(8) , dimension(3,3,nlines), intent(in) :: F0, Fp0
    real(8) , dimension(3,3,nlines), intent(out) :: F0int, Fp0int
    real(8),  dimension(nlines,12) :: s0,S0in
    real(8) , dimension(12) :: tautr,s1,tau, consis, s0int
    real(8) , dimension(3) :: m,n
    logical, dimension(12) :: PA, Active
    logical :: consise
    integer :: i,countera,ij, numact,j,q
    double precision, dimension(12):: x
    
    Tag = 0
    Dp = 0
    gammatoti = gammatot
    
    grains: do j = 1,nlines   !Iterate over all grains                                                                                                                                                                                   
    !integrate velocity gradient to form deformation gradient
     ! Rotate velocity gradient to crystal coordinate 
    Lc(1:3,1:3,j) = matmul(R(1:3,1:3,j),matmul(La,transpose(R(1:3,1:3,j))))
    Dpc = 0
    
    F1 = matmul(F0(1:3,1:3,j),(id + Lc(1:3,1:3,j)*dt0))
    !Step 1, Calculate trial elastic strain
    
    call minv3(Fp0(1:3,1:3,j),Fp0inv) ! Calculate inverse of plastic deformation gradient at time 0
    Fetr = matmul(F1,Fp0inv)
    Ctr = matmul(transpose(Fetr),Fetr)
    Etr = (Ctr-id)/2
    
    !Step 2, calculate trial stress
    
    call voigt(Etr,Ttr) !Calculate C[Etr]=Ttr
    
    !Step 3 and 4, Calcualte resolved shear stress on each slip system and detrermine potentially active slip systems (P.A)
    
    do i = 1,12
    call slipsys(slip,Schmid,m,n,i)
    TS = matmul(transpose(Schmid),Ttr)
    tautr(i) = TS(1,1)+TS(2,2)+TS(3,3)
    PA = abs(tautr)-s0(j,1:12) > 0 
    end do
    
    numact = count(PA)
    !Step 5, Calculate shear increments by determination of A matrix
    
    fiveten: do ij = 1,20    ! Iteration between step 5 and 10
    ! Calculate the slip increments on the potentially active slip systems. 
    cons = 0
    x = 0 
    Active = PA
    s0int = s0(j,1:12)
    if (count(PA) > 0) then
    call sincr(Active,x,tautr,s0int,Ctr) 
    end if 
    
    !step 6, Update plastic deformation gradient
    Fp1 = 0
    countera = 1
    do i = 1,12
        if (Active(i)) then 
            call slipsys(slip,Schmid,m,n,i)
            Fpint =tautr(i)/abs(tautr(i))*x(countera)*Schmid
        Fp1 = Fp1+ Fpint
        countera = countera+1
        end if 
    end do
    Fp1 = matmul(id,Fp0(1:3,1:3,j))+ matmul(Fp1,Fp0(1:3,1:3,j))
    
    ! step 7, Check if determinant is 1, if not normalize
    
    call deter(Fp1,det)
    if (abs(det-1.0) > 0.000001) then 
        Fp1 = det**(-1.0/3.0)*Fp1
    end if
    
    !Step 8, calculate elastic deformation gradient and updated stress in reference state.
    
    !Elastic deformation gradient
    call minv3(Fp1,Fp1inv)
    Fe1 = matmul(F1,Fp1inv)
    
    !Calculate tau star
    Tst = Ttr
    countera = 1
    do i = 1,12
        if (Active(i)) then 
            call slipsys(slip,Schmid,m,n,i)
            CS = matmul(Ctr,Schmid)
            call voigt((CS + transpose(CS))/2,Tint)
        Tst = Tst-Tint*tautr(i)/abs(tautr(i))*x(countera)
        countera = countera+1
        end if 
    end do
    
    !Step 9, Update stress in rotated system and hardening parameter
    
    !Update stress
    call deter(Fe1,det)
    
    T1 = matmul(Fe1, matmul(Tst/det,transpose(Fe1)))
    
    
    
    !Calculate increase in critical shear rate.
    s1 = s0(j,1:12)
    
     do i = 1,12
        countera = 1
        do q = 1,12
            if (Active(q)) then 
            call hparam(i,q,s0int,h,slip)
            s1(i) = s1(i) + h*x(countera)
            countera = countera+1
            end if
        end do
    end do
    
    !Step 10, check consistency
    
    !Calculates resolved shear stress using updated stress in addition ti the vector for conistency and the final active slip systems
    do i = 1,12
        call slipsys(slip,Schmid,m,n,i)
        TS = matmul(transpose(Tst),Schmid)
        tau(i) = TS(1,1)+TS(2,2)+TS(3,3)
       consis(i) = abs(tau(i))-s1(i)
       if (Active(i)) then
        cons = cons + consis(i)**2
       end if  
     
    end do
    
    cons = sqrt(cons)
    
    !If the potentially active and final active slip sys are different the algorithm goes back to step 5, if not it will accept the solution. 
    consise = .true.
    if (count(PA) /= 0) then
    do i = 1,12
            if (.not.PA(i) .and. abs(tau(i)) > s1(i) -0.0000001 )   Then ! Checks if non potenially active system is activated.
                    PA(i) = .true.
                    consise = .false.
                    cycle fiveten
                end if 
     end do
    
     do i = 1,12
        if (Active(i) .and. abs(tau(i)) < s1(i) -0.0000001 ) then ! If the active slip systems is not consistent 
            consise = .false.
           PA(i) = .false.
           cycle fiveten
        end if 
    end do
    
    end if 
    
    if (cons > 0.00001) then
        write(*,*) 'timestep to large, consistency not achieved 1'
        write(*,*) consis
        write(*,*) x
        dt0 = dt0/10
        cycle fiveten
    end if
    
    if (consise) then
    exit fiveten
    end if
    !write(7,*) ij
    
    
    end do fiveten! loop ij 
    !Step 11, compute updated texture. Preferentially write subroutine to update euler angles. new R = (I-Fe1)*R0 ? 
    
    !call decomp(Fe1,Rtest)
    call POLAR(Fe1,Rtest)
    
    Rn = matmul(transpose(Rtest),R(1:3,1:3,j))
    call eulang(Rn,phi1,Phi,phi2)
    
    countera = 1
    do i = 1,12
        if (Active(i)) then 
            call slipsys(slip,Schmid,m,n,i)
            Dpc = Dpc+tautr(i)/abs(tautr(i))*x(countera)*1/2*(Schmid+transpose(schmid))/dt0/nlines
        countera = countera+1
        end if 
    end do
    Dpc =  matmul(transpose(Rtest),matmul(Dpc,Rtest))
    Dp = Dp + matmul(transpose(R(1:3,1:3,j)),matmul(Dpc,R(1:3,1:3,j)))
    
    gammatoti = gammatoti + sum(x)/nlines
    

    
    F0int(1:3,1:3,j) = F1
    Fp0int(1:3,1:3,j) = Fp1
    s0in(j,1:12) = s1
    Tagcint(1:3,1:3,j) = matmul(transpose(R(1:3,1:3,j)),matmul(T1,R(1:3,1:3,j)))
    Tag = Tag+Tagcint(1:3,1:3,j)/nlines
    
    !write(7,*) 'cauchy stress'
    
    !write(7,*) 'updated resolved shear stress'
    !write(7,*) tau
    !write(7,*) 'deltagamma'
    !write(7,*) x
    !write(7,*) PA
    !write(7,*)
    
    
    
end do grains 

end subroutine timestep









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

subroutine slipsys(slip,Schmid,m,n,snum)

!Declear all returning variables
real(8), dimension(3,3),intent(out) :: Schmid
real(8), dimension(6,12), intent(in) :: slip
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

subroutine hparam(alpha, beta, s0, h,slip) !subroutine for determination of hardening moduli
    !Declear varialbles 
        integer, intent(in) :: alpha, beta
        real(8), dimension(12), intent(in) :: s0
        real(8) , intent(out) :: h
        real(8) :: q, h0, s, a, hbeta
        real(8) , dimension(6,12) :: slip
    !Parameters
        
        h0 = 180
        s = 148
        a = 2.25
    !Calculation
        if (DOT_PRODUCT(slip(1:3,alpha),slip(1:3,beta)) == 1) then 
            q = 1
        else 
            q = 1.4
        end if 

        hbeta = h0*(1-s0(beta)/s)**a
        if (alpha == beta ) then
            h = (q+(1-q))*hbeta
        else 
            h = q*hbeta  
        end if 
return
end subroutine hparam



subroutine sincr(PA,x,tautr,s0,Ctr)    !Calculates the slip increments
    use global
    implicit none
    
    real(8) ::  coeffA, switch, sgn,h 
    real(8) , dimension(3,3) :: Ctr
    real(8) , dimension(12) :: tautr, s0
    logical, dimension(12) :: PA
    double precision, dimension(12,12) :: A, U, VT, Eps, Aplus, Atest, Adisp, Eps1
    integer :: i,j,countera,LWMAX , sizeA,counterb, INFO,lda, sw
    double precision, dimension(12):: x,b,Sing,test
    double precision, dimension(:), allocatable:: WORK
    integer :: minpos,teller

 lda = 12
 sw = 0
    LWMAX = 1000
Allocate(WORK(LWMAX))
switch = -1
teller = 1
! The process must be repeated until all calculated slip increments are positive
outer: do
    A = 0
    U = 0
    VT = 0
    Aplus = 0
    x = 0
    b =0
    Sing = 0
    Eps = 0
    test =0
    Atest =0
    Adisp = 0 
    Eps1 = 0
  

!Allocate size of matrices used for SVD 
sizeA = COUNT(PA)



countera = 1   ! Counter used to index the position of each element of A
do i = 1,12
    if (PA(i) ) then
        b(countera) = abs(tautr(i))-s0(i)
        counterb = 1
        do j = 1,12
            if (PA(j) )then
                call hparam(i,j,s0,h,slip)
                !write(7,*) h
                call Acoeff(i,j,Ctr,coeffA)
               
                sgn = tautr(i)*tautr(j)/abs(tautr(i)*tautr(j))
                
                A(countera,counterb) = sgn*coeffA+h
                counterb = counterb+1
            end if 
        end do
        countera = countera+1
    end if 

end do

Adisp = A

! Compute singular value decomposition of A

!call dgesvd('A','A',sizeA,sizeA,A(1:sizeA,1:sizeA),sizeA,Sing(1:sizeA),U(1:sizeA,1:sizeA),sizeA,VT(1:sizeA,1:sizeA),sizeA, &
!             WORK, LWMAX ,INFO)
call dgesvd('A','A',sizeA,sizeA,Adisp,lda,Sing,U,lda,VT,lda, &
            WORK, LWMAX ,INFO)

!Calculate sigma+  
do i = 1,12   
if (Sing(i) > 0.000000001) then
    eps(i,i) = 1/Sing(i)
    Eps1(i,i) = Sing(i)
Else
    Eps(i,i) = 0
    Eps1(i,i) = 0
end if
end do
!write(7,*) 'Sing'
!write(7,*) Sing
!write(7,*) INFO


Aplus = matmul(transpose(VT),matmul(Eps,transpose(U)))   !pseudoinverse

Atest = matmul(U,matmul(Eps1,VT))


x(1:sizeA) = matmul(Aplus(1:sizeA,1:sizeA),b(1:sizeA))

test(1:sizeA) = matmul(A(1:sizeA,1:sizeA),x(1:sizeA))

!write(7,*) 'x'
!write(7,*) x


!Loop through the calculated deltagammas in order to remove the negative ones from the logical array PA
minpos = minloc(x, DIM=1)
countera = 1
do i = 1,12
    if (countera == sizeA+1) then
        exit
    else if (PA(i) .and. x(countera) <= 0 .and. countera == minpos) then
        PA(i) = .FALSE.
        countera = countera+1
    else if (PA(i)) then
        countera =countera+1
    end if
end do



!Criterion to determine final set of active slip systems, if all deltagamma > 0 the iteration is ended
if (count(x(1:sizeA) <= 0 ) == 0) then
    switch = 1
   
    exit outer
end if
end do outer

return
end subroutine sincr




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
    real(8) :: A, B, C, F, G, H, dI2,dI3, dtheta, I2, I3, theta, m = 8.8 , dsum,sum, sFi
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
    theta = acos(I3/I2**(3./2.))
    sum = ((2*cos((2*theta+pi)/6.))**(m))+((2*cos((2*theta -  3*pi)/6.))**(m)) + ((-2*cos((2*theta +  5*pi)/6.))**(m))
    sFi = (3*I2)**(m/2.)*sum
    
    
    do i = 1,6
        dI2 = 2./3.*(F*partials(i,4)+G*partials(i,5)+H*partials(i,6)) &
        +2./54.*((A-C)*(partials(i,1)-partials(i,3))) &
        +2./54.*((C-B)*(partials(i,3)-partials(i,2))) &
        +2./54.*((B-A)*(partials(i,2)-partials(i,1)))
        

        dI3 = 1./54.*((partials(i,3)-partials(i,2))*(A-C)*(B-A) & 
        + ((partials(i,1)-partials(i,3))*(B-A)+(A-C)*(partials(i,2)-partials(i,1)))*(C-B)) &
        + partials(i,4)*G*H + F*partials(i,5)*H + F*H*partials(i,6) &
        - (2*F/6*partials(i,4)*(C-B)+F**2/6*(partials(i,3)-partials(i,2))) &
        - (2*G/6*partials(i,5)*(A-C)+G**2/6*(partials(i,1)-partials(i,3))) &
        - (2*H/6*partials(i,6)*(B-A)+H**2/6*(partials(i,2)-partials(i,1))) 

        dtheta = -1/sqrt(1-I3**2/I2**3)*(dI3*I2**(3./2.)-3./2.*I2**(1./2.)*dI2*I3)/I2**3
        
        dsum = m*dtheta* &
        (  ( 2*cos((2*theta +   pi)/6.))**(m-1) * (-4./6.*sin((2*theta +   pi)/6)) &
         + ( 2*cos((2*theta - 3*pi)/6.))**(m-1) * (-4./6.*sin((2*theta - 3*pi)/6.)) &
         + (-2*cos((2*theta + 5*pi)/6.))**(m-1) * ( 4./6.*sin((2*theta + 5*pi)/6)))
    
        n(i) = m/2*3**(m/2.)*I2**(m/2.-1.)*dI2*sum+(3*I2)**(m/2.)*dsum


    end do
       !n = n/norm2(n)
!write(*,*) n/norm2(n)
!write(8,*) tag(1,1), tag(2,2), n(1), n(2)
call vec2tens(grad,n)



grad = grad*1/(2*m)*(sFi/2)**(1/m-1)
!write(*,*) grad

return
end subroutine hoshfordnormal

!! Subroutine for calculation of the isotropic elastic constant
subroutine Elasticconstant(Chook, mu)
 use global
    implicit none
    real(8), dimension(3,3) :: La,Dp,tag,F,E
    real(8) :: gammatot,gammatoti, dt0, mu , lambda
    real(8) , dimension(3,3,nlines)  :: Fp0,F0,Fp0i,F0i
    real(8),  dimension(nlines,12) :: S0,S0i
    real(8) , dimension(6,6)              :: Chook
    integer :: n, i
    
            La(1,1) = 1.
            La(2,2) = -1./2.
            La(3,3) =-1./2.
            La(1,2) =0
            La(2,1) = 0
            La(1,3) = 0
            La(3,1) = 0
            La(2,3) = 0
            La(3,2) = 0
            
            do i = 1,nlines
            F0(1:3,1:3,i) = id
            Fp0(1:3,1:3,i) = id
            end do
            s0(:,1:12) = 10
            gammatot = 0
     
           ! do while (gammatot == 0)    
           !     do i = 1,nlines
           !         F0(1:3,1:3,i) = id
           !         Fp0(1:3,1:3,i) = id
           !         end do
           !     dt0 = dt *n 
           !     gammatot = 0
           !     s0(:,1:12) = 10
           !     call timestep(Tag, Dp, La, gammatot, gammatoti , Fp0, Fp0i, F0, F0i,S0i,s0,dt0)
           !     gammatot = gammatoti
           !      n = n + 1
           ! end do
            gammatoti = 0   
            n = 1
            dt0 = dt *n 
            call timestep(Tag, Dp, La, gammatot, gammatoti , Fp0, Fp0i, F0, F0i,S0i,s0,dt0)
            F = id + La*dt0
            E = 1./2.*(matmul(transpose(F),F)-id)
            mu =  1./2.*(tag(3,3)-tag(1,1))/(E(3,3)-E(1,1))
            lambda = (tag(1,1)-2*mu*E(1,1))/(E(1,1)+E(2,2)+E(3,3))
           
           
            !Elasticity tensor
            Chook = 0
            forall (i = 1:3) Chook(i,i)= 2*mu +lambda
            forall (i = 4:6) Chook(i,i)= lambda
            forall (i = 1:2) Chook(i+1,i) = lambda
            forall (i = 1:2)  Chook(i,i+1) = lambda
            Chook(3,1) = lambda
            Chook(1,3) = lambda
            
end subroutine

subroutine Yoshidamodel(Tag,D,Dp)
    use global
    implicit none

    real(8) , dimension(3,3) :: Dp,Dpt,Dpn,Tag,D,N, DdevT, tensprod,Ddev,Nnorm
    integer :: i,j,k,l
    real(8) :: kronecker,lambdadot,lambdadottemp, theta,alpha, G, theta0, c1 = 0.3
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

subroutine elastoplasticmoduli(Cep,coeff,tag,Active)
    use global
    implicit none


    real(8), dimension(3,3,3,3) :: Cep,Cel4, dyadic
    real(8), dimension(3,3) :: P,W,tag, CP, schmid_a,schmid_b
    integer :: i,j,k,l,m,n,coeff
    logical, dimension(12) :: Active


    do m = 1,12
    
    !if (Active(m)) then   
    !call slipsys(Schmid_a,m)
    P = 1/2*(Schmid_a+transpose(schmid_a))
    W = 1/2*(Schmid_a-transpose(schmid_a))
    
    do i = 1,3
        do j = 1,3
        end do
        end do
    !    end if 
    end do
    end subroutine elastoplasticmoduli