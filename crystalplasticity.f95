module crystalplasticity
use global

contains
subroutine newton(k,part,bryter,F0i,Fp0i,S0i,pw)
    
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
    
        implicit none
    
     
    !Declear all variables
 
    real(8) :: phi1, Phi, phi2,dt0, sigmaeq, gammatot,pw,gammatoti,gammaskrank,dot, dl, nor1,nor2,epsp
    real(8) , dimension(3,3)  :: grad,T0 , Dp, Lb, Tagb, La0,Dp2
                        
    real(8) , dimension(3,3) :: La
    real(8) , dimension(3,3), intent(out)  :: Tag
    real(8) , dimension(3,3,nlines) ::  F0, Fp0,Fp0i,Fp0int,F0i,F0int, Tagc, Tagcint
    real(8),  dimension(nlines,12) :: s0,S0i,S0in
    integer :: i, switch , o,p,k,h, bryter,secit
    real(8) , dimension(4,4)  :: Jacob, Jinv
    real(8) , dimension(4)  :: sigma, offdl
    integer :: LDA = 4,NRHS = 1, Info,  minl, maxl,nit,bcond=2
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
        epsp = 0
    !Copies of the input variables, in order not to update the initial condition when calculating instantaneous yield surface.
        S0 = s0i 
        Fp0 = Fp0i  
        F0 = F0i  

        propconst = (/0.0, 0.0, 0.0, 0.0, 0.0/)
    
    
 
   
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
if (bcond == 2) then
    La(1,1) = 1/sqrt(1.0+propconst(1)**2)
    La(2,2) = propconst(1)**2/sqrt(1.0+propconst(1)**2)
    La(1,2) = 0
    La(2,1) = 0
    La(1,3) = 0
    La(3,1) = 0
    La(2,3) = 0
    La(3,2) = 0
    La(3,3) = -1.0/3.0*(La(1,1)+La(2,2))

end if   

    gammaskrank = 0.000000001
    

 
    secit = 0
    dt0 = dt
    switch = 1
    if (bryter == 3) then
        goto 13
    end if


    iter: do    
    
       
   
   
Select case (bcond)
case (1)    
!!! Iteration to ensure boundary condition is satisfied at througout all timesteps. 
    nit = 0   
    boundarycond: do  while (nit < 10)  
       ! write(*,*) La
    nor1 = norm2(La)
       call timestep(Tag, Dp, La, gammatot, gammatoti, Fp0, Fp0int, F0, F0int,S0in,s0,dt0)
       !write(*,*) sigma , Tag(1,1), Tag(2,2), gammatoti
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
           !write(*,*) sigma2 , Tag(1,1), Tag(2,2), gammatoti
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
    !call eqvstr(Tag,sigmaeq)
    !write(10,*) sigmaeq

    if (bryter == 1 .or. bryter == 5 .or. bryter == 4) then
       
        
       
        if (gammatoti > pw .and. abs((gammatoti - pw)/pw) > pwpercision) then
            if (pw == 0) then
                dt0 = dt0/2
            else
            dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot)*0.5
            end if
            !write(*,*) dt0, sigma, gammatoti-pw, nor1, nor2
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
        epsp = epsp+norm2(Dp)*sqrt(2./3.)*dt
        if (bryter == 5) then
            if (gammatot > gammaskrank) then
                write(8,*) bryter, gammatoti
                call contract2(La,Dp,dot)
                dot = dot/norm2(La)/norm2(Dp)
            write(11,*) Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3) , dot , acos(dot), gammatot
            
            gammaskrank = gammaskrank + 0.0000001
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
            
            write(3,*) tag(1,1), epsp
            !
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
          
            gammaskrank = gammaskrank + 0.00000005
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
            write(3,*) tag(1,1) , gammatot
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
    call slipsys(Schmid,m,n,i)
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
            call slipsys(Schmid,m,n,i)
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
            call slipsys(Schmid,m,n,i)
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
    ! The logical parameter hardening is set in global module to determine if a hardening material is to be considered.
    if (hardening .eqv. .true.) then
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
    end if
    !Step 10, check consistency
    
    !Calculates resolved shear stress using updated stress in addition ti the vector for conistency and the final active slip systems
    do i = 1,12
        call slipsys(Schmid,m,n,i)
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
            call slipsys(Schmid,m,n,i)
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
    
end do grains 

end subroutine timestep









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
                !!! The logical hardening is set in the global module if a hardening material is to be considered.
                if (hardening .eqv. .true.) then
                call hparam(i,j,s0,h,slip)
                else
                h = 0
                end if    
                
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


!! Subroutine for calculation of the isotropic elastic constant
subroutine elasticconstant(Chook, mu)
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
               
end subroutine Elasticconstant


   subroutine pseudoinv(PA,tautr,s0,Ctr,Aplus)    !Calculates the slip increments
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
                !!! The logical hardening is set in the global module if a hardening material is to be considered.
                if (hardening .eqv. .true.) then
                call hparam(i,j,s0,h,slip)
                else
                h = 0
                end if    
                
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
end subroutine pseudoinv


subroutine elastoplasticmoduli(Cep,tag,Active,tautr,s0,Ctr)
    implicit none


    real(8), dimension(3,3,3,3) :: Cep,Cel4, dyadic
    real(8), dimension(3,3) :: P,W,tag, CP, schmid_a,schmid_b,Ctr,Leftdy,PC
    integer :: i,j,k,l,t,q,coeff,countera,counterb
    logical, dimension(12) :: Active
    real(8), dimension(3) :: m,n
    real(8), dimension(12,12) :: Aplus
    real(8), dimension(12) :: s0, tautr

    do i = 1,3
        Cel4(i,i,i,i) = Cel(1,1)
        do j = 1,3
                   Cel4(i,i,j,j) = Cel(1,2)
                   Cel4(i,j,i,j) = Cel(6,6)
        end do
    end do

    Do i = 1,3
    do j = 1,3
        do l = 1,3
            do k = 1,3
                CP(i,j) = CP(i,j)+Cel4(i,j,l,k)*P(k,l)
            end do
        end do
    end do
end do
        


    call pseudoinv(Active, tautr,s0,Ctr, Aplus)
    countera = 1 
    do t = 1,12

    
    if (Active(t) .eqv. .true.) then  
        
            call slipsys(Schmid_a,m,n,t)
            P = 1./2.*(Schmid_a+transpose(schmid_a))
            W = 1./2.*(Schmid_a-transpose(schmid_a)) 
            counterb = 1
            
            
            ! calculate second order tensor (Ce:P+W_a*Sigma-Sigma*W_a)
            leftdy = 0

            Do i = 1,3
            do j = 1,3
                do l = 1,3
                    do k = 1,3
                        CP(i,j) = CP(i,j)+Cel4(i,j,l,k)*P(k,l)
                    end do
                end do
            end do
             end do
             write(*,*) CP
             CP = 0
            call voigt(P,CP) 
            write(*,*) CP
            do i = 1,3
                do j = 1,3
                   do k = 1,3
                     leftdy = leftdy(i,j) + W(i,k)*Tag(k,j)- Tag(i,k)*W(k,j)
                   end do
                end do
            end do
            
            do q = 1,12
            if (Active(q)) then

                call slipsys(Schmid_a,m,n,q)
              P = 1./2.*(Schmid_a+transpose(schmid_a))
            ! Calculate right hand side of dyadic product (P:Ce)
            
              do i = 1,3
                do j = 1,3

                end do
            end do
            end if
            end do
    
        
       end if 
    end do
    end subroutine elastoplasticmoduli

end module 