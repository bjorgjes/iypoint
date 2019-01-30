program ypoint
implicit none
integer :: k, part,nlines,i,bryter,chunck
real(8) :: t1,t2,omp_get_wtime,pw1,pw2
real(8) , dimension(:,:), Allocatable ::eul
real(8) , dimension(:,:,:), Allocatable  :: F0i, Fp0i,F0,Fp0
real(8),  dimension(:,:), Allocatable :: s0i,S0
t1 = omp_get_wtime()
open(unit=11,file="result2.txt",status='replace')
open(unit=3,file="Stresss.txt",status='replace')
open(unit=8,file="eulerangles.txt",status='replace')
open(unit=13,file="Dp.txt",status='replace')
call countlin('euleranglesin.txt',nlines)
write(*,*) nlines

open(action='read',unit=15,file="euleranglesin.txt",status='old')

allocate(eul(nlines,3))
do i = 1,nlines
    read(15,*) eul(i,1:3)
    write(*,*) eul(i,1:3)
end do
write(*,*)

Allocate(F0i(3,3,nlines),Fp0i(3,3,nlines))
Allocate(s0i(nlines,12))

pw1 = 0.001
bryter = 5
!call newton(1,5,nlines,eul,bryter,F0i,Fp0i,S0i,pw1)   
write(*,*) 'check1'
!F0 = F0i
!S0 = S0i
!Fp0 = Fp0i
pw2 = 0.01
part = 200
call OMP_SET_NUM_THREADS(7)
!$OMP PARALLEL PRIVATE(k,F0i,S0i,Fp0i,bryter)
!$OMP DO
do k = 0,2*part
    bryter = 5
    call newton(k,part,nlines,eul,bryter,F0i,Fp0i,S0i,pw2)
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



pw2 = 0.000001
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



subroutine newton(k,part,nlines,eul,bryter,F0i,Fp0i,S0i,pw)
implicit none
real(8) , dimension(3,3)  :: La, Tag, Dp
real(8) , dimension(4,4)  :: Jacob, Jinv
real(8) , dimension(4)  :: sigma, offdl
real(8) :: epsilon, jacobidel, l, pi = 3.14865,dt,pw
integer :: i, j , LDA = 4,NRHS = 1, Info, k, minl, maxl, part,nlines,omp_get_thread_num,teller,bry
integer  :: bryter
integer , dimension(4) ::pos1, pos2, IPIV
real(8), dimension(nlines,3), intent(in) :: eul
real(8) , dimension(3,3,nlines)  :: Fp0i,F0i
real(8),  dimension(nlines,12) :: S0i

bry = 0
if (bryter == 4) then 
    bry = 1
else if (bryter == 3) then
    bryter = 1
    bry = 1
else if (bryter == 5) then
    bryter = 1
    bry = 2
else if (bryter == 6) then
    bryter = 4
    bry = 5
end if



offdl = 0
epsilon = 0.00001
pos1 = (/1, 1, 2, 3/)
pos2 =(/2, 3, 3, 3/)

dt = 0.00001

10 CONTINUE
teller = 0

!Initialize La
l =1 !strainrate
!do h = (k-1)*chunck,k*chunck-1
!write(*,*) cos(pi*k/part), sin(pi*k/part)
!write(*,*)
La = 0 
La(1,1) = cos(pi*k/part)
La(2,2) = sin(pi*k/part)
La(3,3) =-0.3*(La(1,1)+La(2,2))
La(1,2) = offdl(1)
La(2,1) = offdl(1)
La(1,3) = offdl(2)
La(3,1) = offdl(2)
La(2,3) = offdl(3)
La(3,2) = offdl(3)

    call taylor(La,Tag,nlines,eul,bryter,F0i,Fp0i,S0i,dt,pw,Dp)

    if (bry == 1) then ! If one want relaxed state. 
        !Performes relaxation
        call taylor(La,Tag,nlines,eul,3,F0i,Fp0i,S0i,dt,pw,Dp)
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

subroutine centraldiff(La,jacobidel,epsilon, pos1l, pos2l,pos1T,pos2T,nlines,eul,bryter,F0i,Fp0i,s0i,dt,pw,Dp)
    implicit none
    real(8) , dimension(3,3)  :: La, Tag1, Tag2, Laint1, Laint2, Dp
    real(8), intent(in) :: epsilon , dt,pw
    real(8), intent(out) :: Jacobidel
    integer, intent(in) :: pos1l,pos2l,pos1T,pos2T,nlines,bryter
    real(8), dimension(nlines,3), intent(in) :: eul
    real(8) , dimension(3,3,nlines)  :: Fp0i,F0i
    real(8),  dimension(nlines,12) :: S0i

    Laint1 = La
    Laint2 = La
    if (pos1l /= pos2l) then
    Laint1(pos1l,pos2l) = Laint1(pos1l,pos2l)+epsilon
    Laint2(pos1l,pos2l) = Laint2(pos1l,pos2l)-epsilon
    Laint1(pos2l,pos1l) = Laint1(pos1l,pos2l)+epsilon
    Laint2(pos2l,pos1l) = Laint2(pos1l,pos2l)-epsilon
    else if (pos1l == pos2l) then 
        Laint1(pos1l,pos2l) = Laint1(pos1l,pos2l)+epsilon 
        Laint2(pos1l,pos2l) = Laint2(pos1l,pos2l)-epsilon
    end if      
    Call taylor(Laint1,Tag1,nlines,eul,bryter,F0i,Fp0i,S0i,dt,pw,Dp)
    call taylor(Laint2,Tag2,nlines,eul,bryter,F0i,Fp0i,S0i,dt,pw,Dp)

    jacobidel = (Tag1(pos1T,pos2T)-Tag2(pos1T,pos2T))/(2*epsilon)
    return
end subroutine centraldiff






subroutine taylor(La,Tag,nlines,eul,bryter,F0i,Fp0i,S0i,dt,pw,Dp)
    !use setparam
        implicit none
    
     
    !Declear all variables
 
    real(8) :: phi1, Phi, phi2, det, dt,dt0, sigmaeq, cons, gammatot,pw,gammatoti,gammaskrank,dot, dl
    real(8) , dimension(3,3)  :: F1, Fp1, Fp1inv, Fe1, Fetr,Fp0inv, Ctr,T0,T1,Ttr,Etr , & 
                        id, Schmid,TS,Tint, CS, Tst, Fpint ,Rn,  Rtest,Dpc,Dp, Lb,Tagb,La0
    real(8) , dimension(3,3) :: La
    real(8) , dimension(3,3), intent(out)  :: Tag
    real(8) , dimension(3,3,nlines) :: R, F0, Fp0,Fp0i,Fp0int,F0i,F0int, Tagc, Tagcint, Lc
    real(8),  dimension(nlines,12) :: s0,S0i,S0in
    real(8) , dimension(6,12) :: slip
    real(8) , dimension(6,6) :: Cel
    real(8) , dimension(12) :: tautr,s1,tau, consis, s0int
    real(8) , dimension(3) :: m,n
    logical, dimension(12) :: PA, Active
    logical :: consise
    integer :: i,countera,ij, switch , numact, o,p,k,h, bryter,bakit,nlines,j,e,secit
    double precision, dimension(12):: x
    real(8), dimension(nlines,3) , intent(in) :: eul
    real(8) , dimension(4,4)  :: Jacob, Jinv
    real(8) , dimension(4)  :: sigma, offdl
    integer :: LDA = 4,NRHS = 1, Info,  minl, maxl,nit
    integer , dimension(4) ::pos1, pos2, IPIV
    !call countlin('euleranglesin2.txt',nlines)
    !open(unit=7,file="result.txt",status='replace')  ! outputfile used for testing. 
    !open(unit=8,file="eulerangles.txt",status='replace')
    !open(unit=10,file="Stress.txt",status='replace')
    !open(unit=9,file="euleranglesin2.txt",status='old')

    !Timeincrement
    !dt0 = 0.0000001

    pos1 = (/1, 1, 2, 3/)
    pos2 =(/2, 3, 3, 3/)
    dl = 0.000001
    La0 = La
    !Define velocity gradient
    !strainrate
    gammatot = 0
!Copies of the input variables, in order not to update the initial condition when calculating instantaneous yield surface.
  S0 = s0i 
  Fp0 = Fp0i  
  F0 = F0i  
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
    
     bakit = 0
    !dt = dt0 
   
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
    
    iter: do while (switch < 5000000)     
    if (bryter == 3) then
        goto 13
    end if
       
   
   
   !!! Iteration to ensure boundary condition is satisfied at througout all timesteps. 
    nit = 0   
    boundarycond: do  while (nit < 10)  
       ! write(*,*) La
       call timestep(Tag, Dp, La, gammatot, gammatoti , nlines, Fp0, Fp0int, F0, F0int,R,S0in,s0,dt0, slip , Cel,x)
      ! write(*,*) Tag
       do h = 1,4
        sigma(h) = Tag(pos1(h),pos2(h))
       end do
      ! write(*,*) sigma
         minl = minloc(sigma, DIM = 1)
         maxl = maxloc(sigma, DIM = 1)
         if (abs(sigma(minl))<0.00001 .and. abs(sigma(maxl)) < 0.00001) then
            
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

                call timestep(Tagb, Dp, Lb, gammatot, gammatoti , nlines, Fp0, Fp0int, F0, F0int,R,S0in,s0,dt0, slip , Cel,x)
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
        
        !write(*,*)  
        !write(*,*) jacob(1,1:4)
        !write(*,*) jacob(2,1:4)
        !write(*,*) jacob(3,1:4)
        !write(*,*) jacob(4,1:4)
        !write(*,*) 
        minl = minloc(sigma, DIM = 1)
        maxl = maxloc(sigma, DIM = 1)
       
       
       
        if (abs(sigma(minl))<0.01 .and. abs(sigma(maxl)) < 0.01) then
           ! write(*,*) tag
            exit boundarycond
            
        end if 

   nit = nit+1
   !write(*,*) nit
    end do boundarycond
   
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
       
        
       
        if (gammatoti > pw .and. abs((gammatoti - pw)/pw) > 0.000000001) then
            if (pw == 0) then
                dt0 = dt0/2
            else
            dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot) 
            end if
            
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
                call contract2(La,Dp,dot)
                dot = dot/norm2(La)/norm2(Dp)
            write(11,*) Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3) , dot , acos(dot), gammatot
            
            gammaskrank = gammaskrank + 0.0001
            end if
        end if 

        if (abs((gammatot -pw)/pw) <= 0.000000001) then
            Fp0i = Fp0
            F0i  = F0
            S0i = S0   
            exit iter
        end if
    
    else if (bryter == 2 .or. bryter == 6) then
        
        if (pw /= 0) then
        if (gammatoti > pw .and. abs((gammatoti - pw)/pw) > 0.000000001) then
            if (pw == 0) then
                dt0 = dt0/2
            else
            dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot) 
            end if
            switch = switch +1
            secit = secit +1
            if (secit > 15) then 
                exit iter
            end if 
            cycle iter
        end if 
        end if

        if (pw == 0) then
            if (gammatoti > pw .and. abs(gammatoti) > 0.000000001) then
            dt0 = dt0/2
            !dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot) 
            switch = switch +1
            secit = secit +1
            if (secit > 30) then 
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
                call contract2(La,Dp,dot)
                dot = dot/norm2(La)/norm2(Dp)
            write(11,*) Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3) , dot , acos(dot), gammatot
            gammaskrank = gammaskrank + 0.00001
            end if
        end if 
        
        if (pw /= 0 .and. abs((gammatot -pw)/pw)<= 0.000000001) then
            exit iter
        else if (pw == 0 .and. gammatot <= 0.000000001 .and. gammatot > 0) then
           ! write(*,*) 'check'
            exit iter

        end if
    else if (bryter == 3) then
        13 continue
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
        write(*,*) La
        do i = 1,3
            
            call timestep(Tag, Dp, -La, gammatot, gammatoti , nlines, Fp0, Fp0int, F0, F0int,R,S0in,s0,dt0, slip , Cel,x)
            write(*,*) tag
            write(*,*) x
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
   
        exit iter
    bakit = bakit+1  
    write(*,*) bakit 
    write(*,*) x
    if (bakit == 10) then
        !write(*,*) gammatot
        write(*,*) bakit
        exit iter
    end if 
!else if (bryter == 4) then 
!    if (gammatoti > pw .and. abs((gammatoti - pw)/pw) > 0.000000001) then
!        if (pw == 0) then
!            dt0 = dt0/2
!        else
!        dt0 = (pw-gammatot)*dt0/(gammatoti-gammatot) 
!        end if
!        switch = switch +1
!        secit = secit +1
!        if (secit > 15) then 
!            exit iter
!        end if 
!        cycle iter
!    end if 
    
    
!    F0 = F0int
!    Fp0 = Fp0int
!    Tagc = Tagcint
!    s0 = s0in
!    gammatot = gammatoti 
!    
!    if (abs((gammatot -pw)/pw)<= 0.000000001) then
!        Fp0i = Fp0
!        F0i  = F0
!        S0i = S0   
!       
 !       exit iter
!
!    end if

    end if


    switch = switch +1
    end do iter
    !write(7,*) PA
    !write(7,*) phi1, Phi, phi2 
    !end do !l 
    !end do !p
    !close(9,STATUS='keep')
    !close(7,status='keep')
    !close(8,status='keep')
    return
    end subroutine taylor





subroutine timestep(Tag, Dp, La, gammatot, gammatoti , nlines, Fp0, Fp0int, F0, F0int,R,S0in,s0,dt0, slip,Cel,x )
    
      !use setparam
    implicit none
    
     
    !Declear all variables
 
    real(8) :: phi1, Phi, phi2, det,dt0, cons, gammatot,gammatoti
    real(8) , dimension(3,3)  :: F1, Fp1, Fp1inv, Fe1, Fetr,Fp0inv, Ctr,T1,Ttr,Etr , & 
                        id, Schmid,TS,Tint, CS, Tst, Fpint ,Rn,  Rtest,Dpc,Dp
    real(8) , dimension(3,3), intent(in)  :: La
    real(8) , dimension(3,3), intent(out)  :: Tag
    real(8) , dimension(3,3,nlines) :: R, F0, Fp0,Fp0int,F0int, Tagcint, Lc
    real(8),  dimension(nlines,12) :: s0,S0in
    real(8) , dimension(6,12) :: slip
    real(8) , dimension(6,6) :: Cel
    real(8) , dimension(12) :: tautr,s1,tau, consis, s0int
    real(8) , dimension(3) :: m,n
    logical, dimension(12) :: PA, Active
    logical :: consise
    integer :: i,countera,ij, switch , numact, o,nlines,j,secit
    double precision, dimension(12):: x
 
    !Create identity matrix
    id = 0
    forall(i = 1:3) id(i,i) = 1
    
    
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
    call sincr(Active,x,tautr,s0int,Ctr,Cel,slip,phi1,Phi,phi2) 
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
    ! do i = 1,12
    !    do j = 1,12
    !        call hparam(i,j,s0,h,slip)
    !        s1(i) = s1(i) + h*x(j)
    !    end do
    !end do
    
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
real(8) :: l = 2.0, k = 3.0
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

subroutine deter(A,det)  ! Subroutine to calculate determinant of a 3x3 matrix. 
!Declear variables
    real(8), dimension(3,3),intent(in) :: A
    real(8), intent(out) :: det

    det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
return
end subroutine deter


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

subroutine Acoeff(alpha,beta, Ctr, coeff,Cel,slip)
!Declear variables
    integer, intent(in) :: alpha, beta
    real(8), dimension(3,3),intent(in) :: Ctr
    real(8), dimension(3,3) :: Salpha, Sbeta, rm, csv, tot
    real(8), dimension(3) :: m,n
    real(8), dimension(6,6),intent(in) :: Cel
    real(8), dimension(6,12),intent(in) :: slip
    real(8), dimension(6) :: CS
    real(8), intent(out) :: coeff
!Computation
    
    Call slipsys(slip,Salpha,m,n,alpha)
    Call slipsys(slip,Sbeta,m,n,beta)
    csv = matmul(Ctr,Sbeta)
    csv = (csv+transpose(csv))/2
    CS = (/csv(1,1),csv(2,2),csv(3,3),csv(2,3),csv(1,3),csv(1,2)/)
    rm(1,1) = DOT_PRODUCT(Cel(1,1:6),CS)
    rm(2,2) = DOT_PRODUCT(Cel(2,1:6),CS)
    rm(3,3) = DOT_PRODUCT(Cel(3,1:6),CS)
    rm(1,2) = DOT_PRODUCT(Cel(6,1:6),CS)
    rm(2,1) = DOT_PRODUCT(Cel(6,1:6),CS)
    rm(1,3) = DOT_PRODUCT(Cel(5,1:6),CS)
    rm(3,1) = DOT_PRODUCT(Cel(5,1:6),CS)
    rm(3,2) = DOT_PRODUCT(Cel(4,1:6),CS)
    rm(2,3) = DOT_PRODUCT(Cel(4,1:6),CS)
tot = matmul(transpose(Salpha),rm)
coeff = tot(1,1)+tot(2,2)+tot(3,3)
return
end subroutine Acoeff

subroutine sincr(PA,x,tautr,s0,Ctr,Cel,slip,phi1,Phi,phi2)    !Calculates the slip increments
    
    
    real(8) :: h, coeffA, switch, sgn,phi1,Phi,phi2
    real(8) , dimension(3,3) :: Ctr
    real(8) , dimension(6,12) :: slip
    real(8) , dimension(6,6) :: Cel
    real(8) , dimension(12) :: tautr, s0
    logical, dimension(12) :: PA
    double precision, dimension(12,12) :: A, U, VT, Eps, Aplus, Atest, Adisp, Eps1
    integer :: i,j,countera,LWMAX , sizeA,counterb, INFO,lda, sw
    double precision, dimension(12):: x,b,Sing,test
    double precision, dimension(:), allocatable:: WORK
    integer :: minpos

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
                !call hparam(i,j,s0,h,slip)
                !write(7,*) h
                call Acoeff(i,j,Ctr,coeffA,Cel,slip)
               
                sgn = tautr(i)*tautr(j)/abs(tautr(i)*tautr(j))
                
                A(countera,counterb) = sgn*coeffA
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

subroutine voigt(A,B)
    real(8) , dimension(6,6) :: Cel
    real(8) , dimension(6) :: vec
    real(8) , dimension(3,3) :: A,B
    
    !Elasticity tensor
    Cel = 0
    forall (i = 1:3) Cel(i,i)= 170*1000
    forall (i = 4:6) Cel(i,i)= 75*1000
    forall (i = 1:2) Cel(i+1,i) = 124*1000
    forall (i = 1:2)  Cel(i,i+1) = 124*1000
    Cel(3,1) = 124*1000
    Cel(1,3) = 124*1000

    vec = (/A(1,1),A(2,2),A(3,3),A(2,3),A(1,3),A(1,2)/)

    B(1,1) = DOT_PRODUCT(Cel(1,1:6),vec)
    B(2,2) = DOT_PRODUCT(Cel(2,1:6),vec)
    B(3,3) = DOT_PRODUCT(Cel(3,1:6),vec)
    B(1,2) = DOT_PRODUCT(Cel(6,1:6),vec)
    B(2,1) = DOT_PRODUCT(Cel(6,1:6),vec)
    B(1,3) = DOT_PRODUCT(Cel(5,1:6),vec)
    B(3,1) = DOT_PRODUCT(Cel(5,1:6),vec)
    B(3,2) = DOT_PRODUCT(Cel(4,1:6),vec)
    B(2,3) = DOT_PRODUCT(Cel(4,1:6),vec)

end subroutine voigt

subroutine textur(Fe1,slip,nslip)
    
    real(8), dimension(6,12) :: slip, nslip
    real(8), dimension(3,3) :: Fe1, Feinv
    integer :: i
    
call minv3(Fe1,Feinv)

    do i = 1,12
        nslip(1:3,i) = matmul(transpose(Fe1),slip(1:3,i))
        nslip(1:3,i) = nslip(1:3,i)/SQRT(nslip(1,i)**2+nslip(2,i)**2+nslip(3,i)**2)

        nslip(4:6,i) = matmul(transpose(Feinv),slip(4:6,i))
        nslip(4:6,i) = nslip(4:6,i)/SQRT(nslip(4,i)**2+nslip(5,i)**2+nslip(6,i)**2)
    end do
    return 
    
end subroutine textur

subroutine textr(R,slip,nslip)
    
    real(8), dimension(6,12) :: slip, nslip
    real(8), dimension(3,3) :: R, Rinv
    integer :: i
    
    call minv3(R,Rinv)

    do i = 1,12
        nslip(1:3,i) = matmul(R,slip(1:3,i))
        

        nslip(4:6,i) = matmul(transpose(Rinv),slip(4:6,i))
        
    end do
    return 
    
end subroutine textr

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
subroutine decomp(Fe1,R)
    real(8), dimension(3,3) :: R, Fe1, Ev, U2,VL,Eval, U, Uinv
    integer :: n = 3, LWORK = 1000
    real(8), dimension(3) :: WR, WI, INFO
    real(8), dimension(1000) :: WORK

U2 = matmul(transpose(Fe1),Fe1)
    call dgeev('N','V',n,U2,n,WR,WI,VL,n,Ev,n,WORK,LWORK,INFO)
do i = 1,3
Eval (i,i) = sqrt(WR(i))
end do

U = matmul(transpose(Ev),matmul(Eval,Ev))

call minv3(U,Uinv)

R = matmul(Fe1,Uinv)

return
end subroutine decomp


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

subroutine eqvstr(T,sigmaeq)
    implicit none
    real(8) , intent(in), dimension(3,3) :: T
    real(8) , intent(out) :: sigmaeq
    integer :: i, j 

sigmaeq = 1.00/2.00*((T(1,1)-T(2,2))**2+(T(2,2)-T(3,3))**2+(T(3,3)-T(1,1))**2+6.00*(T(1,2)**2+T(2,3)**2+T(3,1)**2))

sigmaeq = sqrt(sigmaeq)
return
end subroutine eqvstr

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

subroutine gausseid(A,b,x0)
    implicit none
    real(8) , dimension (4,4) :: A, L , Aint, U
    real(8) , dimension(4) :: x0, x1, b,error,bint
    integer :: i

! Normalize the diagonal elements
Aint(1,1:4) = A(1,1:4)/A(1,1) 
Aint(2,1:4) = A(2,1:4)/A(2,2) 
Aint(3,1:4) = A(3,1:4)/A(3,3) 
Aint(4,1:4) = A(4,1:4)/A(4,4) 
Bint(1) = B(1)/A(1,1)
Bint(2) = B(2)/A(2,2)
Bint(3) = B(3)/A(3,3)
Bint(4) = B(4)/A(4,4)


!Define L And U 
L = 0
L(2,1) = Aint(2,1)
L(3,1:2) = Aint(3,1:2)
L(4,1:3) = Aint(4,1:3)

! Define U
U = 0
U(1,2:4) = Aint(1,2:4)
U(2,3:4) = Aint(2,3:4)
U(3,4) = Aint(3,4)

error = 1
   do while (sum(error)> 0.0000000001)
    x1 = 0
    do i = 1,4 
        x1(i) = bint(i)- dot_product(L(i,1:4),x1) - dot_product(U(i,1:4),x0)
    end do
   
    x0 = x1
    error = abs(x1-x0)
    
end do
write(11,*) x1
return
end subroutine


!! PERFORMS DOBLE DOT PRODUCT ON 3X3 TENSORS 
!!
!!
subroutine contract2(T,S,dprod)
    implicit none

    real(8), dimension(3,3), Intent(in) :: T, S
    real(8), Intent(out) :: dprod

    dprod = T(1,1)*S(1,1) + T(2,2)*S(2,2) + T(3,3)*S(3,3) +&
            T(1,2)*S(1,2) + T(2,1)*S(2,1) + T(1,3)*S(1,3) + &
            T(3,1)*S(3,1) + T(2,3)*S(2,3) + T(3,2)*S(3,2)
    
return
end subroutine contract2
