module constitutive
    !use global
    use mathmod
    use crystalplasticity

    




    contains
    subroutine constexpr(l,part,bryter,bcond,strain,tag,epsp,propconst,fid)

        use crystalplasticity
        implicit none

        real(8), dimension(3,3) :: tag, tagi
        real(8) :: epsp, epspi, strain
        real(8) , dimension(4)  :: sigma, offdl
        real(8), dimension(5) ::  offdl2, sigma2, IPIV2
        integer, dimension(6) :: pos1, pos2
        real(8) :: dt0,gammaskrank, dl
        real(8) , dimension(3,3)  :: Lb, Tagb, La, Dp,tagc, N
        real(8), dimension(6) :: propconst
        integer ::  switch , p,k,h,l, bryter,secit,part
        real(8) , dimension(4,4)  :: Jacob, Jinv
        
        integer :: LDA = 4,NRHS = 1, Info,  minl, maxl,nit,bcond,fid
        integer , dimension(4) :: IPIV
        real(8) , dimension(5,5)  :: Jacob2, Jinv2
        real(8) :: pwpercision
        logical :: consistent, consistentcontroll
        character*16 :: filename
        character*19 :: filename2
        gammaskrank = epsp
        pwpercision = 0.0000000001
        secit = 0
        
        dl = 0.0000001
        switch = 0
        consistent = .false.
        consistentcontroll = .false.
        if (bcond == 1 .and. bryter == 5 .or. bryter == 6) then
        write(filename,'("Dp_con_",I2,"_",I2,".txt")') fid , l   
        open(unit=fid+l+600, file=filename, status='unknown')
        write(filename2,'("Dp_con_ang",I2,"_",I2,".txt")') fid , l   
        open(unit=fid+800+l, file=filename2, status='unknown')
        else   
        write(filename,'("Dp_con_",I2,"_",I2,".txt")') fid , 99 
        open(unit=fid+l+600, file=filename, status='unknown')
        write(filename2,'("Dp_con_ang",I2,"_",I2,".txt")') fid , 99
        open(unit=fid+800+l, file=filename2, status='unknown')
        end if

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
          !  consistentcontroll = .true.
            pos1 = (/1, 1, 2, 3, 2, 1/)
                pos2 =(/2, 3, 3, 3, 2, 1/)
    case(2) 
        if (bcond == 2) then
    
            La(1,1) = propconst(6)
            La(1,2) = 0
            La(2,2) = propconst(5)
            La(2,1) = 0
            La(1,3) = 0
            La(3,1) = 0
            La(2,3) = 0
            La(3,2) = 0
            La(3,3) = -0.35*(La(1,1)+La(2,2))
            La = La/norm2(La)
            
            if ( abs(propconst(6)) >= abs(propconst(5))) then
                pos1 = (/1, 1, 2, 3, 2, 1/)
                pos2 =(/2, 3, 3, 3, 2, 1/)
                propconst = propconst/propconst(6)
            else if ( abs(propconst(6)) < abs(propconst(5))) then
                pos1 = (/1, 1, 2, 3, 1, 2/)
                pos2 =(/2, 3, 3, 3, 1, 2/)
                propconst = (/0.d+1, 0.d+1, 0.d+1, 0.d+1, propconst(6), propconst(5)/)/propconst(5)
            end if 
        
        end if   
    end select    
    
            dt0 = dt/10.
iter: do while (switch < 10000000)



Select case (bcond)

case (1)    
  consistent = .true.
consistentcontroll = .true.
!!! Iteration to ensure boundary condition is satisfied at througout all timesteps. 
    nit = 0   
    boundarycond: do  while (nit < 100)  
    !write(*,*) nit, norm2(La) 
    tagi = tag
    epspi = epsp
    consistent = consistentcontroll
       call yoshi3(Tagi,La,epspi,dt0,consistent,Dp)
 
       do h = 1,4
        sigma(h) = Tagi(pos1(h),pos2(h))
       end do
       !call sleep(1)
       !write(*,*) sigma ,epspi, consistent, consistentcontroll
         minl = minloc(sigma, DIM = 1)
         maxl = maxloc(sigma, DIM = 1)
         
         if (abs(sigma(minl)) < 0.00000000001 .and. abs(sigma(maxl)) < 0.00000000001) then
            consistentcontroll = consistent
            
            !write(*,*) sigma , Tagi(1,1), Tagi(2,2), epspi
            exit boundarycond
         end if 
!!!!!! Calculate Jacobian matrix in order to satisfy boundary condition
       do p = 1,4
            
                consistent = consistentcontroll
                tagb = tag
                epspi = epsp
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) + dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl 
                end if

                call yoshi3(Tagb,Lb,epspi,dt0,consistent,Dp)
            do k = 1,4    
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
   !call sleep(2)
  ! call sleep(2)
  ! write(*,*) 
   !write(*,*) 
   tagi = tag
   epspi = epsp
   consistent = consistentcontroll
   !write(*,*) tagi(1,1), tagi(2,2),epspi, consistent, consistentcontroll
  
   call yoshi3(Tagi,La,epspi,dt0,consistent,Dp)
   

    do h = 1,5
        sigma2(h) = Tagi(pos1(h),pos2(h))-propconst(h)*Tagi(pos1(6),pos2(6))
    end do
    !   write(*,*) La
       !write(*,*) epspi
       
         minl = minloc(sigma2, DIM = 1)
         maxl = maxloc(sigma2, DIM = 1)
         if (abs(sigma2(minl)) < 0.00000000001 .and. abs(sigma2(maxl)) < 0.00000000001) then
          consistentcontroll = consistent
          !write(*,*) tagi(1,1), tagi(2,2), epspi,consistent, consistentcontroll
            write(*,*) sigma2 , epspi
       !      write(*,*) tagi(1,1), tagi(2,2)
            exit boundary2
         end if 
     do p = 1,5
          
                
                
                consistent = consistentcontroll
                tagb = tag
                epspi = epsp
                Lb = La
                if (pos1(p) /= pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl
                Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) + dl
                else if (pos1(p) == pos2(p)) then
                Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) + dl 
                end if

                call yoshi3(Tagb,Lb,epspi,dt0,consistent,Dp)
               
                
                !consistent = consistentcontroll
                !tagc = tag
                !epspi = epsp
                !Lb = La
                !if (pos1(p) /= pos2(p)) then
                !Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) - dl
                !Lb(pos2(p),pos1(p)) = La(pos2(p),pos1(p)) - dl
                !else if (pos1(p) == pos2(p)) then
                !Lb(pos1(p),pos2(p)) = La(pos1(p),pos2(p)) - dl 
                !end if
                !call yoshi2(Tagc,Lb,epspi,dt0,consistent)
                
                do k = 1,5
            jacob2(k,p) = ((Tagb(pos1(k),pos2(k))-propconst(k)*Tagb(pos1(6),pos2(6)))- &
            (Tagi(pos1(k),pos2(k))-propconst(k)*tagi(pos1(6),pos2(6))))/dl
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
        La(pos1(5),pos2(5)) = La(pos1(5),pos2(5)) + offdl2(5)
      
        nit = nit+1

   end do boundary2
end select

if (bryter == 1 .or. bryter == 5 .or. bryter == 4) then
       
        
       
    if (epspi > strain .and. abs((epspi - strain)/strain) > pwpercision) then
        if (strain == 0) then
            dt0 = dt0/2
        else
        dt0 = (strain-epsp)*dt0/(epspi-epsp)*0.7
        end if
        secit = secit +1
        !write(*,*) 'Tidsskritt redusert'
        !write(*,*) dt0, epspi, epsp
        !write(*,*)
        switch = switch + 1
        if (secit > 30) then 
            write(*,*) epspi
            write(*,*) 'early exit'
            exit iter
        end if 
        cycle iter
    end if 
    
    tag = tagi
   ! write(*,*) epspi, nit, dt0, l, consistentcontroll
    epsp = epspi
    !write(*,*) epsp
    if (bryter == 5) then
        if (epsp > gammaskrank) then
            write(8,*) bryter, epspi
            
        
        write(11,*) Tag(1,3), Tag(1,2), Tag(2,3), Tag(3,3) ,epsp
        !call Yoshidamodel(tag,La,Dp)
        call hoshfordnormal(tag,N)
        
       write(fid+800+l,*) acos(contract2(La,Dp)/norm2(La)/norm2(Dp))*180/pi, acos(contract2(N,Dp)/norm2(N)/norm2(Dp))*180/pi, epsp
        !write(fid+40+l,*) acos(contract2(La,N)/norm2(La)/norm2(N))*180/pi, acos(contract2(N,Dp)/norm2(N)/norm2(Dp))*180/pi, epsp
        write(fid+l+600,*) Tag(1,1), Tag(2,2) , Dp(1,1)/norm2(La), Dp(2,2)/norm2(La) 
        gammaskrank = gammaskrank + dgamma
        
        end if
    end if 

    if (abs((epsp - strain)/strain) <= pwpercision) then
        !call Yoshidamodel(tag,La,Dp)
        write(18,*) Dp
        write(18,*) La-id*(La(1,1)+La(2,2)+La(3,3))/3
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
        !write(*,*) dt0, sigma, epspi-strain
       ! switch = switch +1
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
    ! write(*,*) epsp

    if (bryter == 6) then
        if (epsp > gammaskrank) then
            write(8,*) bryter, epspi
         if (switch  == 0 ) then
            !call Yoshidamodel(tag,La,Dp)
            write(18,*) Dp
            write(18,*) La-id*(La(1,1)+La(2,2)+La(3,3))/3
         end if
        write(11,*) Tag(1,1), tag(2,2), Tag(1,3),Tag(1,2), Tag(2,3), Tag(3,3),  epsp
        !call Yoshidamodel(tag,La,Dp)
        call hoshfordnormal(tag,N)
        !write(fid+800+l,*) acos(contract2(La,Dp)/norm2(La)/norm2(Dp))*180/pi, acos(contract2(N,Dp)/norm2(N)/norm2(Dp))*180/pi, epsp


        gammaskrank = gammaskrank + dgamma
        !write(fid+l+600,*) Tag(1,1), Tag(2,2) , Dp(1,1)/norm2(La), Dp(2,2)/norm2(La) 
        !write(8,*) Tag(1,1),Tag(2,2), Dp(1,1)/sqrt(Dp(1,1)**2+Dp(2,2)**2),Dp(2,2)/sqrt(Dp(1,1)**2+Dp(2,2)**2.)
        
        end if
    end if 
    
    if (strain /= 0 .and. abs((epsp -strain)/strain)<= pwpercision) then
        call hoshfordnormal(tag,N)
        write(fid+800+l,*) acos(contract2(La,Dp)/norm2(La)/norm2(Dp))*180/pi, acos(contract2(N,Dp)/norm2(N)/norm2(Dp))*180/pi, epsp
        write(fid+l+600,*) Tag(1,1), Tag(2,2) , Dp(1,1)/norm2(La), Dp(2,2)/norm2(La) 
        exit iter
    else if (strain == 0 .and. epsp <= pwpercision .and. epsp > 0) then
       ! write(*,*) 'check'
        exit iter

    end if
end if 
switch = switch +1 
    end do iter
!close(unit=fid+l+600)
!close(unit=fid+l+800)
end subroutine constexpr

subroutine elasticsolution(D,tag)
use crystalplasticity
    implicit none
real(8) sigma, sigma0, G, dt0, consistency, sigma1, dt1
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

subroutine yoshi2(tag,D,epsp,dt0,consistent,Dp)
    use crystalplasticity
    use mathmod
      implicit none

      real(8), dimension(3,3,3,3) :: Cep, dyadic, T, CT,Cel4
      real(8), dimension(3,3) :: N, D, tag, dtag,Nnorm, Dtan, tagint, Ddev, tagcheck,tag1, Dp, De
      real(8) :: h ,h2, G,NCN, theta, alpha, sigma, sigma0, dt0, epsp, lambdadot, lambdadot2
      real(8) :: feps,feps2 , sigmacheck, dt1,epsp2, consistency, dl ,modelnum
      real(8) , dimension(6,6) :: Chook
      integer :: i,j,k,l,p,q,Info
      logical :: consistent
      real(8) , dimension(6) :: Dvec, hvec, hvec2
      real(8) , dimension(7) :: newtvec, solution,IPIV
      real(8) , dimension(7,7) :: Jacobi
      integer , dimension(6) :: pos1, pos2
      modelnum = 2
        pos1 = (/1, 2, 3, 2, 1, 1/)
        pos2 = (/1, 2, 3, 3, 3, 2/)
      dl = 0.000001
     
      sigma0 = gaveps(epsp)
      h = haveps(epsp)
      dtag = 0
      Dp = 0
!if (modelnum == 1) then
!    theta0 = pi/18.
!else if (modelnum == 2) then
!    theta0 = pi/18./2.
!end if 

   
      call Elasticconstant(Chook,G)

      call eqvstress(tag,sigma)
      
      Call hoshfordnormal(tag,N)
      
  
 
 
  

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
  if (consistent .eqv. .false.) then
  !call eqvstresslse.(tag,sigma)
  
  
      
      do i = 1,3
          do j = 1,3
              do k =1,3
                  do l = 1,3
                      dtag(i,j) = dtag(i,j) + Cel4(i,j,k,l)*D(l,k)
                  end do
              end do
          end do
      end do

      if (contract2(dtag,N) < 0 ) then
          stop 'Unloading'
      !write(*,*) 'Unloading'
      end if 
      tagcheck = tag+dtag*dt0
      !write(*,*) sigma-sigma0, sigma , sigma0
  call eqvstress(tagcheck,sigmacheck)
  
  if (sigmacheck > sigma0 ) then
      consistency = abs(sigmacheck - sigma0)
      dt1 = dt0
     
      do while (consistency > 0.0000000000001) 
      !do i = 1,20
          tagcheck = tag+dtag*dt1
          !write(*,*) tagcheck - tag
          call eqvstress(tagcheck,sigmacheck)
          consistency = abs(sigmacheck - sigma0)
         
          
          dt1 = (sigma0-sigma)/(sigmacheck-sigma)*dt1
          
         !write(*,*) dt1, sigma, sigmacheck, sigma0
          if (sigmacheck < sigma0) then
          tag = tagcheck
          sigma = sigmacheck
          end if
      end do
      tag = tagcheck
      sigma = sigmacheck
      consistent = .true.
      !write(*,*) consistency
 

  else
      tag = tagcheck
  end if 
  
  !write(*,*) sigmacheck - sigma0
  Dp = 0
else
!function Dp(sigma
    q = 0
!end function
! Initial guess, elastic predictor
    newtvec = 1
Ddev = D - id*(D(1,1)+D(2,2)+D(3,3))/3    
tag1 = tag
epsp2 = epsp
Call hoshfordnormal(tag1,N)
lambdadot = contract2(N,vec2tens(matmul(Chook,tens2vec(D))))&
/(contract2(N,vec2tens(matmul(Chook,tens2vec(N))))+sqrt(2./3.)*norm2(N)*h)
!write(*,*) newtvec
do while (norm2(newtvec) > 0.000000001)

Call hoshfordnormal(tag1,N)
Nnorm = N/norm2(N)
!Calculate theta

theta = acos(contract2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.,N)/norm2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.)/Norm2(N))


if (theta >= 0 .and. theta <= theta0 ) then
        alpha = 1-c1*sqrt(sigma0/G)-c2*macauley(h)/G
else if (theta > theta0 .and. theta < pi/2) then
    if (modelnum == 1) then 
        alpha = (1-c1*sqrt(sigma0/G)-c2*macauley(h)/G)*((pi/2-theta)/(pi/2-theta0))
    else if (modelnum == 2) then
        alpha = (1-c1*sqrt(sigma0/G)-c2*macauley(h)/G)*Tan(c3*(theta-theta0)+theta0)/tan(theta)
    end if
    !alpha = (pi/2-theta)/(pi/2-theta0)
else if ( theta >= pi/2 .and. theta < pi ) then
    alpha = 0
    write(*,*) 'Warning unloading'
end if   
    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                    T(i,j,k,l) = 1./2.*(kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k))-1./3.*(id(i,j)*id(k,l)) &
                                -Nnorm(i,j)*Nnorm(k,l)     
                end do
            end do
        end do
    end do
   Dtan = 0
    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                   Dtan(i,j) = Dtan(i,j)+ T(i,j,k,l)*Ddev(l,k)
                end do
            end do
        end do
    end do
!lambdadot = contract2(N,vec2tens(matmul(Chook,tens2vec(D))))&
!/(contract2(N,vec2tens(matmul(Chook,tens2vec(N))))+sqrt(2./3.)*norm2(N)*h)
!h = haveps(epsp2)
    epsp2 = epsp + lambdadot*sqrt(2./3.)*norm2(N)*dt0
    call eqvstress(tag1,sigmacheck)
feps = sigmacheck -gaveps(epsp2)
Dp = lambdadot*N+alpha*Dtan
De = D-lambdadot*N-alpha*Dtan
hvec = tens2vec(tag1)-tens2vec(tag) - matmul(Chook,(/ De(1,1), De(2,2),De(3,3),2*De(2,3), 2*De(1,3),2*De(1,2) /))*dt0
!feps = epsp2 - epsp - lambdadot*sqrt(2./3.)*norm2(N)*dt0
!write(*,*) hvec

newtvec(1:6) = hvec
newtvec(7) = feps
!write(*,*) newtvec
!!!! Calculate jacobian

do i = 1,7
    tagint = tag1
    epsp2 = epsp
    
    if (i < 7) then
    if (pos1(i) /= pos2(i)) then 
    tagint(pos1(i),pos2(i)) = tagint(pos1(i),pos2(i))+ dl
    tagint(pos2(i),pos1(i)) = tagint(pos2(i),pos1(i))+ dl
    else
    tagint(pos1(i),pos2(i)) = tagint(pos1(i),pos2(i))+ dl
    end if 
            Call hoshfordnormal(tagint,N)
            theta = acos(contract2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.,N)/norm2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.)/Norm2(N))


            if (theta >= 0 .and. theta <= theta0 ) then
                alpha = 1-c1*sqrt(gaveps(epsp2)/G)-c2*macauley(h)/G
            else if (theta > theta0 .and. theta < pi/2) then
                        if (modelnum == 1) then 
                            alpha = (1-c1*sqrt(gaveps(epsp2)/G)-c2*macauley(haveps(epsp2))/G)*((pi/2-theta)/(pi/2-theta0))
                        else if (modelnum == 2) then
                            alpha = (1-c1*sqrt(gaveps(epsp2)/G)-c2*macauley(haveps(epsp2))/G)* &
                            tan(c3*(theta-theta0)+theta0)/tan(theta)
                        end if
            else if ( theta >= pi/2 .and. theta < pi ) then
                        alpha = 0
                        write(*,*) 'Warning unloading'
            end if   

            Nnorm = N/norm2(N)
                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                        T(p,j,k,l) = 1./2.*(kronecker(p,k)*kronecker(j,l)+kronecker(p,l)*kronecker(j,k))-1./3.*(id(p,j)*id(k,l)) &
                                                -Nnorm(p,j)*Nnorm(k,l)     
                                end do
                            end do
                        end do
                    end do
                Dtan = 0
                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                                    Dtan(p,j) = Dtan(p,j)+ T(p,j,k,l)*Ddev(l,k)
                                end do
                            end do
                        end do
                    end do
                  !  h2 = haveps(epsp2)
        !lambdadot2 = contract2(N,vec2tens(matmul(Chook,tens2vec(D)))) &
        !/(contract2(N,vec2tens(matmul(Chook,tens2vec(N))))+sqrt(2./3.)*norm2(N)*h2)
                   ! feps2 = epsp2 - epsp - lambdadot2*sqrt(2./3.)*norm2(N)*dt0
                    epsp2 = epsp + lambdadot*sqrt(2./3.)*norm2(N)*dt0
       call eqvstress(tagint,sigmacheck)
feps2 = sigmacheck -gaveps(epsp2)

De = D-lambdadot*N-alpha*Dtan
hvec2 = tens2vec(tagint)-tens2vec(tag) - matmul(Chook,(/De(1,1), De(2,2),De(3,3),2*De(2,3), 2*De(1,3),2*De(1,2) /))*dt0
        do k = 1,7
            if (k < 7) then
                jacobi(k,i) = (hvec2(k)-hvec(k))/dl
            else 
                jacobi(k,i) = (feps2-feps)/dl 
            end if
        end do

        else 
            Call hoshfordnormal(tag1,N)
            theta = acos(contract2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.,N)/norm2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.)/Norm2(N))
            
            if (theta >= 0 .and. theta <= theta0 ) then
                             alpha = 1-c1*sqrt(gaveps(epsp2)/G)-c2*macauley(haveps(epsp2))/G
            else if (theta > theta0 .and. theta < pi/2) then
                            if (modelnum == 1) then 
                                alpha = (1-c1*sqrt(gaveps(epsp2)/G)-c2*macauley(haveps(epsp2))/G)*((pi/2-theta)/(pi/2-theta0))
                            else if (modelnum == 2) then
                                alpha = (1-c1*sqrt(gaveps(epsp2)/G)-c2*macauley(haveps(epsp2))/G)*&
                                tan(c3*(theta-theta0)+theta0)/tan(theta)
                            end if
            else if ( theta >= pi/2 .and. theta < pi ) then
                            alpha = 0
                            write(*,*) 'Warning unloading'
            end if   
            
            
            Nnorm = N/norm2(N)
                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                    T(p,j,k,l) = 1./2.*(kronecker(p,k)*kronecker(j,l)+kronecker(p,l)*kronecker(j,k))-1./3.*(id(p,j)*id(k,l)) &
                    -Nnorm(p,j)*Nnorm(k,l)     
                                end do
                            end do
                        end do
                    end do
                Dtan = 0
                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                                Dtan(p,j) = Dtan(p,j)+ T(p,j,k,l)*Ddev(l,k)
                                end do
                            end do
                        end do
                    end do
                   ! h2 = haveps(epsp2+dl)
       ! lambdadot2 = contract2(N,vec2tens(matmul(Chook,tens2vec(D)))) &
       ! /(contract2(N,vec2tens(matmul(Chook,tens2vec(N))))+sqrt(2./3.)*norm2(N)*h2)
        
       ! epsp2 = epsp+lambdadot2*sqrt(2./3.)*norm2(N)*dt0

        epsp2 = epsp + (lambdadot+dl)*sqrt(2./3.)*norm2(N)*dt0
        call eqvstress(tag1,sigmacheck)
        feps2 = sigmacheck -gaveps(epsp2)
      
        De = D-(lambdadot+dl)*N-alpha*Dtan
hvec2 = tens2vec(tag1)-tens2vec(tag) - matmul(Chook,(/ De(1,1), De(2,2),De(3,3),2*De(2,3), 2*De(1,3),2*De(1,2) /))*dt0
        do k = 1,7
            if (k < 7) then
                jacobi(k,i) = (hvec2(k)-hvec(k))/dl
            else
                jacobi(k,i) = (feps2-feps)/dl 
            end if
        end do
    end if 
end do
solution = -newtvec
    call dgesv(7,1,jacobi,7,IPIV,solution,7 , Info)
   
    tag1 = tag1+vec2tens(solution(1:6))
    lambdadot = lambdadot+solution(7)
!write(*,*) solution
    if (q > 10) then
       ! write(*,*) 'not fully converged solution'
       ! write(*,*) newtvec
    end if
    if (q > 15 .and. norm2(newtvec) < 1e-5) then
        write(*,*) 'not fully converged solution'
        exit 
    end if
    q = q+1
end do
tag = tag1
epsp = epsp2
!call hoshfordnormal(tag,N)
!Dp = lambdadot*N+alpha*Dtan
end if 

 ! call eqvstress(tag,sigma)
 ! write(*,*) sigma,gaveps(epsp), sigma-gaveps(epsp)
 !
 ! write(*,*) dtag
 ! write(*,*) 


  return
  end subroutine yoshi2


  subroutine yoshi3(tag,D,epsp,dt0,consistent,Dp)
    use crystalplasticity
    use mathmod
    use global
      implicit none

      real(8), dimension(3,3,3,3) :: Cep, dyadic, T, CT,Cel4
      real(8), dimension(3,3) :: N, tag, dtag,Nnorm, Dtan, tagint, Ddev, tagcheck,tag1, Dp, De, Dcorr
      real(8) :: h ,h2, G,epspelpred, theta, alpha, sigma, sigma0, dt0, epsp, lambdadot, lambdadot2
      real(8) :: feps,feps2 , sigmacheck, dt1,epsp2, sigmaelpred, dl ,modelnum
      real(8) , dimension(6,6) :: Chook
      integer :: i,j,k,l,p,q,Info
      logical :: consistent
      real(8) , dimension(6) :: Dvec, hvec, hvec2
      real(8) , dimension(7) :: newtvec, solution,IPIV
      real(8) , dimension(7,7) :: Jacobi
      integer , dimension(6) :: pos1, pos2
      real(8), intent(in), dimension(3,3) :: D
      !real(8) :: c1= 0.3, c2 = 0.5, c3 = 0.4 ,theta0 

      modelnum = 2
        pos1 = (/1, 2, 3, 2, 1, 1/)
        pos2 = (/1, 2, 3, 3, 3, 2/)
      dl = 0.000001
     
      sigma0 = gaveps(epsp)
      h = haveps(epsp)
      dtag = 0
      Dp = 0


   
      call Elasticconstant(Chook,G)
      !write(*,*) G
      call eqvstress(tag,sigma)
      
      Call hoshfordnormal(tag,N)
      
  
 
 
  

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
  
!!!! Calculate elastic predictor    
  do i = 1,3
    do j = 1,3
        do k =1,3
            do l = 1,3
                dtag(i,j) = dtag(i,j) + Cel4(i,j,k,l)*D(l,k)
            end do
        end do
    end do
end do
!!! Check for unloading
if (contract2(dtag,N) < 0 ) then
   write(*,*) 'Unloading'
end if 
!!! update stress and calculate equivalent stress
tagcheck = tag+dtag*dt0
call eqvstress(tagcheck,sigmaelpred)


!!!! Check if yield strength is surpassed, if the yield stress is reached plastic correction is needed.
if (sigmaelpred < gaveps(epsp)) then
tag = tagcheck



else if (sigmaelpred >= gaveps(epsp)) then
q = 0

! Initial guess, elastic predictor
newtvec = 1
Ddev = D - id*(D(1,1)+D(2,2)+D(3,3))/3    
tag1 = tag
epsp2 = epsp
Call hoshfordnormal(tag1,N)
!!! This equation is a rough estimate, believe voigt notation is not appropriate for the unsymmetric N
lambdadot = contract2(N,vec2tens(matmul(Chook,(/ D(1,1), D(2,2),D(3,3),2*D(2,3), 2*D(1,3),2*D(1,2) /))))&
/(contract2(N,vec2tens(matmul(Chook,tens2vec(N))))+sqrt(2./3.)*norm2(N)*h)


!!! performs newton-raphson until solution is sufficently converged. 
do while (norm2(newtvec) > 0.000000001)

!!!! First calculate using current stress point+initial guess lambdadot, and checks if solution is sufficiently converged. 
!!!! The same step is performed with the updated stress and lambdadot after the solution has been updated from the newton raphson iteration
Call hoshfordnormal(tag1,N)
Nnorm = N/norm2(N)

    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                    T(i,j,k,l) = 1./2.*(kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k))-1./3.*(id(i,j)*id(k,l)) &
                                -Nnorm(i,j)*Nnorm(k,l)     
                end do
            end do
        end do
    end do
   Dtan = 0
    do i = 1,3
        do j = 1,3
            do k =1,3
                do l = 1,3
                   Dtan(i,j) = Dtan(i,j)+ T(i,j,k,l)*Ddev(l,k)
                end do
            end do
        end do
    end do

    theta = acos(contract2(Ddev,N)/norm2(Ddev)/Norm2(N))
    epsp2 = epsp + lambdadot*sqrt(2./3.)*norm2(N)*dt0
    alpha = alphacoeff(theta,epsp2,modelnum,G)
    call eqvstress(tag1,sigmacheck)
feps = sigmacheck -gaveps(epsp2)
Dp = lambdadot*N+alpha*Dtan
Dcorr = Dp
hvec = tens2vec(tag1)-tens2vec(tagcheck) + &
       matmul(Chook,(/ Dcorr(1,1), Dcorr(2,2),Dcorr(3,3),2*Dcorr(2,3), 2*Dcorr(1,3),2*Dcorr(1,2) /))*dt0

newtvec(1:6) = hvec
newtvec(7) = feps
!write(*,*) newtvec
!!!! Calculate jacobian
!!!  The six stress components are used, toghether with lambda dot

do i = 1,7
    tagint = tag1
    epsp2 = epsp
    
    !!  First is to calculate derivatives of the six stress components. 
    if (i < 7) then
    if (pos1(i) /= pos2(i)) then 
    tagint(pos1(i),pos2(i)) = tagint(pos1(i),pos2(i))+ dl
    tagint(pos2(i),pos1(i)) = tagint(pos2(i),pos1(i))+ dl
    else
    tagint(pos1(i),pos2(i)) = tagint(pos1(i),pos2(i))+ dl
    end if 
            Call hoshfordnormal(tagint,N)
            Nnorm = N/norm2(N)
            theta = acos(contract2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.,N)/norm2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.)/Norm2(N))

                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                        T(p,j,k,l) = 1./2.*(kronecker(p,k)*kronecker(j,l)+kronecker(p,l)*kronecker(j,k))-1./3.*(id(p,j)*id(k,l)) &
                                                -Nnorm(p,j)*Nnorm(k,l)     
                                end do
                            end do
                        end do
                    end do
                Dtan = 0
                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                                    Dtan(p,j) = Dtan(p,j)+ T(p,j,k,l)*Ddev(l,k)
                                end do
                            end do
                        end do
                    end do
epsp2 = epsp + lambdadot*sqrt(2./3.)*norm2(N)*dt0
alpha = alphacoeff(theta,epsp2,modelnum,G)
call eqvstress(tagint,sigmacheck)
feps2 = sigmacheck -gaveps(epsp2)

    Dcorr = lambdadot*N + alpha*Dtan
    hvec2 = tens2vec(tagint)-tens2vec(tagcheck) + &
            matmul(Chook,(/Dcorr(1,1), Dcorr(2,2),Dcorr(3,3),2*Dcorr(2,3), 2*Dcorr(1,3),2*Dcorr(1,2) /))*dt0

!!!! Calculate the derivatives of each of the stress equations and the consistency equation. 
        do k = 1,7
            if (k < 7) then
                jacobi(k,i) = (hvec2(k)-hvec(k))/dl
            else 
                jacobi(k,i) = (feps2-feps)/dl 
            end if
        end do
!!!!! The last parameter is lambda, which must be calculated separately
        else 
            Call hoshfordnormal(tag1,N)
            Nnorm = N/norm2(N)
            theta = acos(contract2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.,N)/norm2(D-id*(D(1,1)+D(2,2)+D(3,3))/3.)/Norm2(N))
            

                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                    T(p,j,k,l) = 1./2.*(kronecker(p,k)*kronecker(j,l)+kronecker(p,l)*kronecker(j,k))-1./3.*(id(p,j)*id(k,l)) &
                    -Nnorm(p,j)*Nnorm(k,l)     
                                end do
                            end do
                        end do
                    end do
                Dtan = 0
                    do p = 1,3
                        do j = 1,3
                            do k =1,3
                                do l = 1,3
                                Dtan(p,j) = Dtan(p,j)+ T(p,j,k,l)*Ddev(l,k)
                                end do
                            end do
                        end do
                    end do

        epsp2 = epsp + (lambdadot+dl)*sqrt(2./3.)*norm2(N)*dt0
        alpha = alphacoeff(theta,epsp2,modelnum,G)
        call eqvstress(tagint,sigmacheck)
        feps2 = sigmacheck -gaveps(epsp2)
        Dcorr = (lambdadot+dl)*N+alpha*Dtan
        hvec2 = tens2vec(tag1)-tens2vec(tagcheck) + &
                matmul(Chook,(/ Dcorr(1,1), Dcorr(2,2),Dcorr(3,3),2*Dcorr(2,3), 2*Dcorr(1,3),2*Dcorr(1,2) /))*dt0
        do k = 1,7
            if (k < 7) then
                jacobi(k,i) = (hvec2(k)-hvec(k))/dl
            else
                jacobi(k,i) = (feps2-feps)/dl 
            end if
        end do
    end if 
end do
solution = -newtvec
    call dgesv(7,1,jacobi,7,IPIV,solution,7 , Info)
   
    tag1 = tag1+vec2tens(solution(1:6))
    lambdadot = lambdadot+solution(7)
!write(*,*) solution
    if (q > 10) then
       ! write(*,*) 'not fully converged solution'
        !write(*,*) newtvec
    end if
    if (q > 15 ) then
        write(*,*) 'not fully converged solution'
        write(*,*) newtvec
        exit 
    end if
    q = q+1
end do
tag = tag1
epsp = epsp2
!call hoshfordnormal(tag,N)
!Dp = lambdadot*N+alpha*Dtan
end if 

 ! call eqvstress(tag,sigma)
 ! write(*,*) sigma,gaveps(epsp), sigma-gaveps(epsp)
 !
 ! write(*,*) dtag
 ! write(*,*) 


  return
  end subroutine yoshi3



 function modelerror(c_1,c_2,c_3,theta_0)
    implicit none
    real(8)                  :: c_1,c_2,c_3,theta_0, modelerror
    real(8)                  :: epsp, pw, cpstrain
    real(8), dimension(3,3)  :: tag, cptag, devtag, devcptag
    real(8), dimension(6)    :: propconst
    integer                  :: fid , nlines, k, bryter, bcond
    real(8), dimension(9)    :: cptagv
    !character*15 :: filename
    !threadnum = 1
    !!!! Set model parameters
    c1 = c_1
    c2 = c_2
    c3 = c_3
    theta0 = theta_0
    tag = 0
    tag(1,1) = 49.103557586671037
    tag(2,2) = 49.103557586671037
    tag(3,3) = 1.1107815788452554E-012
    epsp = 0.02

    !open(action='read',unit=15,file="stresscon.txt",status='old')
    !read(15,*) tag(1,1:3), tag(2,1:3), tag(3,1:3), epsp
    !close(15)
    !!!! perform experiment pre strain path change
    !epsp = 0
    !tag = 0
    !pw = 0.02
    !bcond = 2
    !bryter = 1
    !propconst = (/ 0.d+0, 0.d+0, 0.d+0, 0.d+0, 1.d+0, 1.d+0 /)
    !fid = 20
    !call constexpr(0,16,bryter, bcond,pw, tag, epsp, propconst,fid)

    
    !threadnum = omp_get_thread_num()+1 
    !write(filename,'("stress_cp_",I1,".txt")') threadnum
    !write(*,*) threadnum, filename
    
    !!!! Strain path change experiment
    !!!! Open file containing cp data
    call countlin('stress_cp.txt',nlines)
    open(action='read',unit=17,file='stress_cp.txt',status='old')
    
    
    modelerror = 0
    bcond = 1
    bryter = 2
    fid = 21
    do k = 1,15
        if (k<=5) then
        pw = 0.02+0.00005*k
        else 
        pw = 0.02+0.0005*(k-5)
        end if     
        call constexpr(0,16,bryter,bcond,pw, tag, epsp,propconst,fid)
        devtag = tag - (tag(1,1)+tag(2,2)+tag(3,3))/3
        
        !!!! Check that strain is equal
        read(17,*) cptagv , cpstrain
       
        if (abs(cpstrain-epsp) > 0.0000001) then
            write(*,*) 'unequal strains'
            write(*,*) cpstrain , epsp
        end if 
        cptag(1,1:3) = cptagv(1:3)
        cptag(2,1:3) = cptagv(4:6)
        cptag(3,1:3) = cptagv(7:9)
        devcptag = cptag - (cptag(1,1)+ cptag(2,2) + cptag(3,3))/3
       ! write(*,*) 
        modelerror = modelerror + abs(acos(contract2(devtag,devcptag)/norm2(devtag)/norm2(devcptag)))/nlines
        !modelerror = modelerror + norm2(cptag-tag)/nlines
        !write(*,*) modelerror, k
    end do 
    close(17) 
    return  
end function modelerror

end module constitutive