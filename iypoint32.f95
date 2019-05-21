program ypoint
use crystalplasticity
use constitutive
    implicit none

integer :: part,bryter, k, i,bcond, fid, allocatestatus, fid0
real(8) :: t1,t2,omp_get_wtime,pw1,pw2,epsp,epsp1,cpstrain, error,cpstrain0
real(8) , dimension(3,3) :: tag, tag1, N
real(8) , dimension(6) :: propconst
!real(8) , dimension(:,:), Allocatable ::eul
real(8) , dimension(:,:,:), Allocatable  :: F0,F01,Fp0,Fp01
real(8),  dimension(:,:), Allocatable :: S0,S01
real(8), dimension(4) :: solution, initial
t1 = omp_get_wtime()
open(unit=11,file="result.txt",status='replace')
open(unit=5,file="alphacp.txt",status='replace')
open(unit=3,file="Stress.txt",status='replace')
open(unit=8,file="eulerangles.txt",status='replace')
open(unit=13,file="Dp_cp.txt",status='replace')
open(unit=18,file="Dp_con.txt",status='replace')
open(unit=14,file="Dp2_cp.txt",status='replace')
open(unit=16,file="Grad_cp.txt",status='replace')
open(unit=4,file="Cep.txt",status='replace')
call init()
write(*,*) nlines

Allocate(F0(3,3,nlines),Fp0(3,3,nlines),F01(3,3,nlines),Fp01(3,3,nlines), stat=allocatestatus)
if (allocatestatus /= 0) stop "Not enough memory"
Allocate(s0(nlines,12),s01(nlines,12), stat=allocatestatus)
if (allocatestatus /= 0) stop "Not enough memory"

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
!!!   ------------------------------
!!!   bryter = 5 or 6 correspond to 1 or 4, only difference is that data is written to file, ie Dp, Gradient, Angles etc, hence these are used to perform path change experiments. 
!!!
!!!   bryter = 5 - Calculates point on the yield surface corresponding to the prescribed L11,L22 and plastic work and updates F0,F0p,S0
!!!                WRITES TO FILE, Initialize an undeformed crystal
!!!
!!!   bryter = 6 - Takes a predeformed crystal (F0,Fp0,S0) and perform given deformation in stress or strain controll to a desired plastic strain. 
!!!
!!!   ---------------------------
!!!   bryter = 7 - Performs deformation in either stress or strain controll from an UNDEFORMED crystal to a prescribed plastic strain. Used to calculate yield surfaces. 
!!!                Writes the calculated stress point to the file "stress.txt"
!!!
!!!   bryter = 8 - Performs deformation in either stress or strain controll from an PREDEFORMED crystal to a prescribed plastic strain. Used to calculate instantaneous yield surfaces. 
!!!                Writes the calculated stress point to the file "stress.txt"



fid0 = 49

initial = (/3.2d0, 0.5d+0, 0.5d+0, pi/9 /)
call steepestdecent(solution, initial)
c1 = 5
c2 = 0.5
c3 = 0.5
theta0 = pi/10.
tag = 0
epsp = 0
!tag(1,1) = 49.103557586671037
!tag(2,2) = 49.103557586671037
!tag(3,3) = 1.1107815788452554E-012
!epsp = 0.02
!error = modelerror(9.d2, 0.6d0, 0.5d0, pi/9)
!write(*,*) error
fid = fid0
!!propconst = (/0.0, 0.0, 0.0, 0.0, 22, 11/)
propconst = (/0.0, 0.0, 0.0, 0.0, 1.0, 1.0/)
bryter = 5
k = 1
cpstrain = 0
pw1 = 0.02
bcond = 2
!call constexpr(k,16,bryter, bcond,pw1, tag, epsp, propconst,fid)
!call newton(0,16,bryter,bcond,F0,Fp0,S0,pw1,cpstrain,propconst,fid)   
tag1 = tag
epsp1 = epsp
cpstrain0 = cpstrain
!write(*,*) cpstrain
dt = 0.00001
part = 100
bryter = 6
F01 = F0
Fp01 = Fp0
S01 = S0
call OMP_SET_NUM_THREADS(7)
!$OMP PARALLEL PRIVATE(propconst,k, tag,epsp,fid,bryter,cpstrain,S0,Fp0,F0)
!$OMP DO
do k = -16,14,2
    tag = tag1
cpstrain = cpstrain0
    epsp = epsp1
    S0 = s01
    Fp0 = Fp01
    F0 = F01
    !if (k<=5 ) then 
    !pw1 = 0.02+0.00005*k
    !else
    !pw1 = 0.02+0.0005*(k-5)
    !end if
   pw1 = 0.03
    bcond = 1
    bryter = 6
    fid = fid0
     !!propconst = (/0.0, 0.0, 0.0, 0.0, 22, 11/)
    propconst = (/0.d+0, 0.d+0, 0.d+0, 0.d+0, sin(2*pi*k/part), cos(2*pi*k/part)/)
    !write(*,*) 'start'
    !write(*,*) k
    !call constexpr(0,64,bryter,bcond,pw1, tag, epsp,propconst,fid)
    write(*,*) tag(1,1), tag(2,2) , k
    fid = fid0
    !call newton(k,64,bryter,bcond,F0,Fp0,S0,pw1,cpstrain,propconst,fid) !
    tag1 = tag
    epsp1 = epsp
    !cpstrain0 = cpstrain
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
i = 4
bcond = 2
!part = 4d

   ! write(*,*) tag(1,1), tag(2,2) , k

call OMP_SET_NUM_THREADS(2)
!!!$OMP PARALLEL DO PRIVATE(F0,S0,Fp0,bryter,propconst,bcond,pw1,pw2,k,F01,Fp01,S01,tag1,epsp1, tag,epsp,fid,part,i)
!!$OMP PARALLEL PRIVATE(F0,S0,Fp0,bryter,propconst,bcond,pw1,pw2,k &
!!$OMP                 ,F01,Fp01,S01,tag1,epsp1, tag,epsp,fid,part,i, &
!!$OMP                   slip,Cel,eul,nlines,id,R, hardening, pi, dt )
!!$OMP DO
!do k = -4,part
!    bryter = 6
!    tag = tag1
!    epsp = epsp1
!    pw1 = 0.008
!    bcond = 1
!    fid = 20
!    call constexpr(k,4,bryter,bcond,pw1, tag, epsp,propconst,fid)
!    write(*,*) tag(1,1), tag(2,2) , k
!    
!   Fp0 = Fp01
!   S0 = S01
!   F0 = F01
!    bryter = 6
!    fid = 20
!    pw1 = 0.008
!    pw2 = 0.006
!    bcond = 1
!!write(*,*) Fp0(1:3,1,1)

  ! propconst = (/0.0, 0.0, 0.0, 0.0, -0.2*k/)
!    call newton(k,i,bryter,bcond,F0,Fp0,S0,pw2,propconst,fid)   
!end do
!!$OMP END DO
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




subroutine steepestdecent(solution, initial)
    use constitutive
    use crystalplasticity
    implicit none
    real(8) :: error1, error2, tol, error, dl, gamma, beta, alpha
    real(8), dimension(4) :: gradient, gradient1, initial, solution, tempsol,solution1, p
    integer :: k, i , counter
    
solution = initial
dl = 0.000001
    !!! calculate gradient 
error = modelerror(solution(1),solution(2),solution(3), solution(4))
write(*,*) error
do i = 1,100


!write(*,*) error
if ( norm2(gradient) < 0.000000001 ) then
  !  exit
end if
gradient = 0
!call OMP_SET_NUM_THREADS(1)
!!$OMP PARALLEL PRIVATE(c1,c2,c3, theta0,tempsol)
!!$OMP DO
do k = 1,4
    tempsol= solution
    !write(*,*) k
   tempsol(k) = tempsol(k) + dl 
if (k == 3 .and. tempsol(k) > 1) then
    tempsol(k) = solution(k) - dl
    gradient(k) =-(modelerror(tempsol(1), tempsol(2), tempsol(3), tempsol(4))-error)/dl
   else 
   write(*,*) modelerror(tempsol(1), tempsol(2), tempsol(3), tempsol(4))-error
   gradient(k) = (modelerror(tempsol(1), tempsol(2), tempsol(3), tempsol(4))-error)/dl
   end if
end do
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL
write(*,*) 'gradient'
write(*,*) gradient
!!!! Perform line search
p = -gradient/norm2(gradient)

alpha = 0.5
beta = 0.5
error2 = error
!write(*,*) error2 - error + alpha*beta*dot_product(p,gradient)
!write(*,*) p
counter = 0
linesearch: do 
    tempsol = solution + alpha*p
    write(*,*) p, alpha
    !write(*,*) alpha
    !!!!! check if solution is within range.

    if (tempsol(1) < 0 ) then
        if (solution(1) > 0.00000001 ) then
            alpha = -solution(1)/p(1)-0.00000001
            cycle  linesearch
                else
             p(1) = 0
             p = p/norm2(p)
             cycle  linesearch
        end if
    end if
    if (tempsol(2) < 0 ) then
        if (solution(2) > 0.00000001 ) then
            alpha = -solution(2)/p(2)-0.00000001
            cycle  linesearch
                else
            p(2) = 0
            p = p/norm2(p)
            cycle  linesearch
        end if
    end if
    if (tempsol(3) < 0) then
        if (solution(3) > 0.00000001 ) then
            alpha = -solution(3)/p(3)-0.00000001
            cycle linesearch
            else
                p(3) = 0
                p = p/norm2(p)
                cycle linesearch
        end if
    end if
    if (tempsol(3) > 1) then
        
        if (solution(3) < 1-0.00001 ) then
        alpha = (1-solution(3))/p(3)-0.00000001
        cycle linesearch
        else 
            p(3) = 0   
            p = p/norm2(p)
            cycle linesearch
        end if
    end if
    if (tempsol(4) < 0) then
        if (solution(4) > 0.00000001 ) then
            alpha = -solution(4)/p(4)-0.00000001
            cycle  linesearch
                else
             p(4) = 0
             p = p/norm2(p)
             cycle  linesearch
        end if
    end if


    error2 = modelerror(tempsol(1), tempsol(2), tempsol(3), tempsol(4))
  
if (counter <= 10) then
    if (error2 < error + alpha*beta*dot_product(p,gradient) ) then
        !write(*,*) error2 - error - alpha*beta*dot_product(p,gradient)
        exit linesearch
    end if
    alpha = alpha/2
else if (counter > 10 .and. counter <= 20 ) then
    if (error2 < error + alpha*beta*dot_product(p,gradient) ) then
        !write(*,*) error2 - error - alpha*beta*dot_product(p,gradient)
        exit linesearch
    end if
    alpha = alpha/4

else if (counter > 20) then 
    exit linesearch

end if
counter = counter + 1
end do linesearch


solution = tempsol
error = error2
write(*,*) solution, error2

end do

end subroutine steepestdecent



   

