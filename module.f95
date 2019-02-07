module global
    implicit none

    real(8) , dimension(6,12)             :: slip
    real(8) , dimension(6,6)              :: Cel
    real(8) , dimension(:,:), Allocatable ::eul
    integer                               :: nlines
    real(8), dimension(3,3)               :: id
    real(8)                               ::  dt, pi 
                                          
    

    contains
    subroutine init()
    integer :: i
        pi = 4.D0*DATAN(1.D0)
        dt = 0.00001
    
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
     do i = 1,nlines
         read(15,*) eul(i,1:3)
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


end module