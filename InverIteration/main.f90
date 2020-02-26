program inviteration
  ! program3:inverse iteration to calculate the corresponding eigenvectors
  ! compile:gfortran -fdefault-real-8 LinAl.f90 inviteration.f90-o iter

  ! before compile, please input: make clean
  ! compile as: make
	


  use LinAl
  real,dimension(:,:), allocatable :: A,As,B,Bs,A1 !,A2,A3
  real,dimension(:), allocatable :: x1 !,x2,x3
  !real,dimension(:,:), allocatable :: Ide
  !real,dimension(:),allocatable :: x, b, bs, error, E 
  integer :: msize, nsize
  real :: lamb !,lamb2,lamb3
  character*100 filenameA
  if(iargc().ne.1) then
     write(*,*) 'Wrong number of arguments. Please input the filds'
     stop
  endif
  call getarg(1,filenameA)
  call readmat(A,msize,nsize,filenameA)
write(*,*) 'please input estimated eigvalues(please make clean every time for use):' 
  read(*,*) lamb
  allocate(As(msize,nsize))
allocate(A1(msize,nsize))
!allocate(A2(msize,nsize))
!allocate(A3(msize,nsize))
allocate(x1(msize))
!allocate(x2(msize))
!allocate(x3(msize))
  As=A
  A1=A 
!A2=A
!A3=A
  write(*,*) 'The matrix A is:'
  call writemat(As,msize,nsize)
  allocate(B(msize,msize))
  allocate(Bs(msize,msize))
  call inversemat(As,Bs,msize)
write(*,*) '*********************************'
   write(*,*) 'The inverse matrix of A:'
  call writemat(Bs,msize,msize)
write(*,*) '*************please be patient,it need some time to run to get result********************'
call iteration(A,x1,lamb,msize)
write(*,*) 'eigenvector corresponding to estimated eigvalue you input:'
call writemat(x1,msize,1)
deallocate(A,As,B,Bs,x1,A1)
end program inviteration
