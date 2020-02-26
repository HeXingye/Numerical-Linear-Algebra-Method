program hs
  ! program2:QR algorithm without shift and with shift
  ! when A is  symmetric it become be tridiagonal
  ! compile:gfortran -fdefault-real-8 LinAl.f90 hs.f90 -o hs

  ! before compile, please input: make clean
  ! compile as: make
  ! or compile as:
  !gfortran -fdefault-real-8 LinAl.f90 hs.f90 -o hs
  !./hs Amat.dat	


  use LinAl
  real,dimension(:,:), allocatable :: A,As,Q,Qs
  real,dimension(:,:), allocatable :: Ide
  !real,dimension(:),allocatable :: x, b, bs, error, E 
  integer :: msize, nsize
  character*100 filenameA
  if(iargc().ne.1) then
     write(*,*) 'Wrong number of arguments. Please input the filds'
     stop
  endif
  call getarg(1,filenameA)
  call readmat(A,msize,nsize,filenameA)
  allocate(As(msize,nsize))
  As=A
  write(*,*) 'The matrix A is:'
  call writemat(A,msize,nsize)
write(*,*) '**********************************************'
  allocate(Q(msize,msize))
allocate(Qs(msize,msize))
Q=0.
Qs=0.
  !allocate(R(msize,nsize))
  call qrnoshift(A,Q,msize,nsize)
   write(*,*) 'QR of without shift of R:'
  call writemat(A,msize,nsize)
write(*,*) 'QR of without shift of Q:'
  call writemat(Q,msize,msize)

 
write(*,*) '**********************************************'
call qrshift(As,Qs,msize,nsize)
   write(*,*) 'QR of with shift of R:'
  call writemat(As,msize,nsize)
write(*,*) 'QR of with shift of Q:'
  call writemat(Qs,msize,msize)
 
 deallocate(A,As,Q,Qs)
end program hs
