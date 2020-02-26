program hs
  ! program1:Reduction to Hessenberg form of a square matrix A
  ! when A is  symmetric it become be tridiagonal
  ! compile:gfortran -fdefault-real-8 LinAl.f90 hs.f90 -o hs

  ! before compile, please input: make clean
  ! compile as: make
  ! or compile as:
  ! gfortran -fdefault-real-8 LinAl.f90 hs.f90 -o hs
  !./hs Amat.dat	


  use LinAl
  real,dimension(:,:), allocatable :: A,As
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
  !allocate(Q(msize,msize))
!allocate(R(nsize,nsize))
  call qrfactor(A,msize,nsize)
   write(*,*) 'Reduction to Hessenberg form of a symetric matrix A:'
  call writemat(A,msize,nsize)
  deallocate(A,As)
end program hs
