program qr
  ! program2:QR Factorization Using Householder 
  ! compile:gfortran LinAl.f90 qr.f90 -o qr.exe
  ! before compile, please input: make clean
  ! compile as: make
  ! or compile as:
  ! gfortran -fdefault-real-8 LinAl.f90 qr.f90 -o qr
  !./qr atkinson.dat	

  use LinAl

  real,dimension(:,:), allocatable :: A,As,Q,Qt,R,Pa
  real,dimension(:,:), allocatable :: Ide
  real,dimension(:),allocatable :: x, b, bs, error, E ! vectors to store the x and y values to solve the linear problem Ax=b, and find error in solution
  logical :: ising=.false.
  integer :: msize, nsize, i, j
  real :: err, f ! to store error between f value at given x and computed solution value
  character*100 filename,outfilename
! Tis program read x,b from fiel atkinson.dat
! and produce Vandermonde matrices according to 
! the order user input
! then do QR factorization
! the error is given by |A*x-b|

if(iargc().ne.1) then 
  write(*,*) 'Wrong number of arguments. Please input the fild atkinson'
  stop
  endif
  call getarg(1,filename)

write(*,*) 'This is a 3-order polynomial' 
  nsize=4

! count how many lines of x,y
call countlines(filename,msize) 
  allocate(A(msize,nsize))
  allocate(As(msize,nsize))
  allocate(Ide(msize,msize))
call redvectors(x,b,msize,filename)
!read the vector of x,b from dat
write(*,*) 'There are',msize,'columns'
write(*,*) '************************************************' 

!produce the vandermonde matrix A according to polynomial order
call vandermonde(A,x,msize,nsize)
! Copy A,b into As and bs
  As = A
  bs = b 
  allocate(Q(msize,msize))
allocate(R(nsize,nsize))
  call qrfactor(A,Q,msize,nsize)
  !write(*,*)'Q is given by:'
  !call writemat(Q,msize,msize)
  !output A in triangular form
write(*,*) '****************QR Factorization****************' 
write(*,*) '************************************************' 
  write(*,*) 'Factorization A=QR:'
  call writemat(A,msize,nsize)
  write(*,*) 'Factorization A=QR with part R:'
call Rpart(A,R,msize,nsize)
call writemat(R,nsize,nsize)
write(*,*) '************************************************'
write(*,*) '****************compute A-QR********************'
 write(*,*) 'The matrix A-QR is:'
 call writemat(As-matmul(Q,A),msize,nsize)
write(*,*)'2-norm of A-QR=',normF(As-matmul(Q,A),msize,nsize)
  !Construct Qt containing first n columns of Q
  allocate(Qt(msize,nsize))
  do i=1,msize
     Qt(i,1:nsize)=Q(i,1:nsize)
  end do
  !Construct R containing first n rows of A

  !allocate(R(nsize,nsize))
  do i=1,nsize
     R(i,1:nsize)=A(i,1:nsize)
  end do
write(*,*) '************************************************'
write(*,*) '****************compute Q^TQ-I********************'
call imatrix(Ide,msize)
write(*,*) 'The matrix Q^T*Q-I is:'
call writemat(matmul(transpose(Q),Q)-Ide,msize,msize)
write(*,*)'2-norm of (Q^T*Q-I) =',normF(matmul(transpose(Q),Q)-Ide,msize,msize)

write(*,*) '************************************************'
write(*,*) '****************Least square equation********************'
allocate(pa(msize,msize))
  Pa=matmul(transpose(Q),Q)
write(*,*) 'The project matrix is:'
call writemat(pa,msize,msize)

  !solution of Rx=Q^T*b by backsustitution
  b = matmul(transpose(Qt),b)
write(*,*) '************************************************'
write(*,*) '****************solution vector x********************'
  write(*,*) ' The vector before solution'
  write(*,*) (b(i), i=1,nsize)
  call backsub(R,b,nsize)

  !output the solution onto screen
  write(*,*) 'The solution vector x is:'
  !call writemat(B,msize,nsize)
  write(*,*) (b(i), i=1,nsize)
  ! compute error matrix
  allocate(E(nsize))
  E = bs - matmul(As,b)
  write(*,*)'2-norm of error of solution is',norm2(E,nsize)
  !Error^2=sum(f(x(i),bs)-b(i))^2
  err = 0.0
call fitcurve(As,bs,b,x,err,msize,nsize)
  write(*,*)'the rms error between fitted curve and data is:', sqrt(err/msize)
  write(*,*)'Frobenius norm of Q^T*Q - I =',normF(matmul(transpose(Q),Q)-Ide,msize,msize)
 
deallocate(A,As,Q,Qt,R,b,Ide,E)

end program qr
