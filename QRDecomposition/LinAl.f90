module LinAl
  implicit none

contains

  !********************************************************

  subroutine readmat(mat,msize,nsize,filename)
    character*100 filename
    real, dimension(:,:), allocatable, intent(in out) :: mat
    integer :: msize, nsize
    integer :: i,j
    ! Reads a file containing the matrix A 
    ! Sample file:
    ! 3  4
    ! 1.2     1.3     -1.4    3.31
    ! 31.1    0.1     5.411   -1.23
    ! -5.4    7.42    10      -17.4
    ! Note that the first 2 lines are the matrix dimensions, 
    ! then the next msize lines are the matrix entries
    ! Note that entries must be separated by a tab.
    ! Then allocates an array of size (msize,nsize), populates the matrix,
    ! and returns the array. 

    ! This routine takes as INPUTS:
    ! filename = a character string with the name of file to read
    ! This routine returns as OUTPUTS:
    ! msize = first dimension of read matrix
    ! nsize = second dimension of read matrix
    ! mat = read matrix. Note that mat type is real.
    open(10,file=filename)
    ! Read the matrix dimensions
    read(10,*) msize,nsize
    ! Allocate matrix
    allocate(mat(msize,nsize))
    mat=0.0
    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo
    close(10)
  end subroutine readmat

  subroutine redvectors(mat1,mat2,msize,filename)
! Reads a vector x,b
    character*100 filename
    real, dimension(:), allocatable, intent(in out) :: mat1,mat2
    integer :: msize
    integer :: i
  allocate(mat1(msize))
allocate(mat2(msize))
  mat1=0.0
 mat2=0.0
     open(10,file=filename)
   do i = 1,msize 
     read(10,*) mat1(i),mat2(i)
  enddo
  close(10)
  end subroutine redvectors


  subroutine countlines(filename,msize)
!count how many lines in dat fiels,this is m size
    character*100 filename
    character (len=100) :: temp
    integer :: msize,i,io
  i = 0
  io = 0
  open(10,file=filename)
  do while(io.eq.0)  !if reading is correct
     i = i + 1
     read(10,*,iostat=io) temp,temp
  enddo
  close(10)
  msize = i-1
  end subroutine countlines

  subroutine vandermonde(A,x,msize,nsize)
    real,dimension(:,:), allocatable :: A
    real,dimension(:),allocatable :: x
    integer :: msize,nsize,i,j
  do i=1,msize
     do j = 1,nsize
        A(i,j) = x(i)**(j-1)
     enddo
  enddo
  end subroutine vandermonde

  subroutine fitcurve(As,bs,b,x,err,msize,nsize)
!what the 2-norm error on the fit is.
!f is the polynomial value
    real :: err, f
    real,dimension(:,:), allocatable :: A,As
    real,dimension(:),allocatable :: bs, error,b,x
    integer :: msize,nsize,i,j
allocate(error(msize))
 open(20,file='fitcurve.dat')
  do i=1,msize
     f = 0.0
     do j = 1, nsize
        f = f + As(i,j)*b(j)
     end do
     write(20,*) x(i), f
     error(i) = bs(i) - f
     err = err + error(i)**2
  enddo
  close(20)
  end subroutine fitcurve

!calculate norms of the columns of vector Error(m*n)
subroutine colnorm2(mat,msize,nsize,colnom) 
real,dimension(:,:), allocatable :: mat
real,dimension(:), allocatable :: colnom
integer :: msize, nsize, i, j
do j=1,nsize
     do i=1,msize
        colnom(i)= mat(i,j)
     end do
     write(*,*) 'The norm of column',j,'is',norm2(colnom,msize)
  end do
  end subroutine colnorm2

real function norm2(vec,sizev)
    !vector vec and its size sizev
    !compute the 2-orm of vec
    integer :: sizev, i
    real, dimension(sizev) :: vec
    norm2 = 0.
    do i=1,sizev
       norm2 = vec(i)**2+norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

real function normF(A,msize,nsize)
    !function to compute the Frobenius norm of a matrix of size dimensions (msize,nsize) 
    integer :: msize, nsize,i,j
    real, dimension(msize,nsize) :: A
    normF = 0.
    do i=1,msize
       do j=1,nsize
          normF = normF + A(i,j)**2
          !write(*,*) norm2
       enddo
    enddo
    normF = sqrt(normF)
  end function normF

subroutine imatrix(Ide,msize)
    !function to compute the Frobenius norm of a matrix of size dimensions (msize,nsize) 
    integer :: msize,i,j
    real, dimension(msize,msize) :: Ide
    do j=1,msize
     do i=1,msize
      if(i==j) then
        Ide(i,j) = 1    
       endif
     end do
  end do
  end subroutine imatrix

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in rows
    integer :: msize,nsize,i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat


    !A=QR part R
subroutine Rpart(A,R,msize,nsize)
real, dimension(msize,nsize) :: A
real, dimension(nsize,nsize) :: R
integer :: msize,nsize,i,j
R(nsize,nsize)=0.0
    do j=nsize,1,-1
       do i=1,j,1
          R(i,j)=R(i,j)+A(i,j)
       end do
  end do
end subroutine Rpart


subroutine qrfactor(A,Q,msize,nsize)
    !do qr factorization
integer :: msize,nsize,i, j ,k
real, dimension(msize,nsize) :: A
real, dimension(msize,msize) :: Q,H,Ide,u
real, dimension(msize) :: v !,u,ut
real :: sum
real :: s

! first, we need to initial the identity matrix
!do i=1, msize
  !do j=1, msize
   !if(i.eq.j) then
   ! Ide(i,j)=1.0
  ! end if
 ! end do
!end do
call imatrix(Ide,msize)
Q=Ide

!perform a housholder transformation of QR
!H=I-v^Tv=I-2*(uu^T)/(u^Tu) where u=v/|v|
!Q=H1H2.....Hm-1
!A is an m*n matrix

do i=1,nsize
sum=0.0
v=0.0
do j=i,msize
  sum=sum+A(j,i)**2
end do
! compute sign(a11)*|a1|*e1
!This is to find sing(a11)
s=sign(sqrt(sum),A(i,i)) !this equal to sign(a11)*||a1||
do k=i,msize
  if(k.eq.i) then
    v(k)=A(i,i)+s
!vk = [0, · · · 0, ajj + sj , aj+1,j , · · · , amj ]
  else
    v(k)=A(k,i)
  endif
end do
!compute u=v/|v| where |v| is norm if v
!u=v/norm2(v,msize)
!ut=transpose(u)
!update A. A = A − 2*uj*uj^t*A=(I-2*uj*uj^t)*A=H*A
!H=Ide-2*matmul(u,ut)
v = v/norm2(v,msize)
H = Ide - 2.0*spread(v(1:msize),dim=2,ncopies=msize)*spread(v(1:msize),dim=1,ncopies=msize)
A=matmul(H,A)
Q=matmul(Q,H)
end do
end subroutine qrfactor


  !subroutine backsub(A,Q,B,msize,nsize)
subroutine backsub(A,b,nsize)
!compute y=Q^T*b=Rx
!SOLVE Rx=y BY back substitution
!backsub Rx=Q^T*b
    integer :: i,j,k,nsize !,msize
    !real, dimension(msize,msize) :: Q
    real, dimension(nsize,nsize) :: A 
    real, dimension(nsize) :: b, x
    real :: sum
    do i=nsize,1,-1
       sum = 0.     
       do j=i+1,nsize
          sum = sum + A(i,j)*b(j)
       end do
       b(i) = (b(i) - sum)/A(i,i)
    end do
end subroutine backsub


end module LinAl
