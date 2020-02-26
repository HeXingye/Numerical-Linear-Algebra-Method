module LinAl
  implicit none

contains

  !********************************************************

  subroutine readmat(mat,msize,nsize,filename)
  ! This subroutine is for reading and allocating matrix
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

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)

  end subroutine readmat

  subroutine writemat(mat,msize,nsize)
    integer :: msize, nsize, i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat

  real function norm2(vec,v_size)
  !compute 2-norm 
    integer :: v_size, i
    real, dimension(v_size) :: vec
    norm2 = 0.
    do i=1,v_size
       norm2 = norm2 + vec(i)**2
       !write(*,*) norm2
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


 subroutine qrfactor(A,Q,msize,nsize)
    !do qr factorization
integer :: msize,nsize,i, j ,k
real, dimension(msize,nsize) :: A,R
real, dimension(msize,msize) :: Q,H,Ide,u
real, dimension(msize) :: v
real :: sum
real :: s
Ide=0.
! first, we need to initial the identity matrix
do i=1, msize
  do j=1, msize
   if(i.eq.j) then
    Ide(i,j)=1.0
   end if
  end do
end do

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
!update A. 
!A = A − 2*uj*uj^t*A=(I-2*uj*uj^t)*A=H*A
!H=Ide-2*matmul(u,ut)
v = v/norm2(v,msize)
H = Ide - 2.0*spread(v(1:msize),dim=2,ncopies=msize)*spread(v(1:msize),dim=1,ncopies=msize)
A=matmul(H,A)
Q=matmul(Q,H)
end do
end subroutine qrfactor

subroutine qrnoshift(A,Q,msize,nsize)
    !do qr factorization without shift
integer :: msize,nsize,i,j,k
real, dimension(msize,nsize) :: A,Ide,Z
real, dimension(msize,msize) :: Q

Ide=0.0
do i=1, msize
  do j=1, msize
   if(i.eq.j) then
    Ide(i,j)=1.0
   end if
  end do
end do
Q=Ide
do i=1,20
  Z=matmul(A,Q)
call qrfactor(Z,Q,msize,nsize)
end do
end subroutine qrnoshift

subroutine qrshift(As,Qs,msize,nsize)
    !do qr factorization with shift
integer :: msize,nsize,i,j,k
real, dimension(msize,nsize) :: As,Ide
real, dimension(msize,msize) :: Qs
real :: u,error
Ide=0.0
do i=1, msize
  do j=1, msize
   if(i.eq.j) then
    Ide(i,j)=1.0
   end if
  end do
end do
! Qs=Ide
error=100
do i=1,5
u=As(msize,msize)
As=As-u*Ide
call qrfactor(As,Qs,msize,nsize)
As=matmul(As,Qs)+u*Ide
end do
end subroutine qrshift

end module LinAl
