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
    ! mat = read matrix. Note that mat type is double precision.

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

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in rows
    integer :: msize, nsize, i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat


subroutine inversemat(A,B,n)
    !outputs matrix mat of dimension (msize,nize) in rows
    integer :: i,j,k,l,m,n,irow
    real, dimension(n,n) :: A,B
    real :: big, dum
!identity matrix B
do i=1,n
 do j=1,n
   if(i.eq.j) then
     B(i,j)=1.0
   else
     B(i,j)=0.0
   end if
end do
end do
 
do i=1,n
  big=A(i,i) !find pivot!
 do j=1,n
   if (abs(A(j,i)).gt.big) then
      big=abs(A(j,i))
      irow=j
   end if
end do
!exchange rows
if (big.gt.abs(A(i,i))) then
 do k=1,n
   dum=A(i,k)
   A(i,k)=A(irow,k)
   A(irow,k)= dum
   dum=B(i,k)
   B(i,k)=B(irow,k)
   B(irow,k)=dum
 end do
end if
!do operation to let rows below A(i,i) become zero
dum= A(i,i)
do j=1,n
   A(i,j)=A(i,j)/dum
   B(i,j)=B(i,j)/dum
end do

do j=i+1,n
 dum=A(j,i)
 do k=1,n
  A(j,k)=A(j,k)-dum*A(i,k)
  B(j,k)=B(j,k)-dum*B(i,k)
 end do
end do
end do

do i=1,n-1
 do j=i+1,n
  dum=A(i,j)
  do l=1,n
   A(i,l)=A(i,l)-dum*A(j,l)
    B(i,l)=B(i,l)-dum*B(j,l)
  end do
end do
end do
  
!end do
end subroutine inversemat

subroutine  iteration(A,x,lamb,n)
    !outputs matrix mat of dimension (msize,nize) in rows
    integer :: i,j,k,l,m,n,irow
    real, dimension(n,n) :: A,B,Ide,Bs
    real,dimension(n) :: y,x,r,z
    real :: u,nor,error,lamb
!three eigenvalues λ1 = −8.0286, λ2 = 7.9329, λ3 = 5.6689.
u=lamb
!identity matrix B
do i=1,n
 do j=1,n
   if(i.eq.j) then
     Ide(i,j)=1.0
   else
     Ide(i,j)=0.0
   end if
end do
end do
!initial x
do i=1,n
 if (i.eq.1) then 
x(i)=1.0
 else
x(i)=0.0
end if
end do
Bs=A-u*Ide
call inversemat(Bs,B,n) !B=Bs^-1=(A − µI)^−1
error=1.0
do while(error.gt.0.1)
r=matmul(B,x)
nor=norm2(r,n)
y=r/nor
z=y-x
x=y
error=norm2(z,n)
end do
end subroutine iteration
end module LinAl
