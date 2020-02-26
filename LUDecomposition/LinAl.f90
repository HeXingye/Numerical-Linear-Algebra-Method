module LinAl
  implicit none

contains

  !********************************************************

  subroutine readmat(mat,msize,nsize,filename)

    character*100 filename
    double precision, dimension(:,:), allocatable, intent(in out) :: mat
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

!calculate trace of a matrix
character function trace(mat,msize,nsize)
    integer :: msize,nsize, i
    real :: tr
    double precision, dimension(msize,msize) :: mat   
   if(nsize.ne.msize) then 
     write(*,*) 'cannot calculate trace, it is not a square matrix'
     return
   else
       tr=0.
       do i=1,msize
       tr = tr + mat(i,i) 
   enddo   
  end if
    write(*,*) tr
end function trace


!calculate norms of the columns of vector Error(m*n)
  subroutine colnorm2(mat,msize,nsize,colnom) 
double precision,dimension(:,:), allocatable :: mat
double precision,dimension(:), allocatable :: colnom
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
    double precision, dimension(sizev) :: vec
    norm2 = 0.
    do i=1,sizev
       norm2 = vec(i)**2+norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in rows
    integer :: msize, nsize, i, j
    double precision, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat



    !A=LU part L
subroutine Lpart(A,L,msize)
double precision, dimension(msize,msize) :: A
double precision, dimension(msize,msize) :: L
integer :: msize,i,j
L(msize,msize)=0.0
    do i=1,msize
       do j=1,i
       if (i==j) then
          L(i,j)=1.0
       else
          L(i,j)=L(i,j)+A(i,j)
       end if
       end do
  end do
end subroutine Lpart

    !A=LU part U
subroutine Upart(A,U,msize)
double precision, dimension(msize,msize) :: A
double precision, dimension(msize,msize) :: U
integer :: msize,i,j
U(msize,msize)=0.0
    do j=msize,1,-1
       do i=1,j,1
          U(i,j)=U(i,j)+A(i,j)
       end do
  end do
end subroutine Upart


  subroutine lufactor(A,msize,s,ising)
    !Factors the matrix A of dimension (msize,msize) into L,U 
    !L is the Lower triangular
    !U is the upper triangular
    !L-I and U are stored in A

    integer :: msize,i,j,k,colmax,srow,index,temp
    double precision, dimension(msize,msize) :: A
    double precision :: p,dum
    integer, dimension(msize) :: s
    logical :: ising

    do j=1,msize
       s(j)=j
    enddo
    !Loop over the columns
    do j=1,msize
       !Find the pivot row index and pivot p=max|A(k,j)|
       p=0.
       do k=j,msize
          srow=dabs(A(k,j))
          if(p.lt.srow)then
             p=srow
             index=k
          end if
       enddo
       !exchange rows if index is not same as j
       if(index.ne.j) then
          do k=1,msize   ! swap rows in A 
             dum = A(index,k)
             A(index,k) = A(j,k)
             A(j,k) = dum
          enddo
          !exchange entries of s
          temp=s(j)
          s(j)=index
          s(index)=temp
       endif
       !check for zero entry and stop with flag
        if(A(j,j).eq.0.0) then
          write(*,*) 'A is singular'
          ising = .TRUE.
          goto 10
       endif
       !Calculate the i-th column of L
       do i=j+1,msize
          A(i,j)= A(i,j)/A(j,j)
          do k=j+1,msize
             A(i,k)=A(i,k) - A(i,j)*A(j,k)
          end do
       end do
    end do
10 continue

  end subroutine lufactor

  subroutine backsub(A,B,msize,nsize,s)
    double precision, dimension(msize,msize) :: A
    double precision, dimension(msize,nsize) :: B,Bperm
    integer, dimension(msize) :: s
    real :: sum
    integer :: i,j,k,l,msize,nsize


    !perfoms backsubstitution for Ax=B
    !we get new B=PB and store them in B
    !iteration
    do j=1,nsize
       l=0
       do i=1,msize
          Bperm(i,j) = B(s(i),j)
          !write(*,*) s(i)
       end do
       B(:,j) = Bperm(:,j)

       !do L^-1*P*B and store in B 
       do i=1,msize
          sum = B(i,j)
          if (l.ne.0) then
             do k=l,i-1
                sum=sum-A(i,k)*B(k,j) 
             enddo
          else if (sum.ne.0.) then
             l=i
          endif
          B(i,j)=sum         
       enddo

       !backsubtitution
       do i=msize,1,-1
          sum = 0.0
          if(A(i,i).eq.0.0)then
10         write(*,*)'singular matrix'
             goto  20
          end if
          do k=i+1,msize
             sum = sum + A(i,k)*B(k,j)
          end do
          B(i,j) = -(sum-B(i,j))/A(i,i)
       end do
    end do

20 return
  end subroutine backsub


end module LinAl
