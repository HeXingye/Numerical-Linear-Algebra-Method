module LinAl
   implicit none
 
 contains
 
   !********************************************************
 
   subroutine readandallocatemat(mat,msize,nsize,filename)
 
     character*100 filename
     real, dimension(:,:), allocatable, intent(in out) :: mat
     integer :: msize, nsize
 
     integer :: i,j
 
     ! Reads a file containing the matrix mat 
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
 
   end subroutine readandallocatemat
 
 
   real function trace(mat,msize)
      !function to calculate trace of a square matrix mat of size msize
     integer :: msize, i
     real, dimension(msize,msize) :: mat
     trace = 0.
     do i=1,msize
        trace = trace + mat(i,i)
     enddo
   end function trace
 
   real function norm2(vec,vsize)
   !function to compute the (Euclidean)norm of a vector vec of size vsize  !outputs 
     integer :: vsize, i
     real, dimension(vsize) :: vec
     norm2 = 0.
     do i=1,vsize
        norm2 = norm2 + vec(i)**2
        !write(*,*) norm2
     enddo
     norm2 = sqrt(norm2)
   end function norm2
 
   subroutine writemat(mat,msize,nsize)
     !outputs matrix mat of dimension (msize,nize) in row-wise split format
     integer :: msize, nsize, i, j
     real, dimension(msize,nsize) :: mat
     do i=1,msize
        write(*,*) ( mat(i,j), j=1,nsize )
     enddo
   end subroutine writemat
 
   subroutine gausse(A,B,msize,nsize,ising)
 
      integer :: msize, nsize
      real, dimension(msize,msize) :: A
      real, dimension(msize,nsize) :: B
      real, dimension(msize) :: scale
  
      logical :: ising
      integer :: i,j,k,l       ! locally defined variables
      real :: p,temp 
 
     !***** Performs  Gaussian Elimination with partial pivoting*********! 
     ! **** Inputs are:
     ! A = matrix to be pivoted, dimensions are A(msize,msize)
     ! B = matrix of RHS vectors, dimensions are B(msize,nsize) 
     ! ising = logical flag for checking if A is singular
 
     ! loop over rows to find largest element, store as the scale
 
     ising=.FALSE.
 ! Gaussian Elimiaton with pivoting
        do k=1, msize
            do i=k,msize
                scale(i) = 0.0
                do j=k,msize
                    scale(i) = max(scale(i),abs(A(i,j)))
                end do
            end do
            p = abs(A(k,k)/scale(k)) ! find largest p
            l = k
            do j=k+1,msize
                if((A(j,k)/scale(j)) > p) then
                    p = (A(j,k)/scale(j))
                    l = j
                end if
            end do
            if(p .eq. 0.0) then ! Check whether it is sigular
                iSing = .true.
                return
            end if

            if (l .ne. k) then
                do j=k,msize        !exchange rows of A
                    temp = A(k,j)
                    A(k,j) = A(l,j)
                    A(l,j) = temp
                end do
                do j=k,nsize
                    temp = B(k,j) !exchange rows of B
                    B(k,j) = B(l,j)
                    B(l,j) = temp
                enddo
            end if

            ! gaussian elimation
     
            do i=k+1,msize
               if(A(i,k).ne.0.0) then
                temp= A(i,k)/A(k,k)
                do j = 1, nsize !j
                    B(i,j)=B(i,j)-temp*B(k,j)
                enddo
                A(i,k) = 0.0
            !make all elements below diagonal=zero
                do j=k+1,msize
                    A(i,j) = A(i,j)-temp*A(k,j)
                end do
               end if
            end do
        enddo

    return

 write(*,*) ' matrix is singular'
    ising=.TRUE.

return
   end subroutine gausse
 
   subroutine backsub(A,B,msize,nsize)
     real, dimension(msize,msize) :: A
     real, dimension(msize,nsize) :: B
     real :: sum
     integer :: i,j,k,msize,nsize
     do j=1,nsize
        do i=msize,1,-1
           sum = 0.
           if(A(i,i).eq.0.0)then
 1000         write(*,*)'singular matrix'
              goto  2000
           end if
           do k=i+1,msize
              sum = sum + A(i,k)*B(k,j)
           end do
           B(i,j) = (B(i,j) - sum)/A(i,i)
        end do
     end do
 
 2000 return
   end subroutine backsub
 
 end module LinAl
 