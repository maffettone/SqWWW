! LAST EDIT: Phil Maffettone 2016-05-01
module maths
  ! Useful mathematics functions for double precision reals
  ! Some complex functionality is built in
  ! Depends on LAPACK library, http://www.netlib.org/lapack/explore-html/index.html
  ! Compile with -llapack

  use types

  implicit none

  public

  interface det
     module procedure det_cmplx,det_real
  end interface det

contains
  ! --------------------------Public Contents-----------------------------------
  ! FUNCTION: inv(A)
  ! FUNCTION: eigenvalues(A)
  ! FUNCTION: det(A)
  ! FUNCTION: cross(A,B)
  ! SUBROUTINE: RotMatrix(angle,vector,matrix)
  ! FUNCTION: factorial(n)
  ! FUNCTION: linear_interp(x,y,x_targ)
  ! FUNCTION: bin_data(x,y,n_bin)
  ! FUNCTION: convolve(x,h)
  ! ----------------------------------------------------------------------------

  function inv(A) result(Ainv)
    ! Returns the inverse of a matrix calculated by finding the LU decomposition
    ! decomposition.
    real(dp), dimension(:,:), intent(in) :: A               !Matrix to be inverted
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv        !Returned inverse matrix
    real(dp), dimension(size(A,1)) :: work                  !Work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv                   !pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv
  ! ----------------------------------------------------------------------------

  function eigenvalues(A) result(w)
    ! Calculates Matrix Eigenvalues
    ! Real symetric
    ! Very simple, very specific
    implicit none
    real(dp), dimension(:,:), intent(in) :: A               !Input matrix
    real(dp), allocatable :: w(:)                           !Return array of eigenvalues
    character(1) :: JOBZ                                    !Informaiton for Lapack
    character(1) :: UPLO                                    !Info for Lapack
    integer :: N,LDA,info,lwork                             !Integers for lapack
    real(dp), allocatable :: work(:)                        !Work array

    info = 0
    JOBZ = 'N'
    UPLO = 'U'
    N = size(A(1,:))
    allocate(work(max(1,lwork)), w(N))
    LDA = N

    lwork = -1
    call DSYEV(jobz,uplo,n,A,lda,w,work,lwork,info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(max(1,lwork)))

    call DSYEV(jobz,uplo,n,A,lda,w,work,lwork,info)
    if (info /= 0) then
       write(*,*) 'Eigenvalue error'
    end if
  end function eigenvalues
  ! ----------------------------------------------------------------------------

  function det_cmplx(N,mat) result(r)
    ! Calculates Matrix Determinant by Guassian Elimination then diagonal mult
    ! Complex, then real
    implicit none

    integer, intent(in) :: N                                !Dimension of square matrix
    complex(dp), intent(in), dimension(:,:) :: mat          !Matrix for determinant
    complex(dp),allocatable:: m(:,:)                        !Local copy to be manipulated
    complex(dp) :: r                                        !Return value
    integer :: i, info                                      !Dummy integer and error flag
    integer, allocatable :: ipiv(:)                         !Pivot indecies for row interchange
    real(dp) :: sgn                                         !Value for sign of determinant

    allocate(ipiv(N), m(size(mat,1),size(mat,2)))
    ipiv = 0
    m=mat
    call zgetrf(N, N, m, N, ipiv, info)
    r=1
    do i=1,N
       r=r*m(i,i)
    end do
    sgn=1
    do i =1,N
       if(ipiv(i)/=i)then
          sgn=-sgn
       end if
    end do
    r=sgn*r
  end function det_cmplx
  function det_real(N,mat) result(r)
    implicit none
    integer, intent(in) :: N                                !Dimension of square matrix
    real(dp), intent(in), dimension(:,:) :: mat             !Matrix for determinant
    complex(dp),allocatable:: m(:,:)                        !Local copy to be manipulated
    complex(dp) :: det                                      !Return value
    integer :: i, info                                      !Dummy integer and error flag
    integer, allocatable :: ipiv(:)                         !Pivot indecies for row interchange
    real(dp) :: sgn,r                                       !Value for sign of determinant

    allocate(ipiv(N), m(size(mat,1),size(mat,2)))
    ipiv = 0
    m=cmplx(mat,kind=dp)
    call zgetrf(N, N, m, N, ipiv, info)
    det=1
    do i=1,N
       det=det*m(i,i)
    end do
    sgn=1
    do i =1,N
       if(ipiv(i)/=i)then
          sgn=-sgn
       end if
    end do
    r=real(sgn*det)
  end function det_real
  ! ---------------------------------------------------------------------------------

  function cross(a,b)
    ! Generic cross product function in double precision
    real(dp), dimension(3) :: cross
    real(dp), dimension(3), intent(in) :: a,b
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
  ! ----------------------------------------------------------------------------

  subroutine RotMatrix(angle,vector,matrix)
    ! Generate rotation matrix for <angle> about <vector>
    ! Uses Rodrigues rotation formula
    real(dp), intent(in) :: angle,vector(3)                 !Angle and axis
    real(dp), intent(out) :: matrix(3,3)                    !Rotaiton Matrix
    real(dp) :: uvec(3)                                     !Unit vector||axis
    real(dp) :: identity(3,3)                               !Identity matrix
    real(dp) :: W(3,3)                                      !Internal matrix for calculation

    uvec = vector/norm2(vector)
    identity =reshape(dble([1,0,0,0,1,0,0,0,1]),shape(identity))
    W = reshape([0._dp,uvec(3),-1*uvec(2),&
         -1*uvec(3),0._dp,uvec(1),&
         uvec(2),-1*uvec(1),0._dp],shape(W))
    matrix = identity + W*sin(angle)+matmul(W,W)*(1-cos(angle))
  end subroutine RotMatrix
  ! -------------------------------------------------------------------------------

  recursive function factorial(n) result(r)
    ! Calculates N!
    integer,intent(in) :: n                                 !Number N
    integer :: r                                            !N factorial
    if (n == 0) then
       r = 1
    else if (n>0) then
       r = n*factorial(n-1)
    else
       r=0
       write(*,"(A)") 'Negative number passed to factorial function. Zero retuned.'
    end if
  end function factorial
  ! ------------------------------------------------------------------------------

  function linear_interp(x,y,x0) result(r)
    ! Linearly interpolates a funciton y(x) for a value of y0 at the point x0
    ! Assumes the x,y data has been sorted
    real(dp), intent(in) :: x(:),y(:)                       !Discrete function y(x)
    real(dp) :: x0                                          !Desired point x0
    real(dp) :: r                                           !Desired value y(x0)
    integer :: i                                            !Dummy integer

    do i=1,size(x)
       if (x(i)<=x0 .and. x(i+1) >= x0) then
          r = y(i) + (y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i))
          exit
       end if
    end do
  end function linear_interp
  ! -------------------------------------------------------------------------------

  function bin_data(x,y,n_bin) result(r)
    ! Recasts discretized function on evenly spaced bins using linear interpolation
    ! Result is the same size as x,y data
    real(dp), intent(in) :: x(:),y(:)                       !Discrete function y(x)
    real(dp) :: r(size(x),2)                                !Return funciton
    integer, intent(in) :: n_bin                            !Number of desired bins
    real(dp) :: x_min                                       !Minimum x value
    real(dp) :: x_max                                       !Maximum x value
    real(dp) :: bin_size                                    !Bin width
    integer :: i                                            !Dummy integer

    r=0.
    x_max = maxval(x)
    x_min = minval(x)
    bin_size = (x_max-x_min)/n_bin

    do i=1,size(x)
       r(i,1) = x_min+bin_size/2.0_dp
       x_min = x_min + bin_size
       r(i,2) = linear_interp(x,y,r(i,1))
    end do
  end function bin_data
  ! -----------------------------------------------------------------------------

  function convolve(x,h) result(r)
    ! Convolves 1D function x with 1D kernel h
    ! Slow convolution algorithm without using an FFT

    real(dp), allocatable :: r(:)
    real(dp) :: x(:)                                        !Sgnal array, f(t)
    real(dp) :: h(:)                                        !Noise/impulse array, g(t)
    integer :: kernelsize                                   !Size of g(t) array
    integer :: datasize                                     !Size of f(t) array
    integer :: i,j,k                                        !Dummy integers

    datasize = size(x)
    kernelsize = size(h)

    allocate(r(datasize))

    ! Last part
    do i=kernelsize,datasize
       r(i) = 0.0_dp
       j = i
       do k=1,kernelsize
          r(i) = r(i) + x(j)*h(k)
          j = j-1
       end do
    end do

    ! First part, ignores the leading component of the convolution
    ! This way r is of the same size as x
    do i=1,kernelsize
       r(i) = 0.0_dp
       j=i
       k=1
       do while(j>0)
          r(i) = r(i) +x(j)*h(k)
          j = j-1
          k = k+1
       end do
    end do

  end function convolve
  !-----------------------------------------------------------------------------

end module maths
