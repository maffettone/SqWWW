! LAST EDIT: Phil Maffettone 2016-04-15
module types

  ! The following contains global constants and paramters for use in programs
  ! It also contains some simple functions, which would be universally useful


  implicit none

  public
  ! Double Precision
  integer, parameter, public :: dp = kind(1.d0)
  ! Quadrouple Precision
  integer, parameter, public :: qp = selected_real_kind (32)

  ! Fundemental Constants
  ! Pi, double precision
  real(dp), parameter, public :: pi = acos(-1.d0)
  ! N_av, Avagadros number
  real(dp), parameter, public :: N_av = 6.0221413d23
  ! Boltzmans constant [J/K]
  real(dp), parameter, public :: k_b = 1.3806488d-23
  ! Planks constant, h/2pi [Js]
  real(dp), parameter, public :: hbar = 1.05459d-34
  ! Electron charge [C]
  real(dp), parameter, public :: electron = 1.60219d-19
  ! Dielectric permittivity of free space [F/m]
  real(dp), parameter, public :: epsilon0 = 8.8542d-12
  ! Imaginary unit
  complex :: imag = (0,1)


  ! Useful Classes
  type sym_op
     ! Symmetry operations for lattice points Rx+T
     real(dp), dimension (3,3) :: r                         !Rotation Matrix
     Real(dp), dimension (3) :: t                           !Translation Vector
  end type sym_op


contains

  ! --------------------Public Contents-----------------------------------------
  ! strlen(string_with_trailing_blanks)
  ! upper_case(string)
  ! lower_case(string)
  ! alphabetical(character)
  ! ----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------


  function strlen(st) result (r)
    ! Calculates number of occupied elements in a string with trailing blanks
    integer :: r, i                                         !Return and dummy
    character(len=*), intent(in) :: st                      !String
    character(len=1) :: c
    i=len(st)
    c=st(i:i)
    do while(c .eq. ' ')
       i = i-1
       c=st(i:i)
    end do
    r = i
  end function strlen
  ! ----------------------------------------------------------------------------

  elemental subroutine lower_case(word)
    ! convert entire word to lower case
    character (len=*) , intent(in out) :: word
    integer                            :: i,ic,nlen
    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
    end do
  end subroutine lower_case

  elemental subroutine upper_case(word)
    ! convert a word to upper case
    character (len=*), intent (in out) :: word
    integer :: i, ic, nlen
    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= iachar("a") .and. ic <=iachar("z")) word(i:i) = char(ic-32)
    end do
  end subroutine upper_case
  ! ---------------------------------------------------------------------------


  logical function alphabetical(a) result(r)
    ! Determines if a character is a memeber of the alphabet
    character(len=1), intent(in) :: a
    character(len=26) :: alpha, alpha2
    integer :: i
    alpha = 'abcdefghijklmnopqrstuvwxyz'
    alpha2 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    r= .false.
    do i=1,26
       if (a==alpha(i:i) .or. a==alpha2(i:i)) then
          r = .true.
       end if
    end do
  end function alphabetical
  ! -----------------------------------------------------------------------------


  function sym_op_(i) result(r)
    ! Collected Space groups as needed
    ! Constructor for sym_op object
    Integer, intent(in) :: i                                !International Space Group Number
    type(sym_op),allocatable :: r(:)                        !Symmetry operation
    type(sym_op),allocatable :: sym(:)                      !Temporary sym_op
    integer :: j                                            !Dummy variable

    SELECT case(i)
    case(19)
       ! sym19  P212121
       allocate(sym(4))
       sym(1)%r = reshape(dble((/1.,0.,0.,0.,1.,0.,0.,0.,1./)),(/3,3/))
       sym(1)%t = (/0.,0.,0./)
       sym(2)%r = reshape(dble((/1.,0.,0.,0.,-1.,0.,0.,0.,-1./)),(/3,3/))
       sym(2)%t = (/0.5,0.5,0./)
       sym(3)%r = reshape(dble((/-1.,0.,0.,0.,1.,0.,0.,0.,-1./)),(/3,3/))
       sym(3)%t = (/0.,0.5,0.5/)
       sym(4)%r = reshape(dble((/-1.,0.,0.,0.,-1.,0.,0.,0.,1./)),(/3,3/))
       sym(4)%t = (/0.5,0.,0.5/)
    case(198)
       allocate(sym(12))
       !  sym 198 P212121
       ! --------
       sym(1)%r = reshape(dble((/1.,0.,0.,0.,1.,0.,0.,0.,1./)),(/3,3/))
       sym(1)%t = (/0.,0.,0./)
       sym(2)%r = reshape(dble((/1.,0.,0.,0.,-1.,0.,0.,0.,-1./)),(/3,3/))
       sym(2)%t = (/0.5,0.5,0./)
       sym(3)%r = reshape(dble((/-1.,0.,0.,0.,1.,0.,0.,0.,-1./)),(/3,3/))
       sym(3)%t = (/0.,0.5,0.5/)
       sym(4)%r = reshape(dble((/-1.,0.,0.,0.,-1.,0.,0.,0.,1./)),(/3,3/))
       sym(4)%t = (/0.5,0.,0.5/)
       sym(5)%r = reshape(dble((/0.,1.,0.,0.,0.,1.,1.,0.,0./)),(/3,3/))
       sym(5)%t = (/0.,0.,0./)
       sym(6)%r = reshape(dble((/0.,-1.,0.,0.,0.,1.,-1.,0.,0./)),(/3,3/))
       sym(6)%t = (/0.5,0.,0.5/)
       sym(7)%r = reshape(dble((/0.,-1.,0.,0.,0.,-1.,1.,0.,0./)),(/3,3/))
       sym(7)%t = (/0.5,0.5,0./)
       sym(8)%r = reshape(dble((/0.,1.,0.,0.,0.,-1.,-1.,0.,0./)),(/3,3/))
       sym(8)%t = (/0.,0.5,0.5/)
       sym(9)%r = reshape(dble((/0.,0.,1.,1.,0.,0.,0.,1.,0./)),(/3,3/))
       sym(9)%t = (/0.,0.,0./)
       sym(10)%r = reshape(dble((/0.,0.,-1.,-1.,0.,0.,0.,1.,0./)),(/3,3/))
       sym(10)%t = (/0.,0.5,0.5/)
       sym(11)%r = reshape(dble((/0.,0.,1.,-1.,0.,0.,0.,-1.,0./)),(/3,3/))
       sym(11)%t = (/0.5,0.,0.5/)
       sym(12)%r = reshape(dble((/0.,0.,-1.,1.,0.,0.,0.,-1.,0./)),(/3,3/))
       sym(12)%t = (/0.5,0.5,0./)
    case(205)
       allocate(sym(24))
       ! sym 205 Pa-3
       ! ----------
       sym(1)%r = reshape(dble((/1.,0.,0.,0.,1.,0.,0.,0.,1./)),(/3,3/))
       sym(1)%t = (/0.,0.,0./)
       sym(2)%r = reshape(dble((/1.,0.,0.,0.,-1.,0.,0.,0.,-1./)),(/3,3/))
       sym(2)%t = (/0.5,0.5,0./)
       sym(3)%r = reshape(dble((/-1.,0.,0.,0.,1.,0.,0.,0.,-1./)),(/3,3/))
       sym(3)%t = (/0.,0.5,0.5/)
       sym(4)%r = reshape(dble((/-1.,0.,0.,0.,-1.,0.,0.,0.,1./)),(/3,3/))
       sym(4)%t = (/0.5,0.,0.5/)
       sym(5)%r = reshape(dble((/0.,0.,1.,1.,0.,0.,0.,1.,0./)),(/3,3/))
       sym(5)%t = (/0.,0.,0./)
       sym(6)%r = reshape(dble((/0.,0.,-1.,-1.,0.,0.,0.,1.,0./)),(/3,3/))
       sym(6)%t = (/0.,0.5,0.5/)
       sym(7)%r = reshape(dble((/0.,0.,1.,-1.,0.,0.,0.,-1.,0./)),(/3,3/))
       sym(7)%t = (/0.5,0.,0.5/)
       sym(8)%r = reshape(dble((/0.,0.,-1.,1.,0.,0.,0.,-1.,0./)),(/3,3/))
       sym(8)%t = (/0.5,0.5,0./)
       sym(9)%r = reshape(dble((/0.,1.,0.,0.,0.,1.,1.,0.,0./)),(/3,3/))
       sym(9)%t = (/0.,0.,0./)
       sym(10)%r = reshape(dble((/0.,-1.,0.,0.,0.,1.,-1.,0.,0./)),(/3,3/))
       sym(10)%t = (/0.5,0.,0.5/)
       sym(11)%r = reshape(dble((/0.,-1.,0.,0.,0.,-1.,1.,0.,0./)),(/3,3/))
       sym(11)%t = (/0.5,0.5,0./)
       sym(12)%r = reshape(dble((/0.,1.,0.,0.,0.,-1.,-1.,0.,0./)),(/3,3/))
       sym(12)%t = (/0.,0.5,0.5/)
       do j=1,12
          sym(j+12)%r = sym(j)%r*(-1)
          sym(j+12)%t = sym(j)%t
       end do
    case default
       write(*,*) 'Space group', i, 'is not in the types library'
       stop
    end select
    r=sym
  end function sym_op_
  ! ----------------------------------------------------------------------------



end module types
