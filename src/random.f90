! LAST EDIT: Phil Maffettone 2016-04-15
module random
  ! Module for Psuedo-Random-Number-Generation
  use types

  implicit none

contains
  ! --------------------Public Contents-----------------------------------------
  ! init_random_seed()
  ! rand_norm_dist(radom[0,1])
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  !
  ! ----------------------------------------------------------------------------

  subroutine init_random_seed()
    ! Utlization of hardware to generate seed for non-redundant random numbers
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))

    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)

  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed
  ! ----------------------------------------------------------------------------


  subroutine rand_normal_dist(r)
    ! Takes a random number [0,1] and returns a rand on the standard normal dist
    ! Follows Box-Muller method
    real(dp), intent(inout) :: r
    real(dp), save :: u
    real(dp) :: v
    logical,save :: FirstCall = .true.

    if (FirstCall) then
       FirstCall = .false.
       call random_number(u)
    end if
    v=u
    u=r
    r=sqrt(-2.0_dp*log(u))*cos(2*pi*v)
  end subroutine rand_normal_dist
  ! ----------------------------------------------------------------------------

 
end module random
