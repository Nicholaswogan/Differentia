module forwarddiff_dual
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use forwarddiff_const, only: wp
  implicit none
  private
  
  public :: dual
  public :: initialize_dual_array

  public :: assignment (=)
  public :: operator (+), operator (-)
  public :: operator (*), operator (/)
  public :: operator (**)
  public :: operator (==)
  public :: operator (<=), operator (<), operator (>=), operator (>)
  public :: operator (/=)

  public :: abs, acos, asin, atan, atan2
  public :: cos, exp, int, log, log10
  public :: max, maxval, min, minval
  public :: sin, tan, sqrt, sum, maxloc
  
  type :: dual
    real(wp) :: val
    real(wp), allocatable :: der(:)
  end type
  interface dual
    module procedure :: create_dual1
  end interface

  interface initialize_dual_array
    module procedure :: initialize_dual_array_1
    module procedure :: initialize_dual_array_2
    module procedure :: initialize_dual_array_3
  end interface
      
  !~~~ operator overloading ~~~!

  interface assignment (=)
    module procedure :: assign_di
    module procedure :: assign_dr
    module procedure :: assign_id
  end interface

  interface operator (+)
    module procedure :: add_d
    module procedure :: add_dd
    module procedure :: add_di
    module procedure :: add_dr
    module procedure :: add_id
    module procedure :: add_rd
  end interface

  interface operator (-)
    module procedure :: minus_d  
    module procedure :: minus_dd
    module procedure :: minus_di
    module procedure :: minus_dr
    module procedure :: minus_id
    module procedure :: minus_rd
  end interface
  
  interface operator (*)
    module procedure :: mult_dd
    module procedure :: mult_di
    module procedure :: mult_dr
    module procedure :: mult_id
    module procedure :: mult_rd
  end interface
  
  interface operator (/)
    module procedure :: div_dd
    module procedure :: div_di
    module procedure :: div_dr
    module procedure :: div_id
    module procedure :: div_rd
  end interface

  interface operator (**)
    module procedure :: pow_i
    module procedure :: pow_dr
    module procedure :: pow_rd
    module procedure :: pow_dd
  end interface

  interface operator (==)
    module procedure :: eq_dd
    module procedure :: eq_di
    module procedure :: eq_dr
    module procedure :: eq_id
    module procedure :: eq_rd
  end interface

  interface operator (<=)
    module procedure :: le_dd
    module procedure :: le_di
    module procedure :: le_dr
    module procedure :: le_id
    module procedure :: le_rd
  end interface

  interface operator (<)
    module procedure :: lt_dd
    module procedure :: lt_di
    module procedure :: lt_dr
    module procedure :: lt_id
    module procedure :: lt_rd
  end interface

  interface operator (>=)
    module procedure :: ge_dd
    module procedure :: ge_di
    module procedure :: ge_dr
    module procedure :: ge_id
    module procedure :: ge_rd
  end interface

  interface operator (>)
    module procedure :: gt_dd
    module procedure :: gt_di
    module procedure :: gt_dr
    module procedure :: gt_id
    module procedure :: gt_rd
  end interface

  interface operator (/=)
    module procedure :: ne_dd
    module procedure :: ne_di
    module procedure :: ne_dr
    module procedure :: ne_id
    module procedure :: ne_rd
  end interface

  !~~~ overloading intrinsic functions ~~~!

  interface abs
    module procedure :: abs_d
  end interface
  
  interface acos
    module procedure :: acos_d
  end interface
  
  interface asin
    module procedure :: asin_d
  end interface
  
  interface atan
    module procedure :: atan_d
  end interface
  
  interface atan2
    module procedure :: atan2_d
  end interface
  
  interface cos
    module procedure :: cos_d
  end interface
  
  interface exp
    module procedure :: exp_d
  end interface
  
  interface int
    module procedure :: int_d
  end interface
  
  interface log
    module procedure :: log_d
  end interface
  
  interface log10
    module procedure :: log10_d
  end interface

  interface max
    module procedure :: max_dd
    module procedure :: max_di
    module procedure :: max_dr
    module procedure :: max_rd
  end interface
  
  interface maxval
    module procedure :: maxval_d
  end interface

  interface min
    module procedure :: min_dd
    module procedure :: min_dr
    module procedure :: min_rd
  end interface
  
  interface minval
    module procedure :: minval_d
  end interface
  
  interface sin
    module procedure :: sin_d
  end interface
  
  interface tan
    module procedure :: tan_d
  end interface
  
  interface sqrt
    module procedure :: sqrt_d
  end interface
  
  interface sum
    module procedure :: sum_d
  end interface
  
  interface maxloc
    module procedure :: maxloc_d
  end interface

contains

  !~~~ initialize ~~~!

  function create_dual1(ndv) result(b)
    integer, intent(in) :: ndv
    type(dual) :: b
    allocate(b%der(ndv))
  end function

  subroutine initialize_dual_array_1(arr, ndv)
    type(dual), intent(inout) :: arr(:)
    integer, intent(in) :: ndv
    integer :: i
    do i = 1,size(arr)
      allocate(arr(i)%der(ndv))
    enddo
  end subroutine

  subroutine initialize_dual_array_2(arr, ndv)
    type(dual), intent(inout) :: arr(:,:)
    integer, intent(in) :: ndv
    integer :: i, j
    do i = 1,size(arr,2)
      do j = 1,size(arr,1)
        allocate(arr(j,i)%der(ndv))
      enddo
    enddo
  end subroutine

  subroutine initialize_dual_array_3(arr, ndv)
    type(dual), intent(inout) :: arr(:,:,:)
    integer, intent(in) :: ndv
    integer :: i, j, k
    do i = 1,size(arr,3)
      do j = 1,size(arr,2)
        do k = 1,size(arr,1)
          allocate(arr(k,j,i)%der(ndv))
        enddo
      enddo
    enddo
  end subroutine

  !~~~ assignment ~~~!

  elemental subroutine assign_di(u, i)
    type(dual), intent(inout) :: u
    integer, intent(in) :: i
    if (.not. allocated(u%der)) then
      error stop "Uninitialized dual"
    endif
    u%val = real(i,wp)
    u%der = 0.0_wp
  end subroutine

  elemental subroutine assign_dr(u, r)
    type(dual), intent(inout) :: u
    real(wp), intent(in) :: r
    if (.not. allocated(u%der)) then
      error stop "Uninitialized dual"
    endif
    u%val = r
    u%der = 0.0_wp
  end subroutine

  elemental subroutine assign_id(i, v)
    type(dual), intent(in) :: v
    integer, intent(out) :: i
    i = int(v%val)
  end subroutine

  !~~~ addition ~~~!

  elemental function add_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res
    res = u
  end function

  elemental function add_dd(u, v) result(res)
    type(dual), intent(in) :: u
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = u%val + v%val
    res%der = u%der + v%der
  end function

  elemental function add_di(u, i) result(res)
    type(dual), intent(in) :: u
    integer, intent(in) :: i
    type(dual) :: res
    res%val = real(i,wp) + u%val
    res%der = u%der
  end function

  elemental function add_dr(u, r) result(res)
    type(dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res
    res%val = r + u%val
    res%der = u%der
  end function

  elemental function add_id(i, v) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = real(i,wp) + v%val
    res%der = v%der
  end function
  
  elemental function add_rd(r, v) result(res)
    real(wp), intent(in) :: r
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = r + v%val
    res%der = v%der
  end function

  !~~~ subtraction ~~~!

  elemental function minus_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res
    res%val = - u%val
    res%der = - u%der
  end function

  elemental function minus_dd(u, v) result(res)
    type(dual), intent(in) :: u, v
    type(dual) :: res
    res%val = u%val - v%val
    res%der = u%der - v%der
  end function

  elemental function minus_di(u, i) result(res)
    type(dual), intent(in) :: u
    integer, intent(in) :: i
    type(dual) :: res
    res%val = u%val - real(i,wp)
    res%der = u%der
  end function

  elemental function minus_dr(u, r) result(res)
    type (dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res
    res%val = u%val - r
    res%der = u%der
  end function

  elemental function minus_id(i, v) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = real(i,wp) - v%val
    res%der = -v%der
  end function

  elemental function minus_rd(r, v) result(res)
    real(wp), intent(in) :: r
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = r - v%val
    res%der = - v%der
  end function

  !~~~ multiplication ~~~!
  
  elemental function mult_dd(u, v) result(res)
    type(dual), intent(in) :: u
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = u%val*v%val
    res%der = u%val*v%der + v%val*u%der
  end function

  elemental function mult_di(u, i) result(res)
    type(dual), intent(in) :: u
    integer, intent(in) :: i
    type(dual) :: res
    real(wp) :: r
    r = real(i,wp)
    res%val = r * u%val
    res%der = r * u%der
  end function

  elemental function mult_dr(u, r) result(res)
    type(dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res
    res%val = u%val*r
    res%der = u%der*r
  end function

  elemental function mult_id(i, v) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: v
    type(dual) :: res
    real(wp) :: r
    r = real(i,wp)
    res%val = r * v%val
    res%der = r * v%der
  end function
  
  elemental function mult_rd(r, v) result(res)
    real(wp), intent(in) :: r
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = r*v%val
    res%der = r*v%der
  end function

  !~~~ division ~~~!
  
  elemental function div_dd(u, v) result(res)
    type(dual), intent(in) :: u
    type(dual), intent(in) :: v
    type(dual) :: res
    real(wp) :: inv
    inv = 1.0_wp/v%val
    res%val = u%val*inv
    res%der = (u%der - res%val*v%der)*inv
  end function

  elemental function div_di(u, i) result(res)
    type(dual), intent(in) :: u
    integer, intent(in) :: i
    type(dual) :: res
    real(wp) :: inv
    inv = 1.0_wp / real(i,wp)
    res%val = u%val * inv
    res%der = u%der * inv
  end function
  
  elemental function div_dr(u, r) result(res)
    type(dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res
    real(wp) :: inv
    inv = 1.0_wp/r
    res%val = u%val*inv
    res%der = u%der*inv
  end function

  elemental function div_id(i, v) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: v
    type(dual) :: res
    real(wp) :: inv
    inv = 1.0_wp/v%val
    res%val = real(i,wp) * inv
    res%der = -res%val * inv * v%der
  end function
  
  elemental function div_rd(r, v) result(res)
    real(wp), intent(in) :: r
    type(dual), intent(in) :: v
    type(dual) :: res
    real(wp) :: inv
    inv = 1.0_wp / v%val
    res%val = r * inv
    res%der = -res%val*inv*v%der
  end function

  !~~~ power ~~~!

  elemental function pow_i(u, i) result(res)
    type(dual), intent(in) :: u
    integer, intent(in) :: i
    type(dual) :: res
    real(wp) :: pow_x
    pow_x = u%val ** (i - 1)
    res%val = u%val * pow_x
    res%der = real(i,wp) * pow_x * u%der
  end function

  elemental function pow_dr(u, r) result(res)
    type(dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res
    real(wp) :: pow_x
    pow_x = u%val**(r - 1.0_wp)
    res%val = u%val*pow_x
    res%der = r*pow_x*u%der
  end function

  elemental function pow_rd(r, u) result(res)
    real(wp), intent(in) :: r
    type(dual), intent(in) :: u
    type(dual) :: res
    res%val = r**u%val
    res%der = log(r)*r**u%val*u%der
  end function

  elemental function pow_dd(u, v) result(res)
    type(dual), intent(in) :: u
    type(dual), intent(in) :: v
    type(dual) :: res
    res%val = u%val**v%val
    res%der = res%val * (v%val / u%val * u%der + log(u%val) * v%der)
  end function

  !~~~ comparison (==) ~~~!

  elemental function eq_dd(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%val == rhs%val)
  end function

  elemental function eq_di(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    integer, intent(in) :: rhs
    logical :: res
    res = (lhs%val == rhs)
  end function eq_di

  elemental function eq_dr(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    real(wp), intent(in) :: rhs
    logical :: res
    res = (lhs%val == rhs)
  end function

  elemental function eq_id(lhs, rhs) result(res)
    integer, intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs == rhs%val)
  end function eq_id

  elemental function eq_rd(lhs, rhs) result(res)
    real(wp), intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs == rhs%val)
  end function

  !~~~ comparison (<=) ~~~!

  elemental function le_dd(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%val <= rhs%val)
  end function le_dd

  elemental function le_di(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    integer, intent(in) :: rhs
    logical :: res
    res = (lhs%val <= rhs)
  end function le_di

  elemental function le_dr(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    real(wp), intent(in) :: rhs
    logical :: res
    res = (lhs%val <= rhs)
  end function le_dr

  elemental function le_id(i, rhs) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: rhs
    logical :: res
    res = (i <= rhs%val)
  end function le_id

  elemental function le_rd(lhs, rhs) result(res)
    real(wp), intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs <= rhs%val)
  end function le_rd

  !~~~ comparison (<) ~~~!

  elemental function lt_dd(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%val < rhs%val)
  end function

  elemental function lt_di(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    integer, intent(in) :: rhs
    logical :: res
    res = (lhs%val < rhs)
  end function

  elemental function lt_dr(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    real(wp), intent(in) :: rhs
    logical :: res
    res = (lhs%val < rhs)
  end function

  elemental function lt_id(i, rhs) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: rhs
    logical :: res
    res = (i < rhs%val)
  end function

  elemental function lt_rd(lhs, rhs) result(res)
    real(wp), intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs < rhs%val)
  end function

  !~~~ comparison (>=) ~~~!

  elemental function ge_dd(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%val >= rhs%val)
  end function

  elemental function ge_di(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    integer, intent(in) :: rhs
    logical :: res
    res = (lhs%val >= rhs)
  end function

  elemental function ge_dr(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    real(wp), intent(in) :: rhs
    logical :: res
    res = (lhs%val >= rhs)
  end function

  elemental function ge_id(i, rhs) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: rhs
    logical :: res
    res = (i >= rhs%val)
  end function

  elemental function ge_rd(lhs, rhs) result(res)
    real(wp), intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs >= rhs%val)
  end function

  !~~~ comparison (>) ~~~!

  elemental function gt_dd(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%val > rhs%val)
  end function

  elemental function gt_di(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    integer, intent(in) :: rhs
    logical :: res
    res = (lhs%val > rhs)
  end function

  elemental function gt_dr(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    real(wp), intent(in) :: rhs
    logical :: res
    res = (lhs%val > rhs)
  end function

  elemental function gt_id(i, rhs) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: rhs
    logical :: res
    res = (i > rhs%val)
  end function

  elemental function gt_rd(lhs, rhs) result(res)
    real(wp), intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs > rhs%val)
  end function

  !~~~ comparison (/=) ~~~!

  elemental function ne_dd(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%val /= rhs%val)
  end function

  elemental function ne_di(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    integer, intent(in) :: rhs
    logical :: res
    res = (lhs%val /= rhs)
  end function

  elemental function ne_dr(lhs, rhs) result(res)
    type(dual), intent(in) :: lhs
    real(wp), intent(in) :: rhs
    logical :: res
    res = (lhs%val /= rhs)
  end function

  elemental function ne_id(i, rhs) result(res)
    integer, intent(in) :: i
    type(dual), intent(in) :: rhs
    logical :: res
    res = (i /= rhs%val)
  end function

  elemental function ne_rd(lhs, rhs) result(res)
    real(wp), intent(in) :: lhs
    type(dual), intent(in) :: rhs
    logical :: res
    res = (lhs /= rhs%val)
  end function

  !~~~ overloading intrinsics ~~~!

  elemental function abs_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res
    integer :: i

    if(u%val > 0) then
      res%val = u%val
      res%der = u%der
    else if (u%val < 0) then
      res%val = -u%val
      res%der = -u%der
    else
      res%val = 0.0_wp
      allocate(res%der(size(u%der)))
      do i = 1, size(res%der)
        if (u%der(i) == 0.0_wp) then
          res%der(i) = 0.0_wp
        else
          res%der(i) = ieee_value(1.0_wp, ieee_quiet_nan)
        end if
      end do
    endif

  end function

  elemental function acos_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    res%val = acos(u%val)
    if (u%val == 1.0_wp .or. u%val == -1.0_wp) then
      allocate(res%der(size(u%der)))
      res%der = ieee_value(1.0_wp, ieee_quiet_nan)
    else
      res%der = -u%der / sqrt(1.0_wp - u%val**2)
    end if

  end function

  elemental function asin_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    res%val = asin(u%val)
    if (u%val == 1.0_wp .or. u%val == -1.0_wp) then
      allocate(res%der(size(u%der)))
      res%der = ieee_value(1.0_wp, ieee_quiet_nan)
    else
      res%der = u%der / sqrt(1.0_wp - u%val**2)
    end if

  end function

  elemental function atan_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    res%val = atan(u%val)
    res%der = u%der / (1.0_wp + u%val**2)

  end function

  elemental function atan2_d(u, v) result(res)
    type(dual), intent(in) :: u, v
    type(dual) :: res

    real(wp) :: usq_plus_vsq

    res%val = atan2(u%val, v%val)

    usq_plus_vsq = u%val**2 + v%val**2
    res%der = v%val / usq_plus_vsq * u%der - u%val / usq_plus_vsq * v%der

  end function

  elemental function cos_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    res%val = cos(u%val)
    res%der = -sin(u%val) * u%der

  end function

  elemental function exp_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    real(wp) :: exp_x

    exp_x = exp(u%val)
    res%val = exp_x
    res%der = u%der * exp_x

  end function

  elemental function int_d(u) result(res)
    type(dual), intent(in) :: u
    integer :: res

    res = int(u%val)

  end function

  elemental function log_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    real(wp) :: inv

    inv = 1.0_wp / u%val
    res%val = log(u%val)
    res%der = u%der * inv

  end function log_d

  elemental function log10_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    real(wp) :: inv

    inv = 1.0_wp / (u%val * log(10.0_wp))
    res%val = log10(u%val)
    res%der = u%der * inv

  end function

  elemental function max_dd(val1, val2, val3, val4, val5) result(res)
    type(dual), intent(in) :: val1, val2
    type(dual), intent(in), optional :: val3, val4, val5
    type(dual) :: res

    if (val1%val > val2%val) then
      res = val1
    else
      res = val2
    endif
    if (present(val3)) then
      if (res%val < val3%val) res = val3
    endif
    if (present(val4)) then
      if (res%val < val4%val) res = val4
    endif
    if(present(val5))then
      if (res%val < val5%val) res = val5
    endif

  end function

  elemental function max_di(u, i) result(res)
    type(dual), intent(in) :: u
    integer, intent(in) :: i
    type(dual) :: res

    if (u%val > i) then
      res = u
    else
      allocate(res%der(size(u%der)))
      res = i
    endif

  end function

  elemental function max_dr(u, r) result(res)
    type(dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res

    if (u%val > r) then
      res = u
    else
      allocate(res%der(size(u%der)))
      res = r
    endif

  end function

  elemental function max_rd(n, u) result(res)
    real(wp), intent(in) :: n
    type(dual), intent(in) :: u
    type(dual) :: res

    if (u%val > n) then
      res = u
    else
      allocate(res%der(size(u%der)))
      res = n
    endif

  end function

  function maxval_d(u) result(res)
    type(dual), intent(in) :: u(:)
    type(dual) :: res
    integer :: i

    i = maxloc(u%val, 1)
    res = u(i)

  end function

  elemental function min_dd(val1, val2, val3, val4) result(res)
    type(dual), intent(in) :: val1, val2
    type(dual), intent(in), optional :: val3, val4
    type(dual) :: res

    if (val1%val < val2%val) then
      res = val1
    else
      res = val2
    endif
    if(present(val3))then
      if (res%val > val3%val) res = val3
    endif
    if(present(val4))then
      if (res%val > val4%val) res = val4
    endif

  end function

  elemental function min_dr(u, r) result(res)
    type(dual), intent(in) :: u
    real(wp), intent(in) :: r
    type(dual) :: res

    if (u%val < r) then
      res = u
    else
      allocate(res%der(size(u%der)))
      res = r
    endif

  end function

  elemental function min_rd(n, u) result(res)
    real(wp), intent(in) :: n
    type(dual), intent(in) :: u
    type(dual) :: res

    if (u%val < n) then
      res = u
    else
      allocate(res%der(size(u%der)))
      res = n
    endif

  end function

  function minval_d(u) result(res)
    type(dual), intent(in) :: u(:)
    type(dual) :: res
    integer :: i

    i = minloc(u%val, 1)
    res = u(i)

  end function

  elemental function sin_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    res%val = sin(u%val)
    res%der = cos(u%val) * u%der

  end function

  elemental function tan_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res

    res%val = tan(u%val)
    res%der = u%der / cos(u%val)**2

  end function

  elemental function sqrt_d(u) result(res)
    type(dual), intent(in) :: u
    type(dual) :: res
    integer :: i

    res%val = sqrt(u%val)

    if (res%val /= 0.0_wp) then
      res%der = 0.5_wp * u%der / res%val
    else
      allocate(res%der(size(u%der)))
      do i = 1, size(u%der)
        if (u%der(i) == 0.0_wp) then
          res%der(i) = 0.0_wp
        else
          res%der(i) = ieee_value(1.0_wp, ieee_quiet_nan)
        end if
      end do
    end if

  end function

  function sum_d(u) result(res)
    type(dual), intent(in) :: u(:)
    type(dual) :: res
    integer :: i, j

    res%val = sum(u%val)
    if (size(u) == 0) then
      error stop "Can not sum a zero-length dual array."
    endif
    allocate(res%der(size(u(1)%der)))
    res%der = 0.0_wp
    do j = 1,size(u)
      do i = 1,size(u(1)%der)
        res%der(i) = res%der(i) + u(j)%der(i)
      end do
    enddo

  end function

  function maxloc_d(array) result(ind)
    type(dual), intent(in) :: array(:)
    integer :: ind(1)

    ind = maxloc(array%val)

  end function

end module