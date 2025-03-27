module forwarddiff_derivative
  use forwarddiff_const, only: wp
  use forwarddiff_dual, only: dual
  implicit none
  private

  !> Work memory for the `jacobian` routine.
  type :: JacobianWorkMemory
    integer :: jt
    type(dual), allocatable :: xx(:), ff(:)
  end type
  interface JacobianWorkMemory
    module procedure :: create_JacobianWorkMemory
  end interface

  public :: JacobianWorkMemory
  public :: derivative, derivative_sig
  public :: gradient, gradient_sig
  public :: jacobian, jacobian_sig

  abstract interface
    !> Interface for the function input to the `derivative`
    !> routine for computing the derivative of scalar 
    !> functions.
    function derivative_sig(x) result(res)
      import :: dual
      type(dual), intent(in) :: x !! Input scalar
      type(dual) :: res !! Resulting scalar
    end function

    !> Interface for the function input to the `gradient`
    !> subroutine for computing the gradient of a function
    !> mapping a vector to a scalar.
    function gradient_sig(x) result(res)
      import :: dual
      type(dual), target, intent(in) :: x(:) !! Input vector
      type(dual) :: res !! Resulting scalar
    end function

    !> Interface for the function input to the `jacobian`
    !> subroutine for computing the jacobian of a function
    !> mapping a vector to a vector.
    subroutine jacobian_sig(x, f)
      import :: dual
      type(dual), target, intent(in) :: x(:) !! Input vector
      type(dual), target, intent(out) :: f(:) !! Resulting vector
    end subroutine
  end interface

contains

  !> Computes the derivative of input scalar function
  !> `fcn` at the input `x`.
  subroutine derivative(fcn, x, f, dfdx)
    procedure(derivative_sig) :: fcn !! Input scalar function
    real(wp), intent(in) :: x
    real(wp), intent(out) :: f !! `fcn` evaluated at `x`
    real(wp), intent(out) :: dfdx !! The derivative of `fcn` at `x`
    type(dual) :: ff
    ff = fcn(dual(x, [1.0_wp]))
    f = ff%val
    dfdx = ff%der(1)
  end subroutine

  !> Computes the gradient of the input function
  !> `fcn` at the input `x`.
  subroutine gradient(fcn, x, f, dfdx, err)
    procedure(gradient_sig) :: fcn !! Input function mapping a vector to a scalar
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f !! `fcn` evaluated at `x`
    real(wp), intent(out) :: dfdx(:) !! The gradient of `fcn` at `x`
    !> If an error occurs, `err` will be allocated with an error message.
    character(:), allocatable, intent(out) :: err

    type(dual) :: xx(size(x))
    type(dual) :: ff
    integer :: i

    if (size(x) /= size(dfdx)) then
      err = 'Output `dfdx` array is not the right size.'
      return
    endif

    ! Set x
    xx%val = x

    ! Seed the dual number
    do i = 1,size(x)
      allocate(xx(i)%der(size(x)))
      xx(i)%der = 0.0_wp
      xx(i)%der(i) = 1.0_wp
    enddo

    ! Do differentiation
    ff = fcn(xx)

    ! Unpack f(x)
    f = ff%val

    ! Unpack gradient
    dfdx(:) = ff%der(:)

  end subroutine

  function create_JacobianWorkMemory(n, jt, bandwidth, blocksize, err) result(wrk)
    use forwarddiff_const, only: DenseJacobian, BandedJacobian, BlockDiagonalJacobian
    integer, intent(in) :: n
    integer, optional, intent(in) :: jt
    integer, optional, intent(in) :: bandwidth
    integer, optional, intent(in) :: blocksize
    character(:), allocatable, intent(out) :: err
    type(JacobianWorkMemory) :: wrk

    integer :: i, jt_

    ! Determine type of jacobian
    if (present(jt)) then
      jt_ = jt
    else
      jt_ = DenseJacobian
    endif

    wrk%jt = jt_
    allocate(wrk%xx(n))
    allocate(wrk%ff(n))
    if (jt_ == DenseJacobian) then
      do i = 1,n
        allocate(wrk%xx(i)%der(n))
        allocate(wrk%ff(i)%der(n))
      enddo
    elseif (jt_ == BandedJacobian) then
      if (.not.present(bandwidth)) then
        err = '`bandwidth` must be an argument when computing a banded jacobian.'
        return
      endif
      do i = 1,n
        allocate(wrk%xx(i)%der(bandwidth))
        allocate(wrk%ff(i)%der(bandwidth))
      enddo
    elseif (jt_ == BlockDiagonalJacobian) then
      if (.not.present(blocksize)) then
        err = '`blocksize` must be an argument when computing a block diagonal jacobian.'
        return
      endif
      do i = 1,n
        allocate(wrk%xx(i)%der(blocksize))
        allocate(wrk%ff(i)%der(blocksize))
      enddo
    else
      err = 'Invalid value for the Jacobian type indicator `jt`.'
      return
    endif

  end function

  !> Computes the Jacobian of the input function `fcn
  !> at the input `x`, with support for some sparse Jacobians. Only can
  !> consider square Jacobians.
  subroutine jacobian(fcn, x, f, dfdx, wrk, jt, bandwidth, blocksize, err)
    use forwarddiff_const, only: DenseJacobian, BandedJacobian, BlockDiagonalJacobian
    procedure(jacobian_sig) :: fcn !! Input function mapping a vector to a vector
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:) !! `fcn` evaluated at `x`
    real(wp), intent(out) :: dfdx(:,:) 
    !! The Jacobian of `fcn` evaluated at `x`. For dense Jacobians (`jt == DenseJacobian`),
    !! `dfdx` has shape (n,n) where n is `size(x)`. The first dimension indexes rows 
    !! of the Jacobian, while the second dimension indexes columns. For example,
    !!
    !!     | df(1,1) df(1,2) df(1,3) |
    !! J = | df(2,1) df(2,2) df(2,3) |
    !!     | df(3,1) df(3,2) df(3,3) |
    !!
    !! If the Jacobian is banded (`jt == BandedJacobian`), then `dfdx` has shape 
    !! (bandwidth,n), where bandwidth is the Jacobian bandwidth. In this case, the
    !! The diagonals of the Jacobian are loaded into the rows of `dfdx`
    !!
    !!     | df(1,1) df(1,2) 0       0       0       |    
    !!     | df(2,1) df(2,2) df(2,3) 0       0       |  ->  | 0       df(1,2) df(2,3) df(3,4) df(4,5) |
    !! J = | 0       df(3,2) df(3,3) df(3,4) 0       |  ->  | df(1,1) df(2,2) df(3,3) df(4,4) df(5,5) |
    !!     | 0       0       df(4,3) df(4,4) df(4,5) |  ->  | df(2,1) df(3,2) df(4,3) df(5,4) 0       |
    !!     | 0       0       0       df(5,4) df(5,5) |
    !!
    !! If the Jacobian is block-diagonal (`jt == BlockDiagonalJacobian`), then `dfdx` has
    !! shape (blocksize,n) where blocksize is the blocksize of the Jacobian. In this case,
    !! `dfdx` stores the Jacobian entries in the following way:
    !!
    !!     | df(1,1) df(1,2) 0       0       |      
    !! J = | df(2,1) df(2,2) 0       0       |  ->  | df(1,1) df(1,2) df(3,3) df(3,4) |
    !!     | 0       0       df(3,3) df(3,4) |  ->  | df(2,1) df(2,2) df(4,3) df(4,4) |
    !!     | 0       0       df(4,3) df(4,4) |
    type(JacobianWorkMemory), target, optional, intent(inout) :: wrk
    !! Work memory for the calculation. If not provided, memory will by dynamically allocated
    !! and destroyed during the calculation (slightly slower).
    integer, optional, intent(in) :: jt
    !! Jacobian sparsity indicator. The default is `jt = DenseJacobian`.
    !! - If `jt == DenseJacobian`, then the algorithm assumes the Jacobian is dense.
    !! - If `jt == BandedJacobian`, then the algorithm assumes the Jacobian is banded
    !!   with `bandwidth`.
    !! - If `jt == BlockDiagonalJacobian`, then the algorithm assumes the Jacobian is
    !!   block diagonal with `blocksize`. Note that size(x) must be an integer multiple
    !!   of `blocksize`.
    integer, optional, intent(in) :: bandwidth
    !! The Jacobian bandwidth for `jt == BandedJacobian`
    integer, optional, intent(in) :: blocksize
    !! The Jacobian blocksize for `jt == BlockDiagonalJacobian`
    character(:), allocatable, intent(out) :: err
    !! If an error occurs, `err` will be allocated with an error message.

    integer :: jt_
    type(JacobianWorkMemory), target :: wrk_tmp
    type(JacobianWorkMemory), pointer :: wrk_ptr

    ! Check dimensions work out
    if (size(x) /= size(f)) then
      err = 'Output `f` array is not the right size.'
      return
    endif

    ! Determine type of jacobian
    if (present(jt)) then
      jt_ = jt
    else
      jt_ = DenseJacobian
    endif

    if (present(wrk)) then
      wrk_ptr => wrk
      if (wrk_ptr%jt /= jt_) then
        err = 'The work memory has a Jacobian type with the `jt` input for the subroutine `jacobian`'
        return
      endif
    else
      wrk_tmp = JacobianWorkMemory(size(x), jt, bandwidth, blocksize, err)
      if (allocated(err)) return
      wrk_ptr => wrk_tmp
    endif

    if (jt_ == DenseJacobian) then
      call jacobian_dense(fcn, x, f, dfdx, wrk_ptr, err)
      if (allocated(err)) return
    elseif (jt_ == BandedJacobian) then
      if (.not.present(bandwidth)) then
        err = '`bandwidth` must be an argument when computing a banded jacobian.'
        return
      endif
      call jacobian_banded(fcn, x, f, dfdx, wrk_ptr, bandwidth, err)
      if (allocated(err)) return
    elseif (jt_ == BlockDiagonalJacobian) then
      if (.not.present(blocksize)) then
        err = '`blocksize` must be an argument when computing a block diagonal jacobian.'
        return
      endif
      call jacobian_blockdiagonal(fcn, x, f, dfdx, wrk_ptr, blocksize, err)
      if (allocated(err)) return
    else
      err = 'Invalid value for the Jacobian type indicator `jt`.'
      return
    endif

  end subroutine

  subroutine jacobian_dense(fcn, x, f, dfdx, wrk, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    type(JacobianWorkMemory), target, intent(inout) :: wrk
    character(:), allocatable, intent(out) :: err

    type(dual), pointer :: xx(:)
    type(dual), pointer :: ff(:)
    integer :: i, j

    if (size(x) /= size(dfdx,1) .or. size(x) /= size(dfdx,2)) then
      err = 'Output `dfdx` array is not the right size.'
      return
    endif

    xx => wrk%xx
    ff => wrk%ff
    
    ! Set x
    xx%val = x

    ! Seed dual
    do i = 1,size(x)
      xx(i)%der = 0.0_wp
      xx(i)%der(i) = 1.0_wp
    enddo

    ! Do differentiation
    call fcn(xx, ff)

    ! Unpack f(x)
    f = ff%val

    ! Unpack jacobian.
    ! | df(1,1) df(1,2) df(1,3) |
    ! | df(2,1) df(2,2) df(2,3) |
    ! | df(3,1) df(3,2) df(3,3) |

    do j = 1,size(x)
      do i = 1,size(x)
        dfdx(i,j) = ff(i)%der(j)
      enddo
    enddo
  
  end subroutine

  subroutine jacobian_banded(fcn, x, f, dfdx, wrk, bandwidth, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    type(JacobianWorkMemory), target, intent(inout) :: wrk
    integer, intent(in) :: bandwidth
    character(:), allocatable, intent(out) :: err

    type(dual), pointer :: xx(:)
    type(dual), pointer :: ff(:)
    integer :: hbw
    integer :: i, j, ii, jj, kk

    ! Check bandwidth
    if (bandwidth > size(x)) then
      err = '`bandwidth` can not be > size(x).'
      return
    endif
    if (bandwidth < 1) then
      err = '`bandwidth` can not be < 1.'
      return
    endif
    if (mod(bandwidth,2) == 0) then
      err = '`bandwidth` must be odd.'
      return
    endif

    if (bandwidth /= size(dfdx,1) .or. size(x) /= size(dfdx,2)) then
      err = 'Output `dfdx` array is not the right size.'
      return
    endif

    xx => wrk%xx
    ff => wrk%ff

    ! Set x
    xx%val = x

    ! Allocate dual number
    do i = 1,size(x)
      xx(i)%der = 0.0_wp
    enddo

    ! Seed the dual number
    j = 1
    outer : do
      do i = 1,bandwidth
        if (j > size(x)) exit outer 
        xx(j)%der(i) = 1.0_wp
        j = j + 1
      enddo
    enddo outer

    ! Do differentiation
    call fcn(xx, ff)

    ! Unpack f(x)
    f = ff%val

    ! Unpack banded jacobian
    ! In this case, we load diagonals of jacobian into each row of dfdx.
    ! So, dfdx(1,:) is the "highest" diagonal. Illustration:
    !
    ! | df(1,1) df(1,2) 0       0       0       |      
    ! | df(2,1) df(2,2) df(2,3) 0       0       |  ->  | 0       df(1,2) df(2,3) df(3,4) df(4,5) |
    ! | 0       df(3,2) df(3,3) df(3,4) 0       |  ->  | df(1,1) df(2,2) df(3,3) df(4,4) df(5,5) |
    ! | 0       0       df(4,3) df(4,4) df(4,5) |  ->  | df(2,1) df(3,2) df(4,3) df(5,4) 0       |
    ! | 0       0       0       df(5,4) df(5,5) |
    !

    hbw = (bandwidth - 1)/2 ! halfbandwidth 
    j = 1
    outer1 : do
      do i = 1,bandwidth
        if (j > size(x)) exit outer1
        do jj = -hbw,hbw
          kk = jj + hbw + 1
          ii = j + jj
          if (ii < 1) then
            dfdx(kk,j) = 0.0_wp
          elseif (ii > size(x)) then
            dfdx(kk,j) = 0.0_wp
          else
            dfdx(kk,j) = ff(ii)%der(i)
          endif
        enddo
        j = j + 1
      enddo
    enddo outer1

  end subroutine

  subroutine jacobian_blockdiagonal(fcn, x, f, dfdx, wrk, blocksize, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    type(JacobianWorkMemory), target, intent(inout) :: wrk
    integer, intent(in) :: blocksize
    character(:), allocatable, intent(out) :: err

    type(dual), pointer :: xx(:)
    type(dual), pointer :: ff(:)
    integer :: i, j, ii, jj

    ! Check blocksize
    if (blocksize > size(x)) then
      err = '`blocksize` can not be > size(x).'
      return
    endif
    if (blocksize < 1) then
      err = '`blocksize` can not be < 1.'
      return
    endif

    if (mod(size(x),blocksize) /= 0) then
      err = 'size(x) must be an integer multiple of `blocksize`.'
      return
    endif

    if (blocksize /= size(dfdx,1) .or. size(x) /= size(dfdx,2)) then
      err = 'Output `dfdx` array is not the right size.'
      return
    endif

    xx => wrk%xx
    ff => wrk%ff

    ! Set x
    xx%val = x

    ! Allocate dual number
    do i = 1,size(x)
      xx(i)%der = 0.0_wp
    enddo

    ! Seed the dual number
    j = 1
    outer : do
      do i = 1,blocksize
        if (j > size(x)) exit outer 
        xx(j)%der(i) = 1.0_wp
        j = j + 1
      enddo
    enddo outer

    ! Do differentiation
    call fcn(xx, ff)

    ! Unpack f(x)
    f = ff%val

    ! Unpack block jacobian
    !
    ! | df(1,1) df(1,2) 0       0       |      
    ! | df(2,1) df(2,2) 0       0       |  ->  | df(1,1) df(1,2) df(3,3) df(3,4) |
    ! | 0       0       df(3,3) df(3,4) |  ->  | df(2,1) df(2,2) df(4,3) df(4,4) |
    ! | 0       0       df(4,3) df(4,4) |
    do ii = 1,size(x)/blocksize
      jj = (ii - 1)*blocksize
      do i = 1,blocksize
        do j = 1,blocksize
          dfdx(i,j+jj) = ff(i+jj)%der(j)
        enddo
      enddo
    enddo

  end subroutine

end module