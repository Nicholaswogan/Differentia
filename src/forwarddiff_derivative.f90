module forwarddiff_derivative
  use forwarddiff_const, only: wp
  use forwarddiff_dual, only: dual
  implicit none
  private

  public :: derivative, derivative_sig
  public :: gradient, gradient_sig
  public :: jacobian, jacobian_sig

  abstract interface
    function derivative_sig(x) result(res)
      import :: dual
      type(dual), intent(in) :: x
      type(dual) :: res
    end function

    function gradient_sig(x) result(res)
      import :: dual
      type(dual), target, intent(in) :: x(:)
      type(dual) :: res
    end function

    subroutine jacobian_sig(x, f)
      import :: dual
      type(dual), target, intent(in) :: x(:)
      type(dual), target, intent(out) :: f(:)
    end subroutine
  end interface

contains

  subroutine derivative(fcn, x, f, dfdx)
    procedure(derivative_sig) :: fcn
    real(wp), intent(in) :: x
    real(wp), intent(out) :: f
    real(wp), intent(out) :: dfdx
    type(dual) :: ff
    ff = fcn(dual(x, [1.0_wp]))
    f = ff%val
    dfdx = ff%der(1)
  end subroutine

  subroutine gradient(fcn, x, f, dfdx, err)
    procedure(gradient_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f
    real(wp), intent(out) :: dfdx(:)
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

  subroutine jacobian(fcn, x, f, dfdx, jt, bandwidth, blocksize, err)
    use forwarddiff_const, only: DenseJacobian, BandedJacobian, BlockDiagonalJacobian
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    integer, optional, intent(in) :: jt
    integer, optional, intent(in) :: bandwidth
    integer, optional, intent(in) :: blocksize
    character(:), allocatable, intent(out) :: err

    integer :: jt_

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

    if (jt_ == DenseJacobian) then
      call jacobian_dense(fcn, x, f, dfdx, err)
      if (allocated(err)) return
    elseif (jt_ == BandedJacobian) then
      if (.not.present(bandwidth)) then
        err = '`bandwidth` must be an argument when computing a banded jacobian.'
        return
      endif
      call jacobian_banded(fcn, x, f, dfdx, bandwidth, err)
      if (allocated(err)) return
    elseif (jt_ == BlockDiagonalJacobian) then
      if (.not.present(blocksize)) then
        err = '`blocksize` must be an argument when computing a block diagonal jacobian.'
        return
      endif
      call jacobian_blockdiagonal(fcn, x, f, dfdx, blocksize, err)
      if (allocated(err)) return
    endif

  end subroutine

  subroutine jacobian_dense(fcn, x, f, dfdx, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    character(:), allocatable, intent(out) :: err

    type(dual) :: xx(size(x))
    type(dual) :: ff(size(x))
    integer :: i, j

    if (size(x) /= size(dfdx,1) .or. size(x) /= size(dfdx,2)) then
      err = 'Output `dfdx` array is not the right size.'
      return
    endif

    ! Set x
    xx%val = x

    ! Allocate and seed dual
    do i = 1,size(x)
      allocate(xx(i)%der(size(x)))
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

  subroutine jacobian_banded(fcn, x, f, dfdx, bandwidth, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    integer, intent(in) :: bandwidth
    character(:), allocatable, intent(out) :: err

    type(dual) :: xx(size(x))
    type(dual) :: ff(size(x))
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

    ! Set x
    xx%val = x

    ! Allocate dual number
    do i = 1,size(x)
      allocate(xx(i)%der(bandwidth))
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

  subroutine jacobian_blockdiagonal(fcn, x, f, dfdx, blocksize, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    integer, intent(in) :: blocksize
    character(:), allocatable, intent(out) :: err

    type(dual) :: xx(size(x))
    type(dual) :: ff(size(x))
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

    ! Set x
    xx%val = x

    ! Allocate dual number
    do i = 1,size(x)
      allocate(xx(i)%der(blocksize))
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