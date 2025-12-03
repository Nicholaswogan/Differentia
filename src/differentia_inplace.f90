module differentia_inplace
  use differentia_const, only: wp
  use differentia_dual,  only: dual
  implicit none
  private

  public :: zero_dual
  public :: add_dual, add_dual_real
  public :: mult_dual, mult_dual_real
  public :: sum_dual
  public :: matvec_dual

contains

  subroutine zero_dual(y)
    type(dual), intent(inout) :: y(:)
    integer :: i

    if (size(y) == 0) return
    if (.not. allocated(y(1)%der)) error stop "dual_zero: y%der not allocated"
    do i = 1, size(y)
      if (.not. allocated(y(i)%der)) error stop "zero_dual: y%der not allocated"
      y(i)%val = 0.0_wp
      y(i)%der = 0.0_wp
    end do
  end subroutine

  subroutine add_dual(a, b, res)
    type(dual), intent(in)    :: a, b
    type(dual), intent(inout) :: res

    if (.not. allocated(a%der) .or. .not. allocated(b%der)) error stop "add_dual: input der not allocated"
    if (allocated(res%der)) then
      if (size(res%der) /= size(a%der)) error stop "add_dual: res der size mismatch"
    else
      allocate(res%der(size(a%der)))
    end if

    res%val = a%val + b%val
    res%der = a%der + b%der
  end subroutine

  subroutine add_dual_real(a, r, res)
    type(dual), intent(in)    :: a
    real(wp), intent(in)      :: r
    type(dual), intent(inout) :: res

    if (.not. allocated(a%der)) error stop "add_dual_real: input der not allocated"
    if (allocated(res%der)) then
      if (size(res%der) /= size(a%der)) error stop "add_dual_real: res der size mismatch"
    else
      allocate(res%der(size(a%der)))
    end if

    res%val = a%val + r
    res%der = a%der
  end subroutine

  subroutine mult_dual(a, b, res)
    type(dual), intent(in)    :: a, b
    type(dual), intent(inout) :: res

    if (.not. allocated(a%der) .or. .not. allocated(b%der)) error stop "mul_dual: input der not allocated"
    if (allocated(res%der)) then
      if (size(res%der) /= size(a%der)) error stop "mul_dual: res der size mismatch"
    else
      allocate(res%der(size(a%der)))
    end if

    res%val = a%val * b%val
    res%der = a%der * b%val + b%der * a%val
  end subroutine

  subroutine mult_dual_real(a, r, res)
    type(dual), intent(in)    :: a
    real(wp), intent(in)      :: r
    type(dual), intent(inout) :: res

    if (.not. allocated(a%der)) error stop "mul_dual_real: input der not allocated"
    if (allocated(res%der)) then
      if (size(res%der) /= size(a%der)) error stop "mul_dual_real: res der size mismatch"
    else
      allocate(res%der(size(a%der)))
    end if

    res%val = a%val * r
    res%der = a%der * r
  end subroutine

  subroutine sum_dual(vec, res)
    type(dual), intent(in)    :: vec(:)
    type(dual), intent(inout) :: res
    integer :: i

    if (size(vec) == 0) error stop "sum_dual: zero-length input"
    if (.not. allocated(vec(1)%der)) error stop "sum_dual: input der not allocated"
    if (allocated(res%der)) then
      if (size(res%der) /= size(vec(1)%der)) error stop "sum_dual: res der size mismatch"
    else
      allocate(res%der(size(vec(1)%der)))
    end if

    res%val = 0.0_wp
    res%der = 0.0_wp
    do i = 1, size(vec)
      res%val = res%val + vec(i)%val
      res%der = res%der + vec(i)%der
    end do
  end subroutine

  subroutine matvec_dual(A, x, y)
    real(wp), intent(in)      :: A(:,:)
    type(dual), intent(in)    :: x(:)
    type(dual), intent(inout) :: y(:)
    integer :: i, j, ndv

    if (size(A,1) /= size(y) .or. size(A,2) /= size(x)) error stop "matvec_dual: dimension mismatch"
    if (size(y) == 0) return
    if (.not. allocated(x(1)%der)) error stop "matvec_dual: x%der not allocated"
    ndv = size(x(1)%der)

    do i = 1, size(y)
      if (.not. allocated(y(i)%der)) error stop "matvec_dual: y%der not allocated"
      if (size(y(i)%der) /= ndv) error stop "matvec_dual: y der size mismatch"
      y(i)%val = 0.0_wp
      y(i)%der = 0.0_wp
    end do

    do j = 1, size(x)
      do i = 1, size(y)
        y(i)%val = y(i)%val + A(i,j) * x(j)%val
        y(i)%der = y(i)%der + A(i,j) * x(j)%der
      end do
    end do
  end subroutine

end module differentia_inplace
