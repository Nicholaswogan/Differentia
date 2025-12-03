program bench_dense_jacobian_inplace
  use differentia, only: wp, jacobian, JacobianWorkMemory
  use differentia_inplace, only: matvec_dual
  implicit none

  integer, parameter :: n = 100
  integer, parameter :: iterations = 1000

  real(wp), allocatable :: A(:,:)
  real(wp) :: x(n), f(n), J(n,n)
  real(wp) :: total, average
  integer :: r
  real(wp) :: t_start, t_end
  character(:), allocatable :: err
  type(JacobianWorkMemory) :: wrk

  allocate(A(n,n))
  call build_inputs(A, x)

  wrk = JacobianWorkMemory(n, err=err)
  if (allocated(err)) error stop err

  ! warm-up
  call jacobian(dense_fun, x, f, J, wrk=wrk, err=err)
  if (allocated(err)) error stop err

  call cpu_time(t_start)
  do r = 1, iterations
    call jacobian(dense_fun, x, f, J, wrk=wrk, err=err)
    if (allocated(err)) error stop err
  end do
  call cpu_time(t_end)

  total = t_end - t_start
  average = total / real(iterations, wp)

  write(*,'("n=",I0,", iterations=",I0,", avg_time_s=",ES12.5)') n, iterations, average

contains

  subroutine build_inputs(A, x)
    real(wp), intent(out) :: A(:,:), x(:)
    integer :: i, j, nloc
    nloc = size(x)
    do j = 1, nloc
      x(j) = 0.01_wp*real(j, wp)
    end do
    do j = 1, nloc
      do i = 1, nloc
        A(i,j) = 0.001_wp*sin(0.01_wp*real(i + j, wp)) + &
                 0.001_wp*cos(0.02_wp*real(i*j, wp))
      end do
    end do
  end subroutine build_inputs

  subroutine dense_fun(x, f)
    use differentia, only: dual
    type(dual), target, intent(in) :: x(:)
    type(dual), target, intent(inout) :: f(:)
    integer :: i, ndv

    ndv = size(x)

    call matvec_dual(A, x, f)

    do i = 1, ndv
      f(i)%val = f(i)%val + sin(x(i)%val)
      f(i)%der = f(i)%der + cos(x(i)%val)*x(i)%der
    end do
  end subroutine dense_fun

end program bench_dense_jacobian_inplace
