program test_forwarddiff
  use forwarddiff, only: wp, derivative, gradient, jacobian
  implicit none

  call test_dual()
  call test_jacobian()

contains

  subroutine test_dual()
    real(wp) :: x, f, dfdx
    real(wp) :: xx(2), dfdx1(2)
    character(:), allocatable :: err

    open(unit=2,file='test.dat',status='replace',form='unformatted')

    x = 3.0_wp
    call derivative(func_operators, x, f, dfdx)
    print*,f, dfdx
    write(2) f, dfdx

    x = 2.0_wp
    call derivative(func_intrinsics1, x, f, dfdx)
    print*,f, dfdx
    write(2) f, dfdx

    x = 0.1_wp
    call derivative(func_intrinsics2, x, f, dfdx)
    print*,f, dfdx
    write(2) f, dfdx

    xx = [1.0_wp, 2.0_wp]
    call gradient(func_grad1, xx, f, dfdx1, err)
    print*,f, dfdx1
    write(2) f, dfdx1

    xx = [3.0_wp, 4.0_wp]
    call gradient(func_grad2, xx, f, dfdx1, err)
    print*,f, dfdx1
    write(2) f, dfdx1

    close(2)

  end subroutine

  function func_operators(x) result(res)
    use forwarddiff
    type(dual), intent(in) :: x
    type(dual) :: res

    res = x + x + 3.0_wp + x
    res = res - x - 3.0_wp - x
    res = res*x + 2.0_wp*x + x*(-5.0_wp)
    res = res/x + 2.0_wp/x + x/5.0_wp
    res = res**x + res**1.5_wp
    res = 10.0_wp**res

  end function

  function func_intrinsics1(x) result(res)
    use forwarddiff
    type(dual), intent(in) :: x
    type(dual) :: res

    res = abs(x)
    res = res + cos(x)
    res = res + exp(x)
    res = res + log(x)
    res = res + log10(x)
    res = res + sin(x)
    res = res + tan(x)
    res = res + sqrt(x)

  end function

  function func_intrinsics2(x) result(res)
    use forwarddiff
    type(dual), intent(in) :: x
    type(dual) :: res

    res = acos(x)
    res = res + asin(x)
    res = res + atan(x)
    res = max(res, x)
    res = max(res, 1.0_wp)
    res = max(1.0_wp, res)
    res = min(res, res)
    res = min(res, 2.0_wp)
    res = min(2.0_wp, res)

  end function

  function func_grad1(x) result(res)
    use forwarddiff
    type(dual), target, intent(in) :: x(:)
    type(dual) :: res
    res = x(1)*x(1)*x(2) + x(1) + x(2)
  end function

  function func_grad2(x) result(res)
    use forwarddiff
    type(dual), target, intent(in) :: x(:)
    type(dual) :: res
    res = sum(x*3.14_wp)
  end function

  subroutine test_jacobian()
    use forwarddiff, only: JacobianWorkMemory
    real(wp) :: u(3), f(3), dfdu(3,3)
    real(wp) :: f1(3), dfdu1(3,3)
    character(:), allocatable :: err
    type(JacobianWorkMemory) :: wrk

    u = [1.0_wp, 2.0_wp, 3.0_wp]
    call jacobian(rhs_rober_dual, u, f, dfdu, err=err)

    ! Try with work memory preallocated as well
    wrk = JacobianWorkMemory(3, err=err)
    call jacobian(rhs_rober_dual, u, f, dfdu, wrk=wrk, err=err)

    ! print*,f
    print*,''
    print*,dfdu(1,:)
    print*,dfdu(2,:)
    print*,dfdu(3,:)

    ! Check against analytical solution
    call rhs_rober(u, f1)
    call jac_rober(u, dfdu1)

    if (.not. all(f == f1)) then
      error stop "test_jacobian failed"
    endif
    if (.not. all(dfdu == dfdu1)) then
      error stop "test_jacobian failed"
    endif

  end subroutine

  subroutine rhs_rober_dual(u, du)
    use forwarddiff
    type(dual), target, intent(in) :: u(:)
    type(dual), target, intent(out) :: du(:)

    real(wp), parameter :: k1 = 0.04_wp, &
                           k2 = 3.0e7_wp, &
                           k3 = 1.0e4_wp
    
    du(1) = -k1*u(1) + k3*u(2)*u(3)
    du(2) =  k1*u(1) - k2*u(2)**2.0_wp - k3*u(2)*u(3)
    du(3) =  k2*u(2)**2.0_wp

  end subroutine

  subroutine rhs_rober(u, du)
    real(wp), target, intent(in) :: u(:)
    real(wp), target, intent(out) :: du(:)

    real(wp), parameter :: k1 = 0.04_wp, &
                           k2 = 3.0e7_wp, &
                           k3 = 1.0e4_wp
    
    du(1) = -k1*u(1) + k3*u(2)*u(3)
    du(2) =  k1*u(1) - k2*u(2)**2.0_wp - k3*u(2)*u(3)
    du(3) =  k2*u(2)**2.0_wp

  end subroutine
  
  subroutine jac_rober(u, pd)
    real(wp), intent(in) :: u(:)
    real(wp), intent(out) :: pd(:,:)

    real(wp), parameter :: k1 = 0.04_wp, &
                           k2 = 3.0e7_wp, &
                           k3 = 1.0e4_wp
    
    pd(:,1) = [-k1, k1, 0.0_wp] ! column 1
    pd(:,2) = [k3*u(3), -2.0_wp*k2*u(2) - k3*u(3), 2.0_wp*k2*u(2)] ! column 2
    pd(:,3) = [k3*u(2), -k3*u(2), 0.0_wp] ! column 3
    
  end subroutine

end program