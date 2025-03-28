
program test_sparse
  use differentia, only: wp, jacobian
  implicit none

  call test_banded()
  call test_blockdiagonal1()
  call test_blockdiagonal2()

contains

  subroutine test_banded()
    use differentia, only: BandedJacobian
    integer, parameter :: nz = 6
    integer, parameter :: bandwidth = 3
    real(wp) :: u(nz), f(nz), dfdu(bandwidth,nz)
    real(wp) :: f1(nz), dfdu1(nz,nz)
    character(:), allocatable :: err
    integer :: i

    print*,'test_banded'
    
    do i = 1,nz
      u(i) = i
    enddo
    call jacobian(rhs_banded_dual, u, f, dfdu, jt=BandedJacobian, bandwidth=bandwidth, err=err)
    if (allocated(err)) then
      print*,err
      stop
    endif

    do i = 1,bandwidth
      print*,dfdu(i,:)
    enddo
    print*,''

    call jacobian(rhs_banded_dual, u, f1, dfdu1, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    do i = 1,nz
      print*,dfdu1(i,:)
    enddo
    print*,''

  end subroutine

  subroutine rhs_banded_dual(u, du)
    use differentia
    type(dual), target, intent(in) :: u(:)
    type(dual), target, intent(out) :: du(:)
    integer :: i

    du(1) = 3.0_wp*u(2) - u(1)
    do i = 2,size(u)-1
      du(i) = 3.0_wp*u(i+1) - 2.0_wp*u(i) + u(i-1)
    enddo
    du(size(u)) = - u(size(u)) + u(size(u)-1)

  end subroutine

  subroutine test_blockdiagonal1()
    use differentia, only: BlockDiagonalJacobian
    integer, parameter :: nz = 6
    integer, parameter :: blocksize = 2
    real(wp) :: u(nz), f(nz), dfdu(blocksize,nz)
    real(wp) :: f1(nz), dfdu1(nz,nz)
    character(:), allocatable :: err
    integer :: i

    print*,'test_blockdiagonal1'

    do i = 1,nz
      u(i) = i
    enddo
    call jacobian(rhs_blocked1_dual, u, f, dfdu, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    do i = 1,blocksize
      print*,dfdu(i,:) 
    enddo
    print*,''

    call jacobian(rhs_blocked1_dual, u, f1, dfdu1, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    do i = 1,nz
      print*,dfdu1(i,:) 
    enddo
    print*,''

  end subroutine

  subroutine test_blockdiagonal2()
    use differentia, only: BlockDiagonalJacobian
    integer, parameter :: nz = 6
    integer, parameter :: blocksize = 3
    real(wp) :: u(nz), f(nz), dfdu(blocksize,nz)
    real(wp) :: f1(nz), dfdu1(nz,nz)
    character(:), allocatable :: err
    integer :: i

    print*,'test_blockdiagonal2'

    do i = 1,nz
      u(i) = i
    enddo
    call jacobian(rhs_blocked2_dual, u, f, dfdu, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    do i = 1,blocksize
      print*,dfdu(i,:) 
    enddo
    print*,''

    call jacobian(rhs_blocked2_dual, u, f1, dfdu1, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    do i = 1,nz
      print*,dfdu1(i,:) 
    enddo
    print*,''

  end subroutine

  subroutine rhs_blocked1_dual(u, du)
    use differentia
    type(dual), target, intent(in) :: u(:)
    type(dual), target, intent(out) :: du(:)
    integer :: i

    du(1) = u(1) + u(2)*u(1)
    du(2) = u(2) + u(2)*u(1)

    du(3) = u(3)*u(4) + u(3)
    du(4) = u(4) + u(3)*u(3)

    du(5) = 2.0_wp*u(5) + 0.5_wp*u(6)
    du(6) = u(5) + u(6)*u(6)

  end subroutine

  subroutine rhs_blocked2_dual(u, du)
    use differentia
    type(dual), target, intent(in) :: u(:)
    type(dual), target, intent(out) :: du(:)
    integer :: i

    du(1) = u(1) + u(2)*u(1)
    du(2) = u(2) + u(2)*u(1)
    du(3) = u(3) + u(3)*u(1)

    du(4) = u(4) + u(5)*u(6)
    du(5) = 2.0_wp*u(5) + 0.5_wp*u(6)
    du(6) = u(5) + u(6)*u(6)

  end subroutine

end program