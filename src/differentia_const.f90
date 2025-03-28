module differentia_const
  use iso_fortran_env, only: wp => real64
  implicit none
  public

  enum, bind(c)
    ! Jacobian Type
    enumerator :: &
      DenseJacobian = 1, &
      BandedJacobian = 2, &
      BlockDiagonalJacobian = 3
  end enum

end module