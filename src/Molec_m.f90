module Molec_m
  USE NumParameters_m
  implicit none
  private

  real(kind=Rk) :: mass = ONE

  public :: Calc_pot,mass

contains
  FUNCTION Calc_pot(Q)
    real(kind=Rk)             :: Calc_pot

    real(kind=Rk), intent(in) :: Q

  Calc_pot = HALF * Q*Q

  END FUNCTION Calc_pot

end module Molec_m
