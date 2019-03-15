module potentials
  implicit none

contains

  pure function springPot(k, r)
    real(kind = 8), intent(in) :: k, r
    real(kind = 8) :: springPot
    springPot = 0.5d0 * k * r**2
  end function springPot

  pure function springForce(k, r)
    real(kind = 8), intent(in) :: k, r
    real(kind = 8) :: springForce
    springForce = k
  end function springForce

  pure function lennardJones(r) result(potential)
    real(kind = 8), intent(in) :: r
    real(kind = 8) :: potential
    
    potential = 4 * (r**(-12) - r**(-6))
  end function

  pure function ljForce(r)
    real(kind = 8), intent(in) :: r
    real(kind = 8) :: ljForce
    ljForce = 48.0d0 * (r**(-14) - 0.5d0*r**(-8))
  end function

  pure function mcmdEnergy(fitparams)
    real(kind = 8), intent(in) :: fitparams
    real(kind = 8) :: mcmdEnergy
  end function


end module potentials
