# Molecular Dynamics Software based on Pure Functions in Fortran
Author: Shinji Iida

This MD software aims to simulate simple physical systems, such as lattice models.
1-d, 2-d and 3-d system can be handled, so it can be used for testing very simple system's behaviour.

Most functions in this software are made to pure and thus adding new functions would be easy.

## Implimentation
* System:
    - Face-Centered Cubic lattice, 2d and 1d simple lattice

* Potential:
    - Lenard Jones

* Periodic boundary

* Integrator:
    - Velocity Verlet

* Thermostat:
    - velocity rescaling

Following functionality will be implemented:
  * Input of PDB file for a simulated system
  * Calculation of long-range interaction via Ewald, Particle Mesh Ewald, and Zero- multipole summation.
  * Amber Force Field, Harmonic potential,
  * Book Keeping
  * Leap-frog integrator
  * Nose-Hoover thermostat
