# Enhanced-Sampling Molecular Dynamics Software based on Pure Functions in Fortran
Author: Shinji Iida

This MD software aims to simulate simple physical systems, such as lattice models.
1-d, 2-d and 3-d system can be handled, enabling testing very simple system's behaviour.

Most functions in this software are made to pure and thus adding new functions would be easy.

## Implimentation
* System:
    - Face-Centered Cubic lattice

* Potential:
    - Lenard Jones

* Periodic boundary

* Integrator:
    - Velocity Verlet

* Thermostat:
    - Velocity rescaling with heating setting

Following functionalities will be implemented:
  * Usability of PDB file as an input of a simulated system
  * 1D and 2D simple model (Harmonic system) simulator with and without periodic boundary condition
  * Calculation of long-range interaction via Ewald, Particle Mesh Ewald, and Zero-multipole summation.
  * Amber Force Field: bond, angle, torsion, improer tosion.
  * Book Keeping
  * Leap-frog integrator
  * Nose-Hoover thermostat
  * Multicanonical Molecular Dynamics
  * Metadynamics
  * GO-like model
  * CUDA Fortran
