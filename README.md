# Fast Multipole Method for Electromagnetics
![uniform_quad_n10000](https://github.com/user-attachments/assets/8824a917-6583-4daf-a728-9f4f3ce0b224)

C++ implementation of 2D / 3D Fast Multipole Method for computing N-particle electrostatic interactions.

## Compilation

## Execution
Specify the following configuration parameters in config.txt (a sample file is included under /config):

* Mode      : Simulation mode ( 0 - READ, 1 - WRITE )
* Dist      : Particle distribution ( 0 - Uniform, 1 - Gaussian, 2 - Grid )
* Qdist     : Charge distribution ( 0 - Plus, 1 - Minus, 2 - Dipole, 3 - Quadrupole, 4 - Octopole, 5 - Random )
* Exp prec  : Exponential expansion precision ( 0 - order 8, 1 - order 17, 2 - order 26)
* Order     : Multipole expansion order
* Nsrcs     : Number of particles to simulate
* Root leng : Length of root node
* Max parts : Maximum number of particles contained in a leaf node
* Do direct : Do direct (particle to particle) computation after FMM 