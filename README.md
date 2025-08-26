## Fast Multipole Method for Electromagnetics
![uniform_quad_n10000](https://github.com/user-attachments/assets/8824a917-6583-4daf-a728-9f4f3ce0b224)

C++ implementation of 2D / 3D Fast Multipole Method (FMM) for computing `N`-body field-mediated interactions.

### Algorithm

The FMM algorithm accelerates computation of `N`-particle interactions represented by Laplace/Helmholtz-like kernels,
effecting a speedup from order `O(N^2)` to `O(N log N)`. It recursively divides the computational domain into successively 
finer boxes, aggregates source contributions via their multipole expansions, and transforms these into local expansions 
in each box from which solutions due to distant sources can be quickly evaluated.

This implementation of the 3D FMM further uses plane wave expansions to accelerate the multipole to local (M2L)
transforms, traditionally the performance bottleneck for the whole algorithm. Future work is planned to extend the 3D FMM
into frequency and time-domain simulations.

### Build

##### Requirements

* C++ compiler with C++ 17 or newer
* Eigen (library included)

### Run
Specify the following configuration parameters in `/config/config.txt` (a sample file is included):

* `Mode`      : Input mode ( `0` - READ, `1` - WRITE )
* `Pdist`     : Particle distribution ( `0` - Uniform, `1` - Gaussian, `2` - Sphere, `3` - Cylinder )
* `Qdist`     : Charge distribution ( `0` - Plus, `1` - Minus, `2` - Dipole, `3` - Quadrupole, `4` - Octupole, `5` - Random )
* `Exp prec`  : Exponential expansion precision ( `0` - order 8, `1` - order 17, `2` - order 26)
* `Order`     : Multipole expansion order
* `Nsrcs`     : Number of particles to simulate
* `Root leng` : Length of root node
* `Max parts` : Maximum number of particles in a leaf node
* `Do direct` : Do direct (particle to particle) computation after FMM 

Run the executable, which outputs potential and field data to `/out`. Then plot with your choice of software.