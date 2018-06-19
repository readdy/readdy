# Changelog

- simplified c++ potentials interface (#125)
- python method to include potentials defined in a separate extension (#125)
- added headers to conda package (#124) 

## 1.0.1
- added assertions and tests for units (#101)
- reworked cpu kernel (#103)
- compile with gcc7 and c++17 (#103, #107)
- added energy observable (#104)
- added `readdy.__version__` (#107)
- optimized neighbor list (#107)
- added virial and pressure observables (#109, #110)
- added topologies observable (#111, #112)
- fixed topology reactions (#114, #115, #119)
- updated periodic boundary handling (#120)

## 1.0.0
- reactions:
    - fusion, fission, conversion, decay, enzymatic type
    - uncontrolled approximation
    - Gillespie reaction order
- particle propagation:
    - brownian motion
    - external and pair potentials:
        - box potential
        - spherical potential
        - harmonic repulsion
        - lennard-jones
        - piecewise harmonic
        - screened electrostatics
- particle complexes (topologies)
- spatial and temporal changes of these complexes (topology reactions)
- observables:
    - radial distribution function
    - occurred reactions
    - particle positions
    - particles
    - number of particles (per type or in total)
    - forces per particle
    - reaction counts
- kernels:
    - singlecpu
    - cpu (pthread parallelization)
