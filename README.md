# PyMIEX
PyMIEX is a python wrapped FORTRAN code MIEX &mdash; simulation of Mie scattering in case of arbitrarily large size parameters. MIEX has a number of applications in astrophysics, namely observations of dust and ice in protoplanetary discs. It is also applicable for calculating of optical parameters of the granular surfaces, for example snow or ice with large grains.

Orginal code can be found at: https://data.mendeley.com/datasets/5g7s2crbyf/1; original paper: https://doi.org/10.1016/j.cpc.2004.06.070.

## Getting started

Prerequisites:
- unix (Ubuntu, Debian, Mac etc.) (Windows may also work, but was not tested)
- gfortran (other compilers may also work, were not tested)
- numpy 

Download code:\
`git clone github.com/teimy/pymiex`\
Change directory:\
`cd pymiex`\
Run f2py:\
`python -m numpy.f2py -c -m miex_function miex_function.f90`\
This creates a `.so` file in pymiex directory which you can use to interact with FORTRAN wrapped functions. A few python wrapped functions are located in pymiex, which you can use by importing it: \
`import pymiex as pm` \
Examples of usage are located in `examples.ipynb`.

## How to cite
Please be sure to check original paper by S. Wolf and N.V. Voshchinnikov. Please cite [this](https://doi.org/10.1016/j.cpc.2004.06.070) paper if you use this code in your paper.

## Futher reading
Routine `aa2` is based on [Voshchinnikov, N. (2004). Optics of cosmic dust I.](https://ui.adsabs.harvard.edu/abs/2004ASPRv..12....1V/abstract). \
Another python wrap of Mie scattering by Will Martin: https://github.com/wgm2111/wgm-mie-scattering. It is based on Michael Mishchenko FORTRAN code, highly recommend checking it (and giving a star!).
