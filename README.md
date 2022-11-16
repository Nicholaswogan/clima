# Clima

**NOTE: Clima has not been published in a peer-reviewed journal yet. Please do not use this model in a paper that will be published. Wait for me to publish the model first, so you are able to cite it. Email me if you have questions (wogan@uw.edu)**

`Clima` is a radiative transfer code and a climate model. It is currently under construction.

## Installation

You need a Fortran compiler (`gfortran>=9.30`, [install instructions here](https://fortran-lang.org/learn/os_setup/install_gfortran)) and C compiler (e.g. install with `conda install -c conda-forge clang`)

Create a `conda` environment with all dependencies

```sh
conda create -n clima -c conda-forge python numpy scipy pyyaml scikit-build cython h5py
```

Clone this Gitub repository: 

```sh
git clone --depth=1 --recursive https://github.com/Nicholaswogan/clima.git
```

Navigate to the root directory with a terminal, activate your new `conda` environment, then install with pip:

```sh
conda activate clima
python -m pip install --no-deps --no-build-isolation . -v
# or `python setup.py install`, if the above command doesn't work for some reason
```

If your installation fails, and it is not clear why, please raise an issue.

## Contact

If you have questions email me: wogan@uw.edu

## Funding and Acknowledgements

Funding for the development of Clima comes from
- [The Virtual Planetary Laboratory](https://depts.washington.edu/naivpl/content/welcome-virtual-planetary-laboratory)
- [Simons Collaboration on the Origins of Life](https://www.simonsfoundation.org/life-sciences/origins-of-life/simons-collaboration-on-the-origins-of-life/)

This model was build in collaboration with
- Josh Krissansen-Totton
- Maggie Thompson
- Eric Wolf
- David Catling
- Kevin Zahnle
- Sandra Bastelberger
- Giada Arney
- Shawn Domagal-Goldman