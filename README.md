# GPU-Accelerated Energy-Conserving Methods for the Hyperbolized Serre-Green-Naghdi Equations in 2D

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.svg)](https://doi.org/10.5281/zenodo.18038219)

This repository contains information and code to reproduce the results presented
in the article
```bibtex
@online{wittenstein2026gpu,
  title={{GPU}-Accelerated Energy-Conserving Methods for the Hyperbolized {S}erre-{G}reen-{N}aghdi Equations in {2D}},
  author={Wittenstein, Collin and Marks, Vincent and Ricchiuto, Mario and Ranocha, Hendrik},
  year={2026},
  eprint={2601.02540},
  archivePrefix={arXiv},
  primaryClass={math.NA},
  url={https://arxiv.org/abs/2601.02540}, 
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{wittenstein2025hyperbolizedRepro,
  title={Reproducibility repository for
         "{GPU}-Accelerated Energy-Conserving Methods for the Hyperbolized {S}erre-{G}reen-{N}aghdi Equations in {2D}"},
  author={Wittenstein, Collin and Marks, Vincent and Ricchiuto, Mario and Ranocha, Hendrik},
  year={2026},
  howpublished={\url{https://github.com/cwittens/2025_Hyp_SGN_2D}},
  doi={10.5281/zenodo.18038219}
}
```

## Abstract

We develop energy-conserving numerical methods for a two-dimensional hyperbolic
approximation of the Serre–Green–Naghdi equations with variable bathymetry
for both periodic and reflecting boundary conditions. The hyperbolic
formulation avoids the costly inversion of an elliptic operator present in the
classical model. Our schemes combine split forms with summation-by-parts (SBP)
operators to construct semidiscretizations that conserve the total water mass
and the total energy. We provide analytical proofs of these conservation
properties and also verify them numerically. While the framework is general,
our implementation focuses on second-order finite-difference SBP operators. The
methods are implemented in Julia for CPU and GPU architectures (AMD and NVIDIA)
and achieve substantial speedups on modern accelerators. We validate the
approach through convergence studies based on solitary-wave and
manufactured-solution tests, and by comparisons to analytical, experimental,
and existing numerical results. All source code to reproduce our results is
available online.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/).
The numerical experiments presented in this article were performed using
Julia v1.12.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors

- [Collin Wittenstein](https://github.com/cwittens) (Johannes Gutenberg University Mainz, Germany)
- [Vincent Marks](https://github.com/vincmarks) (Johannes Gutenberg University Mainz, Germany)
- [Mario Ricchiuto](https://team.inria.fr/cardamom/marioricchiuto) (Inria at University of Bordeaux, France)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
