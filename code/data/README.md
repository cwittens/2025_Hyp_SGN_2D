# Experimental Data

This directory contains experimental and numerical reference datasets used in this paper.

## Available Datasets

### Dingemans Experiment Data (`Dingemans.csv`)

The data are taken from the experiments of Maarten Dingemans:

```bibtex
@techreport{dingemans1994comparison,
  title={Comparison of computations with {B}oussinesq-like models and laboratory measurements},
  author={Maarten W. Dingemans},
  institution={Delft Hydraulics},
  year={1994},
  number={H1684.12},
  url={http://resolver.tudelft.nl/uuid:c2091d53-f455-48af-a84b-ac86680455e9}
}

@book{dingemans1997water,
  title={Water Wave Propagation Over Uneven Bottoms},
  author={Maarten W. Dingemans},
  year={1997},
  volume={13},
  doi={10.1142/1241},
  publisher={World Scientific}
}
```

#### Data Format

The CSV file contains time series data with the first column being time values and subsequent columns representing wave height measurements at different spatial locations (six wave gauges).

#### Usage

In order to access the data, execute:

```julia
t_values, x_values, experimental_data = data_dingemans()
```

The function returns:
- `t_values`: Time values at which measurements were recorded
- `x_values`: Spatial coordinates of the six wave gauges (3.04, 9.44, 20.04, 26.04, 30.44, 37.04)
- `experimental_data`: Matrix of wave heights at each gauge location (columns) over time (rows)

### Henderson Colliding Waves Data (`Henderson/`)

The data are originally from experimental measurements by Professor Diane Henderson for colliding solitary waves, but were digitally scraped from the paper:

```bibtex
@article{DUTYKH2013,
  title = {Finite volume and pseudo-spectral schemes for the fully nonlinear 1D Serre equations},
  volume = {24},
  ISSN = {1469-4425},
  url = {http://dx.doi.org/10.1017/S0956792513000168},
  DOI = {10.1017/s0956792513000168},
  number = {5},
  journal = {European Journal of Applied Mathematics},
  publisher = {Cambridge University Press (CUP)},
  author = {DUTYKH,  DENYS and CLAMOND,  DIDIER and MILEWSKI,  PAUL and MITSOTAKIS,  DIMITRIOS},
  year = {2013},
  month = may,
  pages = {761–787}
}

```

#### Data Format

The directory contains multiple CSV files (`1850.csv`, `1860.csv`, `1870.csv`, `1880.csv`, `1892.csv`, `1900.csv`, `1905.csv`, `1910.csv`, `1915.csv`, `1919.csv`, `1933.csv`, `1950.csv`, `1985.csv`, `2000.csv`) representing wave height measurements at different time instances. Each file contains two columns: x-coordinate and wave height η.

#### Usage

In order to access the data, execute:

```julia
x_data, h_data = data_henderson(i)
```

Where `i` ranges from 1 to 14 corresponding to the different time snapshots. The function returns:
- `x_data`: Spatial coordinates along the experimental domain
- `h_data`: Wave heights at the corresponding spatial locations

### Busto et al. Reference Data (`Busto_et_al.csv`)

Numerical reference data for wave propagation over Gaussian bathymetry from:

```bibtex
@article{Busto2021,
  title = {On High Order ADER Discontinuous Galerkin Schemes for First Order Hyperbolic Reformulations of Nonlinear Dispersive Systems},
  volume = {87},
  ISSN = {1573-7691},
  url = {http://dx.doi.org/10.1007/s10915-021-01429-8},
  DOI = {10.1007/s10915-021-01429-8},
  number = {2},
  journal = {Journal of Scientific Computing},
  publisher = {Springer Science and Business Media LLC},
  author = {Busto,  Saray and Dumbser,  Michael and Escalante,  Cipriano and Favrie,  Nicolas and Gavrilyuk,  Sergey},
  year = {2021},
  month = mar 
}
```

#### Data Format

The CSV file contains two columns: x-coordinate and wave height for comparison with wave-over-Gaussian bathymetry simulations.

#### Usage

In order to access the data, execute:

```julia
x_data, h_data = data_busto()
```

The function returns:
- `x_data`: Spatial coordinates along the domain
- `h_data`: Wave heights at the corresponding spatial locations

### Mitsotakis et al. Reference Data (`Mitsotakis_et_al/`)

Numerical reference data for wave reflection studies from:

```bibtex
@article{Mitsotakis2016,
  title = {A modified Galerkin/finite element method for the numerical solution of the Serre‐Green‐Naghdi system},
  volume = {83},
  ISSN = {1097-0363},
  url = {http://dx.doi.org/10.1002/fld.4293},
  DOI = {10.1002/fld.4293},
  number = {10},
  journal = {International Journal for Numerical Methods in Fluids},
  publisher = {Wiley},
  author = {Mitsotakis,  D. and Synolakis,  C. and McGuinness,  M.},
  year = {2016},
  month = sep,
  pages = {755–778}
}
```

#### Data Format

The directory contains different CSV files for Figure 6a and 6b from the Mitsotakis et al paper at different times. Each CSV file contains x,y coordinate pairs representing wave profiles at different times and conditions.

#### Usage

For figure 6 data, execute:

```julia
times, all_data = data_mitsotakis(i)
```

Where `i` is either 1 or 2 for different experimental scenarios (either A = 0.075 or A = 0.65). The function returns:
- `times`: Array of time values for the different datasets
- `all_data`: Array of datasets, each containing x,y coordinate pairs for wave profiles

### Favre Wave Data (`Favre/`)

Experimental data for Favre waves taken from:

```bibtex
@article{Chassagne_Filippini_Ricchiuto_Bonneton_2019,
  title = {Dispersive and dispersive-like bores in channels with sloping banks},
  volume = {870},
  DOI = {10.1017/jfm.2019.287},
  journal = {Journal of Fluid Mechanics},
  author = {Chassagne, R. and Filippini, A. G. and Ricchiuto, M. and Bonneton, P.},
  year = {2019},
  pages = {595–616}
}
```

#### Data Format

The directory contains multiple text files (`dh{ε}t{t}_FNPF.txt`) representing wave height measurements at different amplitudes (ε = 0.1, 0.2, 0.3) and time instances (t = 50, 60, 70). Each file contains two columns: x-coordinate and wave height.

#### Usage

In order to access the data, execute:

```julia
all_data = data_favre()
```

The function returns a dictionary where keys are tuples `(ε, t)` and values are the corresponding data arrays with x-coordinates and wave heights.

### Froude Number Data (`Froude/`)

Experimental data for Froude number studies from Favre and Treske:

```bibtex
@book{Favre1935,
  author = {Favre, Henri},
  title = {Étude théorique et expérimentale des ondes de translation dans les canaux découverts},
  publisher = {Dunod},
  address = {Paris},
  year = {1935},
   url={https://api.semanticscholar.org/CorpusID:126909361}
}

@article{Treske1994,
  author = {Treske, A.},
  title = {Undular bores (Favre-waves) in open channels -- Experimental studies},
  journal = {Journal of Hydraulic Research},
  volume = {32},
  number = {3},
  pages = {355--370},
  year = {1994},
  month = {May},
  doi = {10.1080/00221689409498738}
}
```

#### Data Format

The directory contains four text files:
- `Favre_amplmax_100.txt` and `Favre_amplmax_200.txt`: Favre experimental data
- `Treske_amplmax_80.txt` and `Treske_amplmax_160.txt`: Treske experimental data

Each file contains two columns: Froude number and maximum amplitude.

#### Usage

In order to access the data, execute:

```julia
data_favre_100, data_favre_200, data_treske_80, data_treske_160 = data_froude()
```

The function returns four arrays containing Froude number and maximum amplitude data from the respective experiments.
