---
title: "Constants and datasets"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-datasets-constants.html
summary: pyshtools provides easy access to many research-grade datasets and common physical constants.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

## Constants

The *constants* subpackage defines physical constants related to the terrestrial planets and moons. Each constant is an instance of an [astropy](http://docs.astropy.org/en/stable/constants/index.html) `Constant` class, which has the attributes `name`, `value`, `uncertainty`, `unit`, and `reference`.

The pyshtools constants are organized primarily by planet: Mercury, Venus, Earth, Moon, and Mars. Each planet has several attributes, such as the the mean planetary radius `r`, `mass`, and flattening `f`. To see all information about an individual constant, it is only necessary to use the print function:
```python
In [1]: print(pysh.constants.Mars.r)
  Name   = Mean radius of Mars
  Value  = 3389500.0
  Uncertainty  = 0.0
  Unit  = m
  Reference = MarsTopo2600: Wieczorek, M. A. (2015). Gravity and topography of the terrestrial planets. In T. Spohn & G. Schubert (Eds.), Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153-193). Oxford, Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.
```
To use the value of a constant in a calculation, it is only necessary to access its `value` attribute:

```python
In [2]: 2 * pysh.constants.Mars.r.value
6779000.0
```
Physical constants from the *Committee on Data for Science and Technology* are provided in the submodule `codata`. A few of these (such as `G` and `mu0`) are referenced in the main constants namespace.

## Datasets

pyshtools provides easy access to many research-grade gravity, topography, and magnetic field datasets of the terrestrial planets. To load a dataset, it is only necessary to call the relevant method from the *datasets* submodule as in these examples:
```
    hlm = pysh.datasets.Venus.VenusTopo719()  # Venus shape
    clm = pysh.datasets.Earth.EGM2008()  # Earth gravity
    glm = pysh.datasets.Earth.WDMAM2_800()  # Earth magnetic field
    clm = pysh.datasets.Moon.GRGM1200B()  # Gravity of the Moon
```

When accessing a dataset, the file will first be downloaded from the original source using [pooch](https://www.fatiando.org/pooch/latest/) and then stored in the pyshtools subdirectory of the user's cache directory (if it had not been done previously). The file hash will be verified to ensure that it has not been modified, and the file will then be used to initialize and return an `SHCoeffs`, `SHGravCoeffs` or `SHMagCoeffs` class instance. In most cases the files are stored in their original form: If they need to be decompressed or unzipped, this is done on the fly when they are used. Only when zip archives contain several files are the files stored in unzipped form.

The coefficients can be read up to a maximum specified degree by providing the optional variable `lmax`. For IGRF magnetic field coefficients, the year of the output coefficients can be specified by the optional argument `year` (the default is 2020).

The following is the list of implemented datasets:

### Mercury

| Dataset | Description |
| ---------- | ----------- |
| GTMES150 | GSFC 150 degree and order spherical harmonic model of the shape of the planet Mercury. |
| JGMESS160A | JPL 160 degree and order spherical harmonic model of the gravitational potential of Mercury (Konopliv et al. 2020). This model applies a Kaula law constraint to all degrees. |
| JGMESS160A_ACCEL | JPL 160 degree and order spherical harmonic model of the gravitational potential of Mercury (Konopliv et al. 2020). This model applies a surface acceleration constraint that is based on the degree strength given by the coefficient covariance. |
| JGMESS160A_TOPOSIG | JPL 160 degree and order spherical harmonic model of the gravitational potential of Mercury (Konopliv et al. 2020). This model applies a constraint similar to the Kaula constraint except that the uncertainty is given by the corresponding magnitude of the gravity derived from topography. |
| GGMES100V08 | GSFC 100 degree and order spherical harmonic model of the gravitational potential of Mercury (Genova et al. 2019). This model applies a Kaula constraint to all degrees. |


### Venus

| Dataset | Description |
| ---------- | ----------- |
| VenusTopo719 | 719 degree and order spherical harmonic model of the shape of the planet Venus (Wieczorek 2015). |
| MGNP180U | JPL 180 degree and order spherical harmonic model of the gravitational potential of Venus (Konopliv et al. 1999). |


### Earth

| Dataset | Description |
| ---------- | ----------- |
| Earth2012 | Shape and topography (with respect to mean sea level) of Earth expanded to degree and order 2160 (Hirt et al. 2012). |
| Earth2014 | Topography of Earth with respect to mean sea level, expanded to degree and order 2160 (Hirt and Rexer 2015). |
| EGM2008 | Degree 2190 model of the Earth's gravity field in a tide-free system. This model is based on data from altimetry, ground-based measurements, and the satellite GRACE (Pavlis et al. 2012). |
| EIGEN_6C4 | Degree 2190 model of the Earth's gravity field in a tide-free system (Förste et al. 2014). This model is based on data from altimetry, ground-based measurements, and the satellites GOCE, GRACE, and LAGEOS. |
| GGM05C | Degree 360 model of the Earth's gravity field in a zero-tide system (Ries et al. 2016). This model is based on data from altimetry, ground-based measurements, and the satellites GOCE and GRACE. |
| GOCO06S | Degree 300 model of the Earth's gravity field in a zero-tide system (Kvas et al. 2019). This model is based solely on satellite data. |
| EIGEN_GRGS_RL04_MEAN_FIELD | Degree 300 model of the Earth's gravity field in a tide-free system (Lemoine et al. 2019). This model is based solely on satellite data. |
| XGM2019E | Degree 2190 model of the Earth's gravity field in a zero-tide system (Zingerle et al. 2020). This combined model is based on data from altimetry, ground-based measurements, topography, and satellites. |
| IGRF_13 | Degree 13 time-variable model of the Earth's main magnetic field that is valid between 1900 and 2020 (Thébault et al. 2015). |
| SWARM_MLI_2D_0501 | Degree 133 magnetic field model of the Earth's lithsophere that is based largely on satellite data (Thébault et al. 2013). Though this model is based largely on data from the SWARM mission, data from the CHAMP mission  and some ground-based measurements are also used. |
| NGDC_720_V3 | Degree 740 magnetic field model of the Earth's lithsophere that was compiled from satellite, marine, aeromagnetic and ground-based magnetic surveys (Maus 2010). |
| WDMAM2_800 | Degree 800 model of the Earth's lithospheric magnetic field that is based on a worldwide compilation of near-surface magnetic data. This is the second version of the World Digital Magnetic Anomaly Map (WDMAM) (Lesur et al. 2016). |


### The Moon

| Dataset | Description |
| ---------- | ----------- |
| MoonTopo2600p | 2600 degree and order spherical harmonic model of the shape of Earth's Moon in a principal axis coordinate system (Wieczorek 2015). |
| GRGM900C | GSFC 900 degree and order spherical harmonic model of the gravitational potential of the Moon. This model applies a Kaula constraint for degrees greater than 600 (Lemoine et al. 2014). |
| GRGM1200B | GSFC 1200 degree and order spherical harmonic model of the gravitational potential of the Moon (Goossens et al. 2020). This model applies a Kaula constraint for degrees greater than 600. |
| GRGM1200B_RM1_1E0 | GSFC 1200 degree and order spherical harmonic model of the gravitational potential of the Moon (Goossens et al. 2020). This model uses a rank-minus-1 constraint based on gravity from surface topography for degrees greater than 600 with a value of lambda equal to 1. |
| GL0900D | JPL 900 degree and order spherical harmonic model of the gravitational potential of the Moon (Konopliv et al. 2014). This model applies a Kaula constraint for degrees greater than 700. |
| GL1500E | JPL 1500 degree and order spherical harmonic model of the gravitational potential of the Moon (Konopliv et al. 2014). This model applies a Kaula constraint for degrees greater than 700. |
| T2015_449 | 449 degree and order spherical harmonic model of the magnetic potential of the Moon. This model was used in Wieczorek (2018) and is a spherical harmonic expansion of the global magnetic field model of Tsunakawa et al. (2015). |
| Ravat2020 | 450 degree and order spherical harmonic model of the magnetic potential of the Moon. This model is based on using magnetic monopoles with an L1-norm regularisation. |


### Mars

| Dataset | Description |
| ---------- | ----------- |
| MarsTopo2600 | 2600 degree and order spherical harmonic model of the shape of the planet Mars (Wieczorek 2015). |
| GMM3 | GSFC 120 degree and order spherical harmonic model of the gravitational potential of Mars (Genova et al. 2016). This model applies a Kaula constraint for degrees greater than 90. |
| GMM3_RM1_1E0 | GSFC 150 degree and order spherical harmonic model of the gravitational potential of Mars (Goossens et al. 2017). This model uses the same data as GMM3, but with a rank-minus-1 constraint based on gravity from surface topography for degrees greater than 50 with a value of lambda equal to 1.|
| MRO120D |JPL 120 degree and order spherical harmonic model of the gravitational potential of Mars (Konopliv et al. 2016). This model applies a Kaula constraint for degrees greater than 80. |
| Langlais2019 | 134 degree and order spherical harmonic model of the magnetic potential of Mars (Langlais et al. 2019). This model makes use of data from MGS MAG, MGS ER and MAVEN MAG. |
| Morschhauser2014 | 110 degree and order spherical harmonic model of the magnetic potential of Mars (Morschhauser et al. 2014). |


### Vesta

| Dataset | Description |
| ---------- | ----------- |
| VESTA20H | JPL 20 degree and order spherical harmonic model of the gravitational potential of (4) Vesta (Konopliv et al. 2014). |


### Ceres

| Dataset | Description |
| ---------- | ----------- |
| CERES18D | JPL 18 degree and order spherical harmonic model of the gravitational potential of (1) Ceres (Konopliv et al. 2018). |
