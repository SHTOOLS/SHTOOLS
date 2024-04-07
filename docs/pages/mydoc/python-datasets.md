---
title: "Datasets"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-datasets.html
summary: pyshtools provides easy access to many research-grade and historical datasets.
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

## Datasets

pyshtools provides easy access to many research-grade gravity, shape, and magnetic field datasets of the terrestrial planets, moons, and asteroids. To load a dataset, it is only necessary to call the relevant method from the *datasets* submodule as in these examples:
```
    hlm = pysh.datasets.Venus.VenusTopo719()  # Venus shape
    clm = pysh.datasets.Earth.EGM2008()  # Earth gravity
    glm = pysh.datasets.Earth.WDMAM2_800()  # Earth magnetic field
    clm = pysh.datasets.Moon.GRGM1200B()  # Gravity of the Moon
```

When accessing a dataset, the file will first be downloaded from the original source using [pooch](https://www.fatiando.org/pooch/latest/) and then stored in the pyshtools subdirectory of the user's cache directory (if it had not been done previously). The file hash will be verified to ensure that it has not been modified, and the file will then be used to initialize and return an `SHCoeffs`, `SHGravCoeffs` or `SHMagCoeffs` class instance. In most cases the files are stored in their original form: If they need to be decompressed or unzipped, this is done on the fly when they are used. Only when zip archives contain several files are the files stored in unzipped form.

The coefficients can be read up to a maximum specified degree by providing the optional variable `lmax`. For some datasets (such as the shape of Mars), several different models of differing resolution exist, and only the smallest file that is necessary based on the user input value of `lmax` will be downloaded. For IGRF magnetic field coefficients, the year of the output coefficients can be specified by the optional argument `year` (the default is 2020).

The following is the list of implemented datasets. Additional older or deprecated datasets can be found in the `historical` module for each body.

### Mercury

| Dataset | Description |
| ---------- | ----------- |
| USGS_SPG_shape | 5759 degree and order spherical harmonic shape model of Mercury based on the USGS stereo-photogrammetric DTM (Maia 2024). |
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
| SWARM_MLI_2D_0501 | Degree 133 magnetic field model of the Earth's lithosphere that is based largely on satellite data (Thébault et al. 2013). Though this model is based largely on data from the SWARM mission, data from the CHAMP mission  and some ground-based measurements are also used. |
| NGDC_720_V3 | Degree 740 magnetic field model of the Earth's lithosphere that was compiled from satellite, marine, aeromagnetic and ground-based magnetic surveys (Maus 2010). |
| WDMAM2_800 | Degree 800 model of the Earth's lithospheric magnetic field that is based on a worldwide compilation of near-surface magnetic data. This is the second version of the World Digital Magnetic Anomaly Map (WDMAM) (Lesur et al. 2016). |


### The Moon

| Dataset | Description |
| ---------- | ----------- |
| LOLA_shape_pa | 5759 degree and order spherical harmonic model of the shape of Earth's Moon in a principal axis coordinate system (Wieczorek 2024). |
| LOLA_shape | 5759 degree and order spherical harmonic model of the shape of Earth's Moon in mean Earth/polar axis coordinate system (Wieczorek 2024). |
| GRGM900C | GSFC 900 degree and order spherical harmonic model of the gravitational potential of the Moon. This model applies a Kaula constraint for degrees greater than 600 (Lemoine et al. 2014). |
| GRGM1200B | GSFC 1200 degree and order spherical harmonic model of the gravitational potential of the Moon (Goossens et al. 2020). This model applies a Kaula constraint for degrees greater than 600. |
| GRGM1200B_RM1_1E0 | GSFC 1200 degree and order spherical harmonic model of the gravitational potential of the Moon (Goossens et al. 2020). This model uses a rank-minus-1 constraint based on gravity from surface topography for degrees greater than 600 with a value of lambda equal to 1. |
| GL0900D | JPL 900 degree and order spherical harmonic model of the gravitational potential of the Moon (Konopliv et al. 2014). This model applies a Kaula constraint for degrees greater than 700. |
| GL1500E | JPL 1500 degree and order spherical harmonic model of the gravitational potential of the Moon (Konopliv et al. 2014). This model applies a Kaula constraint for degrees greater than 700. |
| T2015_449 | 449 degree and order spherical harmonic model of the magnetic potential of the Moon. This model was used in Wieczorek (2018) and is a spherical harmonic expansion of the global magnetic field model of Tsunakawa et al. (2015). |
| Ravat2020 | 450 degree and order spherical harmonic model of the magnetic potential of the Moon. This model is based on using magnetic monopoles with an L1-norm regularization. |


### Mars

| Dataset | Description |
| ---------- | ----------- |
| MOLA_shape | 5759 degree and order spherical harmonic model of the shape of the planet Mars (Wieczorek 2024). |
| GMM3 | GSFC 120 degree and order spherical harmonic model of the gravitational potential of Mars (Genova et al. 2016). This model applies a Kaula constraint for degrees greater than 90. |
| GMM3_RM1_1E0 | GSFC 150 degree and order spherical harmonic model of the gravitational potential of Mars (Goossens et al. 2017). This model uses the same data as GMM3, but with a rank-minus-1 constraint based on gravity from surface topography for degrees greater than 50 with a value of lambda equal to 1. |
| MRO120F |JPL 120 degree and order spherical harmonic model of the gravitational potential of Mars (Konopliv et al. 2020). This model applies a Kaula constraint for degrees greater than 80. |
| Langlais2019 | 134 degree and order spherical harmonic model of the magnetic potential of Mars (Langlais et al. 2019). This model makes use of data from MGS MAG, MGS ER and MAVEN MAG. |
| Morschhauser2014 | 110 degree and order spherical harmonic model of the magnetic potential of Mars (Morschhauser et al. 2014). |


### (1) Ceres

| Dataset | Description |
| ---------- | ----------- |
| DLR_SPG_shape | 5399 degree and order spherical harmonic shape model of asteroid (1) Ceres based on the DLR stereo-photogrammetric DTM (Wieczorek 2024). |
| JPL_SPC_shape | 1023 degree and order spherical harmonic shape model of asteroid (1) Ceres based on the JPL stereo-photoclinometric DTM (Wieczorek 2024). |
| CERES18D | JPL 18 degree and order spherical harmonic model of the gravitational potential of asteroid (1) Ceres (Konopliv et al. 2018). |


### (4) Vesta

| Dataset | Description |
| ---------- | ----------- |
| DLR_SPG_shape | 5759 degree and order spherical harmonic shape model of asteroid (4) Vesta based on the DLR stereo-photogrammetric DTM (Wieczorek 2024). |
| VESTA20H | JPL 20 degree and order spherical harmonic model of the gravitational potential of asteroid (4) Vesta (Konopliv et al. 2014). |


### (433) Eros

| Dataset | Description |
| ---------- | ----------- |
| NLR_shape | 719 degree and order spherical harmonic shape model of asteroid (433) Eros based on laser altimeter data (Wieczorek 2024). |
| SPC_shape | 511 degree and order spherical harmonic shape model of asteroid (433) Eros based on a stereo-photoclinometry shape model (Wieczorek 2024). |
| JGE15A01 | JPL 15 degree and order spherical harmonic model of the gravitational potential of asteroid (433) Eros (Miller et al. 2002). |


### Jupiter

| Dataset | Description |
| ---------- | ----------- |
| Kaspi2023_gravity | Degree 40 and order 0 spherical harmonic model of the gravitational potential of Jupiter (Kaspi et al. 2023). |


### Io (Jupiter)

| Dataset | Description |
| ---------- | ----------- |
| Anderson2001 | Degree and order 2 spherical harmonic model of the gravitational potential of Io (Anderson et al. 2001). |


### Europa (Jupiter)

| Dataset | Description |
| ---------- | ----------- |
| Anderson1998 | Degree and order 2 spherical harmonic model of the gravitational potential of Europa (Anderson et al. 1998). |


### Ganymede (Jupiter)

| Dataset | Description |
| ---------- | ----------- |
| Ganymede2022 | Degree and order 5 spherical harmonic model of the gravitational potential of Ganymede (Gomez Casajus et al. 2022). |
| Anderson1996_1 | Degree and order 2 spherical harmonic model of the gravitational potential of Ganymede from the first Galileo encounter (Anderson et al. 1996). |
| Anderson1996_2 | Degree and order 2 spherical harmonic model of the gravitational potential of Ganymede from the second Galileo encounter (Anderson et al. 1996). |


### Callisto (Jupiter)

| Dataset | Description |
| ---------- | ----------- |
| Anderson2001 | Degree and order 2 spherical harmonic model of the gravitational potential of Callisto (Anderson et al. 2001). |


### Saturn

| Dataset | Description |
| ---------- | ----------- |
| Jacobson2022_gravity | Degree 12 and order 0 spherical harmonic model of the gravitational potential of Saturn (Jacobson 2022). |


### Titan (Saturn)

| Dataset | Description |
| ---------- | ----------- |
| Corlies2017_shape | Degree and order 8 spherical harmonic shape model of Saturn's moon Titan based on radar altimeter, SAR topography, and radar stereo-photogrammetric data from the Cassini mission (Corlies et al. 2014). |
| Mitri2014_shape | Degree and order 6 spherical harmonic shape model of Saturn's moon Titan based on radar altimeter and SAR topography data from the Cassini mission (Mitri et al. 2014). |
| Durante2019_gravity | Degree and order 5 spherical harmonic model of the gravity field of Saturn's moon Titan (Durante et al. 2019). |


### Enceladus (Saturn)

| Dataset | Description |
| ---------- | ----------- |
| JPL_SPC_shape | Degree and order 1023 spherical harmonic model of the shape of Saturn's moon Enceladus based on the JPL stereo-photoclinometric DTM (Wieczorek 2024). |
| Iess2014_gravity | Degree and order 3 spherical harmonic model of the gravitational field of Saturn's moon Enceladus (Iess et al. 2014). |
| Park2024_gravity | Degree and order 3 spherical harmonic model of the gravitational field of Saturn's moon Enceladus (Park et al. 2024). |


### Uranus

| Dataset | Description |
| ---------- | ----------- |
| Jacobson2014_gravity | Degree 6 and order 0 spherical harmonic model of the gravitational potential of Uranus (Jacobson 2014). |


### Neptune

| Dataset | Description |
| ---------- | ----------- |
| Jacobson2009_gravity | Degree 4 and order 0 spherical harmonic model of the gravitational potential of Neptune (Jacobson 2009). |
