# SHAdmitCorr 

Calculate the admittance and correlation spectra of two real functions.

# Usage

call SHAdmitCorr (`gilm`, `tilm`, `lmax`, `admit`, `corr`, `admit_error`)

# Parameters

`gilm` : input, real\*8, dimension (2, `lmaxg`+1, `lmaxg`+1)
:   The real spherical harmonic coefficients of the function `G`.

`tilm` : input, real\*8, dimension (2, `lmaxt`+1, `lmaxt`+1)
:   The real spherical harmonic coefficients of the function `T`.

`lmax` : input, integer
:   The maximum spherical harmonic degree that will be calculated for the admittance and correlation spectra. This must be less than or equal to the minimum of `lmaxg` and `lmaxt`.

`admit` : output, real\*8, dimension (`lmax`+1)
:   The admittance function, which is equal to `Sgt/Stt`.

`corr` : output, real\*8, dimension (`lmax`+1)
:   The degree correlation function, which is equal to `Sgt/sqrt(Sgg Stt)`.

`admit_error` : output, optional, real\*8, dimension (`lmax`+1)
:   The uncertainty of the admittance function, assuming that `gilm` and `tilm` are related by a linear isotropic transfer function, and that the lack of correlation is a result of uncorrelated noise.

# Description

`SHAdmitCorr` will calculate the admittance and correlation spectra associated with two real functions expressed in real spherical harmonics. The admittance is defined as `Sgt/Stt`, where `Sgt` is the cross-power spectrum of two functions `G` and `T`. The degree-correlation spectrum is defined as `Sgt/sqrt(Sgg Stt)`, which can possess values between -1 and 1. 

If the optional argument `admit_error` is specified, then the error of the admittance will be calculated by assuming that `G` and `T` are related by a linear isotropic transfer function:` Gilm = Ql Tilm + Nilm`, where `N` is noise that is uncorrelated with the topography. It is important to note that the relationship between two fields is often not described by such an isotropic expression.

# See also

[shpowerspectrum](shpowerspectrum.html), [shcrosspowerspectrum](shcrosspowerspectrum.html)
