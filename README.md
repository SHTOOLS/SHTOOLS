	<h2>What is SHTOOLS?</h2>
		
	<p>SHTOOLS is an archive of fortran 95 based software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere. </p>

	<h2>What makes SHTOOLS different?</h2>

	<p class="nomarginbot">While several collections of code exist for working with data expressed in spherical harmonics, this one is unique for several reasons:</p>

	<ul>	
		<li>It can accommodate any standard normalization of the spherical harmonic functions ("geodesy" 4&pi; normalized,  Schmidt semi-normalized, orthonormalized, and unnormalized).</li>
		
		<li>Both real and complex spherical harmonics are supported.</li>

		<li>Spherical harmonic transforms are calculated by exact quadrature rules using either (1) the sampling theorem of <i>Driscoll and Healy</i> (1994) where data are equally sampled (or spaced) in latitude and longitude, or (2) Gauss-Legendre quadrature. A least squares inversion routine for irregularly sampled data is included.</li>

		<li>One can choose to use or exclude the Condon-Shortley phase factor of (-1)<sup>m</sup> with the associated Legendre functions.</li>

		<li>The spherical harmonic transforms are proven to be accurate to approximately degree 2800, corresponding to a spatial resolution of better than 4 arc minutes.</li>

		<li>Routines are included for performing localized multitaper spectral analyses.</li>

		<li>Routines are included for performing standard gravity and magnetic field calculations, such as computation of the geoid and the determination of the potential associated with finite-amplitude topography.</li>

		<li>The routines are fast. Spherical harmonic transforms and reconstructions take on the order of 1 second for bandwidths less than 600 and about 3 minutes for bandwidths close to 2800.</li>
	</ul>
		
	<h2>How do I use SHTOOLS?</h2>

	<p>SHTOOLS is invoked by standard subroutine and function calls in your fortran 90/95 program. For some routines, it will be necessary to have installed the freely available Fourier transform package <a href="http://www.fftw.org">FFTW</a>, and the linear algebra packages <a href="http://www.netlib.org/lapack/">LAPACK</a> and <a href="http://www.netlib.org/blas/">BLAS</a>.</p>

	<h2>How much does it cost?</h2>

	<p class="extramarginbot">The SHTOOLS software package is free. Two free Fortran 95 compilers exist (gfortran as part of <a href="http://gcc.gnu.org/">GCC</a> and <a href="http://www.g95.org/">g95</a>). SHTOOLS can be modified and distributed according to the revised BSD license.</p>
	