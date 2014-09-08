real*8 function SHConfidence(l_conf, r)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will calculate the probability that two sets of spherical
!	harmonic coefficients, which possess a correlation coefficients r, are 
!	linearly correlated at a given degree. This is calucated according to
!	equation A7 in Pauer et al. (2007), which is from Eckhard 1984, and MathWorld 
!	http://mathworld.wolfram.com/CorrelationCoefficientBivariateNormalDistribution.html.
!
!	CALLING PARAMETERS
!		INPUT
!			l_conf		degree to calculate confidence levels.
!			r		the correlation coefficient of two sets
!					of spherical harmonic coeficients at degree l_conf.
!		OUTPUT
!			cl		The confidence level.
!
!	Written by Mark Wieczorek (April 18, 2007)
!
!	Copyright (c) 2007, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: r
	integer, intent(in) :: l_conf
	real*8 :: prod
	integer:: l, i

	SHConfidence = abs(r)
	prod = 1.0d0
	
	do l=2, l_conf, 1
		i = l-1
		prod = prod * dble(2*i-1) / dble(2*i)
		SHConfidence = SHConfidence + prod * abs(r) * (1.0d0 - r**2)**(l-1)
	enddo

end function SHConfidence
