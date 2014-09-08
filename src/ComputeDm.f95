subroutine ComputeDm(dllm, lmax, m, theta0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will compute the kernel D for the space-concentration 
!	problem of a spherical cap for a given spherical harmonic order. 
!	(When m /= 0, the eigenfunctions will not be isotropic). All terms 
!	are computed exactly using Gauss-Legendre quadrature. The eigenfunctions 
!	of this kernel are spherical harmonic coefficients normalized according
!	to the "geodesy" convention. The diagonal elements of Dllm approaches 
!	unity when theta0 appoaches 180 degrees.
!
!	Calling Parameters
!		IN
!			lmax:		Maximum spherical harmonic degree.
!			theta0:		Angular radius of spherical cap IN RADIANS.
!			m:		Angular order used in the concentration problem.
!		OUT
!			dllm:		Symmetric kernel of size lmax+1 by lmax+1.
!
!	Dependencies:	PreGLQ, PlmBar, NGLQ
!
!	Written by Mark Wieczorek (June 2004)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	use SHTOOLS, only: PreGLQ, PlmBar, NGLQ, PlmIndex

	implicit none
	
	real*8, intent(out) ::	dllm(:,:)
	real*8, intent(in) ::	theta0
	integer, intent(in) :: 	lmax, m
	real*8 ::	upper, zero(2*lmax+1), w(2*lmax+1), x
	real*8, allocatable ::	plm(:)
	integer ::	l, lp, i, n, astat

	if (abs(m) > lmax) then
		print*, "Error --- ComputeDm"
		print*, "M must be less than or equal to LMAX."
		print*, "Input values of LMAX and M are", lmax, m
		stop
	elseif (size(dllm(1,:)) < lmax+1 .or. size(dllm(:,1)) < lmax+1) then
		print*, "Error --- ComputeDm"
		print*, "DLLM must be dimensioned as (LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(dllm(1,:)), size(dllm(:,1))
		stop
	endif
	
	allocate(plm( (lmax+1)*(lmax+2)/2 ), stat = astat)
	if (astat /= 0) then
		print*, "Error --- ComputeDM"
		print*, "Problem allocating array PLM", astat
		stop
	endif

	dllm = 0.0d0

	n = NGLQ(2*lmax)
	
	upper = 1.0d0
	x = cos(theta0)
				
	Call PreGLQ(x, upper, n, zero, w)

	do i=1, n
		call PlmBar(plm, lmax, zero(i))
		
		do l=abs(m), lmax
			do lp = l, lmax
			
				dllm(l+1, lp+1) = dllm(l+1,lp+1) + w(i) * plm(PlmIndex(l,abs(m))) * plm(PlmIndex(lp,abs(m)))	
				
			enddo
		enddo
	enddo
	
	
	do l=abs(m), lmax
		do lp=l+1, lmax, 1
			dllm(lp+1,l+1) = dllm(l+1,lp+1)
		enddo
	enddo
	
	if (m == 0) then
		dllm = dllm / 2.0d0
	else
		dllm = dllm / 4.0d0
	endif 
	
	call PlmBar(plm, -1, zero(1)) ! deallocate memory
	
	deallocate(plm)

end subroutine ComputeDm

