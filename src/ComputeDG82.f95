subroutine ComputeDG82(DG82, lmax, m, theta0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will compute the kernel of Grunbaum et al. 1982
!	that commutes with the space concentration kernel. This 
!	kernel is tridiagonal and has simple expressions for 
!	its elements. Note that this kernel is multiplied by -1 in
!	comparison to Grunbaum et al. in order that the eigenvectors 
!	will be listed in the same order as the eigenvectors of the
!	equivalent space concentration matrix Dllm. While the eigenfunctions
!	of this kernel correspond to the eigenvalues of Dllm, the eigenvalues
!	do NOT!
!
!	Calling Parameters
!		IN
!			lmax:		Maximum spherical harmonic degree.
!			theta0:		Angular radius of spherical cap IN RADIANS.
!			m:		Angular order to concentration problem.
!		OUT
!			D0G82:		Symmetric tridiagonal kernel with a maximum size of lmax+1 by lmax+1.
!
!	Dependencies:	none
!
!	Written by Mark Wieczorek (June 2004)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none

	real*8, intent(out) ::	DG82(:,:)
	real*8, intent(in) ::	theta0
	integer, intent(in) :: 	lmax, m
	real*8 ::		x
	integer ::		i, n


	n = lmax + 1 - abs(m)
	
	if(n<1) then
		print*, "Error --- ComputeDG82"
		print*, "abs(M) must be less than or equal to LMAX"
		print*, "Input values of l and m are ", lmax, m
		stop
	elseif (size(DG82(1,:)) < n .or. size(DG82(:,1)) < n) then
		print*, "Error --- ComputeDG82"
		print*, "DG82 must be dimensioned as (LMAX-abs(M)+1,LMAX-abs(M)+1) where LMAX and M are ", lmax, m
		print*, "Input array is dimensioned as ", size(DG82(1,:)), size(DG82(:,1))
		stop
	endif

	DG82 = 0.0d0

	x = cos(theta0)
	
		
	DG82(1,1) = x * dble(1+m) * dble(m)
	
	do i=2, n, 1
		DG82(i,i) = x * dble(i+m) * dble(i+m-1)
		DG82(i,i-1) = - sqrt(dble(i-1+m)**2 - dble(m)**2) * &
			( dble(i-1+m)**2 - dble(lmax+1)**2 ) / &
			sqrt(4.0d0*dble(i-1+m)**2 - 1.0d0)
		DG82(i-1,i) = DG82(i,i-1)
	enddo
		
end subroutine ComputeDG82

