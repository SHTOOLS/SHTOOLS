subroutine SHCilmToVector(cilm, vector, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert an array of spherical harmonic
!	coefficients Cilm to an ordered 1D vector.
!
!	Calling Parameters
!		IN
!			cilm	Input spherical harmonic coefficients with dimension
!				cilm(2, lmax+1, lmax+1).
!			lmax	Maximum spherical harmonic degree of the input coefficients
!		OUT
!			vector	1D vector of ordered spherical harmonic coefficients with
!				dimension (lmax+1)**2. The ordering is described in YilmIndex.
!
!
!	Written by Mark Wieczorek (August 2009)
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) ::	cilm(:,:,:)
	real*8, intent(out) ::	vector(:)
	integer, intent(in) ::	lmax
	integer ::		i, l, m
	
	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- SHCilmToVector"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX + 1)."
		print*, "LMAX = ", lmax
		print*, "Dimension of CILM = ", size(cilm(:,1,1)),  size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if (size(vector(:)) < (lmax+1)**2) then
		print*, "Error --- SHCilmToVector"
		print*, "VECTOR must be have dimension (LMAX+1)**2."
		print*, "LMAX = ", lmax
		print*, "Dimension of VECTOR = ", size(vector(:))
		stop
	endif
	
	if (lmax < 0) then
		print*, "Error --- SHCilmToVector"
		print*, "LMAX must be positive."
		print*, "LMAX = ", lmax
		stop
	endif
	
	i=0
	do l=0, lmax
		do m=0, l
			i = i + 1
			vector(i) = cilm(1, l+1, m+1)
		enddo
		
		do m=1, l
			i = i + 1
			vector(i) = cilm(2, l+1, m+1)
		enddo
	enddo
	
end subroutine SHCilmToVector


subroutine SHVectorToCilm(vector, cilm, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert a 1D ordered vector of spherical harmonic
!	coefficients to a 3D array Cilm(i, l, m).
!
!	Calling Parameters
!		IN
!			vector	1D vector of ordered spherical harmonic coefficients with
!				dimension (lmax+1)**2. The ordering is described in YilmIndex.
!			lmax	Maximum spherical harmonic degree of the input coefficients
!		OUT
!			cilm	Output spherical harmonic coefficients with dimension
!				cilm(2, lmax+1, lmax+1).
!			
!
!	Written by Mark Wieczorek (August 2009)
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(out) ::	cilm(:,:,:)
	real*8, intent(in) ::	vector(:)
	integer, intent(in) ::	lmax
	integer ::		k, i, l, m
	
	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- SHVectorToCilm"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX + 1)."
		print*, "LMAX = ", lmax
		print*, "Dimension of CILM = ", size(cilm(:,1,1)),  size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if (size(vector(:)) < (lmax+1)**2) then
		print*, "Error --- SHVectorToCilm"
		print*, "VECTOR must be have dimension (LMAX+1)**2."
		print*, "LMAX = ", lmax
		print*, "Dimension of VECTOR = ", size(vector(:))
		stop
	endif
	
	if (lmax < 0) then
		print*, "Error --- SHVectorToCilm"
		print*, "LMAX must be positive."
		print*, "LMAX = ", lmax
		stop
	endif
	
	l = 0
	m = 0
	i = 1
	
	cilm(1,1,1) = vector(1)
	
	i = 2
	
	do k=2, (lmax+1)**2
		m = m + 1
		if (m > l .and. i == 1) then
			i = 2
			m = 1
		elseif (m > l .and. i == 2) then
			l = l + 1
			m = 0
			i = 1
		endif
		
		cilm(i,l+1,m+1) = vector(k)
	
	enddo
	
end subroutine SHVectorToCilm

