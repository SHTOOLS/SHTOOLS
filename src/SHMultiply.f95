subroutine SHMultiply(shout, sh1, lmax1, sh2, lmax2, precomp, norm, csphase)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will multiply two spherical harmonic fields which are expressed
!	up to maximum spherical harmonic degrees lmax1 and lmax2. The output spherical 
!	harmonic coefficients will have a maximum spherical harmonic degree equal 
!	to lmax1 + lmax2.
!
!	Calling Parameters:
!		IN
!			sh1		Spherical harmonic field with maximum spherical harmonic
!					degree lmax1.
!			sh2		Spherical harmonic field with maximun spherical harmonic
!					degree lmax2.
!			lmax1		Maximum spherical harmonic degree of sh1.
!			lmax2		Maximum spherical harmonic degree of sh2. 
!		OUT	
!			shout		Spherical harmonic expansion of spatial multiplication of
!					sh1 and sh2, with a maximum spherical harmonic degree of 
!					lmax1 + lmax2.
!		OPTIONAL (IN)
!			precomp		If 1, the array plx will be precomputed when calling
!					the subroutine PreCompute. If 0 (default), then this array will 
!					not be precomputed.
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			norm:		Normalization to be used when calculating Legendre functions
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!
!	Dependencies:	PreCompute, MakeGridGLQ, SHExpandGLQ, CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek 2004
!	September 3, 2005. Modified to include the optional argument precomp which determines whether
!	the array plx will be precomputed in the subroutine PreCompute.
!	July 25, 2012. Modified to save arrays ZERO and W to reduce redundancy with multiple calls.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PreCompute, MakeGridGLQ, SHExpandGLQ, CSPHASE_DEFAULT

	implicit none
		
	real*8, intent(out) ::	shout(:,:,:)
	real*8, intent(in) ::	sh1(:,:,:), sh2(:,:,:)
	integer, intent(in) ::	lmax1, lmax2
	integer, intent(in), optional ::	precomp, norm, csphase
	integer :: 	lmaxout, phase, mnorm, astat(2), nlat, nlong
	real*8, allocatable, save ::	zero(:), w(:)
	integer, save ::	first = 1, lmaxout_last = 0
	real*8, allocatable :: 	grid1glq(:,:), grid2glq(:,:), plx(:,:)


	if (size(sh1(:,1,1)) < 2 .or. size(sh1(1,:,1)) < lmax1+1 .or. size(sh1(1,1,:)) < lmax1+1) then
		print*, "Error --- SHMultiply"
		print*, "SHI must be dimensioned as (2, LMAX1+1, LMAX1+1) where LMAX1 is", lmax1
		print*, "Input array is dimensioned ", size(sh1(:,1,1)),  size(sh1(1,:,1)), size(sh1(1,1,:)) 
		stop
	elseif (size(sh2(:,1,1)) < 2 .or. size(sh2(1,:,1)) < lmax2+1 .or. size(sh2(1,1,:)) < lmax2+1) then
		print*, "Error --- SHMultiply"
		print*, "SH2 must be dimensioned as (2,LMAX2+1, LMAX2+1) where LMAX2 is", lmax2
		print*, "Input array is dimensioned ", size(sh2(:,1,1)),  size(sh2(1,:,1)), size(sh2(1,1,:)) 
		stop
	elseif (size(shout(:,1,1)) < 2 .or. size(shout(1,:,1)) < lmax1+lmax2+1 .or. size(shout(1,1,:)) < lmax1+lmax2+1) then
		print*, "Error --- SHMultiply"
		print*, "SHOUT must be dimensioned as (2, LMAX1+LMAX2+1, LMAX1+LMAX2+1) where LMAX1 and LMAX2 are", lmax1, lmax2
		print*, "Input array is dimensioned ", size(shout(:,1,1)),  size(shout(1,:,1)), size(shout(1,1,:)) 
		stop
	endif

	if (present(csphase)) then
     		if (csphase /= -1 .and. csphase /= 1) then
     			print*, "Error --- SHMultiply"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif
     	
     	if (present(precomp)) then
     		if (precomp /= 1 .and. precomp /=0) then
     			print*, "Error --- SHMultiply"
     			print*, "PRECOMP must be either 0 (do not precompute PLX) or 1 (precompute PLX)."
     			print*, "Input value is ", precomp
     			stop
     		endif
     	endif
     	
	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error - SHExpandGLQ"
			print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
			print*, "Input value is ", norm
			stop
		endif
		mnorm = norm
	else
		mnorm = 1
	endif
	
	lmaxout = lmax1 + lmax2
	nlat = lmax1+lmax2+1
	nlong = 2*(lmax1+lmax2)+1
	
	if (first == 1) then
		first = 0
		lmaxout_last = lmaxout
		allocate(zero(lmaxout+1), stat = astat(1))
		allocate(w(lmaxout+1), stat = astat(2))
		
		if (sum(astat(1:2)) /= 0) then
			print*, "Error --- SHMultiply"
			print*, "Problem allocating arrays ZERO and W", astat(1), astat(2)
			stop
		endif
		
		call PreCompute(lmaxout, zero, w, csphase = phase, norm = mnorm)
	endif
	
	if (lmaxout /= lmaxout_last) then
		lmaxout_last = lmaxout
		deallocate(zero)
		deallocate(w)
		allocate(zero(lmaxout+1), stat = astat(1))
		allocate(w(lmaxout+1), stat = astat(2))
		
		if (sum(astat(1:2)) /= 0) then
			print*, "Error --- SHMultiply"
			print*, "Problem allocating arrays ZERO and W", astat(1), astat(2)
			stop
		endif
		
		call PreCompute(lmaxout, zero, w, csphase = phase, norm = mnorm)
		
	endif
	
	allocate(grid1glq(nlat, nlong), stat = astat(1))
	allocate(grid2glq(nlat, nlong), stat = astat(2))
	if (sum(astat(1:2)) /= 0) then
		print*, "Error --- SHMultiply"
		print*, "Problem allocating arrays GRID1GLQ and GRID2GLQ", astat(1), astat(2)
		stop
	endif
	
	if (present(precomp)) then
		if (precomp == 0) then
	
			call MakeGridGLQ(grid1glq, sh1(1:2,1:lmax1+1, 1:lmax1+1), lmaxout, zero = zero, csphase = phase, norm = mnorm)

			call MakeGridGLQ(grid2glq, sh2(1:2,1:lmax2+1, 1:lmax2+1), lmaxout, zero = zero, csphase = phase, norm = mnorm)
		
			grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong)*grid2glq(1:nlat,1:nlong)
	
			call SHExpandGLQ(shout, lmaxout, grid1glq, w, zero = zero, csphase = phase, norm = mnorm)
			
		else
			allocate ( plx(lmax1+lmax2+1, (lmax1+lmax2+1)*(lmax1+lmax2+2)/2 ), stat=astat(1) )
			if (astat(1) /=0) then 
				print*, "Error --- SHMultiply"
				print*, "Problem allocating array PLX", astat(1)
				stop
			endif
			
			call PreCompute(lmaxout, zero, w, plx = plx, csphase = phase, norm = mnorm)
	
			call MakeGridGLQ(grid1glq, sh1(1:2,1:lmax1+1, 1:lmax1+1), lmaxout, plx = plx, csphase = phase, norm = mnorm)

			call MakeGridGLQ(grid2glq, sh2(1:2,1:lmax2+1, 1:lmax2+1), lmaxout, plx = plx, csphase = phase, norm = mnorm)
		
			grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong)*grid2glq(1:nlat,1:nlong)
	
			call SHExpandGLQ(shout, lmaxout, grid1glq, w, plx = plx, csphase = phase, norm = mnorm)
		
			deallocate (plx)
		
		endif

	else
	
		call MakeGridGLQ(grid1glq, sh1(1:2,1:lmax1+1, 1:lmax1+1), lmaxout, zero = zero, csphase = phase, norm = mnorm)

		call MakeGridGLQ(grid2glq, sh2(1:2,1:lmax2+1, 1:lmax2+1), lmaxout, zero = zero, csphase = phase, norm = mnorm)
		
		grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong)*grid2glq(1:nlat,1:nlong)
			
		call SHExpandGLQ(shout, lmaxout, grid1glq, w, zero = zero, csphase = phase, norm = mnorm)
		
	endif	
	
	deallocate(grid1glq)
	deallocate(grid2glq)

end subroutine SHMultiply

