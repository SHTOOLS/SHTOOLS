subroutine PreCompute(lmax, zero, w, plx, wisdom_file, norm, csphase, cnorm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will precompute several arrays and matrixes that are 
!	necessary for the Gauss-Legendre quadrature sherical harmonic
!	expansion routines. This is done only to speed up computations
!	in the case where more than one expansion is performed.
!	In addition, system-wide FFTW wisdom will be loaded, as well 
!	as an optional user wisdom file. Note, if memory is an issue
!	because of large spherical harmonic degrees, DO NOT PRECOMPUTE
!	THE ARRAY PLX!
!
!	Calling Parameters:
!		IN
!			lmax:		Maximum spherical harmonic degree to
!					be used in all following expansions.
!		OUT
!			zero:		An array of length (lmax + 1) which contains the Gauss points 
!					(i.e., zeros of P(l+1,m=0)).
!			w:		An array of length lmax + 1 which contains the weights used in the 
!					Gauss-Legendre quadrature routines.
!		OPTIONAL
!			wisdom_file	User wisdom file name.
!			plx:		Array of normalized associated Legendre Polnomials computed
!					at the Gauss points. The size of this array must be
!					(1:lmax+1, 1:(lmax+1)*(lmax+2)/2).
!			norm:		Normalization to be used when calculating Legendre functions
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			cnorm:		1: compute complex normalized Legendre functions.
!					0 (default): compute real normalized Legendre functions.
!
!
!	Dependencies:	PreGLQ, NGLQSH, PlmBar, PlmSchmidt, PLegendreA, PlmON, CSPHASE_DEFAULT, PlmIndex
!
!	Written by Mark Wieczorek 2003
!	September 4, 2005. The array plx was turned into an optional argument. This
!	is because when the spherical harmonic degree is very large (>360), it takes
!	up a lot of memory. When doing the spherical harmonic transforms, if this array
!	is not sent to the subroutine, it will be computed for each latitude. This 
!	takes place at the same speed as before. The only reason for precomputing this
!	is in case it it reused several times.
!	May 1, 2008. Added optional arguement CNORM to calculate complex normalized Legendre functions.
!
!	Copyright (c) 2005-2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PreGLQ, NGLQSH, PlmBar, PlmSchmidt, PLegendreA, PlmON, CSPHASE_DEFAULT, PlmIndex

	implicit none
	       	
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	zero(:), w(:)
	real*8, intent(out), optional ::	plx(:,:)
	integer, intent(in), optional ::	norm, csphase, cnorm
	character(*), intent(in), optional ::	wisdom_file
	character ::	system_wisdom(80)
	integer ::		n, i, isuccess, fn, istat, astat, phase, l, m, i_s, cnormin
	integer, save ::	first_call = 1
	real*8  :: 		upper, lower, pi
	real*8, allocatable :: 	 pl(:)
	
	
	if (size(zero) < lmax+1) then 
		print*, "Error --- PreCompute"
		print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(zero)
		stop
	elseif (size(w) < lmax+1) then 
		print*, "Error --- PreCompute"
		print*, "W must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(w)
		stop
	endif
	
	if (present(plx)) then
		if (size(plx(:,1)) < lmax +1 .or. size(plx(1,:)) < (lmax+1)*(lmax+2)/2) then
			print*, "Error --- PreCompute"
			print*, "PLX must be dimensioned as (LMAX+1, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
			print*, "Input array is dimensioned as ", size(plx(:,1)), size(plx(1,:))
			stop
		endif
	endif
	
	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error - PreCompute"
			print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
			print*, "Input value is ", norm
			stop
		endif
	endif
	
	if (present(csphase)) then
     		if (csphase /= -1 .and. csphase /= 1) then
     			print*, "PreCompute --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif
     	
     	allocate(pl(((lmax+2)*(lmax+1))/2), stat = astat)
     	if (astat /= 0) then
		print*, "Error --- PreCompute"
		print*, "Problem allocating array PL", astat
		stop
	endif

	if (present(cnorm)) then
		cnormin = cnorm
	else
		cnormin = 0
	endif
	
	pi = acos(-1.0d0)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Load FFTW wisdom files if this is the first time this routine is being called
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (first_call == 1) then
	
		system_wisdom = "/etc/fftw/wisdom"
		call dfftw_import_system_wisdom(isuccess)
	
		if (isuccess /= 0) then
			print*, "Using system-wide FFTW wisdom."
		endif

		if (present(wisdom_file)) then

			fn = 12
			open(unit=fn, file=wisdom_file, iostat = istat)
			
			if (istat /= 0 ) then
				print*, "Error --- PreCompute"
				print*, "Problem opening user defined FFTW wisdom file."
			endif
			
			call import_wisdom_from_file(isuccess, fn)
			if (isuccess == 0) then
				print*, "Local Wisdom file successfully read."
			else
				print*, "Local Wisdom file not read."
			endif

		endif
		
		first_call = 0
		
	endif
	
	
	upper = 1.0d0
	lower = -1.0d0			
	
	n = NGLQSH(lmax)	! lmax is the maxium degree of the equation we
				! are going to integrate. 
						
	call PreGLQ(lower, upper, n, zero, w)	! Determine Gauss Points and Weights.
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Next, fill the array plx with all of the associated legendre polynomials
	! 	evaluated at the Gauss points. The first index of this array corresponds
	! 	to the Gauss point, whereas the (lmax+2)*(lmax+1)/2 elements in the second column
	! 	correspond to all of the assiated legendre polynomials, ordered according to
	! 	index = (l+1)*l/2 + m + 1.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (present(plx)) then
		
		do i=1, (n+1)/2
			if (present(norm)) then
				if (norm == 1) call PlmBar(pl, lmax, zero(i), csphase = phase, cnorm=cnormin)
				if (norm == 2) call PlmSchmidt(pl, lmax, zero(i), csphase = phase, cnorm=cnormin)
				if (norm == 3) call PLegendreA(pl, lmax, zero(i), csphase = phase)
				if (norm == 4) call PlmON(pl, lmax, zero(i), csphase = phase, cnorm=cnormin)
			else
				call PlmBar(pl, lmax, zero(i), csphase = phase, cnorm=cnormin)
			endif
			
			if (i == (n+1)/2 .and. mod(n,2) /= 0) then
				plx(i,1:((lmax+1)*(lmax+2))/2) = pl(1:((lmax+1)*(lmax+2))/2)
			else
				i_s = n+1-i
				plx(i,1:((lmax+1)*(lmax+2))/2) = pl(1:((lmax+1)*(lmax+2))/2)
				
				do l=0, lmax
					do m=0, l
						plx(i_s,PlmIndex(l,m)) = plx(i,PlmIndex(l,m)) * (-1.0d0)**(l-m)
					enddo
				enddo
				
			endif
					
		enddo
		
		if (present(norm)) then
			if (norm == 1) call PlmBar(pl, -1, zero(1), csphase = phase, cnorm=cnormin)
			if (norm == 2) call PlmSchmidt(pl, -1, zero(1), csphase = phase, cnorm=cnormin)
			if (norm == 4) call PlmON(pl, -1, zero(1), csphase = phase, cnorm=cnormin)
		else
			call PlmBar(pl, -1, zero(1), csphase = phase, cnorm=cnormin)
		endif
	
	endif
	
	deallocate(pl)

end subroutine PreCompute

