subroutine HilmRhoH(cilm, ba, grid, lmax, nmax, mass, r0, rho, gridtype, w, plx, zero, filter_type, filter_deg, lmax_calc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will compute the next estimate of Moho coefficients 
!	given an initial estimate of the moho in a gridded data file. This
!	is simply Equation 18 in Wieczorek and Phillips (1998), modified to 
!	take into account lateral variations in density as from eq. 30 of 
!	Wieczorek (2007). Note that the degree-0 topography must be included 
!	in the gridded relief. Note that the array plx is optional, and should 
!	not be precomputed when memory is an issue (i.e., lmax>360).
!
!	Calling Parameters:
!		IN
!			ba:		Bouguer Anomaly spherical harmonic coefficients.
!			grid:		Initial estimate of Moho refief, gridded according
!					to a call to Makegrid.
!			lmax:		Maxmimum spherical harmonic degree to compute.
!			nmax:		Order of potential coefficient expansion.
!			mass:		Mass of planet.	
!			R0:		The radius that the coefficients are referenced
!					to. 		
!			rho:		Laterally varying density contrast between the mantle and crust.
!			gridtype	1 = Gauss-Legendre quadrature grid corresponding to LMAX.
!					2 = N by N Driscoll and Healy grid corresponding to LMAX.
!					3 = N by 2N Driscoll and Healy grid corresponding to LMAX.
!		OUT
!			cilm:		Estimate of Moho relief spherical harmonic coefficients, 
!					with dimensions (2, lmax+1, lmax+1).
!			
!		OPTIONAL		
!			w:		Gauss-Legendre points used in integrations (for GRIDTYPE=1, 
!					determined from a call to PreCompute).
!			zero:		Array of dimension lmax+1 that contains the latitudinal
!					gridpoints used in the Gauss-Legendre quadrature integration
!					scheme (For GRIDTYPE=1). Only needed if plx is not included.
!			plx:		Input array of Associated Legendre Polnomials computed
!					at the Gauss points (for GRIDTYPE=1, determined from a call to
!					PreCompute). If this is not included, then the optional
!					array zero MUST be inlcuded.
!			filter_deg:	If specified, each interation will be filtered according to
!					equation 18 in Wieczorek and Phillips (1998), where the value
!					of filter_deg corresponds to the spherical harmonic degree where 
!					the filter is 1/2.
!			filter_type:	If filter_deg is specified, this must be as well. 
!					A value of (0) corresponds to none, (1) corresponds to the minimum amplitude filter in 
!					Wieczorek and Phillips (1998), whereas as value of (2) corresponds
!					to a minimum curvature filter.
!			lmax_calc:	Maximum degree to compute spherical harmonic expansions up to.
!
!
!	All units assumed to be SI.
!
!	Dependencies:		NGLQSH, SHExpandGLQ, SHExpandDH, wl, wl_curv, MakeGridGLQ, MakeGridDH
!
!	Written by Mark Wieczorek 2003
!		September 3, 2005. Modifed so that the array plx is now optional.
!		May 19, 2010. Modified to use Dricoll and Healy grids.
!		May 15, 2012. Modified to include lateral variations in density. Also modified
!				first spherical harmonic expansion to only compute the degree 0 term.
!
!	Copyright (c) 2005-2010, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: NGLQSH, SHExpandGLQ, SHExpandDH, Wl, WlCurv, MakeGridGLQ, MakeGridDH

	implicit none
		
	real*8, intent(out) :: 			cilm(:,:,:)
	real*8, intent(in) ::			ba(:,:,:), grid(:,:), mass, r0, rho(:,:)
	real*8, intent(in), optional ::		plx(:,:), zero(:), w(:)
	integer, intent(in) ::			lmax, nmax, gridtype
	integer, intent(in), optional :: 	filter_type, filter_deg, lmax_calc
	real*8 ::				prod, pi, d, depth, filter(lmax+1)
	real*8, allocatable ::			cilmn(:, :, :), grid2(:,:)
	integer  ::				j, l, n, nlong, nlat, astat(2), lmax_out, lmax_calc2, n_out

	pi = acos(-1.0d0)
	
	if (present(lmax_calc)) then
		if (lmax_calc > lmax .or. lmax_calc < 0) then
			print*, "Error -- HilmRhoH"
			print*, "LMAX_CALC must be less than or equal to LMAX."
			print*, "LMAX = ", lmax
			print*, "LMAX_CALC = ", lmax_calc
			stop
		endif
		lmax_calc2 = lmax_calc
	else
		lmax_calc2 = lmax
	endif

	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax_calc2+1 .or. size(cilm(1,1,:)) < lmax_calc2+1) then
		print*, "Error --- HilmRhoH"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax_calc2
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	elseif (size(ba(:,1,1)) < 2 .or. size(ba(1,:,1)) < lmax_calc2+1 .or. size(ba(1,1,:)) < lmax_calc2+1) then
		print*, "Error --- HilmRhoH"
		print*, "BA must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax_calc2
		print*, "Input dimension is ", size(ba(:,1,1)), size(ba(1,:,1)), size(ba(1,1,:))
		stop
	endif
	
	if (gridtype == 1) then
		nlat = NGLQSH(lmax) ! lmax+1
		nlong = 2*lmax+1
	elseif (gridtype == 2) then 
		nlat = 2*lmax+2
		nlong = nlat
	elseif (gridtype == 3) then
		nlat = 2*lmax+2
		nlong = 2*nlat
	else
		print*, "Error --- HilmRhoH"
		print*, "GRIDTYPE must be 1 (Gauss-Legendre Quadrature), 2 (Driscoll-Healy NxN) or 3 (Driscoll-Healy Nx2N)"
		stop
	endif
	
	
	if(size(grid(1,:)) < nlong .or. size(grid(:,1)) <  nlat .or. size(rho(1,:)) < nlong .or. size(rho(:,1)) <  nlat) then
		print*, "Error --- HilmRhoH"
		print*, "GRID and RHO must be dimensioned as (NLAT, NLONG) where"
		print*, "NLAT = ", nlat
		print*, "NLONG = ", nlong
		print*, "Input dimension of GRID is ", size(grid(1,:)), size(grid(:,1))
		print*, "Input dimension of RHO is ", size(rho(1,:)), size(rho(:,1))
		stop
	endif
	
	if (present(w)) then
		if (gridtype /= 1) then
			print*, "Error --- HilmRhoH"
			print*, "W can only be present when GRIDTYPE is 1 (Gauss-Legendre Quadrature)"
			print*, "GRIDTYPE = ", gridtype
			stop
		elseif(size(w) < lmax+1) then
			print*, "Error --- HilmRhoH"
			print*, "W must be dimensioned as (LMAX+1) where LMAX is ", lmax
			print*, "Input dimension is ", size(w)
			stop
		endif
	endif

	if (present(zero)) then
		if (gridtype /= 1) then
			print*, "Error --- HilmRhoH"
			print*, "ZERO can only be present when GRIDTYPE is 1 (Gauss-Legendre Quadrature)"
			print*, "GRIDTYPE = ", gridtype
			stop
		elseif (size(zero) < lmax + 1) then
			print*, "Error --- HilmRhoH"
			print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
			print*, "Input dimension is ", size(zero)
			stop
		endif
	endif
	
	if (present(plx)) then
		if (gridtype /= 1) then
			print*, "Error --- HilmRhoH"
			print*, "PLX can only be present when GRIDTYPE is 1 (Gauss-Legendre Quadrature)"
			print*, "GRIDTYPE = ", gridtype
			stop
		elseif (size(plx(:,1)) < lmax+1 .or. size(plx(1,:)) < (lmax+1)*(lmax+2)/2) then
			print*, "Error --- HilmRhoH"
			print*, "PLX must be dimensioned as (LMAX+1, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
			print*, "Input dimension is ", size(plx(:,1)), size(plx(1,:))
			stop
		endif
	endif

	if (present(w)) then
		if (present(plx) .and. present(zero)) then
			print*, "Error --- HilmRhoH"
			print*, "If the optional parameter W is present, then only one of PLX and ZERO must also be present."
			stop
		endif
		if (.not.present(plx) .and. .not.present(zero)) then
			print*, "Error --- HilmRhoH"
			print*, "If the optional parameter W is present, then one of PLX or ZERO must also be present."
			stop
		endif
	endif
	
	if(present(filter_type)) then
		if (filter_type /=0 .and. filter_type /=1 .and. filter_type /=2) then
			print*, "Error --- HilmRhoH"
			print*, "FILTER_TYPE must be either 0 (none), 1 (minimum amplitude), or 2 (minimum curvature)"
			print*, "Input value is ", filter_type
			stop
		endif
		if (filter_type == 1 .or. filter_type == 2) then
			if (.not.present(filter_deg)) then
				print*, "Error --- HilmRhoH"
				print*, "If FILTER_TYPE is present and equal to 1 or 2, so mush FILTER_DEG."
				stop
			endif
		endif
	endif
	
	allocate(cilmn(2, lmax_calc2+1, lmax_calc2+1), stat = astat(1))
	allocate(grid2(nlat, nlong), stat = astat(2))
	if (astat(1) /= 0 .or. astat(2) /= 0) then
		print*, "Error --- HilmRhoH"
		print*, "Problem allocating arrays CILMN and GRID2", astat(1), astat(2)
		stop
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Do the expansions
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cilm = 0.0d0
	cilmn = 0.0d0
	
	! Determine reference radius of Moho relief 
	grid2(1:nlat,1:nlong) = grid(1:nlat,1:nlong)
	if (gridtype == 1) then
		if (present(plx)) then
			call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
				plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), norm = 1, csphase = 1, lmax_calc = 0)
		else 
			call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
				zero = zero(1:lmax+1), norm = 1, csphase = 1, lmax_calc = 0)
		endif
	elseif (gridtype == 2) then
		call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilmn, lmax_out, norm = 1, &
			sampling = 1, csphase = 1, lmax_calc = 0)
	elseif (gridtype == 3) then
		call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilmn, lmax_out, norm = 1, &
			sampling = 2, csphase = 1, lmax_calc = 0)
	endif

	d = cilmn(1,1,1)
	depth = r0-d
	print*, "Average depth of Moho (km) = ", depth/1.d3
	
		
	! calculate (rho*h)_00
	grid2(1:nlat,1:nlong) = (grid(1:nlat,1:nlong)-d) * rho(1:nlat,1:nlong)
	if (gridtype == 1) then
		if (present(plx)) then
			call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
				plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), norm = 1, csphase = 1, lmax_calc = 0)
		else 
			call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
				zero = zero(1:lmax+1), norm = 1, csphase = 1, lmax_calc = 0)
		endif
	elseif (gridtype == 2) then
		call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilmn, lmax_out, norm = 1, &
			sampling = 1, csphase = 1, lmax_calc = 0)
	elseif (gridtype == 3) then
		call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilmn, lmax_out, norm = 1, &
			sampling = 2, csphase = 1, lmax_calc = 0)
	endif

	cilm(1,1,1) = cilmn(1,1,1)
		
	! Calculate first term resulting from Bouguer anomaly	
		
	filter = 1.0d0
	if (present(filter_type) .and. present(filter_deg)) then
		do l=1, lmax_calc2
			if (filter_type == 1) then
				filter(l+1) = Wl(l, filter_deg, r0, d)
			elseif(filter_type == 2) then
				filter(l+1) = WlCurv(l, filter_deg, r0, d)
			endif
		enddo
	endif
	
	do l=1, lmax_calc2
		cilm(1:2,l+1,1:l+1) = filter(l+1)*ba(1:2,l+1,1:l+1) * mass * dble(2*l+1) * ( (r0/d)**l) / &
			(4.0d0*pi*d**2)
	enddo
		
	! calculate higher order terms
	
	do n=2, nmax
		
		grid2(1:nlat,1:nlong) = rho(1:nlat,1:nlong) * ((grid(1:nlat,1:nlong) - d)/d)**n

		if (gridtype == 1) then
			if (present(plx)) then
				call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
					plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), norm = 1, csphase = 1, lmax_calc = lmax_calc2)
			else
				call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
					zero = zero(1:lmax+1), norm = 1, csphase = 1, lmax_calc = lmax_calc2)
			endif
		elseif (gridtype == 2) then
			call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilmn, lmax_out, norm = 1, &
				sampling = 1, csphase = 1, lmax_calc = lmax_calc2)
		elseif (gridtype == 3) then
			call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilmn, lmax_out, norm = 1, &
				sampling = 2, csphase = 1, lmax_calc = lmax_calc2)
		endif
		
		do l = 1, lmax_calc2
			prod = 1.0d0
			do j=1, n
				prod = prod * dble(l+4-j)
			enddo
			prod = d * prod/( dble(l+3) * dble(fact(n)) )
			
			cilm(1:2,l+1,1:l+1) = cilm(1:2,l+1,1:l+1) - filter(l+1)*cilmn(1:2,l+1,1:l+1)*prod
		enddo
	enddo
	
	! make grid of (h*rho)_lm
	if (gridtype == 1) then
		if (present(plx)) then
			call MakeGridGLQ(grid2(1:nlat,1:nlong), cilm, lmax,  &
				plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), norm = 1, csphase = 1, lmax_calc = lmax_calc2)
		else
			call MakeGridGLQ(grid2(1:nlat,1:nlong), cilm, lmax,  &
				zero = zero(1:lmax+1), norm = 1, csphase = 1, lmax_calc = lmax_calc2)
		endif
	elseif (gridtype == 2) then
		call MakeGridDH(grid2(1:nlat,1:nlong), n_out, cilm, lmax, norm = 1, &
			sampling = 1, csphase = 1, lmax_calc = lmax_calc2)
	elseif (gridtype == 3) then
		call MakeGridDH(grid2(1:nlat,1:nlong), n_out, cilm, lmax, norm = 1, &
			sampling = 2, csphase = 1, lmax_calc = lmax_calc2)
	endif
	
	! convert grid
	
	grid2(1:nlat,1:nlong) = grid2(1:nlat,1:nlong) / rho(1:nlat,1:nlong)
	
	if (gridtype == 1) then
		if (present(plx)) then
			call SHExpandGLQ(cilm, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
				plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), norm = 1, csphase = 1, lmax_calc = lmax_calc2)
		else
			call SHExpandGLQ(cilm, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), &
				zero = zero(1:lmax+1), norm = 1, csphase = 1, lmax_calc = lmax_calc2)
		endif
	elseif (gridtype == 2) then
		call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilm, lmax_out, norm = 1, &
			sampling = 1, csphase = 1, lmax_calc = lmax_calc2)
	elseif (gridtype == 3) then
		call SHExpandDH(grid2(1:nlat,1:nlong), nlat, cilm, lmax_out, norm = 1, &
			sampling = 2, csphase = 1, lmax_calc = lmax_calc2)
	endif
	
	cilm(1,1,1) = d
	
	deallocate(cilmn)
	deallocate(grid2)	
	
	contains
	
		function fact(i)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	This function computes the factorial of an integer
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			implicit none
			integer ::	i, j, fact
	
			if (i == 0) then
				fact = 1	
			elseif (i .lt. 0) then
				print*, "Argument to FACT must be positive"
				stop
			else
				fact = 1
				do j = 1, i
					fact = fact * j
				enddo
			endif
			
		end function fact
	
end subroutine HilmRhoH

