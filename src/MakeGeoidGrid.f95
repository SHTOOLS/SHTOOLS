subroutine MakeGeoidGrid(geoid, cilm, lmax, r0pot, GM, PotRef, omega, r, gridtype, &
			order, nlat, nlong, interval, lmax_calc, a, f)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subrouine will calculate the height to an equipotential surface with
!	respect to a spherical reference radius R for a body rotating with angular rotation
!	rate OMEGA. This method uses a first-, second-, or third-order Taylor series expansion 
!	of the potential, evaluated at radius R. The output gridype is specificied by GRIDTYPE, 
!	which can be either Cartesian 2D, Guass-Legendre quadrature, or Driscoll and Healy. 
!	If the optional parameters A and F are specified, the geoid will be referenced to a 
!	flattened ellispoid with semimajor axis A and flattening F.
!	
!	Note that this routine is only strictly valid when the geoid is above the surface!
!	To calculated the height of the geoid when it is below the surface, one would
!	need to know the density structure of the planet. Furthermore, the calculation of the 
!	potential (and its derivatives) is only strictly valid when R lies above the maximum radius
!	of the planet.
!
!	Latitude is geocentric latitude.
!
!	Calling Parameters
!		IN
!			cilm		Spherical harmonic potential coefficients.
!			lmax		Maximum spherical harmonic degree. For gridtype=1, this
!					is simply the maximum degree used in evaluating the 
!					spherical harmonic coefficients. For the other grids,
!					this sets the spacing of the output grid.
!			r0pot		Reference radius of the potential coefficients.
!			GM		GM associated with the potential coefficients.
!			PotRef		Reference potential of the geoid.
!			omega		Angular rotation rate used when calculating the geoid
!					on in the rotating coordinate system.
!			r		Reference sperical radius to calculate the geoid heights.
!			gridtype	1 = Gauss-Legendre quadrature grid corresponding to LMAX.
!					2 = N by N Driscoll and Healy grid corresponding to LMAX.
!					3 = N by 2N Driscoll and Healy grid corresponding to LMAX.
!					4 = 2D Cartesian using MakeGrid2D.
!			order 		Order of the Taylor expansion. Either 1, 2, or 3.
!		OUT
!			geoid		Gridded values of the geoid, in meters, referenced to the spherical
!					radius R.
!			nlat, nlong	Number of latitude and longitude points for the grid.
!		OPTIONAL, IN
!			interval	Grid spacing of the output grid in DEGREES. Used only when
!					GRIDTYPE = 1.
!			lmax_calc	For GRIDTYPE 2, 3, and 4, this specifies the maximum spherical
!					harmonic degree to evaluate the function to.
!			a		Semimajor axis of the reference flattened ellipsoid.
!			f		Flattening of the reference ellipsoid.
!			
!			
!
!	Dependencies :  MakeGrid2D, MakeGridGLQ, MakeGridDH, NGLQSH, PreCompute
!
!	Written by Mark Wieczorek (October 2005)
!
!	May 31, 2006 -  Modified so that the Condon-Shortley phase is never used, regardless of 
!			what the default value of CSPHASE_DEFAULT is.
!	June 23, 2007 - Completely rewritten. Can now use 4 different types of grids, and expanded
!			to third order.
!	April 2009 -	Changed the order of the gridtypes
!	June 2013 - 	Fixed bug that could give rise to runtime errors when the dimensions of CILM were less than lmax.
!
!	Copyright (c) 2005-2007, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: MakeGrid2D, MakeGridGLQ, MakeGridDH, NGLQSH, PreCompute

	implicit none
	
	real*8, intent(out) ::	geoid(:,:)
	real*8, intent(in) ::	cilm(:,:,:), r0pot, GM, r, PotRef, omega
	integer, intent(in) :: lmax, order, gridtype
	integer, intent(in), optional :: lmax_calc
	integer, intent(out) :: nlat, nlong
	real*8, intent(in), optional ::	interval, a, f
	real*8 ::	pi, r_ex, lat
	integer :: 	l, nlat1, nlong1, lmax_comp, astat, n, i, astat1, astat2
	real*8, allocatable :: grida(:,:), gridb(:,:), gridc(:,:), gridd(:,:), zero(:), w(:), qq(:,:), pp(:,:), uu(:,:), &
				cilm1(:,:,:), cilm2(:,:,:)

	if ( (present(f) .and. .not. present(a)) .or. (present(a) .and. .not. present(f)) ) then
		print*, "Error --- MakeGeoidGrid"
		print*, "Both F and A must be specified."
		print*, "A, ", present(a)
		print*, "F, ", present(f)
		stop
	endif
	
	if (present(lmax_calc)) then
		if (lmax_calc > lmax) then
			print*, "Error --- MakeGeoidGrid"
			print*, "LMAX_CALC must be less than or equal to LMAX."
			print*, "LMAX_CALC = ", lmax_calc
			print*, "LMAX = ", lmax
			stop
		endif
	endif

	if (gridtype < 1 .or. gridtype > 4) then
		print*, "Error --- MakeGeoidGrid"
		print*, "GRIDTYPE must be either (1) Cartesian 2D, (2) Gauss-Legendre quadrature,"
		print*, "(3) N by N Driscoll and Healy, or (4) N by 2N Driscoll and Healy."
		print*, "Input value = ", gridtype
		stop
	endif
	
	if (gridtype == 4 .and. .not. present(interval)) then
		print*, "Error --- MakeGeoidGrid"
		print*, "If GRIDTYPE = 4 (2D Cartesian), the optional parameter INTERVAL"
		print*, "must be specified."
		stop
	endif
	
	if (order > 3 .or. order < 1) then
		print*, "Error --- MakeGeoidGrid"
		print*, "ORDER must be 1, 2, or 3."
		stop
	endif

	
	if (gridtype == 1) then
		nlong = 2*lmax + 1
		nlat = NGLQSH(lmax)	! nlat = lmax+1
	elseif (gridtype == 2) then
		nlat = 2*lmax+2
		nlong = nlat
	elseif (gridtype == 3) then
		nlat = 2*lmax+2
		nlong = 2*nlat
	elseif (gridtype == 4) then
		nlat = 180.0/interval + 1
		nlong = 360./interval + 1
	endif
		
	if (size(geoid(:,1)) < nlat .or. size(geoid(1,:)) < nlong ) then
		print*, "Error --- MakeGeoidGrid"
		print*, "GEOID must be dimensioned as (180/INTERVAL+1, 360/INTERVAL+1) where INTERVAL is ", interval
		print*, "Input array has dimension ", size(geoid(:,1)), size(geoid(1,:))
		stop
	elseif (size(cilm(:,1,1)) < 2) then
		print*, "Error --- MakeGeoidGrid"
		print*, "CILM must be dimensioned as (2, *, *)."
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif

	if (present(lmax_calc)) then
		lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, lmax_calc)
	else
		lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)
	endif
	
	if (gridtype == 1) then
		allocate(zero(lmax+1), stat = astat)
		if (astat /=0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array ZERO, ", astat
			stop
		endif
		allocate(w(lmax+1), stat = astat)
		if (astat /=0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array W, ", astat
			stop
		endif
		call PreCompute(lmax, zero, w, norm=1, csphase=1)
	endif
	
	allocate(cilm1(2,lmax+1,lmax+1), stat = astat1)
	allocate(cilm2(2,lmax+1,lmax+1), stat = astat2)
	if (astat1 /= 0 .or. astat2 /= 0) then
		print*, "Error --- MakeGeoidGrid"
		print*, "Problem allocating arrays CILM1 and CILM2", astat1, astat2
		stop
	endif

	pi = acos(-1.0d0)
	
	cilm1(1:2,1:lmax_comp+1, 1:lmax_comp+1) = cilm(1:2,1:lmax_comp+1, 1:lmax_comp+1)
	cilm1(1,1,1) = 1.0d0	! Make sure that the degree-0 term is included

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 
	!	Solve the equation a + b x + c x**2 + d x**3 = 0
	!	for each grid point.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
	! Create grid A.
	
	do l=0, lmax_comp
		cilm2(1:2,l+1,1:l+1) = (GM/r) * cilm1(1:2,l+1,1:l+1) * (r0pot/r)**l	
	enddo
	! add rotation and subtract reference field
	cilm2(1,1,1) = cilm2(1,1,1) + (omega*r)**2 / 3.0d0 - PotRef
	cilm2(1,3,1) = cilm2(1,3,1) - (omega*r)**2 / (3.0d0*sqrt(5.0d0))
	
	allocate(grida(nlat, nlong), stat = astat)
	if (astat /= 0) then
		print*, "Error --- MakeGeoidGrid"
		print*, "Problem allocating array GRIDA", astat
		stop
	endif

	select case(gridtype)
		case(1); call MakeGridGLQ(grida, cilm2, lmax, zero=zero, norm=1, csphase=1, lmax_calc=lmax_comp)
		case(2); call MakeGridDH(grida, n, cilm2, lmax, norm=1, csphase=1, sampling=1, lmax_calc = lmax_comp)
		case(3); call MakeGridDH(grida, n, cilm2, lmax, norm=1, csphase=1, sampling=2, lmax_calc = lmax_comp)
		case(4); call MakeGrid2D(grida, cilm2, lmax_comp, interval, nlat1, nlong1, norm=1, csphase=1)
	end select

	! Create Grid B	
	
	do l=0, lmax_comp
		cilm2(1:2,l+1,1:l+1) = -(GM/r**2) * dble(l+1) * cilm1(1:2,l+1,1:l+1) * (r0pot/r)**l
	enddo
	! add rotational terms
	cilm2(1,1,1) = cilm2(1,1,1) + 2.0d0 * r * omega**2 / 3.0d0
	cilm2(1,3,1) = cilm2(1,3,1) - 2.0d0 * r * omega**2 / (3.0d0*sqrt(5.0d0))
	
	allocate(gridb(nlat, nlong), stat = astat)
	if (astat /= 0) then
		print*, "Error --- MakeGeoidGrid"
		print*, "Problem allocating array GRIDB", astat
		stop
	endif
	select case(gridtype)
		case(1); call MakeGridGLQ(gridb, cilm2, lmax, zero=zero, norm=1, csphase=1, lmax_calc=lmax_comp)
		case(2); call MakeGridDH(gridb, n, cilm2, lmax, norm=1, csphase=1, sampling=1, lmax_calc = lmax_comp)
		case(3); call MakeGridDH(gridb, n, cilm2, lmax, norm=1, csphase=1, sampling=2, lmax_calc = lmax_comp)
		case(4); call MakeGrid2D(gridb, cilm2, lmax_comp, interval, nlat1, nlong1, norm=1, csphase=1)
	end select


	! Create Grid C
	
	if (order == 2 .or. order == 3) then
		
		do l=0, lmax_comp
			cilm2(1:2,l+1,1:l+1) = GM/(2.0d0*r**3) * dble(l+1) * dble(l+2) * cilm1(1:2,l+1,1:l+1) * (r0pot/r)**l
		enddo
		! add rotational terms
		cilm2(1,1,1) = cilm2(1,1,1) + omega**2 / 3.0d0
		cilm2(1,3,1) = cilm2(1,3,1) - omega**2 / (3.0d0*sqrt(5.0d0))

		allocate(gridc(nlat, nlong), stat = astat)
		if (astat /= 0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array GRIDC", astat
			stop
		endif
	
		select case(gridtype)
			case(1); call MakeGridGLQ(gridc, cilm2, lmax, zero=zero, norm=1, csphase=1, lmax_calc=lmax_comp)
			case(2); call MakeGridDH(gridc, n, cilm2, lmax, norm=1, csphase=1, sampling=1, lmax_calc = lmax_comp)
			case(3); call MakeGridDH(gridc, n, cilm2, lmax, norm=1, csphase=1, sampling=2, lmax_calc = lmax_comp)
			case(4); call MakeGrid2D(gridc, cilm2, lmax_comp, interval, nlat1, nlong1, norm=1, csphase=1)
		end select
	
	endif
	
	
	! Create Grid D
	
	if (order == 3) then
		
		do l=0, lmax_comp
			cilm2(1:2,l+1,1:l+1) = -GM/(6.0d0*r**4) * dble(l+1) * dble(l+2) * dble(l+3) * &
						cilm1(1:2,l+1,1:l+1) * (r0pot/r)**l
		enddo

		allocate(gridd(nlat, nlong), stat = astat)
		if (astat /= 0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array GRIDD", astat
			stop
		endif

		select case(gridtype)
			case(1); call MakeGridGLQ(gridd, cilm2, lmax, zero=zero, norm=1, csphase=1, lmax_calc=lmax_comp)
			case(2); call MakeGridDH(gridd, n, cilm2, lmax, norm=1, csphase=1, sampling=1, lmax_calc = lmax_comp)
			case(3); call MakeGridDH(gridd, n, cilm2, lmax, norm=1, csphase=1, sampling=2, lmax_calc = lmax_comp)
			case(4); call MakeGrid2D(gridd, cilm2, lmax_comp, interval, nlat1, nlong1, norm=1, csphase=1)
		end select
	
	endif
	
	! solve equation
	
	if (order == 1) then
	
		geoid(1:nlat,1:nlong) = -grida(1:nlat,1:nlong) / gridb(1:nlat,1:nlong)
		
	elseif (order == 2) then
		
		geoid(1:nlat, 1:nlong) = (-gridb(1:nlat, 1:nlong) - &
		 		sqrt(gridb(1:nlat, 1:nlong)**2 - 4.0d0*gridc(1:nlat,1:nlong) * &
				grida(1:nlat, 1:nlong)) ) / (2.0d0*gridc(1:nlat,1:nlong))
		
	elseif (order == 3) then
		allocate(qq(nlat, nlong), stat = astat)
		if (astat /= 0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array QQ", astat
			stop
		endif
		allocate(pp(nlat, nlong), stat = astat)
		if (astat /= 0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array PP", astat
			stop
		endif
		allocate(uu(nlat, nlong), stat = astat)
		if (astat /= 0) then
			print*, "Error --- MakeGeoidGrid"
			print*, "Problem allocating array UU", astat
			stop
		endif
		
		pp(1:nlat,1:nlong) = gridb(1:nlat,1:nlong)/gridd(1:nlat,1:nlong) - &
			((gridc(1:nlat,1:nlong)/gridd(1:nlat,1:nlong))**2) / 3.0d0
		qq = grida(1:nlat,1:nlong)/gridd(1:nlat,1:nlong) + &
			2.0d0*((gridc(1:nlat,1:nlong)/gridd(1:nlat,1:nlong))**3) / 27.0d0 - &
			9.0d0*gridc(1:nlat,1:nlong)*gridb(1:nlat,1:nlong)/(gridd(1:nlat,1:nlong)**2) / 27.0d0
		uu(1:nlat,1:nlong) = (qq(1:nlat,1:nlong)/2.0d0 + &
			sqrt((qq(1:nlat,1:nlong)**2)/4.0d0 + (pp(1:nlat,1:nlong)**3)/27.0d0) )**(1.0d0/3.0d0)
		geoid(1:nlat,1:nlong) = pp(1:nlat,1:nlong)/3.0d0/uu(1:nlat,1:nlong) - &
			uu(1:nlat,1:nlong) - gridc(1:nlat,1:nlong)/gridd(1:nlat,1:nlong)/3.0d0
		
		deallocate(qq)
		deallocate(pp)
		deallocate(uu)
		
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Reference geoid to a flattened ellipsoid
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (present(a) .and. present(f)) then
		do i = 1, nlat
		
			if (gridtype == 4) then
				lat = 90.0d0 - dble(i-1)*interval
			elseif (gridtype == 1) then
				lat = asin(zero(i)) * 180.0d0 / pi
			elseif (gridtype == 2 .or. gridtype == 3) then
				lat = 90.0d0 - 180.0d0 * dble(i-1)/dble(nlat)
			endif
			
			r_ex = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
			(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
			r_ex = sqrt(r_ex)

			geoid(i,1:nlong) = geoid(i,1:nlong) + r - r_ex
		enddo
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Clean up
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (gridtype == 2) then
		deallocate(zero)
		deallocate(w)
	endif
	deallocate(grida)
	deallocate(gridb)
	if (order == 2) deallocate(gridc)
	if (order == 3) then
		deallocate(gridc)
		deallocate(gridd)
	endif
	
	deallocate(cilm1)
	deallocate(cilm2)
	

end subroutine MakeGeoidGrid
