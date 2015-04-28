program MarsCrustalThickness
!-------------------------------------------------------------------------------
!
!	This program will compute the relief along the crust-mantle interface that
!	is needed to explain the martian gravity field. This crustal thickness model 
!	is "anchored" by chosing a minimum crustal thickness.
!	
!	In this program, the maximum spherical harmonic degree that is used
!	when calculating the relief raised to the nth power is the input spherical 
!   harmonic degree. This is not entirely correct, as the relief raised to the 
!   nth power generates a spherical harmonic field up to lmax*n. However, as it 
!   is also not correct to assume that the topo and gravity field are 
!   bandlimited with a maximum spherical harmonic degree of lmax, this 
!   approximation is probably ok, espeically considering that the iterations
!	are filtered.
!
!	In order to improve stability when iterating for the Moho relief, the following iterative
!	scheme is used:
!
!		h3 = (h2+h1)/2
!		h4 = f(h3)
!
!	where "h" represents the moho relife, and "f" represents the SHTOOLS 
!   funcion "BAtoHilm".
!
!	Copyright (c) 2015, Mark A. Wieczorek
!	All rights reserved.
!
!-------------------------------------------------------------------------------
	use SHTOOLS
	use PLANETSCONSTANTS

	implicit none
	real*8 :: 	r0, mass, rho_crust, rho_mantle, pi, gm, thinnest(2), &
			d, t0, delta, delta_max, thick_delta, param(8), grav, &
			r_grav, interval, rref, timein, timeout, max_thick, min_thick
	real*8, allocatable ::	topo_grid(:,:), moho_grid(:,:), moho_grid2(:,:), &
	        moho_grid3(:,:), temp_grid(:,:), topo_c(:,:,:), moho_c(:,:,:), &
	        bc(:,:,:), ba(:,:,:), cilm(:,:,:), pot(:,:,:), misfit(:,:,:)
	integer ::	l, m, lmax, i, nmax, nlat, nlong, gridtype, astat(12), n_out, &
			iter, j, r1, lmaxp, lmaxt, filter_type, half, degmax, sampling
	character*120 ::	grav_file, moho_out, thick_grid_out, topo_file, &
	        misfit_file

	print*,  "rho_crust (kg/m3) > "
	read(*,*) rho_crust
	print*, "rho_mantle (kg/m3) > "
	read(*,*) rho_mantle
	
	pi = acos(-1.0d0)
	delta_max = 5.0d0
	grav = Grav_constant

	gridtype = 3
	if (gridtype == 2) then
		sampling = 1
	elseif (gridtype == 3) then
		sampling = 2
	else
		stop
	endif

	nmax = 6	! nmax of Wieczorek and Phillips (1998)
	
	print*, "Input filter type (1) Minimum amplitude, (2) minimum curvature, (0) no filter "
	read(*,*) filter_type
	if (filter_type /= 0 ) then
		print*, "Degree at which the filter is 1/2 " 
		read(*,*) half
	endif
		
	grav_file = "../../ExampleDataFiles/jgmro_110b_sha.tab"
	topo_file = "../../ExampleDataFiles/MarsTopo719.shape"
	
	print*, "Remove degree 1 topo coefficients from Bouguer Correction? (0:no, 1:yes) > "
	read(*,*) r1
	
	print*, "maximum degree to compute Moho relief to >"
	read(*,*) degmax
	
	lmax = 2 * degmax ! this partially takes into account aliasing problems. Technically, it should be nmax*degmax
	
	nlat = 2*lmax+2
	nlong = 2*nlat

	print*, "Minimum assumed crustal thickness (km) > "
	read(*,*) t0
	t0 = t0 * 1.0d3
		
	print*, "Moho spherical harmonic coeficient output filename > "
	read(*,*) moho_out
	print*, "Grid spacing for output crustal thickness map (degrees) > "
	read(*,*) interval
	print*, "gridded crustal thickness output filename >"
	read(*,*) thick_grid_out
	print*, "Gravity misfit spherical harmonic filename >"
	read(*,*) misfit_file
	
	call cpu_time(timein)
	
	allocate(topo_grid(nlat, nlong), stat = astat(1))
	allocate(moho_grid(nlat, nlong), stat = astat(2))
	allocate(moho_grid2(nlat, nlong), stat = astat(3))
	allocate(moho_grid3(nlat, nlong), stat = astat(4))
	allocate(temp_grid(nlat, nlong), stat = astat(5))
	allocate(topo_c(2,degmax+1,degmax+1), stat = astat(6))
	allocate(moho_c(2,degmax+1,degmax+1), stat = astat(7))
	allocate(bc(2,degmax+1,degmax+1), stat = astat(8))
	allocate(ba(2,degmax+1,degmax+1), stat = astat(9))
	allocate(cilm(2,degmax+1,degmax+1), stat = astat(10))
	allocate(pot(2,degmax+1,degmax+1), stat = astat(11))
	allocate(misfit(2,degmax+1,degmax+1), stat = astat(12))
	
	if (sum(astat(1:12)) /= 0) then
		print*, "Problem allocating arrays."
		stop
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Read topo and grav files 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print*, "Reading data from ", grav_file
	call SHRead(grav_file, pot, lmaxp, header=param(1:2))
	
	gm = param(2) * 1.d9
	r_grav = param(1) * 1.d3
	mass = gm/grav
	
	print*, "Mass (kg) = ", mass
	print*, "Average surface gravity (m/s2) = ", gm/r_grav**2
	print*, "Lmax of gravitational potential file = ", lmaxp

	print*, "Reading data from ", topo_file
	call SHRead(topo_file, topo_c, lmaxt)
	print*, "Lmax of topography file = ", lmaxt
		
	r0 = topo_c(1,1,1)
	print*, "r0 (km) = ", r0/1.d3
		
	! Downward continue gravity coefficients from 3397 to MPR
	
	do l=2, lmaxp
		pot(1:2,l+1,1:l+1) = pot(1:2,l+1,1:l+1) * (r_grav/r0)**l
	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Create Bouger anomaly up to degree nmax
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	print*, "Creating Bouger anomaly"
	
	call MakeGridDH(topo_grid, n_out, topo_c, lmax, norm = 1, sampling = sampling, csphase = 1, lmax_calc = degmax)
	
	print*, "Maximum Topo (km) = ", maxval(topo_grid(1:nlat, 1:nlong))/1.d3
	print*, "Minimum Topo (km) = ", minval(topo_grid(1:nlat, 1:nlong))/1.d3
	
	call CilmPlus(bc, topo_grid, degmax, nmax, mass, rref, rho_crust, gridtype, n = nlat)
		
	ba = pot - bc	! This is the bouguer anomaly
	
	if (r1 == 1) ba(1:2,2,1:2) = 0.0d0  
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Compute crustal thickness by iterating at each
	!	reference moho depth. Gravity anomalies are calculated
	!	using n*lmax coefficients, but only the first lmax are
	!	used in the crustal thickness calculations.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	d = r0 - 44.0d3		! initial reference moho depth

	thick_delta = 1.d9
	
	do while (abs(thick_delta) > delta_max)
	
		write(*,*)
		print*, "Reference depth (km) = ", (r0-d)/1.d3
		
		moho_c(1,1,1) = d
		
		! first iteration
		do l=1, degmax
			if (filter_type == 0) then
				moho_c(1:2,l+1,1:l+1) = ba(1:2,l+1,1:l+1) * mass * dble(2*l+1) * ((r0/d)**l) &
				/(4.0d0*pi*(rho_mantle-rho_crust)*d**2)
			elseif (filter_type == 1) then
				moho_c(1:2,l+1,1:l+1) = downcontfilterma(l, half, r0, d) * ba(1:2,l+1,1:l+1) * mass * dble(2*l+1) * ((r0/d)**l) &
				/(4.0d0*pi*(rho_mantle-rho_crust)*d**2)
			elseif (filter_type == 2) then
				moho_c(1:2,l+1,1:l+1) = downcontfiltermc(l, half, r0, d) * ba(1:2,l+1,1:l+1) * mass * dble(2*l+1) * ((r0/d)**l) &
				/(4.0d0*pi*(rho_mantle-rho_crust)*d**2)
			endif
		enddo
        
		call MakeGridDH(moho_grid3, n_out, moho_c, lmax, norm = 1, sampling = sampling, csphase = 1, lmax_calc = degmax)
	
		max_thick = maxval(topo_grid(1:nlat, 1:nlong)-moho_grid3(1:nlat, 1:nlong))
		min_thick = minval(topo_grid(1:nlat, 1:nlong)-moho_grid3(1:nlat, 1:nlong))
		print*, "Maximum Crustal thickness (km) = ", max_thick/1.d3
		print*, "Minimum Crustal thickness (km) = ", min_thick/1.d3

		! second iteration
		if (filter_type == 0) then
			call BAtoHilm(moho_c, ba, moho_grid3, lmax, nmax, mass, r0, rho_mantle-rho_crust, gridtype, &
				lmax_calc = degmax)
		else
			call BAtoHilm(moho_c, ba, moho_grid3, lmax, nmax, mass, r0, rho_mantle-rho_crust, gridtype, &
				filter_type=filter_type, filter_deg=half, lmax_calc = degmax)
		endif
		call MakeGridDH(moho_grid2, n_out, moho_c, lmax, norm = 1, sampling = sampling, csphase = 1, lmax_calc = degmax)
	
		delta = maxval(abs(moho_grid3(1:nlat, 1:nlong)-moho_grid2(1:nlat, 1:nlong)))
		print*, "Delta (km) = ", delta/1.d3
		
		temp_grid(1:nlat, 1:nlong) = topo_grid(1:nlat, 1:nlong)-moho_grid2(1:nlat, 1:nlong)
		max_thick = maxval(temp_grid(1:nlat, 1:nlong))
		min_thick = minval(temp_grid(1:nlat, 1:nlong))
		print*, "Maximum Crustal thickness (km) = ", max_thick/1.d3
		print*, "Minimum Crustal thickness (km) = ", min_thick/1.d3
	
		iter = 0
		delta = 1.0d9
	
		do while(delta > delta_max)
		
			iter = iter + 1	
			print*, "Iteration ", iter
			
			moho_grid(1:nlat,1:nlong) = (moho_grid2(1:nlat,1:nlong) + moho_grid3(1:nlat,1:nlong)) / 2.0d0
	
			delta = maxval(abs(moho_grid(1:nlat, 1:nlong)-moho_grid2(1:nlat, 1:nlong)))
			print*, "Delta (km) = ", delta/1.d3
		
			temp_grid(1:nlat, 1:nlong) = topo_grid(1:nlat, 1:nlong)-moho_grid(1:nlat, 1:nlong)
			max_thick = maxval(temp_grid(1:nlat, 1:nlong))
			min_thick = minval(temp_grid(1:nlat, 1:nlong))
			print*, "Maximum Crustal thickness (km) = ", max_thick/1.d3
			print*, "Minimum Crustal thickness (km) = ", min_thick/1.d3
		
			moho_grid3(1:nlat, 1:nlong) = moho_grid2(1:nlat, 1:nlong)
			moho_grid2(1:nlat, 1:nlong) = moho_grid(1:nlat, 1:nlong)
		
			iter = iter +1	
			print*, "Iteration ", iter

			if (filter_type == 0) then
				call BAtoHilm(moho_c, ba, moho_grid2, lmax, nmax, mass, r0, rho_mantle-rho_crust, gridtype, &
					lmax_calc = degmax)
			else
				call BAtoHilm(moho_c, ba, moho_grid2, lmax, nmax, mass, r0, rho_mantle-rho_crust, gridtype, &
					filter_type=filter_type, filter_deg=half, lmax_calc = degmax)
			endif
			call MakeGridDH(moho_grid, n_out, moho_c, lmax, norm = 1, sampling = sampling, csphase = 1, lmax_calc = degmax)
				
			delta = maxval(abs(moho_grid(1:nlat, 1:nlong)-moho_grid2(1:nlat, 1:nlong)))
			print*, "Delta (km) = ", delta/1.d3
		
			temp_grid(1:nlat, 1:nlong) = topo_grid(1:nlat, 1:nlong)-moho_grid(1:nlat, 1:nlong)
			max_thick = maxval(temp_grid(1:nlat, 1:nlong))
			min_thick = minval(temp_grid(1:nlat, 1:nlong))
			print*, "Maximum Crustal thickness (km) = ", max_thick/1.d3
			print*, "Minimum Crustal thickness (km) = ", min_thick/1.d3
		
			moho_grid3(1:nlat, 1:nlong) = moho_grid2(1:nlat, 1:nlong)
			moho_grid2(1:nlat, 1:nlong) = moho_grid(1:nlat, 1:nlong)
		
			if (max_thick > 100.d3) then
				print*, "Not converging"
				stop
			endif
				
		enddo
	
		d = d + (min_thick-t0)
		thick_delta  = min_thick  - t0
	
	enddo
	
	thinnest(1:2) = minloc(temp_grid(1:nlat, 1:nlong))
	if (gridtype == 2) then
		print*, "Location of thinnest crust (lat, long) = ", &
			90.0d0 - dble((thinnest(1) - 1))*180.0/dble(nlat), dble(thinnest(2)-1)*360.0/dble(nlat)
	elseif(gridtype == 3) then			
		print*, "Location of thinnest crust (lat, long) = ", &
			90.0d0 - dble((thinnest(1) - 1))*180.0/dble(nlat), dble(thinnest(2)-1)*180.0/dble(nlat)
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine misfit between observed and calculated gravity, and write
	!	data to external files. Note that here, only coefficients up to lmax
	!	are considered.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call CilmPlus(cilm, moho_grid, degmax, nmax, mass, rref, rho_mantle-rho_crust, gridtype, n=nlat)
	
	! upward continue moho coefficients
	do l=0, degmax
		cilm(1:2,l+1,1:l+1) = cilm(1:2,l+1,1:l+1) * (d/r0)**l
	enddo
	
	misfit = pot - (bc + cilm) 	! this is the misfit
	misfit(1,1,1) = 0.0d0		! ignore degree-0 misfit
	
	do l=0, degmax
		cilm(1:2,l+1,1:l+1) = misfit(1:2,l+1,1:l+1)*dble(l+1.0d0)*gm*(1.0d5)/r0**2	
	enddo
	
	call MakeGridDH(temp_grid, n_out, cilm, lmax, norm = 1, sampling = sampling, csphase = 1, lmax_calc = degmax)
	
	print*, "Maximum misfit (mgals) = ", maxval(temp_grid(1:nlat, 1:nlong))
	print*, "Minimum misfit (mgals) = ", minval(temp_grid(1:nlat, 1:nlong))

	print*, "Mean Crustal Thickness (km) =", (r0-moho_c(1,1,1))/1.d3
	
	print*, "Writing output data"
	
	open(12,file=moho_out)
	
	do l=0,degmax
		do m=0,l
			write(12,*) l, m, moho_c(1,l+1,m+1), moho_c(2,l+1,m+1)
		enddo
	enddo
	
	close(12)
	
	call MakeGrid2d(topo_grid, topo_c - moho_c, degmax, interval, nlat, nlong)
	
	open(12, file=thick_grid_out)
	write(12,*) nlat, nlong
	do i=1, nlat
		do j=1,nlong
			write(12,*) topo_grid(i,j)/1.d3
		enddo
	enddo
	
	close(12)
	
	open(12, file=misfit_file)
	do l=0,degmax
		do m=0,l
			write(12,*) l, m, misfit(1,l+1,m+1), misfit(2,l+1,m+1)
		enddo
	enddo
	
	close(12)
	
	deallocate(topo_grid)
	deallocate(moho_grid)
	deallocate(moho_grid2)
	deallocate(moho_grid3)
	deallocate(temp_grid)
	deallocate(topo_c)
	deallocate(moho_c)
	deallocate(bc)
	deallocate(ba)
	deallocate(cilm)
	deallocate(pot)
	deallocate(misfit)
	
	call cpu_time(timeout)
	print*, "time (sec) = ", timeout-timein

end program MarsCrustalThickness
