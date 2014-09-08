subroutine MakeCircleCoord(coord, lat, lon, theta0, cinterval, cnum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will return the latitude and longitude
!	coordinates of a circle of radius theta from a given 
!	point. The first index in the output vectors corresponds
!	to the point directly north of the central cirlce coordinates,
!	and subsequent points arranged in a clockwise manner.
!
!	Calling Parameters:
!		IN
!			lat, lon:	Latitude and longitude in degrees.
!			theta0:		Angular radius of the circle in degrees.
!		OUT
!			coord:		360/interval (latitude, longitude) coordinates.
!		OPTIONAL, IN
!			cinterval:	Angular spacing of latitude and longitude points
!					in degrees (default=1).
!		OPTIONAL, OUT
!			cnum:		Number of points in the output vectors.
!		
!
!	Dependencies: None
!
!	Written by Mark Wieczorek December 2005.
!
!	October 29, 2012. Modified the order of the output vectors so that the first
!	point is directly north, and the following points are arranged clockwise.
!
!	Copyright (c) 2005-2012, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) ::	lat, lon, theta0
	real*8, intent(out) :: coord(:,:)
	real*8, intent(in), optional ::	cinterval
	integer, intent(out), optional :: cnum
	
	real*8 ::  pi, interval, xold, yold, zold, x, y, z, x1, phi
	integer ::   k, num
	
	if (theta0 == 0.0d0) then
		coord(1,1) = lat
		coord(1,2) = lon
		if (present(cnum)) then
			cnum = 1
		endif
		return
	endif
	
	if (present(cinterval)) then
		interval = cinterval
	else
		interval = 1.0d0
	endif
	
	num = 360.0d0/interval
	
	if (present(cnum)) then
		cnum = num
	endif
	
	if (size(coord(:,1)) < num .or. size(coord(1,:)) < 2) then
		print*, "Error --- MakeCircleCoord"
		print*, "COORD must be dimensioned as (NUM, 2) where NUM is ", NUM
		print*, "Input array is dimensioned as ", size(coord(:,1)), size(coord(1,:))
		stop
	endif
	
	pi = acos(-1.0d0)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate grid points. First create a cirlce, then rotate these
	!	points.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	zold = cos(theta0*pi/180.0d0)
		
	do k=1, num
		phi = pi - dble(k-1)*(2.0d0*pi/dble(num))
		xold = sin(theta0*pi/180.0d0)*cos(phi)
		yold = sin(theta0*pi/180.0d0)*sin(phi)
			
		! rotate coordinate system 90-lat degrees about y axis
			
		x1 = xold*cos(pi/2.0-lat*pi/180.0d0) + zold*sin(pi/2.0-lat*pi/180.0d0)
		z = -xold*sin(pi/2.0-lat*pi/180.0d0) + zold*cos(pi/2.0-lat*pi/180.0d0)

		! rotate coordinate system lon degrees about z axis
			
		x = x1*cos(lon*pi/180.0d0) - yold*sin(lon*pi/180.0d0)
		y = x1*sin(lon*pi/180.0d0) + yold*cos(lon*pi/180.0d0)
		
		coord(k,1) = (pi/2.0d0 - acos(z/sqrt(x**2+y**2+z**2)) ) * 180.0d0/pi
		coord(k,2) = atan2(y, x) * 180.0d0/pi

	enddo

end subroutine MakeCircleCoord
