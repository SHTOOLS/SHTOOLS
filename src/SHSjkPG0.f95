real*8 function SHSjkPG0(incspectra, j, k, l, m, evec, lwin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the expected windowed cross-power spectra for a two fields,
!	which are windowed by tapers j and k. Note that this is only valid for zonal tapers. 
!	This corresponds to the variable (this is eq. D8 of Wieczoek and Simons 2005):
!
!	   (j,k)
!	< S      (l,m) >
!	   PG
!
!	It can be shown that SPGjk(l,-m) = SPGjk(l,m).
!
!	Calling Parameteters
!		IN
!			incspectra	Knonw input (cross) power spectra
!					as a function of degree.
!			j, k		Taper numbers
!			l, m 		Angular degree and order.
!			evec		Matrix containing the spectral coefficients of the 
!					zonal tapers.
!			lwin		Spherical harmonic bandwidth of the windows.
!
!	Dependencies: Wigner3j.
!
!	Written by Mark Wieczorek, May 2006.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	use SHTOOLS, only: Wigner3j
	implicit none
	
	real*8, intent(in) ::	incspectra(:), evec(:,:)
	integer, intent(in) ::	lwin, l, m, j, k
	integer ::	i, l1, l2, imin, imax, &
			l10min, l10max, l1min, l1max, l20min, l20max, l2min, l2max
	real*8 :: 	wl10(lwin+l+1), wl20(lwin+l+1), &
			wl1(lwin+l+1), wl2(lwin+l+1)
	
		
	if (size(evec(1,:)) < lwin+1 .or. size(evec(:,1)) < lwin+1) then
		print*, "Error --- SHSjkPG"
		print*, "EVEC must be dimensioned (LWIN+1, LWIN+1) where LWIN is ", lwin
		print*, "Input array is dimensioned", size(evec(:,1)), size(evec(1,:))
		stop
	elseif(size(incspectra(:)) < l+lwin + 1) then
		print*, "Error --- SHSjkPG"
		print*, "INCSPECTRA must be dimensioned as (L+LWIN+1) where L and LWIN are ", l, lwin
		print*, "Input arrays is dimensioned ", size(incspectra)
		stop
	elseif(j>lwin+1 .or. k>lwin+1) then
		print*, "Error --- SHSjkPG"
		print*, "J and K must be less then LWIN+1."
		print*, "J = ", j
		print*, "K = ", k
		stop
	endif
	
	SHSjkPG0 = 0.0d0
	
	if (l<abs(m)) return
	
	do l1=0, lwin, 1
		call Wigner3j(wl10, l10min, l10max, l1, l, 0, 0, 0)
		call Wigner3j(wl1, l1min, l1max, l1, l, m, 0, -m)
		
		do l2 = 0, lwin, 1
			call Wigner3j(wl20, l20min, l20max, l2, l, 0, 0, 0)
			call Wigner3j(wl2, l2min, l2max, l2, l, m, 0, -m)
			imin = max(l1min, l2min)
			imax = min(l1max, l2max)
		
			do i=imin, imax, 1	
				if (mod(i+l1+l,2) == 0 .and. mod(i+l2+l,2) == 0 ) then
					SHSjkPG0 = SHSjkPG0 + evec(l1+1,j)*evec(l2+1,k) * &
						incspectra(i+1) * sqrt(2.0d0*l1+1.0d0) * sqrt(2.0d0*l2+1.0d0) * &
						(2.0d0*l+1.0d0) *  wl10(i-l10min+1) * wl20(i-l20min+1) * &
						wl1(i-l1min+1) * wl2(i-l2min+1)
				endif	
			enddo
		enddo
	enddo	
	
end function SHSjkPG0


