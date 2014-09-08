Program TimingAccuracyDHC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will test the accuracy of the spherical harmonic complex 
!	Driscoll and Healy tranformation routines by expanding a field in the 
!	space domain, transforming this to spherical harmonics, and comparing 
!	the relative error of the coefficients.
!
!
!	Written by Mark Wieczorek (July 2006)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS
	
	implicit none	
		
	integer, parameter ::	maxdeg = 2800
	character*200 ::		outfile1, outfile2, outfile3, outfile4, outfile
	complex*16, allocatable ::	cilm(:,:,:), cilm2(:,:,:), griddh(:,:)
	real*8 ::		huge8, maxerror, err1, err2, beta, rms, timein(3), timeout(3)
	integer ::		lmax, l, m, seed, n, sampling, lmaxout, astat(3)
	
	allocate(cilm(2,maxdeg+1,maxdeg+1), stat=astat(1))
	allocate(cilm2(2,maxdeg+1,maxdeg+1), stat=astat(2))
	allocate(griddh(2*maxdeg+2,4*maxdeg+4), stat=astat(3))
	
	if (sum(astat(1:3)) /= 0) then
		print*, "Error --- TimingAccuracyDHC"
		print*, "Problem allocating arrays CILM, CILM2, GRIDDH"
		stop
	endif
	
	
	sampling = 1
	print*, "Value of beta for power law Sff = l^(-beta) > "
	read(*,*) beta
	
	print*, "output file name > "
	read(*,*) outfile
	
	outfile1 = trim(outfile) // ".timef"
	outfile2 = trim(outfile) // ".timei"
	outfile3 = trim(outfile) // ".maxerror"
	outfile4 = trim(outfile) // ".rmserror"
		
	huge8 = huge(maxerror)
	
	seed = -1053253
	
	cilm = (0.0d0, 0.0d0)
	
	do l=1, maxdeg
		
		do m=0, l
			if (m==0) then
				cilm(1,l+1,m+1) = dcmplx(RandomGaussian(seed), RandomGaussian(seed))
			else
				cilm(1,l+1,m+1) = dcmplx(RandomGaussian(seed), RandomGaussian(seed))
				cilm(2,l+1,m+1) = dcmplx(RandomGaussian(seed), RandomGaussian(seed))
			endif
		enddo
		cilm(1:2, l+1, 1:l+1) = cilm(1:2, l+1, 1:l+1)*sqrt(dble(l)**beta)/sqrt(2.0d0*l+1)
	enddo	
	
	print*, "Lmax, Maximum rel. error of Cilm, RMS relative error, Time inverse (sec), Time forward (sec)"

	lmax = 1
	
	do 
		
		lmax = lmax*2
		if (lmax > maxdeg) lmax = maxdeg
		
		n = 2*lmax + 2
		
		if (lmax==2) then
			open(12,file=outfile1, status="replace")
			open(13,file=outfile2, status="replace")
			open(14,file=outfile3, status="replace")
			open(15,file=outfile4, status="replace")
		else
			open(12,file=outfile1, position="append")
			open(13,file=outfile2, position="append")
			open(14,file=outfile3, position="append")
			open(15,file=outfile4, position="append")
		endif
		
		call cpu_time(timein(2))	
		if (sampling == 1) then
			call MakeGridDHC(griddh(1:n, 1:n), n, cilm(1:2,1:lmax+1, 1:lmax+1), lmax, sampling=sampling, norm=1)
		else
			call MakeGridDHC(griddh(1:n, 1:2*n), n, cilm(1:2,1:lmax+1, 1:lmax+1), lmax, sampling=sampling)
		endif
		call cpu_time(timeout(2))
		
		
		call cpu_time(timein(3))
		if (sampling == 1) then
			call SHExpandDHC(griddh(1:n, 1:n), n, cilm2(1:2, 1:lmax+1, 1:lmax+1), lmaxout, sampling=sampling, norm=1)
		else
			call SHExpandDHC(griddh(1:n, 1:2*n), n, cilm2(1:2, 1:lmax+1, 1:lmax+1), lmaxout, sampling=sampling)
		endif
		call cpu_time(timeout(3))
		
		maxerror = 0.0d0
		rms = 0.0d0
		
		do l=1, lmax
			
			do m=0, l
				if (m==0) then
					err1 = abs( (cilm(1,l+1,m+1)-cilm2(1,l+1,m+1)) ) / abs( cilm(1,l+1,m+1) )
					if (err1 >= maxerror) maxerror = err1
					rms = rms + err1**2
				else
					err1 = abs( (cilm(1,l+1,m+1)-cilm2(1,l+1,m+1)) ) / abs( cilm(1,l+1,m+1) )
					err2 = abs( (cilm(2,l+1,m+1)-cilm2(2,l+1,m+1)) ) / abs( cilm(2,l+1,m+1) )
					if (err1 >= maxerror) maxerror = err1
					if (err2 >= maxerror) maxerror = err2
					rms = rms + err1**2 + err2**2
				endif
			enddo
		enddo
		rms = sqrt(rms/dble(l+1)**2)
		
		! elasped time in seconds!
		print*, lmax, maxerror, rms, timeout(2)-timein(2), timeout(3)-timein(3)
		write(12,*) lmax, timeout(2)-timein(2)
		write(13,*) lmax, timeout(3)-timein(3)
		write(14,*) lmax, maxerror
		write(15,*) lmax, rms
		
		if (maxerror > huge8) then
			print*, "Overflow problems"
			print*, "Overflow problems"
			close(12)
			close(13)
			close(14)
			close(15)
			stop
		endif
		
		if (maxerror > 10.0d0) then
			close(12)
			close(13)
			close(14)
			close(15)
			stop
		endif
		
		close(12)
		close(13)
		close(14)
		close(15)
		
		if (lmax == maxdeg) exit
			
	enddo

end program TimingAccuracyDHC
