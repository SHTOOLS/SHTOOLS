subroutine SHRead(filename, cilm, lmax, skip, header, error)
!-------------------------------------------------------------------------------
!
!   This subroutine will open an ascii file, and read in all of the 
!   spherical harmonic coefficients. If the option "header"
!   is specified, then the first Length(header) records of the 
!   spherical harmonic file will be retuned in the array header.
!   
!   The file will be read until the end of the file is encountered
!   or until the maximum length of cilm (as dimensioned in the calling 
!   program) is reached.
!
!   Calling Parameters
!
!       IN
!           filename        Character name of the ascii file.
!
!       OUT
!           cilm        Spherical harmonic coeficients with dimensions
!                       (2, lmax+1, lmax+1).
!           lmax        Maximum spherical harmonic degree of cilm.
!
!       OPTIONAL
!           header      Array of the first length(header) data records
!                       in the file. Called as header = header, or 
!                       header = header(1:number_of_header_records).
!           skip        Number of lines to skip
!           error       Error of the spherical harmonic coefficients, assumed 
!                       to in the format (l, m, c1lm, c2lm, error1lm, error2lm).
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    character(*), intent(in) :: filename
    integer, intent(out) :: lmax
    real*8, intent(out) :: cilm(:,:,:)
    real*8, intent(out), optional :: header(:), error(:,:,:)
    integer, intent(in), optional :: skip
    integer :: l, m, stat, ll, mm, lmax2, lstart, headlen, fu
    
    lmax = 0
    cilm = 0.0d0
    fu = 101
    
    if (size(cilm(:,1,1)) < 2 ) then
        print*, "Error --- SHRead"
        print*, "CILM must be dimensioned (2, *, *)."
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        stop
        
    end if
    
    lmax2 = min(size(cilm(1,1,:) ) - 1, size(cilm(1,:,1) ) - 1)
    
    open(fu, file=filename, status="old")
        
    if (present(header)) headlen = size(header)
    
    !---------------------------------------------------------------------------
    !
    !   Skip lines and read header information
    !
    !---------------------------------------------------------------------------
    if (present(skip) ) then
    
        do l = 1, skip, 1
            read(fu,*, iostat=stat)
            
            if (stat /= 0 ) then
                print*, "Error --- SHRead"
                print*, "Problem skipping first lines of ", filename
                print*, "Line number = ", l
                print*, "Number of lines to skip = ", skip
                stop
                
            end if
            
        end do
        
    end if
    
    if (present(header) ) then
        read(fu,*, iostat=stat) (header(l), l=1, headlen,1)
        
        if (stat /= 0 ) then
            print*, "Error --- SHRead"
            print*, "Problem reading header line ", filename 
            stop
            
        end if
        
    end if
    
    !---------------------------------------------------------------------------
    !
    !   Determine first l value
    !
    !---------------------------------------------------------------------------
    read(fu,*, iostat=stat) ll
    
    if (stat /= 0 ) then
        print*, "Error --- SHRead"
        print*, "Problem reading first line of ", filename
        stop
        
    end if

    lstart = ll
    
    rewind (fu)
    
    if ( present(skip) ) then
        do l = 1, skip, 1
            read(fu,*, iostat=stat)
            
        end do
        
    endif
    
    if (present(header)) then
        read(fu,*, iostat=stat) (header(l), l=1, headlen)
        
    end if
    
    !---------------------------------------------------------------------------
    !
    !   Read coefficients
    !
    !---------------------------------------------------------------------------
    do l = lstart, lmax2, 1
        do m = 0,l
            if ( present(error) ) then
                read(fu,*,iostat=stat) ll, mm, cilm(1,l+1,m+1), &
                            cilm(2,l+1,m+1), error(1,l+1,m+1), error(2,l+1,m+1)
                            
            else
                if (m == 0) then
                    read(fu,*,iostat=stat) ll, mm, cilm(1,l+1,m+1)
                    
                else
                    read(fu,*,iostat=stat) ll, mm, cilm(1,l+1,m+1), &
                            cilm(2,l+1,m+1)
                            
                end if
                
            end if
                if (stat < 0) then
                exit
                
            else if (stat > 0) then
                print*, "Error --- SHRead"
                print*, "Problem reading file ", filename
                stop
                
            else if (ll /=l .or. mm /= m) then
                print*, "Error --- SHRead"
                print*, "Problem reading file ", filename
                print*, "Expected indices (l,m) = ", l, m
                print*, "Read indices (l,m) = ", ll, mm
                stop
                
            end if
            
        end do
        
        if (stat < 0) exit
    
        lmax = l    
        
    end do
    
    close (fu)
        
end subroutine SHRead
