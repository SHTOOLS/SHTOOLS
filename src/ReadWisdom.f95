!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!	Notes from the FFTW3.0.1 f77_wisdom.f file:
!
!     	This is an example implementation of Fortran wisdom export/import
!     	to/from a Fortran unit (file), exploiting the generic
!     	dfftw_export_wisdom/dfftw_import_wisdom functions.
!     
!     	We cannot compile this file into the FFTW library itself, lest all
!     	FFTW-calling programs be required to link to the Fortran I/O
!     	libraries.
!
!	History
!		Based on file f77_wisdom.f and turned into f95 code by Mark Wieczorek
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_char(c, iunit)
      	character :: 	c
      	integer	::	iunit
      	
      	write(iunit, fmt="(a)", advance="no") c
      	
end subroutine write_char      


subroutine export_wisdom_to_file(iunit)
      	integer :: 	iunit
      	external 	write_char
      	
     	call dfftw_export_wisdom(write_char, iunit)

end subroutine export_wisdom_to_file


subroutine read_char(c, iunit)
      	integer :: 	c
      	integer :: 	iunit
      
        read(iunit, fmt="(a)", advance="no") c
      
end subroutine read_char
      

subroutine import_wisdom_from_file(isuccess, iunit)
      integer :: 	isuccess
      integer ::	iunit
      external		read_char
      
      call dfftw_import_wisdom(isuccess, read_char, iunit)

end subroutine import_wisdom_from_file

