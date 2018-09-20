!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: demo_setup.f90,v 1.0 19-03-2018, IBU
!
!                This source code is part of
!
!   Symbolic Information Flow Measure Code for styding the
!   information flow in dynamical systems
!
!                        VERSION 1.0
!
! Written by Hiqmet Kamberaj.
! Copyright (C) 2018 Hiqmet Kamberaj.
! Check out h.kamberaj@gmail.com for more information.
!
! This program is free software; you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the Free Software Foundation; 
! GPL-3.0
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this program; 
! if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
! Boston, MA 02111-1307 USA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program demo_setup
!
      USE sifm_kinds
      use SYSTEM_CLASS, only : SETUP_DRIVER
      use TEUTILS_CLASS, only : addStrings, new_unit
      implicit none
!
      Integer :: Natoms=2
      Integer :: offset = 1
      integer :: Nframes=5000
      integer :: Ndim=1
      integer :: iStart=0
      integer :: iSkip =1
      integer :: iStop =5000
      integer :: debug=1
      integer :: benchmark = 1
      character(len=10) :: buffer
      real(sifm_real) :: Ax, Cxy, Sigma_x, Sigma_y
!		  
      IF (iargc() < 8) THEN
         CALL getarg(1, buffer)
          IF (buffer(1:4) == "help") THEN
              call help()
              return
          ELSE
              write(*,'("Wrong command line")')		
              call help()
              STOP "Restart the run!"
          ENDIF
      ENDIF
      CALL getarg(1, buffer)
      read(buffer,*, ERR=999) iStart
      write(*,*) "Time to start trajectory - iStart:          ", iStart

      CALL getarg(2, buffer)
      read(buffer,*, ERR=999) iSkip
      write(*,*) "Number of steps to skip - iSkip:            ", iSkip

      CALL getarg(3, buffer)
      read(buffer,*, ERR=999) iStop
      write(*,*) "Time to stop trajectory - iStop:            ", iStop

      CALL getarg(4, buffer)
      read(buffer,*, ERR=999) Offset
      write(*,*) "Offset residue numbering - Offset:          ", Offset

      CALL getarg(5, buffer)
      read(buffer,*, ERR=999) Natoms
      write(*,*) "Number of time series - Natoms:             ", Natoms
	  
      CALL getarg(6, buffer)
      read(buffer,*, ERR=999) Ndim
      write(*,*) "Dimensionality of the problem - Ndim:       ", Ndim

      CALL getarg(7, buffer)
      read(buffer,*, ERR=999) Debug
      write(*,*) "Set debuging flag value - Debug:            ", Debug

      CALL getarg(8, buffer)
      read(buffer,*, ERR=999) benchmark
      write(*,*) "Set benchmark flag value - Benchmark:       ", benchmark
	  
      IF (benchmark == 1) THEN
          CALL getarg(9, buffer)
          read(buffer,*, ERR=999) Ax
          write(*,*) "Set Ax value - Ax:                      ", Ax    
          CALL getarg(10, buffer) 
          read(buffer,*, ERR=999) Cxy
          write(*,*) "Set Cxy value - Cxy:                    ", Cxy   
          CALL getarg(11, buffer)
          read(buffer,*, ERR=999) Sigma_x
          write(*,*) "Set SigmaX value - SigmaX:              ", Sigma_x		  
          CALL getarg(12, buffer)
          read(buffer,*, ERR=999) Sigma_y
          write(*,*) "Set SigmaY value - SigmaY:              ", Sigma_y    		  
      ENDIF
      goto 1000 

999   write(*,*) "Bad input???"
1000  continue   
!
      nframes = (iStop - iStart + 1) / iSkip
!	  

      IF (benchmark == 1) then
          CALL SETUP_DRIVER(Nframes, Ndim, Natoms, Offset, iStart, iSkip, iStop,     &
	                    debug, benchmark, Ax, Cxy, Sigma_x, Sigma_y)
      else
          CALL SETUP_DRIVER(Nframes, Ndim, Natoms, Offset, iStart, iSkip, iStop,     &
	                debug, benchmark)

      endif
!
      STOP 'Program Finished Successfully'
! 
Contains
      subroutine help()
!     	  
            write(*,'("Enter:")')
            write(*,*) "1  - Time to start trajectory     - iStart"
            write(*,*) "2  - Number of steps to skip      - iSkip"
            write(*,*) "3  - Time to stop trajectory      - iStop"
            write(*,*) "4  - Offset residue numbering     - Offset"
            write(*,*) "5  - Number of time series        - Natoms"
            write(*,*) "6  - Dimenionality of the problem - Ndim"
            write(*,*) "7  - Set debuging flag value      - Debug"				
            write(*,*) "8  - Set Bechmark flag value      - Benchmark"			
!					
          return
     end subroutine 
end program
