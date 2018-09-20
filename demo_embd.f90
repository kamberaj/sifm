!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: demo_embd.f90,v 1.0 19-03-2018, IBU
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
program demo_embd
!
     USE sifm_kinds
     use EMBD_CLASS, only : EMBD_DRIVER
     implicit none
!
      Integer :: Natoms=2
      integer :: Nframes=5000
      integer :: Ndim=1
      integer :: M1=2
      integer :: M2=20
      integer :: T1=1
      integer :: T2=1
      integer :: debug=1
      character(len=10) :: buffer
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
      read(buffer,*, ERR=999) Nframes
      write(*,*) "Number of Time Frames - Nframes:            ", Nframes

      CALL getarg(2, buffer)
      read(buffer,*, ERR=999) Natoms
      write(*,*) "Number of time series - Natoms:             ", Natoms
	  
      CALL getarg(3, buffer)
      read(buffer,*, ERR=999) Ndim
      write(*,*) "Dimensionality of the problem - Ndim:       ", Ndim

      CALL getarg(4, buffer)
      read(buffer,*, ERR=999) T1
      write(*,*) "Time lag minimum value - T1:                ", T1
	  
      CALL getarg(5, buffer)
      read(buffer,*, ERR=999) T2
      write(*,*) "Time lag maximum  value - T2:               ", T2

      CALL getarg(6, buffer)
      read(buffer,*, ERR=999) M1
      write(*,*) "Embedded dimension minimum value - M1:      ", M1

      CALL getarg(7, buffer)
      read(buffer,*, ERR=999) M2
      write(*,*) "Embedded dimension maximum value - M2:      ", M2

      CALL getarg(8, buffer)
      read(buffer,*, ERR=999) Debug
      write(*,*) "Set debuging flag value - Debug:            ", Debug
	  
      goto 1000 

999   write(*,*) "Bad input???"
1000  continue   
!
      CALL EMBD_DRIVER(Nframes, Ndim, Natoms, M1, M2, T1, T2, debug)
!
      STOP 'Program Finished Successfully'
!     
contains 
!
      subroutine help()
         write(*,'("Enter:")')
         write(*,*) "1  - Number of Time Frames            - Nframes"
         write(*,*) "2  - Number of time series            - Natoms"
         write(*,*) "3  - Dimenionality of the problem     - Ndim"
         write(*,*) "4  - Time lag minimum value           - T1"
         write(*,*) "5  - Time lag maximum value           - T2"
         write(*,*) "6  - Embedded dimension minimum value - M1"
         write(*,*) "7  - Embedded dimension maximum value - M2"
         write(*,*) "8  - Set debuging flag value          - Debug"			
         return
      end subroutine help
! 	  
end program
