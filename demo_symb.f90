!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: demo_symb.f90,v 1.0 19-03-2018, IBU
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
program demo_symb
!
      USE sifm_kinds
      use SYMB_CLASS, only : SYMB_DRIVER
      implicit none
!
      Integer :: Natoms=2
      integer :: Nframes=5000
      integer :: Ndim=1
      integer :: debug=1
      integer :: qSymbolic = 3
      integer :: Nmc = 100
      character(len=10) :: buffer 
!
      IF (iargc() < 6) THEN
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
      read(buffer,*, ERR=999) Debug
      write(*,*) "Set debuging flag value - Debug:            ", Debug

      CALL getarg(5, buffer)
      read(buffer,*, ERR=999) qSymbolic
      write(*,*) "Set Symbolic method flag - qSymbolic:       ", qSymbolic

      CALL getarg(6, buffer)
      read(buffer,*, ERR=999) Nmc
      write(*,*) "Get Nr of Monte Carlo Steps:                ", Nmc
	  
      goto 1000 

999   write(*,*) "Bad input???"
1000  continue   
!	  
      CALL SYMB_DRIVER(Nframes, Ndim, Natoms, debug, qSymbolic, Nmc)
!
      STOP 'Program Finished Successfully'
!    
contains
     subroutine help()
!
          write(*,'("Enter:")')
          write(*,*) "1  - Number of Time Frames                                                - Nframes"
          write(*,*) "2  - Number of time series                                                - Natoms"
          write(*,*) "3  - Dimenionality of the problem                                         - Ndim"
          write(*,*) "4  - Set debuging flag value                                              - Debug"			
          write(*,*) "5-   Symbolic method flag (1: 0s, 1s; 2: 0s, 1s, 2s, etc; 3: Monte Carlo) - qSymbolic"			
          write(*,*) "6  - Set Nr of Monte Carlo Steps                                          - Nmc"			
!
          return
    end subroutine help  	  
end program
