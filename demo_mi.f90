! Id: demo_mi.f90,v 1.0 19-03-2018, IBU
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
program demo_mi
!
      USE sifm_kinds
      use MI_CLASS, only : MI_DRIVER, SMI_DRIVER
      implicit none
!
      integer :: Nframes = 5000
      integer :: Ndim = 1
      integer :: Natoms = 2
      integer :: debug = 1
      integer :: qMIMethod = 1
      integer :: qMIShuffle = 0
      integer :: Nshuffles = 0
      real(sifm_real) :: Rcut=0.02_sifm_real, StatP=0.95_sifm_real
      character(len=10) :: buffer
!
      IF (iargc() < 9) THEN
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
      write(*,*) "Length of time series - Nframes:                    ", Nframes

      CALL getarg(2, buffer)
      read(buffer,*, ERR=999) Natoms
      write(*,*) "Number of time series - Natoms:                     ", Natoms
	  
      CALL getarg(3, buffer)
      read(buffer,*, ERR=999) Ndim
      write(*,*) "Dimensionality of the problem - Ndim:               ", Ndim

      CALL getarg(4, buffer)
      read(buffer,*, ERR=999) qMIMethod
      write(*,*) "Flag for Method of MI calculation - qMIMethod:      ", qMIMethod
	  
      CALL getarg(5, buffer)
      read(buffer,*, ERR=999) qMIShuffle
      write(*,*) "Flag for Method of Shuffling of MI - qMIShuffle:    ", qMIShuffle
	  
      CALL getarg(6, buffer)
      read(buffer,*, ERR=999) Nshuffles
      write(*,*) "Number of Shuffling of time series - Nshuffles:     ", Nshuffles
	  
      CALL getarg(7, buffer)
      read(buffer,*, ERR=999) Rcut
      write(*,*) "Cutoff for MI minimum value - Rcut:                 ", Rcut

      CALL getarg(8, buffer)
      read(buffer,*, ERR=999) statP
      write(*,*) "Confidence level for averages - statP:              ", statP
	  
      CALL getarg(9, buffer)
      read(buffer,*, ERR=999) Debug
      write(*,*) "Set debuging flag value - Debug:                    ", Debug
	 
      goto 1000 

999   write(*,*) "Bad input???", buffer

1000  continue   
! 
      IF (qMIMethod == 1) THEN
          CALL SMI_DRIVER(Nframes, Ndim, Natoms, qMIShuffle, debug, Nshuffles, Rcut, StatP)
      ELSE
          CALL  MI_DRIVER(Nframes, Ndim, Natoms, qMIShuffle, debug, Nshuffles, Rcut, StatP)      
      ENDIF
!
      STOP 'Program Finished Successfully'
!  
contains
     subroutine help()

          write(*,'("Enter:")')
          write(*,*) "1  - Length of time series                                                  - Nframes"
          write(*,*) "2  - Number of time series                                                  - Natoms"
          write(*,*) "3  - Dimenionality of the problem                                           - Ndim"
          write(*,*) "4  - Flag for Method of MI calculation (1: Discrete; 2: Schreiber)          - qMIMethod"
          write(*,*) "5  - Flag for Method of Shuffling of MI (1: Permutation; 2: Block shuffling)- qMIShuffle"
          write(*,*) "6  - Number of Shuffling of time series                                     - Nshuffles"
          write(*,*) "7  - Cutoff for Mutual Information minimum value                            - Rcut"
          write(*,*) "8  - Confidence level for averages                                          - statP"
          write(*,*) "9  - Set debuging flag value                                                - Debug"
!
          return
  
!
  end subroutine help 
!   	  
end program
