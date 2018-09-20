!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: demo_te.f90,v 1.0 19-03-2018, IBU
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
program demo_te
!
      USE sifm_kinds
      use TE_CLASS, only : TE_DRIVER, STE1_DRIVER, STE2_DRIVER, STE3_DRIVER
      implicit none
!
      integer :: Nframes = 5000
      integer :: Ndim = 1
      integer :: Natoms = 2
      integer :: debug = 1
      integer :: qteMethod = 1
      integer :: qteShuffle = 0
      integer :: qteNorm = 0
      integer :: Nshuffles = 0
      real(sifm_real) :: Rcut=0.02_sifm_real, StatP=0.95_sifm_real
      character(len=10) :: buffer
!
      IF (iargc() < 10) THEN
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
      write(*,*) "Length of time series - Nframes:                 ", Nframes

      CALL getarg(2, buffer)
      read(buffer,*, ERR=999) Natoms
      write(*,*) "Number of time series - Natoms:                  ", Natoms
	  
      CALL getarg(3, buffer)
      read(buffer,*, ERR=999) Ndim
      write(*,*) "Dimensionality of the problem - Ndim:            ", Ndim

      CALL getarg(4, buffer)
      read(buffer,*, ERR=999) qTeMethod 
      write(*,*) "Flag for Method of TE calculation - qteMethod:   ", qTeMethod

      CALL getarg(5, buffer)
      read(buffer,*, ERR=999) qTeNorm
      write(*,*) "Flag for Normalization of TE - qteNorm:          ", qTeNorm
	  
      CALL getarg(6, buffer)
      read(buffer,*, ERR=999) qTeShuffle
      write(*,*) "Flag for Method of Shuffling of TE - qteShuffle: ", qTeShuffle
	  
      CALL getarg(7, buffer)
      read(buffer,*, ERR=999) Nshuffles
      write(*,*) "Number of Shuffling of time series - Nshuffles:  ", Nshuffles
	  
      CALL getarg(8, buffer)
      read(buffer,*, ERR=999) Rcut
      write(*,*) "Cutoff for TE minimum value - Rcut:              ", Rcut

      CALL getarg(9, buffer)
      read(buffer,*, ERR=999) statP
      write(*,*) "Confidence level for averages - statP:           ", statP
	  
      CALL getarg(10, buffer)
      read(buffer,*, ERR=999) Debug
      write(*,*) "Set debuging flag value - Debug:                 ", Debug
	 
      goto 1000 

999   write(*,*) "Bad input???"
1000  continue   
! 
      IF (qTEMethod == 1) THEN
          CALL STE1_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)
      ELSEIF (qTEMethod == 2) THEN
          CALL STE2_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)
      ELSEIF (qTEMethod == 3) THEN
          CALL STE3_DRIVER(Nframes, Ndim, Natoms,qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)
      ELSE
          CALL TE_DRIVER(Nframes, Ndim, Natoms, qteShuffle, qteNorm, debug, Nshuffles, Rcut, StatP)      
      ENDIF
!
      STOP 'Program Finished Successfully'
! 
Contains
      subroutine help()
!
          write(*,'("Enter:")')
          write(*,*) "1  - Length of time series                                                                   - Nframes"
          write(*,*) "2  - Number of time series                                                                   - Natoms"
          write(*,*) "3  - Dimenionality of the problem                                                            - Ndim"
          write(*,*) "4  - Flag for Method of TE calculation (1: Discrete (m=1); 2: Discrete (m>1); 3: Schreiber)  - qTEMethod"
          write(*,*) "5  - Flag for Normalization of TE (0: No normalization; 1: normalization                     - qTENorm"
          write(*,*) "6  - Flag for Method of Shuffling of TE (1: Permutation; 2: Block shuffling)                 - qTEShuffle"
          write(*,*) "7  - Number of Shuffling of time series                                                      - Nshuffles"
          write(*,*) "8  - Cutoff for Mutual Information minimum value                                             - Rcut"
          write(*,*) "9  - Confidence level for averages                                                           - statP"
          write(*,*) "10 - Set debuging flag value                                                                 - Debug"
!					
          return
     end subroutine     	  
end program demo_te
