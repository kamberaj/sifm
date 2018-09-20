!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: defs_mpi.f90,v 1.0 19-03-2018, IBU
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
module TEMPI_CLASS
use sifm_kinds
implicit none
!
include 'mpif.h'
!
! Kind constants for system-independent precision specification.
 ! H. Kamberaj

type val_t
    real(sifm_real) :: value
    integer :: rank
end type val_t

      integer, save :: NumProcs = 1
      integer, save :: MyID = 0
      integer, parameter :: master = 0
      integer , parameter :: msgtag2 = 12
      integer, save, dimension ( MPI_STATUS_SIZE ) :: status

contains

subroutine MPI_start()
implicit none

integer :: ierr
!
call MPI_INIT ( ierr )
call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE ( MPI_COMM_WORLD, NumProcs, ierr )
!
Return
end subroutine MPI_start

subroutine MPI_finish()
implicit none

integer :: ierr
!
call MPI_FINALIZE (ierr)
!
Return
end subroutine MPI_finish

end module TEMPI_CLASS

