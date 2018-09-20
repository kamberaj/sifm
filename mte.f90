!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: mte.f90,v 1.0 19-03-2018, IBU
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
module TE_class
use sifm_kinds
use TEutils_class
use TErandom_class
use TELinkList_class
implicit none
!
      integer, save :: nframes, ndim, natoms
      real(sifm_real), save, dimension(:,:), allocatable :: TES,TE
      character(len=1), save, allocatable :: xyzs(:,:,:)
      real(sifm_real), save, allocatable :: xyz(:,:,:)
      integer, save, allocatable :: mopt(:,:), topt(:,:)

!
CONTAINS
!
subroutine STE1_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, debug, &
                      Nshuffles, r0, statP)
  use tempi_class
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!  
      integer :: ierr
      real(sifm_real) :: LocalTime_Start, LocalTime_end
      real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  nframes = nf
  natoms = na
  ndim = nd
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_TransferEntropy(qTeShuffle, 1) 
!
  CALL read_embdparam()
  Mopt = 1
  CALL read_xyzs()
!
! -- Initialize MPI
  call MPI_start()
  call MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, 1)
!
  Localtime_start = MPI_Wtime()
!
  IF (qteShuffle > 0) THEN
      CALL getSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  ELSE
      CALL getSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
  ENDIF
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for TE:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  IF (myID == Master) call WRITE_TE(qTeShuffle,1)
! 
  CALL DEALLOCATE_TransferEntropy()
!
  call MPI_finish()
!
  return
end subroutine STE1_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine STE2_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, debug, &
                      Nshuffles, r0, statP)
  use tempi_class
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!  
      integer :: ierr
      real(sifm_real) :: LocalTime_Start, LocalTime_end
      real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  nframes = nf
  natoms = na
  ndim = nd
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_TransferEntropy(qTeShuffle, 1) 
!
  CALL read_embdparam()
  CALL read_xyzs()
!
! -- Initialize MPI
  call MPI_start()
  call MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, 1)
!
  Localtime_start = MPI_Wtime()
!
  IF (qteShuffle > 0) THEN
      CALL getSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  ELSE
      CALL getSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
  ENDIF
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for TE:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  IF (myID == Master) call WRITE_TE(qTeShuffle,1)
! 
  CALL DEALLOCATE_TransferEntropy()
!
  call MPI_finish()
!
  return
end subroutine STE2_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine STE3_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, debug, &
                      Nshuffles, r0, statP)
  use tempi_class
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!  
      integer :: ierr
      real(sifm_real) :: LocalTime_Start, LocalTime_end
      real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  nframes = nf
  natoms = na
  ndim = nd
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_TransferEntropy(qTeShuffle, 3) 
!
  write(*,*) "Read embd parameters"
  CALL read_embdparam()
  CALL read_xyzs()
  write(*,*) ".................... Done!"
!
! -- Initialize MPI
  call MPI_start()
  call MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, 3)
!
  Localtime_start = MPI_Wtime()
!
  CALL getCrossSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
!
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for TE:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  IF (myID == Master) call WRITE_TE(qTeShuffle,3)
! 
  CALL DEALLOCATE_TransferEntropy()
!
  call MPI_finish()
!
  return
end subroutine STE3_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TE_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, debug, &
                      Nshuffles, r0, statP)
  use tempi_class
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional   :: r0, statP
!  
      integer :: ierr
      real(sifm_real) :: LocalTime_Start, LocalTime_end
      real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  nframes = nf
  natoms = na
  ndim = nd
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_TransferEntropy(qTeShuffle, 4) 
!
  CALL read_embdparam()
  CALL read_xyz()
!
! -- Initialize MPI
  call MPI_start()
  call MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, 2)
!
  Localtime_start = MPI_Wtime()
!
  IF (qteShuffle > 0) THEN
      CALL getTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  ELSE
      CALL getTransferEntropyNDIM_MPIDOF(qteNorm, debug)
  ENDIF
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for TE:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  IF (myID == Master) call WRITE_TE(qTeShuffle,4)
! 
  CALL DEALLOCATE_TransferEntropy()
!
  call MPI_finish()
!
  return
end subroutine te_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the transfer entropy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_TE(qTeShuffle, model)
Implicit None
integer, intent(in) :: qTeShuffle, model
!  
integer             :: i, j, nunit
!  
nunit = new_unit()
open(unit=nunit, file='te.txt', status='unknown', action='write')
if (model == 3) then
DO i=1, Natoms
   write(nunit,*) (TE(i,j), j=1, ndim)
ENDDO
else
DO i=1, Natoms
   write(nunit,*) (TE(i,j), j=1, natoms)
ENDDO
endif
!
IF (qteshuffle > 0) THEN
nunit = new_unit()
open(unit=nunit, file='ste.txt', status='unknown', action='write')
if (model == 3) then
DO i=1, Natoms
   write(nunit,*) (TES(i,j), j=1, ndim)
ENDDO
else
DO i=1, Natoms
   write(nunit,*) (TES(i,j), j=1, natoms)
ENDDO
endif
ENDIF
!
return
!
end subroutine write_te
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Allocate_TransferEntropy(qteShuffle, model)
      implicit none
      integer, intent(in) :: qTeShuffle, model
!
      integer :: ierr
!
      IF (.not. Allocated(TE)) THEN
          if (model == 3) then
	      ALLOCATE( TE(Natoms, Ndim), stat=ierr )
              if (ierr /= 0) stop 'Error allocating TE'
          else
	      ALLOCATE( TE(Natoms, Natoms), stat=ierr )
              if (ierr /= 0) stop 'Error allocating TE'
          endif
      ENDIF
!
      IF (qTEShuffle > 0) THEN
          IF (.not. Allocated(TES)) THEN
	      ALLOCATE( TES(Natoms, Natoms), stat=ierr )
              if (ierr /= 0) stop 'Error allocating TES'
          ENDIF
      ENDIF
!
      IF (model < 4) THEN
          IF (.not. Allocated(xyzs)) THEN
              Allocate(XYZs(Natoms, ndim, Nframes), stat=ierr)
              IF (ierr /= 0) Stop "Error Allocating XYZs"
          ENDIF
      ELSE
          IF (.not. Allocated(xyz)) THEN
              Allocate(XYZ(Natoms, ndim, Nframes), stat=ierr)
              IF (ierr /= 0) Stop "Error Allocating XYZ"
          ENDIF
      ENDIF
      IF (.not. Allocated(topt)) THEN
          Allocate(topt(Natoms, ndim), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating topt"
      ENDIF
      IF (.not. Allocated(mopt)) THEN
          Allocate(mopt(Natoms, ndim), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating mopt"
      ENDIF
!
      return
end subroutine Allocate_TransferEntropy
!
subroutine deAllocate_TransferEntropy()
      implicit none
!
      integer :: ierr
!
      IF (Allocated(TE)) THEN
	  DEALLOCATE( TE, stat=ierr )
          if (ierr /= 0) stop 'Error deallocating TE'
      ENDIF
!
      IF (Allocated(TES)) THEN
	  DEALLOCATE( TES, stat=ierr )
          if (ierr /= 0) stop 'Error deallocating TES'
      ENDIF
!
      IF (Allocated(xyzs)) THEN
	  DEAllocate(xyzs, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating XYZs"
      ENDIF
      IF (Allocated(xyz)) THEN
	  DEAllocate(xyz, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating XYZ"
      ENDIF
      IF (Allocated(mopt)) THEN
	  DEAllocate(mopt, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating mopt"
      ENDIF
      IF (Allocated(topt)) THEN
	  DEAllocate(topt, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating topt"
      ENDIF
!
      return
end subroutine deAllocate_TransferEntropy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xyz()
Implicit None
!  
  integer :: i, d, t, nunit, time, atindex
  character(len=7) :: remark
  character(len=80) :: line  
!
nunit = new_unit()
open(unit=nunit, file='traj.xyz', status='old', action='read')
read(nunit, '(A80)') line
read(nunit, '(A80)') line
DO t = 1, Nframes
   DO i = 1, Natoms
      read(nunit,'(A80)') line
      IF (line(1:5) == 'FRAME') THEN
          read(line,'(A7, 2I10, 3F12.6)') remark, time, atindex, ( XYZ(i,d,t), d = 1, Ndim ) 
      ENDIF
   ENDDO
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit)
!
return
!
end subroutine read_xyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the embedded dimension parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_embdparam()
Implicit None
!  
  integer :: i, d, nunit
!
  nunit = new_unit()
  open(unit=nunit, file='embd.txt', status='old', action='read')
!  
  DO i = 1, Natoms
     read(nunit,*) (topt(i,d), d=1, Ndim), (mopt(i,d), d=1, Ndim)
  ENDDO
  CLOSE(nunit)
!
return
!
end subroutine read_embdparam
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the symbolic Trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xyzs()
Implicit None
!  
  integer :: i, d, t, nunit, time, atindex
  character(len=7) :: remark
  character(len=80) :: line  
!
nunit = new_unit()
open(unit=nunit, file='traj.xyzs', status='old', action='read')
read(nunit, '(A80)') line
read(nunit, '(A80)') line
DO t = 1, Nframes
   DO i = 1, Natoms
      read(nunit,'(A80)') line
      IF (line(1:5) == 'FRAME') THEN
          read(line,'(A7, 2I10, 3A3)') remark, time, atindex, ( XYZS(i,d,t), d = 1, Ndim ) 
      ENDIF
   ENDDO
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit)
!
return
!
end subroutine read_xyzs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MPI_BroadCast(qteShuffle, qteNorm, debug, Nshuffles, r0, statP, model)
use TEMPI_CLASS
implicit none
!
      integer, intent(in)        :: qTeNorm, model
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional   :: r0, statP
!
!
integer :: ierr
!
! --- If running on parallel Broadcast some information
call MPI_BCAST (Natoms,     1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Ndim,       1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Nframes,    1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (qTeNorm,    1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (debug,      1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
if (model == 3) then
call MPI_BCAST (TE,         Natoms*Ndim, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
else
call MPI_BCAST (TE,         Natoms*Natoms, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
endif
call MPI_BCAST (qTEShuffle, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
IF (qTEShuffle > 0) THEN
    call MPI_BCAST (Nshuffles,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (TES,        Natoms*Natoms, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (statP,      1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (r0,         1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
ENDIF
!
if (model < 4) then
    call MPI_BCAST (XYZS,   Natoms*Ndim*Nframes, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr )
else
    call MPI_BCAST (XYZ,    Natoms*Ndim*Nframes, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
endif
!
call MPI_BCAST (Mopt,       Natoms*Ndim, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Topt,       Natoms*Ndim, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
return
END Subroutine MPI_Broadcast
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getCrossSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
use TEMPI_CLASS
implicit none
integer, intent(in)        :: qTeNorm
integer, intent(in)        :: debug
!
!{Local variables}
integer :: i, iList, K
real(sifm_real) :: txy, tyx, hxx, hyy, hx, hy, &
                   time_start, time_end
!
integer :: ierr
!
real(sifm_real), dimension(:), allocatable :: LocalTij, LocalTji
real(sifm_real), dimension(:), allocatable :: GlobalTij,GlobalTji
integer :: iStart, iEnd, Stride
!
IF ( mod(Natoms, NumProcs) == 0 ) THEN
     Stride = Natoms / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LocalTij(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Tij'
     ALLOCATE( LocalTji(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Tji'
     ALLOCATE( GlobalTij(Natoms), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Tij'
     ALLOCATE( GlobalTji(Natoms), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Tji'
!
ELSE
     Stop "Number of tpairs must be Multiple of Number of Processors"
ENDIF
!
! --- Compute the entropy transfer
!-------- call all the pairs
TE = rzero
Pairs_LOOP: DO i = iStart, iEnd
      K = i - iStart + 1
!
      IF (myID == Master) call cpu_time(time_start)
!
   !--- Compute the direct transfer entropy (only on master slave)
      CALL symbolic_TE_entropies1D(nframes,XYZs(i,1,:),XYZs(i,2,:),                   &
	                               Mopt(i,1),Mopt(i,2),Topt(i,1),Topt(i,2),       &
				       Txy,Tyx,hxx,hyy,hx,hy)   
      IF (qTeNorm > 0) THEN
          IF ( abs(Hyy) >  prtiny ) THEN
   	       Tyx = Tyx/Hyy
          ELSE
   	       Tyx = Tyx/(Hx+Hy)
          ENDIF
          IF ( abs(Hxx) >  prtiny ) THEN
   	       Txy = Txy/Hxx
          ELSE
   	       Txy = Txy/(Hx+Hy)
          ENDIF
      ENDIF 
!
      LocalTij(K) = Txy
      LocalTji(K) = Tyx	
!
      IF ( debug > 0 ) THEN
          write(*,'("--------------------------------------------------------")')
          write(*,'("Instantaneous Txy (+bias)  =        ", 2I5, F12.6)') K, Stride, LocalTij(K)
          write(*,'("Instantaneous Tyx (+bias)  =        ", 2I5, F12.6)') K, Stride, LocalTji(K)
      ENDIF
		 
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F20.3)') Time_End - Time_Start
      ENDIF
!
enddo Pairs_LOOP
!
CALL MPI_gather( LocalTij,  Stride, MPI_DOUBLE_PRECISION, &
                 GlobalTij, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
CALL MPI_gather( LocalTji,  Stride, MPI_DOUBLE_PRECISION, &
                 GlobalTji, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )

IF (myID == master) THEN
    DO i = 1, Natoms
       TE(i,1) = GlobalTij(i)
       TE(i,2) = GlobalTji(i)
    ENDDO
ENDIF
!
DEALLOCATE( LocalTij, stat=ierr )
if (ierr /= 0) stop 'Error allocating Local Tij'
DEALLOCATE( LocalTji, stat=ierr )
if (ierr /= 0) stop 'Error allocating Local Tji'
DEALLOCATE( GlobalTij, stat=ierr )
if (ierr /= 0) stop 'Error allocating Global Tij'
DEALLOCATE( GlobalTji, stat=ierr )
if (ierr /= 0) stop 'Error allocating Global Tji'
!
return
end subroutine getCrossSTransferEntropyNDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
use TEMPI_CLASS
implicit none
integer, intent(in)        :: qTeNorm
integer, intent(in)        :: debug
!
!{Local variables}
integer, allocatable :: listPairs(:,:)
integer :: Npairs
integer :: i, j, K, Ilist
real(sifm_real) :: txy, tyx, hxx, hyy, hx, hy, &
                   time_start, time_end
!
integer :: ierr
!
real(sifm_real), dimension(:), allocatable :: LocalTij, LocalTji
real(sifm_real), dimension(:), allocatable :: GlobalTij,GlobalTji
integer :: iStart, iEnd, Stride

!--- Compute total number of possible pairs
Npairs = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      Npairs = Npairs + 1
   ENDDO
ENDDO
IF ( .not. Allocated(listPairs) ) THEN
     ALLOCATE( listPairs(Npairs, 2), stat=ierr )
     if (ierr /= 0) stop 'Error allocating List of Pairs'
ENDIF
k = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      k = k + 1
      listPairs(k, 1) = i
      listPairs(k, 2) = j
   ENDDO
ENDDO
!
IF ( mod(Npairs, NumProcs) == 0 ) THEN
     Stride = Npairs / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LocalTij(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Tij'
     ALLOCATE( LocalTji(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Tji'
     ALLOCATE( GlobalTij(Npairs), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Tij'
     ALLOCATE( GlobalTji(Npairs), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Tji'
!
ELSE
     Stop "Number of tpairs must be Multiple of Number of Processors"
ENDIF
!
! --- Compute the entropy transfer
!-------- call all the pairs
TE = rzero
Pairs_LOOP: DO iList = iStart, iEnd
      K = iList - iStart + 1
!
      IF (myID == Master) call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
   !--- Compute the direct transfer entropy (only on master slave)
      IF (Ndim > 1) THEN
          CALL symbolic_TE_entropiesND(nframes,ndim,XYZs(i,:,:),XYZs(j,:,:),        &
	                               Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),       &
				       Txy,Tyx,hxx,hyy,hx,hy)   
      ELSE
          CALL symbolic_TE_entropies1D(nframes, XYZs(i,1,:),XYZs(j,1,:),            &
	                               Mopt(i,1),Mopt(j,1),Topt(i,1),Topt(j,1),       &
				       Txy,Tyx,hxx,hyy,Hx,Hy)  
      ENDIF
      IF (qTeNorm > 0) THEN
          IF ( abs(Hyy) >  prtiny ) THEN
   	       Tyx = Tyx/Hyy
          ELSE
   	       Tyx = Tyx/(Hx+Hy)
          ENDIF
          IF ( abs(Hxx) >  prtiny ) THEN
   	       Txy = Txy/Hxx
          ELSE
   	       Txy = Txy/(Hx+Hy)
          ENDIF
      ENDIF 
!
      LocalTij(K) = Txy
      LocalTji(K) = Tyx	
!
      IF ( debug > 0 ) THEN
          write(*,'("--------------------------------------------------------")')
          write(*,'("Instantaneous Txy (+bias)  =        ", 2I5, F12.6)') K, Stride, LocalTij(K)
          write(*,'("Instantaneous Tyx (+bias)  =        ", 2I5, F12.6)') K, Stride, LocalTji(K)
      ENDIF
		 
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F20.3)') Time_End - Time_Start
      ENDIF
!
enddo Pairs_LOOP
!
CALL MPI_gather( LocalTij,  Stride, MPI_DOUBLE_PRECISION, &
                 GlobalTij, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
CALL MPI_gather( LocalTji,  Stride, MPI_DOUBLE_PRECISION, &
                 GlobalTji, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )

IF (myID == master) THEN
    DO iList = 1, Npairs
       I = ListPairs(iList, 1)
       J = ListPairs(iList, 2)
       TE(i,j) = GlobalTij(iList)
       TE(j,i) = GlobalTji(iList)
    ENDDO
ENDIF
!
IF ( Allocated(listPairs) ) THEN
     DEALLOCATE( listPairs, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Pairs'
ENDIF
!
DEALLOCATE( LocalTij, stat=ierr )
if (ierr /= 0) stop 'Error allocating Local Tij'
DEALLOCATE( LocalTji, stat=ierr )
if (ierr /= 0) stop 'Error allocating Local Tji'
DEALLOCATE( GlobalTij, stat=ierr )
if (ierr /= 0) stop 'Error allocating Global Tij'
DEALLOCATE( GlobalTji, stat=ierr )
if (ierr /= 0) stop 'Error allocating Global Tji'
!
return
end subroutine getSTransferEntropyNDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
use TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!
!{Local variables}
integer :: ntrials
integer, allocatable :: listPairs(:,:)
integer :: Npairs
integer :: i, j, iTrial, t, iList
character(len=1), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real) :: txy, tyx, hxx, hyy, hx, hy
!
integer :: ierr, Stride, K, iStart, iEnd
real(sifm_real), allocatable :: GlobalTE(:), LocalTE(:)

!--- Compute total number of possible pairs
Npairs = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      Npairs = Npairs + 1
   ENDDO
ENDDO
ALLOCATE( listPairs(Npairs, 2), stat=ierr )
if (ierr /= 0) stop 'Error allocating List of Pairs'
t = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      t = t + 1
      listPairs(t, 1) = i
      listPairs(t, 2) = j
   ENDDO
ENDDO
!
! -- Initialize 2*natoms random number sequences
IF ( qteshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, MyID)
     ntrials = nshuffles
endif
!
IF ( mod(Ntrials, NumProcs) == 0 ) THEN
     Stride = Ntrials / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LocalTE(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local TExy'
     ALLOCATE( GlobalTE(Ntrials), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global TExy'
!
ELSE
     Stop "Number of Trial Shuffles must be Multiple of Number of Processors"
ENDIF

! --- Compute the entropy transfer
TE = rzero; TES = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
   !--- Compute the direct transfer entropy (only on master slave)
      IF (myID == Master) THEN
          IF (Ndim > 1) THEN
              CALL symbolic_TE_entropiesND(nframes,ndim,XYZs(i,:,:),XYZs(j,:,:),        &
	                                     Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),       &
				             Txy,Tyx,hxx,hyy,hx,hy)   
          ELSE
              CALL symbolic_TE_entropies1D(nframes, XYZs(i,1,:),XYZs(j,1,:),            &
	                                     Mopt(i,1),Mopt(j,1),Topt(i,1),Topt(j,1),       &
				             Txy,Tyx,hxx,hyy,Hx,Hy)  
          ENDIF
          TE(i,j) = Txy
          TE(j,i) = Tyx
          write(*,'("--------------------------------------------------------")')
          write(*,'("Instantaneous Txy (+bias)  =        ", F12.6)') TE(i,j)
          write(*,'("Instantaneous Tyx (+bias)  =        ", F12.6)') TE(j,i)
      ENDIF
!---- Remove the bias, by shuffling the trajectories
      IF (qteshuffle > 0) THEN

       ! -- Subtract the bias from T(x->y) and normalize it
       xs_shuffle = XYZS(i,:,:)
       LocalTE = rzero
       DO itrial = iStart, iEnd
          K = iTrial - iStart + 1
          IF (qTEShuffle == 1) THEN
              CALL shuffleND(nframes,ndim,xs_shuffle,I)
          ELSE
              CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(i,:), Topt(i,:), I)
          ENDIF
          CALL symbolic_TE_entropyND(nframes,ndim,xs_shuffle,XYZs(j,:,:),          &
		                     Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),Txy)
          LocalTE(K) = Txy
       ENDDO
       CALL MPI_gather( LocalTE,  Stride, MPI_DOUBLE_PRECISION,  &
                        GlobalTE, Stride, MPI_DOUBLE_PRECISION, &
                        master, MPI_COMM_WORLD, ierr )             
       IF (myID == Master) THEN
           Txy = mean( GlobalTE, Ntrials )
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of Bias Txy (Before Statistical Test)  =  ", F12.6)') Txy
           CALL statTest(Ntrials, R0, Txy, Txy, StatP)	
           write(*,'("Average value of Bias Txy (After Statistical Test)   =  ", F12.6)') Txy		 
           IF (qTENorm == 0) THEN
	       TES(i,j) = TE(i,j) - Txy
           ELSE
               IF ( abs(Hxx) > prtiny ) THEN
   	            TES(i,j) = (TE(i,j) - Txy)/Hxx
               ELSE
   	            TES(i,j) = (TE(i,j) - Txy)/(Hx+Hy)
               ENDIF
           ENDIF
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of true TSxy (Before Statistical Test) =  ", 3F12.6)') TE(i,j), Txy, TES(i,j)
           CALL statTest(Ntrials, R0, TES(i,j), TES(i,j), StatP)	
           write(*,'("Average value of true TSxy (After Statistical Test)  =  ", F12.6)') TES(i,j)
           write(*,'("--------------------------------------------------------")')
           IF ( debug == 1 ) THEN
                write(*,'("Txy (+bias, -bias) =  ", 5F12.6)') TE(i,j), TES(i,j), Hxx,Hx,Hy
           Endif
       ENDIF

       ! -- Subtract the bias from T(y->x) and normalize it
       xs_shuffle = XYZS(j,:,:)
       LocalTE = rzero
       DO itrial = iStart, iEnd
          K = iTrial - iStart + 1
          IF (qTEShuffle == 1) THEN
              CALL shuffleND(nframes,ndim,xs_shuffle,Natoms+J)
          ELSE
              CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
          ENDIF
          CALL symbolic_TE_entropyND(nframes,ndim,xs_shuffle, XYZs(i,:,:), &
                                     Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),Tyx)
          LocalTE(K) = Tyx
       ENDDO	  
       CALL MPI_gather( LocalTE,  Stride, MPI_DOUBLE_PRECISION,  &
                        GlobalTE, Stride, MPI_DOUBLE_PRECISION, &
                        master, MPI_COMM_WORLD, ierr )             
       IF (myID == Master) THEN
           Tyx = mean( GlobalTE, Ntrials )
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of Bias Tyx (Before Statistical Test)  =  ", F12.6)') Tyx			 
           CALL statTest(Ntrials, R0, Tyx, Tyx, STATP)
           write(*,'("Average value of Bias Tyx (After Statistical Test)   =  ", F12.6)') Tyx
           write(*,'("--------------------------------------------------------")')			 
           IF (qTeNorm == 0) THEN
	       TES(j,i) = TE(j,i) - Tyx
           ELSE 
               IF ( abs(Hyy) >  prtiny ) THEN
   	            TES(j,i) = (TE(j,i) - Tyx)/Hyy
               ELSE
   	            TES(j,i) = (TE(j,i) - Tyx)/(Hx+Hy)
               ENDIF
           ENDIF 
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of true TSyx (Before Statistical Test) =  ", F12.6)') TES(j,i)
           CALL statTest(Ntrials, R0, TES(j,i), TES(j,i), StatP)	
           write(*,'("Average value of true TSyx (After Statistical Test)  =  ", F12.6)') TES(j,i)
           write(*,'("--------------------------------------------------------")')
           IF ( debug == 1 ) THEN
                write(*,'("Tyx (+Bias, -Bias) =  ", 5F12.6)') TE(j,i), TES(j,i), Hyy,Hx,Hy
           Endif
       ENDIF
      ENDIF
!
enddo Pairs_LOOP
!
IF ( qteshuffle > 0) then
     CALL rndnum_clear()
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
DEALLOCATE( LocalTE, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Local TE'
DEALLOCATE( GlobalTE, stat=ierr )
if (ierr /= 0) stop 'Error deallocating global TE'
!
return
end subroutine getSTransferEntropyNDIM_MPISFL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getTransferEntropyNDIM_MPIDOF(qteNorm, debug)
use TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: debug
!
!{Local variables}
integer, allocatable :: listPairs(:,:)
integer :: Npairs
integer :: i, j, K, Ilist
real(sifm_real) :: txy, tyx, hxx, hyy, hx, hy, &
                   time_start, time_end
!
integer :: ierr
!
real(sifm_real), dimension(:), allocatable :: LocalTij, LocalTji
real(sifm_real), dimension(:), allocatable :: GlobalTij,GlobalTji
integer :: iStart, iEnd, Stride

!--- Compute total number of possible pairs
Npairs = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      Npairs = Npairs + 1
   ENDDO
ENDDO
IF ( .not. Allocated(listPairs) ) THEN
     ALLOCATE( listPairs(Npairs, 2), stat=ierr )
     if (ierr /= 0) stop 'Error allocating List of Pairs'
ENDIF
k = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      k = k + 1
      listPairs(k, 1) = i
      listPairs(k, 2) = j
   ENDDO
ENDDO
!
IF ( mod(Npairs, NumProcs) == 0 ) THEN
     Stride = Npairs / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LocalTij(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Tij'
     ALLOCATE( LocalTji(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Tji'
     ALLOCATE( GlobalTij(Npairs), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Tij'
     ALLOCATE( GlobalTji(Npairs), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Tji'
!
ELSE
     Stop "Number of tpairs must be Multiple of Number of Processors"
ENDIF
!
! --- Compute the entropy transfer
TE = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = iStart, iEnd
      K = iList - iStart + 1
!
      IF (myID == Master) call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
   !--- Compute the direct transfer entropy (only on master slave)
      IF (Ndim > 1) THEN
          CALL symbolic_TE2_entropiesND(nframes,ndim,XYZ(i,:,:),XYZ(j,:,:),        &
	                               Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),       &
				       Txy,Tyx,hxx,hyy,hx,hy)   
      ELSE
          CALL symbolic_TE2_entropies1D(nframes, XYZ(i,1,:),XYZ(j,1,:),            &
	                               Mopt(i,1),Mopt(j,1),Topt(i,1),Topt(j,1),       &
				       Txy,Tyx,hxx,hyy,Hx,Hy)  
      ENDIF
      IF (qTeNorm > 0) THEN
          IF ( abs(Hyy) >  prtiny ) THEN
   	       Tyx = Tyx/Hyy
          ELSE
   	       Tyx = Tyx/(Hx+Hy)
          ENDIF
          IF ( abs(Hxx) >  prtiny ) THEN
   	       Txy = Txy/Hxx
          ELSE
   	       Txy = Txy/(Hx+Hy)
          ENDIF
      ENDIF 
!
      LocalTij(K) = Txy
      LocalTji(K) = Tyx	
!
      IF (debug > 0) THEN
          write(*,'("--------------------------------------------------------")')
          write(*,'("Instantaneous Txy (+bias)  =        ", 2I5, F12.6)') K, Stride, LocalTij(K)
          write(*,'("Instantaneous Tyx (+bias)  =        ", 2I5, F12.6)') K, Stride, LocalTji(K)
      ENDIF
		 
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F20.3)') Time_End - Time_Start
      ENDIF
!
enddo Pairs_LOOP
!
CALL MPI_gather( LocalTij,  Stride, MPI_DOUBLE_PRECISION, &
                 GlobalTij, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
CALL MPI_gather( LocalTji,  Stride, MPI_DOUBLE_PRECISION, &
                 GlobalTji, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )

IF (myID == master) THEN
    DO iList = 1, Npairs
       I = ListPairs(iList, 1)
       J = ListPairs(iList, 2)
       TE(i,j) = GlobalTij(iList)
       TE(j,i) = GlobalTji(iList)
    ENDDO
ENDIF
!
IF ( Allocated(listPairs) ) THEN
     DEALLOCATE( listPairs, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Pairs'
ENDIF
!
DEALLOCATE( LocalTij, stat=ierr )
if (ierr /= 0) stop 'Error allocating Local Tij'
DEALLOCATE( LocalTji, stat=ierr )
if (ierr /= 0) stop 'Error allocating Local Tji'
DEALLOCATE( GlobalTij, stat=ierr )
if (ierr /= 0) stop 'Error allocating Global Tij'
DEALLOCATE( GlobalTji, stat=ierr )
if (ierr /= 0) stop 'Error allocating Global Tji'
!
return
end subroutine getTransferEntropyNDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
use TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional   :: r0, statP
!
!{Local variables}
integer :: ntrials
integer, allocatable :: listPairs(:,:)
integer :: Npairs
integer :: i, j, iTrial, t, iList
real(sifm_real), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real) :: txy, tyx, hxx, hyy, hx, hy
!
integer :: ierr, Stride, K, iStart, iEnd
real(sifm_real), allocatable :: GlobalTE(:), LocalTE(:)

!--- Compute total number of possible pairs
Npairs = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      Npairs = Npairs + 1
   ENDDO
ENDDO
ALLOCATE( listPairs(Npairs, 2), stat=ierr )
if (ierr /= 0) stop 'Error allocating List of Pairs'
t = 0
DO i = 1, Natoms-1
   !-------- call the second atom
   DO j = i+1, Natoms
      t = t + 1
      listPairs(t, 1) = i
      listPairs(t, 2) = j
   ENDDO
ENDDO
!
! -- Initialize 2*natoms random number sequences
IF ( qteshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, MyID)
     ntrials = nshuffles
endif
!
IF ( mod(Ntrials, NumProcs) == 0 ) THEN
     Stride = Ntrials / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LocalTE(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local TExy'
     ALLOCATE( GlobalTE(Ntrials), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global TExy'
!
ELSE
     Stop "Number of Trial Shuffles must be Multiple of Number of Processors"
ENDIF

! --- Compute the entropy transfer
TE = rzero; TES = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
   !--- Compute the direct transfer entropy (only on master slave)
      IF (myID == Master) THEN
          IF (Ndim > 1) THEN
              CALL symbolic_TE2_entropiesND(nframes,ndim,XYZ(i,:,:),XYZ(j,:,:),        &
	                                    Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),   &
				            Txy,Tyx,hxx,hyy,hx,hy)   
          ELSE
              CALL symbolic_TE2_entropies1D(nframes, XYZ(i,1,:),XYZ(j,1,:),            &
	                                    Mopt(i,1),Mopt(j,1),Topt(i,1),Topt(j,1),   &
				            Txy,Tyx,hxx,hyy,Hx,Hy)  
          ENDIF
          TE(i,j) = Txy
          TE(j,i) = Tyx
          write(*,'("--------------------------------------------------------")')
          write(*,'("Instantaneous Txy (+bias)  =        ", F12.6)') TE(i,j)
          write(*,'("Instantaneous Tyx (+bias)  =        ", F12.6)') TE(j,i)
      ENDIF
!---- Remove the bias, by shuffling the trajectories
      IF (qteshuffle > 0) THEN

       ! -- Subtract the bias from T(x->y) and normalize it
       xs_shuffle = XYZ(i,:,:)
       LocalTE = rzero
       DO itrial = iStart, iEnd
          K = iTrial - iStart + 1
          IF (qTEShuffle == 1) THEN
              CALL shuffleND2(nframes,ndim,xs_shuffle,I)
          ELSE
              CALL BlockShuffleND2(nframes, ndim, xs_shuffle, Mopt(i,:), Topt(i,:), I)
          ENDIF
          CALL symbolic_TE2_entropyND(nframes,ndim,xs_shuffle,XYZ(j,:,:),          &
		                     Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),Txy)
          LocalTE(K) = Txy
       ENDDO
       CALL MPI_gather( LocalTE,  Stride, MPI_DOUBLE_PRECISION,  &
                        GlobalTE, Stride, MPI_DOUBLE_PRECISION, &
                        master, MPI_COMM_WORLD, ierr )             
       IF (myID == Master) THEN
           Txy = mean( GlobalTE, Ntrials )
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of Bias Txy (Before Statistical Test)  =  ", F12.6)') Txy
           CALL statTest(Ntrials, R0, Txy, Txy, StatP)	
           write(*,'("Average value of Bias Txy (After Statistical Test)   =  ", F12.6)') Txy		 
           IF (qTENorm == 0) THEN
	       TES(i,j) = TE(i,j) - Txy
           ELSE
               IF ( abs(Hxx) > prtiny ) THEN
   	            TES(i,j) = (TE(i,j) - Txy)/Hxx
               ELSE
   	            TES(i,j) = (TE(i,j) - Txy)/(Hx+Hy)
               ENDIF
           ENDIF
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of true TSxy (Before Statistical Test) =  ", 3F12.6)') TE(i,j), Txy, TES(i,j)
           CALL statTest(Ntrials, R0, TES(i,j), TES(i,j), StatP)	
           write(*,'("Average value of true TSxy (After Statistical Test)  =  ", F12.6)') TES(i,j)
           write(*,'("--------------------------------------------------------")')
           IF ( debug == 1 ) THEN
                write(*,'("Txy (+bias, -bias) =  ", 5F12.6)') TE(i,j), TES(i,j), Hxx,Hx,Hy
           Endif
       ENDIF

       ! -- Subtract the bias from T(y->x) and normalize it
       xs_shuffle = XYZ(j,:,:)
       LocalTE = rzero
       DO itrial = iStart, iEnd
          K = iTrial - iStart + 1
          IF (qTEShuffle == 1) THEN
              CALL shuffleND2(nframes,ndim,xs_shuffle,Natoms+J)
          ELSE
              CALL BlockShuffleND2(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
          ENDIF
          CALL symbolic_TE2_entropyND(nframes,ndim,xs_shuffle, XYZ(i,:,:),           &
                                      Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),Tyx)
          LocalTE(K) = Tyx
       ENDDO	  
       CALL MPI_gather( LocalTE,  Stride, MPI_DOUBLE_PRECISION,  &
                        GlobalTE, Stride, MPI_DOUBLE_PRECISION, &
                        master, MPI_COMM_WORLD, ierr )             
       IF (myID == Master) THEN
           Tyx = mean( GlobalTE, Ntrials )
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of Bias Tyx (Before Statistical Test)  =  ", F12.6)') Tyx			 
           CALL statTest(Ntrials, R0, Tyx, Tyx, STATP)
           write(*,'("Average value of Bias Tyx (After Statistical Test)   =  ", F12.6)') Tyx
           write(*,'("--------------------------------------------------------")')			 
           IF (qTeNorm == 0) THEN
	       TES(j,i) = TE(j,i) - Tyx
           ELSE 
               IF ( abs(Hyy) >  prtiny ) THEN
   	            TES(j,i) = (TE(j,i) - Tyx)/Hyy
               ELSE
   	            TES(j,i) = (TE(j,i) - Tyx)/(Hx+Hy)
               ENDIF
           ENDIF 
           write(*,'("--------------------------------------------------------")')
           write(*,'("Average value of true TSyx (Before Statistical Test) =  ", F12.6)') TES(j,i)
           CALL statTest(Ntrials, R0, TES(j,i), TES(j,i), StatP)	
           write(*,'("Average value of true TSyx (After Statistical Test)  =  ", F12.6)') TES(j,i)
           write(*,'("--------------------------------------------------------")')
           IF ( debug == 1 ) THEN
                write(*,'("Tyx (+Bias, -Bias) =  ", 5F12.6)') TE(j,i), TES(j,i), Hyy,Hx,Hy
           Endif
       ENDIF
      ENDIF
!
enddo Pairs_LOOP
!
IF ( qteshuffle > 0) then
     CALL rndnum_clear()
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
DEALLOCATE( LocalTE, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Local TE'
DEALLOCATE( GlobalTE, stat=ierr )
if (ierr /= 0) stop 'Error deallocating global TE'
!
return
end subroutine getTransferEntropyNDIM_MPISFL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropies1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)           :: ndata,m1,m2,tau1,tau2
character(len=1), intent(in)  :: xs(:),ys(:)
real(sifm_real), intent(out)  :: Txy,Tyx
real(sifm_real), intent(out)  :: Hxx,Hyy,hx,hy
!
type (StateElem), pointer  :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer  :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter         :: Delta=1
integer                    :: istat,L,k,m,tau
character(len=1)           :: xsc,ysc
character(len=250)         :: xL,x1L,yL,y1L
character(len=250)         :: xyL, xy1L, yx1L
real(sifm_real)            :: hx1,hy1,hxy,hxy1,hyx1
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m=max(m1,m2);tau=max(tau1,tau2)
!!! --- Loop over all data points
time_Loop: DO L=((M-1)*Tau+1),(Ndata-Delta)
  xL=""; yL=""
  x1L=""; y1L=""
  xyL=""; xy1L=""; yx1L=""
!
  x_Loop: DO k=1,m1
          xsc  = xs(L-(k-1)*tau1)
	  call addStrings(xL,xsc)
          call addStrings(x1L,xsc)
          call addStrings(xyL,xsc)
   	  call addStrings(xy1L,xsc)
          call addStrings(yx1L,xsc)
  end do x_Loop
!
  y_Loop: DO k=1,m2
          ysc  = ys(L-(k-1)*tau2)
	  call addStrings(yL,ysc)
  	  call addStrings(y1L,ysc)
	  call addStrings(xyL,ysc)
	  call addStrings(xy1L,ysc)
   	  call addStrings(yx1L,ysc)
  end do y_Loop
!
  ysc  = ys(L+Delta)
  call addStrings(y1L,ysc)
  call addStrings(xy1L, ysc)
!
  xsc  = xs(L+Delta)
  call addStrings(x1L,xsc)
  call addStrings(yx1L, xsc)
!
  call Linklist(headx,tailx,xL)   
  call LinkList(heady,taily,yL)
  call Linklist(headx1,tailx1,x1L)   
  call LinkList(heady1,taily1,y1L)
  call Linklist(headxy,tailxy,xyL)
  call LinkList(headxy1,tailxy1,xy1L)
  call Linklist(headyx1,tailyx1,yx1L)   
!
end do time_Loop
!!! --- calculate the Shannon entropy
Hx = SymbolicShannonEntropy(Headx)
call freeList(headx)
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hx1 = SymbolicShannonEntropy(Headx1)
call freeList(headx1)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
Hyx1 = SymbolicShannonEntropy(Headyx1)
call freeList(headyx1)
!
Txy=(hy1-hy) - (hxy1-hxy)
Tyx=(hx1-hx) - (hyx1-hxy)
Hxx=(hy1-hy)
Hyy=(hx1-hx)
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
return
!
end subroutine symbolic_TE_entropies1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropiesND(ndata,ndim,xs,ys,m1,m2,tau1,tau2, &
                                   Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)           :: ndata,ndim
integer, dimension(ndim)      :: m1,m2,tau1,tau2
character(len=1), intent(in)  :: xs(:,:),ys(:,:)
real(sifm_real), intent(out)  :: Txy,Tyx
real(sifm_real), intent(out)  :: Hxx,Hyy,hx,hy
!
type (StateElem), pointer  :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer  :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter         :: Delta = 1
integer                    :: d,istat,L,k,m,tau
character(len=1)           :: sc
character(len=250)         :: xL,x1L,yL,y1L
character(len=250)         :: xyL, xy1L, yx1L
real(sifm_real)            :: hx1,hy1,hxy,hxy1,hyx1
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m = max(imax1d(m1),imax1d(m2))
tau = max(imax1d(tau1),imax1d(tau2))
!!! --- Loop over all data points
time_Loop: DO L=((M-1)*Tau+1),(Ndata-Delta)
  xL=""; yL=""
  x1L=""; y1L=""
  xyL=""; xy1L=""; yx1L=""
!  
  xdim_Loop1: DO d=1,Ndim 
      x_Loop1: DO k=1,m1(d)
            sc   = xs(d, L-(k-1)*tau1(d))
	    call addStrings(xL,sc)
	    call addStrings(x1L,sc)
	    call addStrings(xyL,sc)
	    call addStrings(xy1L,sc)
            call addStrings(yx1L, sc)
      end do x_Loop1
      sc   = xs(d,L+Delta)
      call addStrings(x1L,sc)
      call addStrings(yx1L, sc)
  end do xdim_Loop1
!  
  ydim_Loop1: DO d=1, Ndim
       y_Loop1: DO k=1,m2(d)
            sc   = ys(d,L-(k-1)*tau2(d))
	    call addStrings(yL, sc)
	    call addStrings(y1L, sc)
	    call addStrings(xyL, sc)
	    call addStrings(xy1L,sc)
	    call addStrings(yx1L, sc)
       end do y_Loop1
       sc   = ys(d,L+Delta)
       call addStrings(y1L,sc)
       call addStrings(xy1L, sc)
  end do ydim_Loop1
!  
  call Linklist(headx,tailx,xL)   
  call LinkList(heady,taily,yL)
  call Linklist(headx1,tailx1,x1L)   
  call LinkList(heady1,taily1,y1L)
  call Linklist(headxy,tailxy,xyL)
  call LinkList(headxy1,tailxy1,xy1L)
  call Linklist(headyx1,tailyx1,yx1L)   

end do time_Loop
!!! --- calculate the Shannon entropy
Hx = SymbolicShannonEntropy(Headx)
call freeList(headx)
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hx1 = SymbolicShannonEntropy(Headx1)
call freeList(headx1)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
Hyx1 = SymbolicShannonEntropy(Headyx1)
call freeList(headyx1)
!
Txy=(hy1-hy) - (hxy1-hxy)
Tyx=(hx1-hx) - (hyx1-hxy)
Hxx=(hy1-hy)
Hyy=(hx1-hx)
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
return
!
end subroutine symbolic_TE_entropiesND   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropy1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)           :: ndata,m1,m2,tau1,tau2
character(len=1), intent(in)  :: xs(:),ys(:)
real(sifm_real), intent(out)  :: Txy
!
type (StateElem), pointer :: heady,taily,heady1,taily1
type (StateElem), pointer :: headxy,tailxy,headxy1,tailxy1
integer, parameter        :: Delta=1
integer                   :: L,k,m,tau
character(len=1)          :: sc
character(len=250)        :: yL,yy1L
character(len=250)        :: xyL, xyy1L
real(sifm_real)           :: hy,hy1,hxy,hxy1
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(m1,m2)
tau = max(tau1,tau2)
!
Txy = rzero
time_Loop: DO L=((M-1)*Tau+1),(Ndata-Delta)
!
   yL ="";  yy1L =""
   xyL="";  xyy1L=""
!
   x_Loop: DO k=1,M1
           sc    = xs(L-(k-1)*Tau1)
           call addStrings(xyL, sc)
	   call addStrings(xyy1L, sc)
   end do x_Loop
!
   y_Loop: DO k=1,M2
           sc      = ys(L-(k-1)*Tau2)
	   call addStrings(yL,sc)
	   call addStrings(xyL, sc)
	   call addStrings(yy1L,sc)
           call addStrings(xyy1L, sc)
   end do y_Loop
!
   sc    = ys(L+Delta)
   call addStrings(yy1L,sc)
   call addStrings(xyy1L, sc)
! 
   call LinkList(heady,taily,yL)  
   call LinkList(heady1,taily1,yy1L)
   call Linklist(headxy,tailxy,xyL)
   call LinkList(headxy1,tailxy1,xyy1L)
!
end do time_Loop
!
!!! --- calculate the Shannon entropy
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
!
Txy=(hy1-hy) - (hxy1-hxy)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
return
!
end subroutine symbolic_TE_entropy1D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropyND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)        :: ndata,ndim
integer, dimension(Ndim), intent(in) :: m1,m2,tau1,tau2
character(len=1), dimension(Ndim,Ndata), intent(in) :: xs,ys
real(sifm_real), intent(out) :: Txy
!
type (StateElem), pointer :: heady,taily,heady1,taily1
type (StateElem), pointer :: headxy,tailxy,headxy1,tailxy1
integer, parameter        :: Delta = 1
integer                   :: d,L,k,m,tau
character(len=1)          :: sc
character(len=250)        :: yL,y1L
character(len=250)        :: xyL, xy1L
real(sifm_real)           :: hy,hy1,hxy,hxy1
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(imax1d(m1), imax1d(m2))
tau = max(imax1d(tau1), imax1d(tau2))
!
Txy = rzero
time_Loop: DO L=((M-1)*Tau+1),(Ndata-Delta)
   yL=""; y1L=""
   xyL=""; xy1L=""
!   
   xdim_Loop: DO d=1, Ndim
   x_Loop: DO k=1,m1(d)
             sc   = xs(d,L-(k-1)*tau1(d))
	     call addStrings(xyL, sc)
	     call addStrings(xy1L, sc)
   end do x_Loop
   end do xdim_Loop
!   
   ydim_Loop: DO d=1, Ndim
   y_Loop: DO k=1,m2(d)
           sc  = ys(d,L-(k-1)*tau2(d))
	   call addStrings(yL, sc)
           call addStrings(xyL, sc)
           call addStrings(y1L, sc)
           call addStrings(xy1L, sc)
   end do y_Loop
   sc    = ys(d,L+Delta)
   call addStrings(y1L, sc)
   call addStrings(xy1L, sc)
   end do ydim_Loop
!   
   call LinkList(heady,taily,yL)  
   call LinkList(heady1,taily1,y1L)
   call Linklist(headxy,tailxy,xyL)
   call LinkList(headxy1,tailxy1,xy1L)

end do time_Loop
!!! --- calculate the Shannon entropy
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
!
Txy=(hy1-hy) - (hxy1-hxy)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
return
!
end subroutine symbolic_TE_entropyND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE2_entropies1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)           :: ndata,m1,m2,tau1,tau2
real(sifm_real), intent(in)    :: xs(:),ys(:)
real(sifm_real), intent(out)   :: Txy,Tyx
real(sifm_real), intent(out)   :: Hxx,Hyy,hx,hy
!
type (StateElem), pointer  :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer  :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter         :: Delta = 1
integer                    :: istat,L,k,m,tau
character(len=1)           :: xsc,ysc
character(len=250)         :: xL,x1L,yL,y1L
character(len=250)         :: xyL, xy1L, yx1L
real(sifm_real)             :: hx1,hy1,hxy,hxy1,hyx1
integer, dimension(m1)     :: xsL
integer, dimension(m1+1)   :: xsL1
integer, dimension(m2)     :: ysL
integer, dimension(m2+1)   :: ysL1
real(sifm_real), dimension(M1)    :: X
real(sifm_real), dimension(M2)    :: Y
real(sifm_real), dimension(M1+1)  :: X1
real(sifm_real), dimension(M2+1)  :: Y1
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m=max(m1,m2);tau=max(tau1,tau2)

!!! --- Loop over all data points
time_Loop: DO L=(M-1)*Tau+1,(Ndata-Delta)
  xL=""; yL=""
  x1L=""; y1L=""
  xyL=""; xy1L=""; yx1L=""
!
  X = vcopy(XS,L,M1,Tau1,-1)
  CALL sort(M1, X, XSL)
  x_Loop: DO k=1,m1
          xsc  = toString( xsL(k) )
	  call addStrings(xL,xsc)
          call addStrings(xyL,xsc)
   	  call addStrings(xy1L,xsc)
  end do x_Loop
!
  X1(1:M1) = X(1:M1)
  X1((m1+1):(m1+1)) = XS(L+Delta) 
  CALL sort(M1+1, X1, XSL1)
  x1_Loop: DO k=1, M1+1
            xsc  = ToString( XSL1(k) )
 	    call addStrings(x1L, xsc)
            call addStrings(yx1L, xsc)
  end do x1_Loop 
!
  Y = vcopy(YS,L,M2,Tau2,-1)
  CALL sort(M2, Y, YSL) 
  y_Loop: DO k=1,m2
          ysc  = toString( ysL(k) )
	  call addStrings(yL, ysc)
	  call addStrings(xyL, ysc)
   	  call addStrings(yx1L, ysc)
  end do y_Loop
!
   Y1(1:M2) = Y(1:M2)
   Y1((m2+1):(m2+1)) = YS(L+Delta) 
   CALL sort(M2+1, Y1, YSL1)
   y1_Loop: DO k=1, M2+1
            ysc = ToString( YSL1(k) )
 	    call addStrings(y1L, ysc)
            call addStrings(xy1L, ysc)
   end do y1_Loop 
!  
  call Linklist(headx,tailx,xL)   
  call LinkList(heady,taily,yL)
  call Linklist(headx1,tailx1,x1L)   
  call LinkList(heady1,taily1,y1L)
  call Linklist(headxy,tailxy,xyL)
  call LinkList(headxy1,tailxy1,xy1L)
  call Linklist(headyx1,tailyx1,yx1L)   
!
end do time_Loop
!!! --- calculate the Shannon entropy
Hx = SymbolicShannonEntropy(Headx)
call freeList(headx)
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hx1 = SymbolicShannonEntropy(Headx1)
call freeList(headx1)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
Hyx1 = SymbolicShannonEntropy(Headyx1)
call freeList(headyx1)
!
Txy=(hy1-hy) - (hxy1-hxy)
Tyx=(hx1-hx) - (hyx1-hxy)
Hxx=(hy1-hy)
Hyy=(hx1-hx)
!
nullify(headx,tailx,heady,taily,heady1,taily1,headx1,tailx1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
return
!
end subroutine symbolic_TE2_entropies1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE2_entropiesND(ndata,ndim,xs,ys,m1,m2,tau1,tau2, &
                                    Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)           :: ndata,ndim
integer, dimension(ndim)      :: m1,m2,tau1,tau2
real(sifm_real), intent(in)    :: xs(:,:),ys(:,:)
real(sifm_real), intent(out)   :: Txy,Tyx
real(sifm_real), intent(out)   :: Hxx,Hyy,hx,hy
!
type (StateElem), pointer  :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer  :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter         :: Delta = 1
integer                    :: d,istat,L,k,m,tau
character(len=1)           :: sc
character(len=250)         :: xL,x1L,yL,y1L
character(len=250)         :: xyL, xy1L, yx1L
real(sifm_real)             :: hx1,hy1,hxy,hxy1,hyx1
integer, allocatable       :: XSL(:),YSL(:)
integer, allocatable       :: XSL1(:),YSL1(:)
real(sifm_real), allocatable :: X1(:),Y1(:)
real(sifm_real), allocatable :: X(:),Y(:)
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m = max(imax1d(m1),imax1d(m2))
tau = max(imax1d(tau1),imax1d(tau2))
!
!!! --- Loop over all data points
time_Loop: DO L=(M-1)*Tau+1,(Ndata-Delta)
  xL=""; yL=""
  x1L=""; y1L=""
  xyL=""; xy1L=""; yx1L=""
!  
  xdim_Loop: DO d=1,Ndim 
     Allocate( XSL(M1(d)) )
     Allocate( X(M1(d)) )
     X = vcopy(XS(d,:),L,M1(d),Tau1(d),-1)
     CALL sort(M1(d), X, XSL)
     x_Loop: DO k=1,m1(d)
           sc = toString(xsL(k))
	   call addStrings(xL,sc)
	   call addStrings(xyL,sc)
	   call addStrings(xy1L,sc)
     end do x_Loop
     Allocate( XSL1(M1(d)+1) )
     Allocate( X1(M1(d)+1) )
     X1(1:M1(d)) = X(1:M1(d))
     X1((M1(d)+1):(M1(d)+1)) = XS(d,L+Delta) 
     CALL sort(M1(d)+1, X1, XSL1)
     x1_Loop: DO k=1, M1(d)+1
            sc  = ToString( XSL1(k) )
 	    call addStrings(x1L, sc)
            call addStrings(yx1L, sc)
     end do x1_Loop 
     Deallocate(XSL)
     Deallocate(X)
     Deallocate(XSL1)
     Deallocate(X1)
  end do xdim_Loop
!  
  ydim_Loop: DO d=1, Ndim
     Allocate( YSL(m2(d)) )
     Allocate( Y(m2(d)) )
     Y = vcopy(YS(d,:),L,M2(d),Tau2(d),-1)
     CALL sort(M2(d), Y, YSL)
     y_Loop: DO k=1,m2(d)
           sc = toString(ysL(k))
	   call addStrings(yL, sc)
	   call addStrings(xyL, sc)
	   call addStrings(yx1L, sc)
     end do y_Loop
     Allocate( YSL1(m2(d)+1) )
     Allocate( Y1(m2(d)+1) )
     Y1(1:M2(d)) = Y(1:M2(d))
     Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta) 
     CALL sort(M2(d)+1, Y1, YSL1)
     y1_Loop: DO k=1, M2(d)+1
            sc = ToString( YSL1(k) )
 	    call addStrings(y1L, sc)
            call addStrings(xy1L, sc)
     end do y1_Loop 
     Deallocate( YSL )
     Deallocate( Y )
     Deallocate( YSL1 )
     Deallocate( Y1 )
  end do ydim_Loop
! 
  call Linklist(headx,tailx,xL)   
  call LinkList(heady,taily,yL)
  call Linklist(headx1,tailx1,x1L)   
  call LinkList(heady1,taily1,y1L)
  call Linklist(headxy,tailxy,xyL)
  call LinkList(headxy1,tailxy1,xy1L)
  call Linklist(headyx1,tailyx1,yx1L)   
!
end do time_Loop
!!! --- calculate the Shannon entropy
Hx = SymbolicShannonEntropy(Headx)
call freeList(headx)
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hx1 = SymbolicShannonEntropy(Headx1)
call freeList(headx1)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
Hyx1 = SymbolicShannonEntropy(Headyx1)
call freeList(headyx1)
!
Txy=(hy1-hy) - (hxy1-hxy)
Tyx=(hx1-hx) - (hyx1-hxy)
Hxx=(hy1-hy)
Hyy=(hx1-hx)
!
nullify(headx,tailx,heady,taily,heady1,taily1,headx1,tailx1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
return
!
end subroutine symbolic_TE2_entropiesND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE2_entropy1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)           :: ndata,m1,m2,tau1,tau2
real(sifm_real), intent(in)    :: xs(:),ys(:)
real(sifm_real), intent(out) :: Txy
!
type (StateElem), pointer :: heady,taily,heady1,taily1
type (StateElem), pointer :: headxy,tailxy,headxy1,tailxy1
integer, parameter        :: Delta = 1
integer                   :: L,k,m,tau
character(len=1)          :: sc
character(len=250)        :: yL,y1L
character(len=250)        :: xyL, xy1L
real(sifm_real)            :: hy,hy1,hxy,hxy1
integer, dimension(M1)    :: XSL
integer, dimension(M2)    :: YSL
integer, dimension(M2+1)  :: YSL1
real(sifm_real), dimension(M1)    :: X
real(sifm_real), dimension(M2)    :: Y
real(sifm_real), dimension(M2+1)  :: Y1
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(m1,m2)
tau = max(tau1,tau2)
!
Txy = rzero
time_Loop: DO L=(M-1)*Tau+1,(Ndata-Delta)
   yL=""; y1L=""
   xyL=""; xy1L=""
!
   X = vcopy(XS,L,M1,Tau1,-1)
   CALL sort(M1, X, XSL)
   x_Loop: DO k=1, M1
           sc = ToString( XSL(k) )
           call addStrings(xyL, sc)
	   call addStrings(xy1L, sc)
   end do x_Loop
!
   Y = vcopy(YS,L,M2,Tau2,-1)
   CALL sort(M2, Y, YSL)
   y_Loop: DO k=1, M2
           sc = ToString( YSL(k) )
	   call addStrings(yL, sc)
	   call addStrings(xyL, sc)
   end do y_Loop
!
   Y1(1:M2) = Y(1:M2)
   Y1((m2+1):(m2+1)) = YS(L+Delta) 
   CALL sort(M2+1, Y1, YSL1)
   y1_Loop: DO k=1, M2+1
            sc = ToString( YSL1(k) )
 	    call addStrings(y1L, sc)
            call addStrings(xy1L, sc)
   end do y1_Loop 
! 
   call LinkList(heady,taily,yL)  
   call LinkList(heady1,taily1,y1L)
   call Linklist(headxy,tailxy,xyL)
   call LinkList(headxy1,tailxy1,xy1L)
!
end do time_Loop
!
!!! --- calculate the Shannon entropy
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
!
Txy=(hy1-hy) - (hxy1-hxy)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
return
!
end subroutine symbolic_TE2_entropy1D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE2_entropyND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)        :: ndata,ndim
integer, dimension(Ndim), intent(in) :: m1,m2,tau1,tau2
real(sifm_real), dimension(Ndim,Ndata), intent(in) :: xs,ys
real(sifm_real), intent(out) :: Txy
!
type (StateElem), pointer :: heady,taily,heady1,taily1
type (StateElem), pointer :: headxy,tailxy,headxy1,tailxy1
integer, parameter          :: Delta = 1
integer                     :: d,L,k,m,tau
character(len=1)            :: sc
character(len=250)          :: yL,y1L
character(len=250)          :: xyL, xy1L
real(sifm_real)              :: hy,hy1,hxy,hxy1
integer, allocatable        :: XSL(:), YSL(:), YSL1(:)
real(sifm_real), allocatable :: X(:), Y(:), Y1(:)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(imax1d(m1), imax1d(m2))
tau = max(imax1d(tau1), imax1d(tau2))
!
Txy = rzero
time_Loop: DO L=(M-1)*Tau+1,(Ndata-Delta)
   yL=""; y1L=""
   xyL=""; xy1L=""
!   
   xdim_Loop: DO d=1, Ndim
      ALLOCATE( XSL(M1(d)), X(M1(d)) )
      X = vcopy(XS(d,:),L,M1(d),Tau1(d),-1)
      CALL sort(M1(d), X, XSL)
      x_Loop: DO k=1,M1(d)
           sc = toString( XSL(k) )
	   call addStrings(xyL, sc)
	   call addStrings(xy1L, sc)
      end do x_Loop
      DEALLOCATE( XSL,X )
   end do xdim_Loop
!   
   ydim_Loop: DO d=1, Ndim
      ALLOCATE( YSL(M2(d)), Y(M2(d)) )
      Y = vcopy(YS(d,:),L,M2(d),Tau2(d),-1)
      CALL sort(M2(d), Y, YSL)
      y_Loop: DO k=1,m2(d)
           sc = toString( YSL(k) )
	   call addStrings(yL,sc)
           call addStrings(xyL, sc)
      end do y_Loop
!
      Allocate( YSL1(m2(d)+1), Y1(m2(d)+1) )
      Y1(1:M2(d)) = Y(1:M2(d))
      Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta) 
      CALL sort(M2(d)+1, Y1, YSL1)
      y1_Loop: DO k=1, M2(d)+1
            sc = toString( YSL1(k) )
	    call addStrings(y1L,sc)
	    call addStrings(xy1L, sc)
      end do y1_Loop  
      DEALLOCATE( YSL, YSL1 )
      DEALLOCATE( Y, Y1 )
   end do ydim_Loop
!   
   call LinkList(heady,taily,yL)  
   call LinkList(heady1,taily1,y1L)
   call Linklist(headxy,tailxy,xyL)
   call LinkList(headxy1,tailxy1,xy1L)

end do time_Loop
!!! --- calculate the Shannon entropy
Hy = SymbolicShannonEntropy(Heady)
call freeList(heady)
Hy1 = SymbolicShannonEntropy(Heady1)
call freeList(heady1)
Hxy = SymbolicShannonEntropy(Headxy)
call freeList(headxy)
Hxy1 = SymbolicShannonEntropy(Headxy1)
call freeList(headxy1)
!
Txy=(hy1-hy) - (hxy1-hxy)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
return
!
end subroutine symbolic_TE2_entropyND
!
end module TE_class
