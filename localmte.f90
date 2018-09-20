!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: localmte.f90,v 1.0 19-03-2018, IBU
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
module LTE_class
use sifm_kinds
use TEutils_class
use TErandom_class
use TELinkList_class
implicit none
!
      integer, save :: Natoms, Nframes
      integer, save :: Ndim, Ndof
      real(sifm_real), save, allocatable :: XYZ(:,:,:)
      character(len=1), save, allocatable :: XYZS(:,:,:)
      integer, save, allocatable :: mopt(:,:), topt(:,:)
      real(sifm_real), save, dimension(:,:,:), allocatable :: LocalTES,LocalTE
      integer, save :: debug

!{ Thermostat type class }
type nhc_t
      integer                      :: len=3
      real(sifm_real), dimension(3) :: mass,pos,force,vel
end type nhc_t


!
CONTAINS
!
subroutine LTE_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, db, Nshuffles, r0, statP)
  use tempi_class
  implicit none
  integer, intent(in)        :: qTeNorm
  integer, intent(in)        :: qTeShuffle
  integer, intent(in)        :: Na, Nf
  integer, intent(in)        :: Nd, db
  integer, optional          :: Nshuffles
  real(sifm_real), optional  :: r0, statP
!  
  integer :: ierr
  real(sifm_real) :: LocalTime_Start, LocalTime_end
  real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  Nframes = Nf
  Ndim = Nd
  Natoms = Na
  Ndof = Ndim * Natoms
  Debug = db
!
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_LocalTransferEntropy(qTeShuffle,2) 
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
  IF (qTEShuffle > 0) THEN
      CALL getLocalTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  ELSE
      CALL getLocalTransferEntropyNDIM_MPIDOF(qteNorm, debug)
  ENDIF
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for Local TE:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  call WRITE_LTE(qTeShuffle)
!
  CALL DEALLOCATE_LocalTransferEntropy()
!
  CALL MPI_finish()
!
  STOP 'Program finished successfully'    
!
  return
end subroutine LTE_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LSTE_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, db, Nshuffles, r0, statP)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: qTeNorm
  integer, intent(in)        :: qTeShuffle
  integer, intent(in)        :: Na, Nf
  integer, intent(in)        :: Nd, db
  integer, optional          :: Nshuffles
  real(sifm_real), optional  :: r0, statP
!  
  integer :: ierr
  real(sifm_real) :: LocalTime_Start, LocalTime_end
  real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  Nframes = Nf
  Ndim = Nd
  Natoms = Na
  Ndof = Ndim * Natoms
  Debug = db
!
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_LocalTransferEntropy(qTeShuffle,1) 
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
  IF (qTEShuffle > 0) THEN
      CALL getLocalSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  ELSE
      CALL getLocalSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
  ENDIF
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for Local TE:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  call WRITE_LTE(qTeShuffle)
!
  CALL DEALLOCATE_LocalTransferEntropy()
!
  CALL MPI_finish()
!
  STOP 'Program finished successfully'    
!
  return
end subroutine LSTE_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WRITE_LTE(qTEShuffle)
!
Integer, intent(in) :: qTESHUFFLE
!
integer             :: t, i, j, nunit
! 
nunit = new_unit()
open(unit=nunit, file='localte.txt', status='unknown', action='write')
DO t=1, Nframes
   DO i=1, Natoms
      write(nunit,*) (LocalTE(t,i,j), j=1, natoms)
   ENDDO
ENDDO
!
IF (qteshuffle > 0) THEN
nunit = new_unit()
open(unit=nunit, file='localste.txt', status='unknown', action='write')
DO t=1, Nframes
   DO i=1, Natoms
      write(nunit,*) (LocalTES(t,i,j), j=1, natoms)
   ENDDO
ENDDO
!
ENDIF
!
return
!
end subroutine WRITE_LTE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Allocate_LocalTransferEntropy(qTeShuffle, model)
      implicit none
      integer, intent(in) :: qTeShuffle, model
!
      integer :: ierr
!
      IF (.not. Allocated(LocalTE)) THEN
	  ALLOCATE( LocalTE(Nframes, Natoms, Natoms), stat=ierr )
          if (ierr /= 0) stop 'Error allocating Local TE'
      ENDIF
!
      IF (qTEShuffle > 0) THEN
          IF (.not. Allocated(LocalTES)) THEN
	      ALLOCATE( LocalTES(Nframes, Natoms, Natoms), stat=ierr )
              if (ierr /= 0) stop 'Error allocating Local TES'
          ENDIF
      ENDIF
!
      IF (model == 1) THEN
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
!
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
end subroutine Allocate_LocalTransferEntropy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine deAllocate_LocalTransferEntropy()
      implicit none
!
      integer :: ierr
!
      IF (Allocated(LocalTE)) THEN
	  DEALLOCATE( LocalTE, stat=ierr )
          if (ierr /= 0) stop 'Error deallocating Local TE'
      ENDIF
!
      IF (Allocated(LocalTES)) THEN
	  DEALLOCATE( LocalTES, stat=ierr )
          if (ierr /= 0) stop 'Error deallocating Local TES'
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
end subroutine deAllocate_LocalTransferEntropy
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
call MPI_BCAST (debug,      1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (LocalTE,    Natoms*Natoms*Nframes, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (qTeNorm,    1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (qTEShuffle, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
IF (qTEShuffle > 0) THEN
    call MPI_BCAST (Nshuffles,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (r0,         1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (statP,      1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (LocalTES,   Natoms*Natoms*Nframes, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
ENDIF
!
if (model == 1) then
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLocalSTransferEntropyNDIM_MPIDOF(qteNorm, debug)
USE TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: debug
!
!{Local variables}
integer, allocatable :: listPairs(:,:)
integer :: Npairs, ilist
integer :: i, j, t, ierr, K
real(sifm_real), dimension(Nframes) :: Ltxy,Ltyx
real(sifm_real), dimension(:,:), allocatable :: LocalTxy,LocalTyx
real(sifm_real), dimension(:,:), allocatable :: GlobalTxy,GlobalTyx
real(sifm_real) :: time_start, time_end
integer :: iStart, iEnd, Stride
!
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
IF ( mod(Npairs, NumProcs) == 0 ) THEN
     Stride = Npairs / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd
ELSE
     write(*,*) Npairs, NumProcs
     Stop "Number of trials must be Multiple of Number of Processors"
ENDIF
!
allocate(LocalTxy(Stride,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-Txy"
allocate(LocalTyx(Stride,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating Local-Tyx"
allocate(GlobalTxy(Npairs,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-Txy"
allocate(GlobalTyx(Npairs,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-Tyx"
!
! --- Compute the entropy transfer
LocalTE  = rzero
LocalTxy = rzero
LocalTyx = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = iStart, iEnd
      K = iList - iStart + 1
!
      IF (myID == Master) call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
!
      !--- Compute the direct transfer entropy
      CALL symbolic_LTE_entropies(nframes, ndim, XYZs(i,:,:), XYZs(j,:,:),               &
                                  Mopt(i,:), Mopt(j,:), Topt(i,:), Topt(j,:),            &
                                  LTxy, LTyx)
!
      LocalTxy(K,:) = LTxy
      LocalTyx(K,:) = LTyx
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Local Txy and Tyx:   ", F10.6)') Time_End - Time_Start
      ENDIF
!
enddo Pairs_LOOP
!
CALL MPI_gather( LocalTxy,  Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 GlobalTxy, Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
CALL MPI_gather( LocalTyx,  Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 GlobalTyx, Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )

IF (myID == master) THEN
    DO iList = 1, Npairs
       I = ListPairs(iList, 1)
       J = ListPairs(iList, 2)
       LocalTE(:,i,j) = GlobalTxy(iList,:)
       LocalTE(:,j,i) = GlobalTyx(iList,:)
    ENDDO
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
deallocate(LocalTxy, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-Txy"
deallocate(LocalTyx, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating Local-Tyx"
deallocate(GlobalTxy, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-Txy"
deallocate(GlobalTyx, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-Tyx"
!
return
end subroutine getLocalSTransferEntropyNDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLocalSTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
USE TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!
!{Local variables}
integer, allocatable :: listPairs(:,:)
integer :: Npairs, ilist
integer :: ntrials
integer :: i, j, itrial, t, ierr, K
character(len=1), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real), dimension(Nframes) :: Ltxy,Ltyx
real(sifm_real), dimension(:,:), allocatable :: LocalTxy,LocalTyx
real(sifm_real), dimension(:,:), allocatable :: GlobalTxy,GlobalTyx
real(sifm_real) :: time_start, time_end
integer :: iStart, iEnd, Stride
!
! -- Initialize 2*natoms random number sequences
IF ( qteshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, 0)
     ntrials = nshuffles
endif
!
IF ( mod(Ntrials, NumProcs) == 0 ) THEN
     Stride = Ntrials / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd
ELSE
     write(*,*) Ntrials, NumProcs
     Stop "Number of trials must be Multiple of Number of Processors"
ENDIF

IF ( qteshuffle > 0 ) then
     allocate(LocalTxy(Stride,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating local-Txy"
     allocate(LocalTyx(Stride,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating Local-Tyx"
     allocate(GlobalTxy(Ntrials,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating global-Txy"
     allocate(GlobalTyx(Ntrials,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating global-Tyx"
ENDIF
!

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
! --- Compute the entropy transfer
localTE = RZERO
IF ( qteshuffle > 0 ) THEN
     localTES = RZERO
ENDIF
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      IF (myID == Master) call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
!
      !--- Compute the direct transfer entropy
      IF (myID == Master) THEN
          CALL symbolic_LTE_entropies(nframes, ndim, XYZs(i,:,:), XYZs(j,:,:),                   &
                                      Mopt(i,:), Mopt(j,:), Topt(i,:), Topt(j,:),            &
                                      LTxy, LTyx)
          LocalTE(:,i,j) = LTxy
          LocalTE(:,j,i) = LTyx
      ENDIF
      
      IF (qteshuffle > 0) THEN
           ! -- Subtract the bias from T(x->y) and normalize it
           LocalTxy = rzero
           xs_shuffle = XYZS(i,:,:)
           DO iTrial = iStart, iEnd
               K = iTrial - iStart + 1
               IF (qTEShuffle == 1) THEN
                   CALL shuffleND(nframes,ndim,xs_shuffle,I)
               ELSE
                   CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(i,:), Topt(i,:), I)
               ENDIF
               CALL symbolic_LTE_entropy(nframes,ndim,xs_shuffle,XYZs(j,:,:),                &
                                         Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),            &
                                         LocalTxy(K,:))
           ENDDO
           CALL MPI_gather( LocalTxy,  Nframes*Stride, MPI_INTEGER, &
                            GlobalTxy, Nframes*Stride, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierr )
           IF (myID == Master) THEN
               LTxy = sum( GlobalTxy, dim=1 ) / real(Ntrials, sifm_real)
               LocalTES(:,i,j) = LocalTE(:,i,j) - LTxy
           ENDIF
 
      ! -- Subtract the bias from T(y->x) and normalize it
           LocalTyx = rzero
           xs_shuffle = XYZs(j,:,:)
           DO iTrial = iStart, iEnd
              K = iTrial - iStart + 1
              IF (qTEShuffle == 1) THEN
                  CALL shuffleND(nframes,ndim,xs_shuffle,Natoms+J)
              ELSE
                  CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
              ENDIF
              CALL symbolic_LTE_entropy(nframes,ndim,xs_shuffle, XYZs(i,:,:),            &
                                        Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),         &
                                        LocalTyx(K,:))
           ENDDO
           CALL MPI_gather( LocalTyx,  Nframes*Stride, MPI_INTEGER, &
                            GlobalTyx, Nframes*Stride, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierr )
           IF (myID == Master) THEN
               LTyx = sum( GlobalTyx, dim=1 ) / real(Ntrials, sifm_real)
               LocalTES(:,j,i) = LocalTE(:,j,i) - LTyx
           ENDIF
      ENDIF 
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Local Txy and Tyx:   ", F10.6)') Time_End - Time_Start
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

IF ( qteshuffle > 0 ) then
     deallocate(LocalTxy, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating local-Txy"
     deallocate(LocalTyx, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating Local-Tyx"
     deallocate(GlobalTxy, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating global-Txy"
     deallocate(GlobalTyx, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating global-Tyx"
ENDIF
!
return
end subroutine getLocalSTransferEntropyNDIM_MPISFL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLocalTransferEntropyNDIM_MPIDOF(qteNorm, debug)
USE TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: debug
!
!{Local variables}
integer, allocatable :: listPairs(:,:)
integer :: Npairs, ilist
integer :: i, j, t, ierr, K
real(sifm_real), dimension(Nframes) :: Ltxy,Ltyx
real(sifm_real), dimension(:,:), allocatable :: LocalTxy,LocalTyx
real(sifm_real), dimension(:,:), allocatable :: GlobalTxy,GlobalTyx
real(sifm_real) :: time_start, time_end
integer :: iStart, iEnd, Stride
!
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
IF ( mod(Npairs, NumProcs) == 0 ) THEN
     Stride = Npairs / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd
ELSE
     write(*,*) Npairs, NumProcs
     Stop "Number of trials must be Multiple of Number of Processors"
ENDIF

allocate(LocalTxy(Stride,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-Txy"
allocate(LocalTyx(Stride,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating Local-Tyx"
allocate(GlobalTxy(Npairs,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-Txy"
allocate(GlobalTyx(Npairs,Nframes), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-Tyx"
!
! --- Compute the entropy transfer
LocalTE = RZERO
LocalTxy = rzero
LocalTyx = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = iStart, iEnd
      K = iList - iStart + 1
!
      IF (myID == Master) call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
!
      !--- Compute the direct transfer entropy
      CALL symbolic_LTE2_entropies(nframes, ndim, XYZ(i,:,:), XYZ(j,:,:),                    &
                                      Mopt(i,:), Mopt(j,:), Topt(i,:), Topt(j,:),            &
                                      LTxy, LTyx)
!
      LocalTxy(K,:) = LTxy
      LocalTyx(K,:) = LTyx
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Local Txy and Tyx:   ", F10.6)') Time_End - Time_Start
      ENDIF
!
enddo Pairs_LOOP
!
CALL MPI_gather( LocalTxy,  Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 GlobalTxy, Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
CALL MPI_gather( LocalTyx,  Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 GlobalTyx, Stride*Nframes, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )

IF (myID == master) THEN
    DO iList = 1, Npairs
       I = ListPairs(iList, 1)
       J = ListPairs(iList, 2)
       LocalTE(:,i,j) = GlobalTxy(iList,:)
       LocalTE(:,j,i) = GlobalTyx(iList,:)
    ENDDO
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
deallocate(LocalTxy, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-Txy"
deallocate(LocalTyx, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating Local-Tyx"
deallocate(GlobalTxy, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-Txy"
deallocate(GlobalTyx, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-Tyx"
!
return
end subroutine getLocalTransferEntropyNDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLocalTransferEntropyNDIM_MPISFL(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
USE TEMPI_CLASS
implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!
!{Local variables}
integer, allocatable :: listPairs(:,:)
integer :: Npairs, ilist
integer :: ntrials
integer :: i, j, itrial, t, ierr, K
real(sifm_real), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real), dimension(Nframes) :: Ltxy,Ltyx
real(sifm_real), dimension(:,:), allocatable :: LocalTxy,LocalTyx
real(sifm_real), dimension(:,:), allocatable :: GlobalTxy,GlobalTyx
real(sifm_real) :: time_start, time_end
integer :: iStart, iEnd, Stride
!
! -- Initialize 2*natoms random number sequences
IF ( qteshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, 0)
     ntrials = nshuffles
endif
!
IF ( mod(Ntrials, NumProcs) == 0 ) THEN
     Stride = Ntrials / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd
ELSE
     write(*,*) Ntrials, NumProcs
     Stop "Number of DoF must be Multiple of Number of Processors"
ENDIF

IF ( qteshuffle > 0 ) then
     allocate(LocalTxy(Stride,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating local-Txy"
     allocate(LocalTyx(Stride,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating Local-Tyx"
     allocate(GlobalTxy(Ntrials,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating global-Txy"
     allocate(GlobalTyx(Ntrials,Nframes), STAT=ierr )
     if (ierr /= 0) STOP "Error allocating global-Tyx"
ENDIF
!

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
! --- Compute the entropy transfer
localTE = rzero
IF ( qteshuffle > 0 ) THEN
     localTES = rzero
ENDIF
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      IF (myID == Master) call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)
!
      !--- Compute the direct transfer entropy
      IF (myID == Master) THEN
          CALL symbolic_LTE2_entropies(nframes, ndim, XYZ(i,:,:), XYZ(j,:,:),                 &
                                       Mopt(i,:), Mopt(j,:), Topt(i,:), Topt(j,:),            &
                                       LTxy, LTyx)
          LocalTE(:,i,j) = LTxy
          LocalTE(:,j,i) = LTyx
      ENDIF
      
      IF (qteshuffle > 0) THEN
           ! -- Subtract the bias from T(x->y) and normalize it
           LocalTxy = rzero
           xs_shuffle = XYZ(i,:,:)
           DO iTrial = iStart, iEnd
               K = iTrial - iStart + 1
               IF (qTEShuffle == 1) THEN
                   CALL shuffleND2(nframes,ndim,xs_shuffle,I)
               ELSE
                   CALL BlockShuffleND2(nframes, ndim, xs_shuffle, Mopt(I,:), Topt(I,:), I)
               ENDIF
               CALL symbolic_LTE2_entropy(nframes,ndim,xs_shuffle,XYZ(j,:,:),                &
                                          Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),            &
                                          LocalTxy(K,:))
           ENDDO
           CALL MPI_gather( LocalTxy,  Nframes*Stride, MPI_INTEGER, &
                            GlobalTxy, Nframes*Stride, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierr )
           IF (myID == Master) THEN
               LTxy = sum( GlobalTxy, dim=1 ) / real(Ntrials, sifm_real)
               LocalTES(:,i,j) = LocalTE(:,i,j) - LTxy
           ENDIF

      ! -- Subtract the bias from T(y->x) and normalize it
           LocalTyx = rzero
           xs_shuffle = XYZ(j,:,:)
           DO iTrial = iStart, iEnd
              K = iTrial - iStart + 1
              IF (qTEShuffle == 1) THEN
                  CALL shuffleND2(nframes,ndim,xs_shuffle,Natoms+J)
              ELSE
                  CALL BlockShuffleND2(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
              ENDIF
              CALL symbolic_LTE2_entropy(nframes,ndim,xs_shuffle, XYZ(i,:,:),            &
                                        Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),         &
                                        LocalTyx(K,:))
           ENDDO
           CALL MPI_gather( LocalTyx,  Nframes*Stride, MPI_INTEGER, &
                            GlobalTyx, Nframes*Stride, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierr )
           IF (myID == Master) THEN
               LTyx = sum( GlobalTyx, dim=1 ) / real(Ntrials, sifm_real)
               LocalTES(:,j,i) = LocalTE(:,j,i) - LTyx
           ENDIF
      ENDIF 
!
      IF (myID == Master) THEN
          CALL CPU_time( Time_End )
          write(*,'("CPU Elapsed time for computation of Local Txy and Tyx:   ", F10.6)') Time_End - Time_Start
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

IF ( qteshuffle > 0 ) then
     deallocate(LocalTxy, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating local-Txy"
     deallocate(LocalTyx, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating Local-Tyx"
     deallocate(GlobalTxy, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating global-Txy"
     deallocate(GlobalTyx, STAT=ierr )
     if (ierr /= 0) STOP "Error deallocating global-Tyx"
ENDIF
!
return
end subroutine getLocalTransferEntropyNDIM_MPISFL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_LTE2_entropies(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy,LTyx)
implicit none
!
integer, intent(in)                               :: ndata,ndim
integer, dimension(ndim), intent(in)              :: m1,m2,tau1,tau2
real(sifm_real), dimension(ndim,ndata), intent(in) :: xs,ys
real(sifm_real), dimension(ndata), intent(out)     :: LTxy,LTyx
!
type (StateElem), pointer                 :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer                 :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter                        :: Delta=1
integer                                   :: ncounts_x,ncounts_y,ncounts_x1,ncounts_y1, &
                                             ncounts_xy,ncounts_xy1,ncounts_yx1
integer                                   :: nstates_x,nstates_y,nstates_x1,nstates_y1, &
                                             nstates_xy,nstates_xy1,nstates_yx1
integer                                   :: d,L,k,m,tau,t,ndims
character(len=1)                          :: xsc,ysc
character(len=250)                        :: xL,x1L,yL,y1L
character(len=250)                        :: xyL, xy1L, yx1L
real(sifm_real), dimension(:), allocatable :: hx1,hy1,hxy,hyx,hxy1,hyx1,hx,hy
integer, allocatable                      :: XSL(:),YSL(:)
integer, allocatable                      :: XSL1(:),YSL1(:)
real(sifm_real), allocatable               :: X1(:),Y1(:)
real(sifm_real), allocatable               :: X(:),Y(:)
real(sifm_real), dimension(:), allocatable :: prob_x, prob_x1, prob_y1, Prob_y
real(sifm_real), dimension(:), allocatable :: Prob_xy, Prob_xy1, prob_yx1
integer, dimension(:), allocatable        :: frequency_x, frequency_x1, frequency_y1, frequency_y
integer, dimension(:), allocatable        :: frequency_xy, frequency_xy1, frequency_yx1
integer, dimension(:), allocatable        :: states_y1,states_y,states_x,states_x1, &
                                             states_xy,states_xy1,states_yx1
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m = max(imax1d(m1),imax1d(m2))
tau = max(imax1d(tau1),imax1d(tau2))
!
ndims = Ndata-Delta-(M-1)*Tau
call mteallocate_r(hx,   ndims,  "hx",   "symbolic-LTE-entropies")
call mteallocate_r(hx1,  ndims,  "hx1",  "symbolic-LTE-entropies")
call mteallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropies")
call mteallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropies")
call mteallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropies")
call mteallocate_r(hyx,  ndims,  "hyx",  "symbolic-LTE-entropies")
call mteallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropies")
call mteallocate_r(hyx1, ndims,  "hyx1", "symbolic-LTE-entropies")
!
call mteallocate_i(states_x,   ndims, "states-x",   "symbolic-LTE-entropies")
call mteallocate_i(states_x1,  ndims, "states-x1",  "symbolic-LTE-entropies")
call mteallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropies")
call mteallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropies")
call mteallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropies")
call mteallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropies")
call mteallocate_i(states_yx1, ndims, "states-yx1", "symbolic-LTE-entropies")

!!! --- Loop over all data points
t = 0
time_Loop1: DO L=((M-1)*Tau+1),(Ndata-Delta)
!
  t = t + 1
  xL=""; yL=""
  x1L=""; y1L=""
  xyL=""; xy1L=""; yx1L=""
!
  xdim_loop: DO d=1, ndim
!
  Allocate( XSL(M1(d)) )
  Allocate( X(M1(d)) )
  X = vcopy(XS(d,:),L,M1(d),Tau1(d),-1)
  CALL sort(M1(d), X, XSL)
!
  x_Loop1: DO k=1,m1(d)
          xsc  = toString( xsL(k) )
          call addStrings(xL,xsc)
          call addStrings(xyL,xsc)
          call addStrings(xy1L,xsc)
  end do x_Loop1
!
  Allocate( XSL1(M1(d)+1) )
  Allocate( X1(M1(d)+1) )
  X1(1:M1(d)) = X(1:M1(d))
  X1((M1(d)+1):(M1(d)+1)) = XS(d,L+Delta)
  CALL sort(M1(d)+1, X1, XSL1)
  x1_Loop: DO k=1, M1(d)+1
           xsc  = ToString( XSL1(k) )
           call addStrings(x1L, xsc)
           call addStrings(yx1L, xsc)
  end do x1_Loop
  Deallocate(XSL)
  Deallocate(X)
  Deallocate(XSL1)
  Deallocate(X1)
  end do xdim_loop
!
  ydim_loop: do d=1, ndim
  Allocate( YSL(m2(d)) )
  Allocate( Y(m2(d)) )
  Y = vcopy(YS(d,:),L,M2(d),Tau2(d),-1)
  CALL sort(M2(d), Y, YSL)
  y_Loop1: DO k=1,m2(d)
         ysc  = toString( ysl(k) )
         call addStrings(yL,ysc)
         call addStrings(xyL,ysc)
         call addStrings(yx1L,ysc)
  end do y_Loop1
!
  Allocate( YSL1(m2(d)+1) )
  Allocate( Y1(m2(d)+1) )
  Y1(1:M2(d)) = Y(1:M2(d))
  Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta)
  CALL sort(M2(d)+1, Y1, YSL1)
  y1_Loop: DO k=1, M2(d)+1
           ysc = ToString( YSL1(k) )
           call addStrings(y1L, ysc)
           call addStrings(xy1L, ysc)
  end do y1_Loop
  Deallocate( YSL )
  Deallocate( Y )
  Deallocate( YSL1 )
  Deallocate( Y1 )
  end do ydim_loop
!
  call LocalLinkList(headx,tailx,xL,states_x(t))  
  call LocalLinkList(headx1,tailx1,x1L,states_x1(t))
  call LocalLinkList(heady,taily,yL,states_y(t))  
  call LocalLinkList(heady1,taily1,y1L,states_y1(t))
  call LocalLinklist(headxy,tailxy,xyL,states_xy(t))
  call LocalLinkList(headxy1,tailxy1,xy1L,states_xy1(t))  
  call LocalLinkList(headyx1,tailyx1,yx1L,states_yx1(t))  
!
end do time_Loop1
!
ncounts_x   = getTotalStates(Headx)
ncounts_y   = getTotalStates(Heady)
ncounts_x1  = getTotalStates(Headx1)
ncounts_y1  = getTotalStates(Heady1)
ncounts_xy  = getTotalStates(Headxy)
ncounts_xy1 = getTotalStates(Headxy1)
ncounts_yx1 = getTotalStates(Headyx1)
!
nstates_x = CountStates(headx)
nstates_y = CountStates(heady)
nstates_x1 = CountStates(headx1)
nstates_y1 = CountStates(heady1)
nstates_xy = CountStates(headxy)
nstates_xy1 = CountStates(headxy1)
nstates_yx1 = CountStates(headyx1)
!
call mteallocate_r(prob_x,   ncounts_x,   "prob-x",   "symbolic-LTE-entropies")
call mteallocate_r(prob_y,   ncounts_y,   "prob-y",   "symbolic-LTE-entropies")
call mteallocate_r(prob_x1,  ncounts_x1,  "prob-x1",  "symbolic-LTE-entropies")
call mteallocate_r(prob_y1,  ncounts_y1,  "prob-y1",  "symbolic-LTE-entropies")
call mteallocate_r(prob_xy,  ncounts_xy,  "prob-xy",  "symbolic-LTE-entropies")
call mteallocate_r(prob_xy1, ncounts_xy1, "prob-xy1", "symbolic-LTE-entropies")
call mteallocate_r(prob_yx1, ncounts_yx1, "prob-yx1", "symbolic-LTE-entropies")
!
call mteallocate_i(frequency_x,   ncounts_x,   "frequency-x",   "symbolic-LTE-entropies")
call mteallocate_i(frequency_y,   ncounts_y,   "frequency-y",   "symbolic-LTE-entropies")
call mteallocate_i(frequency_x1,  ncounts_x1,  "frequency-x1",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_y1,  ncounts_y1,  "frequency-y1",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_xy,  ncounts_xy,  "frequency-xy",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_xy1, ncounts_xy1, "frequency-xy1", "symbolic-LTE-entropies")
call mteallocate_i(frequency_yx1, ncounts_yx1, "frequency-yx1", "symbolic-LTE-entropies")
!
call SymbolicJointProbabilityTraj(Headx, ncounts_x, Prob_x, frequency_x)
call SymbolicJointProbabilityTraj(Headx1, ncounts_x1, Prob_x1, frequency_x1)
call SymbolicJointProbabilityTraj(Heady, ncounts_y, Prob_y, frequency_y)
call SymbolicJointProbabilityTraj(Heady1, ncounts_y1, Prob_y1, frequency_y1)
call SymbolicJointProbabilityTraj(Headxy, ncounts_xy, Prob_xy, frequency_xy)
call SymbolicJointProbabilityTraj(Headxy1, ncounts_xy1, Prob_xy1, frequency_xy1)
call SymbolicJointProbabilityTraj(Headyx1, ncounts_yx1, Prob_yx1, frequency_yx1)
!
t = 0
time_Loop2: DO L=((M-1)*Tau+1),(Ndata-Delta)
  t = t + 1
  Hx(t)   = -log( Prob_x(states_x(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)
  Hx1(t)  = -log( Prob_x1(states_x1(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)
  Hyx(t)  = -log( Prob_xy(states_xy(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)
  Hyx1(t) = -log( Prob_yx1(states_yx1(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)

  Hy(t)   = -log( Prob_y(states_y(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  Hy1(t)  = -log( Prob_y1(states_y1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  Hxy(t)  = -log( Prob_xy(states_xy(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  Hxy1(t) = -log( Prob_xy1(states_xy1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
end do time_Loop2
!
do t = 1, nframes
   if (t <= ndims) then
   LTxy(t) = (hy1(t)-hy(t)) - (hxy1(t)-hxy(t))
   LTyx(t) = (hx1(t)-hx(t)) - (hyx1(t)-hyx(t))
   else
   LTxy(t) = rzero
   LTyx(t) = rzero
   endif
enddo
!
call freeList(headx)
call freeList(heady)
call freeList(headx1)
call freeList(heady1)
call freeList(headxy)
call freeList(headxy1)
call freeList(headyx1)
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
call mtedeallocate_r(hx,   ndims,  "hx",   "symbolic-LTE-entropies")
call mtedeallocate_r(hx1,  ndims,  "hx1",  "symbolic-LTE-entropies")
call mtedeallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropies")
call mtedeallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropies")
call mtedeallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropies")
call mtedeallocate_r(hyx,  ndims,  "hyx",  "symbolic-LTE-entropies")
call mtedeallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropies")
call mtedeallocate_r(hyx1, ndims,  "hyx1", "symbolic-LTE-entropies")
!
call mtedeallocate_i(states_x,   ndims, "states-x",   "symbolic-LTE-entropies")
call mtedeallocate_i(states_x1,  ndims, "states-x1",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropies")
call mtedeallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropies")
call mtedeallocate_i(states_yx1, ndims, "states-yx1", "symbolic-LTE-entropies")
!
call mtedeallocate_r(prob_x,   ncounts_x,   "prob-x",   "symbolic-LTE-entropies")
call mtedeallocate_r(prob_y,   ncounts_y,   "prob-y",   "symbolic-LTE-entropies")
call mtedeallocate_r(prob_x1,  ncounts_x1,  "prob-x1",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_y1,  ncounts_y1,  "prob-y1",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_xy,  ncounts_xy,  "prob-xy",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_xy1, ncounts_xy1, "prob-xy1", "symbolic-LTE-entropies")
call mtedeallocate_r(prob_yx1, ncounts_yx1, "prob-yx1", "symbolic-LTE-entropies")
!
call mtedeallocate_i(frequency_x,   ncounts_x,   "frequency-x",   "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_y,   ncounts_y,   "frequency-y",   "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_x1,  ncounts_x1,  "frequency-x1",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_y1,  ncounts_y1,  "frequency-y1",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_xy,  ncounts_xy,  "frequency-xy",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_xy1, ncounts_xy1, "frequency-xy1", "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_yx1, ncounts_yx1, "frequency-yx1", "symbolic-LTE-entropies")
!
return
end subroutine symbolic_LTE2_entropies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_LTE2_entropy(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy)
implicit none
!
integer, intent(in)                               :: ndata,ndim
integer, dimension(ndim), intent(in)              :: m1,m2,tau1,tau2
real(sifm_real), dimension(ndim,ndata), intent(in) :: xs,ys
real(sifm_real), dimension(ndata), intent(out)     :: LTxy
!
type (StateElem), pointer                  :: heady,taily,heady1,taily1
type (StateElem), pointer                  :: headxy,tailxy,headxy1,tailxy1
integer, parameter                         :: Delta=1
integer                                    :: ncounts_y,ncounts_y1, &
                                              ncounts_xy,ncounts_xy1
integer                                    :: nstates_y,nstates_y1, &
                                              nstates_xy,nstates_xy1
integer                                    :: d,L,k,m,tau,t,ndims
character(len=1)                           :: xsc,ysc
character(len=250)                         :: xL, yL,y1L
character(len=250)                         :: xyL, xy1L
real(sifm_real), dimension(:), allocatable :: hy1,hxy,hxy1,hy
integer, allocatable                       :: XSL(:),YSL(:)
integer, allocatable                       :: YSL1(:)
real(sifm_real), allocatable               :: Y1(:)
real(sifm_real), allocatable               :: X(:),Y(:)
real(sifm_real), dimension(:), allocatable :: prob_y1, Prob_y
real(sifm_real), dimension(:), allocatable :: Prob_xy, Prob_xy1
integer, dimension(:), allocatable        :: frequency_y1, frequency_y
integer, dimension(:), allocatable        :: frequency_xy, frequency_xy1
integer, dimension(:), allocatable        :: states_y1, states_y, states_xy, states_xy1
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(imax1d(m1),imax1d(m2))
tau = max(imax1d(tau1),imax1d(tau2))
!
ndims = Ndata-Delta-(M-1)*Tau
call mteallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropies")
call mteallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropies")
call mteallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropies")
call mteallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropies")
!
call mteallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropies")
call mteallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropies")
call mteallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropies")
call mteallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropies")

!!! --- Loop over all data points
t = 0
time_Loop1: DO L=((M-1)*Tau+1),(Ndata-Delta)
!
  t = t + 1
  yL=""; y1L=""
  xyL=""; xy1L=""
!
  xdim_loop: DO d=1, ndim
!
  Allocate( XSL(M1(d)) )
  Allocate( X(M1(d)) )
  X = vcopy(XS(d,:),L,M1(d),Tau1(d),-1)
  CALL sort(M1(d), X, XSL)
!
  x_Loop1: DO k=1,m1(d)
          xsc  = toString( xsL(k) )
          call addStrings(xyL,xsc)
          call addStrings(xy1L,xsc)
  end do x_Loop1
  Deallocate(XSL)
  Deallocate(X)
  end do xdim_loop
!
  ydim_loop: do d=1, ndim
  Allocate( YSL(m2(d)) )
  Allocate( Y(m2(d)) )
  Y = vcopy(YS(d,:),L,M2(d),Tau2(d),-1)
  CALL sort(M2(d), Y, YSL)
  y_Loop1: DO k=1,m2(d)
         ysc  = toString( ysl(k) )
         call addStrings(yL,ysc)
         call addStrings(xyL,ysc)
  end do y_Loop1
!
  Allocate( YSL1(m2(d)+1) )
  Allocate( Y1(m2(d)+1) )
  Y1(1:M2(d)) = Y(1:M2(d))
  Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta)
  CALL sort(M2(d)+1, Y1, YSL1)
  y1_Loop: DO k=1, M2(d)+1
           ysc = ToString( YSL1(k) )
           call addStrings(y1L, ysc)
           call addStrings(xy1L, ysc)
  end do y1_Loop
  Deallocate( YSL )
  Deallocate( Y )
  Deallocate( YSL1 )
  Deallocate( Y1 )
  end do ydim_loop
!
  call LocalLinkList(heady,taily,yL,states_y(t))  
  call LocalLinkList(heady1,taily1,y1L,states_y1(t))
  call LocalLinklist(headxy,tailxy,xyL,states_xy(t))
  call LocalLinkList(headxy1,tailxy1,xy1L,states_xy1(t))   
!
end do time_Loop1
!
ncounts_y   = getTotalStates(Heady)
ncounts_y1  = getTotalStates(Heady1)
ncounts_xy  = getTotalStates(Headxy)
ncounts_xy1 = getTotalStates(Headxy1)
!
nstates_y = CountStates(heady)
nstates_y1 = CountStates(heady1)
nstates_xy = CountStates(headxy)
nstates_xy1 = CountStates(headxy1)

!
call mteallocate_r(prob_y,   ncounts_y,   "prob-y",   "symbolic-LTE-entropies")
call mteallocate_r(prob_y1,  ncounts_y1,  "prob-y1",  "symbolic-LTE-entropies")
call mteallocate_r(prob_xy,  ncounts_xy,  "prob-xy",  "symbolic-LTE-entropies")
call mteallocate_r(prob_xy1, ncounts_xy1, "prob-xy1", "symbolic-LTE-entropies")
!
call mteallocate_i(frequency_y,   ncounts_y,   "frequency-y",   "symbolic-LTE-entropies")
call mteallocate_i(frequency_y1,  ncounts_y1,  "frequency-y1",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_xy,  ncounts_xy,  "frequency-xy",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_xy1, ncounts_xy1, "frequency-xy1", "symbolic-LTE-entropies")
!
call SymbolicJointProbabilityTraj(Heady, ncounts_y, Prob_y, frequency_y)
call SymbolicJointProbabilityTraj(Heady1, ncounts_y1, Prob_y1, frequency_y1)
call SymbolicJointProbabilityTraj(Headxy, ncounts_xy, Prob_xy, frequency_xy)
call SymbolicJointProbabilityTraj(Headxy1, ncounts_xy1, Prob_xy1, frequency_xy1)
!
t = 0
time_Loop2: DO L=((M-1)*Tau+1),(Ndata-Delta)
  t = t + 1
  Hy(t)   = -log( Prob_y(states_y(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  Hy1(t)  = -log( Prob_y1(states_y1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  Hxy(t)  = -log( Prob_xy(states_xy(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  Hxy1(t) = -log( Prob_xy1(states_xy1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
end do time_Loop2
!
do t = 1, nframes
   if (t <= ndims) then
   LTxy(t) = (hy1(t)-hy(t)) - (hxy1(t)-hxy(t))
   else
   LTxy(t) = rzero
   endif
enddo
!
call freeList(heady)
call freeList(heady1)
call freeList(headxy)
call freeList(headxy1)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
call mtedeallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropies")
call mtedeallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropies")
call mtedeallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropies")
call mtedeallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropies")
!
call mtedeallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropies")
call mtedeallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropies")
!
call mtedeallocate_r(prob_y,   ncounts_y,   "prob-y",   "symbolic-LTE-entropies")
call mtedeallocate_r(prob_y1,  ncounts_y1,  "prob-y1",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_xy,  ncounts_xy,  "prob-xy",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_xy1, ncounts_xy1, "prob-xy1", "symbolic-LTE-entropies")
!
call mtedeallocate_i(frequency_y,   ncounts_y,   "frequency-y",   "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_y1,  ncounts_y1,  "frequency-y1",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_xy,  ncounts_xy,  "frequency-xy",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_xy1, ncounts_xy1, "frequency-xy1", "symbolic-LTE-entropies")
!
return
end subroutine symbolic_LTE2_entropy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_LTE_entropies(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy,LTyx)
implicit none
!
integer, intent(in)                           :: ndata,ndim
integer, dimension(ndim), intent(in)          :: m1,m2,tau1,tau2
character(len=1), dimension(ndim,ndata), intent(in)    :: xs,ys
real(sifm_real), dimension(ndata), intent(out) :: LTxy,LTyx
!
type (StateElem), pointer                 :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer                 :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter                        :: Delta=1
integer                                   :: ncounts_x,ncounts_y,ncounts_x1,ncounts_y1, &
                                             ncounts_xy,ncounts_xy1,ncounts_yx1
integer                                   :: nstates_x,nstates_y,nstates_x1,nstates_y1, &
                                             nstates_xy,nstates_xy1,nstates_yx1
integer                                   :: d,L,k,m,tau,t,ndims
character(len=1)                          :: xsc,ysc
character(len=250)                        :: xL,x1L,yL,y1L
character(len=250)                        :: xyL, xy1L, yx1L
real(sifm_real), dimension(:), allocatable :: hx1,hy1,hxy,hyx,hxy1,hyx1,hx,hy
real(sifm_real), dimension(:), allocatable :: prob_x, prob_x1, prob_y1, Prob_y
real(sifm_real), dimension(:), allocatable :: Prob_xy, Prob_xy1, prob_yx1
integer, dimension(:), allocatable        :: frequency_x, frequency_x1, frequency_y1, frequency_y
integer, dimension(:), allocatable        :: frequency_xy, frequency_xy1, frequency_yx1
integer, dimension(:), allocatable        :: states_y1,states_y,states_x,states_x1, &
                                             states_xy,states_xy1,states_yx1
real(sifm_real)                            :: Txy
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m = max(imax1d(m1),imax1d(m2))
tau = max(imax1d(tau1),imax1d(tau2))
!
ndims = Ndata-Delta-(M-1)*Tau
call mteallocate_r(hx,   ndims,  "hx",   "symbolic-LTE-entropies")
call mteallocate_r(hx1,  ndims,  "hx1",  "symbolic-LTE-entropies")
call mteallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropies")
call mteallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropies")
call mteallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropies")
call mteallocate_r(hyx,  ndims,  "hyx",  "symbolic-LTE-entropies")
call mteallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropies")
call mteallocate_r(hyx1, ndims,  "hyx1", "symbolic-LTE-entropies")
!
call mteallocate_i(states_x,   ndims, "states-x",   "symbolic-LTE-entropies")
call mteallocate_i(states_x1,  ndims, "states-x1",  "symbolic-LTE-entropies")
call mteallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropies")
call mteallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropies")
call mteallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropies")
call mteallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropies")
call mteallocate_i(states_yx1, ndims, "states-yx1", "symbolic-LTE-entropies")

!!! --- Loop over all data points
t = 0
time_Loop1: DO L=((M-1)*Tau+1),(Ndata-Delta)
!
   t = t + 1
   xL=""; yL=""
   x1L=""; y1L=""
   xyL=""; xy1L=""; yx1L=""
!
   xdim_loop: DO d=1, ndim
      x_Loop1: DO k=1,m1(d)
        xsc  = xs(d,L-(k-1)*tau1(d))
        call addStrings(xL,xsc)
        call addStrings(xyL,xsc)
        call addStrings(x1L,xsc)
        call addStrings(yx1L,xsc)
        call addStrings(xy1L,xsc)
      end do x_Loop1
      xsc  = xs(d,L+Delta)
      call addStrings(x1L,xsc)
      call addStrings(yx1L,xsc)
   end do xdim_loop
!
   ydim_loop: do d=1, ndim
      y_Loop1: DO k=1,m2(d)
        ysc  = ys(d,L-(k-1)*tau2(d))
        call addStrings(yL,ysc)
        call addStrings(xyL,ysc)
        call addStrings(y1L,ysc)
        call addStrings(xy1L,ysc)
        call addStrings(yx1L,ysc)
      end do y_Loop1
      ysc  = ys(d,L+Delta)
      call addStrings(y1L,ysc)
      call addStrings(xy1L,ysc)
   end do ydim_loop
!
   call LocalLinkList(headx,tailx,xL,states_x(t))
   call LocalLinkList(headx1,tailx1,x1L,states_x1(t))
   call LocalLinkList(heady,taily,yL,states_y(t))
   call LocalLinkList(heady1,taily1,y1L,states_y1(t))
   call LocalLinklist(headxy,tailxy,xyL,states_xy(t))
   call LocalLinkList(headxy1,tailxy1,xy1L,states_xy1(t))
   call LocalLinkList(headyx1,tailyx1,yx1L,states_yx1(t))
!
end do time_Loop1
!
ncounts_x   = getTotalStates(Headx)
ncounts_y   = getTotalStates(Heady)
ncounts_x1  = getTotalStates(Headx1)
ncounts_y1  = getTotalStates(Heady1)
ncounts_xy  = getTotalStates(Headxy)
ncounts_xy1 = getTotalStates(Headxy1)
ncounts_yx1 = getTotalStates(Headyx1)
!
nstates_x = CountStates(headx)
nstates_y = CountStates(heady)
nstates_x1 = CountStates(headx1)
nstates_y1 = CountStates(heady1)
nstates_xy = CountStates(headxy)
nstates_xy1 = CountStates(headxy1)
nstates_yx1 = CountStates(headyx1)
!
call mteallocate_r(prob_x,   ncounts_x,   "prob-x",   "symbolic-LTE-entropies")
call mteallocate_r(prob_y,   ncounts_y,   "prob-y",   "symbolic-LTE-entropies")
call mteallocate_r(prob_x1,  ncounts_x1,  "prob-x1",  "symbolic-LTE-entropies")
call mteallocate_r(prob_y1,  ncounts_y1,  "prob-y1",  "symbolic-LTE-entropies")
call mteallocate_r(prob_xy,  ncounts_xy,  "prob-xy",  "symbolic-LTE-entropies")
call mteallocate_r(prob_xy1, ncounts_xy1, "prob-xy1", "symbolic-LTE-entropies")
call mteallocate_r(prob_yx1, ncounts_yx1, "prob-yx1", "symbolic-LTE-entropies")
!
call mteallocate_i(frequency_x,   ncounts_x,   "frequency-x",   "symbolic-LTE-entropies")
call mteallocate_i(frequency_y,   ncounts_y,   "frequency-y",   "symbolic-LTE-entropies")
call mteallocate_i(frequency_x1,  ncounts_x1,  "frequency-x1",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_y1,  ncounts_y1,  "frequency-y1",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_xy,  ncounts_xy,  "frequency-xy",  "symbolic-LTE-entropies")
call mteallocate_i(frequency_xy1, ncounts_xy1, "frequency-xy1", "symbolic-LTE-entropies")
call mteallocate_i(frequency_yx1, ncounts_yx1, "frequency-yx1", "symbolic-LTE-entropies")
!
call SymbolicJointProbabilityTraj(Headx, ncounts_x, Prob_x, frequency_x)
call SymbolicJointProbabilityTraj(Headx1, ncounts_x1, Prob_x1, frequency_x1)
call SymbolicJointProbabilityTraj(Heady, ncounts_y, Prob_y, frequency_y)
call SymbolicJointProbabilityTraj(Heady1, ncounts_y1, Prob_y1, frequency_y1)
call SymbolicJointProbabilityTraj(Headxy, ncounts_xy, Prob_xy, frequency_xy)
call SymbolicJointProbabilityTraj(Headxy1, ncounts_xy1, Prob_xy1, frequency_xy1)
call SymbolicJointProbabilityTraj(Headyx1, ncounts_yx1, Prob_yx1, frequency_yx1)
!
t = 0
time_Loop2: DO L=((M-1)*Tau+1),(Ndata-Delta)
   t = t + 1
   Hx(t)   = -log( Prob_x(states_x(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)
   Hx1(t)  = -log( Prob_x1(states_x1(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)
   Hyx(t)  = -log( Prob_xy(states_xy(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)
   Hyx1(t) = -log( Prob_yx1(states_yx1(t)) ) / real(frequency_yx1(states_yx1(t)), sifm_real)

   Hy(t)   = -log( Prob_y(states_y(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
   Hy1(t)  = -log( Prob_y1(states_y1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
   Hxy(t)  = -log( Prob_xy(states_xy(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
   Hxy1(t) = -log( Prob_xy1(states_xy1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
end do time_Loop2
!
do t = 1, nframes
   if (t <= ndims) then
       LTxy(t) = (hy1(t)-hy(t)) - (hxy1(t)-hxy(t))
       LTyx(t) = (hx1(t)-hx(t)) - (hyx1(t)-hyx(t))
   else
       LTxy(t) = rzero
       LTyx(t) = rzero
   endif
enddo
!
!
IF (debug == 1) THEN
    Txy = dot_product( LTxy(1:ndims), Prob_xy1(states_xy1(1:ndims)) )
    write(*,'("Total transfer entropy from local:   :", F12.6)') Txy
    Txy = rzero
    call symbolic_TE_entropynD(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy)
    write(*,'("Total transfer entropy from global:  :", F12.6)') Txy 
ENDIF
!
call freeList(headx)
call freeList(heady)
call freeList(headx1)
call freeList(heady1)
call freeList(headxy)
call freeList(headxy1)
call freeList(headyx1)
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
call mtedeallocate_r(hx,   ndims,  "hx",   "symbolic-LTE-entropies")
call mtedeallocate_r(hx1,  ndims,  "hx1",  "symbolic-LTE-entropies")
call mtedeallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropies")
call mtedeallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropies")
call mtedeallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropies")
call mtedeallocate_r(hyx,  ndims,  "hyx",  "symbolic-LTE-entropies")
call mtedeallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropies")
call mtedeallocate_r(hyx1, ndims,  "hyx1", "symbolic-LTE-entropies")
!
call mtedeallocate_i(states_x,   ndims, "states-x",   "symbolic-LTE-entropies")
call mtedeallocate_i(states_x1,  ndims, "states-x1",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropies")
call mtedeallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropies")
call mtedeallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropies")
call mtedeallocate_i(states_yx1, ndims, "states-yx1", "symbolic-LTE-entropies")
!
call mtedeallocate_r(prob_x,   ncounts_x,   "prob-x",   "symbolic-LTE-entropies")
call mtedeallocate_r(prob_y,   ncounts_y,   "prob-y",   "symbolic-LTE-entropies")
call mtedeallocate_r(prob_x1,  ncounts_x1,  "prob-x1",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_y1,  ncounts_y1,  "prob-y1",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_xy,  ncounts_xy,  "prob-xy",  "symbolic-LTE-entropies")
call mtedeallocate_r(prob_xy1, ncounts_xy1, "prob-xy1", "symbolic-LTE-entropies")
call mtedeallocate_r(prob_yx1, ncounts_yx1, "prob-yx1", "symbolic-LTE-entropies")
!
call mtedeallocate_i(frequency_x,   ncounts_x,   "frequency-x",   "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_y,   ncounts_y,   "frequency-y",   "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_x1,  ncounts_x1,  "frequency-x1",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_y1,  ncounts_y1,  "frequency-y1",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_xy,  ncounts_xy,  "frequency-xy",  "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_xy1, ncounts_xy1, "frequency-xy1", "symbolic-LTE-entropies")
call mtedeallocate_i(frequency_yx1, ncounts_yx1, "frequency-yx1", "symbolic-LTE-entropies")
!
return
end subroutine symbolic_LTE_entropies
!
subroutine symbolic_LTE_entropy(ndata,ndim,xs,ys,m1,m2,tau1,tau2,LTxy)
implicit none
!
integer, intent(in)            :: ndata,ndim
integer, intent(in)            :: m1(ndim),m2(ndim),tau1(ndim),tau2(ndim)
character(len=1), intent(in)   :: xs(ndim,ndata),ys(ndim,ndata)
real(sifm_real), intent(inout) :: LTxy(ndata)
!
type (StateElem), pointer                 :: heady,taily,heady1,taily1
type (StateElem), pointer                 :: headxy,tailxy,headxy1,tailxy1
integer, parameter                        :: Delta=1
integer                                   :: ncounts_y,ncounts_y1,ncounts_xy,ncounts_xy1
integer                                   :: nstates_y,nstates_y1,nstates_xy,nstates_xy1
integer                                   :: d,L,k,m,tau,t,ierror,ndims
character(len=1)                          :: xsc,ysc
character(len=250)                        :: yL,yy1L
character(len=250)                        :: xyL,xyy1L
real(sifm_real), dimension(:), allocatable :: hy1,hxy,hxy1,hy
real(sifm_real), dimension(:), allocatable :: prob_y1,Prob_xy,Prob_xy1,Prob_y
integer, dimension(:), allocatable        :: frequency_y1,frequency_xy,frequency_xy1,frequency_y
integer, dimension(:), allocatable        :: states_y1,states_xy,states_xy1,states_y
real(sifm_real)                            :: Txy
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(imax1d(m1),imax1d(m2))
tau = max(imax1d(tau1),imax1d(tau2))
!
ndims = Ndata-Delta-(M-1)*Tau
call mteallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropy")
call mteallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropy")
call mteallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropy")
call mteallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropy")
!
call mteallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropy")
call mteallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropy")
call mteallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropy")
call mteallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropy")
!
!!! --- Loop over all data points
t = 0
time_Loop1: DO L=((M-1)*Tau+1),(Ndata-Delta)
!
  t = t + 1
  yL=""; yy1L=""
  xyL=""; xyy1L=""
! 
  xdim_loop: DO d = 1, ndim
     x_Loop1: DO k=1,m1(d)
         xsc   = xs(d,L-(k-1)*tau1(d))
         call addStrings(xyL,xsc)
         call addStrings(xyy1L,xsc)
     end do x_Loop1
  end do xdim_loop
!
  ydim_loop: DO d = 1, ndim
     y_Loop1: DO k=1,m2(d)
         ysc   = ys(d,L-(k-1)*tau2(d))
         call addStrings(yL,ysc)
         call addStrings(xyL,ysc)
         call addStrings(yy1L,ysc)
         call addStrings(xyy1L,ysc)
     end do y_Loop1
     ysc   = ys(d,L+Delta)
     call addStrings(yy1L,ysc)
     call addStrings(xyy1L,ysc)
  end do ydim_loop
!
  call LocalLinkList(heady,taily,yL,states_y(t))  
  call LocalLinkList(heady1,taily1,yy1L,states_y1(t))
  call LocalLinklist(headxy,tailxy,xyL,states_xy(t))
  call LocalLinkList(headxy1,tailxy1,xyy1L,states_xy1(t))  
!
end do time_Loop1
!
ncounts_y   = getTotalStates(Heady)
ncounts_y1  = getTotalStates(Heady1)
ncounts_xy  = getTotalStates(Headxy)
ncounts_xy1 = getTotalStates(Headxy1)
!
nstates_y   = CountStates(heady)
nstates_y1  = CountStates(heady1)
nstates_xy  = CountStates(headxy)
nstates_xy1 = CountStates(headxy1)
!
call mteallocate_r(prob_y, nstates_y, "prob-y", "symbolic-LTE-entropy")
call mteallocate_r(prob_y1, nstates_y1, "prob-y1", "symbolic-LTE-entropy")
call mteallocate_r(prob_xy, nstates_xy, "prob-xy", "symbolic-LTE-entropy")
call mteallocate_r(prob_xy1, nstates_xy1, "prob-xy1", "symbolic-LTE-entropy")
!
call mteallocate_i(frequency_y, nstates_y, "frequency-y", "symbolic-LTE-entropy")
call mteallocate_i(frequency_y1, nstates_y1, "frequency-y1", "symbolic-LTE-entropy")
call mteallocate_i(frequency_xy, nstates_xy, "frequency-xy", "symbolic-LTE-entropy")
call mteallocate_i(frequency_xy1, nstates_xy1, "frequency-xy1", "symbolic-LTE-entropy")
!
call SymbolicJointProbabilityTraj(Heady, ncounts_y, Prob_y, frequency_y)
call SymbolicJointProbabilityTraj(Heady1, ncounts_y1, Prob_y1, frequency_y1)
call SymbolicJointProbabilityTraj(Headxy, ncounts_xy, Prob_xy, frequency_xy)
call SymbolicJointProbabilityTraj(Headxy1, ncounts_xy1, Prob_xy1, frequency_xy1)
!
t = 0
time_Loop2: DO L=((M-1)*Tau+1),(Ndata-Delta)
  t = t + 1
  hy(t)   = -dlog( Prob_y(states_y(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  hy1(t)  = -dlog( Prob_y1(states_y1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  hxy(t)  = -dlog( Prob_xy(states_xy(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
  hxy1(t) = -dlog( Prob_xy1(states_xy1(t)) ) / real(frequency_xy1(states_xy1(t)), sifm_real)
end do time_Loop2
!
do t = 1, nframes
   IF (t <= ndims) THEN
       LTxy(t) = LTxy(t) + (hy1(t)-hy(t)) - (hxy1(t)-hxy(t))
   ELSE
       LTxy(t) = rzero
   ENDIF
enddo
!
IF (debug == 1) THEN
    Txy = dot_product( LTxy(1:ndims), Prob_xy1(states_xy1(1:ndims)) )
    write(*,'("Total transfer entropy from local:   :", F12.6)') Txy
    Txy = rzero
    call symbolic_TE_entropynD(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy)
    write(*,'("Total transfer entropy from global:  :", F12.6)') Txy 
ENDIF

call freeList(heady)
call freeList(heady1)
call freeList(headxy)
call freeList(headxy1)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
call mtedeallocate_r(hy1,  ndims,  "hy1",  "symbolic-LTE-entropy")
call mtedeallocate_r(hxy,  ndims,  "hxy",  "symbolic-LTE-entropy")
call mtedeallocate_r(hxy1, ndims,  "hxy1", "symbolic-LTE-entropy")
call mtedeallocate_r(hy,   ndims,  "hy",   "symbolic-LTE-entropy")
!
call mtedeallocate_i(states_y1,  ndims, "states-y1",  "symbolic-LTE-entropy")
call mtedeallocate_i(states_xy,  ndims, "states-xy",  "symbolic-LTE-entropy")
call mtedeallocate_i(states_xy1, ndims, "states-xy1", "symbolic-LTE-entropy")
call mtedeallocate_i(states_y,   ndims, "states-y",   "symbolic-LTE-entropy")
!
call mtedeallocate_r(prob_y, nstates_y, "prob-y", "symbolic-LTE-entropy")
call mtedeallocate_r(prob_y1, nstates_y1, "prob-y1", "symbolic-LTE-entropy")
call mtedeallocate_r(prob_xy, nstates_xy, "prob-xy", "symbolic-LTE-entropy")
call mtedeallocate_r(prob_xy1, nstates_xy1, "prob-xy1", "symbolic-LTE-entropy")
!
call mtedeallocate_i(frequency_y, nstates_y, "frequency-y", "symbolic-LTE-entropy")
call mtedeallocate_i(frequency_y1, nstates_y1, "frequency-y1", "symbolic-LTE-entropy")
call mtedeallocate_i(frequency_xy, nstates_xy, "frequency-xy", "symbolic-LTE-entropy")
call mtedeallocate_i(frequency_xy1, nstates_xy1, "frequency-xy1", "symbolic-LTE-entropy")
!
return
end subroutine symbolic_LTE_entropy
!
subroutine symbolic_Local_entropy1D(ndata,xs,m,tau,LH,H)
implicit none
!
integer, intent(in)          :: ndata,m,tau
character(len=1), intent(in) :: xs(ndata)
real(sifm_real), intent(out) :: LH(:)
real(sifm_real), intent(out) :: H
!
type (StateElem), pointer :: head, tail
integer :: nstates,ncounts
integer :: istat,L,k,t,ierror,ndims
character(len=1) :: xsc
character(len=250) :: xL
real(sifm_real), dimension(:), allocatable :: prob
integer, dimension(:), allocatable :: frequency
integer, dimension(:), allocatable :: states
!
nullify(head,tail)
!
ndims = Ndata-(M-1)*Tau
!
Allocate( states(ndims), stat=ierror)
IF (ierror /= 0) STOP "Error allocating states"
!
!!! --- Loop over all data points
t = 0
time_Loop1: DO L=((M-1)*Tau+1),Ndata
!
  t = t + 1
  xL=""
! 
  x_Loop1: DO k=1,m
         xsc = xs(L-(k-1)*tau)
         call addStrings(xL,xsc)
  end do x_Loop1
!
  call LocalLinkList(head,tail,xL,states(t))  
!
end do time_Loop1
!
ncounts = getTotalStates(Head)
write(*,*) "Ntotal=   ", ncounts
nstates = CountStates(Head)
!
Allocate( Prob(nstates), stat=ierror )
IF (ierror /= 0) STOP "Error allocating Prob"
Allocate( Frequency(nstates), stat=ierror )
IF (ierror /= 0) STOP "Error allocating frequency"
!
call SymbolicJointProbabilityTraj(Head, ncounts, Prob, frequency)
!
t = 0
H = rzero
time_Loop2: DO L=((M-1)*Tau+1),Ndata
  t = t + 1
  LH(t) = -dlog( Prob(states(t)) ) / real(Frequency(states(t)), sifm_real)
  H = H -Prob(states(t)) * dlog( Prob(states(t)) ) / real(Frequency(states(t)), sifm_real)
end do time_Loop2
write(*,'(F12.6)') H
!
H = dot_product( Prob(states(1:ndims)), LH(1:ndims) ) !/ real(M, sifm_real)
write(*,'(F12.6)') H 
!
H = SymbolicShannonEntropy(Head) !/ real(M, sifm_real)
write(*,'(F12.6)') H
!
call freeList(head)
!
nullify(head,tail)
!
DeAllocate( Prob, stat=ierror)
IF (ierror /= 0) STOP "Error deallocating Prob"
!
end subroutine symbolic_Local_entropy1d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropyND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)        :: ndata,ndim
integer, dimension(Ndim), intent(in) :: m1,m2,tau1,tau2
character(len=1), dimension(Ndim,Ndata), intent(in) :: xs,ys
real(sifm_real), intent(inout) :: Txy
!
type (StateElem), pointer :: heady,taily,heady1,taily1
type (StateElem), pointer :: headxy,tailxy,headxy1,tailxy1
integer, parameter        :: Delta = 1
integer                   :: d,L,k,m,tau
character(len=1)          :: sc
character(len=250)        :: yL,y1L
character(len=250)        :: xyL, xy1L
real(sifm_real)            :: hy,hy1,hxy,hxy1
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
m = max(imax1d(m1), imax1d(m2))
tau = max(imax1d(tau1), imax1d(tau2))
!
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
call addStrings(yL,sc)
call addStrings(xyL, sc)
call addStrings(y1L,sc)
call addStrings(xy1L, sc)
end do y_Loop
sc    = ys(d,L+Delta)
call addStrings(y1L,sc)
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
Txy=Txy + (hy1-hy) - (hxy1-hxy)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
end subroutine symbolic_TE_entropyND
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE2_entropyND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)        :: ndata,ndim
integer, dimension(Ndim), intent(in) :: m1,m2,tau1,tau2
real(sifm_real), dimension(Ndim,Ndata), intent(in) :: xs,ys
real(sifm_real), intent(inout) :: Txy
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
Txy=Txy + (hy1-hy) - (hxy1-hxy)
!
nullify(heady,taily,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1)
!
end subroutine symbolic_TE2_entropyND
!
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
!
end module LTE_class
