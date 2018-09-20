!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: minfo.f90,v 1.0 19-03-2018, IBU
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
module MI_class
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
      integer, save, allocatable :: mopt(:,:), tau(:,:)
      real(sifm_real), save, dimension(:,:), allocatable :: MI
      real(sifm_real), save, dimension(:,:), allocatable :: MIS
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMI_DRIVER(Nf, Nd, Na, qMIShuffle, debug, Nshuffles, r0, statP)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: qMIShuffle
  integer, intent(in)        :: Na, Nf
  integer, intent(in)        :: Nd, debug
  integer, optional          :: Nshuffles
  real(sifm_real), optional  :: r0, statP
!
  Integer :: I
  integer :: ierr
  real(sifm_real) :: LocalTime_Start, LocalTime_end
  real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  Nframes = Nf
!
  Ndim = Nd
  Natoms = Na
  Ndof = Natoms * Ndim
!
  IF (debug == 1) THEN
      write(*,'("# of frames to be investigated                             :  ", I5)') nframes
      write(*,'("Dimension of the fluctuations                              :  ", I5)') ndim
      write(*,'("# of degrees of freedom                                    :  ", I5)') ndof
  ENDIF
!
  CALL Allocate_MI(1, qMIShuffle)
!
  CALL read_embdparam()
  CALL read_xyzs()
!
! -- Initialize MPI
  call MPI_start()
  call MPI_BroadCast(qMIShuffle, debug, Nshuffles, r0, statP, 1)
!
  Localtime_start = MPI_Wtime()
!
  IF (qMIShuffle > 0) THEN
      CALL getSMINDIM_MPISFL(qMIShuffle, debug, Nshuffles, r0, statp)
  ELSE
      CALL getSMINDIM_MPIDOF(debug)
  ENDIF
!
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for MI:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  IF (myID == Master) CALL write_mi(qMIShuffle) 
!
  CALL DEAllocate_MI()
!
  call MPI_finish()
!
  RETURN
!
end subroutine SMI_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MI_DRIVER(Nf, Nd, Na, qMIShuffle, debug, Nshuffles, r0, statP)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: qMIShuffle
  integer, intent(in)        :: Na, Nf
  integer, intent(in)        :: Nd, debug
  integer, optional          :: Nshuffles
  real(sifm_real), optional  :: r0, statP
!
  Integer :: I
  integer :: ierr
  real(sifm_real) :: LocalTime_Start, LocalTime_end
  real(sifm_real) :: GlobalTime_Start, GlobalTime_end
!
  Nframes = Nf
!
  Ndim = Nd
  Natoms = Na
  Ndof = Natoms * Ndim
!
  IF (debug == 1) THEN
      write(*,'("# of frames to be investigated                             :  ", I5)') nframes
      write(*,'("Dimension of the fluctuations                              :  ", I5)') ndim
      write(*,'("# of degrees of freedom                                    :  ", I5)') ndof
  ENDIF
!
  CALL Allocate_MI(2,qMIShuffle)
!
  CALL read_embdparam()
  CALL read_xyz()
!
!
! -- Initialize MPI
  call MPI_start()
  call MPI_BroadCast(qMIShuffle, debug, Nshuffles, r0, statP, 2)
!
  Localtime_start = MPI_Wtime()
!
  IF (qMIShuffle > 0) THEN
      CALL getMINDIM_MPISFL(qMIShuffle, debug, Nshuffles, r0, statp)
  ELSE
      CALL getMINDIM_MPIDOF(debug)
  ENDIF
!
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for MI:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  IF (myID == Master) CALL write_mi(qMIShuffle) 
!
  CALL DEAllocate_MI()
!
  call MPI_finish()
!
  RETURN
!
end subroutine MI_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the Mutual Information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_MI(qMIShuffle)
  implicit none
!
  integer, intent(in) :: qMIShuffle
!  
  integer             :: i, j, nunit
!  
  nunit = new_unit()
  open(unit=nunit, file = 'mi.txt', status='unknown', action='write')

  Write(*,'("Print Mutual Information: -----------------------------")')
  DO i=1, Natoms
     DO j=1, Natoms
        write(nunit,'(2I10, F12.6)') i,j, MI(i,j)
     ENDDO
  ENDDO
  close(nunit)
!
  IF (qMIShuffle > 0) THEN
  nunit = new_unit()
  open(unit=nunit, file = 'mis.txt', status='unknown', action='write')

  Write(*,'("Print Mutual Information: -----------------------------")')
  DO i=1, Natoms
     DO j=1, Natoms
        write(nunit,'(2I10, F12.6)') i,j, MIS(i,j)
     ENDDO
  ENDDO
  close(nunit)
  ENDIF
!
return
end subroutine write_MI
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Allocate_MI(qmimethod, qMIShuffle)
      implicit none
!
      integer, intent(in) :: qmimethod, qMIShuffle
!
      integer :: ierr
!
      IF ( qMIMETHOD == 2 ) THEN
           IF (.not. Allocated(xyz)) THEN
               Allocate(XYZ(Natoms,ndim,nframes), stat=ierr)
               IF (ierr /= 0) Stop "Error Allocating X"
           ENDIF
      ELSE
!
           IF (.not. Allocated(xyzs)) THEN
               Allocate(XYZs(Natoms,ndim,nframes), stat=ierr)
               IF (ierr /= 0) Stop "Error Allocating Xs"
           ENDIF
      ENDIF
!  
      IF (.not. Allocated(tau)) THEN
      Allocate(Tau(Natoms, ndim), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating Tau"
      ENDIF
!
      IF (.not. Allocated(Mopt)) THEN
          Allocate(mopt(Natoms, ndim), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating Mopt"
      ENDIF
!
      IF (.not. Allocated(MI)) THEN
          ALLOCATE( MI(natoms, Natoms), stat=ierr )
          if (ierr /= 0) stop 'Error allocating Mij'
      ENDIF
!
      IF (qMIShuffle > 0) THEN
          IF (.not. Allocated(MIS)) THEN
              ALLOCATE( MIS(natoms, Natoms), stat=ierr )
              if (ierr /= 0) stop 'Error allocating MSij'
          ENDIF
      ENDIF
!
      return
end subroutine Allocate_MI
!
subroutine deAllocate_MI()
      implicit none
!
      integer :: ierr
!
      IF (Allocated(xyz)) THEN
          DEAllocate(XYZ, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating X"
      ENDIF
!
      IF (Allocated(xyzs)) THEN
          DEAllocate(XYZS, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating Xs"
      ENDIF
!  
      IF (Allocated(tau)) THEN
          DEAllocate(Tau, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating Tau"
      ENDIF
!
      IF (Allocated(Mopt)) THEN
          DEAllocate(mopt, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating Mopt"
      ENDIF
!
      IF (Allocated(MI)) THEN
          DEALLOCATE( MI, stat=ierr )
          if (ierr /= 0) stop 'Error deallocating Mij'
      ENDIF
!
      IF (Allocated(MIS)) THEN
          DEALLOCATE( MIS, stat=ierr )
          if (ierr /= 0) stop 'Error deallocating MSij'
      ENDIF
!
      return
end subroutine deAllocate_MI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the real trajectory
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
     read(nunit,*) (tau(i,d), d=1, Ndim), (mopt(i,d), d=1, Ndim)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MPI_BroadCast(qMIShuffle, debug, Nshuffles, r0, statP, model)
use TEMPI_CLASS
implicit none
!
      integer, intent(in)        :: model
      integer, intent(in)        :: qMIShuffle
      integer, intent(in)        :: debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!
!
integer :: ierr
!
! --- If running on parallel Broadcast some information
call MPI_BCAST (Natoms,     1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Ndim,       1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Nframes,    1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (debug,      1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (MI,         Natoms*Natoms, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (qMIShuffle, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
IF (qMIShuffle > 0) THEN
    call MPI_BCAST (Nshuffles,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (MIS,        Natoms*Natoms, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (statP,      1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
    call MPI_BCAST (r0,         1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
ENDIF
!
if (model == 1) then
    call MPI_BCAST (XYZS,   Natoms*Ndim*Nframes, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr )
else
    call MPI_BCAST (XYZ,    Natoms*Ndim*Nframes, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
endif
!
call MPI_BCAST (Mopt,       Natoms*Ndim, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Tau,        Natoms*Ndim, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
return
END Subroutine MPI_Broadcast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSMINDIM_MPIDOF(debug)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: debug
!
!{Local variables}
integer :: ntrials
integer :: i, j, t
!
integer, allocatable :: listPairs(:,:)
integer :: Npairs, iList, K
integer :: ierr, Stride, iStart, iEnd
real(sifm_real) :: tMxy
!
real(sifm_real), allocatable, dimension(:) :: gMxy, lMxy
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
IF ( mod(Npairs, NumProcs) == 0 ) THEN
     Stride = Npairs / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!       
     ALLOCATE( LMxy(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local  Mxy'
     ALLOCATE( GMxy(Npairs), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Mxy'
!
ELSE
     Stop "Number of Trial Shuffles must be Multiple of Number of Processors"
ENDIF
! --- Compute the Mutual Information
!-------- call all the pairs
Pairs_LOOP: DO iList = iStart, iEnd
        K = iList - iStart + 1
!
        I = ListPairs(iList, 1)
        J = ListPairs(iList, 2)
!        
        CALL symbolic_MI_ND(nframes,ndim,XYZs(i,:,:),XYZs(j,:,:),            &
                            Mopt(i,:),Mopt(j,:),Tau(i,:),Tau(j,:),           &
                            tMxy)
        LMxy(K) = tMxy
!
ENDDO PAIRS_LOOP
!
CALL MPI_gather( LMxy, Stride, MPI_DOUBLE_PRECISION, &
                 GMxy, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
IF (myID == master) THEN
    DO iList = 1, Npairs
       I = ListPairs(iList, 1)
       J = ListPairs(iList, 2)
       MI(i,j) = GMxy(iList)
       MI(j,i) = MI(i,j)
    ENDDO
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
DEALLOCATE( LMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Local Mxy'
DEALLOCATE( GMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Global Mxy'

!
RETURN
!
END subroutine getSMINDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSMINDIM_MPISFL(qMIShuffle, debug, Nshuffles, r0, statP)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: qMIShuffle
  integer, intent(in)        :: debug
  integer, optional          :: Nshuffles
  real(sifm_real), optional  :: r0, statP
!
!{Local variables}
integer :: ntrials
integer :: i, j, t, itrial
character(len=1), dimension(ndim,nframes) :: xis_shuffle
character(len=1), dimension(ndim,nframes) :: xjs_shuffle
!
integer, allocatable :: listPairs(:,:)
integer :: Npairs, iList
integer :: ierr, Stride, iStart, iEnd
real(sifm_real) :: tMxy
!
real(sifm_real), allocatable, dimension(:) :: gMxy, lMxy
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

! -- Initialize 2*natoms random number sequences
IF ( qMIshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, 0)
     ntrials = nshuffles
     MIS = rzero
endif

IF ( mod(Ntrials, NumProcs) == 0 ) THEN
     Stride = Ntrials / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LMxy(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Mxy'
     ALLOCATE( GMxy(Ntrials), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Mxy'
!
ELSE
     Stop "Number of Trial Shuffles must be Multiple of Number of Processors"
ENDIF
! --- Compute the Mutual Information
MI = rzero

!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
        I = ListPairs(iList, 1)
        J = ListPairs(iList, 2)
!        
        IF (myID == Master) THEN
            CALL symbolic_MI_ND(nframes,ndim,XYZs(i,:,:),XYZs(j,:,:),            &
                                Mopt(i,:),Mopt(j,:),Tau(i,:),Tau(j,:),           &
                                MI(i,j))
        ENDIF
!
        !--- Compute the direct mutual information
        IF (qMIshuffle > 0) THEN
             xis_shuffle = XYZS(i,:,:)
             xjs_shuffle = XYZS(j,:,:)
             Do iTrial = iStart, iEnd
                t = iTrial - iStart + 1
                IF (QMIShuffle == 1) THEN
                    CALL shuffleND(nframes,ndim,xis_shuffle,I)
                    CALL shuffleND(nframes,ndim,xjs_shuffle,J+Natoms)
                ELSEIF (QMIShuffle == 2) THEN
                    CALL BlockShuffleND(nframes, ndim, xis_shuffle, Mopt(i,:), Tau(i,:), I)
                    CALL BlockShuffleND(nframes, ndim, xjs_shuffle, Mopt(j,:), Tau(j,:), J+Natoms)
                ELSE
                    STOP "Error: Unknown shuffling method; QMISHUFFLE has to be 1 or 2"
                ENDIF
                CALL symbolic_MI_ND(nframes,ndim,xis_shuffle,xjs_shuffle,            &
                                    Mopt(i,:),Mopt(j,:),Tau(i,:),Tau(j,:),           &
                                    LMxy(t))
             ENDDO
             CALL MPI_gather( LMxy, Stride, MPI_DOUBLE_PRECISION, &
                              GMxy, Stride, MPI_DOUBLE_PRECISION, &
                              master, MPI_COMM_WORLD, ierr )             
             IF (myID == Master) THEN
                 tMxy = mean( gMxy, Ntrials )
!
                 write(*,'("Average value of Bias Mxy (Before Statistical Test)  =  ", F12.6)') tMxy
                 CALL statTest(Ntrials, R0, tMxy, tMxy, StatP)	
                 write(*,'("Average value of Bias Mxy (After Statistical Test)   =  ", F12.6)') tMxy
!
                 MIS(i,j) = MI(i,j) - tMxy
             ENDIF
        ENDIF
        IF (debug >= 1 .and. myID == Master) THEN
            IF ( qMIshuffle > 0 ) THEN
                 write(*,'("Symbolic Mutual Information:   ", I5, I5, 3F12.6)') I, J, tMxy, MI(i,j), MIS(i,j)
            ELSE
                 write(*,'("Symbolic Mutual Information:   ", I5, I5, F12.6)') I, J, MI(i,j)
            ENDIF
        ENDIF
!
        IF (myID == Master) THEN
            IF ( qMIshuffle > 0 ) MIS(j,i) = MIS(i,j)
            MI(j,i) = MI(i,j)
        ENDIF
!
ENDDO PAIRS_LOOP
!
IF ( qMIshuffle > 0) then
     CALL rndnum_clear()
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
DEALLOCATE( LMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Local Mxy'
DEALLOCATE( GMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Global Mxy'
!
RETURN
!
END subroutine getSMINDIM_MPISFL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getMINDIM_MPIDOF(debug)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: debug
!
!{Local variables}
integer :: ntrials
integer :: i, j, t
!
integer, allocatable :: listPairs(:,:)
integer :: Npairs, iList, K
integer :: ierr, Stride, iStart, iEnd
real(sifm_real) :: tMxy
!
real(sifm_real), allocatable, dimension(:) :: gMxy, lMxy
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
IF ( mod(Npairs, NumProcs) == 0 ) THEN
     Stride = Npairs / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!       
     ALLOCATE( LMxy(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local  Mxy'
     ALLOCATE( GMxy(Npairs), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Mxy'
!
ELSE
     Stop "Number of Trial Shuffles must be Multiple of Number of Processors"
ENDIF
! --- Compute the Mutual Information
!-------- call all the pairs
Pairs_LOOP: DO iList = iStart, iEnd
        K = iList - iStart + 1
!
        I = ListPairs(iList, 1)
        J = ListPairs(iList, 2)
!        
        CALL symbolic_MI2_ND(nframes,ndim,XYZ(i,:,:),XYZ(j,:,:),              &
                             Mopt(i,:),Mopt(j,:),Tau(i,:),Tau(j,:),           &
                             tMxy)
        LMxy(K) = tMxy
!
ENDDO PAIRS_LOOP
!
CALL MPI_gather( LMxy, Stride, MPI_DOUBLE_PRECISION, &
                 GMxy, Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
IF (myID == master) THEN
    DO iList = 1, Npairs
       I = ListPairs(iList, 1)
       J = ListPairs(iList, 2)
       MI(i,j) = GMxy(iList)
       MI(j,i) = MI(i,j)
    ENDDO
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
DEALLOCATE( LMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Local Mxy'
DEALLOCATE( GMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Global Mxy'
!
RETURN
!
END subroutine getMINDIM_MPIDOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getMINDIM_MPISFL(qMIShuffle, debug, Nshuffles, r0, statP)
  USE TEMPI_CLASS
  implicit none
  integer, intent(in)        :: qMIShuffle
  integer, intent(in)        :: debug
  integer, optional          :: Nshuffles
  real(sifm_real), optional  :: r0, statP
!
!{Local variables}
integer :: ntrials
integer :: i, j, t, itrial
real(sifm_real), dimension(ndim,nframes) :: xis_shuffle
real(sifm_real), dimension(ndim,nframes) :: xjs_shuffle
!
integer, allocatable :: listPairs(:,:)
integer :: Npairs, iList
integer :: ierr, Stride, iStart, iEnd
real(sifm_real) :: tMxy
!
real(sifm_real), allocatable, dimension(:) :: gMxy, lMxy
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

! -- Initialize 2*natoms random number sequences
IF ( qMIshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, 0)
     ntrials = nshuffles
     MIS = rzero
endif

IF ( mod(Ntrials, NumProcs) == 0 ) THEN
     Stride = Ntrials / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd 
!    
     ALLOCATE( LMxy(Stride), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Local Mxy'
     ALLOCATE( GMxy(Ntrials), stat=ierr )
     if (ierr /= 0) stop 'Error allocating Global Mxy'
!
ELSE
     Stop "Number of Trial Shuffles must be Multiple of Number of Processors"
ENDIF
! --- Compute the Mutual Information
MI = rzero

!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
        I = ListPairs(iList, 1)
        J = ListPairs(iList, 2)
!        
        IF (myID == Master) THEN
            CALL symbolic_MI2_ND(nframes,ndim,XYZ(i,:,:),XYZ(j,:,:),              &
                                 Mopt(i,:),Mopt(j,:),Tau(i,:),Tau(j,:),           &
                                 MI(i,j))
        ENDIF
!
        !--- Compute the direct mutual information
        IF (qMIshuffle > 0) THEN
             xis_shuffle = XYZ(i,:,:)
             xjs_shuffle = XYZ(j,:,:)
             Do iTrial = iStart, iEnd
                t = iTrial - iStart + 1
                IF (QMIShuffle == 1) THEN
                    CALL shuffleND2(nframes,ndim,xis_shuffle,I)
                    CALL shuffleND2(nframes,ndim,xjs_shuffle,J+Natoms)
                ELSEIF (QMIShuffle == 2) THEN
                    CALL BlockShuffleND2(nframes, ndim, xis_shuffle, Mopt(i,:), Tau(i,:), I)
                    CALL BlockShuffleND2(nframes, ndim, xjs_shuffle, Mopt(j,:), Tau(j,:), J+Natoms)
                ELSE
                    STOP "Error: Unknown shuffling method; QMISHUFFLE has to be 1 or 2"
                ENDIF
                CALL symbolic_MI2_ND(nframes,ndim,xis_shuffle,xjs_shuffle,            &
                                     Mopt(i,:),Mopt(j,:),Tau(i,:),Tau(j,:),           &
                                     LMxy(t))
             ENDDO
             CALL MPI_gather( LMxy, Stride, MPI_DOUBLE_PRECISION, &
                              GMxy, Stride, MPI_DOUBLE_PRECISION, &
                              master, MPI_COMM_WORLD, ierr )             
             IF (myID == Master) THEN
                 tMxy = mean( gMxy, Ntrials )
!
                 CALL statTest(Ntrials, R0, tMxy, tMxy, StatP)	
!
                 MIS(i,j) = MI(i,j) - tMxy
             ENDIF
        ENDIF
        IF (debug >= 1 .and. myID == Master) THEN
            IF ( qMIshuffle > 0 ) THEN
                 write(*,'("Symbolic Mutual Information:   ", I5, I5, 3F12.6)') I, J, tMxy, MI(i,j), MIS(i,j)
            ELSE
                 write(*,'("Symbolic Mutual Information:   ", I5, I5, F12.6)') I, J, MI(i,j)
            ENDIF
        ENDIF
!
        IF (myID == Master) THEN
            IF ( qMIshuffle > 0 ) MIS(j,i) = MIS(i,j)
            MI(j,i) = MI(i,j)
        ENDIF
!
ENDDO PAIRS_LOOP
!
IF ( qMIshuffle > 0) then
     CALL rndnum_clear()
ENDIF
!
DEALLOCATE( listPairs, stat=ierr )
if (ierr /= 0) stop 'Error deallocating List of Pairs'
DEALLOCATE( LMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Local Mxy'
DEALLOCATE( GMxy, stat=ierr )
if (ierr /= 0) stop 'Error deallocating Global Mxy'

!
RETURN
!
END subroutine getMINDIM_MPISFL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_MI_ND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,mi)
implicit none
!
integer, intent(in)           :: ndata,ndim
integer, dimension(ndim)      :: m1,m2,tau1,tau2
character(len=1), intent(in)  :: xs(:,:),ys(:,:)
real(sifm_real), intent(out)  :: mi
!
type (StateElem), pointer  :: headx, tailx, heady, taily
type (StateElem), pointer  :: headxy, tailxy
integer                    :: d, L, k, mx, my, M, mxy, tau, L0
character(len=1)           :: sc
character(len=250)         :: xL, x1L, yL
character(len=250)         :: xyL
real(sifm_real)            :: hx, hy, hxy, NormX, NormY, NormXY
!
nullify(headx,tailx,heady,taily)
nullify(headxy,tailxy)
!
mx  = imax1d(m1)
my  = imax1d(m2)
mxy = max(mx, my)
M   = max(imax1d(m1), imax1d(m2)) 
tau = max(imax1d(tau1), imax1d(tau2))
L0  = (M - 1) * Tau + 1
!
!!! --- Loop over all data points
time_Loop: DO L=L0, Ndata
  xL=""; yL=""; xyL=""
!  
  xdim_Loop1: DO d=1,Ndim 
  x_Loop1: DO k=1,m1(d)
           sc  = xs(d, L-(k-1)*tau1(d))
           call addStrings(xL, sc)
           call addStrings(xyL,sc)
  end do x_Loop1
  end do xdim_Loop1
!  
  ydim_Loop1: DO d=1, Ndim
  y_Loop1: DO k=1,m2(d)
           sc  = ys(d, L-(k-1)*tau2(d))
           call addStrings(yL, sc)
           call addStrings(xyL,sc)
  end do y_Loop1
  end do ydim_Loop1
!  
  call Linklist(headx,tailx,xL)   
  call LinkList(heady,taily,yL)
  call Linklist(headxy,tailxy,xyL)

end do time_Loop
!
!!! --- calculate the Shannon entropy
NormX  = toBits / real(MX,  sifm_real)
NormY  = toBits / real(MY,  sifm_real)
NormXY = toBits / real(MXY, sifm_real)
Hx = SymbolicShannonEntropy(Headx) * NormX
call freeList(headx)
Hy = SymbolicShannonEntropy(Heady) * NormY
call freeList(heady)
Hxy = SymbolicShannonEntropy(Headxy) * NormXY
call freeList(headxy)
!
MI = Hx + Hy - Hxy
Write(*,'("entropies:   ", 6F12.5)') Hx, Hy, Hxy, hx+hy - hxy, Hxy - Hy, Hxy - Hx
IF (MI < rzero) MI = rzero
!
nullify(headx,tailx,heady,taily)
nullify(headxy,tailxy)
!
end subroutine symbolic_MI_ND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_MI2_ND(ndata,ndim,xs,ys,m1,m2,tau1,tau2,MI)
implicit none
!
integer, intent(in)          :: ndata,ndim
integer, dimension(ndim)     :: m1,m2,tau1,tau2
real(sifm_real), intent(in)  :: xs(:,:),ys(:,:)
real(sifm_real), intent(out) :: MI
!
type (StateElem), pointer    :: headx,tailx,heady,taily
type (StateElem), pointer    :: headxy,tailxy
integer                      :: d,istat,L,k,m,mx,my,mxy,L0,tau
character(len=1)             :: sc
character(len=250)           :: xL,yL
character(len=250)           :: xyL
real(sifm_real)              :: hx, hy, hxy, NormX, NormY, NormXY
integer, allocatable         :: XSL(:),YSL(:)
real(sifm_real), allocatable :: X(:),Y(:)
!
nullify(headx,tailx,heady,taily)
nullify(headxy,tailxy)
!
mx  = sum(m1)  !imax1d(m1)
my  = sum(m2)  !imax1d(m2)
mxy = Mx + My  !max(mx, my)
M   = max(imax1d(m1), imax1d(m2)) 
tau = max(imax1d(tau1), imax1d(tau2))
L0  = (M - 1) * Tau + 1
tau = max(imax1d(tau1),imax1d(tau2))
!
!!! --- Loop over all data points
time_Loop: DO L=L0, Ndata
  xL=""; yL=""; xyL=""
!  
  xdim_Loop: DO d=1,Ndim 
     Allocate( XSL(M1(d)) )
     Allocate( X(M1(d)) )
     X = vcopy(XS(d,:),L,M1(d),Tau1(d),-1)
     CALL sort(M1(d), X, XSL)
     x_Loop: DO k=1,m1(d)
           sc  = toString(xsL(k))
           call addStrings(xL,sc)
           call addStrings(xyL,sc)
     end do x_Loop
     Deallocate(XSL)
     Deallocate(X)
  end do xdim_Loop
!  
  ydim_Loop: DO d=1, Ndim
     Allocate( YSL(m2(d)) )
     Allocate( Y(m2(d)) )
     Y = vcopy(YS(d,:),L,M2(d),Tau2(d),-1)
     CALL sort(M2(d), Y, YSL)
     y_Loop: DO k=1,m2(d)
           sc  = toString(ysL(k))
           call addStrings(yL, sc)
           call addStrings(xyL, sc)
     end do y_Loop
     Deallocate( YSL )
     Deallocate( Y )
  end do ydim_Loop
! 
  call Linklist(headx,tailx,xL)   
  call LinkList(heady,taily,yL)
  call Linklist(headxy,tailxy,xyL) 
!
end do time_Loop
!
!!! --- calculate the Shannon entropy
NormX  = toBits / real(Mx,  sifm_real)
NormY  = toBits / real(My,  sifm_real)
NormXY = toBits / real(Mxy, sifm_real)
Hx = SymbolicShannonEntropy(Headx) * NormX
call freeList(headx)
Hy = SymbolicShannonEntropy(Heady) * NormY
call freeList(heady)
Hxy = SymbolicShannonEntropy(Headxy) * NormXY
call freeList(headxy)
!
MI = Hx + Hy - Hxy
!
Write(*,'(4F12.5)') Hx, Hy, Hxy, Hxy-Hx
!
nullify(headx,tailx,heady,taily)
nullify(headxy,tailxy)
!
end subroutine symbolic_MI2_ND
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function symbolicEntropyRate1D(ndata,xs,m,tau) result(Hrate)
implicit none
!
integer, intent(in)           :: ndata
integer                       :: m,tau
integer, intent(in)           :: xs(:)
real(sifm_real)                :: Hrate
!
type (StateElem), pointer  :: headx,tailx,heady,taily
integer                    :: L,k
character(len=1)           :: sc
character(len=250)         :: xL,yL
real(sifm_real)             :: hx,hy
!
nullify(headx,tailx,heady,taily)
!
!!! --- Loop over all data points
time_Loop: DO L=((M-1)*Tau+1), Ndata
xL=""; yL=""
!
x_Loop: DO k=1,m
  sc  = toString(xs(L-(k-1)*tau))
  call addStrings(xL, sc)
end do x_Loop
!
y_Loop: DO k=1,(m+1)
  sc  = toString(xs(L-(k-1)*tau))
  call addStrings(yL, sc)
end do y_Loop
!
call Linklist(headx,tailx,xL)
call LinkList(heady,taily,yL)
!
end do time_Loop
!!! --- calculate the Shannon entropy
Hx = SymbolicShannonEntropy(Headx)/real(m,sifm_real)/log(2.0_sifm_real)
call FreeList(headx)
Hy = SymbolicShannonEntropy(Heady)/real(m,sifm_real)/log(2.0_sifm_real)
call FreeList(heady)
!
Hrate = Hy - Hx
!
nullify(headx,tailx,heady,taily)
!
end function symbolicEntropyRate1D
!
end module MI_class

