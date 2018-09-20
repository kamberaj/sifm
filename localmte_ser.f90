!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: localmte_ser.f90,v 1.0 19-03-2018, IBU
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
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, db
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!  
  real(sifm_real)      :: Time_Start, Time_end
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
  CALL Allocate_LocalTransferEntropy(qTeShuffle, 2) 
!
  CALL read_embdparam()
  CALL read_xyz()
!
  CALL CPU_time(Time_start)
  CALL getLocalTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  CALL CPU_time( Time_End )
  write(*,'("CPU Elapsed time for TE:   ", F10.6)') Time_End - Time_Start
!
  call WRITE_LTE(qTeShuffle)
!
  CALL DEALLOCATE_LocalTransferEntropy()
!
  STOP 'Program finished successfully'    
!
  return
end subroutine LTE_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LSTE_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, db, Nshuffles, r0, statP)
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, db
      integer, optional          :: Nshuffles
      real(sifm_real), optional  :: r0, statP
!  
  real(sifm_real)      :: Time_Start, Time_end
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
  CALL Allocate_LocalTransferEntropy(qTeShuffle, 1) 
!
  CALL read_embdparam()
  CALL read_xyzs()
!
  CALL CPU_time(Time_start)
  CALL getLocalSTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  CALL CPU_time( Time_End )
  write(*,'("CPU Elapsed time for TE:   ", F10.6)') Time_End - Time_Start
!
  call WRITE_LTE(qTeShuffle)
!
  CALL DEALLOCATE_LocalTransferEntropy()
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
          if (ierr /= 0) stop 'Error allocating local TE'
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

!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLocalSTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
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
integer :: i, j, itrial, t, ierr
character(len=1), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real), dimension(Nframes) :: Ltxy,Ltyx
real(sifm_real), dimension(:,:), allocatable :: Localtxy,Localtyx
real(sifm_real) :: time_start, time_end
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
! -- Initialize 2*natoms random number sequences
IF ( qteshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, 0)
     ntrials = nshuffles
     ALLOCATE( LocalTxy(Ntrials, Nframes), stat=ierr )
     if (ierr /= 0) stop 'Error allocating List of Local-Txy'
     ALLOCATE( LocalTyx(Ntrials, Nframes), stat=ierr )
     if (ierr /= 0) stop 'Error allocating List of Local-Tyx'
ENDIF

! --- Compute the entropy transfer
localTE  = rzero
IF (qTEShuffle > 0) THEN
    localTES = rzero
ENDIF
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)

      !--- Compute the direct transfer entropy
      CALL symbolic_LTE_entropies(nframes, ndim, XYZs(i,:,:), XYZs(j,:,:),                   &
                                      Mopt(i,:), Mopt(j,:), Topt(i,:), Topt(j,:),            &
                                      LTxy, LTyx)
      localTE(:,i,j) = LTxy
      localTE(:,j,i) = LTyx
      IF (qteshuffle > 0) THEN
           ! -- Subtract the bias from T(x->y) and normalize it
           LocalTxy = rzero        
           xs_shuffle = XYZS(i,:,:)
           DO itrial=1, Ntrials
               IF (qTEShuffle == 1) THEN
                   CALL shuffleND(nframes,ndim,xs_shuffle,I)
               ELSE
                   CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(i,:), Topt(i,:), I)
               ENDIF               
               CALL symbolic_LTE_entropy(nframes,ndim,xs_shuffle,XYZs(j,:,:),                &
                                         Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),            &
                                         LocalTxy(itrial,:))
           ENDDO
           LTxy = sum( LocalTxy, dim=1 ) / real(Ntrials, sifm_real)
           LocalTES(:,i,j) = LocalTE(:,i,j) - LTxy

      ! -- Subtract the bias from T(y->x) and normalize it
           LocalTyx = rzero
           xs_shuffle = XYZs(j,:,:)
           DO itrial=1, Ntrials
              IF (qTEShuffle == 1) THEN
                  CALL shuffleND(nframes,ndim,xs_shuffle,Natoms+J)
              ELSE
                  CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
              ENDIF
              CALL symbolic_LTE_entropy(nframes,ndim,xs_shuffle, XYZs(i,:,:),            &
                                        Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),         &
                                        LocalTyx(Itrial,:))
           ENDDO
           LTyx = sum( LocalTyx, dim=1 ) / real(Ntrials, sifm_real)
           LocalTES(:,j,i) = LocalTE(:,j,i) - LTyx
      ENDIF   
!
      CALL CPU_time( Time_End )
      write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F10.6)') Time_End - Time_Start
!
enddo Pairs_LOOP
!
IF ( qteshuffle > 0) then
     CALL rndnum_clear()
     DEALLOCATE( LocalTxy, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Local-Txy'
     DEALLOCATE( LocalTyx, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Local-Tyx'
ENDIF
!
IF ( Allocated(listPairs) ) THEN
     DEALLOCATE( listPairs, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Pairs'
ENDIF
!
return
end subroutine getLocalSTransferEntropyNDIM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLocalTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
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
integer :: i, j, itrial, t, ierr
real(sifm_real), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real), dimension(Nframes) :: Ltxy,Ltyx
real(sifm_real), dimension(:,:), allocatable :: Localtxy,Localtyx
real(sifm_real) :: time_start, time_end
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
! -- Initialize 2*natoms random number sequences
IF ( qteshuffle > 0 ) then
     CALL rndnum_iniall(2*Natoms, 1, 0)
     ntrials = nshuffles
     ALLOCATE( LocalTxy(Ntrials, Nframes), stat=ierr )
     if (ierr /= 0) stop 'Error allocating List of Local-Txy'
     ALLOCATE( LocalTyx(Ntrials, Nframes), stat=ierr )
     if (ierr /= 0) stop 'Error allocating List of Local-Tyx'
ENDIF

! --- Compute the entropy transfer
localTE = rzero
IF ( qteshuffle > 0 ) THEN
     localTES = rzero
ENDIF
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      call cpu_time(time_start)
!
      I = ListPairs(iList, 1)
      J = ListPairs(iList, 2)

      !--- Compute the direct transfer entropy
      CALL symbolic_LTE2_entropies(nframes, ndim, XYZ(i,:,:), XYZ(j,:,:),                   &
                                      Mopt(i,:), Mopt(j,:), Topt(i,:), Topt(j,:),            &
                                      LTxy, LTyx)
      LocalTE(:,i,j) = LTxy
      LocalTE(:,j,i) = LTyx
      IF (qteshuffle > 0) THEN
           ! -- Subtract the bias from T(x->y) and normalize it
           LocalTxy = rzero       
           xs_shuffle = XYZ(i,:,:)
           DO itrial=1, Ntrials
               IF (qTEShuffle == 1) THEN
                   CALL shuffleND2(nframes,ndim,xs_shuffle,I)
               ELSE
                   CALL BlockShuffleND2(nframes, ndim, xs_shuffle, Mopt(I,:), Topt(I,:), I)
               ENDIF
               CALL symbolic_LTE2_entropy(nframes,ndim,xs_shuffle,XYZ(j,:,:),                &
                                          Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),            &
                                          LocalTxy(itrial,:))
           ENDDO
           LTxy = sum( LocalTxy, dim=1 ) / real(Ntrials, sifm_real)
           LocalTES(:,i,j) = LocalTE(:,i,j) - LTxy

      ! -- Subtract the bias from T(y->x) and normalize it
           LocalTyx = rzero
           xs_shuffle = XYZ(j,:,:)
           DO itrial=1, Ntrials
              IF (qTEShuffle == 1) THEN
                  CALL shuffleND2(nframes,ndim,xs_shuffle,Natoms+J)
              ELSE
                  CALL BlockShuffleND2(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
              ENDIF
              CALL symbolic_LTE2_entropy(nframes,ndim,xs_shuffle, XYZ(i,:,:),            &
                                         Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),         &
                                         LocalTyx(Itrial,:))
           ENDDO
           LTyx = sum( LocalTyx, dim=1 ) / real(Ntrials, sifm_real)
           LocalTES(:,j,i) = LocalTE(:,j,i) - LTyx
      ENDIF   
!
      CALL CPU_time( Time_End )
      write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F10.6)') Time_End - Time_Start
!
enddo Pairs_LOOP
!
IF ( qteshuffle > 0) then
     CALL rndnum_clear()
     DEALLOCATE( LocalTxy, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Local-Txy'
     DEALLOCATE( LocalTyx, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Local-Tyx'
ENDIF
!
IF ( Allocated(listPairs) ) THEN
     DEALLOCATE( listPairs, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Pairs'
ENDIF
!
return
end subroutine getLocalTransferEntropyNDIM
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
          xL   = addStrings(xL,xsc)
          xyL  = addStrings(xyL,xsc)
          xy1L = addStrings(xy1L,xsc)
  end do x_Loop1
!
  Allocate( XSL1(M1(d)+1) )
  Allocate( X1(M1(d)+1) )
  X1(1:M1(d)) = X(1:M1(d))
  X1((M1(d)+1):(M1(d)+1)) = XS(d,L+Delta)
  CALL sort(M1(d)+1, X1, XSL1)
  x1_Loop: DO k=1, M1(d)+1
           xsc  = ToString( XSL1(k) )
           x1L = addStrings(x1L, xsc)
           yx1L= addStrings(yx1L, xsc)
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
         yL   = addStrings(yL,ysc)
         xyL  = addStrings(xyL,ysc)
         yx1L = addStrings(yx1L,ysc)
  end do y_Loop1
!
  Allocate( YSL1(m2(d)+1) )
  Allocate( Y1(m2(d)+1) )
  Y1(1:M2(d)) = Y(1:M2(d))
  Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta)
  CALL sort(M2(d)+1, Y1, YSL1)
  y1_Loop: DO k=1, M2(d)+1
           ysc = ToString( YSL1(k) )
           y1L = addStrings(y1L, ysc)
           xy1L = addStrings(xy1L, ysc)
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
          xyL  = addStrings(xyL,xsc)
          xy1L = addStrings(xy1L,xsc)
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
         yL   = addStrings(yL,ysc)
         xyL  = addStrings(xyL,ysc)
  end do y_Loop1
!
  Allocate( YSL1(m2(d)+1) )
  Allocate( Y1(m2(d)+1) )
  Y1(1:M2(d)) = Y(1:M2(d))
  Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta)
  CALL sort(M2(d)+1, Y1, YSL1)
  y1_Loop: DO k=1, M2(d)+1
           ysc = ToString( YSL1(k) )
           y1L = addStrings(y1L, ysc)
           xy1L = addStrings(xy1L, ysc)
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
        xL   = addStrings(xL,xsc)
        xyL  = addStrings(xyL,xsc)
        x1L  = addStrings(x1L,xsc)
        yx1L = addStrings(yx1L,xsc)
        xy1L = addStrings(xy1L,xsc)
      end do x_Loop1
      xsc  = xs(d,L+Delta)
      x1L  = addStrings(x1L,xsc)
      yx1L = addStrings(yx1L,xsc)
   end do xdim_loop
!
   ydim_loop: do d=1, ndim
      y_Loop1: DO k=1,m2(d)
        ysc  = ys(d,L-(k-1)*tau2(d))
        yL   = addStrings(yL,ysc)
        xyL  = addStrings(xyL,ysc)
        y1L  = addStrings(y1L,ysc)
        xy1L = addStrings(xy1L,ysc)
        yx1L = addStrings(yx1L,ysc)
      end do y_Loop1
      ysc  = ys(d,L+Delta)
      y1L  = addStrings(y1L,ysc)
      xy1L = addStrings(xy1L,ysc)
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
         xyL   = addStrings(xyL,xsc)
         xyy1L = addStrings(xyy1L,xsc)
     end do x_Loop1
  end do xdim_loop
!
  ydim_loop: DO d = 1, ndim
     y_Loop1: DO k=1,m2(d)
         ysc   = ys(d,L-(k-1)*tau2(d))
         yL    = addStrings(yL,ysc)
         xyL   = addStrings(xyL,ysc)
         yy1L  = addStrings(yy1L,ysc)
         xyy1L = addStrings(xyy1L,ysc)
     end do y_Loop1
     ysc   = ys(d,L+Delta)
     yy1L  = addStrings(yy1L,ysc)
     xyy1L = addStrings(xyy1L,ysc)
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
integer, intent(in)         :: ndata,m,tau
character(len=1), intent(in)         :: xs(ndata)
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
         xL = addStrings(xL,xsc)
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
xyL  = addStrings(xyL, sc)
xy1L = addStrings(xy1L, sc)
end do x_Loop
end do xdim_Loop
!
ydim_Loop: DO d=1, Ndim
y_Loop: DO k=1,m2(d)
sc  = ys(d,L-(k-1)*tau2(d))
yL  = addStrings(yL,sc)
xyL = addStrings(xyL, sc)
y1L = addStrings(y1L,sc)
xy1L= addStrings(xy1L, sc)
end do y_Loop
sc    = ys(d,L+Delta)
y1L   = addStrings(y1L,sc)
xy1L  = addStrings(xy1L, sc)
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
xyL = addStrings(xyL, sc)
xy1L = addStrings(xy1L, sc)
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
yL = addStrings(yL,sc)
xyL = addStrings(xyL, sc)
end do y_Loop
!
Allocate( YSL1(m2(d)+1), Y1(m2(d)+1) )
Y1(1:M2(d)) = Y(1:M2(d))
Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta)
CALL sort(M2(d)+1, Y1, YSL1)
y1_Loop: DO k=1, M2(d)+1
sc = toString( YSL1(k) )
y1L = addStrings(y1L,sc)
xy1L = addStrings(xy1L, sc)
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
xL = addStrings(xL,sc)
xyL = addStrings(xyL,sc)
xy1L = addStrings(xy1L,sc)
end do x_Loop
Allocate( XSL1(M1(d)+1) )
Allocate( X1(M1(d)+1) )
X1(1:M1(d)) = X(1:M1(d))
X1((M1(d)+1):(M1(d)+1)) = XS(d,L+Delta)
CALL sort(M1(d)+1, X1, XSL1)
x1_Loop: DO k=1, M1(d)+1
sc  = ToString( XSL1(k) )
x1L = addStrings(x1L, sc)
yx1L= addStrings(yx1L, sc)
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
yL = addStrings(yL, sc)
xyL = addStrings(xyL, sc)
yx1L = addStrings(yx1L, sc)
end do y_Loop
Allocate( YSL1(m2(d)+1) )
Allocate( Y1(m2(d)+1) )
Y1(1:M2(d)) = Y(1:M2(d))
Y1((m2(d)+1):(m2(d)+1)) = YS(d,L+Delta)
CALL sort(M2(d)+1, Y1, YSL1)
y1_Loop: DO k=1, M2(d)+1
sc = ToString( YSL1(k) )
y1L = addStrings(y1L, sc)
xy1L = addStrings(xy1L, sc)
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
