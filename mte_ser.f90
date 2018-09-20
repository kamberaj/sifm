!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: mte_ser.f90,v 1.0 19-03-2018, IBU
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
      real(sifm_real), save, dimension(:,:), allocatable :: TE, TES
      character(len=1), save, allocatable :: xyzs(:,:,:)
      real(sifm_real), save, allocatable :: xyz(:,:,:)
      integer, save, allocatable :: mopt(:,:), topt(:,:)

!
CONTAINS
!
subroutine STE_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, debug, &
                      Nshuffles, r0, statP)
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional   :: r0, statP
!  
      real(sifm_real)      :: Time_Start, Time_end
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
  CALL CPU_time(Time_start)
  CALL getSTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  CALL CPU_time( Time_End )
  write(*,'("CPU Elapsed time for TE:   ", F20.3)') Time_End - Time_Start
!
  call WRITE_TE(qTeShuffle)
! 
  CALL DEALLOCATE_TransferEntropy()
!
  return
end subroutine STE_DRIVER
!
subroutine TE_DRIVER(Nf, Nd, Na, qteShuffle, qteNorm, debug, &
                      Nshuffles, r0, statP)
  implicit none
      integer, intent(in)        :: qTeNorm
      integer, intent(in)        :: qTeShuffle
      integer, intent(in)        :: Na, Nf
      integer, intent(in)        :: Nd, debug
      integer, optional          :: Nshuffles
      real(sifm_real), optional   :: r0, statP
!  
  integer             :: i, j, d, t, ia, ja
  real(sifm_real)      :: Time_Start, Time_end
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
  CALL Allocate_TransferEntropy(qTeShuffle, 2) 
!
  CALL read_embdparam()
  CALL read_xyz()
!
  CALL CPU_time(Time_start)
  CALL getTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statp)
  CALL CPU_time( Time_End )
  write(*,'("CPU Elapsed time for TE:   ", F20.3)') Time_End - Time_Start
!
  call WRITE_TE(qTeShuffle)
!  
  CALL DEALLOCATE_TransferEntropy()
!
  return
end subroutine te_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the transfer entropy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_TE(qTeShuffle)
Implicit None
integer, intent(in) :: qTeShuffle
!  
integer             :: i, j, nunit
!  
nunit = new_unit()
open(unit=nunit, file='te.txt', status='unknown', action='write')
DO i=1, Natoms
   write(nunit,*) (TE(i,j), j=1, natoms)
ENDDO
!
IF (qteshuffle > 0) THEN
nunit = new_unit()
open(unit=nunit, file='ste.txt', status='unknown', action='write')
DO i=1, Natoms
   write(nunit,*) (TES(i,j), j=1, natoms)
ENDDO
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
	  ALLOCATE( TE(Natoms, Natoms), stat=ierr )
          if (ierr /= 0) stop 'Error allocating TE'
      ENDIF
!
      IF (qTEShuffle > 0) THEN
          IF (.not. Allocated(TES)) THEN
	      ALLOCATE( TES(Natoms, Natoms), stat=ierr )
              if (ierr /= 0) stop 'Error allocating TES'
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
subroutine getSTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
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
integer :: i, j, itrial,t, Ilist
character(len=1), dimension(ndim,nframes) :: xs_shuffle
real(sifm_real) :: txy, tyx, avr_Txy, avr_Tyx, hxx, hyy, hx, hy, &
                  time_start, time_end
!
integer :: ierr
real(sifm_real) :: gtxy,gtyx
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
endif

! --- Compute the entropy transfer
TE = rzero; TES = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      call cpu_time(time_start)
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
      TE(i,j) = Txy; TE(j,i) = Tyx
      write(*,'("--------------------------------------------------------")')
      write(*,'("Instantaneous Txy (+bias)  =        ", F12.6)') Txy
      write(*,'("Instantaneous Tyx (+bias)  =        ", F12.6)') Tyx
!---- Remove the bias, by shuffling the trajectories
      IF (qteshuffle > 0) THEN

       ! -- Subtract the bias from T(x->y) and normalize it
       gTxy=0.0_sifm_real
       xs_shuffle = XYZS(i,:,:)
       DO itrial = 1, Ntrials
          IF (qTEShuffle == 1) THEN
              CALL shuffleND(nframes,ndim,xs_shuffle,I)
          ELSE
              CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(i,:), Topt(i,:), I)
          ENDIF
          CALL symbolic_TE_entropyND(nframes,ndim,xs_shuffle,XYZs(j,:,:),          &
		                     Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),gTxy)
       ENDDO
       Txy = gTxy / real(Ntrials, sifm_real)
       write(*,'("--------------------------------------------------------")')
       write(*,'("Average value of Bias Txy (Before Statistical Test)  =  ", F12.6)') Txy
       CALL statTest(Ntrials, R0, Txy, Txy, StatP)	
       write(*,'("Average value of Bias Txy (After Statistical Test)   =  ", F12.6)') Txy		 
       IF (qTENorm == 0) THEN
	   TES(i,j) = (TE(i,j) - Txy)
       ELSE
           IF ( abs(Hxx) > prtiny ) THEN
   	        TES(i,j) = (TE(i,j) - Txy)/Hxx
           ELSE
   	        TES(i,j) = (TE(i,j) - Txy)/(Hx+Hy)
           ENDIF
       ENDIF
       write(*,'("--------------------------------------------------------")')
       write(*,'("Average value of true TSxy (Before Statistical Test) =  ", F12.6)') TES(i,j)
       CALL statTest(Ntrials, R0, TES(i,j), TES(i,j), StatP)	
       write(*,'("Average value of true TSxy (After Statistical Test)  =  ", F12.6)') TES(i,j)
       write(*,'("--------------------------------------------------------")')
       IF ( debug == 1 ) THEN
            write(*,'("Txy (+bias, -bias) =  ", 5F12.6)') TE(i,j), TES(i,j), Hxx,Hx,Hy
       Endif
 
       ! -- Subtract the bias from T(y->x) and normalize it
       gTyx = 0.0_sifm_real
       xs_shuffle = XYZs(j,:,:)
       DO itrial = 1, Ntrials
           IF (qTEShuffle == 1) THEN
              CALL shuffleND(nframes,ndim,xs_shuffle,Natoms+J)
          ELSE
              CALL BlockShuffleND(nframes, ndim, xs_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
          ENDIF
          CALL symbolic_TE_entropyND(nframes,ndim,xs_shuffle, XYZs(i,:,:), &
                                          Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),gTyx)
       ENDDO	  
       Tyx = gTyx / real(Ntrials, sifm_real)
       write(*,'("--------------------------------------------------------")')
       write(*,'("Average value of Bias Tyx (Before Statistical Test)  =  ", F12.6)') Tyx			 
       CALL statTest(Ntrials, R0, TYx, Tyx, STATP)
       write(*,'("Average value of Bias Tyx (After Statistical Test)   =  ", F12.6)') Tyx
       write(*,'("--------------------------------------------------------")')			 
       IF (qTeNorm == 0) THEN
	   TES(j,i) = (TE(j,i) - Tyx)
       ELSE 
           IF ( abs(Hyy) > prtiny ) THEN
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
!
      CALL CPU_time( Time_End )
      write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F20.3)') Time_End - Time_Start
!
enddo Pairs_LOOP
!
IF ( qteshuffle > 0) then
     CALL rndnum_clear()
ENDIF
!
IF ( Allocated(listPairs) ) THEN
     DEALLOCATE( listPairs, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Pairs'
ENDIF
!
return
end subroutine getSTransferEntropyNDIM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getTransferEntropyNDIM(qteShuffle, qteNorm, debug, Nshuffles, r0, statP)
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
integer :: i, j, itrial,t, Ilist
real(sifm_real), dimension(ndim,nframes) :: x_shuffle
real(sifm_real) :: txy, tyx, avr_Txy, avr_Tyx, hxx, hyy, hx, hy, &
                  time_start, time_end
!
integer :: ierr
real(sifm_real) :: gtxy,gtyx

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
endif

! --- Compute the entropy transfer
TE = rzero; TES = rzero
!-------- call all the pairs
Pairs_LOOP: DO iList = 1, Npairs
!
      call cpu_time(time_start)
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
      TE(i,j) = Txy; TE(j,i) = Tyx
      write(*,'("--------------------------------------------------------")')
      write(*,'("Instantaneous Txy (+bias)  =        ", F12.6)') Txy
      write(*,'("Instantaneous Tyx (+bias)  =        ", F12.6)') Tyx
!---- Remove the bias, by shuffling the trajectories
      IF (qteshuffle > 0) THEN

       ! -- Subtract the bias from T(x->y) and normalize it
       gTxy=0.0_sifm_real
       x_shuffle = XYZ(i,:,:)
       DO itrial=1, Ntrials
          IF (qTEShuffle == 1) THEN
              CALL shuffleND2(nframes,ndim,x_shuffle,I)
          ELSE
              CALL BlockShuffleND2(nframes, ndim, x_shuffle, Mopt(i,:), Topt(i,:), I)
          ENDIF
          CALL symbolic_TE2_entropyND(nframes,ndim,x_shuffle,XYZ(j,:,:),          &
		                      Mopt(i,:),Mopt(j,:),Topt(i,:),Topt(j,:),gTxy)
       ENDDO
       Txy = gTxy / real(Ntrials, sifm_real)
       write(*,'("--------------------------------------------------------")')
       write(*,'("Average value of Bias Txy (Before Statistical Test)  =  ", F12.6)') Txy
       CALL statTest(Ntrials, R0, Txy, Txy, StatP)	
       write(*,'("Average value of Bias Txy (After Statistical Test)   =  ", F12.6)') Txy		 
       IF (qTENorm == 0) THEN
	   TES(i,j) = (TE(i,j) - Txy)
       ELSE
           IF ( abs(Hxx) > prtiny ) THEN
   	        TES(i,j) = (TE(i,j) - Txy)/Hxx
           ELSE
   	        TES(i,j) = (TE(i,j) - Txy)/(Hx+Hy)
           ENDIF
       ENDIF
       write(*,'("--------------------------------------------------------")')
       write(*,'("Average value of true TSxy (Before Statistical Test) =  ", F12.6)') TES(i,j)
       CALL statTest(Ntrials, R0, TES(i,j), TES(i,j), StatP)	
       write(*,'("Average value of true TSxy (After Statistical Test)  =  ", F12.6)') TES(i,j)
       write(*,'("--------------------------------------------------------")')
       IF ( debug == 1 ) THEN
            write(*,'("Txy (+bias, -bias) =  ", 5F12.6)') TE(i,j), TES(i,j), Hxx,Hx,Hy
       Endif

       ! -- Subtract the bias from T(y->x) and normalize it
       gTyx = 0.0_sifm_real
       x_shuffle = XYZ(j,:,:)
       DO itrial = 1, Ntrials
          IF (qTEShuffle == 1) THEN
              CALL shuffleND2(nframes,ndim,x_shuffle,Natoms+J)
          ELSE
              CALL BlockShuffleND2(nframes, ndim, x_shuffle, Mopt(j,:), Topt(j,:), Natoms+J)
          ENDIF
          CALL symbolic_TE2_entropyND(nframes,ndim,x_shuffle, XYZ(i,:,:), &
                                      Mopt(j,:),Mopt(i,:),Topt(j,:),Topt(i,:),gTyx)
       ENDDO	  
       Tyx = gTyx / real(Ntrials, sifm_real)
       write(*,'("--------------------------------------------------------")')
       write(*,'("Average value of Bias Tyx (Before Statistical Test)  =  ", F12.6)') Tyx			 
       CALL statTest(Ntrials, R0, TYx, Tyx, STATP)
       write(*,'("Average value of Bias Tyx (After Statistical Test)   =  ", F12.6)') Tyx
       write(*,'("--------------------------------------------------------")')			 
       IF (qTeNorm == 0) THEN
	   TES(j,i) = (TE(j,i) - Tyx)
       ELSE 
           IF ( abs(Hyy) > prtiny ) THEN
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
!
      CALL CPU_time( Time_End )
      write(*,'("CPU Elapsed time for computation of Txy and Tyx:   ", F20.3)') Time_End - Time_Start
!
enddo Pairs_LOOP
!
IF ( qteshuffle > 0) then
     CALL rndnum_clear()
ENDIF
!
IF ( Allocated(listPairs) ) THEN
     DEALLOCATE( listPairs, stat=ierr )
     if (ierr /= 0) stop 'Error deallocating List of Pairs'
ENDIF
!
return
end subroutine getTransferEntropyNDIM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropies1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)            :: ndata,m1,m2,tau1,tau2
character(len=1), intent(in)   :: xs(:),ys(:)
real(sifm_real), intent(out)   :: Txy,Tyx
real(sifm_real), intent(out)   :: Hxx,Hyy,hx,hy
!
type (StateElem), pointer  :: headx,tailx,heady,taily,headx1,tailx1,heady1,taily1
type (StateElem), pointer  :: headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1
integer, parameter         :: Delta=1
integer                    :: istat,L,k,m,tau
character(len=1)           :: xsc,ysc
character(len=250)         :: xL,x1L,yL,y1L
character(len=250)         :: xyL, xy1L, yx1L
real(sifm_real)             :: hx1,hy1,hxy,hxy1,hyx1
!
nullify(headx,tailx,heady,taily,headx1,tailx1,heady1,taily1)
nullify(headxy,tailxy,headxy1,tailxy1,headyx1,tailyx1)
!
m=max(m1,m2);tau=max(tau1,tau2)
!!! --- Loop over all data points
write(*,*) "Here"
time_Loop: DO L=((M-1)*Tau+1),(Ndata-Delta)
  xL=""; yL=""
  x1L=""; y1L=""
  xyL=""; xy1L=""; yx1L=""
!
  x_Loop: DO k=1,m1
          xsc  = xs(L-(k-1)*tau1)
	  xL   = addStrings(xL,xsc)
          x1L  = addStrings(x1L,xsc)
          xyL  = addStrings(xyL,xsc)
   	  xy1L = addStrings(xy1L,xsc)
          yx1L = addStrings(yx1L,xsc)
  end do x_Loop
!
  y_Loop: DO k=1,m2
          ysc  = ys(L-(k-1)*tau2)
	  yL   = addStrings(yL,ysc)
  	  y1L  = addStrings(y1L,ysc)
	  xyL  = addStrings(xyL,ysc)
	  xy1L = addStrings(xy1L,ysc)
   	  yx1L = addStrings(yx1L,ysc)
  end do y_Loop
!
  ysc  = ys(L+Delta)
  y1L  = addStrings(y1L,ysc)
  xy1L = addStrings(xy1L, ysc)
!
  xsc  = xs(L+Delta)
  x1L  = addStrings(x1L,xsc)
  yx1L = addStrings(yx1L, xsc)
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
end subroutine symbolic_TE_entropies1d
!
subroutine symbolic_TE_entropiesND2(ndata,ndim,xs,ys,m1,m2,tau1,tau2, &
                                    Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)           :: ndata,ndim
integer, dimension(ndim)      :: m1,m2,tau1,tau2
character(len=1), intent(in)  :: xs(:,:),ys(:,:)
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
	    xL   = addStrings(xL,sc)
	    x1L  = addStrings(x1L,sc)
	    xyL  = addStrings(xyL,sc)
	    xy1L = addStrings(xy1L,sc)
            yx1L = addStrings(yx1L, sc)
  end do x_Loop1
  end do xdim_Loop1

  xdim_Loop2: DO d=1,Ndim 
      sc   = xs(d,L+Delta)
      x1L  = addStrings(x1L,sc)
      yx1L = addStrings(yx1L, sc)
  end do xdim_Loop2
!  
  ydim_Loop1: DO d=1, Ndim
  y_Loop1: DO k=1,m2(d)
            sc   = ys(d,L-(k-1)*tau2(d))
	    yL   = addStrings(yL, sc)
	    y1L  = addStrings(y1L, sc)
	    xyL  = addStrings(xyL, sc)
	    xy1L = addStrings(xy1L,sc)
	    yx1L = addStrings(yx1L, sc)
  end do y_Loop1
  end do ydim_Loop1

  ydim_Loop2: DO d=1, Ndim
       sc   = ys(d,L+Delta)
       y1L  = addStrings(y1L,sc)
       xy1L = addStrings(xy1L, sc)
  end do ydim_Loop2
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
end subroutine symbolic_TE_entropiesND2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropiesND(ndata,ndim,xs,ys,m1,m2,tau1,tau2, &
                                    Txy,Tyx,Hxx,Hyy,hx,hy)
implicit none
!
integer, intent(in)            :: ndata,ndim
integer, dimension(ndim)       :: m1,m2,tau1,tau2
character(len=1), intent(in)   :: xs(:,:),ys(:,:)
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
            sc = xs(d, L-(k-1)*tau1(d))
	    xL = addStrings(xL,sc)
	    x1L = addStrings(x1L,sc)
	    xyL = addStrings(xyL,sc)
	    xy1L = addStrings(xy1L,sc)
            yx1L = addStrings(yx1L, sc)
      end do x_Loop1
      sc   = xs(d,L+Delta)
      x1L  = addStrings(x1L,sc)
      yx1L = addStrings(yx1L, sc)
  end do xdim_Loop1
!  
  ydim_Loop1: DO d=1, Ndim
       y_Loop1: DO k=1,m2(d)
            sc = ys(d,L-(k-1)*tau2(d))
	    yL = addStrings(yL, sc)
	    y1L = addStrings(y1L, sc)
	    xyL = addStrings(xyL, sc)
	    xy1L = addStrings(xy1L,sc)
	    yx1L = addStrings(yx1L, sc)
       end do y_Loop1
       sc   = ys(d,L+Delta)
       y1L  = addStrings(y1L,sc)
       xy1L = addStrings(xy1L, sc)
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
end subroutine symbolic_TE_entropiesND   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE_entropy1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)            :: ndata,m1,m2,tau1,tau2
character(len=1), intent(in)   :: xs(:),ys(:)
real(sifm_real), intent(inout) :: Txy
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
time_Loop: DO L=((M-1)*Tau+1),(Ndata-Delta)
!
   yL ="";  yy1L =""
   xyL="";  xyy1L=""
!
   x_Loop: DO k=1,M1
           sc    = xs(L-(k-1)*Tau1)
           xyL   = addStrings(xyL, sc)
	   xyy1L = addStrings(xyy1L, sc)
   end do x_Loop
!
   y_Loop: DO k=1,M2
             sc    = ys(L-(k-1)*Tau2)
	     yL    = addStrings(yL,sc)
	     xyL   = addStrings(xyL, sc)
	     yy1L  = addStrings(yy1L,sc)
             xyy1L = addStrings(xyy1L, sc)
   end do y_Loop
!
   sc    = ys(L+Delta)
   yy1L  = addStrings(yy1L,sc)
   xyy1L = addStrings(xyy1L, sc)
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
Txy=Txy + (hy1-hy) - (hxy1-hxy)
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
real(sifm_real), intent(inout) :: Txy
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
	    xL   = addStrings(xL,xsc)
          xyL  = addStrings(xyL,xsc)
   	    xy1L = addStrings(xy1L,xsc)
  end do x_Loop
!
  X1(1:M1) = X(1:M1)
  X1((m1+1):(m1+1)) = XS(L+Delta) 
  CALL sort(M1+1, X1, XSL1)
  x1_Loop: DO k=1, M1+1
            xsc  = ToString( XSL1(k) )
 	      x1L  = addStrings(x1L, xsc)
            yx1L = addStrings(yx1L, xsc)
  end do x1_Loop 
!
  Y = vcopy(YS,L,M2,Tau2,-1)
  CALL sort(M2, Y, YSL) 
  y_Loop: DO k=1,m2
          ysc  = toString( ysL(k) )
	    yL   = addStrings(yL, ysc)
	    xyL  = addStrings(xyL, ysc)
   	    yx1L = addStrings(yx1L, ysc)
  end do y_Loop
!
   Y1(1:M2) = Y(1:M2)
   Y1((m2+1):(m2+1)) = YS(L+Delta) 
   CALL sort(M2+1, Y1, YSL1)
   y1_Loop: DO k=1, M2+1
            ysc = ToString( YSL1(k) )
 	      y1L = addStrings(y1L, ysc)
            xy1L = addStrings(xy1L, ysc)
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
end subroutine symbolic_TE2_entropiesND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_TE2_entropy1D(ndata,xs,ys,m1,m2,tau1,tau2,Txy)
implicit none
!
integer, intent(in)           :: ndata,m1,m2,tau1,tau2
real(sifm_real), intent(in)    :: xs(:),ys(:)
real(sifm_real), intent(inout) :: Txy
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
time_Loop: DO L=(M-1)*Tau+1,(Ndata-Delta)
   yL=""; y1L=""
   xyL=""; xy1L=""
!
   X = vcopy(XS,L,M1,Tau1,-1)
   CALL sort(M1, X, XSL)
   x_Loop: DO k=1, M1
           sc = ToString( XSL(k) )
           xyL = addStrings(xyL, sc)
	     xy1L = addStrings(xy1L, sc)
   end do x_Loop
!
   Y = vcopy(YS,L,M2,Tau2,-1)
   CALL sort(M2, Y, YSL)
   y_Loop: DO k=1, M2
           sc = ToString( YSL(k) )
	     yL = addStrings(yL, sc)
	     xyL = addStrings(xyL, sc)
   end do y_Loop
!
   Y1(1:M2) = Y(1:M2)
   Y1((m2+1):(m2+1)) = YS(L+Delta) 
   CALL sort(M2+1, Y1, YSL1)
   y1_Loop: DO k=1, M2+1
            sc = ToString( YSL1(k) )
 	      y1L = addStrings(y1L, sc)
            xy1L = addStrings(xy1L, sc)
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
Txy=Txy + (hy1-hy) - (hxy1-hxy)
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
end module TE_class
