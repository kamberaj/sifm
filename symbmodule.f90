!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: symbmodule.f90,v 1.0 19-03-2018, IBU
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
module SYMB_class
use sifm_kinds
use TEutils_class
use TErandom_class
implicit none
!
      integer, save :: Natoms, Nframes
      integer, save :: Ndim
      character(len=1), save, allocatable :: xyzs(:,:,:)
      real(sifm_real), save, allocatable :: xyz(:,:,:)
      integer, save, allocatable :: mopt(:,:), topt(:,:)
!
CONTAINS
!
subroutine SYMB_DRIVER(Nf, Nd, Na, debug, qSymbolic, Nmc)
  implicit none
  integer, intent(in) :: Nf, Nd, Na, Nmc
  Integer, intent(in) :: debug, qSymbolic
!  
  integer :: i, j, d, t, ierr, ndof
  real(sifm_real) :: Time_Start, Time_end
  real(sifm_real), dimension(:,:), allocatable :: X
  character(len=1), dimension(:,:), allocatable :: XS
  integer, dimension(:), allocatable :: Tau, M
!
  Nframes = Nf
  Ndim = Nd
  Natoms = Na
  Ndof = Nd * Na
!
  ALLOCATE( X(Ndof, Nframes), stat=ierr )
  if (ierr /= 0) stop 'Error allocating X'
  ALLOCATE( Xs(Ndof, Nframes), stat=ierr )
  if (ierr /= 0) stop 'Error allocating Xs'
  ALLOCATE( Tau(Ndof), stat=ierr )
  if (ierr /= 0) stop 'Error allocating Tau'
  ALLOCATE( M(Ndof), stat=ierr )
  if (ierr /= 0) stop 'Error allocating M'
!
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_Symbolics() 
!
  call read_xyz()
  call read_embdparam()
!
  J = 0
  DO i = 1, Natoms
     DO d = 1, Ndim
        J = J + 1
        X(J,:) = xyz(i,d,:)
        Tau(j) = Topt(i,d)
        M(j)   = Mopt(i,d)
     ENDDO
  ENDDO
!
  CALL CPU_time(Time_start)
  write(*,*) "Symbolic method:   ", qSymbolic
  IF (qsymbolic == 3) THEN
      CALL symbolize_trajectoryMC(Ndof, Nframes, X, XS, tau, M, Nmc)
  ELSEIF (qSymbolic == 2) THEN
      CALL symbolize_trajectory(Natoms, Ndim, Nframes, Xyz, XYZS, topt, mopt)
  ELSEIF (qsymbolic == 1) THEN
      CALL symbolize(Natoms, Ndim, Nframes, Xyz, XYZS)
  ELSE
      STOP 'Method has not been emplemented yet!'  
  ENDIF
!
  CALL CPU_time( Time_End )
  write(*,'("CPU Elapsed time for SYMB Module:   ", F10.6)') Time_End - Time_Start
!
  IF (qsymbolic == 3) THEN
  J = 0
  DO i = 1, Natoms
     DO d = 1, Ndim
        J = J + 1
        xyzs(i,d,:) = XS(J,:)
     ENDDO
  ENDDO
  ENDIF
!
  Call write_xyzs()    
!
  CALL DEALLOCATE_Symbolics()
!
  DEALLOCATE( X, stat=ierr )
  if (ierr /= 0) stop 'Error deallocating X'
  DEALLOCATE( Xs, stat=ierr )
  if (ierr /= 0) stop 'Error deallocating Xs'
  DEALLOCATE( Tau, stat=ierr )
  if (ierr /= 0) stop 'Error deallocating Tau'
  DEALLOCATE( M, stat=ierr )
  if (ierr /= 0) stop 'Error deallocating M'

!
  return
end subroutine SYMB_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the trajectories
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
!!  Print out the symbolic Trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_xyzs()
Implicit None
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='traj.xyzs', status='unknown', action='write')
write(nunit, '("SIFM: XYZS Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   DO i = 1, Natoms
      write(nunit, '(A7, 2I10, 3A3)') "FRAME  ", t, i, (xyzs(i,d,t), d=1, Ndim) 
   ENDDO
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
!
return
!
end subroutine write_xyzs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Allocate_Symbolics()
      implicit none
!
      integer :: ierr
!
!
      IF (.not. Allocated(xyzs)) THEN
          Allocate(XYZs(Natoms, ndim, Nframes), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating XYZs"
      ENDIF
      IF (.not. Allocated(xyz)) THEN
          Allocate(XYZ(Natoms, ndim, Nframes), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating XYZ"
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
end subroutine Allocate_Symbolics
!
subroutine deAllocate_Symbolics()
      implicit none
!
      integer :: ierr
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
end subroutine deAllocate_Symbolics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Symbolise trajectory using Method 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolize_trajectoryMC(Ndof, nframes, X, global_XS, tau, mopt, Nmc)
use TEMPI_CLASS
implicit none
!
integer, intent(in) :: Ndof, Nframes, Nmc
integer, dimension(:), intent(in) :: tau, mopt
real(sifm_real), dimension(:,:), intent(in) :: X
character(len=1), dimension(:,:), intent(inout) :: global_Xs
!
integer :: Ndata
integer :: i, ibin, nbin, ii, ierr, t, Istep
real(sifm_real) :: delta, Local_H, x_max, x_min
character(len=1), dimension(nframes) :: XS_temp
real(sifm_real), allocatable, dimension(:) :: Xc
character(len=1), allocatable, dimension(:) :: Local_Xs
real(sifm_real) :: time_start, time_end
!
real(sifm_real) :: rnd
type(val_t) :: Global_H_Max, Local_H_Max
integer :: Local_Nbin, Global_Nbin
!
CALL MPI_start()
!
call MPI_BCAST (Mopt,    Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Tau,     Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (NMC,     1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Ndof,    1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Nframes, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (X,       Ndof*NFrames, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
!
allocate( Local_Xs(Nframes), stat=ierr )
if (ierr /= 0) STOP "Error allocating Local-Xs"
!
if (myID == Master) then
    write(*,*) "Optimize nbin and Critical points:"
    CALL CPU_time( Time_start )
endif
!
DOF_LOOP: DO i = 1, Ndof
      IF (myID == Master) &
          write(*,'("Collect Info for IDIM:  ", I5)') I
      !Ndata    = Nframes - (Mopt(i) - 1) * Tau(i)
      !Nbin_min = floor( (one*Ndata)**(one/Mopt(i)) ) + 1
      Nbin = Mopt(i)
      IF (Nbin < 2) Nbin = 2
      IF (Nbin > 26) Nbin = 26

      IF (myID == Master) &
          write(*,'("Nbin:  ", 2I10)') Nbin
      X_min = minval( X(i, :) )
      X_max = maxval( X(i, :) )

! -- Initialize a set of random numbers
      CALL rndnum_iniall(1, 1, MyID)
!
      Local_H_Max%Value = rzero
      Local_Nbin = 0
!
      allocate( Xc(nbin+1), stat=ierr )
      IF (ierr /= 0) STOP "Error allocating Xc"
!
! -- Set XC(1) = min(X) and XC(nbin+1) = max(X)
      Xc(1) = X_min
      Xc(nbin+1) = X_max
      MC_LOOP: Do istep = myID + 1, Nmc, NumProcs
!
         DO ibin=2, Nbin
            Xc(ibin) = ( Xc(ibin-1) - Xc(Nbin+1) ) * randf(ii) + Xc(Nbin+1)
         ENDDO
!
         TIME_LOOP: Do t=1, nframes
            XS_temp(t) = ''
            Do ibin = 1, Nbin
               IF ( X(i,t)>=(Xc(ibin)-PRTINY)  .and.  X(i,t)<=(Xc(ibin+1)+PRTINY) ) THEN
                    XS_temp(t) = Symbol(ibin)
                    goto 9999
               ENDIF
            ENDDO
            IF ( XS_temp(t) == '' ) THEN
                 write(*,*) t, i, XS_temp(t)
                 stop 'Symbolization: Empty symbol'
            ENDIF
            9999 CONTINUE
         ENDDO TIME_LOOP
!
         CALL symbolic_entropy1D(Nframes, XS_temp, mopt(i), tau(i), Local_H)
!
         IF ( Local_H > Local_H_Max%Value) THEN
              Local_H_Max%Value = Local_H
              Local_H_Max%Rank  = myID
              Local_XS = XS_temp
              Local_Nbin  = Nbin
         ENDIF
!
      ENDDO MC_LOOP
!
      deallocate( Xc, stat=ierr )
      IF (ierr /= 0) STOP "Error deallocating Xc in Symbolize Traj using MC method"
!
      CALL MPI_Reduce( Local_H_Max, Global_H_Max, 1, MPI_2DOUBLE_PRECISION,    & 
                       MPI_MAXLOC, master, MPI_COMM_WORLD, ierr )
      IF (myID == Master) write(*,*) Global_H_Max%value, Global_H_Max%rank
      IF ( myID /= master ) THEN ! Send Info to Master Processor
           call MPI_Send( Local_XS, Nframes, MPI_CHARACTER, master, &
                          msgtag2, MPI_COMM_WORLD, ierr )
           call MPI_Send( Local_Nbin, 1, MPI_INTEGER, master, &
                          msgtag2, MPI_COMM_WORLD, ierr )      
      ELSE
           IF ( Global_H_max%rank /= master) THEN
                call MPI_Recv( Local_XS, Nframes, MPI_CHARACTER, Global_H_max%rank, msgtag2, &
                               MPI_COMM_WORLD, status, ierr )
                call MPI_Recv( Local_Nbin, 1, MPI_INTEGER, Global_H_max%rank, msgtag2, &
                               MPI_COMM_WORLD, status, ierr )
                Global_Nbin = Local_Nbin
                Global_XS(i,:) = Local_XS(:)
           ELSE
                Global_Nbin = Local_Nbin
                Global_XS(i,:) = Local_XS(:)
           ENDIF
      ENDIF
      if (myID == Master) &
          Write(*,'("Idim=   # of bins along the direction:  ", I5, 2X, I10)') i, global_Nbin
!
      CALL rndnum_clear()
!
ENDDO DOF_LOOP
!
call MPI_BCAST (Global_XS, Ndof*NFrames, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr )
!
deallocate( Local_Xs, stat=ierr )
if (ierr /= 0) STOP "Error deallocating Local-Xs"
!
if (myID == Master) then
    CALL CPU_time( Time_End )
    write(*,'("CPU Elapsed time for Symbolize Traj - Monte Carlo:   ", F10.6)') Time_End - Time_Start
endif
!
call MPI_finish()
!
RETURN
!
end subroutine symbolize_trajectoryMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Symbolise trajectory using Method 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolize_trajectory(Natoms, Ndim, nframes, X, XS, topt, mopt)
implicit none
integer, intent(in) :: Natoms, Ndim, Nframes
integer, dimension(:,:), intent(in) :: topt, mopt
real(sifm_real), dimension(:,:,:), intent(in) :: x
character(len=1), dimension(:,:,:), intent(out) :: xs
!
integer :: ndata
integer :: i, ia, ibin, nbin, ierr, t
real(sifm_real) :: delta, h, x_max, x_min
real(sifm_real), allocatable, dimension(:) :: xc
!
write(*,*) "Optimize nbin and Critical points:"
!
ATOMS_LOOP: DO ia = 1, Natoms
DIM_LOOP: do i=1, ndim
!
  !Ndata    = Nframes - (Mopt(ia,i) - 1) * Topt(ia,i)
  !Nbin_min = floor( (one*Ndata)**(one/Mopt(ia,i)) ) + 1
  Nbin = Mopt(ia,i)
  IF (Nbin < 2) Nbin = 2
  IF (Nbin > 26) Nbin = 26
  write(*,'("Nbin:  ", I10)') Nbin
!
  X_min = minval( X(ia,i,:) )
  X_max = maxval( X(ia,i,:) )
!
  allocate( xc(nbin+1), stat=ierr )
  if (ierr /= 0) STOP "Error allocating Xc"
!
! -- Set XC(1) = min(X) and XC(nbin+1) = max(X)
  Delta = (X_max - X_min)/real(Nbin,sifm_real)
  Xc(1) = X_min
  Xc(nbin+1) = X_max
  DO ibin = 1, Nbin-1
     Xc(ibin+1) = Xc(ibin) + Delta
  ENDDO
  DO t=1, nframes
     XS(ia,i,t) = ''
     DO ibin = 1, Nbin
        IF ( X(ia,i,t) >= (Xc(ibin)-prtiny)  .and.  X(ia,i,t) <= (Xc(ibin+1)+prtiny) ) THEN
             XS(ia,i,t) = Symbol(ibin)
             GOTO 9999
        ENDIF
     ENDDO
     IF ( XS(ia,i,t) == '' ) THEN
          write(*,*) t, ia, i, XS(ia,i,t)
          stop 'Symbolization: Empty symbol'
     ENDIF
     9999 continue
  ENDDO
!
  CALL symbolic_entropy1D(nframes, XS(ia,i,:), mopt(ia,i), topt(ia,i), h)
!
  deallocate( xc, stat=ierr )
  IF (ierr /= 0) STOP "Error deallocating xc"
!
  Write(*,'("Idim=    # of bins along the direction and Max entropy:  ", 3I5, F12.6)') i, Nbin, mopt(ia,i), H
!
enddo DIM_LOOP
enddo ATOMS_LOOP
!
end subroutine symbolize_trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Symbolise trajectory using Method 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolize(Natoms, Ndim, nframes, X, XS)
implicit none
integer, intent(in) :: Natoms, Ndim, Nframes
real(sifm_real), dimension(:,:,:), intent(in) :: x
character(len=1), dimension(:,:,:), intent(out) :: xs
!
integer :: i, ia, t
!
write(*,*) "Optimize nbin and Critical points:"
!
ATOMS_LOOP: DO ia = 1, Natoms
    DIM_LOOP: do i=1, ndim
      TIME_LOOP: DO t = 1, Nframes
         IF (X(ia,i,t) > rzero) THEN
             XS(ia,i,t) = Symbol(1)
         ELSE
             XS(ia,i,t) = Symbol(2)
         ENDIF
      ENDDO TIME_LOOP
    enddo DIM_LOOP
enddo ATOMS_LOOP
!
end subroutine symbolize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolic_entropy1D(ndata,xs,m,tau,H)
use telinklist_class
implicit none
!
integer, intent(in) :: ndata
integer,  intent(in) :: m,tau
character(len=1), dimension(:), intent(in) :: xs
real(sifm_real), intent(out) :: H
!
type (StateElem), pointer :: head,tail
integer                     :: L, L0, k
character(len=250)          :: xL
character(len=1)            :: sc
!
nullify(head,tail)
!
L0 = (M - 1) * Tau + 1
time_Loop: DO L = L0, ndata
xL=""
DO k = 1, M 
   call addStrings( XL, XS(L-(k-1)*Tau) )
ENDDO
call LinkList(head,tail,xL)
end do time_Loop
!
H = SymbolicShannonEntropy(Head)
call freeList(head)
!
nullify(head,tail)
!
return
!
end subroutine symbolic_entropy1D
!
end module Symb_class
