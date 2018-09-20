!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: symbmodule_ser.f90,v 1.0 19-03-2018, IBU
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
  write(*,'("CPU Elapsed time for TE:   ", F10.6)') Time_End - Time_Start
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
implicit none
!
integer, intent(in) :: Ndof, Nframes, Nmc
integer, dimension(:), intent(in) :: tau, mopt
real(sifm_real), dimension(:,:), intent(in) :: X
character(len=1), dimension(:,:), intent(inout) :: global_Xs
!
integer :: Nbin_min, Nbin_max, Ndata
integer :: i, ibin, nbin, ii, ierr, t, Istep
real(sifm_real) :: delta, x_max, x_min
character(len=1), dimension(nframes) :: XS_temp
real(sifm_real), allocatable, dimension(:) :: Xc
character(len=1), allocatable, dimension(:) :: Local_Xs
real(sifm_real) :: time_start, time_end
!
real(sifm_real) :: rnd, Global_H_Max, Local_H_Max, Local_H
integer :: Local_Nbin, Global_Nbin
!
allocate( Local_Xs(Nframes), stat=ierr )
if (ierr /= 0) STOP "Error allocating Local-Xs"
!
CALL CPU_time( Time_start )
!
DOF_LOOP: DO i = 1, Ndof
      write(*,'("Collect Info for IDIM:  ", I5)') I
      Ndata    = Nframes - (Mopt(i) - 1) * Tau(i)
      Nbin_min = floor( (1.0_sifm_real*Ndata)**(1.0_sifm_real/Mopt(i)) ) + 1
      IF (Nbin_min < 2) Nbin_min = 2
      Nbin_max = Nbin_min
      write(*,'("Nbin (Min,Max):  ", 2I10)') Nbin_min, Nbin_Max
      X_min = minval( X(i, :) )
      X_max = maxval( X(i, :) )

! -- Initialize a set of random numbers
      CALL rndnum_iniall(nbin_max-nbin_min + 1, 1, 0)
!
      ii = 0
      Local_H_Max = rzero
      Local_Nbin = 0
      BIN_LOOP: do nbin = nbin_min, nbin_max
        ii=ii+1
!
        allocate( Xc(nbin+1), stat=ierr )
        if (ierr /= 0) STOP "Error allocating Xc"
!
! -- Set XC(1) = min(X) and XC(nbin+1) = max(X)
        Xc(1) = X_min
        Xc(nbin+1) = X_max
        MC_LOOP: Do istep = 1, Nmc
!
         DO ibin=2, Nbin
            Xc(ibin) = ( Xc(ibin-1) - Xc(Nbin+1) ) * randf(ii) + Xc(Nbin+1)
         ENDDO
!
         TIME_LOOP: Do t=1, nframes
            XS_temp(t) = ""
            Do ibin = 1, Nbin
               IF ( X(i,t)>=(Xc(ibin)-prtiny)  .and.  X(i,t)<=(Xc(ibin+1)+prtiny) ) THEN
                    XS_temp(t) = Symbol(ibin)
               ENDIF
            ENDDO
         ENDDO TIME_LOOP
!
         CALL symbolic_entropy1D(Nframes, XS_temp, mopt(i), tau(i), Local_H)
!
         IF ( Local_H > Local_H_Max) THEN
              Local_H_Max = Local_H
              Local_XS = XS_temp
              Local_Nbin  = Nbin
         ENDIF
!
        ENDDO MC_LOOP
!
        deallocate( Xc, stat=ierr )
        IF (ierr /= 0) STOP "Error deallocating Xc in Symbolize Traj using MC method"
!
      enddo BIN_LOOP
!
      Global_Nbin = Local_Nbin
      Global_XS(i,:) = Local_XS(:)
      Write(*,'("Idim=   # of bins along the direction:  ", I5, 2X, I10)') i, global_Nbin
!
      CALL rndnum_clear()
!
ENDDO DOF_LOOP
!
deallocate( Local_Xs, stat=ierr )
if (ierr /= 0) STOP "Error deallocating Local-Xs"
!
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for Symbolize Traj - Monte Carlo:   ", F10.6)') Time_End - Time_Start
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
integer :: Nbin_min, Nbin_max, ndata
integer :: i, ia, ibin, nbin, ii, ierr, t, nbin_arg
real(sifm_real) :: delta, h, h_max, x_max, x_min
character(len=1), dimension(nframes) :: XS_temp
real(sifm_real), allocatable, dimension(:) :: xc
!
write(*,*) "Optimize nbin and Critical points:"
!
ATOMS_LOOP: DO ia = 1, Natoms
DIM_LOOP: do i=1, ndim
!
  X_min = minval( X(ia,i,:) )
  X_max = maxval( X(ia,i,:) )
!
  !Ndata = Nframes - (Mopt(ia,i) - 1) * Topt(ia,i)
  !Nbin_min = floor( (1.0_sifm_real*ndata)**(1.0_sifm_real/Mopt(ia,i)) ) + 1
  !Nbin_max = Ndata
!
  Nbin_min = 2
  Nbin_max = mopt(ia,i)
!
  ii=0
  H_max = 0.0_sifm_real
  BIN_LOOP: do nbin=nbin_min, nbin_max
    ii=ii+1
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
       XS_temp(t) = ""
       DO ibin = 1, Nbin
          IF ( X(ia,i,t) >= (Xc(ibin)-prtiny)  .and.  X(ia,i,t) <= (Xc(ibin+1)+prtiny) ) THEN
               XS_temp(t) = symbol(ibin)
          ENDIF
       ENDDO
    ENDDO
    CALL symbolic_entropy1D(nframes, xs_temp, mopt(ia,i), topt(ia,i), h)
    IF ( h > h_max ) THEN
         nbin_arg = Nbin
         h_max = h
         XS(ia,i,:) = XS_temp
    ENDIF
!
    deallocate( xc, stat=ierr )
    IF (ierr /= 0) STOP "Error deallocating xc"
!
  enddo BIN_LOOP
!
  Write(*,'("Idim=    # of bins along the direction and Max entropy:  ", 3I5, F12.6)') i, Nbin_arg, mopt(ia,i), H_Max
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
             XS(ia,i,t) = symbol(1)
         ELSE
             XS(ia,i,t) = symbol(2)
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
   XL = addStrings( XL, XS(L-(k-1)*Tau) )
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
