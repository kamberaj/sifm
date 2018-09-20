!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: embdmodule_ser.f90,v 1.0 19-03-2018, IBU
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
module EMBD_class
use sifm_kinds
use TEutils_class
use TErandom_class
implicit none
!
      integer, save :: Natoms, Nframes
      integer, save :: Ndim, Ndof
      integer, save :: Mopt_min, Mopt_max, Topt_min, Topt_max
      integer, save, allocatable :: mopt(:), tau(:)
      integer, save :: debug
      real(sifm_real), dimension(:,:), allocatable :: X

!
CONTAINS
!
subroutine EMBD_DRIVER(Nf, Nd, Na, M1, M2, Tau1, tau2, db)
  implicit none
  integer, intent(in) :: Nf, Nd, Na
  Integer, intent(in) :: db, M1, M2, tau1, tau2
!  
  real(sifm_real)      :: Time_Start, Time_end
!
  Nframes = Nf
  Ndim = Nd
  Natoms = Na
  Ndof = Ndim * Natoms
  Debug = db
  Mopt_min = M1
  Mopt_max = M2
  Topt_min = Tau1
  Topt_max = Tau2
!
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
  endif
!
  CALL Allocate_EmbeddedDim() 
!
  Call Read_xyz()
!
  CALL CPU_time(Time_start)
  CALL CompEmbeddedDimensions(tau, mopt)
!
  CALL CPU_time( Time_End )
  write(*,'("CPU Elapsed time for EMBD module:   ", F10.6)') Time_End - Time_Start
!
  write(*,*) "Write Embedded Parameters to a file"
  Call Print_embdparam()
  write(*,*) "... Done!"
!
  CALL DEALLOCATE_EmbeddedDim()
!
  return
end subroutine EMBD_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the trajectories
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xyz()
Implicit None
!  
  integer :: i, d, j, t, nunit, time, atindex
  character(len=7) :: remark
  character(len=80) :: line
  real(sifm_real), dimension(Ndim) :: xyz
  
!
nunit = new_unit()
open(unit=nunit, file='traj.xyz', status='old', action='read')
read(nunit, '(A80)') line
write(*, *) line
read(nunit, '(A80)') line
write(*, *) line
DO t = 1, Nframes
   j = 0
   DO i = 1, Natoms
      read(nunit,'(A80)') line
      IF (line(1:5) == 'FRAME') THEN
          read(line,'(A7, 2I10, 3F12.6)') remark, time, atindex, ( XYZ(d), d = 1, Ndim ) 
          write(*,'(A7, 2I10, 3F12.6)') remark, time, atindex, ( XYZ(d), d = 1, Ndim ) 
          DO d = 1, Ndim
             j = j + 1
             X(j,t) = xyz(d)
          ENDDO
      ENDIF
   ENDDO
ENDDO
read(nunit, '(A80)') line
CLOSE(nunit)
!
return
!
end subroutine read_xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the embedded dimension parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_embdparam()
Implicit None
!  
  integer :: i, j, d, nunit
  integer, dimension (Ndim) :: T, M
!
  nunit = new_unit()
  open(unit=nunit, file='embd.txt', status='unknown', action='write')
!  
  j=0
  DO i = 1, Natoms
     DO d = 1, Ndim
        J = J + 1
        M(d) = Mopt(J)
        T(d) = Tau(J)
     ENDDO
     write(nunit,*) (t(d), d=1, Ndim), (m(d), d=1, Ndim)
  ENDDO
  CLOSE(nunit)
!
return
!
end subroutine Print_embdparam
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Allocate_embeddedDim()
      implicit none
!
      integer :: ierr
!
!
      IF (.not. Allocated(X)) THEN
          Allocate(X(Ndof,Nframes), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating X"
      ENDIF
      IF (.not. Allocated(tau)) THEN
          Allocate(Tau(Ndof), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating Tau"
      ENDIF
!
      IF (.not. Allocated(Mopt)) THEN
	  Allocate(mopt(Ndof), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating Mopt"
      ENDIF   
!
      return
end subroutine Allocate_EmbeddedDim
!
subroutine deAllocate_EmbeddedDim()
      implicit none
!
      integer :: ierr
!  
      IF (Allocated(X)) THEN
	  DEAllocate(X, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating X"
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
      return
end subroutine deAllocate_EmbeddedDim
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Compute the time shift and state vector dimension of time series
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CompEmbeddedDimensions(gtau, gmopt)
Implicit none
!
integer, dimension(:), intent(inout) :: gtau,  gmopt
!
Integer :: nunit1, nunit2, i, j, m, im, Nsize, Ntsize
!
integer, allocatable :: gcount_fnn(:,:)
real(sifm_real), allocatable :: gMI(:,:)
integer :: ierr
real(sifm_real) :: time_start, time_end
!
!
nunit1 = new_unit()
open(unit=nunit1, file='tau.txt', status='unknown', action='write')
nunit2 = new_unit()
open(unit=nunit2, file='mopt.txt', status='unknown', action='write')

!
Nsize  = mopt_max - mopt_min + 1
Ntsize = topt_max - topt_min + 1
!
allocate(gcount_fnn(NDof,Nsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-count-fnn"
allocate(gMI(NDof,Ntsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-MI"
!
call cpu_time( time_start )
Do I = 1, Ndof
   CALL getTimeLag(Nframes, X(i,:), topt_min, topt_max, gtau(I), gMI(I,:))
ENDDO
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for TOPT:   ", F10.6)') Time_End - Time_Start
!
call cpu_time( time_start )
Do I = 1, Ndof
   CALL getOptEmbDim(Nframes, X(i,:), mopt_min, mopt_max, gtau(I), gMopt(I), gcount_fnn(I,:))
ENDDO
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for MOPT:   ", F10.6)') Time_End - Time_Start
!
DO i = 1, Ndof
   im=0
   DO m = mopt_min, mopt_max
      im = im + 1
      write(nunit2,*) I, m, gcount_fnn(I,im), Gmopt(i)
   ENDDO
ENDDO
!
DO i = 1, Ndof
   im = 0
   DO m = topt_min, topt_max
      im = im + 1
      write(nunit1,*) I, m, gMi(I,im), gTau(i)
   ENDDO
ENDDO
!
close(nunit1)
close(nunit2)
!
deallocate(Gcount_fnn, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-count-fnn"
deallocate(Gmi, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-mi"
!
CONTAINS
!
subroutine getTimeLag(N, x, tau_min, tau_max, tau_arg, MI)
implicit none
integer, intent(in) :: N
real(sifm_real), intent(in) :: x(:)
real(sifm_real), intent(inout) :: MI(:)
integer, intent(in) :: tau_max, tau_min
integer, intent(inout) :: tau_arg

! -- Local variables
integer :: i, ibin,i1, i2, norm, it, Tau, Nb, ierr
integer, allocatable, dimension(:) :: n1,n2
integer, allocatable, dimension(:,:) :: n12
real(sifm_real) :: xmin,xmax,delta_x,s1,s2,s12
!
    Xmin = minval(X)
    Xmax = maxval(X)
    call getOptimalNb_entropy(xmax-xmin, xmin, X, N, Nb)
    delta_x = (xmax-xmin)/real(Nb-1,sifm_real)
!
    allocate(n1(nb), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating n1"
    allocate(n2(nb), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating n2"
    allocate(n12(nb,nb), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating n12"
!
    It=0
    TAU_LOOP: DO Tau = tau_min, Tau_max
        It = It + 1
        N1=0; N2=0; N12=0
        DO I = 1, N-Tau
           I1 = floor( (x(i) - xmin) / delta_x ) + 1
           IF (I1 <= NB) THEN
               N1( I1 ) = N1( I1 ) + 1
           ENDIF
           I2 = floor( (x(i+tau) - xmin) / delta_x ) + 1
           IF (I2 <= NB) THEN
               N2( I2 ) = N2( I2 ) + 1
           ENDIF
           IF (I1 <= NB .and. I2 <= NB) THEN
               N12( I1, I2 ) = N12( I1, I2 ) + 1
           ENDIF
        ENDDO
        S1 = 0.0_sifm_real
        norm=sum( N1 )
        do i1 = 1, NB
           if ( n1(i1) > 0 ) then
                S1 = S1 + PSI( real(n1(i1),sifm_real) )*real(n1(i1),sifm_real)
           endif
        enddo
        S1 = log( real(norm,sifm_real) ) - S1 / real(norm,sifm_real)

        S2 = 0.0_sifm_real
        norm=sum( N2 )
        do i2 = 1, NB
           if ( n2(i2) > 0 ) then
                S2 = S2 + PSI( real(n2(i2),sifm_real) )*real(n2(i2),sifm_real)
           endif
        enddo
        S2 = log( real(norm,sifm_real) ) - S2 / real(norm,sifm_real)
   
        S12 = 0.0_sifm_real
        norm = 0
        do i1 = 1, NB
           do i2 = 1, NB
              if ( n12(i1,i2) > 0 ) then
                   S12 = S12 + PSI( real(n12(i1,i2),sifm_real) )*real(n12(i1,i2),sifm_real)
                   norm = norm + n12(i1,i2)
              endif
           enddo
        enddo
        S12 = log( real(norm,sifm_real) ) - S12 / real(norm,sifm_real)
        MI(It) = S1 + S2 - S12
    ENDDO TAU_LOOP

    ! -- Get the first minima
    Tau_arg=0
    DO i = 3, It
       IF ( MI(i) > MI(i-1) ) THEN
            Tau_arg = i-1
            GOTO 1
       ENDIF
    ENDDO
1   continue
    IF (Tau_arg == 0) Tau_arg = IT
!
    deallocate(n1, STAT=ierr )
    if (ierr /= 0) STOP "Error deallocating n1"
    deallocate(n2, STAT=ierr )
    if (ierr /= 0) STOP "Error deallocating n2"
    deallocate(n12, STAT=ierr )
    if (ierr /= 0) STOP "Error deallocating n12"
!
END SUBROUTINE getTimeLag
!
subroutine getOptimalNb_entropy(R, xmin, x, n, m_arg)
real(sifm_real), intent(in) :: R, xmin
real(sifm_real), dimension(:), intent(in) :: x
integer, intent(in) :: n
integer, intent(out) :: m_arg

! --- Local variables
integer :: m, i, i1, m1, m2, norm
real(sifm_real) :: f_prev, delta, f
integer, allocatable :: n1(:)
real(sifm_real), parameter :: tol=0.0001_sifm_real

m1 = 10
m2 = 5000 
m=m1
BIN_LOOP: DO WHILE (m <= m2)
  delta = R / real(m-1, sifm_real)
  allocate( n1(m) )
  n1 = 0
  DO I = 1, N
     I1 = floor( (x(i) - xmin) / delta ) + 1
     IF (I1 <= m) THEN
         N1( I1 ) = N1( I1 ) + 1
     ENDIF
  ENDDO
  ! ---- Compute entropy
  f = 0.0_sifm_real
  norm = sum( n1 )
  do i1 = 1, m
     if (n1(i1) > 0) then
         f = f + PSI( real(n1(i1),sifm_real) )*real(n1(i1),sifm_real)
     endif
  enddo
  f = log( real(norm,sifm_real) ) - f/real(norm,sifm_real)
  IF (m > m1) THEN 
  IF (abs(f - f_prev) <= tol) THEN
      m_arg=m
      exit
  ENDIF
  ENDIF
  f_prev = f
  m=m+1
  deallocate( n1 )
ENDDO BIN_LOOP
return
end subroutine getOptimalNb_entropy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Use the global false nearest neighbors to determine
!! the optimal embedding dimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getOptEmbDim(n, X, m1, m2, tau, Mopt, count_fnn)
!
integer, intent(in) :: n, m1, m2
real(sifm_real), dimension(:), intent(in) :: X
integer, intent(in) :: tau
integer, intent(inout) :: Mopt
integer, intent(inout) :: count_fnn(:)
!
integer :: im, m, i, k, j, ierr, ipos, T, Neff
real(sifm_real) :: R, dist1, dist2, min_dist1
real(sifm_real), parameter :: Rtol=10.0_sifm_real
!
count_fnn = 0
T = Tau
im=0 
!
M_LOOP: DO m = m1, m2
   im = im + 1
   Neff = N - m*T
   DO I = 1, Neff-2
      min_DIST1 = rzero
      Do k = 0, m-1
         min_DIST1 = min_DIST1 + ( X(I+k*T) - X(I+1+k*T) )**(2.0_sifm_real)
      ENDDO
      ipos=I+1
      DO J = I+2, Neff
         DIST1 = rzero
         DO k = 0, m-1
            DIST1 = DIST1 + ( X(I+k*T) - X(J+k*T) )**(2.0_sifm_real)
         ENDDO
         IF (min_dist1 > dist1) THEN
             ipos = J
             min_dist1 = dist1
         ENDIF
      ENDDO  
      dist2 = (  X(I+M*T) - X(ipos+M*T) )**(2.0_sifm_real)
      IF (min_dist1 <= rzero) THEN
          min_dist1=PRTINY
      ENDIF
      R = sqrt( dist2 / min_dist1 )
      IF (R >= Rtol) count_fnn(im) = count_fnn(im) + 1
   ENDDO
ENDDO M_LOOP
!
Mopt=m2
im=0
DO m = m1, m2
     im = im+1
     IF (count_fnn(im) <= 0) THEN
        Mopt = m
        GOTO 1
     ENDIF
ENDDO
1 CONTINUE
!
return
end subroutine getOptEmbDim
!
end subroutine CompEmbeddedDimensions
!
end module EMBD_CLASS
