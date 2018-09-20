!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: embdmodule.f90,v 1.0 19-03-2018, IBU
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
  use TEMPI_CLASS
  implicit none
  integer, intent(in) :: Nf, Nd, Na
  Integer, intent(in) :: db, M1, M2, tau1, tau2
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
  call MPI_start()
  call MPI_Broadcast()
!
  Localtime_start = MPI_Wtime()
  if (Debug == 1) then
      CALL CompEmbeddedDimensions_slow(tau, mopt)
  else
      CALL CompEmbeddedDimensions_fast()
  endif
! -- Timing Analysis  
  LocalTime_end = MPI_WTIME()
  call MPI_REDUCE ( LocalTime_end,   GlobalTime_end,   1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  call MPI_REDUCE ( LocalTime_start, GlobalTime_start, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr )      
  IF (myID == Master) THEN
      write(*,'("CPU Elapsed time for EMBD parameters:   ", F20.3)') &
                 (GlobalTime_End - GlobalTime_Start)/real(NumProcs, sifm_real)
  ENDIF
!
  write(*,*) "Write Embedded Parameters to a file"
  Call Print_embdparam()
  write(*,*) "... Done!"
!
  CALL DEALLOCATE_EmbeddedDim()
!
  call MPI_Finish()
!
  return
end subroutine EMBD_DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Broadcast to every node
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MPI_Broadcast()
use TEMPI_CLASS
implicit none
!
integer :: ierr
!
call MPI_BCAST (Mopt,    Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Mopt_min,1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Mopt_max,1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

call MPI_BCAST (Topt_min,1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Topt_max,1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Tau,     Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

call MPI_BCAST (Natoms,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Ndof,    1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST (Nframes, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

call MPI_BCAST (X,       Ndof*NFrames, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr )
!
return
end subroutine MPI_Broadcast

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
subroutine CompEmbeddedDimensions_fast()
use TEMPI_CLASS
Implicit none
!
Integer :: nunit1, nunit2, i, m, im, Nsize, Ntsize, K
integer, allocatable :: Ltau(:), Lmopt(:)
real(sifm_real) :: time_end, time_start
integer :: ierr, Stride, Istart, iEnd
!
IF ( mod(Ndof, NumProcs) == 0 ) THEN
     Stride = Ndof / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd
ELSE
     Stop "Number of DoF must be Multiple of Number of Processors"
ENDIF
allocate(Ltau(Stride), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-tau"
allocate(Lmopt(Stride), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-M"
!
IF (myID == master) THEN
    nunit1 = new_unit()
    open(unit=nunit1, file='tau.txt', status='unknown', action='write')
    nunit2 = new_unit()
    open(unit=nunit2, file='mopt.txt', status='unknown', action='write')
ENDIF
!
Nsize  = mopt_max - mopt_min + 1
Ntsize = topt_max - topt_min + 1
!
if (myID == master) &
    call cpu_time(time_start)
Do I = iStart, iEnd
   K = I - iStart + 1
   CALL getTimeLag_fast(I, Ltau(K))
ENDDO
if (myID == master) then
    CALL CPU_time( Time_End )
    write(*,'("CPU Elapsed time for TOPT:   ", I5, F10.6)') myID, Time_End - Time_Start
endif
!
CALL MPI_gather( Ltau, Stride, MPI_INTEGER, &
                 tau,  Stride, MPI_INTEGER, &
                 master, MPI_COMM_WORLD, ierr )
!
call MPI_BCAST (tau, Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
if (myID == master) &
     call cpu_time(time_start)
Do I = iStart, iEnd
   K = I - iStart + 1
   call getOptEmbDim_fast(I, LMopt(K))
ENDDO
if (myID == master) then
    CALL CPU_time( Time_End )
    write(*,'("CPU Elapsed time for MOPT:   ", F10.6)') Time_End - Time_Start
endif
!
CALL MPI_gather( Lmopt, Stride, MPI_INTEGER, &
                 mopt, Stride, MPI_INTEGER, &
                 master, MPI_COMM_WORLD, ierr )
!
call MPI_BCAST (Mopt, Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
IF ( myID == master ) THEN
     im = 0
     DO i = 1, Ndof
        DO m = mopt_min, mopt_max
           im = im + 1
           write(nunit2,*) I, m, mopt(i)
        ENDDO
    ENDDO
ENDIF
!
IF ( myID == master ) THEN
     im = 0
     DO i = 1, Ndof
        DO m = topt_min, topt_max
           im = im + 1
           write(nunit1,*) I, m, Tau(i)
        ENDDO
    ENDDO
ENDIF
!
IF (myID == master) THEN
    close(nunit1)
    close(nunit2)
ENDIF
!
deallocate(Ltau, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-tau"
deallocate(Lmopt, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-M"
!
return
!
END subroutine CompEmbeddedDimensions_fast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the time Lag using the fast routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getTimeLag_fast(Idof, tau_arg)
use sifm_kinds
implicit none
integer, intent(in) :: Idof
integer, intent(out) :: tau_arg

! -- Local variables
integer :: i, ibin,i1, i2, norm, it, T, Nb, ierr, Nsize
integer, allocatable, dimension(:) :: n1,n2
integer, allocatable, dimension(:,:) :: n12
real(sifm_real) :: xmin,xmax,delta_x,s1,s2,s12
real(sifm_real), allocatable :: MI(:)
!
    Nsize = Topt_Max - Topt_Min + 1
    Xmin = minval(X(idof,:))
    Xmax = maxval(X(idof,:))
    call getOptimalNb_entropy_fast(xmax - xmin, xmin, idof, Nb)
    delta_x = (xmax-xmin)/real(Nb-1,sifm_real)
!
    allocate(n1(nb), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating n1"
    allocate(n2(nb), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating n2"
    allocate(n12(nb,nb), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating n12"
    allocate(MI(nsize), STAT=ierr )
    if (ierr /= 0) STOP "Error allocating MI"
!
    It=0
    TAU_LOOP: DO T = Topt_Min, Topt_Max
        It = It + 1
        N1=0; N2=0; N12=0
        DO I = 1, Nframes-T
           I1 = floor( (x(idof,i) - xmin) / delta_x ) + 1
           IF (I1 <= NB) THEN
               N1( I1 ) = N1( I1 ) + 1
           ENDIF
           I2 = floor( (x(idof,i+t) - xmin) / delta_x ) + 1
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
END SUBROUTINE getTimeLag_fast
!
subroutine getOptimalNb_entropy_fast(R, xmin, idof, m_arg)
use sifm_kinds
implicit none
real(sifm_real), intent(in) :: R, xmin
integer, intent(in) :: idof
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
  DO I = 1, Nframes
     I1 = floor( (x(idof,i) - xmin) / delta ) + 1
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
end subroutine getOptimalNb_entropy_fast
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Use the global false nearest neighbors to determine
!! the optimal embedding dimension - Fast routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getOptEmbDim_fast(idof, Mopt_arg)
use sifm_kinds
implicit none
!
integer, intent(in) :: idof
integer, intent(out) :: Mopt_arg
!
integer :: im, m, i, k, j, ierr, ipos, T, Neff, nsize
real(sifm_real) :: R, dist1, dist2, min_dist1
real(sifm_real), parameter :: Rtol=10.0_sifm_real
integer, allocatable :: count_fnn(:)
!
nsize = Mopt_max - Mopt_min + 1
allocate( count_fnn(nsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating count-fnn"
count_fnn = 0
!
T = Tau(idof)
im=0 
!
M_LOOP: DO m = Mopt_Min, Mopt_Max
   im = im + 1
   Neff = Nframes - m * T
   DO I = 1, Neff-2
      min_DIST1 = rzero
      Do k = 0, m-1
         min_DIST1 = min_DIST1 + ( X(idof,I+k*T) - X(idof,I+1+k*T) )**(2.0_sifm_real)
      ENDDO
      ipos=I+1
      DO J = I+2, Neff
         DIST1 = rzero
         DO k = 0, m-1
            DIST1 = DIST1 + ( X(idof,I+k*T) - X(idof,J+k*T) )**(2.0_sifm_real)
         ENDDO
         IF (min_dist1 > dist1) THEN
             ipos = J
             min_dist1 = dist1
         ENDIF
      ENDDO  
      dist2 = (  X(idof,I+M*T) - X(idof,ipos+M*T) )**(2.0_sifm_real)
      IF (min_dist1 <= rzero) THEN
          min_dist1=PRTINY
      ENDIF
      R = sqrt( dist2 / min_dist1 )
      IF (R >= Rtol) count_fnn(im) = count_fnn(im) + 1
   ENDDO
ENDDO M_LOOP
!
Mopt_arg=Mopt_Max
im=0
DO m = Mopt_Min, Mopt_Max
     im = im+1
     IF (count_fnn(im) <= 0) THEN
        Mopt_arg = m
        GOTO 1
     ENDIF
ENDDO
1 CONTINUE
!
deallocate( count_fnn, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating count-fnn"
!
return
end subroutine getOptEmbDim_fast
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Compute the time shift and state vector dimension of time series
!!  Slow routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CompEmbeddedDimensions_slow(gtau, gmopt)
use TEMPI_CLASS
Implicit none
!
integer, dimension(:), intent(inout) :: gtau,  gmopt
!
Integer :: nunit1, nunit2, i, m, im, Nsize, Ntsize, K
integer, allocatable :: Ltau(:), Lmopt(:)
!
integer, allocatable :: Lcount_fnn(:)
integer, allocatable :: gcount_fnn(:)
real(sifm_real), allocatable :: LMI(:)
real(sifm_real), allocatable :: gMI(:)
real(sifm_real) :: time_end, time_start
!
integer :: ierr, Stride, Istart, iEnd
!
IF ( mod(Ndof, NumProcs) == 0 ) THEN
     Stride = Ndof / NumProcs
     iStart = (myID * Stride) + 1
     iEnd   = iStart + Stride - 1
     write(*,*) myID, Stride, iStart, iEnd
ELSE
     Stop "Number of DoF must be Multiple of Number of Processors"
ENDIF
allocate(Ltau(Stride), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-tau"
allocate(Lmopt(Stride), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-M"
!
IF (myID == master) THEN
    nunit1 = new_unit()
    open(unit=nunit1, file='tau.txt', status='unknown', action='write')
    nunit2 = new_unit()
    open(unit=nunit2, file='mopt.txt', status='unknown', action='write')
ENDIF
!
Nsize  = mopt_max - mopt_min + 1
Ntsize = topt_max - topt_min + 1
!
allocate(Lcount_fnn(Stride*Nsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-count-fnn"
allocate(gcount_fnn(NDof*Nsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-count-fnn"
allocate(Lmi(Stride*Ntsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating local-MI"
allocate(gMI(NDof*Ntsize), STAT=ierr )
if (ierr /= 0) STOP "Error allocating global-MI"
!
if (myID == master) &
    call cpu_time(time_start)
Do I = iStart, iEnd
   K = I - iStart + 1
   CALL getTimeLag_slow(Nframes, X(I,:), topt_min, topt_max, Ltau(K), LMI( ((K-1)*Ntsize+1):(Ntsize*K) ) )
ENDDO
if (myID == master) then
    CALL CPU_time( Time_End )
    write(*,'("CPU Elapsed time for TOPT:   ", I5, F10.6)') myID, Time_End - Time_Start
endif
!
CALL MPI_gather( Ltau, Stride, MPI_INTEGER, &
                 Gtau, Stride, MPI_INTEGER, &
                 master, MPI_COMM_WORLD, ierr )
!
call MPI_BCAST (gTau, Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
if (myID == master) call cpu_time(time_start)
Do I = iStart, iEnd
   K = I - iStart + 1
   CALL getOptEmbDim_slow(Nframes, X(i,:), mopt_min, mopt_max, gtau(I), LMopt(K), &
                     Lcount_fnn( ((K-1)*Nsize+1):(Nsize*K) ) )
ENDDO
if (myID == master) then
    CALL CPU_time( Time_End )
    write(*,'("CPU Elapsed time for MOPT:   ", F10.6)') Time_End - Time_Start
endif
!
CALL MPI_gather( Lmopt, Stride, MPI_INTEGER, &
                 Gmopt, Stride, MPI_INTEGER, &
                 master, MPI_COMM_WORLD, ierr )
!
call MPI_BCAST (gMopt, Ndof, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
!
CALL MPI_gather( Lcount_fnn, Nsize*Stride, MPI_INTEGER, &
                 Gcount_fnn, Nsize*Stride, MPI_INTEGER, &
                 master, MPI_COMM_WORLD, ierr )

CALL MPI_gather( Lmi, Ntsize*Stride, MPI_DOUBLE_PRECISION, &
                 Gmi, Ntsize*Stride, MPI_DOUBLE_PRECISION, &
                 master, MPI_COMM_WORLD, ierr )
!
IF ( myID == master ) THEN
     im = 0
     DO i = 1, Ndof
        DO m = mopt_min, mopt_max
           im = im + 1
           write(nunit2,*) I, m, gcount_fnn(im), Gmopt(i)
        ENDDO
    ENDDO
ENDIF
!
IF ( myID == master ) THEN
     im = 0
     DO i = 1, Ndof
        DO m = topt_min, topt_max
           im = im + 1
           write(nunit1,*) I, m, gMi(im), gTau(i)
        ENDDO
    ENDDO
ENDIF
!
IF (myID == master) THEN
    close(nunit1)
    close(nunit2)
ENDIF
!
deallocate(Ltau, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-tau"
deallocate(Lmopt, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-M"
deallocate(Lcount_fnn, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-count-fnn"
deallocate(Gcount_fnn, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-count-fnn"
deallocate(Lmi, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating local-mi"
deallocate(Gmi, STAT=ierr )
if (ierr /= 0) STOP "Error deallocating global-mi"
!
return
!
END subroutine CompEmbeddedDimensions_slow
!
subroutine getTimeLag_slow(N, x, tau_min, tau_max, tau_arg, MI)
use sifm_kinds
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
    call getOptimalNb_entropy_slow(xmax-xmin, xmin, X, N, Nb)
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
END SUBROUTINE getTimeLag_slow
!
subroutine getOptimalNb_entropy_slow(R, xmin, x, n, m_arg)
use sifm_kinds
implicit none
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
end subroutine getOptimalNb_entropy_slow
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Use the global false nearest neighbors to determine
!! the optimal embedding dimension - Slow routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getOptEmbDim_slow(n, X, m1, m2, tau, Mopt, count_fnn)
use sifm_kinds
implicit none
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
   Neff = N - m * T
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
end subroutine getOptEmbDim_slow
!
end module EMBD_CLASS
