!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: syssetup.f90,v 1.0 19-03-2018, IBU
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
module System_class
use sifm_kinds
use TEutils_class
use TErandom_class
implicit none
!
      integer, save :: Natoms, Nframes
      integer, save :: Ndim, Offset
      real(sifm_real), save, allocatable :: xyz(:,:,:)
      integer, save, allocatable :: xyzs(:,:,:)
      integer, save :: debug

!
CONTAINS
!
subroutine Setup_DRIVER(Nf, Nd, Na, Off,                            &
                       iStart, iSkip, iStop,                        &
                       db, benchmark,                               &
                       Ax, Cxy, Sigma_x, Sigma_y)
  implicit none
  integer, intent(in) :: Nf, Nd, Na, iSkip, iStop, iStart, benchmark, db, off 
  real(sifm_real), intent(in), optional :: Ax, Cxy, Sigma_x, Sigma_y
!  
  real(sifm_real)      :: Time_Start, Time_end
!
  Nframes = Nf
  Ndim = Nd
  Natoms = Na
  Debug = db
  Offset = off
!
  if (debug == 1) then
      write(*,'("Nframes=   ", I5)') Nframes
      write(*,'("Ndim=      ", I5)') Ndim
      write(*,'("Natoms=    ", I5)') Natoms
      write(*,'("Benchmark= ", I5)') benchmark
  endif
!
  write(*,*) "Memory ..."
  CALL Allocate_system(benchmark) 
  write(*,*) "... allocated"
!
  IF (benchmark == 1) THEN
      call Benchmark1(Ax, Cxy, Sigma_x, Sigma_y, iStart, iSkip, iStop)
  ELSEIF (benchmark == 2) THEN
      call Benchmark2(iStart, iSkip, iStop)
  ELSEIF (benchmark == 3) THEN
      call Benchmark3(iStart, iSkip, iStop)
  ELSEIF (benchmark == 4) THEN
      call Benchmark4(iStart, iSkip, iStop)
  ELSEIF (benchmark == 5) THEN
      call Benchmark5(iStart, iSkip, iStop)
  ELSEIF (benchmark == 6) THEN  ! Lorenz system
      call read_xyz()
  ELSEIF (benchmark == 7) THEN  
      call mdTrajPDBVmd(Istart, iSkip, iStop)
  ELSE
      call mdTrajEnergy(Istart, iSkip, iStop)
  ENDIF
!
  call write_xyz()
  IF ( (benchmark == 4) .or. (benchmark == 5) ) THEN
      call write_xyzs()      
  ENDIF
!
  CALL DEALLOCATE_System()
!
  return
end subroutine setup_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mdTrajEnergy(Istart, iSkip, iStop)
implicit none
integer, intent(in) :: iStart, iSkip, iStop
!
integer           :: t, tt, nunit
integer           :: i, j, k, ii
real(sifm_real), dimension(Natoms) :: xyz_t
!   
nunit = new_unit()
open(unit=nunit, file='fox-rna_mmenergy.txt', status='old', action='read')
!
tt = 0
DO t = 1, iStop
   Read(nunit,*) i, xyz_t
   Write(*,*) i, xyz_t
   IF ( (t >= iStart) .and. (mod(t, iSkip) == 0) .and. (tt < Nframes) ) THEN
         tt = tt + 1
         do j = 1, Natoms
            xyz(j,1,tt) = xyz_t(j)
            xyz(j,2,tt) = 0.0_sifm_real
            do k = 1, Natoms
               IF (k /= j) THEN
                   xyz(j,2,tt) = xyz(j,2,tt) + xyz_t(k)
               ENDIF
            enddo
         enddo
   ENDIF
ENDDO 
DO i = 1, Natoms
   DO j = 1, Ndim
      Xyz(i,j,:) = Xyz(i,j,:) - mean( Xyz(i,j,:), Nframes )
   ENDDO
ENDDO
return
end subroutine MDTrajEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mdTrajPDBVmd(Istart, iSkip, iStop)
      implicit none
      integer, intent(in) :: iStart, iSkip, iStop
!
      integer           :: ia, j, t, tt, iend, ires_prev, nunit
      integer           :: t_ires, ires, t_atindex
      character(len=1)  :: t_chain
      character(len=3)  :: t_resname
      character(len=4)  :: t_segname
      character(len=5)  :: t_atname
      character(len=6)  :: atom
      character(len=80) :: Remark
      real(sifm_real)    :: t_x, t_y, t_z, t_a, t_b
      real(sifm_real), dimension(Natoms, Ndim) :: xyz_t
!   
      nunit = new_unit()
      open(unit=nunit, file='traj.pdb', status='old', action='read')
!
xyz=999.0_sifm_real
read(nunit,*) Remark
tt = 0
DO t = 1, iStop
   ia = 0
   ires_prev = Offset
   ires = 1
   iend = 0
   xyz_t = rzero
   do while (iend == 0) 
         Read(nunit,100,end=90) atom, t_atindex, t_atname, t_resname, t_chain, t_ires, &
              t_x, t_y, t_z, t_a, t_b, t_segname
         write(*,100) atom, t_atindex, t_atname, t_resname, t_chain, t_ires, &
              t_x, t_y, t_z, t_a, t_b, t_segname
         t_segname = adjustL(t_segname)
         If (atom(1:4) == 'ATOM') then
             if (t_ires /= ires_prev) then
                xyz_t(ires,1) = xyz_t(ires,1) / real(Ia, sifm_real)
                xyz_t(ires,2) = xyz_t(ires,2) / real(Ia, sifm_real)
                xyz_t(ires,3) = xyz_t(ires,3) / real(Ia, sifm_real)
                ires_prev = t_ires
                ires = ires + 1
                Ia = 0
            endif
            ia = ia + 1
            t_atname = adjustL(t_atname)
            IF (ires > Natoms) THEN
                write(*,*) Ires, Natoms
                stop 'Number of atoms exceeds!!!'
            ENDIF
            xyz_t(ires,1) = xyz_t(ires,1) + t_x
            xyz_t(ires,2) = xyz_t(ires,2) + t_y
            xyz_t(ires,3) = xyz_t(ires,3) + t_z
         Elseif (atom(1:3) == 'END') then
            xyz_t(ires,1) = xyz_t(ires,1) / real(Ia, sifm_real)
            xyz_t(ires,2) = xyz_t(ires,2) / real(Ia, sifm_real)
            xyz_t(ires,3) = xyz_t(ires,3) / real(Ia, sifm_real)
            goto 92
         endif
         GOTO 91
90       Continue
         iend = 1
91       continue
   ENDDO
92 continue 
IF ( (t >= iStart) .and. (mod(t, iSkip) == 0) .and. (tt < Nframes) ) THEN
      tt = tt + 1
      xyz(:, :, tt) = xyz_t
ENDIF
ENDDO 
DO ia = 1, Natoms
   DO j = 1, Ndim
      Xyz(ia,j,:) = Xyz(ia,j,:) - mean( Xyz(ia,j,:), Nframes )
   ENDDO
ENDDO
100   format(a6, i5, a5, 1x, a3, 1x, A, i4, 4x, 3F8.3, 2F6.2, 4x, a4)    
      return
end subroutine MDTrajPdbVMD
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  The Multivariate Gaussian Process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Benchmark1(Ax,Cxy,Sigmax,Sigmay,iStart,iSkip,iStop)
implicit none
!
real(sifm_real), intent(in) :: Ax, Cxy, Sigmax, Sigmay
integer, intent(in) :: iStart, iStop, iSkip
!
! -- Local variables
integer :: i,d,t,ierr,tt,Nt
real(sifm_real), dimension(:,:,:), allocatable :: xyz_t
real(sifm_real) :: Time_start, Time_end
integer, parameter :: t0 = 1
!
! --- Some initial parameters
!
CALL CPU_time(Time_start)
!
Nt = iStop
ALLOCATE( xyz_t(Natoms, Ndim, Nt), stat=ierr )
if (ierr /= 0) stop 'Error allocating xyz-t'
!
! -- Initialize Some random number sequences
CALL rndnum_iniall(2*Ndim, 1, 0)
!
t = 1
DO d=1, Ndim
   xyz_t(1,d,t) = sigmaX * gauss_rand(d)
   xyz_t(2,d,t) = sigmaY * gauss_rand(Ndim + d)
enddo
Do d = 1, Ndim
   tt = 0
   DO t = 2, Nt
      xyz_t(1,d,t)     = Ax  * xyz_t(1,d,t-t0) + sigmaX * gauss_rand(d)
      IF (t > t0) THEN
          xyz_t(2,d,t) = Cxy * xyz_t(1,d,t-t0) + sigmaY * gauss_rand(Ndim+d)
      ELSE
          xyz_t(2,d,t) = sigmaY * gauss_rand(Ndim+d)
      ENDIF
      IF ( t >= iStart .and. mod(t - 1, iSkip) == 0 ) THEN
           tt = tt + 1
           IF (tt <= Nframes) THEN
               xyz(1,d,tt) = xyz_t(1,d,t)
               xyz(2,d,tt) = xyz_t(2,d,t)
           ENDIF
      ENDIF
   ENDDO
ENDDO
!
DO i = 1, Natoms
   DO d = 1, Ndim
      Xyz(i,d,:) = Xyz(i,d,:) - mean( Xyz(i,d,:), Nframes )
   ENDDO
ENDDO
!
CALL rndnum_clear()
!
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for EMBD module:   ", F10.6)') Time_End - Time_Start
!
DEALLOCATE( xyz_t, stat=ierr )
if (ierr /= 0) stop 'Error allocating xyz-t'
!
RETURN
!
end subroutine Benchmark1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Mutual information in this Experiment:
!!  Y = X+Z
!!  X uniformly distributed in [-1/2:1/2]
!!  Z uniformly distributed in [-a/2:a/2]
!!  f_y(y) = convolution(X,Z) 
!!  I(X;Y) = h(Y) - h(Z)
!!  In Nats: I(X;Y) = a/2  - ln(a) if a <= 1; otherwise I(X;Y)=1/(2*a) if a>=1.
!!  In Bits: I(X;Y) = [a/2 - ln(a)]/ln(2) if a <= 1; otherwise I(X;Y)=[1/(2*a)]/ln(2) if a>=1.
!!  Example: a=0.5, then I(X;Y) = 1.36067 (bits); h(Y)=0.49075 (bits)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Benchmark2(Istart, Iskip, Istop)
implicit none
Integer, intent(in) :: istart, iskip, istop
!
! -- Local variables
integer :: i,d,t,tt
real(sifm_real):: u, x_t, y_t, z_t, a2
real(sifm_real), parameter :: A = 0.5_sifm_real
real(sifm_real) :: Time_start, Time_end
!
CALL CPU_time( Time_Start )
! --- Some initial parameters
a2 = A * 0.5_sifm_real

!
! -- Initialize Some random number sequences
CALL rndnum_iniall(Ndim, 1, 0)
!
Do d = 1, Ndim
   tt = 0
   t  = 0
   DO while ( tt < Nframes )
      t = t + 1
      IF ( (t >= iStart) .and. (mod(t,iSkip) == 0) ) THEN
           IF ( tt < Nframes ) Then
                u = randf(d)
                x_t = u - 0.5_sifm_real
                IF ( u <= A2 ) THEN
                    y_t = -0.5_sifm_real * (one + A) + sqrt(2.0_sifm_real * A * u)
                ELSEIF ( u <= (one - A2) ) THEN
                    y_t = u - 0.5_sifm_real
                ELSE
                    y_t =  0.5_sifm_real * (one + A) - sqrt( 2.0_sifm_real * A * (one - u) )
                ENDIF
                tt = tt + 1
                xyz(1,d,tt) = x_t
                xyz(2,d,tt) = y_t
           ENDIF
      ENDIF
   ENDDO
ENDDO
!
CALL rndnum_clear()
!
DO i = 1, Natoms
   DO d = 1, Ndim
      XYZ(i,d,:) = XYZ(i,d,:) - mean(XYZ(i,d,:), Nframes)
   ENDDO
ENDDO
!
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for Benchmark2:   ", F10.6)') Time_End - Time_Start
!
return
end subroutine Benchmark2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Mutual information in the this Experiment:
!!  Y = X+Z
!!  X  uniformly distributed in [-1/2:1/2]
!!  Z1 uniformly distributed in [-a/2:a/2]
!!  Z2 uniformly distributed in [-a/4:a/4]
!!  f_y(y) = convolution(X,Z1,Z2)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Benchmark3(Istart, Iskip, Istop)
implicit none
Integer, intent(in) :: istart, iskip, istop
!
! -- Local variables
integer :: i, d, t, tt, idof
real(sifm_real) :: Time_start, Time_end
!
CALL CPU_time( Time_Start )
!
! -- Initialize Some random number sequences
CALL rndnum_iniall(Natoms*Ndim, 1, 0)
!
tt = 0
DO t = 1, Istop
   IF ( t >= iStart .and. mod(t,iSkip) == 0 ) THEN
        tt = tt + 1
        IF ( tt <= Nframes ) Then
           idof = 0
           do i = 1, Natoms
                do d = 1, Ndim
                   idof = idof + 1
                   xyz(i,d,tt) = gauss_rand(idof)
                enddo
           enddo
        ENDIF
      ENDIF
ENDDO
!
CALL rndnum_clear()
!
DO i = 1, Natoms
   DO d = 1, Ndim
      XYZ(i,d,:) = XYZ(i,d,:) - mean(XYZ(i,d,:), Nframes)
   ENDDO
ENDDO
!
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for Benchmark3:   ", F10.6)') Time_End - Time_Start
!
return
end subroutine Benchmark3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Mutual information in the this Experiment:
!!  Throwing a fair die:
!!  X -> Variable representing the value of the face closeer to us
!!  Y -> Variable representing the value of the top of the die
!!  Possible values are 1, 2, 3, 4, 5, 6
!!  Theoretical value:
!!  I(T,B) = H(B) - H(B|T) = log(6) - log(4) = log(3) - 1
!!  H(B) = 2.584962501; H(B|T) = 2.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Benchmark4(Istart, Iskip, Istop)
implicit none
Integer, intent(in) :: istart, iskip, istop
! -- Local variables
integer :: i, d, t, tt
real(sifm_real):: x_t, y_t
real(sifm_real) :: Time_start, Time_end
!
CALL CPU_time( Time_Start )
!
! -- Initialize Some random number sequences
CALL rndnum_iniall(Natoms*Ndim, 1, 0)
!
Do d=1, Ndim
   tt = 0
   DO t = 1, Istop
      x_t = randf(d)
      y_t = randf(d)
      IF ( t >= iStart .and. mod(t,iSkip) == 0 ) THEN
           tt = tt + 1
           IF ( tt <= Nframes) Then
                xyz(1,d,tt) = x_t
                xyz(2,d,tt) = y_t
           ENDIF
      ENDIF
   ENDDO
ENDDO
!
DO d = 1, Ndim
   DO t = 1, Nframes
           IF ( xyz(1,d,t) < (1.0_sifm_real/6.0_sifm_real) ) THEN
                xyzs(1,d,t) = 1
                IF ( xyz(2,d,t) < (1.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 3
                ELSEIF ( xyz(2,d,t) < (2.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 4
                ELSEIF ( xyz(2,d,t) < (3.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 5
                ELSE
                     xyzs(2,d,t) = 6
                ENDIF
            ELSEIF ( xyz(1,d,t) < (2.0_sifm_real/6.0_sifm_real) ) THEN
                xyzs(1,d,t) = 2
                IF ( xyz(2,d,t) < (1.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 3
                ELSEIF ( xyz(2,d,t) < (2.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 4
                ELSEIF ( xyz(2,d,t) < (3.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 5
                ELSE
                     xyzs(2,d,t) = 6
                ENDIF
            ELSEIF ( xyz(1,d,t) < (3.0_sifm_real/6.0_sifm_real) ) THEN
                xyzs(1,d,t) = 3
                IF ( xyz(2,d,t) < (1.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 1
                ELSEIF ( xyz(2,d,t) < (2.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 2
                ELSEIF ( xyz(2,d,t) < (3.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 5
                ELSE
                     xyzs(2,d,t) = 6
                ENDIF
            ELSEIF ( xyz(1,d,t) < (4.0_sifm_real/6.0_sifm_real) ) THEN
                xyzs(1,d,t) = 4
                IF ( xyz(2,d,t) < (1.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 1
                ELSEIF ( xyz(2,d,t) < (2.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 2
                ELSEIF ( xyz(2,d,t) < (3.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 5
                ELSE
                     xyzs(2,d,t) = 6
                ENDIF
            ELSEIF ( xyz(1,d,t) < (5.0_sifm_real/6.0_sifm_real) ) THEN
                xyzs(1,d,t) = 5
                IF ( xyz(2,d,t) < (1.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 1
                ELSEIF ( xyz(2,d,t) < (2.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 2
                ELSEIF ( xyz(2,d,t) < (3.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 3
                ELSE
                     xyzs(2,d,t) = 4
                ENDIF
            ELSE
                xyzs(1,d,t) = 6
                IF ( xyz(2,d,t) < (1.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 1
                ELSEIF ( xyz(2,d,t) < (2.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 2
                ELSEIF ( xyz(2,d,t) < (3.0_sifm_real/4.0_sifm_real) ) THEN
                     xyzs(2,d,t) = 3
                ELSE
                     xyzs(2,d,t) = 4
                ENDIF

            ENDIF
   ENDDO
ENDDO
!
CALL rndnum_clear()
!
CALL CPU_time( Time_End )
write(*,'("CPU Elapsed time for Benchmark4:   ", F10.6)') Time_End - Time_Start
!
return
end subroutine Benchmark4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Mutual information in the this Experiment:
!!  Flipping a fair Coin:
!!  X -> Variable representing the value of the bottom of the coin
!!  Y -> Variable representing the value of the top of the coin
!!  Head = 1 and Tail = 0
!!  Theoretical value:
!!  I(T,B) = H(B) - H(B|T) = log(2) - 0 = 1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Benchmark5(Istart, Iskip, Istop)
      implicit none
      Integer, intent(in) :: istart, iskip, istop
!
! -- Local variables
      integer :: i, d, t, tt
      real(sifm_real):: x_t, y_t
      real(sifm_real) :: Time_start, Time_end
!
      CALL CPU_time( Time_Start )
!
! -- Initialize Some random number sequences
      CALL rndnum_iniall(2*Ndim, 1, 0)
!
      Do d=1, Ndim
         tt = 0    
         DO t = 1, Istop
            x_t = randf(d)
            y_t = x_t
            IF ( t >= iStart .and. mod(t,iSkip) == 0 ) THEN
                 tt = tt + 1
                 IF ( tt <= Nframes) Then
                      xyz(1,d,tt) = x_t
                      xyz(2,d,tt) = y_t
                 ENDIF 
            ENDIF
         ENDDO
      ENDDO
!
      DO d = 1, Ndim
         DO t = 1, Nframes
            IF ( xyz(1,d,t) > 0.5_sifm_real ) THEN
                 xyzs(1,d,t) = 1
                 xyzs(2,d,t) = 0
            ELSE
                 xyzs(1,d,t) = 0
                 xyzs(2,d,t) = 1
            ENDIF
         ENDDO
      ENDDO
!
      CALL rndnum_clear()
!
      CALL CPU_time( Time_End )
      write(*,'("CPU Elapsed time for Benchmark5:   ", F10.6)') Time_End - Time_Start
!
return
!
end subroutine Benchmark5
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Print out the trajectories
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_xyz()
Implicit None
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='traj.xyz', status='unknown', action='write')
write(nunit, '("SIFM: XYZ Coordinates Format")')
write(nunit, '("SIFM: Copyright Hiqmet Kamberaj")')
DO t = 1, Nframes
   DO i = 1, Natoms
      write(nunit, '(A7, 2I10, 3F12.6)') "FRAME  ", t, i, (xyz(i,d,t), d=1, Ndim) 
   ENDDO
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
!
return
!
end subroutine write_xyz
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
      write(nunit, '(A7, 2I10, 3I3)') "FRAME  ", t, i, (xyzs(i,d,t), d=1, Ndim) 
   ENDDO
ENDDO
write(nunit, '("SIFM: END")')
CLOSE(nunit)
!
return
!
end subroutine write_xyzs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read the trajectories
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xyz()
Implicit None
!  
  integer :: i, d, t, nunit
!
nunit = new_unit()
open(unit=nunit, file='traj.xyz', status='old', action='read')
DO t = 1, Nframes
   read(nunit, *) ( (xyz(i,d,t), d=1, Ndim), i = 1, Natoms ) 
ENDDO
CLOSE(nunit)
!
return
!
end subroutine read_xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Allocate_System(Benchmark)
      implicit none
!
      integer, intent(in) :: benchmark
!
      integer :: ierr
!
!
      IF (.not. Allocated(XYZ)) THEN
          Allocate(XYZ(Natoms, Ndim, Nframes), stat=ierr)
          IF (ierr /= 0) Stop "Error Allocating XYZ"
      ENDIF 
!
      IF (benchmark == 4 .or. benchmark == 5) THEN 
          IF (.not. Allocated(XYZs)) THEN
               Allocate(XYZs(Natoms, Ndim, Nframes), stat=ierr)
               IF (ierr /= 0) Stop "Error Allocating XYZs"
          ENDIF 
      ENDIF
!
      return
end subroutine Allocate_System
!
subroutine deAllocate_System()
      implicit none
!
      integer :: ierr
!  
      IF (Allocated(XYZ)) THEN
	  DEAllocate(XYZ, stat=ierr)
          IF (ierr /= 0) Stop "Error deAllocating XYZ"
      ENDIF
!
      IF (Allocated(XYZs)) THEN
          DEAllocate(XYZs, stat=ierr)
          IF (ierr /= 0) Stop "Error DEAllocating XYZs" 
      ENDIF
!
      return
end subroutine deAllocate_System
!
end module System_class
