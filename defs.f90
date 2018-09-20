!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: defs.f90,v 1.0 19-03-2018, IBU
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


module SIFM_KINDS 
implicit none

 ! Kind constants for system-independent precision specification.
 ! H. Kamberaj
    integer, parameter :: sifm_real4 = kind(0.)   !-- real single precision
    integer, parameter :: sifm_real8 = kind(0.d0) !-- real double precision

    !-- Kinds for specified real precisions.
    !integer, parameter :: sifm_real4 = selected_real_kind(6, 30)   !-- 8 digits
    !integer, parameter :: sifm_real8 = selected_real_kind(12,30)  !-- 16 digits

    !-- Kinds for specified integer ranges.
    integer, parameter :: sifm_int1 = selected_int_kind(2)  !-- 99 max
    integer, parameter :: sifm_int2 = selected_int_kind(4)  !-- 9,999 max
    integer, parameter :: sifm_int4 = selected_int_kind(8)  !-- 99,999,999 max
    integer, parameter :: sifm_int8 = selected_int_kind(16) !-- 9,999,999,999,999,999 max

    !-- Kind for working real precision.
    integer, parameter :: sifm_real = sifm_real8

    !-- Big constants. These are less than huge so that we have some
    !-- room for comparisons and other arithmetic without overflowing.
    real(sifm_real), parameter :: rzero = 0.0_sifm_real
    real(sifm_real), parameter :: one  = 1.0_sifm_real
    real(sifm_real), parameter :: toBits  = one/log(2.0_sifm_real)
    real(sifm_real), parameter :: half = 0.5_sifm_real
    real(sifm_real), parameter :: M_PI = 3.141592653589793238462643383279502884197_sifm_real
    real(sifm_real), parameter :: M_2PI = 6.283185307179586476925286766559005768394_sifm_real
    real(sifm_real), parameter :: ToDeg = 180.0_sifm_real/M_PI
    real(sifm_real), parameter :: LOG_ZERO=-30.0_sifm_real
    real(sifm_real), parameter :: max_exp = maxexponent(rzero)
    real(sifm_real), parameter :: min_exp = minexponent(rzero)
    real(sifm_real), parameter :: kB = 0.00198614_sifm_real
    REAL(sifm_real), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sifm_real
    real(sifm_real), parameter :: PRTINY = tiny(rzero) * 16.0_sifm_real
    real(sifm_real), parameter :: PRBIG  = huge(sifm_real) / 16.0_sifm_real
    integer, parameter         :: PIBIG  = huge(sifm_int8)/16

    character(len=1), dimension(26) :: Symbol = (/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',  &
                                                  'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',  &
                                                  's', 't', 'u', 'v', 'w', 'x', 'y', 'z' /)

end module SIFM_KINDS

