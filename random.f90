!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: random.f90,v 1.0 19-03-2018, IBU
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
MODULE TERANDOM_CLASS
  use sifm_kinds
  implicit none
  !   varibles for the new random generators
  !   some common variables:
  !
  !  rngseeds           storage for array of seeds
  !
  !  The following two have to be initialized in the:
  !  where we parse the command line and not in iniall routines!!!
  !
  !  rngchoice    1     the default fixed RNG (clcg)
  !  rngchoice    2     get RNG from new fortran standard:
  !                     call random_number()
  !
  !  rngdistrchoice 0   uniform
  !                 1   gaussian
  !
  !  rngseries   1-100  which CLCG series we want to use (default 1, or is it 2?)
  !
  !  nrand              number of seeds (1,4,8,...)
  !

  !
  !------------------------------------------------------------------CC
  !  Header file for RNG: CLCG (source/util/clcg.src) X. Qian 12/00   C
  !  The variables are:                                               C
  !      Maxgen     512 maximum number of independent streams of      C 
  !                 random number sequences  (This value can be       C
  !                 increased as needed.)                             C
  !      IniSD      1                                                 C
  !      LstSD      2                                                 C
  !      NewSD      3                                                 C
  !                 These three option flags are used to initialize   C
  !                 the RNG with different initial conditions.        C
  !                 By default, initial seeds lcgIg{(i=1, ..., 4),g}  C
  !                 and last seeds  lcgLg{(i=1, ..., 4),g}            C
  !                 (for g=1...Maxgen) are set to the original        C
  !                 seeds, and previous seeds, respectively.          C
  !                 Calls to IniGen(g,stype)                          C
  !                 (where stype is IniSD or LstSD or NewSD)          C
  !                 can reset the seeds for stream g to the initial   C
  !                 values (stype = IniSD), previous values           C
  !                 (stype = LstSD) or new values (stype = NewSD).    C
  !      lcgIg      Initial seed values, dimension 4 by Maxgen        C
  !                 for the four LCGs.                                C
  !      lcgLg      Last seed values, dimension 4 by Maxgen           C
  !      lcgCg      Current seed values, dimension 4 by Maxgen        C
  !                                                                   C
  !     lcgmul(4)   The multipliers for the 4 LCGs                    C
  !     lcgmod(4)   The moduli for the 4 LCGs                         C
  !                 THE MULTIPLIER AND MODULI VALUES                  C
  !                     MUST NOT BE CHANGED                           C
  !                                                                   C
  !    lcgaw(4)     lcgmul{j}^{2^w}      w=41, j=1, ..., 4.           C
  !    lcgavw(4)    lcgmul{j}^{2^(v+w)}, v=31, j=1, ..., 4.           C
  !                 These two arrays are used to generate initial     C
  !                 seeds for the specified Maxgen number of the      C
  !                 streams  with the  initial seeds given by user    C
  !                 or from default values.                           C
  !                                                                   C



  !=======================================================================

  integer,parameter :: Maxgen = 8192, IniSD = 1, LstSD = 2, NewSD = 3
  integer lcgmul(4),lcgmod(4),lcgaw(4),lcgavw(4)
  integer lcgIg(4,Maxgen),lcgLg(4,Maxgen),lcgCg(4,Maxgen)
  data lcgmul/45991,207707,138556,49689/
  data lcgmod/2147483647,2147483543,2147483423,2147483323/ 
  integer,allocatable, dimension(:) :: rngseeds
  integer :: nrand
  integer :: rngchoice

contains
  
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !               Combined Linear Congruential Generator (CLCG)        C
  !  Adapted from Pierre L'Ecuyer & Terry H Andres' C version code.    C
  !  References                                                        C
  ! [1] P. L'Ecuyer and T. H. Andres,                                  C
  !     ``A Random Number Generator Based on the Combination           C
  !     of Four LCGs'', Mathematics and Computers in Simulation,       C
  !     44 (1997), 99--107.                                            C
  ! [2] http://www.iro.umontreal.ca/~lecuyer/                          C
  !                                                                    C
  ! For further information, please contact                            C
  !              Tamar Schlick                                         C  
  !              schlick@nyu.edu                                       C
  ! Converted to FORTRAN by                                            C
  !               Xiaoliang Qian  10/7/99                              C
  !                                                                    C
  !     Fixed the usage of this code by                                C
  !               Milan Hodoscek  6/6/04                               C
  !                                                                    C
  !     NOTE: The code needs to be improved for parallel               C
  !                                                                    C
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !  This pseudorandom number generator combines 4 linear              C
  !  congruential generators (LCGs) to get the long period             C
  !  of about 2^121 for the resulting sequence. This sequence has      C
  !  passed many statistical tests of `randomness'.                    C
  !  Essentially, the four LCGs defined as                             C
  !                                                                    C
  !      X{j,n} = lcgmul{j}X{j,n-1} mod lcgmod{j}                (1)   C
  !                            j=1,2,3,4                               C
  !                            n=1,2,3, ...                            C
  !                                                                    C
  !  with lcgmod{1}= 2147483647, lcgmod{2}=2147483543,                 C
  !       lcgmod{3}=2147483423, lcgmod{4}=2147483323                   C
  !  and lcgmul{1}=45991, lcgmul{2}=207707, lcgmul{3}=138556,          C
  !  lcgmul{4}=49689.                                                  C
  !                                                                    C
  !  The construct                                                     C
  !      Z{n} = (Sum [(-1)^{j+1}*X{j,n}/lcgmod{j}]) mod 1   (2)        C
  !                                                                    C
  ! for n=1,2,3, ... is then a uniformly distributed random sequence   C
  ! in (0,1). It can be proved that the LCG corresponding to the       C
  ! combined generator has modulus, multiplier, and period length of   C
  !      21267641435849934371830464348413044909,                       C
  !      5494569482908719143153333426731027229,                        C
  !      (2^{31}-2)(2^{31}-106)(2^{31}-226)(2^{31}-326) ~ 2^{121},     C
  !  respectively.                                                     C
  !                                                                    C
  !  The default initial seed is the vector {11111111, 22222222,       C
  !  33333333, 44444444} and can be changed by calling SetIniSD        C
  !  after calling CLCGInit.                                           C
  !                                                                    C
  !  This RNG can be used under parallel conditions to give            C 
  !  independent random number sequence when each processor            C
  !  calls with different stream number g (e.g., RANDOM(g)).           C
  !                                                                    C
  !  To use these RNG routines, the user should proceed as follows:    C
  !                                                                    C
  !  1. Call routine CLCGInit() to initialize all Maxgen (100) streams C
  !     using four default initial seeds.                              C
  !  2. [Optional] Call SetiniSD(sd) with desired seed array           C
  !     (4 values) to override the default values specified in Init(). C
  !  3. Call function RANDOM(k), where k is an integer from 1 to 100   C
  !     specifying the stream number. For parallel codes, k can be set C
  !     to a processor id number.                                      C
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
 
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
  function RANDF(ig) result(random_rtn)
    !------------------------------------------------------------------------C
    !   Return a double precision uniformly distributed random number in     C
    !   (0,1) from the gth stream and reset the current seed Cg accordingly  C
    !   (i.e., using one of the 100 initial seed vectors generated in the    C
    !   SetiniSD routine).                                                   C
    ! Converted to FORTRAN by                                                C
    !               Xiaoliang Qian  10/7/99                                  C
    !
    !
    implicit none
!    
    integer, intent(in) :: ig
    real(sifm_real) :: random_rtn
!
    integer :: k,s,j,g
    real(sifm_real) :: u(4), rng1
    integer  :: dv(4),mv(4)
    data dv/46693,10339,15499,43218/ 
    data mv/25884, 870,3979,24121/ 
    data u/4.65661287524579692d-10,-4.65661310075985993d-10, &
         4.65661336096842131d-10,-4.65661357780891134d-10/
    
    !
    ! system provided RNG
    if(rngchoice == 2) then
       call random_number(rng1)
       random_rtn = rng1
       return
    endif
    g = ig
    if (g  <=  1)  g = 1
    g= mod(g-1,Maxgen) + 1
    RANDOM_rtn = rzero
    do j = 1,4
       s = lcgCg(j,g)
       k = s/dv(j)
       s = lcgmul(j) * (s - k * dv(j)) - k * mv(j)
       if (s  <  0) s = s + lcgmod(j)
       lcgCg(j,g) = s
       RANDOM_rtn = RANDOM_RTN + u(j) * s
       if (RANDOM_RTN  <  rzero)  RANDOM_RTN = RANDOM_RTN + one
       if (RANDOM_RTN  >=  one)   RANDOM_RTN = RANDOM_RTN - one
    enddo
    return
  end function RANDF

    !------------------------------------------------------------------------C
    !  This optional routine uses the input seed value s for stream g        C
    !  instead of the default settings (routine SetiniSD).                   C
    ! Converted to FORTRAN by                                                C
    !               Xiaoliang Qian  10/7/99                                  C
  subroutine  SetSeed (ig,s)
!
    integer, intent(in) :: ig, s(4)
!
    integer :: j,g
!    
    g = ig
    if (g  <=  1)  g = 1
    g= mod(g-1,Maxgen) + 1    
    do j = 1,4          
       lcgIg(j,g) = s(j)
    enddo
    call IniGen (g,IniSD)               
  end subroutine SetSeed
  
  subroutine  GetSeed (ig,s)
    !------------------------------------------------------------------------C
    !  This optional routine returns current seed value s for stream g       C
    ! Converted to FORTRAN by                                                C
    !               Xiaoliang Qian  10/7/99                                  C
    
!    
    integer, intent(in) :: ig
    integer, intent(out) :: s(4)
!
    integer :: j,g
!    
    g = ig
    if (g  <=  1)  g = 1
    g= mod(g-1,Maxgen) + 1
!
    do j = 1,4          
       s(j)= lcgCg(j,g)
    enddo
    return
  end subroutine GetSeed

    !------------------------------------------------------------------------
    !  This optional routine resets the gth stream so that the initial seed
    !  is either the original initial seed (if stype = IniSD) or the last
    !  seed (if stype = 3).
    ! Converted to FORTRAN by
    !               Xiaoliang Qian  10/7/99
  subroutine IniGen (g,stype)
!    
    integer, intent(in) :: g,stype
!
    integer :: j,ig
!    
    ig=g
    if (ig  <=  1)  ig = 1
    ig= mod(ig-1,Maxgen) + 1
    do j = 1,4
       if (stype == IniSD) then
          lcgLg(j,ig) = lcgIg(j,ig)
       else
          if (stype == NewSD)  &
               lcgLg(j,ig) = MulMod (lcgaw(j),lcgLg(j,ig),lcgmod(j)) 
       endif
       lcgCg(j,ig) = lcgLg(j,ig)
    enddo
    return
  end subroutine IniGen

!--------------------------------------------------------------------C
!  Return s*t mod M. All numbers out of range are truncated          C
!  before the mod operation to avoid overflow. See Park & Miller,    C
!  Comm. ACM 31, 1192 (1988) for this multiplication procedure.      C
! Converted to FORTRAN by                                            C
!               Xiaoliang Qian  10/7/99                              C
!--------------------------------------------------------------------C
  integer function MulMod (s,t,M)
    integer, intent(inout) :: s,t,M
!
    integer, parameter :: H = 32768
    integer :: S0,S1,q,qh,rh,k 
!
    if (s  <  0) s = s + M
    if (t  <  0) t = t + M
    if (s  <  H) then
       S0 = s
       MulMod = 0
    else
       S1 = s / H
       S0 = s - H * S1
       qh = M / H
       rh = M - H * qh

       if (S1  >=  H) then
          S1 = S1 - H
          k = t / qh
          MulMod = H * (t - k * qh) - k * rh

10        if (MulMod  <  0) then
             MulMod = MulMod + M      
             goto 10
          endif

       else
          MulMod = 0
       endif

       if (S1  /=  0) then
          q = M / S1
          k = t / q
          MulMod = MulMod - k * (M - S1 * q)
          if (MulMod  >  0) MulMod = MulMod - M
          MulMod = MulMod + S1 * (t - k * q) 

20        if (MulMod  <  0) then
             MulMod = MulMod + M      
             goto 20
          endif
       endif

       k = MulMod / qh
       MulMod = H * (MulMod - k * qh) - k * rh

30     if (MulMod  <  0) then
          MulMod = MulMod + M      
          goto 30
       endif
    endif

    if (S0  /=  0) then
       Q = M / S0
       k = t / q                         
       MulMod = MulMod - k * (M - S0 * q)
       if (MulMod  >  0) MulMod = MulMod - M
       MulMod = MulMod + S0 * (t - k * q) 
40     if (MulMod  <  0)  then
          MulMod = MulMod + M      
          goto 40
       endif
    endif
    return
  end function MulMod

  !
  !     <*> the following code added by rm venable, FDA Biophysics Lab <*>
  !         Modified by SRD 1/27/91.
  !
  FUNCTION GAUSS_rand(g) result(RND)
    !-----------------------------------------------------------------------
    !     This function routine generates a Gaussian random
    !     deviate of 0.0 mean and standard deviation RF.
    !     The algorithm from Box and Muller.
    !
    !      Bernard R. Brooks   January, 1988
    !
    !     copied from dynlng and renamed <*>  rm venable  <*>  7 dec 1989
    !
    implicit none
    integer, intent(in) :: g
    real(sifm_real) :: RND
!
    real(sifm_real)       :: rng1,rng2
    logical, save        :: set
    real(sifm_real), save :: gset
    real(sifm_real)       :: v1,v2,r,fac
    data set /.TRUE./
!
    IF (set) THEN
       r = 2.0_sifm_real
       do while (r > 1.0_sifm_real)
          if (rngchoice == 1) then
              rng1 = randf(g)
              rng2 = randf(g)
          elseif (rngchoice == 2) then
              call random_number(rng1)
              call random_number(rng2)
          ELSE
              Stop "No choice made for random numbers"
          endif
          v1 = 2.0_sifm_real*rng1 - 1.0_sifm_real
          v2 = 2.0_sifm_real*rng2 - 1.0_sifm_real
          r = v1*v1 + v2*v2
       enddo
       fac = sqrt(-2.0_sifm_real * log(r)/r)
       gset = v1 * fac
       set = .FALSE.
       Rnd = v2 * fac
    ELSE
       set = .TRUE.
       Rnd = gset
    ENDIF
    RETURN
  END FUNCTION GAUSS_Rand

subroutine randf_vec(I,N,r)
  implicit none
  integer, intent(in) :: I,N
  real(sifm_real), intent(out) :: r(N)
  integer :: k
  do k=1, N
     r(k)=randf(I)
  enddo
  return
end subroutine randf_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rndnum_iniall(nrnd, rngc, mynod)
    implicit none
!
    integer, intent(in) :: Nrnd, mynod, rngc
!
    integer(sifm_int8) :: count, rate, maxcount
    integer(sifm_int4) :: count4, rate4, maxcount4
    integer :: ierr
!
    rngchoice = rngc       ! CLCG default now
    nrand = Nrnd  
    IF (Nrand > MaxGen) THEN
        write(*,*) Nrand, MaxGen
        Stop "Increase the Max number of seeds"
    ENDIF    
    IF (.not. Allocated(RNGseeds)) THEN
	ALLOCATE( RNGseeds(Nrand), stat=ierr )
        if (ierr /= 0) stop 'Error allocating RNGSeeds'
    ENDIF
!
    if(rngchoice == 1) then
          ! get the seeds from system time by default
          ! NOTE: no big numbers for CLCG!!
          call system_clock(count4,rate4,maxcount4)
          rngseeds(1:nrand) = count4 + 50000*mynod
          call CLCGInit
    elseif(rngchoice == 2) then ! this is not needed here for now
          call random_seed(nrand)
          ! get the seeds from system time by default
          call system_clock(count,rate,maxcount)
          rngseeds(1:nrand) = count + 50000*mynod
          call random_seed(put=rngseeds)
    endif
    return
  end subroutine rndnum_iniall

  subroutine rndnum_clear()
    implicit none
!
    integer :: ierr
!  
    IF (Allocated(rngseeds)) THEN
	DEAllocate(rngseeds, stat=ierr)
        IF (ierr /= 0) Stop "Error deAllocating RNGSeeds"
    ENDIF
!
    return
  end subroutine rndnum_clear

    !------------------------------------------------------------------------C
    !   Initialize the RNG with seed values in vector sd of dimension 4.     C
    !   Each such initial vector generates one stream of random variates     C
    !   combining the 4 LCG sequences with resulting sequence length         C
    !   2^{121}.                                                             C
    !   The lcgaw and lcgavw arrays of dimension 4 have default values       C
    !       lcgaw{j} = lcgmul{j}^{2^31} mod lcgmod{j}                        C
    !   and lcgavw{j} = lcgmul{j}^{2^41} mod lcgmod{j},                      C
    !   for j=1, ..., 4 corresponding to the 4 LCGs.                         C
    ! Converted to FORTRAN by                                                C
    !               Xiaoliang Qian  10/7/99                                  C
 subroutine CLCGInit
    implicit none
!
!
    integer, dimension(4) :: sd
    integer :: j,i
    integer, parameter :: v = 31, w = 41
    data sd/11111111,22222222,33333333,44444444/
!
    sd(1:4) = rngseeds(1:4)
    !
    do j = 1,4
       lcgaw(j) = lcgmul(j)
       do  i = 1,w
          lcgaw(j) = MulMod (lcgaw(j),lcgaw(j),lcgmod(j))
       enddo
       lcgavw(j) = lcgaw(j)
       do i = 1,v
          lcgavw(j) = MulMod (lcgavw(j),lcgavw(j),lcgmod(j))
       enddo
    enddo
    call SetiniSD (sd)
  end subroutine CLCGInit

    !------------------------------------------------------------------------C
    !  Set initial seed values for all 100 (= Maxgen, defined in clcg.f90)   C
    !  streams using the initial seeds from the first stream.                C
    ! Converted to FORTRAN by                                                C
    !               Xiaoliang Qian  10/7/99                                  C  
  subroutine SetiniSD (s)
    implicit none
!
    integer, intent(in) :: s(4)
!
    integer :: g, j 
    
    do j = 1, 4
       lcgIg(j,1) = s(j)
    enddo
    call IniGen (1,IniSD)
    do g = 2, Maxgen
       do  j = 1,4
          lcgIg(j,g) = MulMod (lcgavw(j),lcgIg(j,g-1),lcgmod(j))
       enddo
       call IniGen (g,IniSD)
    enddo
    RETURN
  end subroutine SetiniSD

end module teRandom_Class
