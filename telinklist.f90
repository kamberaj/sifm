!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Id: telinklist.f90,v 1.0 19-03-2018, IBU
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
Module TELinkList_class
use sifm_kinds
implicit none

type :: StateElem
  integer :: NWCounts = 1
  Integer :: State
  character(len=250) :: word
  type (StateElem), pointer :: Next
end type StateElem

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine LinkList(head,tail,temp)
      implicit none
      type (StateElem), intent(inout), pointer :: head,tail
      character(len=*), intent(in)  :: temp
! 
      type(StateElem), pointer :: ptr
      integer :: addword,istat,current_state
  
      if (.not. associated(head)) then
          allocate(head,stat=istat)
          tail=>head
          nullify(tail%next)
          tail%word = temp
          tail%state = 1
      else
          ptr=>head
          search: DO
             if (.not. associated(ptr%next)) then
                 if (compareStrings(temp, tail%word) == 1) then
                     addword=0
                     tail%NWcounts=tail%NWcounts + 1
                 else
                     addword=1
                 endif
                 exit search
             else if ( compareStrings(ptr%word,temp) == 1) then
                 ptr%NWcounts = ptr%NWcounts + 1
                 addword=0
                 exit search
             endif
             ptr=>ptr%next
          end do search
          if (addword==1) then
              allocate(tail%next,stat=istat)
              current_state = tail%State
              tail=>tail%next
              nullify(tail%next)
              tail%word = temp
              tail%State = current_state + 1
         endif
      endif
      end subroutine LinkList

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!{ Function to compare two strings }
! Output:
!    0 if they are not equal
!    1 if they are equal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compareStrings(a,b) result(r)
character(len=*), intent(in) :: a,b
integer :: r
!{Local variables}
integer :: n1,n2,i
n1 = len_trim(a)
n2 = len_trim(b)
if (n1 /= n2) stop "strings should be equal in length"
r = 1
do i = 1, n1
if ( a(i:i) /= b(i:i) ) then
r=0
return
endif
enddo
return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine LocalLinkList(head,tail,temp,visited_state)
      implicit none
      type (StateElem), intent(inout), pointer :: head,tail
      character(len=*), intent(in)  :: temp
      integer, intent(out) :: visited_state
! 
      type(StateElem), pointer :: ptr
      integer :: addword,istat,previous_state
  
      if (.not. associated(head)) then
          allocate(head,stat=istat)
          tail=>head
          nullify(tail%next)
          tail%word = temp
          tail%state = 1
          visited_state = tail%state
      else
          ptr=>head
          search: DO
             if (.not. associated(ptr%next)) then
                 if (temp == tail%word) then
                     addword=0
                     tail%NWcounts=tail%NWcounts + 1
                     visited_state = tail%state
                 else
                     addword=1
                 endif
                 exit search
             else if ( ptr%word == temp ) then
                 ptr%NWcounts = ptr%NWcounts + 1
                 visited_state = ptr%state
                 addword=0
                 exit search
             endif
             ptr=>ptr%next
          end do search
          if (addword==1) then
              allocate(tail%next,stat=istat)
              previous_state = tail%State
              tail=>tail%next
              nullify(tail%next)
              tail%word = temp
              tail%State = previous_state + 1
              visited_state = tail%state
         endif
      endif
      return
      end subroutine LocalLinkList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine freeList(head)
    IMPLICIT NONE
    type( StateElem ), pointer :: head
    type( StateElem ), pointer :: current, next
    
    current => head
    do while(associated(current))
       next => current%next
       deallocate(current)
       nullify(current)
       current => next    
    enddo    
   END subroutine freeList

   FUNCTION DeleteList(head)
    IMPLICIT NONE
    type( StateElem ), pointer :: head
    type( StateElem ), pointer :: DeleteList
    type( StateElem ), pointer :: h
    IF ( ASSOCIATED(head) ) THEN
       h => head
       head => head%next
       deallocate(h)
    END IF
    DeleteList => head
   END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function SymbolicShannonEntropy(head) result(Entropy)
      type (StateElem), intent(in), pointer :: Head
      real(sifm_real) :: Entropy
!
      type (StateElem), pointer :: ptr
      integer :: Total
      real(sifm_real) :: prob
!  
      ptr => head
      total = 0
      DO WHILE ( associated(ptr) )
         IF ( ptr%NWcounts > 0 ) THEN
              total = total + ptr%NWcounts
         ENDIF
         ptr => ptr%next
      END DO
      IF (total > 0) THEN 
          ptr => head
          Entropy = 0.0_sifm_real
          DO WHILE ( associated(ptr) )
             IF ( ptr%NWcounts > 0 ) THEN
                  prob = real(ptr%NWcounts, sifm_real) / real(total, sifm_real)
                  Entropy = Entropy  - prob * log( prob )
             ENDIF
             ptr => ptr%next
          END DO
      ELSE
          STOP "<SymbolicShannonEntropy> Error in calculating Shannon entropy"
      ENDIF
      return
  end function SymbolicShannonEntropy
!
  subroutine SymbolicJointProbabilityTraj(head, total, prob, frequency)
      type (StateElem), intent(in), pointer :: Head
      real(sifm_real), intent(out) :: prob(:)
      integer, intent(out) :: frequency(:)
      integer, intent(in) :: total
!
      type (StateElem), pointer :: ptr
      integer :: istate
!     
      IF (total > 0) THEN
          ptr => head
          Prob = 0.0_sifm_real
          DO WHILE ( associated(ptr) )
            istate = ptr%state
IF ( ptr%NWcounts > 0 ) THEN
         prob(istate) = real(ptr%NWcounts, sifm_real) / real(total, sifm_real)
                 Frequency(istate) = ptr%NWcounts
ENDIF
            ptr => ptr%next
          END DO
      ENDIF
      return
  end subroutine SymbolicJointProbabilityTraj
!
  function getTotalStates(head) result(total)
      type (StateElem), intent(in), pointer :: Head
      integer :: total      
!
      type (StateElem), pointer :: ptr
!  
      ptr => head
      total = 0
      DO WHILE ( associated(ptr) )
            IF ( ptr%NWcounts > 0 ) THEN
           total = total + ptr%NWcounts
            ENDIF
            ptr => ptr%next
      END DO
      return
  end function
!
  function CountStates(head) result(Nstates)
      type (StateElem), intent(in), pointer :: Head
      integer :: Nstates
!
      type (StateElem), pointer :: ptr
!  
      ptr => head
      Nstates = 0
      DO WHILE ( associated(ptr) )
            IF (ptr%NWcounts > 0) THEN
    Nstates = max(Nstates,ptr%State)
            ENDIF
            ptr => ptr%next
      END DO
      return
  end function
!
end Module TELinkList_class
