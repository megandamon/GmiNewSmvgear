


!  $Id: MAPL_Sort.F90,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

!=============================================================================
!BOP

! !MODULE: MAPL_Sort   -- A utility to sort integers

! !INTERFACE:

module MAPL_SortMod

  implicit none
  private

! !PUBLIC ROUTINES:

public MAPL_Sort

! !DESCRIPTION:
! 
!   {\tt GEOS\_Sort} is a utility to do a quick-sort on integers. The general
!   interface is:
!\bv       
!       subroutine MAPL_Sort(A,B)
!         integer(kind=[4,8]),       intent(INOUT) :: A(:)
!         integer(kind=4), optional, intent(INOUT) :: B(size(A)[,:])
!\ev
!   {\tt GEOS\_Sort} sorts A in ascending order and reorders the columns of B
!   in the same order; i.e., it does the same exchanges to B as were done 
!   to A in sorting it.  If, for example, on input  B contains the ordered integers
!   from 1 to size(A), on output it will contain the positions of the elements of
!   the sorted A in the unsorted A.  

!EOP
!=============================================================================

interface MAPL_Sort
   module procedure SORT0L
   module procedure SORT0S
   module procedure SORT1L
   module procedure SORT1S
   module procedure SORT2L
   module procedure SORT2S
end interface

contains

subroutine SORT0S(A)
  integer(kind=4), intent(INOUT) :: A(:)
  call QSORT0S(A,size(A))
end subroutine SORT0S

subroutine SORT0L(A)
  integer(kind=8), intent(INOUT) :: A(:)
  call QSORT0(A,size(A))
end subroutine SORT0L

subroutine SORT1S(A,B)
  integer(kind=4), intent(INOUT) :: A(:)
  integer(kind=4), intent(INOUT) :: B(:)
  call QSORT2S(A,B,size(A),1)
end subroutine SORT1S

subroutine SORT1L(A,B)
  integer(kind=8), intent(INOUT) :: A(:)
  integer(kind=4), intent(INOUT) :: B(:)
  call QSORT2(A,B,size(A),1)
end subroutine SORT1L

subroutine SORT2S(A,B)
  integer(kind=4), intent(INOUT) :: A(:)
  integer(kind=4), intent(INOUT) :: B(:,:)
  call QSORT2S(A,B,size(A),size(B,2))
end subroutine SORT2S

subroutine SORT2L(A,B)
  integer(kind=8), intent(INOUT) :: A(:)
  integer(kind=4), intent(INOUT) :: B(:,:)
  call QSORT2(A,B,size(A),size(B,2))
end subroutine SORT2L

end module MAPL_SortMod
