
      module GmiFlush_mod

      implicit none

      private
      public  :: GmiFlush

!-----------------------------------------------------------------------------

      CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   GmiFlush
!
! DESCRIPTION
!   This routine is a port-a-potty; it flushes anywhere.
!
! ARGUMENTS
!    unit_no : logical unit number to flush
!
!-----------------------------------------------------------------------------

      subroutine GmiFlush (unit_no)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: unit_no

!     ----------------
!     Begin execution.
!     ----------------

      call Flush (unit_no)

      return

      end subroutine GmiFlush

      end module GmiFlush_mod
