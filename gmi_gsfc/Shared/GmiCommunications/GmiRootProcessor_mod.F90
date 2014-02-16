!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiRootProcessor_mod
!
! !INTERFACE:
!
      module GmiRootProcessor_mod
!
! !USES:
      use GmiMessagePassing_mod, only : stopCode
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: gmiRootProcessor
!

#     include "mpif.h"
!
! !DESCRIPTION:
! Finds out if a given processor is the root.
!
! !AUTHOR:
!  Jules Kouatchou, NASA GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-----------------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gmiRootProcessor
!
! !INTERFACE:
!
      function gmiRootProcessor (commuWorld) result(root)
!
      implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: commuWorld
!
! !RETURN VALUE:
  logical             :: root
!
! !DESCRIPTION:
!   This routine determines which processor is the root processor.
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: procID
!EOP
!-----------------------------------------------------------------------------
!BOC

      root = .false.
      call MPI_COMM_RANK(commuWorld, procID, ierr)

      if (procID == 0) root = .true.

      return

      end function gmiRootProcessor
!EOC
!-----------------------------------------------------------------------------

      end module GmiRootProcessor_mod
