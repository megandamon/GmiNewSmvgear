!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiMessageFinal_mod
!
! !INTERFACE:
!
      module GmiMessageFinal_mod
!
! !USES:
      use GmiMessagePassing_mod, only : synchronizeGroup, stopCode
!
      implicit none
!
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: finalizeMPI
!
! !DESCRIPTION:
! Finalizes MPI.
!
! !AUTHOR:
!  John Tannahill, LLNL , jrt@llnl.gov
!  Jules Kouatchou, NASA GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-----------------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeMPI
!
! !INTERFACE:
!
      subroutine finalizeMPI (msg_comm)
!
      implicit none
!

#     include "mpif.h"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: msg_comm
!
! !DESCRIPTION:
!   This routine closes up MPI, if it is being used.
!
! !LOCAL VARIABLES:
      integer :: errlen
      integer :: ierr1, ierr2
      integer :: ii
      character (len=MPI_MAX_ERROR_STRING) :: errstr
!EOP
!-----------------------------------------------------------------------------
!BOC

!     ===============
      call synchronizeGroup(msg_comm)
!     ===============

!     =================
      call MPI_Finalize (ierr1)
!     =================

      if (ierr1 /= MPI_SUCCESS) then
        call MPI_Error_String (ierr1, errstr, errlen, ierr2)
        Print '(256a)', (errstr(ii:ii), ii=1,errlen)
        call stopCode (msg_comm, " MPI ERROR ")
      end if

      return

      end subroutine finalizeMPI
!EOC
!-----------------------------------------------------------------------------

      end module GmiMessageFinal_mod
