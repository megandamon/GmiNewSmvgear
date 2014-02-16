!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiMessageInit_mod
!
! !INTERFACE:
!
      module GmiMessageInit_mod
!
!USES:
      use GmiMessagePassing_mod, only : stopCode
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: initializeMPI
!
! !DESCRIPTION:
! Initialize MPI.
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
! !IROUTINE: initializeMPI
!
! !INTERFACE:
!
      subroutine initializeMPI (msg_comm)
!
      implicit none
!
!
#     include "mpif.h"
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: msg_comm
!
! !DESCRIPTION:
!   This routine initializes the MPI (Message Passing Interface Standard)
!   programming environment.
!
! !LOCAL VARIABLES:
      integer :: procID
      integer :: errlen
      integer :: ii
      integer :: ierr1, ierr2
      character (len=MPI_MAX_ERROR_STRING) :: errstr
!EOP
!-----------------------------------------------------------------------------
!BOC
      msg_comm = MPI_COMM_WORLD

      call MPI_Init (ierr2)

      if (ierr2 /= MPI_SUCCESS) then
        call MPI_Error_String (ierr2, errstr, errlen, ierr1)
        Print '(256a)', (errstr(ii:ii), ii=1,errlen)
        call stopCode (msg_comm, " MPI ERROR ")
      end if

      call MPI_Comm_Rank (msg_comm, procID, ierr2)

      if (ierr2 /= MPI_SUCCESS) then
        call MPI_Error_String (ierr2, errstr, errlen, ierr1)
        Print '(256a)', (errstr(ii:ii), ii=1,errlen)
        call stopCode (msg_comm, " MPI ERROR ")
      end if

      return

      end subroutine initializeMPI
!EOC
!-----------------------------------------------------------------------------

      end module GmiMessageInit_mod
