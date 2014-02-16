!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiBroadcast_mod
!
! !INTERFACE:
!
      module GmiBroadcast_mod
!
! !USES:
      use GmiMessagePassing_mod, only : broadcastReal8
      use GmiPrintError_mod     , only : GmiPrintError
!
      implicit none
!
! !PUBLIC MEMBERS FUNCTIONS:
      private
      public   :: Gmi_Pole_Broadcast
!
#     include "mpif.h"
!
! !DESCRIPTION:
! Contains a routine for broadcasting information to processors
! having regions around the poles.
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
! !IROUTINE: Gmi_Pole_Broadcast
!
! !INTERFACE:
!
      subroutine Gmi_Pole_Broadcast (count, buffer, &
     &   ju1, j2, ju1_gl, j2_gl, commu_npole, commu_spole)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ju1, j2, ju1_gl, j2_gl
      integer, intent(in) :: commu_spole
      integer, intent(in) :: commu_npole
      integer, intent(in) :: count       ! length of buffer array
      real*8  :: buffer(count)           ! array of size count of data to broadcast
!
! !DESCRIPTION:
!   Broadcasts information to a group using message passing.
!   The information broadcast initially resides in the root process, but
!   resides in all processes of the group at the end of the routine.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: commu
!EOP
!-----------------------------------------------------------------------------
!BOC

      if (((ju1 == ju1_gl) .and. (j2 == j2_gl)) .or.  &
     &    ((ju1 /= ju1_gl) .and. (j2 /= j2_gl))) then

        err_msg = 'Pole problem in Gmi_Pole_Broadcast.'
        call GmiPrintError (err_msg, .true., 2, ju1, j2, 0, 0.0d0, 0.0d0)

      else if (ju1 == ju1_gl) then

        commu = commu_spole

      else if (j2 == j2_gl) then

        commu = commu_npole

      end if

      call broadcastReal8 (buffer, count, 0, commu)

      return

      end subroutine Gmi_Pole_Broadcast
!EOC
!---------------------------------------------------------------------------
      end module GmiBroadcast_mod
