!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiCommunicationMPI2_mod
!
! !INTERFACE:
!
      module GmiCommunicationMPI2_mod
!
! !USES:
      use GmiMessagePassing_mod, only : stopCode
!
      implicit none
!
#     include "mpif.h"


      integer :: mpi2_offset
      integer :: wingrp, winmax, winobj

      integer (KIND=MPI_ADDRESS_KIND) :: windisp
      integer (KIND=MPI_ADDRESS_KIND) :: winsize

!     -------------------------------------------------------------
!     To switch to F90 allocation, activate the statement below and
!     comment out the two statements below that.
!     -------------------------------------------------------------

      real*8, allocatable :: winarray(:)

!EOP
!--------------------------------------------------------------------------
      contains
!--------------------------------------------------------------------------
!BOP
      subroutine Mpi2_Bord_Init &
             (k1, k2, numSpecies, numDomains, gmi_nborder, map2_u, commuWorld)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: k1, k2, gmi_nborder
      integer, intent(in) :: numSpecies, numDomains
      integer, intent(in) :: commuWorld
      integer, intent(in) :: map2_u(1:,1:,1:)
!
! !LOCAL VARIABLES:
      integer :: ierrwin
      integer :: ilmax, ismax, ksmax
      integer :: is1, is2
      integer :: np
!EOP
!--------------------------------------------------------------------------
!BOC
      ilmax = 0

      do np = 1, numDomains

        is1 = map2_u(2,1,np) - map2_u(1,1,np)
        is2 = map2_u(2,2,np) - map2_u(1,2,np)

        ilmax = Max (ilmax, is1)
        ilmax = Max (ilmax, is2)

      end do

      ismax   = 2 * gmi_nborder * (2 * gmi_nborder + ilmax)
      ksmax   = Max (k2-k1+2, numSpecies)

      winmax  = ismax * ksmax
      winsize = 16 * winmax

!     -------------------------------------------------------------
!     To switch to F90 allocation, activate the statement below and
!     comment out the statement below that.
!     -------------------------------------------------------------

      Allocate (winarray(2*winmax))

!c    call Mpi_Alloc_Mem
!c   &  (winsize, mpi_info_null, pwin, ierrwin)


      call Mpi_Win_Create  &
     &  (winarray, winsize, 8, MPI_INFO_NULL, commuWorld,  &
     &   winobj, ierrwin)

      if (ierrwin /= 0) then
        call stopCode (commuWorld, "Error creating MPI2 window.")
      else
        Write (6,*) 'Created MPI2 window of size ', winsize, ' bytes.'
      end if


      call Mpi_Win_Get_Group  &
     &  (winobj, wingrp, ierrwin)

      if (ierrwin /= 0) then
        call stopCode (commuWorld, "Error creating MPI2 group.")
      end if

      mpi2_offset = 1

      return

      end subroutine Mpi2_Bord_Init
!EOC
!--------------------------------------------------------------------------
!BOP
      subroutine Mpi2_Bord_Final (commu )

      implicit none

      integer :: commu
      integer :: ierrwin
!EOP
!--------------------------------------------------------------------------
!BOC
      call Mpi_Group_Free (wingrp, ierrwin)

      if (ierrwin /= 0) then
        call stopCode (commu, "Error freeing MPI2 group.")
      end if

      call Mpi_Win_Free (winobj, ierrwin)

      if (ierrwin /= 0) then
        call stopCode (commu, "Error releasing MPI2 window.")
      end if

      return

      end subroutine Mpi2_Bord_Final
!EOC
!--------------------------------------------------------------------------
!BOP
      subroutine Mpe_Win_Fence  &
     &  (aaa, iii, jjj, kkk)

      real*8, intent(inout) :: aaa(:)
      integer :: iii, jjj, kkk
!EOP
!--------------------------------------------------------------------------
!BOC
      call Mpi_Win_Fence (iii, jjj, kkk)

      return

      end subroutine Mpe_Win_Fence
!EOC
!--------------------------------------------------------------------------
!BOP
      subroutine Mpe_Win_Wait  (aaa, iii, jjj)

      real*8, intent(inout) :: aaa(:)
      integer :: iii, jjj
!EOP
!--------------------------------------------------------------------------
!BOC
      call Mpi_Win_Wait (iii, jjj)

      return

      end subroutine Mpe_Win_Wait
!EOC
!--------------------------------------------------------------------------
      end module GmiCommunicationMPI2_mod

