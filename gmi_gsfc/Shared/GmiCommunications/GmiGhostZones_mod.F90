!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiGhostZones_mod
!
! !INTERFACE:
!
      module GmiGhostZones_mod
!
! !USES:
      use GmiMessagePassing_mod, only : receiveReal8, sendReal8, sendReceiveReal8
      use GmiMessagePassing_mod, only : stopCode
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public   :: Msg3d_Ew
      public   :: Msg3d_Ns

! !DEFINED PARAMETERS:
      integer, parameter :: msg3d_outer_flag = 1
!
! !DESCRIPTION:
! Routine for inter-processor communications around ghost zones.
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
! !IROUTINE: Msg3d_Ew
!
! !INTERFACE:
!
      subroutine Msg3d_Ew  &
     &  (array, ilo, ihi, jlo, jhi, klo, khi, nborder, eb, wb,  &
     &   lpe, lpw, ewperiod, weperiod, msg_id,  &
     &   pass_to_east, pass_to_west, pass_from_east, pass_from_west, &
     &   msg_comm)
!
      implicit none
!
#     include "gem_msg_numbers.h"
#     include "GmiParameters.h"
!
! !INPUT PARAMETERS:
               ! dimensions of array to be passed 
      integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
               ! number of ghost zones defining the border frame around the
               ! subdomain (e.g., nborder=1 for second order, 2 for fourth order)
      integer, intent(in) :: nborder
               ! flag to indicate whether a subdomain is on the east boundary
      integer, intent(in) :: eb
               ! flag to indicate whether a subdomain is on the west boundary
      integer, intent(in) :: wb
               ! local processor for subdomain to the east of this one
      integer, intent(in) :: lpe
               ! local processor for subdomain to the west of this one
      integer, intent(in) :: lpw
               ! flag to indicate east-west periodicity or islands
      integer, intent(in) :: ewperiod
               ! flag to indicate west-east periodicity or islands
      integer, intent(in) :: weperiod
               ! message identifier
      integer, intent(in) :: msg_id
               ! message pasing communicator
      integer, intent(in) :: msg_comm
!
! !INPUT/OUTPUT PARAMETERS
               ! a 3D array, whose borders are to be passed
      real*8 , intent(inOut) ::  array(ilo:ihi, jlo:jhi, klo:khi)
               ! working arrays to carry the messages
      real*8 , intent(inOut) :: pass_to_east  (jlo:jhi, klo:khi, nborder)
      real*8 , intent(inOut) :: pass_to_west  (jlo:jhi, klo:khi, nborder)
      real*8 , intent(inOut) :: pass_from_east(jlo:jhi, klo:khi, nborder)
      real*8 , intent(inOut) :: pass_from_west(jlo:jhi, klo:khi, nborder)
!
! !DESCRIPTION:
! This routine is used to copy, transmit, and receive data from neighboring
! subdomain processors into the ghost zones on the east and west sides.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer :: easterly, westerly
      integer :: ew_x, we_x
      integer :: ierrwin
      integer :: indwin
      integer :: iouter
      integer :: irc
      integer :: j, k
      integer :: lperror
      integer :: sow_we
      integer, save :: ewperiod_save, weperiod_save
      integer, save :: lpe2, lpw2
      integer, save :: lpe_save, lpw_save
      integer, save :: num_send
      integer, save :: send_grp
      integer, save :: rank_send(2)
!
!EOP
!------------------------------------------------------------------------------
!BOC

!     =========================================
      IOUTLOOP: do iouter = 1, msg3d_outer_flag
!     =========================================

!       -------------------------------------------------------------------
!       Generate unique message numbers for this variable.
!       NOTE:  message numbers will be unique if msg_id is chosen to be
!       incremented by at least 4 for all variables which use this routine.
!       -------------------------------------------------------------------

        ew_x = msg_id + EW_ESM
        we_x = msg_id + WE_ESM

!       ----------------------------------------------------------
!       Gather east-west inner borders (3D equations).
!       These loops pack up 3D data from the outermost real zones
!       of this processor's subdomain into a message to be sent to
!       the neighboring subdomain.
!       ----------------------------------------------------------

        easterly = 0
        westerly = 0

        if ((eb == 0) .or. (ewperiod /= 0)) then

          easterly = 1

          do irc = 1, nborder
            do k = klo, khi
              do j = jlo, jhi
                pass_to_east(j,k,irc) = array(ihi-nborder+1-irc,j,k)
              end do
            end do
          end do

        end if

        if ((wb == 0) .or. (weperiod /= 0)) then

          westerly = 1

          do irc = 1, nborder
            do k = klo, khi
              do j = jlo, jhi
                pass_to_west(j,k,irc) = array(ilo+nborder-1+irc,j,k)
              end do
            end do
          end do

        end if

        sow_we = (jhi - jlo + 1) * (khi - klo + 1) * nborder

!         ---------------------------------------------------
!         Transmit from west to east for multiprocessor runs.
!         ---------------------------------------------------

          if ((eb /= 0) .and. (ewperiod == 0)) then

            call receiveReal8 (lpw, sow_we, we_x, pass_from_west, msg_comm)

          else if ((wb /= 0) .and. (weperiod == 0)) then

            call sendReal8 (lpe, sow_we, we_x, pass_to_east,   msg_comm)

          else

            call sendReceiveReal8  &
     &        (lpw, sow_we, we_x, pass_to_east,  &
     &         lpe, sow_we, we_x, pass_from_west, msg_comm)

          end if

!         ---------------------------------------------------
!         Transmit from east to west for multiprocessor runs.
!         ---------------------------------------------------

          if ((wb /= 0) .and. (weperiod == 0)) then

            call receiveReal8 (lpe, sow_we, ew_x, pass_from_east, msg_comm)

          else if ((eb /= 0) .and. (ewperiod == 0)) then

            call sendReal8 (lpw, sow_we, ew_x, pass_to_west,   msg_comm)

          else

            call sendReceiveReal8  &
     &        (lpe, sow_we, ew_x, pass_to_west,  &
     &         lpw, sow_we, ew_x, pass_from_east, msg_comm)

          end if

!       -------------------------------------------------------------
!       Scatter east-west outer borders (3D equations).
!       These loops unpack the message received into the ghost zones.
!       -------------------------------------------------------------

        if ((wb == 0) .or. (weperiod /= 0)) then

            do irc = 1, nborder
              do k = klo, khi
                do j = jlo, jhi
                  array(nborder+ilo-irc,j,k) = pass_from_west(j,k,irc)
                end do
              end do
            end do
        end if

        if ((eb == 0) .or. (ewperiod /= 0)) then

            do irc = 1, nborder
              do k = klo, khi
                do j = jlo, jhi
                  array(ihi-nborder+irc,j,k) = pass_from_east(j,k,irc)
                end do
              end do
            end do
        end if

!     ===============
      end do IOUTLOOP
!     ===============

      return

      end subroutine Msg3d_Ew
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Msg3d_Ns
!
! !INTERFACE:
!
      subroutine Msg3d_Ns  &
     &  (array, ilo, ihi, jlo, jhi, klo, khi, nborder, nb, sb,  &
     &   lpn, lps, nsperiod, snperiod, msg_id,  &
     &   pass_to_north, pass_to_south, pass_from_north, pass_from_south, &
     &   msg_comm)
!
      implicit none
!
#     include "gem_msg_numbers.h"
#     include "GmiParameters.h"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
               ! number of ghost zones defining the border frame around the
               ! subdomain (e.g., nborder=1 for second order, 2 for fourth order)
      integer, intent(in) :: nborder
               ! local processor for subdomain to the north of this one
      integer, intent(in) :: nb
               ! flag to indicate whether a subdomain is on the south boundary
      integer, intent(in) :: sb
               ! local processor for subdomain to the north of this one
      integer, intent(in) :: lpn
               ! local processor for subdomain to the south of this one
      integer, intent(in) :: lps
               ! flag to indicate north-south periodicity or islands
      integer, intent(in) :: nsperiod
               ! flag to indicate south-north periodicity or islands
      integer, intent(in) :: snperiod
               ! message identifier
      integer, intent(in) :: msg_id
               ! message passing communicator
      integer, intent(in) :: msg_comm
!
! !INPUT/OUTPUT PARAMETERS
               ! a 3D array, whose borders are to be passed
      real*8 , intent(inOut) ::  array(ilo:ihi, jlo:jhi, klo:khi)
               ! working arrays to carry the messages
      real*8 , intent(inOut) :: pass_to_north  (ilo:ihi, klo:khi, nborder)
      real*8 , intent(inOut) :: pass_to_south  (ilo:ihi, klo:khi, nborder)
      real*8 , intent(inOut) :: pass_from_north(ilo:ihi, klo:khi, nborder)
      real*8 , intent(inOut) :: pass_from_south(ilo:ihi, klo:khi, nborder)
!
! !DESCRIPTION:
!   This routine is used to copy, transmit, and receive data from neighboring
!   subdomain processors into the ghost zones on the north and south sides.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      integer :: i, k
      integer :: ierrwin
      integer :: indwin
      integer :: iouter
      integer :: irc
      integer :: lperror
      integer :: northerly, southerly
      integer :: ns_x, sn_x
      integer :: sow_sn
      integer, save :: lpn2, lps2
      integer, save :: lpn_save, lps_save
      integer, save :: nsperiod_save, snperiod_save
      integer, save :: num_send
      integer, save :: send_grp
      integer, save :: rank_send(2)
!
!EOP
!------------------------------------------------------------------------------
!BOC


!     =========================================
      IOUTLOOP: do iouter = 1, msg3d_outer_flag
!     =========================================

!       -------------------------------------------------------------------
!       Generate unique message numbers for this variable.
!       NOTE:  message numbers will be unique if msg_id is chosen to be
!       incremented by at least 4 for all variables which use this routine.
!       -------------------------------------------------------------------

        ns_x = msg_id + NS_ESM
        sn_x = msg_id + SN_ESM

!       ----------------------------------------------------------
!       Gather north-south inner borders (3D equations).
!       These loops pack up 3D data from the outermost real zones
!       of this processor's subdomain into a message to be sent to
!       the neighboring subdomain.
!       ----------------------------------------------------------

        northerly = 0
        southerly = 0

        if ((nb == 0) .or. (nsperiod /= 0)) then

          northerly = 1

          do irc = 1, nborder
            do k = klo, khi
              do i = ilo, ihi
                pass_to_north(i,k,irc) = array(i,jhi-nborder+1-irc,k)
              end do
            end do
          end do

        end if

        if ((sb == 0) .or. (snperiod /= 0)) then

          southerly = 1

          do irc = 1, nborder
            do k = klo, khi
              do i = ilo, ihi
                pass_to_south(i,k,irc) = array(i,jlo+nborder-1+irc,k)
              end do
            end do
          end do

        end if

        sow_sn = (ihi - ilo + 1) * (khi - klo + 1) * nborder


!         -----------------------------------------------------
!         Transmit from south to north for multiprocessor runs.
!         -----------------------------------------------------

          if ((nb /= 0) .and. (nsperiod == 0)) then

            call receiveReal8 (lps, sow_sn, sn_x, pass_from_south, msg_comm)

          else if ((sb /= 0) .and. (snperiod == 0)) then

            call sendReal8 (lpn, sow_sn, sn_x, pass_to_north,   msg_comm)

          else

            call sendReceiveReal8  &
     &        (lps, sow_sn, sn_x, pass_to_north,  &
     &         lpn, sow_sn, sn_x, pass_from_south, msg_comm)

          end if

!         -----------------------------------------------------
!         Transmit from north to south for multiprocessor runs.
!         -----------------------------------------------------

          if ((sb /= 0) .and. (snperiod == 0)) then

            call receiveReal8 (lpn, sow_sn, ns_x, pass_from_north, msg_comm)

          else if ((nb /= 0) .and. (nsperiod == 0)) then

            call sendReal8 (lps, sow_sn, ns_x, pass_to_south,   msg_comm)

          else

            call sendReceiveReal8  &
     &        (lpn, sow_sn, ns_x, pass_to_south,  &
     &         lps, sow_sn, ns_x, pass_from_north, msg_comm)

          end if


!       -------------------------------------------------------------
!       Scatter south-north outer borders (3D equations).
!       These loops unpack the message received into the ghost zones.
!       -------------------------------------------------------------

        if ((sb == 0) .or. (snperiod /= 0)) then

            do irc = 1, nborder
              do k = klo, khi
                do i = ilo, ihi
                  array(i,nborder+jlo-irc,k) = pass_from_south(i,k,irc)
                end do
              end do
            end do
        end if

        if ((nb == 0) .or. (nsperiod /= 0)) then

            do irc = 1, nborder
              do k = klo, khi
                do i = ilo, ihi
                  array(i,jhi-nborder+irc,k) = pass_from_north(i,k,irc)
                end do
              end do
            end do
        end if

!     ===============
      end do IOUTLOOP
!     ===============

      return

      end subroutine Msg3d_Ns
!EOC
!------------------------------------------------------------------------

      end module GmiGhostZones_mod
