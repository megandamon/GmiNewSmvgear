!------------------------------------------------------------------------------
! NASA/GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
! !MODULE: GmiMessagePassing_mod
!
! !INTERFACE:
!
    module GmiMessagePassing_mod
!
    implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
    private
    public  :: sendCharacter , receiveCharacter, broadcastCharacter
    public  :: sendInteger   , receiveInteger  , broadcastInteger
    public  :: sendLogical   , receiveLogical  , broadcastLogical
    public  :: sendReal8     , receiveReal8    , broadcastReal8
    public  :: stopCode      , writeMpiError   , sendReceiveReal8
    public  :: writeDebugInfo, synchronizeGroup
!
# include "mpif.h"
!
! !DESCRIPTION:
! This module contains routines performing basic message passing operations
! (send, receive, broadcast, send/receive) and error checking.
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
! ! REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
    contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sendCharcater
!
! !INTERFACE:
!
      subroutine sendCharacter (dest, count, tag, message, commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: dest
      integer          , intent(in) :: count
      integer          , intent(in) :: tag
      integer          , intent(in) :: commu
      character (len=1), intent(in) :: message(count)
!
! !DESCRIPTION:
! Performs a send operation on a charcater array
! and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr
      integer :: rank

!     =============
      call MPI_Send (message, count, MPI_CHARACTER, dest, tag, commu, ierr)
!     =============

      if (ierr /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr)

        call MPI_Comm_Rank (commu, rank, ierr)

        call writeDebugInfo  &
     &    ('send', rank, 'dest', dest, 'tag', tag,  &
     &     'char', 'count', count, 'block', commu)

      end if

      return

      end subroutine sendCharacter
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: receiveCharacter
!
! !INTERFACE:
!
      subroutine receiveCharacter (source, count, tag, message, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: source
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
!
! !OUTPUT PARAMETERS:
      character (len=1), intent(out) :: message(count)
!
! !DESCRIPTION:
! Performs a receive operation on a character array
! and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      logical :: fail
      integer :: count_chk
      integer :: ierr1
      integer :: rank
      integer :: status(MPI_STATUS_SIZE)

      fail = .false.

!     =============
      call MPI_Recv(message, count, MPI_CHARACTER, source, tag, commu, status, ierr1)
!     =============

      if (ierr1 /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr1)

        fail = .true.

      else

!       ==================
        call MPI_Get_Count (status, MPI_CHARACTER, count_chk, ierr1)
!       ==================

        if ((ierr1 /= MPI_SUCCESS) .or. (count /= count_chk)) then

          fail = .true.

          if (ierr1 /= MPI_SUCCESS) then
            call writeMpiError (commu, .false., ierr1)
          end if

        end if

      end if

      if (fail) then

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr1)
!       ==================

        call writeDebugInfo ('recv', rank, 'source', source, 'tag', tag,  &
     &     'char', 'count', count, 'block', commu)

      end if

      return

      end subroutine receiveCharacter
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: broadcastCharacter 
!
! !INTERFACE:
!
      subroutine broadcastCharacter (buffer, count, root, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: count
      integer, intent(in) :: root
      integer, intent(in) :: commu
!
! !INPUT/OUTPUT PARAMETERS:
      character (len=1), intent(inOut) :: buffer(count)
!
! !DESCRIPTION:
! Broadcasts a charcater array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr

!     ==============
      call MPI_Bcast (buffer, count, MPI_CHARACTER, root, commu, ierr)
!     ==============

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if


      return

      end subroutine broadcastCharacter
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sendInteger
!
! !INTERFACE:
!
      subroutine sendInteger (dest, count, tag, message, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: dest
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
      integer, intent(in) :: message(count)
!
! !DESCRIPTION:
! Performs a send operation on an integer array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr
      integer :: rank

!     =============
      call MPI_Send (message, count, MPI_INTEGER, dest, tag, commu, ierr)
!     =============

      if (ierr /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr)

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr)
!       ==================

        call writeDebugInfo ('send', rank, 'dest', dest, 'tag', tag,  &
     &     'integer', 'count', count, 'block', commu)

      end if

      return

      end subroutine sendInteger 
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: receiveInteger
!
! !INTERFACE:
!
      subroutine receiveInteger (source, count, tag, message, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: source
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: message(count)
!
! !DESCRIPTION:
! Performs a receive operation on an integer array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      logical :: fail
      integer :: count_chk
      integer :: ierr1
      integer :: rank
      integer :: status(MPI_STATUS_SIZE)

      fail = .false.

!     =============
      call MPI_Recv (message, count, MPI_INTEGER, source, tag, commu, status, ierr1)
!     =============

      if (ierr1 /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr1)

        fail = .true.

      else

!       ==================
        call MPI_Get_Count (status, MPI_INTEGER, count_chk, ierr1)
!       ==================

        if ((ierr1 /= MPI_SUCCESS) .or. (count /= count_chk)) then

          fail = .true.

          if (ierr1 /= MPI_SUCCESS) then
            call writeMpiError (commu, .false., ierr1)
          end if

        end if

      end if


      if (fail) then

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr1)
!       ==================

        call writeDebugInfo ('recv', rank, 'source', source, 'tag', tag,  &
     &     'integer', 'count', count, 'block', commu)

      end if


      return

      end subroutine receiveInteger
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: broadcastInteger
!
! !INTERFACE:
!
      subroutine broadcastInteger (buffer, count, root, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: count
      integer, intent(in) :: root
      integer, intent(in) :: commu
!
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(inOut) :: buffer(count)
!
! !DESCRIPTION:
! Broadcasts an integer array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr

!     ==============
      call MPI_Bcast (buffer, count, MPI_INTEGER, root, commu, ierr)
!     ==============

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if


      return

      end subroutine broadcastInteger
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sendLogical
!
! !INTERFACE:
!
      subroutine sendLogical (dest, count, tag, message, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: dest
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
      logical, intent(in) :: message(count)
!
! !DESCRIPTION:
! Performs a send operation on a logical array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr
      integer :: rank

!     =============
      call MPI_Send (message, count, MPI_LOGICAL, dest, tag, commu, ierr)
!     =============


      if (ierr /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr)

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr)
!       ==================

        call writeDebugInfo ('send', rank, 'dest', dest, 'tag', tag,  &
     &     'char', 'count', count, 'block', commu)

      end if

      return

      end subroutine sendLogical
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: receiveLogical
!
! !INTERFACE:
!
      subroutine receiveLogical (source, count, tag, message, commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: source
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
!
! !OUTPUT PARAMETERS:
      logical, intent(out) :: message(count)
!
! !DESCRIPTION:
! Performs a receive operation on a logical array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      logical :: fail
      integer :: count_chk
      integer :: ierr1
      integer :: rank
      integer :: status(MPI_STATUS_SIZE)

      fail = .false.


!     =============
      call MPI_Recv (message, count, MPI_LOGICAL, source, tag, commu, status, ierr1)
!     =============

      if (ierr1 /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr1)

        fail = .true.

      else

!       ==================
        call MPI_Get_Count (status, MPI_LOGICAL, count_chk, ierr1)
!       ==================

        if ((ierr1 /= MPI_SUCCESS) .or. (count /= count_chk)) then

          fail = .true.

          if (ierr1 /= MPI_SUCCESS) then
            call writeMpiError (commu, .false., ierr1)
          end if

        end if

      end if

      if (fail) then

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr1)
!       ==================

        call writeDebugInfo ('recv', rank, 'source', source, 'tag', tag,  &
     &     'char', 'count', count, 'block', commu)

      end if

      return

      end subroutine receiveLogical
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: broadcastLogical
!
! !INTERFACE:
!
      subroutine broadcastLogical (buffer, count, root, commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: count
      integer, intent(in) :: root
      integer, intent(in) :: commu
!
! !INPUT/OUTPUT PARAMETERS:
      logical, intent(inOut) :: buffer(count)
!
! !DESCRIPTION:
! Broadcasts a logical array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr

!     ==============
      call MPI_Bcast (buffer, count, MPI_LOGICAL, root, commu, ierr)
!     ==============

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if

      return

      end subroutine broadcastLogical
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sendReal8
!
! !INTERFACE:
!
      subroutine sendReal8 (dest, count, tag, message, commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: dest
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
      real*8 , intent(in) :: message(count)
!
! !DESCRIPTION:
! Performs a send operation on a double precision array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr
      integer :: rank

!     =============
      call MPI_Send (message, count, MPI_REAL8, dest, tag, commu, ierr)
!     =============

      if (ierr /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr)

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr)
!       ==================

        call writeDebugInfo ('send', rank, 'dest', dest, 'tag', tag,  &
     &     'real8', 'count', count, 'block', commu)

      end if

      return

      end subroutine sendReal8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: receiveReal8
!
! !INTERFACE:
!
      subroutine receiveReal8 (source, count, tag, message, commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: source
      integer, intent(in) :: count
      integer, intent(in) :: tag
      integer, intent(in) :: commu
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: message(count)
!
! !DESCRIPTION:
! Performs a receive operation on a double precision array 
! and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      logical :: fail
      integer :: count_chk
      integer :: ierr1
      integer :: rank
      integer :: status(MPI_STATUS_SIZE)

      fail = .false.

!     =============
      call MPI_Recv (message, count, MPI_REAL8, source, tag, commu, status, ierr1)
!     =============

      if (ierr1 /= MPI_SUCCESS) then

        call writeMpiError (commu, .false., ierr1)

        fail = .true.

      else

!       ==================
        call MPI_Get_Count (status, MPI_REAL8, count_chk, ierr1)
!       ==================

        if ((ierr1 /= MPI_SUCCESS) .or. (count /= count_chk)) then

          fail = .true.

          if (ierr1 /= MPI_SUCCESS) then
            call writeMpiError (commu, .false., ierr1)
          end if

        end if

      end if

      if (fail) then

!       ==================
        call MPI_Comm_Rank (commu, rank, ierr1)
!       ==================

        call writeDebugInfo ('recv', rank, 'source', source, 'tag', tag,  &
     &     'real8', 'count', count, 'block', commu)

      end if

      return

      end subroutine receiveReal8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sendReceiveReal8
!
! !INTERFACE:
!
      subroutine sendReceiveReal8 (source, count1, tag1, message1,  &
     &               dest, count2, tag2, message2, commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: source, dest
      integer, intent(in) :: count1, count2
      integer, intent(in) :: tag1, tag2
      integer, intent(in) :: commu
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: message1(count1), message2(count2)
!
! !DESCRIPTION:
! Performs a send/receive operation on a double precision array
! and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      logical :: fail
      integer :: count_chk
      integer :: ierr
      integer :: status(MPI_STATUS_SIZE)

      fail = .false.

!     =================
      call MPI_Sendrecv  &
!     =================
     &  (message1, count1, MPI_REAL8, dest,   tag1,  &
     &   message2, count2, MPI_REAL8, source, tag2,  &
     &   commu, status, ierr)

      if (ierr /= MPI_SUCCESS) then
         call writeMpiError (commu, .false., ierr)
         fail = .true.
      else

!       ==================
        call MPI_Get_Count (status, MPI_REAL8, count_chk, ierr)
!       ==================

        if (ierr /= MPI_SUCCESS) then
           fail = .true.
           call writeMpiError (commu, .false., ierr)
        end if

      end if

      if (fail) then
!       =============
        call stopCode (commu, " SendReceiveReal8 ERROR ")
!       =============
      end if


      return

      end subroutine sendReceiveReal8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: broadcastReal8
!
! !INTERFACE:
!
      subroutine broadcastReal8 (buffer, count, root, commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: count
      integer, intent(in) :: root
      integer, intent(in) :: commu
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: buffer(count)
!
! DESCRIPTION:
! Broadcasts a double precision array and does error checking.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr

!     ==============
      call MPI_Bcast (buffer, count, MPI_REAL8, root, commu, ierr)
!     ==============

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if

      return

      end subroutine broadcastReal8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: synchronizeGroup
!
! !INTERFACE:
!
      subroutine synchronizeGroup (commu)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: commu
!
! !DESCRIPTION:
!   This routine imposes a barrier.  That is, no process in" group" exits
!   until all processes in "group" have at least begun to execute the
!   routine.  No guarantees are made about how close to the same time the
!   processes complete this barrier.

!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ierr

      call MPI_Barrier (commu, ierr)

      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (commu, .true., ierr)
      end if

      return

      end subroutine synchronizeGroup
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeDebugInfo
!
! !INTERFACE:
!
      subroutine writeDebugInfo  &
     &  (op_type, mynode, other_node_str, other_node, tag_str, tag,  &
     &   data_type, buffer_length_str, buffer_length, block_type, &
     &   commu)

      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: commu
      character (len=*), intent(in) :: op_type
      integer          , intent(in) :: mynode
      character (len=*), intent(in) :: other_node_str
      integer          , intent(in) :: other_node
      character (len=*), intent(in) :: tag_str
      integer          , intent(in) :: tag
      character (len=*), intent(in) :: data_type
      character (len=*), intent(in) :: buffer_length_str
      integer          , intent(in) :: buffer_length
      character (len=*), intent(in) :: block_type
!
! !DESCRIPTION:
! Writes debugging information.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      Write (6,900)  &
     &  op_type, mynode,  &
     &  other_node_str, other_node,  &
     &  tag_str, tag, data_type,  &
     &  buffer_length_str, buffer_length,  &
     &  block_type

 900  format (/, ' Error in ', a8, ' macro', /,  &
     &  5x, ' my_node:       ', 32x, ' = ', i16, /,  &
     &  5x, ' other_node:    ', a32, ' = ', i16, /,  &
     &  5x, ' tag:           ', a32, ' = ', i16, /,  &
     &  5x, ' data type:     ', a32, /,  &
     &  5x, ' buffer_length: ', a32, ' = ', i16, /,  &
     &  5x, ' block_type:    ', a32)

!     =============
      call stopCode (commu, " writeDebugInfo ERROR ")
!     =============

      return

      end subroutine writeDebugInfo
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeMpiError
!
! !INTERFACE:
!
      subroutine writeMpiError (commu, do_stop, ierr1)

      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: do_stop
      integer, intent(in) :: commu, ierr1
!
! !DESCRIPTION:
! Determines the error message associated with the error code "ierr1" and
! stop all the processes if "do\_stop" is set to true.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      character (len=MPI_MAX_ERROR_STRING) :: errstr
      integer :: ierr2
      integer :: ierrlen
      integer :: ii

!     =====================
      call MPI_Error_String (ierr1, errstr, ierrlen, ierr2)
!     =====================

      Print '(256a)', (errstr(ii:ii), ii=1,ierrlen)

      if (do_stop) then
!       =============
        call stopCode (commu, " writeMpiError ERROR ")
!       =============
      end if


      return

      end subroutine writeMpiError
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stopCode
!
! !INTERFACE:
!
      subroutine stopCode (commu, msg)

      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: commu
      character (len=*), intent(in) :: msg
!
! !DESCRIPTION:
! Stops all the processes associated with the communicator.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: icode, ierr

      Write (6,*) msg

      call MPI_Abort (commu, icode, ierr)

      Stop

      return

      end subroutine stopCode
!EOC
!------------------------------------------------------------------------------

      end module GmiMessagePassing_mod
