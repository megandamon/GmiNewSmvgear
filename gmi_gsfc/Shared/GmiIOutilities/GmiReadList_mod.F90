!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiReadList_mod
!
! !INTERFACE:
!
   module GmiReadList_mod
!
! !USES:
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   private
   public  Read_List
!
! !DESCRIPTION:
!  Provides a routine to read items from an ASCII file.
!
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Read_List
!
! !INTERFACE:
!
      subroutine Read_List  &
        (list_items, num_list_items, list_to_read, max_list_items)
!
! !USES:
!
      use GmiPrintError_mod,      only : GmiPrintError
      use GmiASCIIoperations_mod, only : AsciiOpenRead

      implicit none
!
! !INPUT PARAMETERS:
!!    max_list_items : maximum number of list items allowed
!!    list_to_read   : file to read the list items from
      integer          , intent(in)   :: max_list_items
      character (len=*), intent(in)   :: list_to_read
!
! !OUTPUT PARAMETERS:
!!    list_items     : array to put the list items into
!!    num_list_items : number of items in the list_items array
      character (len=*), intent(out)  :: list_items(max_list_items)
      integer          , intent(out)  :: num_list_items
!
! !DESCRIPTION:
!  Reads items out of a file and puts them into an array.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      character (len=80) :: tmp_string
      integer            :: il
      integer            :: iost
      integer            :: rllun
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      err_msg(1:2)   = '  '
      num_list_items = 0

      call AsciiOpenRead (rllun, list_to_read)

      READLOOP: do il = 1, max_list_items

        Read (rllun, *, iostat = iost) list_items(il)

!                      =============
        if (iost /= 0) exit READLOOP
!                      =============

        num_list_items = il

      end do READLOOP

      if (iost > 0) then
         err_msg = 'In Read_List:  Error reading file names #1.'
      else if (num_list_items == 0) then
         err_msg = 'In Read_List:  No file names.'
      else if (num_list_items == max_list_items) then
         Read (rllun, *, iostat = iost) tmp_string

         if (iost > 0) then
            err_msg = 'In Read_List:  Error reading file names #2.'
         else if (iost == 0) then
            err_msg = 'In Read_List:  Too many file names.'
         end if
      end if
 
      if (err_msg(1:2) == 'In') then
         call GmiPrintError (err_msg, .true., 2, num_list_items,   &
                         max_list_items, 0, 0.0d0, 0.0d0)
      end if

      Close (rllun)
 
      return
 
      end subroutine Read_List
!EOC
!------------------------------------------------------------------
end module GmiReadList_mod
