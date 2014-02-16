!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiASCIIoperations_mod
!
! !INTERFACE:
!
   module GmiASCIIoperations_mod
!
! !USES:
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   private
   public  AsciiOpenRead
   public  AsciiOpenWrite
!
! !DESCRIPTION:
!  Operations for opening ASCII files for reading and writing.
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
! !IROUTINE: AsciiOpenRead
!
! !INTERFACE:
!
      subroutine AsciiOpenRead (asc_lun, asc_filname)
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
      use GmiFileUnit_mod  , only : GetFileUnitNumber

      implicit none
!
! !INPUT PARAMETERS:
!!    asc_lun     : logical unit for opened ASCII file
!!    asc_filname : name of ASCII file to open for reading
      integer          , intent(out) :: asc_lun
      character (len=*), intent(in) :: asc_filname
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION
!  Opens an ASCII file for reading and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: ierr
      integer :: in1
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
!
      ierr = 0
      call GetFileUnitNumber (asc_lun, ierr)

      if (ierr /= 0) then
         in1 = Len_Trim (asc_filname)

         err_msg = 'Open problem #1 in AsciiOpenRead:  ' // asc_filname(1:in1)

         call GmiPrintError (err_msg, .true., 1, ierr, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = 0
      Open (unit = asc_lun, file = asc_filname, status = 'old', iostat = ierr)

      if (ierr /= 0) then
         in1 = Len_Trim (asc_filname)

         err_msg = 'Open problem #2 in AsciiOpenRead:  ' // asc_filname(1:in1)

         call GmiPrintError (err_msg, .true., 1, ierr, 0, 0, 0.0d0, 0.0d0)
      end if

      return
 
      end subroutine AsciiOpenRead
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AsciiOpenWrite
!
! !INTERFACE:
!
      subroutine AsciiOpenWrite (asc_lun, asc_filname)
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
      use GmiFileUnit_mod  , only : GetFileUnitNumber

      implicit none
!
! !INPUT PARAMETERS:
!!    asc_lun     : logical unit for opened ASCII file
!!    asc_filname : name of ASCII file to open for writing
      integer          , intent(out)  :: asc_lun
      character (len=*), intent(in)  :: asc_filname
!
! !DESCRIPTION:
!   This routine opens an ASCII file for writing and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: ierr
      integer :: in1
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
!
      ierr = 0
      call GetFileUnitNumber (asc_lun, ierr)
 
      if (ierr /= 0) then
         in1 = Len_Trim (asc_filname)
 
         err_msg = 'Open problem #1 in AsciiOpenWrite:  ' // asc_filname(1:in1)
 
         call GmiPrintError (err_msg, .true., 1, ierr, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = 0
      Open (unit = asc_lun, file = asc_filname, status = 'unknown', iostat = ierr)
 
      if (ierr /= 0) then
         in1 = Len_Trim (asc_filname)
 
         err_msg = 'Open problem #2 in AsciiOpenWrite:  ' // asc_filname(1:in1)
 
         call GmiPrintError (err_msg, .true., 1, ierr, 0, 0, 0.0d0, 0.0d0)
      end if

      return
 
      end subroutine AsciiOpenWrite
!EOC
!----------------------------------------------------------------
end module GmiASCIIoperations_mod
