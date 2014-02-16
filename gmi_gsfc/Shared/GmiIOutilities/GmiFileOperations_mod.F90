!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiFileOperations_mod
!
! !INTERFACE:
!
   module GmiFileOperations_mod
!
! !USES:
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   private
   public  doesFileExist
   public  finishFileName
   public  makeOutfileName
   public  makeOutfileNames

#  include "GmiParameters.h"
!
! !DESCRIPTION:
!  Provide routines to find out if a file exists and to manipulate
!  file name.
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
! !IROUTINE: doesFileExist
!
! !INTERFACE:
!
      logical function doesFileExist (fname)
!
! !USES:
!
!
      implicit none
!
! !INPUT PARAMETERS:
!
!!    fname : the name of the file to check
      character (len=*), intent(in) :: fname
!
! !DESCRIPTION:
!   Checks for the existence of a file.
!
! !LOCAL VARIABLES:
      logical :: does_exist
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
      Inquire (file = fname, exist = does_exist)

      doesFileExist = does_exist

      return
 
      end function doesFileExist
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finishFileName
!
! !INTERFACE:
!
      subroutine finishFileName (fname, dpath)
!
! !USES:
!
      use GmiPrintError_mod
!
      implicit none
!
! !INPUT PARAMETERS:
!!    dpath : the directory path
      character (len=*), intent(in)    :: dpath
!
! !INPUT/OUTPUT PARAMETERS:
!!    fname : the file name
      character (len=*), intent(inout) :: fname
!
! !DESCRIPTION:
!  This routine prepends the directory path to the file name
!  if the file name does not begin with a "/".
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG)  :: err_msg
      character (len=MAX_LENGTH_FILE_NAME) :: tmp_string

      integer :: in1, in2
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
      if (fname(1:1) /= '/') then

        in1 = Len_Trim (dpath)
        in2 = Len_Trim (fname)

        if (((in1 + in2) > MAX_LENGTH_FILE_NAME) .or. (Len (fname) < MAX_LENGTH_FILE_NAME)) then
          err_msg = 'Problem in finishFileName.'
          call GmiPrintError (err_msg, .true., 2, in1 + in2, Len (fname), &
                           0, 0.0d0, 0.0d0)
        end if

        if ((in1 /= 0) .and. (in2 /= 0)) then
          tmp_string = dpath(1:in1) // '/' // fname(1:in2)
          fname      = tmp_string
        end if

      end if
 
      return

      end subroutine finishFileName
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: makeOutfileName
!
! !INTERFACE:
!
      subroutine makeOutfileName (fname, suffix, gmi_problem_name)
!
! !USES:
!
!
      implicit none
!
! !INPUT PARAMETERS:
!!   suffix : the suffix to add after the file name is trimmed of blank spaces
!!   gmi_problem_name: problem name
      character (len=*),   intent(in)  :: suffix
      character (len=*), intent(in)  :: gmi_problem_name
!
! !INPUT/OUTPUT PARAMETERS:
!!   fname  : the file name
      character (len=*), intent(inout) :: fname
!
! !DESCRIPTION:
!  This routine uses the problem name and a suffix to create a name for an
!  output file.
!
! !LOCAL VARIABLES:
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
      in1   = Len_Trim (gmi_problem_name)
 
      fname = gmi_problem_name(1:in1) // suffix

      return

      end subroutine makeOutfileName
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: makeOutfileNames
!
! !INTERFACE:
!
      subroutine makeOutfileNames (fname, suffix1, suffix2, gmi_problem_name)
!
! !USES:
!
!
      implicit none
!
! !INPUT PARAMETERS:
!!   suffix1 : First  suffix to add after the file name is trimmed of blank spaces
!!   suffix2 : Second suffix to add after the file name is trimmed of blank spaces
!!   gmi_problem_name: problem name
      character (len=*),   intent(in)  :: suffix1, suffix2
      character (len=*), intent(in)  :: gmi_problem_name
!
! !INPUT/OUTPUT PARAMETERS:
      character (len=*), intent(inout) :: fname ! the file name
!
! !DESCRIPTION:
!  This routine uses the problem name and two suffixes to create a name for an
!  output file.
!
! !LOCAL VARIABLES:
      integer :: in1, in2
!
! !AUTHOR: 
!  Bigyani Das and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      in1   = Len_Trim (gmi_problem_name)
      in2   = len_Trim (suffix1)
      fname = gmi_problem_name(1:in1) // '.' // suffix1(1:in2) // suffix2

      return

      end subroutine makeOutfileNames
!EOC
!-------------------------------------------------------------------------
end module GmiFileOperations_mod
