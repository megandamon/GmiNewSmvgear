
#include "GmiESMF_ErrLog.h"
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: RUTMod - Implements Interface to a minimalistic MAPL GC that
!                   serves as a parent of the MAPL_ExtData GC
!
! !INTERFACE:
!
   MODULE GmiAdvecCoreStubs_mod
!
! !USES:
!
   USE ESMF_Mod
   USE GmiESMF_ErrorChecking_mod

   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices

!
! !DESCRIPTION: 
!
!  {\tt RUT} Root Utility Test is a minimalistic parent grid component that creates 
!  ExtData grid component.
!
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov  Implemented as part of the ExtData 
!                               utility test.
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for the RUT
!
! !INTERFACE:

   SUBROUTINE SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------

      character(len=ESMF_MAXSTR) :: IAm = "SetServices"
      integer                         :: STATUS

!   Local derived type aliases
    type(ESMF_Config)          :: CF
    character(len=ESMF_MAXSTR) :: comp_name

!                              ------------

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, name=comp_name, rc = STATUS )
    Iam = trim(comp_name) // '::' // trim(Iam)

!   Greetings
!   ---------
!    if (MAPL_am_I_root()) then
!         print *, trim(Iam)//': ACTIVE'
!         print *
!    end if

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
!    call MAPL_GridCompSetEntryPoint ( GC, ESMF_SETINIT,  Initialize_, rc = STATUS )
!    call MAPL_GridCompSetEntryPoint ( GC, ESMF_SETRUN,   Run_,        rc = STATUS )
!    call MAPL_GridCompSetEntryPoint ( GC, ESMF_SETFINAL, Finalize_,   rc = STATUS )

!   Generic Set Services
!   --------------------
!    call MAPL_GenericSetServices ( GC, rc = STATUS )

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize RUT
!
! !INTERFACE:
!

   SUBROUTINE Initialize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out)            :: rc         ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------

      character(len=ESMF_MAXSTR) :: IAm = "Initialize_"
      integer :: STATUS

   type(ESMF_Grid)             :: GRID        ! Grid
   type(ESMF_Config)           :: CF          ! Universal Config 

   character(len=ESMF_MAXSTR)     :: comp_name

  
!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, config=CF, rc = STATUS )
   Iam = trim(comp_name) // '::' // trim(Iam)

!  Create grid for this GC
!  ------------------------
!   call MAPL_GridCreate  (GC, rc = STATUS )
!   call ESMF_GridCompGet (GC, grid=GRID, rc = STATUS)

!  Initialize MAPL Generic
!  -----------------------
!   call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, clock,  rc = STATUS )

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs RUT
!
! !INTERFACE:
!

   SUBROUTINE Run_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------
      character(len=ESMF_MAXSTR) :: IAm = "Run_"
      integer :: STATUS

   character(len=ESMF_MAXSTR)    :: comp_name


!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, rc = STATUS )
   Iam = trim(comp_name) // '::' // trim(Iam)

!   Call Run for every Child
!   -------------------------
!    call MAPL_GenericRun ( GC, IMPORT, EXPORT, CLOCK,  rc = STATUS)


!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize RUT
!
! !INTERFACE:
!

   SUBROUTINE Finalize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------

      character(len=ESMF_MAXSTR) :: IAm = "Finalize_"
      integer :: STATUS


   character(len=ESMF_MAXSTR)  :: comp_name

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, rc = STATUS )
   Iam = trim(comp_name) // trim(Iam)

!  Finalize MAPL Generic
!  ---------------------
!   call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  rc = STATUS )

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

 end SUBROUTINE Finalize_
!EOC
!--------------------------------------------------------------------
end module GmiAdvecCoreStubs_mod
