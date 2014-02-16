!------------------------------------------------------------------------------
! NASA GSFC - SSSO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiESMF_ErrorChecking_mod
!
! !INTERFACE:
!
      module GmiESMF_ErrorChecking_mod
!
! !USES:
      USE ESMF_Mod
!
      IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
      PRIVATE
      PUBLIC  :: GMI_VRFY
      PUBLIC  :: GMI_ASRT
      PUBLIC  :: GMI_RTRN
      PUBLIC  :: GmiESMF_Abort
!
! !DESCRIPTION:
! Routines/functions for error checking.
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GMI_VRFY
!
! !INTERFACE:
!
      logical function GMI_VRFY(A,iam,line,rc)
         integer,           intent(IN ) :: A
         character*(*),     intent(IN ) :: iam
         integer,           intent(IN ) :: line
         integer, optional, intent(OUT) :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
           GMI_VRFY = A/=ESMF_SUCCESS
           if(GMI_VRFY)then
             if(present(RC)) then
               print'(A40,I10)',Iam,line
               RC=A
             endif
           endif
      end function GMI_VRFY
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GMI_ASRT
!
! !INTERFACE:
!
      logical function GMI_ASRT(A,iam,line,rc)
         logical,           intent(IN ) :: A
         character*(*),     intent(IN ) :: iam
         integer,           intent(IN ) :: line
         integer, optional, intent(OUT) :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
           GMI_ASRT = .not.A
           if(GMI_ASRT)then
             if(present(RC))then
               print'(A40,I10)',Iam,LINE
               RC=ESMF_FAILURE
             endif
           endif
      end function GMI_ASRT
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GMI_RTRN
!
! !INTERFACE:
!
      logical function GMI_RTRN(A,iam,line,rc)
         integer,           intent(IN ) :: A
         character*(*),     intent(IN ) :: iam
         integer,           intent(IN ) :: line
         integer, optional, intent(OUT) :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
           GMI_RTRN = .true.
           if(A/=ESMF_SUCCESS)print'(A40,I10)',Iam,line
           if(present(RC)) RC=A
      end function GMI_RTRN
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GmiESMF_Abort
!
! !INTERFACE:
!
    subroutine GmiESMF_Abort
!EOP
!------------------------------------------------------------------------------
!BOC
    call ESMF_Finalize(terminationFlag = ESMF_Abort)
    RETURN
    end subroutine GmiESMF_Abort
!EOC
!-------------------------------------------------------------------------


      end module GmiESMF_ErrorChecking_mod
