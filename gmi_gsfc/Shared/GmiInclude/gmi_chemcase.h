
!=============================================================================
!
!
! CODE DEVELOPER
!   Jules Kouatchou
!
! FILE
!   gmi_chemcase.h
!
! DESCRIPTION
!   This include file contains the variable used to store the name of the
!   chemical mechanism selected by the user.
!
!=============================================================================


      character (len=25) :: chem_mecha

!     =====================
      common  / gmichem_mecha /  &
!     =====================
     &  chem_mecha
