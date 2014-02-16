!=======================================================================
!
! FILE
!   setkin_mw.h
!
! DESCRIPTION
!  This include file contains gram molecular weights
!  for the species in the model mechanism.
!
!
!=======================================================================

      real*8, save  :: mw_data(1:NSP)

      data mw_data(1:NSP) / &
     &   1.0000D+00,   2.8960D+01,   2.8960D+01,   2.2200D+02,  &
     &   2.1000D+02,   2.1000D+02,   7.0000D+00,   1.0000D+01,  &
     &   7.0000D+00,   1.0000D+01,   1.4194D+02,   4.4010D+01,  &
     &   4.8000D+01,   4.8000D+01,   1.4605D+02,   2.8960D+01,  &
     &   2.8000D+01,   4.8000D+01  /

!                                  --^--
