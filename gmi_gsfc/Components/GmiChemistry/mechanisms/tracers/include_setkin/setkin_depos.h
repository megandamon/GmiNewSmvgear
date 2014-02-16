!=======================================================================
!
! FILE
!   setkin_depos.h
!
! DESCRIPTION
!   This include file contains physical information about each species which
!   will be used in the deposition routines. It includes an integer declaring
!   aerosol species (aerosol=1 means it is an aerosol, aerosol=0 means it is
!   a gas). There is also a Henry's law constant and oxidizing potential for
!   each species.
!
!
!=======================================================================

      integer      :: aerosol         (1:NSP)

      real*8, save :: hstar           (1:NSP)
      real*8, save :: oxidize         (1:NSP)
      real*8, save :: delH_298_over_R (1:NSP)
      real*8, save :: retention_eff   (1:NSP)

      data aerosol (1:NSP) / &
     &   0,   0,   0,   0,  38,  &
     &  38,  39,  39,  39,  39,  &
     &   0,   0,   0,   0,   0,  &
     &   0,   0,   0 /

      data hstar (1:NSP) / &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   1.000D-02 /

      data delH_298_over_R (1:NSP) / &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00 /

      data oxidize (1:NSP) / &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   1.000D+00 /

      data retention_eff (1:NSP) / &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,   0.000D+00,  &
     &   0.000D+00,   0.000D+00,   0.000D+00 /

!                                  --^--
