!=======================================================================
!
! $Id: setkin_depos.h,v 1.2 2011-08-08 17:46:52 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
!   "setkin_depos.h" was modified by Xiaohong Liu
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
!  Chemistry input file generated: 2/23/00 10:12AM
!  Reaction dictionary:            3/29/99 2:55PM
!  Setkin files generated:         Wed Feb 23 10:30:41 2000
!
!=======================================================================

      integer :: aerosol(1:NSP)

      real*8, save :: hstar          (1:NSP)
      real*8, save :: oxidize        (1:NSP)
      real*8, save :: delH_298_over_R(1:NSP)
      real*8, save :: retention_eff  (1:NSP)

      data aerosol(1:NSP) /  &
     &  0,  0,  0,  0, 21, 20, 21, 22, 21, 20,  &
     & 21, 22,  2,  4,  4,  2,  2, 27, 28, 29,  &
     & 30, 23, 24, 25, 26,  0,  0,  0,  0,  0 /

      data hstar(1:NSP) /  &
     & 1.00d5,    1.00d5,    1.00d5,    0.56d0,    0.00d0,  &
     & 0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     & 0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     & 0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     & 0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     & 1.00d-2,   0.00d0,    4.00d3,    2.00d0,    0.00d0 /

      data oxidize(1:NSP) /  &
     &    1.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     &    0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     &    0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     &    0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     &    0.00d0,    0.00d0,    0.00d0,    0.00d0,    0.00d0,  &
     &    1.00d0,    0.00d0,    0.00d0,    1.00d0,    0.00d0 /

      data delH_298_over_R(1:NSP) /  &
     &-6.80d3, -3.02d3, -3.02d3,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 2.30d3,  0.00d0,  0.00d0,  0.00d0,  0.00d0 /

      data retention_eff(1:NSP) /  &
     & 0.05d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,  &
     & 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0 /

