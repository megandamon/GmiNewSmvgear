!=======================================================================
!
! $Id: setkin_group.h,v 1.1 2006-07-03 02:54:25 kouatch Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_group.h
!
! DESCRIPTION
!   This include file contains information about grouping species for
!   transport purposes.
!
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!
!=======================================================================

      integer, parameter :: NUMGRP      =  0

      integer, parameter :: MAXGRP_ELEM = 10

      character*16, save :: cgrp_nam(1)

      data cgrp_nam(1) /"YOy"/

      integer :: mfe1, nf1

      integer, save :: sgrp_elem_map(MAXGRP_ELEM, 1)

      data  ((sgrp_elem_map(mfe1,nf1), mfe1=1,MAXGRP_ELEM),  &
     &                                 nf1 =1,1) /  &
     &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

      real*8,  save :: sgrp_fac(MAXGRP_ELEM, 1)

      data  ((sgrp_fac(mfe1,nf1), mfe1=1,MAXGRP_ELEM),  &
     &                            nf1 =1,1) /  &
     & 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,  &
     & 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /

!      -----------------------------------------------------------------
!      qq2 : grouped species sum
!      -----------------------------------------------------------------

!cc   integer :: ig ,im ,imsgrp

      real*8,  save :: max_bry_adjust
      real*8,  save :: max_cly_adjust
      real*8,  save :: max_nox_adjust

      real*8  :: qqgrp        (i1:i2 ,ju1:j2 ,k1:k2 ,1)

      real*8  :: group_adjust (i1:i2 ,ju1:j2 ,k1:k2)
      real*8  :: group_factor (i1:i2 ,ju1:j2 ,k1:k2)
      real*8  :: qq2          (i1:i2 ,ju1:j2 ,k1:k2)

!                                  --^--
