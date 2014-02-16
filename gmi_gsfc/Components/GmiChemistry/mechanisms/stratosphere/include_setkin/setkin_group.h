!=======================================================================
!
! $Id: setkin_group.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
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
!  Chemistry input file:    GMIS2 4:59 PM 7/30/2002
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Mon Jan 13 22:23:01 2003
!
!=======================================================================

      integer, parameter :: NUMGRP      =  2

      integer, parameter :: MAXGRP_ELEM = 10

      character*16, save :: cgrp_nam(NUMGRP)

      data cgrp_nam(1) /"Bry"/
      data cgrp_nam(2) /"Cly"/

      integer :: mge1, ng1

      integer, save :: sgrp_elem_map(MAXGRP_ELEM, NUMGRP)

      data  ((sgrp_elem_map(mge1,ng1), mge1=1,MAXGRP_ELEM),  &
     &                                 ng1 =1,NUMGRP) /  &
     & 31,32,33,34,35,36, 0, 0, 0, 0,  &
     & 23,24,25,26,27,28,29,30,32, 0 /

      real*8,  save :: sgrp_fac(MAXGRP_ELEM, NUMGRP)

      data  ((sgrp_fac(mge1,ng1), mge1=1,MAXGRP_ELEM),  &
     &                            ng1 =1,NUMGRP) /  &
     & 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,  &
     & 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,  &
     & 1.0D+00, 2.0D+00, 1.0D+00, 1.0D+00, 2.0D+00,  &
     & 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /

!     ------------------------------------------------------------------------
!     qq2 : grouped species sum
!     ------------------------------------------------------------------------

      real*8,  save :: max_bry_adjust
      real*8,  save :: max_cly_adjust
      real*8,  save :: max_nox_adjust

        real*8  :: qqgrp        (i1:i2, ju1:j2, k1:k2, NUMGRP)

      real*8  :: group_adjust (i1:i2, ju1:j2, k1:k2)
      real*8  :: group_factor (i1:i2, ju1:j2, k1:k2)
      real*8  :: qq2          (i1:i2, ju1:j2, k1:k2)

!                                  --^--
