!=======================================================================
!
! $Id: setkin_group.h,v 1.5 2011-08-09 22:12:58 mrdamon Exp $
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
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!=======================================================================

      integer, parameter :: NUMGRP      =  0

      integer, parameter :: MAXGRP_ELEM = 10

      character*16, save :: cgrp_nam(1)

      data cgrp_nam(1) /"YOy"/

      integer :: mge1, ng1

      integer, save :: sgrp_elem_map(MAXGRP_ELEM, 1)

      data  ((sgrp_elem_map(mge1,ng1), mge1=1,MAXGRP_ELEM), &
     &                                 ng1 =1,1) / &
     &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

      real*8,  save :: sgrp_fac(MAXGRP_ELEM, 1)

      data  ((sgrp_fac(mge1,ng1), mge1=1,MAXGRP_ELEM), &
     &                            ng1 =1,1) / &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00 /

!      -----------------------------------------------------------------
!      qq2 : grouped species sum
!      -----------------------------------------------------------------

      real*8,  save :: max_bry_adjust
      real*8,  save :: max_cly_adjust
      real*8,  save :: max_nox_adjust

      real*8  :: qqgrp        (i1:i2 ,ju1:j2 ,k1:k2 ,1)

      real*8  :: group_adjust (i1:i2 ,ju1:j2 ,k1:k2)
      real*8  :: group_factor (i1:i2 ,ju1:j2 ,k1:k2)
      real*8  :: qq2          (i1:i2 ,ju1:j2 ,k1:k2)

!                                  --^--
