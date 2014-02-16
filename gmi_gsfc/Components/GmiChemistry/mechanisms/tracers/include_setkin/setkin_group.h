!=======================================================================
!
! FILE
!   setkin_group.h
!
! DESCRIPTION
!   This include file contains information about grouping species for
!   transport purposes.
!
!
!=======================================================================

      integer, parameter :: NUMGRP      =  0

      integer, parameter :: MAXGRP_ELEM = 0

      integer, save :: sgrp_elem_map(MAXGRP_ELEM, NUMGRP)

      real*8,  save :: sgrp_fac(MAXGRP_ELEM, NUMGRP)

      real*8  :: qqgrp (i1:i2 ,ju1:j2 ,k1:k2 ,1)

!                                  --^--
