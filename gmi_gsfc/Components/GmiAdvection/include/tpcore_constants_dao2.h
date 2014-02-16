
!=============================================================================
!
! $Id: tpcore_constants_dao2.h,v 1.5 2011-08-09 22:12:57 mrdamon Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   tpcore_constants_dao2.h
!
! DESCRIPTION
!   This include file contains some constants for the DAO2 advection model.
!
!=============================================================================


      logical, parameter ::  &
     &  CROSS   = .true.,    & ! calculate cross terms?
     &  FILL    = .true.   ! fill negatives?


      integer, parameter ::  &
     &  IORD    = 3,         & ! controls various options in E-W      advection
     &  JORD    = 3,         & ! controls various options in N-S      advection
!c   &  KORD    = 5,       ! controls various options in vertical advection
!c   &  KORD    = 3,       ! controls various options in vertical advection
     &  KORD    = 7        ! controls various options in vertical advection


      real*8,  parameter ::  &
     &  UMAX    = 220.0d0  ! estimate (upper limit) of the maximum u-wind
                           ! speed (m/s)

