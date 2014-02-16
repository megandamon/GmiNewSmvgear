!=======================================================================
!
! $Id: setkin_synspc.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_synspc.h
!
! DESCRIPTION
!   This include file contains information about synthetic stratospheric
!   sources for ozone and nitrogen oxides.
!
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
!
!=======================================================================


      logical, parameter :: USE_SYNOZ = .false.
      logical, parameter :: USE_NODOZ = .false.

      integer, parameter :: ISYNOZ = 0

      integer, parameter :: MAXNODOZ_ELEM = 10

      integer, parameter :: NOX_ELEM_MAP(MAXNODOZ_ELEM) =  &
     &  (/ 51, 51, 52, 53, 54,  0,  0,  0,  0,  0 /)

      integer, parameter :: NOY_ELEM_MAP(MAXNODOZ_ELEM) =  &
     &  (/ 23, 24, 25,  0,  0,  0,  0,  0,  0,  0 /)

!                                  --^--

