!=======================================================================
!
! $Id: setkin_synspc.h,v 1.5 2011-08-08 17:46:52 mrdamon Exp $
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
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!
!=======================================================================


      logical, parameter :: USE_SYNOZ = .false.
      logical, parameter :: USE_NODOZ = .false.

      integer, parameter :: MAXNODOZ_ELEM = 10

      integer, parameter :: NOX_ELEM_MAP(MAXNODOZ_ELEM) =  &
     &  (/ 10, 12, 13, 14, 15,  0,  0,  0,  0,  0 /)

      integer, parameter :: NOY_ELEM_MAP(MAXNODOZ_ELEM) =  &
     &  (/ 16, 17, 18,  0,  0,  0,  0,  0,  0,  0 /)

!                                  --^--

