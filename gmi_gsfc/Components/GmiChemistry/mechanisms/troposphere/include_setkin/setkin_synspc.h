!=======================================================================
!
! $Id: setkin_synspc.h,v 1.7 2011-08-09 22:12:58 mrdamon Exp $
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
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!=======================================================================


      logical, parameter :: USE_SYNOZ = .true.
      logical, parameter :: USE_NODOZ = .true.

      integer, parameter :: MAXNODOZ_ELEM = 10

      integer, parameter :: NOX_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/ INO, INO2, INO3, IN2O5, IN2O5,  0,  0,  0,  0,  0 /)

      integer, parameter :: NOY_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/  IHNO2,  IHNO3,  IHNO4,  0,  0,  0,  0,  0,  0,  0 /)

!                                  --^--

