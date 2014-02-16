!=======================================================================
!
! $Id: setkin_synspc.h,v 1.8 2011-08-09 22:52:11 mrdamon Exp $
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
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_SO2_JPL06.db
!  Setkin files generated:  Thu Nov 19 21:21:30 2009
!
!=======================================================================


      logical, parameter :: USE_SYNOZ = .true.
      logical, parameter :: USE_NODOZ = .false.

      integer, parameter :: MAXNODOZ_ELEM = 10

      integer, parameter :: NOX_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/ INO, INO2, INO3, IN2O5, IN2O5,  0,  0,  0,  0,  0 /)

      integer, parameter :: NOY_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/  IHNO2,  IHNO3,  IHNO4,  0,  0,  0,  0,  0,  0,  0 /)

!                                  --^--

