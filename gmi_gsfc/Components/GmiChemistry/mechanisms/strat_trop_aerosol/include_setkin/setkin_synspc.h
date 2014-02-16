!=======================================================================
!
! $Id: setkin_synspc.h,v 1.1 2013-07-31 15:35:22 ssteenro Exp $
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
!  Setkin files generated:  Thu Feb 17 20:08:54 2011
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

