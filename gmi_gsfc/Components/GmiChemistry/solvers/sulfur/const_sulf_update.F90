!=============================================================================
!
! $Id: const_sulf_update.F90,v 1.6 2012-01-05 21:41:30 jkouatch Exp $
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   const_sulf_update.F
!
! ROUTINES
!   Do_Const_Sulf_Update
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Const_Sulf_Update (chem40)
!
! DESCRIPTION
!
!     *****************************************************************
!     * advances the species concentrations and does mixing of clear
!     *   and cloudy regions
!     *****************************************************************
!
!     ********** input:
!
!     *     cp     : concentration (# cm-3 for gas-species, g/g for aerosols)
!     *     dcp    : storage for conc. change
!
!     *     fso4a  : aqueous fossil H2SO4 conc.
!     *     fso4n  : clear sky fossil H2SO4 conc.
!     *     fso2   : fossil SO2 conc.
!
!     *     nso4a  : aqueous natural H2SO4 conc.
!     *     nso4n  : clear sky natural H2SO4 conc.
!     *     nso2   : natural SO2 conc.
!     *     ndms   : natural DMS conc.
!
!     *     h2o2   : H2O2 conc.
!
!     *     dfso4a : aqueous fossil H2SO4 conc.
!     *     dfso4n : clear sky fossil H2SO4 conc.
!     *     dfso2  : fossil SO2 conc.
!
!     *     dnso4a : aqueous natural H2SO4 conc.
!     *     dnso4n : clear sky natural H2SO4 conc.
!     *     dnso2  : natural SO2 conc.
!     *     dndms  : natural DMS conc.
!
!     *     dh2o2  : H2O2 conc.
!
!
!     * For H2SO4 the aqueous and clear sky parts remain seperate.
!       However, the H2O2, DMS and SO2 are a combination of the aqueous
!       and clear sky parts, weighted by the cloud fraction cfp.
!
!     * The clear sky changes were computed and stored in chem11.
!       In chem30 the aqueous changes were computed and,
!       except for H2SO4, added to the clear sky.  e.g., with cfp
!       the cloud fraction at a grid point
!
!       dcp(ns) = (1-cfp) * d[clear sky] + cfp * d[aqueous]
!
!----------------------------------------------------------------------

      subroutine Do_Const_Sulf_Update  &
     &  (itloop, cp, dcp)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

#ifdef nonZeroInd
      integer, parameter :: IH2O2_1    = IH2O2
      integer, parameter :: IFSO2_l    = 1
      integer, parameter :: INSO2_l    = 1
      integer, parameter :: INDMS_l    = 1
      integer, parameter :: IFSO4A_l   = 1
      integer, parameter :: INSO4A_l   = 1
      integer, parameter :: IFSO4N1_l  = 1
      integer, parameter :: IFSO4N2_l  = 1
      integer, parameter :: IFSO4N3_l  = 1
      integer, parameter :: INSO4N1_l  = 1
      integer, parameter :: INSO4N2_l  = 1
      integer, parameter :: INSO4N3_l  = 1
#elif nonZeroInd_tracers
      integer, parameter :: IH2O2_1    = 1
      integer, parameter :: IFSO2_l    = 1
      integer, parameter :: INSO2_l    = 1
      integer, parameter :: INDMS_l    = 1
      integer, parameter :: IFSO4A_l   = 1
      integer, parameter :: INSO4A_l   = 1
      integer, parameter :: IFSO4N1_l  = 1
      integer, parameter :: IFSO4N2_l  = 1
      integer, parameter :: IFSO4N3_l  = 1
      integer, parameter :: INSO4N1_l  = 1
      integer, parameter :: INSO4N2_l  = 1
      integer, parameter :: INSO4N3_l  = 1
#else
      integer, parameter :: IH2O2_1    = IH2O2
      integer, parameter :: IFSO2_l    = IFSO2
      integer, parameter :: INSO2_l    = INSO2
      integer, parameter :: INDMS_l    = INDMS
      integer, parameter :: IFSO4A_l   = IFSO4A
      integer, parameter :: INSO4A_l   = INSO4A
      integer, parameter :: IFSO4N1_l  = IFSO4N1
      integer, parameter :: IFSO4N2_l  = IFSO4N2
      integer, parameter :: IFSO4N3_l  = IFSO4N3
      integer, parameter :: INSO4N1_l  = INSO4N1
      integer, parameter :: INSO4N2_l  = INSO4N2
      integer, parameter :: INSO4N3_l  = INSO4N3
#endif


!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: cp (itloop, NSP)
      real*8  :: dcp(itloop, NSP)


      cp(:, IH2O2_1)     = cp(:, IH2O2_1) + dcp(:,IH2O2_1)
      cp(:, INDMS_l)   = cp(:, INDMS_l) + dcp(:,INDMS_l)
      cp(:, INSO2_l)   = cp(:, INSO2_l) + dcp(:,INSO2_l)
      cp(:, IFSO2_l)   = cp(:, IFSO2_l) + dcp(:,IFSO2_l)

      cp(:, IFSO4A_l)  = cp(:, IFSO4A_l ) + dcp(:,IFSO4A_l)
      cp(:, IFSO4N1_l) = cp(:, IFSO4N1_l) + dcp(:,IFSO4N1_l)
      cp(:, IFSO4N2_l) = cp(:, IFSO4N2_l) + dcp(:,IFSO4N2_l)
      cp(:, IFSO4N3_l) = cp(:, IFSO4N3_l) + dcp(:,IFSO4N3_l)

      cp(:, INSO4A_l)  = cp(:, INSO4A_l ) + dcp(:,INSO4A_l)
      cp(:, INSO4N1_l) = cp(:, INSO4N1_l) + dcp(:,INSO4N1_l)
      cp(:, INSO4N2_l) = cp(:, INSO4N2_l) + dcp(:,INSO4N2_l)
      cp(:, INSO4N3_l) = cp(:, INSO4N3_l) + dcp(:,INSO4N3_l)

      where(cp(:, IH2O2_1)      < 1.0d-30) cp(:, IH2O2_1) = 1.0d-30
      where(cp(:, INDMS_l)    < 1.0d-30) cp(:, INDMS_l) = 1.0d-30
      where(cp(:, INSO2_l)    < 1.0d-30) cp(:, INSO2_l) = 1.0d-30
      where(cp(:, IFSO2_l)    < 1.0d-30) cp(:, IFSO2_l) = 1.0d-30

      where(cp(:, IFSO4A_l)   < 1.0d-30) cp(:, IFSO4A_l ) = 1.0d-30
      where(cp(:, IFSO4N1_l)  < 1.0d-30) cp(:, IFSO4N1_l) = 1.0d-30
      where(cp(:, IFSO4N2_l)  < 1.0d-30) cp(:, IFSO4N2_l) = 1.0d-30
      where(cp(:, IFSO4N3_l)  < 1.0d-30) cp(:, IFSO4N3_l) = 1.0d-30

      where(cp(:, INSO4A_l)   < 1.0d-30) cp(:, INSO4A_l ) = 1.0d-30
      where(cp(:, INSO4N1_l)  < 1.0d-30) cp(:, INSO4N1_l) = 1.0d-30
      where(cp(:, INSO4N2_l)  < 1.0d-30) cp(:, INSO4N2_l) = 1.0d-30
      where(cp(:, INSO4N3_l)  < 1.0d-30) cp(:, INSO4N3_l) = 1.0d-30


      return

      end

