!=============================================================================
!
! $Id: gas_sulfchem.F90,v 1.6 2012-01-05 21:41:30 jkouatch Exp $
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   gas_sulfchem.F
!
! ROUTINES
!   Do_Gas_Sulfchem
!
! HISTORY
!   - December 8, 2005 * Bigyani Das
!     Added variables sfso4, snso4 for primary so4 emission
!     as 2.5% of SO2 emission and also in the SO4 chemical source
!     and sink terms add sfso4, snso4.  used dms_no3 array as the 
!     index for primary so4 emission budget
!     Added controlling variables do_aerocom and do_emiss_dust
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Gas_Sulfchem
!
! DESCRIPTION
!   This is the routine for the gas-chemistry
!
! ARGUMENTS
!   dtclc     : time step for clear air chemistry (s)
!   itloop    : # of zones (ilong * ilat * ivert)
!   cp        : species concentration (molecules/cm3 for gas-species,
!               and g/g for aerosols)
!   dcp       : species concentration change
!   cfac      : total cloud fraction [0 - 1]
!   zmp       : air density (molecules/cm^3)
!   semiss    : array of sulfur emissions (molecules/cm^3/s)
!   oh, ho2, h2o, no3
!             : oh, ho2, h2o, no3 species concentrations (molecules/cm3)
!   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p
!             : gas phase reaction rate coefficient
!
!-----------------------------------------------------------------------------

      subroutine Do_Gas_Sulfchem  &
     &  (dtclc, itloop, cp, dcp, cfac, zmp, semiss,  &
     &   oh, ho2, h2o, no3,  &
     &   massc, pr_sulf_src, do_aerocom,  &
     &   dms_oh, dms_no3, so2_oh,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      real*8,intent(in)  :: dtclc

      integer,intent(in) :: itloop

      real*8,intent(in)  :: cp (itloop, NSP)
      real*8,intent(Out) :: dcp(itloop, NSP)

      real*8,intent(in)  :: cfac(itloop)
      real*8,intent(in)  :: zmp (itloop)
      real*8,intent(in)  :: semiss(itloop, NSP)

      real*8,intent(in)  :: oh (itloop)
      real*8,intent(in)  :: ho2(itloop)
      real*8,intent(in)  :: h2o(itloop)
      real*8,intent(in)  :: no3(itloop)

      real*8 ,intent(in) :: massc(itloop)
      logical,intent(in) :: pr_sulf_src
      logical,intent(in) :: do_aerocom

      real*8,intent(inOut) :: dms_oh (itloop)
      real*8,intent(inOut) :: dms_no3(itloop)
      real*8,intent(inOut) :: so2_oh (itloop)

      real*8,intent(in) :: qh1p(itloop)
      real*8,intent(in) :: qh2p(itloop)
      real*8,intent(in) :: qh3p(itloop)
      real*8,intent(in) :: qj4p(itloop)
      real*8,intent(in) :: qh5p(itloop)
      real*8,intent(in) :: qh8p(itloop)
      real*8,intent(in) :: qh9p(itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8 :: exh2o2(itloop)
      real*8 :: ah2o2 (itloop)
      real*8 :: aldms (itloop)
      real*8 :: exdms (itloop)
      real*8 :: rdms  (itloop)
      real*8 :: also2 (itloop)
      real*8 :: exso2 (itloop)
      real*8 :: rso2  (itloop)

!.... compute the analytic integration terms
!
!          ======
      call chem10  &
!          ======
     &  (dtclc, itloop,  &
     &   oh, ho2, h2o, no3,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p,  &
     &   exh2o2, ah2o2,  &
     &   aldms, exdms, rdms,  &
     &   also2, exso2, rso2)

!.... compute species changes due to clear region chemistry over dtclc
!
!          ======
      call chem11  &
!          ======
     &  (dtclc, itloop, cp, dcp,  &
     &   oh, cfac, zmp, semiss,  &
     &   massc, pr_sulf_src, do_aerocom,  &
     &   dms_oh, dms_no3, so2_oh,  &
     &   exh2o2, ah2o2,  &
     &   aldms, exdms, rdms,  &
     &   also2, exso2, rso2)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chem10
!
! DESCRIPTION
!
!   compute the decay and "asymptotic" terms for the chemisty
!   analytic solutions in clear regions
!
!   dtclc  :  time step for clear air chemistry (s)
!   ohp, ho2p, h2op, no3p
!          :  oh, ho2, h2o, no3 species concentrations (molecules/cm3)
!   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p
!          :  gas phase reaction rate coefficient
!
!.... hydrogen peroxide
!     *     alh2o2 = qh4 + qh5 * OH
!     *     exh2o2 = exp( - alh2o2 * dtclc )
!     *     ah2o2  = (qh3 + qh8 * H2O) * HO2**2 / alh2o2
!
!.... dms
!     *     aldms = qh1 * OH + qh9 * NO3
!     *     exdms = exp( - aldms * dtclc )
!     *     rdms  = qh1 * OH / aldms
!
!.... sulfur dioxide
!     *     also2 = qh2 * OH
!     *     exso2 = exp( - also2 * dtclc )
!     *     rso2  = aldms / ( also2 - aldms )
!     *     aso2  = so2sor / also2                   not used
!
!----------------------------------------------------------------------

      subroutine chem10  &
     &  (dtclc, itloop,  &
     &   ohp, ho2p, h2op, no3p,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p,  &
     &   exh2o2, ah2o2,  &
     &   aldms, exdms, rdms,  &
     &   also2, exso2, rso2)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      real*8,intent(in)  :: dtclc

      integer,intent(in) :: itloop

      real*8,intent(in)  :: ohp (itloop)
      real*8,intent(in)  :: ho2p(itloop)
      real*8,intent(in)  :: h2op(itloop)
      real*8,intent(in)  :: no3p(itloop)

      real*8,intent(in) :: qh1p(itloop)
      real*8,intent(in) :: qh2p(itloop)
      real*8,intent(in) :: qh3p(itloop)
      real*8,intent(in) :: qj4p(itloop)
      real*8,intent(in) :: qh5p(itloop)
      real*8,intent(in) :: qh8p(itloop)
      real*8,intent(in) :: qh9p(itloop)

      real*8,intent(Out)  :: exh2o2(itloop)
      real*8,intent(Out)  :: ah2o2 (itloop)
      real*8,intent(Out)  :: aldms (itloop)
      real*8,intent(Out)  :: exdms (itloop)
      real*8,intent(Out)  :: rdms  (itloop)
      real*8,intent(Out)  :: also2 (itloop)
      real*8,intent(Out)  :: exso2 (itloop)
      real*8,intent(Out)  :: rso2  (itloop)

      real*8  :: alh2o2



      exh2o2 = 0.0d0
      ah2o2  = 0.0d0
      aldms  = 0.0d0
      exdms  = 0.0d0
      rdms   = 0.0d0
      also2  = 0.0d0
      exso2  = 0.0d0
      rso2   = 0.0d0

!     *****************************************************************
!     * (4) H2O2
!     *****************************************************************

      do 100 ijk = 1, itloop

        alh2o2 = qj4p(ijk) + qh5p(ijk) * ohp(ijk)
        exh2o2(ijk) = exp( - alh2o2 * dtclc )
        ah2o2(ijk)  = (qh3p(ijk) + qh8p(ijk) * h2op(ijk))  &
     &              * ho2p(ijk) * ho2p(ijk) / alh2o2

  100 continue

!     *****************************************************************
!     * (3) DMS
!     *****************************************************************

      where( ohp(:) .gt. 10. )
        aldms(:) = qh1p(:) * ohp(:) + qh9p(:) * no3p(:)
        exdms(:) = exp( - aldms(:) * dtclc )
        rdms (:) = qh1p(:) * ohp(:) / aldms(:)
      elsewhere
        aldms(:) = 0.
        exdms(:) = 0.
        rdms (:) = 0.
      end where


!     *****************************************************************
!     * (2) SO2
!     *****************************************************************

      do 300 ijk = 1, itloop

        also2(ijk) = qh2p(ijk) * ohp(ijk)
        exso2(ijk) = exp( - also2(ijk) * dtclc )
        rso2(ijk)  = aldms(ijk) / ( also2(ijk) - aldms(ijk) )

  300 continue

      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chem11
!
! DESCRIPTION
!
!   advances the gas phase species using analytic expressions
!
!     ********** comments:
!
!     * species changes are computed first and stored in dcp
!       these changes are added to those from aqueous chemistry,
!       then applied to the cp
!
!     *     cp  : species concentration (# cm-3 for gas-species)
!     *     dcp : storage for conc. change
!
!     *     The following species are categorized by the type
!           of source driving them.
!
!     *     fso4a      ==  H2SO4,  fossil - aqueous
!     *     fso4n      ==  H2SO4,  fossil - clear sky
!     *     fso2       ==  SO2  ,    "
!
!     *     nso4a      ==  H2SO4,  natural - aqueous
!     *     nso4n      ==  H2SO4,     "    - clear sky
!     *     nso2       ==  SO2  ,     "
!     *     ndms       ==  DMS  ,     "
!
!     *     bso4a      ==  H2SO4,  biomass - aqueous
!     *     bso4n      ==  H2SO4,     "    - clear sky
!     *     bso2       ==  SO2  ,     "
!
!     *     h2o2       ==  H2O2
!
!     *     also2  = qh2 * OH
!     *     exso2  = exp( - also2 * dtclc )
!     *     rso2   = aldms / (  also2 - aldms )
!     *     aldms  = qh1 * OH + qh9 * NO3
!     *     exdms  = exp( - aldms * dtclc )
!     *     rdms   = qh1 * OH / aldms
!     *     alh2o2 = qh4 + qh5 * OH
!     *     exh2o2 = exp( - alh2o2 * dtclc )
!     *     ah2o2  = qh3 * HO2**2 / alh2o2
!
!
!     ********** output:
!
!     *     Species changes due to clear sky reactions.
!
!     *     dfso4n    ==  H2SO4,  fossil  - clear sky
!     *     dfso2     ==  SO2  ,    "
!
!     *     dnso4n    ==  H2SO4,  natural - clear sky
!     *     dnso2     ==  SO2  ,     "
!     *     dndms     ==  DMS  ,     "
!
!     *     dbso4n    ==  H2SO4,  biomass - clear sky
!     *     dbso2     ==  SO2  ,     "
!
!     *     dh2o2     ==  H2O2
!
!     *     dcp       ==  clear air species changes, then
!                         the clear air contribution to
!                         species change, i.e.
!                         ( 1 - cfac ) * dcp.
!                         In chem30 there will be added to this
!                         the contribution due to aqueous chem
!                         in cloud, i.e.
!                         cfac * [ (new aqueous) - (old species) ].
!                         This will be imposed upon the grid points
!                         in chem40.
!
!-----------------------------------------------------------------------------

      subroutine chem11  &
     &  (dtclc, itloop, cp, dcp,  &
     &   ohp, cfac, zmp, semiss,  &
     &   massc, pr_sulf_src, do_aerocom,  &
     &   dms_oh, dms_no3, so2_oh,  &
     &   exh2o2, ah2o2,  &
     &   aldms, exdms, rdms,  &
     &   also2, exso2, rso2)

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

      real*8,intent(in)  :: dtclc

      integer,intent(in) :: itloop

      real*8,intent(in)  :: cp    (itloop, NSP)
      real*8,intent(Out) :: dcp   (itloop, NSP)

      real*8,intent(in)  :: ohp   (itloop)
      real*8,intent(in)  :: cfac   (itloop)
      real*8,intent(in)  :: zmp   (itloop)
      real*8,intent(in)  :: semiss(itloop, NSP)

      real*8,intent(in)  :: massc (itloop)
      logical,intent(in) :: pr_sulf_src
      logical,intent(in) :: do_aerocom

      real*8,intent(inOut) :: dms_oh (itloop)
      real*8,intent(inOut) :: dms_no3(itloop)
      real*8,intent(inOut) :: so2_oh (itloop)

      real*8,intent(in)  :: exh2o2(itloop)
      real*8,intent(in)  :: ah2o2 (itloop)

      real*8,intent(in)  :: aldms (itloop)
      real*8,intent(in)  :: exdms (itloop)
      real*8,intent(in)  :: rdms  (itloop)

      real*8,intent(in)  :: also2 (itloop)
      real*8,intent(in)  :: exso2 (itloop)
      real*8,intent(in)  :: rso2  (itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8 :: h2o2  (itloop)
      real*8 :: dh2o2 (itloop)
      real*8 :: ndms  (itloop)
      real*8 :: dndms (itloop)
      real*8 :: nso2  (itloop)
      real*8 :: dnso2 (itloop)
      real*8 :: dnso4n(itloop)
      real*8 :: fso2  (itloop)
      real*8 :: dfso2 (itloop)
      real*8 :: dfso4n(itloop)

      real*8 :: sfso2(itloop)
      real*8 :: snso2(itloop)
      real*8 :: sndms(itloop)

      real*8 :: sfso4(itloop)    ! primary so4
      real*8 :: snso4(itloop)    ! primary so4

      h2o2   = 0.0d0
      dh2o2  = 0.0d0
      ndms   = 0.0d0
      dndms  = 0.0d0
      nso2   = 0.0d0
      dnso2  = 0.0d0
      dnso4n = 0.0d0
      fso2   = 0.0d0
      dfso2  = 0.0d0
      dfso4n = 0.0d0

      sfso2  = 0.0d0
      snso2  = 0.0d0
      sndms  = 0.0d0

      sfso4  = 0.0d0
      snso4  = 0.0d0

      h2o2(:) = cp(:, IH2O2_1)
      ndms(:) = cp(:, INDMS_l)
      nso2(:) = cp(:, INSO2_l)
      fso2(:) = cp(:, IFSO2_l)

      if (do_aerocom) then
         sfso2(:) = semiss(:, IFSO2_l) * 0.975d0
         snso2(:) = semiss(:, INSO2_l) * 0.975d0
         sndms(:) = semiss(:, INDMS_l)

         sfso4(:) = semiss(:, IFSO2_l) * 0.025d0
         snso4(:) = semiss(:, INSO2_l) * 0.025d0
      else
         sfso2(:) = semiss(:, IFSO2_l)
         snso2(:) = semiss(:, INSO2_l)
         sndms(:) = semiss(:, INDMS_l)   ! no scaling
!     sndms(:) = semiss(:, INDMS_l) /2.0d0   ! scale by 2 to reduce DMS emission
      endif

!     *****************************************************************
!     * (11) H2O2
!     *****************************************************************


      dh2o2(:) = ( ah2o2(:) - h2o2(:) )  &
     &         * ( 1.d0 - exh2o2(:) )
      where ( abs( dh2o2(:) ) < 0.d0 ) dh2o2(:) = 0.d0


!     *****************************************************************
!     * (10) biomass SO2
!     *****************************************************************
!
!           dbso2(:) = ( sbso2(:) / also2(:) - bso2(:) )
!     &              * ( 1.d0 - exso2(:) )
!           where ( abs( dbso2(:) ) < 1. ) dbso2(:) = 0.d0
!
!     *****************************************************************
!     * (9) biomass H2SO4
!     *****************************************************************
!
!           dbso4n(:) = sbso2(:) * dtclc - dbso2(:)
!           where ( abs( dbso4n(:) ) < 1.d0 ) dbso4n(:) = 0.d0
!


!     *****************************************************************
!     * (7) natural DMS
!     *****************************************************************

      where ( ohp(:) .gt. 10.d0 )

        dndms(:) = ( sndms(:) / aldms(:) - ndms(:) )  &
     &             * ( 1.d0 - exdms(:) )

      elsewhere

        dndms(:) = sndms(:) * dtclc

      end where

      where ( abs( dndms(:) ) < 0.d0 ) dndms(:) = 0.d0


!     *****************************************************************
!     * (6) natural SO2
!     *****************************************************************

      where ( ohp(:) .gt. 10.d0 )

        dnso2(:) = ( ( sndms(:) + snso2(:) )  &
     &             / also2(:) - nso2(:) )  &
     &             * ( 1.d0 - exso2(:) )  &
     &             - ( sndms(:) / aldms(:) -  ndms(:) )  &
     &             * ( exdms(:) - exso2(:) ) * rso2(:)

      elsewhere

        dnso2(:) = ( snso2(:) / also2(:) - nso2(:) )  &
     &             * ( 1.d0 - exso2(:) )

      end where

      where ( abs( dnso2(:) ) < 0.d0 ) dnso2(:) = 0.d0


!     *****************************************************************
!     * (5) natural H2SO4
!     *****************************************************************

      if (do_aerocom) then
         where ( ohp(:) .gt. 10.d0 )

            dnso4n(:) = ( sndms(:) + snso2(:) ) * dtclc  &
     &                - dnso2(:) - dndms(:) + snso4(:) * dtclc

         elsewhere

            dnso4n(:) = snso2(:) * dtclc - dnso2(:) + snso4(:) * dtclc

         end where
      else
         where ( ohp(:) .gt. 10.d0 )

              dnso4n(:) = ( sndms(:) + snso2(:) ) * dtclc  &
     &                - dnso2(:) - dndms(:)

          elsewhere

              dnso4n(:) = snso2(:) * dtclc - dnso2(:)

          end where
      endif

      where ( abs( dnso4n(:) ) < 0.d0 ) dnso4n(:) = 0.d0


!     *****************************************************************
!     * (3) fossil SO2
!     *****************************************************************

      dfso2(:) = ( sfso2(:) / also2(:) - fso2(:) )  &
     &           * ( 1.d0 - exso2(:) )

      where ( abs( dfso2(:) ) < 0.d0 ) dfso2(:) = 0.d0

!     *****************************************************************
!     * (4) fossil H2SO4
!     *****************************************************************

      if (do_aerocom) then
        dfso4n(:) = sfso2(:) * dtclc - dfso2(:) + sfso4(:) * dtclc
      else
        dfso4n(:) = sfso2(:) * dtclc - dfso2(:)
      end if

!.sds      where ( abs( dfso4n(:) ) < 1.d0 ) dfso4n(:) = 0.d0
      where ( abs( dfso4n(:) ) < 0.d0 ) dfso4n(:) = 0.d0


!     *****************************************************************
!     * Compute the clear sky contribution to species change.
!     *****************************************************************

!     conv: convert SO4 concentration from # cm-3 to gm gm-1
      conv = wtmair / wtmso4

!.... H2O2
      dcp(:, IH2O2_1)   = ( 1.0d0 - cfac(:) ) * dh2o2(:)

!.... fossil
      dcp(:, IFSO2_l)   = ( 1.0d0 - cfac(:) ) * dfso2(:)

      dcp(:, IFSO4N1_l) = ( 1.0d0 - cfac(:) ) * dfso4n(:)  &
     &                    * ffra_so4(1) / (conv * zmp(:))  ! unit: g/g
      dcp(:, IFSO4N2_l) = ( 1.0d0 - cfac(:) ) * dfso4n(:)  &
     &                    * ffra_so4(2) / (conv * zmp(:))
      dcp(:, IFSO4N3_l) = ( 1.0d0 - cfac(:) ) * dfso4n(:)  &
     &                    * ffra_so4(3) / (conv * zmp(:))

!.... natural
      dcp(:, INDMS_l)   = ( 1.0d0 - cfac(:) ) * dndms(:)
      dcp(:, INSO2_l)   = ( 1.0d0 - cfac(:) ) * dnso2(:)

      dcp(:, INSO4N1_l) = ( 1.0d0 - cfac(:) ) * dnso4n(:)  &
     &                    * nfra_so4(1) / (conv * zmp(:))  ! unit: g/g
      dcp(:, INSO4N2_l) = ( 1.0d0 - cfac(:) ) * dnso4n(:)  &
     &                    * nfra_so4(2) / (conv * zmp(:))
      dcp(:, INSO4N3_l) = ( 1.0d0 - cfac(:) ) * dnso4n(:)  &
     &                    * nfra_so4(3) / (conv * zmp(:))

!.... biomass
!     dbso2(:)  = ( 1.0d0 - cfac(:) ) * dbso2(:)
!     dbso4n(:) = ( 1.0d0 - cfac(:) ) * dbso4n(:)


      if (pr_sulf_src) then

        dms_oh(:) =  &
     &    dms_oh(:) +  &
     &    ( 1.0d0 - cfac(:) ) *  &
     &    max( (sndms(:) * dtclc - dndms(:)), 0.0d0) *  &
     &    rdms(:) *  &
     &    massc(:) / zmp(:) * wtmso2 / wtmair

        if (do_aerocom) then
          dms_no3(:) =  &
     &      dms_no3(:) +  &
     &      ( 1.0d0 - cfac(:) ) *  &
     &      ( (sfso4(:) + snso4(:)) * dtclc ) *  &
     &      massc(:) / zmp(:) * wtmso4 / wtmair

          so2_oh(:) =  &
     &      so2_oh(:) +  &
     &      ( 1.0d0 - cfac(:) ) *  &
     &      (dfso4n(:) + dnso4n(:) - (sfso4(:) + snso4(:)) * dtclc) *  &
     &      massc(:) / zmp(:) * wtmso4 / wtmair

        else
          dms_no3(:) =  &
     &      dms_no3(:) +  &
     &      ( 1.0d0 - cfac(:) ) *  &
     &      max( (sndms(:) * dtclc - dndms(:)), 0.0d0) *  &
     &      ( 1.0d0 - rdms(:)) *  &
     &      massc(:) / zmp(:) * wtmso2 / wtmair

          so2_oh(:) =  &
     &      so2_oh(:) +  &
     &      ( 1.0d0 - cfac(:) ) *  &
     &      (dfso4n(:) + dnso4n(:)) *  &
     &      massc(:) / zmp(:) * wtmso4 / wtmair

        end if

      end if


      return

      end

