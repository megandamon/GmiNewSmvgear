!=======================================================================
!
! $Id: setkin_kcalc.F90,v 1.2 2011-08-08 19:37:04 mrdamon Exp $
!
! ROUTINE
!   kcalc - GMI model (setkin_kcalc.F)
!   1 AUG 02 - PSC
!
! DESCRIPTION
!   Calculates and returns rate constants for thermal
!   reactions for the temperatures and pressures supplied
!
! ARGUMENTS
!  INPUT
!   pressure    : mb - profile of pressures
!   temperature : K - profile of temperatures
!   adcol       : molecules cm-3 - profile of total number density
!   specarr     : molecules cm-3 - profiles of species concentrations
!  OUTPUT
!   rcarr       : cm3 molecule-1 s-1 - rate constant values in units as
!                 appropriate
!
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
!bnd
!bnd Updated to JPL2002 (JPL Publication 02-25) in April 2004.
!bnd Updated only trace gas rate constants.  No new rxns were added
!bnd as this would require the IDL code used to generate setkin.
!bnd The O1D + N2/H2O/O2/N2O/O3 rxns are from Dunlea & Ravishankara (2004)
!bnd and Ravishankara et al. (2004).
!
!=======================================================================
      subroutine kcalc( npres0 ,sadcol  &
     &                 ,pressure ,tropp ,temperature ,lwc ,adcol  &
     &                 ,specarr ,rcarr)

      implicit none

#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "setkin_mw.h"

!.... Argument declarations

      INTEGER, INTENT (IN)  :: npres0

      REAL*8,  INTENT (IN)  :: tropp

      REAL*8,  INTENT (IN)  :: adcol       (       npres0)
      REAL*8,  INTENT (IN)  :: lwc         (       npres0)
      REAL*8,  INTENT (IN)  :: pressure    (       npres0)
      REAL*8,  INTENT (IN)  :: temperature (       npres0)

      REAL*8,  INTENT (IN)  :: sadcol      (NSAD  ,npres0)
      REAL*8,  INTENT (IN)  :: specarr     (NMF   ,npres0)

      REAL*8,  INTENT (OUT) :: rcarr       (NUM_K ,npres0)

!.... Local variable declarations

      real*8  mw(NSP)

      real*8  &
     &  water    (npres0)

      real*8  &
     &  sad_ice  (npres0)  &
     & ,sad_lbs  (npres0)  &
     & ,sad_nat  (npres0)  &
     & ,sad_soot (npres0)  &
     & ,sad_sts  (npres0)

     mw(:) = mw_data(:)

!....          Extract water number density
!
      water(:)    = specarr(12 ,:)

      sad_lbs(:)  = sadcol(1, :)
      sad_sts(:)  = sadcol(2, :)
      sad_nat(:)  = sadcol(3, :)
      sad_ice(:)  = sadcol(4, :)
      sad_soot(:) = sadcol(5, :)
!
!....          Start thermal rate constants
!
!....           O + O2 = O3
!
      rcarr(1,:) = skterlp( 6.000D-34 ,2.40D+00 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           O + O3 = 2 O2
!
      rcarr(2,:) = skarr( 8.000D-12 ,2060.0D+00 ,temperature)
!
!....           N2 + O1D = N2 + O
!
!bnd      rcarr(3,:) = skarr( 1.800D-11 ,-110.0D+00 ,temperature)
      rcarr(3,:) = skarr( 2.100D-11 ,-115.0D+00 ,temperature)
!
!....           O1D + O2 = O + O2
!
!bnd      rcarr(4,:) = skarr( 3.200D-11 ,-70.0D+00 ,temperature)
      rcarr(4,:) = skarr( 3.1200D-11 ,-70.0D+00 ,temperature)
!
!....           O1D + O3 = 2 O2
!
!bnd      rcarr(5,:) = skarr( 1.200D-10 ,0.0D+00 ,temperature)
      rcarr(5,:) = skarr( 2.3700D-10 ,-6.0D+00 ,temperature)
!
!....           H2O + O1D = 2 OH
!
!bnd      rcarr(6,:) = skarr( 2.200D-10 ,0.0D+00 ,temperature)
      rcarr(6,:) = skarr( 1.620D-10 ,-65.0D+00 ,temperature)
!
!....           H2 + O1D = H + OH
!
      rcarr(7,:) = skarr( 1.100D-10 ,0.0D+00 ,temperature)
!
!....           N2O + O1D = N2 + O2
!
!bnd Took percent reaction product pathways from JPL 02-25, but
!bnd reaction rate from Dunlea and Ravishankara (2004).
!bnd      rcarr(8,:) = skarr( 4.900D-11 ,0.0D+00 ,temperature)
      rcarr(8,:) = skarr( 1.110D-10 ,-17.0D+00 ,temperature)
      rcarr(8,:) = rcarr(8,:) * 0.422
!
!....           N2O + O1D = 2 NO
!
!bnd      rcarr(9,:) = skarr( 6.700D-11 ,0.0D+00 ,temperature)
      rcarr(9,:) = skarr( 1.110D-10 ,-17.0D+00 ,temperature)
      rcarr(9,:) = rcarr(9,:) * 0.578
!
!....           CH4 + O1D = CH3O2 + OH
!
      rcarr(10,:) = skarr( 1.125D-10 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H + HO2
!
      rcarr(11,:) = skarr( 3.000D-11 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H2
!
      rcarr(12,:) = skarr( 7.500D-12 ,0.0D+00 ,temperature)
!
!....           CF2Cl2 + O1D = 2 Cl
!
      rcarr(13,:) = skarr( 1.200D-10 ,0.0D+00 ,temperature)
!
!....           CFC113 + O1D = 3 Cl
!
      rcarr(14,:) = skarr( 1.500D-10 ,0.0D+00 ,temperature)
!
!....           CFC114 + O1D = 2 Cl
!
!bnd      rcarr(15,:) = skarr( 1.000D-10 ,0.0D+00 ,temperature)
      rcarr(15,:) = skarr( 1.300D-10 ,0.0D+00 ,temperature)
!
!....           CFC115 + O1D = Cl
!
      rcarr(16,:) = skarr( 1.500D-11 ,0.0D+00 ,temperature)
!
!....           HCFC22 + O1D = Cl
!
      rcarr(17,:) = skarr( 7.200D-11 ,0.0D+00 ,temperature)
!
!....           HCFC141b + O1D = 2 Cl
!
      rcarr(18,:) = skarr( 1.790D-10 ,0.0D+00 ,temperature)
!
!....           HCFC142b + O1D = Cl
!
      rcarr(19,:) = skarr( 1.630D-10 ,0.0D+00 ,temperature)
!
!....           H + O2 = HO2
!
      rcarr(20,:) = sktroe( 5.700D-32 ,1.60D0  &
     &                     ,7.500D-11 ,0.00D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           H + O3 = O2 + OH
!
      rcarr(21,:) = skarr( 1.400D-10 ,470.0D+00 ,temperature)
!
!....           H2 + OH = H + H2O
!
      rcarr(22,:) = skarr( 5.500D-12 ,2000.0D+00 ,temperature)
!
!....           O3 + OH = HO2 + O2
!
!bnd      rcarr(23,:) = skarr( 1.500D-12 ,880.0D+00 ,temperature)
      rcarr(23,:) = skarr( 1.700D-12 ,940.0D+00 ,temperature)
!
!....           O + OH = H + O2
!
      rcarr(24,:) = skarr( 2.200D-11 ,-120.0D+00 ,temperature)
!
!....           OH + OH = H2O + O
!
      rcarr(25,:) = skarr( 4.200D-12 ,240.0D+00 ,temperature)
!
!....           HO2 + O = O2 + OH
!
      rcarr(26,:) = skarr( 3.000D-11 ,-200.0D+00 ,temperature)
!
!....           HO2 + O3 = 2 O2 + OH
!
!bnd      rcarr(27,:) = skarr( 2.000D-14 ,680.0D+00 ,temperature)
      rcarr(27,:) = skarr( 1.000D-14 ,490.0D+00 ,temperature)
!
!....           H + HO2 = 2 OH
!
!bnd      rcarr(28,:) = skarr( 7.000D-11 ,0.0D+00 ,temperature)
      rcarr(28,:) = skarr( 8.100D-11 ,0.0D+00 ,temperature)
!
!....           HO2 + OH = H2O + O2
!
      rcarr(29,:) = skarr( 4.800D-11 ,-250.0D+00 ,temperature)
!
!....           HO2 + HO2 = H2O2 + O2
!
      rcarr(30,:) = skho2dis (temperature ,adcol)
!
!....           H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      rcarr(31,:) = skho2h2o (temperature ,adcol)
!
!....           H2O2 + OH = H2O + HO2
!
      rcarr(32,:) = skarr( 2.900D-12 ,160.0D+00 ,temperature)
!
!....           N + O2 = NO + O
!
      rcarr(33,:) = skarr( 1.500D-11 ,3600.0D+00 ,temperature)
!
!....           N + NO = N2 + O
!
      rcarr(34,:) = skarr( 2.100D-11 ,-100.0D+00 ,temperature)
!
!....           NO + O3 = NO2 + O2
!
      rcarr(35,:) = skarr( 3.000D-12 ,1500.0D+00 ,temperature)
!
!....           NO2 + OH = HNO3
!
!bnd      rcarr(36,:) = sktroe( 2.400D-30 ,3.10D0
!bnd     &                     ,1.700D-11 ,2.10D0 ,0.0D0
!bnd     &                     ,temperature ,adcol)
      rcarr(36,:) = sktroe( 2.000D-30 ,3.00D0  &
     &                     ,2.500D-11 ,0D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           HO2 + NO = NO2 + OH
!
      rcarr(37,:) = skarr( 3.500D-12 ,-250.0D+00 ,temperature)
!
!....           NO2 + O = NO + O2
!
      rcarr(38,:) = skarr( 5.600D-12 ,-180.0D+00 ,temperature)
!
!....           NO2 + O3 = NO3 + O2
!
      rcarr(39,:) = skarr( 1.200D-13 ,2450.0D+00 ,temperature)
!
!....           HO2 + NO2 = HO2NO2
!
      rcarr(40,:) = sktroe( 1.800D-31 ,3.20D0  &
     &                     ,4.700D-12 ,1.40D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           NO3 + O = NO2 + O2
!
      rcarr(41,:) = skarr( 1.000D-11 ,0.0D+00 ,temperature)
!
!....           NO + NO3 = 2 NO2
!
      rcarr(42,:) = skarr( 1.500D-11 ,-170.0D+00 ,temperature)
!
!....           NO2 + NO3 = N2O5
!
      rcarr(43,:) = sktroe( 2.000D-30 ,4.40D0  &
     &                     ,1.400D-12 ,0.70D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           N2O5 = NO2 + NO3
!
      rcarr(44,:) = sktroe( 6.667D-04 ,4.40D0  &
     &                     ,4.667D+14 ,0.70D0 ,10991.0D+00  &
     &                     ,temperature ,adcol)
!
!....           HNO3 + OH = H2O + NO3
!
      rcarr(45,:) = skohhno3 (temperature ,adcol)
!
!....           HO2NO2 = HO2 + NO2
!
      rcarr(46,:) = sktroe( 8.570D-05 ,3.20D0  &
     &                     ,2.240D+15 ,1.40D0 ,10900.0D+00  &
     &                     ,temperature ,adcol)
!
!....           HO2NO2 + OH = H2O + NO2 + O2
!
      rcarr(47,:) = skarr( 1.300D-12 ,-380.0D+00 ,temperature)
!
!....           Cl + O3 = ClO + O2
!
      rcarr(48,:) = skarr( 2.300D-11 ,200.0D+00 ,temperature)
!
!....           Cl + H2 = H + HCl
!
      rcarr(49,:) = skarr( 3.700D-11 ,2300.0D+00 ,temperature)
!
!....           Cl + H2O2 = HCl + HO2
!
      rcarr(50,:) = skarr( 1.100D-11 ,980.0D+00 ,temperature)
!
!....           Cl + HO2 = HCl + O2
!
      rcarr(51,:) = skarr( 1.800D-11 ,-170.0D+00 ,temperature)
!
!....           Cl + HO2 = ClO + OH
!
      rcarr(52,:) = skarr( 4.100D-11 ,450.0D+00 ,temperature)
!
!....           ClO + O = Cl + O2
!
      rcarr(53,:) = skarr( 3.000D-11 ,-70.0D+00 ,temperature)
!
!....           ClO + OH = Cl + HO2
!
      rcarr(54,:) = skarr( 7.400D-12 ,-270.0D+00 ,temperature)
!
!....           ClO + OH = HCl + O2
!
!bnd      rcarr(55,:) = skarr( 3.200D-13 ,-320.0D+00 ,temperature)
      rcarr(55,:) = skarr( 6.000D-13 ,-230.0D+00 ,temperature)
!
!....           ClO + HO2 = HOCl + O2
!
!bnd      rcarr(56,:) = skarr( 4.800D-13 ,-700.0D+00 ,temperature)
      rcarr(56,:) = skarr( 2.700D-12 ,-220.0D+00 ,temperature)
!
!....           ClO + HO2 = HCl + O3
!
      rcarr(57,:) = skarr( 0.000D+00 ,0.0D+00 ,temperature)
!
!....           ClO + NO = Cl + NO2
!
      rcarr(58,:) = skarr( 6.400D-12 ,-290.0D+00 ,temperature)
!
!....           ClO + NO2 = ClONO2
!
      rcarr(59,:) = sktroe( 1.800D-31 ,3.40D0  &
     &                     ,1.500D-11 ,1.90D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           ClO + ClO = 2 Cl + O2
!
      rcarr(60,:) = skarr( 3.000D-11 ,2450.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2 + O2
!
      rcarr(61,:) = skarr( 1.000D-12 ,1590.0D+00 ,temperature)
!
!....           ClO + ClO = Cl + OClO
!
      rcarr(62,:) = skarr( 3.500D-13 ,1370.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2O2
!
!bnd      rcarr(63,:) = sktroe( 2.200D-32 ,3.10D0
!bnd     &                     ,3.400D-12 ,1.00D0 ,0.0D0
!bnd     &                     ,temperature ,adcol)
      rcarr(63,:) = sktroe( 1.600D-32 ,4.50D0  &
     &                     ,2.000D-12 ,2.40D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           Cl2O2 = 2 ClO
!
      rcarr(64,:) = sktroe( 1.730D-05 ,3.10D0  &
     &                     ,2.680D+15 ,1.00D0 ,8744.0D+00  &
     &                     ,temperature ,adcol)
!
!....           HCl + OH = Cl + H2O
!
      rcarr(65,:) = skarr( 2.600D-12 ,350.0D+00 ,temperature)
!
!....           HOCl + OH = ClO + H2O
!
      rcarr(66,:) = skarr( 3.000D-12 ,500.0D+00 ,temperature)
!
!....           ClONO2 + O = ClO + NO3
!
!bnd      rcarr(67,:) = skarr( 4.500D-12 ,900.0D+00 ,temperature)
      rcarr(67,:) = skarr( 2.900D-12 ,800.0D+00 ,temperature)
!
!....           ClONO2 + OH = HOCl + NO3
!
      rcarr(68,:) = skarr( 1.200D-12 ,330.0D+00 ,temperature)
!
!....           Cl + ClONO2 = Cl2 + NO3
!
      rcarr(69,:) = skarr( 6.500D-12 ,-135.0D+00 ,temperature)
!
!....           Br + O3 = BrO + O2
!
      rcarr(70,:) = skarr( 1.700D-11 ,800.0D+00 ,temperature)
!
!....           Br + HO2 = HBr + O2
!
      rcarr(71,:) = skarr( 1.500D-11 ,600.0D+00 ,temperature)
!
!....           Br + CH2O = CO + HBr + HO2
!
      rcarr(72,:) = skarr( 1.700D-11 ,800.0D+00 ,temperature)
!
!....           BrO + O = Br + O2
!
      rcarr(73,:) = skarr( 1.900D-11 ,-230.0D+00 ,temperature)
!
!....           BrO + HO2 = HOBr + O2
!
      rcarr(74,:) = skarr( 3.400D-12 ,-540.0D+00 ,temperature)
!
!....           BrO + NO = Br + NO2
!
      rcarr(75,:) = skarr( 8.800D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + NO2 = BrONO2
!
      rcarr(76,:) = sktroe( 5.200D-31 ,3.20D0  &
     &                     ,6.900D-12 ,2.90D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           BrO + ClO = Br + OClO
!
      rcarr(77,:) = skarr( 9.500D-13 ,-550.0D+00 ,temperature)
!
!....           BrO + ClO = Br + Cl + O2
!
      rcarr(78,:) = skarr( 2.300D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + ClO = BrCl + O2
!
      rcarr(79,:) = skarr( 4.100D-13 ,-290.0D+00 ,temperature)
!
!....           BrO + BrO = 2 Br + O2
!
      rcarr(80,:) = skbrodis (temperature)
!
!....           HBr + OH = Br + H2O
!
      rcarr(81,:) = skarr( 1.100D-11 ,0.0D+00 ,temperature)
!
!....           CO + OH = H
!
      rcarr(82,:) = skcooh (temperature ,pressure)
!
!....           CH4 + OH = CH3O2 + H2O
!
      rcarr(83,:) = skarr( 2.450D-12 ,1775.0D+00 ,temperature)
!
!....           CH2O + OH = CO + H2O + HO2
!
!bnd      rcarr(84,:) = skarr( 1.000D-11 ,0.0D+00 ,temperature)
      rcarr(84,:) = skarr( 9.000D-12 ,0.0D+00 ,temperature)
!
!....           CH2O + O = CO + HO2 + OH
!
      rcarr(85,:) = skarr( 3.400D-11 ,1600.0D+00 ,temperature)
!
!....           CH4 + Cl = CH3O2 + HCl
!
      rcarr(86,:) = skarr( 9.600D-12 ,1360.0D+00 ,temperature)
!
!....           CH2O + Cl = CO + HCl + HO2
!
      rcarr(87,:) = skarr( 8.100D-11 ,30.0D+00 ,temperature)
!
!....           CH3O2 + NO = CH2O + HO2 + NO2
!
!bnd      rcarr(88,:) = skarr( 3.000D-12 ,-280.0D+00 ,temperature)
      rcarr(88,:) = skarr( 2.800D-12 ,-300.0D+00 ,temperature)
!
!....           CH3O2 + HO2 = CH3OOH + O2
!
!bnd      rcarr(89,:) = skarr( 3.800D-13 ,-800.0D+00 ,temperature)
      rcarr(89,:) = skarr( 4.100D-13 ,-750.0D+00 ,temperature)
!
!....           CH3OOH + OH = CH3O2 + H2O
!
!bnd      rcarr(90,:) = skarr( 2.700D-12 ,-200.0D+00 ,temperature)
      rcarr(90,:) = skarr( 3.800D-12 ,-200.0D+00 ,temperature)
!
!....           CH3Cl + OH = Cl + H2O + HO2
!
      rcarr(91,:) = skarr( 2.400D-12 ,1250.0D+00 ,temperature)
!
!....           CH3CCl3 + OH = 3 Cl + H2O
!
      rcarr(92,:) = skarr( 1.600D-12 ,1520.0D+00 ,temperature)
!
!....           HCFC22 + OH = Cl + H2O
!
      rcarr(93,:) = skarr( 1.050D-12 ,1600.0D+00 ,temperature)
!
!....           HCFC141b + OH = 2 Cl + H2O
!
      rcarr(94,:) = skarr( 1.250D-12 ,1600.0D+00 ,temperature)
!
!....           HCFC142b + OH = Cl + H2O
!
      rcarr(95,:) = skarr( 1.300D-12 ,1770.0D+00 ,temperature)
!
!....           CH3Cl + Cl = CO + 2 HCl + HO2
!
      rcarr(96,:) = skarr( 3.200D-11 ,1250.0D+00 ,temperature)
!
!....           CH3Br + OH = Br + H2O + HO2
!
!bnd      rcarr(97,:) = skarr( 2.500D-12 ,1310.0D+00 ,temperature)
      rcarr(97,:) = skarr( 2.350D-12 ,1300.0D+00 ,temperature)
!
!....           N2O5 = 2 HNO3
!
      rcarr(98,:) = sksts_n2o5 (temperature ,adcol ,pressure  &
     &           ,sad_lbs ,water)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(99,:) = sksts_clono2 (temperature, adcol ,pressure  &
     &            ,sad_lbs ,specarr(  28,:) ,specarr(    30,:) ,water)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(100,:) = sksts_brono2 (temperature, adcol ,pressure  &
     &             ,sad_lbs ,water)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(101,:) = sksts_clono2_hcl (temperature ,adcol ,pressure  &
     &           ,sad_lbs ,specarr(    28,:) ,specarr(  30,:) ,water)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(102,:) = sksts_hocl_hcl (temperature ,adcol ,pressure  &
     &           ,sad_lbs ,specarr(  29,:) ,specarr( 28,:)  &
     &           ,specarr( 30,:) ,water)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(103,:) = sksts_hobr_hcl (temperature ,adcol ,pressure  &
     &           ,sad_lbs ,specarr(  35,:) ,specarr( 28,:) ,water)
!
!....           N2O5 = 2 HNO3
!
      rcarr(104,:) = sksts_n2o5 (temperature ,adcol ,pressure  &
     &            ,sad_sts ,water)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(105,:) = sksts_clono2 (temperature, adcol ,pressure  &
     &            ,sad_sts ,specarr(  28,:) ,specarr(    30,:) ,water)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(106,:) = sksts_brono2 (temperature, adcol ,pressure  &
     &             ,sad_sts ,water)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(107,:) = sksts_clono2_hcl (temperature ,adcol ,pressure  &
     &           ,sad_sts ,specarr(    28,:) ,specarr(  30,:) ,water)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(108,:) = sksts_hocl_hcl (temperature ,adcol ,pressure  &
     &           ,sad_sts ,specarr(  29,:) ,specarr( 28,:)  &
     &           ,specarr( 30,:) ,water)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(109,:) = sksts_hobr_hcl (temperature ,adcol ,pressure  &
     &           ,sad_sts ,specarr(  35,:) ,specarr( 28,:) ,water)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(110,:) = sknat_clono2 (temperature ,pressure ,sad_nat ,tropp)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(111,:) = sknat_brono2 (temperature ,pressure ,sad_nat ,tropp)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(112,:) = sknat_hcl_clono2 (temperature ,pressure  & 
     &           ,sad_nat ,specarr(  28,:) ,tropp)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(113,:) = sknat_hcl_hocl (temperature ,pressure  & 
     &           ,sad_nat ,specarr(  28,:) ,tropp)
!
!....           BrONO2 + HCl = BrCl + HNO3
!
      rcarr(114,:) = sknat_hcl_brono2 (temperature ,pressure  & 
     &           ,sad_nat ,specarr(  28,:) ,tropp)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(115,:) = sknat_hcl_hobr (temperature ,pressure  & 
     &           ,sad_nat ,specarr(  28,:) ,tropp)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(116,:) = skice_clono2 (temperature ,sad_ice)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(117,:) = skice_brono2 (temperature ,sad_ice)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(118,:) = skice_hcl_clono2 (temperature  &
     &           ,sad_ice ,specarr(  28,:))
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(119,:) = skice_hcl_hocl (temperature  &
     &           ,sad_ice ,specarr(  28,:))
!
!....           BrONO2 + HCl = BrCl + HNO3
!
      rcarr(120,:) = skice_hcl_brono2 (temperature  &
     &           ,sad_ice ,specarr(  28,:))
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(121,:) = skice_hcl_hobr (temperature  &
     &           ,sad_ice ,specarr(  28,:))
!
!....           HNO3 = NO2 + OH
!
      rcarr(122,:) = sksoot_hno3 (temperature ,sad_soot)
!
!....          End thermal rate constants
!
!....          Begin rate constant functions
!
      CONTAINS
        FUNCTION skarr (af,ae,tk)
          real*8  &
     &      af ,ae ,tk(:)
          real*8, dimension(size(tk)) :: skarr
          skarr(:) = af * exp(-ae / tk(:))
        END FUNCTION skarr
        FUNCTION sklp (af,npwr,tk,ad)
          real*8  &
     &      af ,npwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: sklp
          sklp(:) = ad(:) * af * (300.0d0/tk(:))**npwr
        END FUNCTION sklp
        FUNCTION skhp (ai,mpwr,tk)
          real*8  &
     &      ai ,mpwr ,tk(:)
          real*8, dimension(size(tk)) :: skhp
          skhp(:) = ai * (300.0d0/tk(:))**mpwr
        END FUNCTION skhp
        FUNCTION skfo (af,npwr,ai,mpwr,tk,ad)
          real*8  &
     &      af ,npwr ,ai ,mpwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skfo
          skfo(:) = sklp(af,npwr,tk,ad) / skhp(ai,mpwr,tk)
        END FUNCTION skfo
        FUNCTION skterlp (af,npwr,ae,tk,ad)
          real*8  &
     &      af ,npwr ,ae ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skterlp
          skterlp(:) = sklp(af,npwr,tk,ad) * exp(-ae / tk(:))
        END FUNCTION skterlp
        FUNCTION sktroe (af,npwr,ai,mpwr,ae,tk,ad,fc)
          real*8  &
     &      af ,npwr ,ai ,mpwr ,ae ,tk(:) ,ad(:)
          real*8, OPTIONAL :: fc
          real*8, dimension(size(tk)) :: sktroe
          real*8 fsubc
          real*8 skfo_local(size(tk))
          if (present (fc)) then; fsubc=fc; else; fsubc=0.6d0; end if
          skfo_local(:) = skfo(af,npwr,ai,mpwr,tk,ad)
          sktroe(:) = skterlp(af,npwr,ae,tk,ad) * fsubc**  &
     &                (1.0d0/(1.0d0+log10(skfo_local(:))**2)) /  &
     &                (1.0d0+skfo(af,npwr,ai,mpwr,tk,ad))
        END FUNCTION sktroe
!
!.... skho2dis (temperature ,adcol)
!
!.... JPL 97-4
!
        FUNCTION skho2dis (tk,ad)
!
!....       HO2 + HO2 = H2O2 + O2
!....
!.... PSC - 8/2/2002
!....
!.... This routine returns a bimolecular rate constant that accounts
!.... for the pressure, but not the H2O, dependence of the reaction.
!.... The H2O dependence is treated separately.
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skho2dis
          skho2dis(:) = 2.3d-13 * exp(600.0d0 / tk(:)) +  &
     &                  1.7d-33 * ad(:) * exp(1000.0d0 / tk(:))
        END FUNCTION skho2dis
!
!.... skho2h2o (temperature ,adcol)
!
!.... JPL 97-4
!
        FUNCTION skho2h2o (tk,ad)
!
!....       HO2 + HO2 + H2O = H2O2 + O2 + H2O
!....
!.... PSC - 8/2/2002
!....
!.... JPL 97-4 Note B13
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skho2h2o
          skho2h2o(:) = skho2dis(tk ,ad) *  &
     &                  1.4d-21 * exp(2200.0d0 / tk(:))
        END FUNCTION skho2h2o
!
!.... skohhno3 (temperature ,adcol)
!
!.... JPL 00-003 Note C9
!
        FUNCTION skohhno3 (tk,ad)
!
!                    M
!....      OH + HNO3 = NO3 + H2O
!....
!.... PSC - 8/2/2002
!....
!.... Pressure in hPa
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8  &
     &      xxx1(size(tk)) ,xxx2(size(tk))
          real*8, dimension(size(tk)) :: skohhno3
          xxx1(:)     = ad(:) * 6.5d-34 * exp(1335.0d0 / tk(:))
          xxx2(:)     = 2.7d-17 * exp(2199.0d0 / tk(:))
          skohhno3(:) = 2.4d-14 * exp(460.0d0 / tk(:)) +  &
     &                  xxx1(:) / (1.0d0 + xxx1(:) / xxx2(:))
        END FUNCTION skohhno3
!
!.... skbrodis (temperature)
!
!.... JPL 97-4
!
        FUNCTION skbrodis (tk)
!
!....       BrO + BrO = 2 Br + O2
!....
!.... PSC - 8/23/2002
!....
!.... Combined two product channels into one.
!.... Appears as special function to avoid confusion with
!.... actual elementary reaction in reaction database.
!
          real*8  &
     &      tk(:)
          real*8, dimension(size(tk)) :: skbrodis
          skbrodis(:) = 1.5d-12 * exp(230.0d0 / tk(:))
        END FUNCTION skbrodis
!
!.... skcooh (temperature ,pressure)
!
!.... JPL 97-4
!
        FUNCTION skcooh (tk,pr)
!
!                  M
!....      OH + CO = H|HO2 (+ CO2)
!....
!.... PSC - 8/2/2002
!....
!.... Pressure in hPa
!
          real*8  &
     &      tk(:) ,pr(:)
          real*8, dimension(size(tk)) :: skcooh
          skcooh(:) = 1.5d-13 * (1.0d0 + (pr(:)/1013.25d0)*0.6d0)
        END FUNCTION skcooh
!
!.... sksts_n2o5 (temperature ,adcol ,pressure ,sad_lbs ,water, ptrop)
!
!.... JPL 02-25
!
        FUNCTION sksts_n2o5 (tk,ad,pr,sad,h2o,ptrop)
          implicit none
          real*8 :: tk(:), ad(:) ,pr(:) ,sad(:) ,h2o(:)
          real*8, OPTIONAL :: ptrop
          real*8, DIMENSION (size(tk)) :: sksts_n2o5
!... local work variables
          real*8 :: ak0, ak1, ak2,  &
     &      ph2o, p0h2o, aw, y1, y2, am ,wt, gamma, avgvel
          real*8 :: a1, b1, c1, d1, a2, b2, c2, d2, pi
          integer n
!... coefficients for JPL02 N2O5+H2O wt pct dependence
          real*8, save :: ac0(4), ac1(4), ac2(4)
          DATA ac0 /-25.5265, -0.133188, 0.00930846, -9.0194d-5/
          DATA ac1 / 9283.76,   115.345,   -5.19258,  0.0483464/
          DATA ac2 /-851801.,  -22191.2,    766.916,   -6.85427/
!=======================================================================
!     N2O5 + stratospheric sulfate aerosol = 2 HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!... JPL02 form for N2O5+H2O(l)
!
          pi = acos(-1.0d0)
        do n=1,size(tk)
          if( sad(n) > 0.0 ) then
!... partial pressure of h2o, limit
            ph2o = pr(n) * h2o(n)/ad(n)
            if( ph2o < 5.d-11 ) ph2o = 5.d-11
!... Table 1: H2SO4 wt% from T and pH2O
            p0H2O = EXP( 18.452406985 - 3505.1578807/tk(n)  &
     &        - 330918.55082/tk(n)**2 + 12725068.262/tk(n)**3)
            aw = pH2O / p0H2O
!... limit aw to 1
            aw = min(aw,1.0d0)

            if (aw < 0.05) then
              a1 = 12.37208932
              b1 = -0.16125516114
              c1 = -30.490657554
              d1 = -2.1133114241
              a2 = 13.455394705
              b2 = -0.1921312255
              c2 = -34.285174607
              d2 = -1.7620073078
            elseif (aw >= 0.85) then
              a1 = -180.06541028
              b1 = -0.38601102592
              c1 = -93.317846778
              d1 = 273.88132245
              a2 = -176.95814097
              b2 = -0.3625704815410
              c2 = -90.469744201
              d2 = 267.45509988
            else ! (aw >= 0.05)
              a1 = 11.820654354
              b1 = -0.20786404244
              c1 = -4.807306373
              d1 = -5.1727540348
              a2 = 12.891938068
              b2 = -0.23233847708
              c2 = -6.4261237757
              d2 = -4.9005471319
            endif

            y1 = a1*aw**b1+c1*aw+d1
            y2 = a2*aw**b2+c2*aw+d2
            am = y1+(tk(n)-190.)*(y2-y1)/70.
            wt = 9800.*am/(98.*am+1000.)
!... keep bounds wt pct H2SO4 in check
            wt = min(wt,1.0d+2)
            wt = max(wt,1.0d-7)

!... calculate gamma
            ak0 = ac0(1) + ac0(2)*wt + ac0(3)*wt**2 + ac0(4)*wt**3
            ak1 = ac1(1) + ac1(2)*wt + ac1(3)*wt**2 + ac1(4)*wt**3
            ak2 = ac2(1) + ac2(2)*wt + ac2(3)*wt**2 + ac2(4)*wt**3

            gamma = EXP(ak0 + ak1/tk(n) + ak2/(tk(n)**2))

            avgvel = 100.0*sqrt(8.0*8314.48*tk(n)/(pi*mw(IN2O5)))

            sksts_n2o5(n) = 0.25 * gamma * avgvel * sad(n)
!... make sure rate is (+)ve
            sksts_n2o5(n) = max(sksts_n2o5(n),0.0d0)
!... no aerosol
          else
            sksts_n2o5(n) = 0.0
          endif

          if ( present(ptrop) ) then
            if( pr(n) > ptrop )  sksts_n2o5(n) = 0.0
          endif
         enddo
        END FUNCTION sksts_n2o5
!
!.... sksts_clono2 (temperature ,adcol ,pressure ,sad_lbs ,specarr( HCl,:),specarr( ClONO2,:) ,water, ptrop)
!
!.... JPL 02-25
!
        FUNCTION sksts_clono2 (tk,ad,pr,sad,hcl,cln,h2o,ptrop)
          implicit none
!... passed variables
          real*8 :: tk(:) ,ad(:) ,pr(:), sad(:), h2o(:) ,hcl(:) ,cln(:)
          real*8, OPTIONAL :: ptrop
!... function definition
          real*8, DIMENSION (size(tk)) :: sksts_clono2
!... local work variables
          real*8 :: pi

!          real*8, DIMENSION (size(tk)) :: ph2o, p0H2O, aw, a1, b1, c1,
          real*8 :: ph2o, p0H2O, aw, a1, b1, c1,  &
     &     d1, a2, b2, c2, d2, y1, y2, am, wt, pHCl, pClONO2
           integer n
!          real*8, DIMENSION (size(tk)) :: z1, z2, z3, rho, aMSO4, X, A,
          real*8 :: z1, z2, z3, rho, aMSO4, X, A,  &
     &     To, eval, h, aH, c_ClONO2, SClONO2, HClONO2, DClONO2, akH2O,  &
     &     akH, akhydr, GbH2O, HHCl, aMHCl, akHCl, alClONO2, fClONO2,  &
     &     GClONO2rxn, GbHCl, Gs, FHCl, Gsp, GbHClp, Gb, gClONO2,  &
     &     gClONO2_HCl, gClONO2_H2O, avgvel
!... set aerosol radius
          real*8 :: rad
!... number of particles/cm**3
          real*8, save :: sparts
          DATA sparts /10./

!=======================================================================
!     ClONO2 + stratospheric sulfate aerosol = HOCl + HNO3
!=======================================================================
!
!
! Stratospheric Reaction Coefficients of ClONO2 and HOCl
! Shi, Q., P. Davidovits, Boston College, Chestnut Hill, MA 02167
! D. R. Worsnop, T. J. Jayne, C. E. Kolb, Aerodyne Research, Inc., Billerica,
! MA 01821
! worsnop@aerodyne.com//shiq@aerodyne.com
! INSTRUCTIONS FOR USE:
!
! *INPUT*
!       tk  -- The user must provide Temperature (K, scalar or vector)
! --scalars or vector of same length as tk --
!       pr  -- Pressure (mb)
!       h2o -- water vapor mixing ratio
!       hcl -- HCl mixing ratio
!       cln -- ClONO2 mixing ratio
!       rad -- Aerosol radius (cm), typical value is 1e-5
!
! *OUTPUT*
! Heterogeneous reaction probabilities for:
! ClONO2 reacting with HCl and H2O
! HOCl reacting with HCl.
! These are symbolized in the following code as:
! gClONO2_HCl, gClONO2_H2O, gHOCl
!
! Convert to Water vapor partial press (mbar)
! HCl and ClONO2 partial pressure (atm).
! These are symbolized in the following code, respectively, as:
! tk, r, PH2O, PHCl, PClONO2

        pi = acos(-1.0d0)
        do n=1,size(tk)
          if( sad(n) > 0.0 ) then
!... calculate radius of aerosol particles
            rad = sqrt(sad(n)/(4.*pi*sparts))
!... partial pressure of h2o, limit
            ph2o = pr(n) * h2o(n)/ad(n)
            if( ph2o < 5.d-11 ) ph2o = 5.d-11
!... Table 1: H2SO4 wt% from T and pH2O
            p0H2O = EXP( 18.452406985 - 3505.1578807/tk(n)  &
     &        - 330918.55082/tk(n)**2 + 12725068.262/tk(n)**3)
            aw = pH2O / p0H2O
!... limit aw to 1
            aw = min(aw,1.0d0)

            if (aw < 0.05) then
              a1 = 12.37208932
              b1 = -0.16125516114
              c1 = -30.490657554
              d1 = -2.1133114241
              a2 = 13.455394705
              b2 = -0.1921312255
              c2 = -34.285174607
              d2 = -1.7620073078
            elseif (aw >= 0.85) then
              a1 = -180.06541028
              b1 = -0.38601102592
              c1 = -93.317846778
              d1 = 273.88132245
              a2 = -176.95814097
              b2 = -0.3625704815410
              c2 = -90.469744201
              d2 = 267.45509988
            else !where (aw >= 0.05)
              a1 = 11.820654354
              b1 = -0.20786404244
              c1 = -4.807306373
              d1 = -5.1727540348
              a2 = 12.891938068
              b2 = -0.23233847708
              c2 = -6.4261237757
              d2 = -4.9005471319
            endif

            y1 = a1*aw**b1+c1*aw+d1
            y2 = a2*aw**b2+c2*aw+d2
            am = y1+(tk(n)-190.)*(y2-y1)/70.
            wt = 9800.*am/(98.*am+1000.)
!... keep bounds wt pct H2SO4 in check
            wt = min(wt,1.0d+2)
            wt = max(wt,1.0d-7)

            phcl = pr(n)*hcl(n)/ad(n)/1013.25
            pclono2 = pr(n)*cln(n)/ad(n)/1013.25

!... Table 2: Parameters for H2SO4 Solution
            z1 = 0.12364-5.6e-7*tk(n)**2
            z2 = -0.02954+1.814e-7*tk(n)**2
            z3 = 2.343e-3-1.487e-6*tk(n)-1.324e-8*tk(n)**2
            rho = 1+z1*am+z2*am**1.5+z3*am**2
            aMSO4 = rho*wt/9.8
            X = wt/(wt+(100.-wt)*98./18.)
            A = 169.5+5.18*wt-0.0825*wt**2.+3.27e-3*wt**3.
            To = 144.11+0.166*wt-0.015*wt**2.+2.18e-4*wt**3.
!... This has caused grief because T gets too close to To
            eval = max (tk(n)-To,1.0)

            h = A*tk(n)**(-1.43)*exp(448/(eval))
            aH = exp(60.51-0.095*wt+0.0077*wt**2.-1.61e-5*wt**3.  &
     &        -(1.76+2.52e-4*wt**2)*sqrt(tk(n))  &
     &        +(-805.89+253.05*wt**0.076)/sqrt(tk(n)))

!... Table 3. ClONO2 + H2O and ClONO2+HCl
            c_ClONO2 = 1474*sqrt(tk(n))
            SClONO2 = 0.306+24.0/tk(n)
            HClONO2 = 1.6e-6*exp(4710./tk(n))*exp(-SClONO2*aMSO4)
            DClONO2 = 5e-8*tk(n)/h
            akH2O = 1.95e10*exp(-2800./tk(n))
            akH = 1.22e12*exp(-6200./tk(n))
            akhydr = akH2O*aw + akH*aH*aw
            GbH2O = 4.*HClONO2*0.082*tk(n)  &
     &       *sqrt(DClONO2*akhydr)/c_ClONO2
            HHCl = (0.094-0.61*X+1.2*X**2.)  &
     &        *exp(-8.68+(8515-10718*X**0.7)/tk(n))
            aMHCl = HHCl *PHCl
            akHCl = 7.9e11*aH*DClONO2*amhcl
            alClONO2 = sqrt(DClONO2/(akhydr+akHCl))
            fClONO2 = 1/tanh(rad/alClONO2)- alClONO2/rad
            GClONO2rxn = fClONO2*GbH2O *sqrt(1+akHCl/akhydr)
            GbHCl = GClONO2rxn* akHCl/( akHCl+ akhydr)
            Gs = 66.12*exp(-1374/tk(n))*HClONO2*amhcl
            FHCl = 1
            if(phcl > 0.0) then
               FHCl = 1/(1+0.612*(Gs+GbHCl)* pClONO2/pHCl)
             endif

            Gsp = FHCl*Gs
            GbHClp = FHCl*GbHCl
            Gb = GbHClp + GClONO2rxn* akhydr/( akHCl+ akhydr)
            gClONO2 = 0.0
            gClONO2_HCl = 0.0
            if((Gsp + Gb) > 0.0) then
              gClONO2 = 1./(1.+1./(Gsp + Gb))
              gClONO2_HCl = gClONO2 *(Gsp + GbHClp)/(Gsp + Gb)
             endif
            gClONO2_H2O = 0.0
            if(gClONO2_HCl < gClONO2) then
              gClONO2_H2O = gClONO2 - gClONO2_HCl
             endif

            avgvel = 100.0*sqrt(8.0*8314.48*tk(n)/(pi*mw(ICLONO2)))

            sksts_clono2(n) = 0.25 * gClONO2_H2O * avgvel * sad(n)
!... make sure rate is (+)ve
            sksts_clono2(n) = max(sksts_clono2(n),0.0d0)
!... no aerosol
          else
            sksts_clono2(n) = 0.0
          endif

          if ( present(ptrop) ) then
            if( pr(n) > ptrop )  sksts_clono2(n) = 0.0
          end if
         enddo
        END FUNCTION sksts_clono2
!
!.... sksts_brono2 (temperature ,sad_lbs)
!
!.... JPL 02-25
!
        FUNCTION sksts_brono2 (tk,ad,pr,sad,h2o,ptrop)
          implicit none
!... passed variables
          real*8 :: tk(:), ad(:), pr(:), sad(:), h2o(:)
          real*8, OPTIONAL :: ptrop
!... function definition
          real*8, DIMENSION (size(tk)) :: sksts_brono2
!... local work variables
          real*8 :: pi
          integer :: n
          real*8 :: ph2o, p0H2O, aw, a1, b1, c1,  &
     &     d1, a2, b2, c2, d2, y1, y2, am, wt, gamma, avgvel


!=======================================================================
!     BrONO2 + stratospheric sulfate aerosol = HOBr + HNO3
!=======================================================================
!

!... JPL02-25 for BrONO2+H2O(l)
          pi = acos(-1.0d0)
        do n=1,size(tk)
          if( sad(n) > 0.0 ) then
!... partial pressure of h2o, limit
            ph2o = pr(n) * h2o(n)/ad(n)
            if( ph2o < 5.d-11 ) ph2o = 5.d-11
!... Table 1: H2SO4 wt% from T and pH2O
            p0H2O = EXP( 18.452406985 - 3505.1578807/tk(n)  &
     &        - 330918.55082/tk(n)**2 + 12725068.262/tk(n)**3)
            aw = pH2O / p0H2O
!... limit aw to 1
            aw = min(aw,1.0d0)

            if (aw < 0.05) then
              a1 = 12.37208932
              b1 = -0.16125516114
              c1 = -30.490657554
              d1 = -2.1133114241
              a2 = 13.455394705
              b2 = -0.1921312255
              c2 = -34.285174607
              d2 = -1.7620073078
            elseif (aw >= 0.85) then
              a1 = -180.06541028
              b1 = -0.38601102592
              c1 = -93.317846778
              d1 = 273.88132245
              a2 = -176.95814097
              b2 = -0.3625704815410
              c2 = -90.469744201
              d2 = 267.45509988
            else  ! where (aw >= 0.05)
              a1 = 11.820654354
              b1 = -0.20786404244
              c1 = -4.807306373
              d1 = -5.1727540348
              a2 = 12.891938068
              b2 = -0.23233847708
              c2 = -6.4261237757
              d2 = -4.9005471319
            endif

            y1 = a1*aw**b1+c1*aw+d1
            y2 = a2*aw**b2+c2*aw+d2
            am = y1+(tk(n)-190.)*(y2-y1)/70.
            wt = 9800.*am/(98.*am+1000.)
!... keep bounds wt pct H2SO4 in check
            wt = min(wt,1.0d+2)
            wt = max(wt,1.0d-7)

            gamma = 1./(1./.805+1./(EXP(29.24-.396*wt)+0.114))

            avgvel = 100.0*sqrt(8.0*8314.48*tk(n)/(pi*mw(IBRONO2)))

            sksts_brono2(n) = 0.25 * gamma * avgvel * sad(n)
!... make sure rate is (+)ve
            sksts_brono2(n) = max(sksts_brono2(n),0.0d0)
!... no aerosol
           else
            sksts_brono2(n) = 0.0
           endif

          if ( present(ptrop) ) then
            if( pr(n) > ptrop )  sksts_brono2(n) = 0.0
          end if
         enddo

        END FUNCTION sksts_brono2
!
!.... sksts_clono2_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HCl,:) ,specarr(ClONO2,:), ,water
!
!.... JPL 02-25
!
        FUNCTION sksts_clono2_hcl (tk,ad,pr,sad,hcl,cln,h2o,ptrop)
          implicit none
!... passed variables
          real*8 :: tk(:) ,ad(:) ,pr(:), sad(:), h2o(:) ,hcl(:) ,cln(:)
          real*8, OPTIONAL :: ptrop
!... function definition
          real*8, DIMENSION (size(tk)) :: sksts_clono2_hcl
!... local work variables
          real*8 :: hcle ,pi
          integer :: n
          real*8 :: ph2o, p0H2O, aw, a1, b1, c1,  &
     &     d1, a2, b2, c2, d2, y1, y2, am, wt, pHCl, pClONO2
          real*8 :: z1, z2, z3, rho, aMSO4, X, A,  &
     &     To, eval, h, aH, c_ClONO2, SClONO2, HClONO2, DClONO2, akH2O,  &
     &     akH, akhydr, GbH2O, HHCl, aMHCl, akHCl, alClONO2, fClONO2,  &
     &     GClONO2rxn, GbHCl, Gs, FHCl, Gsp, GbHClp, Gb, gClONO2,  &
     &     gClONO2_HCl, gClONO2_H2O, avgvel
!... set aerosol radius
          real*8 :: rad
!... number of particles/cm**3
          real*8, save :: sparts
          DATA sparts /10./


!=======================================================================
!     ClONO2 + HCl on stratospheric sulfate aerosol = Cl2 + HNO3
!=======================================================================
!
!
! Stratospheric Reaction Coefficients of ClONO2 and HOCl
! Shi, Q., P. Davidovits, Boston College, Chestnut Hill, MA 02167
! D. R. Worsnop, T. J. Jayne, C. E. Kolb, Aerodyne Research, Inc., Billerica,
! MA 01821
! worsnop@aerodyne.com//shiq@aerodyne.com
! INSTRUCTIONS FOR USE:
!
! *INPUT*
!       tk  -- The user must provide Temperature (K, scalar or vector)
! --scalars or vector of same length as tk --
!       pr  -- Pressure (mb)
!       h2o -- water vapor mixing ratio
!       hcl -- HCl mixing ratio
!       cln -- ClONO2 mixing ratio
!       rad -- Aerosol radius (cm), typical value is 1e-5
!
! *OUTPUT*
! Heterogeneous reaction probabilities for:
! ClONO2 reacting with HCl and H2O
! HOCl reacting with HCl.
! These are symbolized in the following code as:
! gClONO2_HCl, gClONO2_H2O, gHOCl
!
! Convert to Water vapor partial press (mbar)
! HCl and ClONO2 partial pressure (atm).
! These are symbolized in the following code, respectively, as:
! tk, r, PH2O, PHCl, PClONO2

        pi = acos(-1.0d0)
        do n=1,size(tk)
          if( sad(n) > 0.0 ) then
!
!... calculate mean aerosol radius
            rad = sqrt(sad(n)/(4.*pi*sparts))
!... partial pressure of h2o, limit
            ph2o = pr(n) * h2o(n)/ad(n)
            if( ph2o < 5.d-11 ) ph2o = 5.d-11
!... Table 1: H2SO4 wt% from T and pH2O
            p0H2O = EXP( 18.452406985 - 3505.1578807/tk(n)  &
     &        - 330918.55082/tk(n)**2 + 12725068.262/tk(n)**3)
            aw = pH2O / p0H2O
!... limit aw to 1
            aw = min(aw,1.0d0)

            if (aw < 0.05) then
              a1 = 12.37208932
              b1 = -0.16125516114
              c1 = -30.490657554
              d1 = -2.1133114241
              a2 = 13.455394705
              b2 = -0.1921312255
              c2 = -34.285174607
              d2 = -1.7620073078
            elseif (aw >= 0.85) then
              a1 = -180.06541028
              b1 = -0.38601102592
              c1 = -93.317846778
              d1 = 273.88132245
              a2 = -176.95814097
              b2 = -0.3625704815410
              c2 = -90.469744201
              d2 = 267.45509988
            else  !where (aw >= 0.05)
              a1 = 11.820654354
              b1 = -0.20786404244
              c1 = -4.807306373
              d1 = -5.1727540348
              a2 = 12.891938068
              b2 = -0.23233847708
              c2 = -6.4261237757
              d2 = -4.9005471319
            endif

            y1 = a1*aw**b1+c1*aw+d1
            y2 = a2*aw**b2+c2*aw+d2
            am = y1+(tk(n)-190.)*(y2-y1)/70.
            wt = 9800.*am/(98.*am+1000.)
!... keep bounds wt pct H2SO4 in check
            wt = min(wt,1.0d+2)
            wt = max(wt,1.0d-7)

            phcl = pr(n)*hcl(n)/ad(n)/1013.25
            pclono2 = pr(n)*cln(n)/ad(n)/1013.25

!... Table 2: Parameters for H2SO4 Solution
            z1 = 0.12364-5.6e-7*tk(n)**2
            z2 = -0.02954+1.814e-7*tk(n)**2
            z3 = 2.343e-3-1.487e-6*tk(n)-1.324e-8*tk(n)**2
            rho = 1+z1*am+z2*am**1.5+z3*am**2
            aMSO4 = rho*wt/9.8
            X = wt/(wt+(100.-wt)*98./18.)
            A = 169.5+5.18*wt-0.0825*wt**2.+3.27e-3*wt**3.
            To = 144.11+0.166*wt-0.015*wt**2.+2.18e-4*wt**3.

!... This has caused grief because T gets too close to To
            eval = max (tk(n)-To,1.0)

            h = A*tk(n)**(-1.43)*exp(448./(eval))
            aH = exp(60.51-0.095*wt+0.0077*wt**2.  &
     &        -1.61e-5*wt**3.-(1.76+2.52e-4*wt**2)*sqrt(tk(n))  &
     &        +(-805.89+253.05*wt**0.076)/sqrt(tk(n)))

!... Table 3. ClONO2 + H2O and ClONO2+HCl
            c_ClONO2 = 1474*sqrt(tk(n))
            SClONO2 = 0.306+24.0/tk(n)
            HClONO2 = 1.6e-6*exp(4710./tk(n))*exp(-SClONO2*aMSO4)
            DClONO2 = 5e-8*tk(n)/h
            akH2O = 1.95e10*exp(-2800./tk(n))
            akH = 1.22e12*exp(-6200./tk(n))
            akhydr = akH2O*aw + akH*aH*aw
            GbH2O = 4.*HClONO2*0.082*tk(n)  &
     &       *sqrt(DClONO2*akhydr)/c_ClONO2
            HHCl = (0.094-0.61*X+1.2*X**2.)  &
     &       *exp(-8.68+(8515-10718*X**0.7)/tk(n))
            aMHCl = HHCl *PHCl
            akHCl = 7.9e11*aH*DClONO2*amhcl
            alClONO2 = sqrt(DClONO2/(akhydr+akHCl))
            fClONO2 = 1/tanh(rad/alClONO2)- alClONO2/rad
            GClONO2rxn = fClONO2*GbH2O*sqrt(1+akHCl/akhydr)
            GbHCl = GClONO2rxn* akHCl/( akHCl+ akhydr)
            Gs = 66.12*exp(-1374/tk(n))*HClONO2*amhcl
            FHCl = 1.
            if(phcl > 0.0)  &
     &         FHCl = 1/(1+0.612*(Gs+GbHCl)* pClONO2/pHCl)

            Gsp = FHCl*Gs
            GbHClp = FHCl*GbHCl
            Gb = GbHClp+GClONO2rxn*akhydr/(akHCl+akhydr)
            gClONO2 = 0.0
            gClONO2_HCl = 0.0
            if((Gsp + Gb) > 0.0) then
              gClONO2 = 1./(1.+1./(Gsp + Gb))
              gClONO2_HCl = gClONO2*(Gsp+GbHClp)/(Gsp+Gb)
             endif

            avgvel = 100.0*sqrt(8.0*8314.48*tk(n)/(pi*mw(ICLONO2)))

            hcle = max(hcl(n),1.0d6*sad(n))
            sksts_clono2_hcl(n) = 0.25*gClONO2_HCl*avgvel*sad(n)/hcle

!... make sure rate is (+)ve
            sksts_clono2_hcl(n) = max(sksts_clono2_hcl(n),0.0d0)
!... no aerosol, no reaction
           else
            sksts_clono2_hcl(n) = 0.0
           endif

          if ( present(ptrop) ) then
            if( pr(n) > ptrop )  sksts_clono2_hcl(n) = 0.0
           endif
         enddo
        END FUNCTION sksts_clono2_hcl
!
!.... sksts_hocl_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOCl,:) ,
!
!.... JPL 02-25
!
        FUNCTION sksts_hocl_hcl (tk,ad,pr,sad,hocl,hcl,cln,h2o,ptrop)
          implicit none

!... passed variables
          real*8 :: tk(:) ,ad(:) , pr(:), sad(:), h2o(:), hocl(:),  &
     &     hcl(:), cln(:)
          real*8, OPTIONAL :: ptrop
!... function definition
          real*8, DIMENSION (size(tk)) :: sksts_hocl_hcl
!... local work variables
          real*8 :: hcle ,pi
          integer :: n
          real*8 :: ph2o, p0H2O, aw, a1, b1, c1,  &
     &     d1, a2, b2, c2, d2, y1, y2, am, wt, pHCl, pClONO2

          real*8 :: z1, z2, z3, rho, aMSO4, X, A,  &
     &     To, eval, h, aH, c_ClONO2, SClONO2, HClONO2, DClONO2, akH2O,  &
     &     akH, akhydr, GbH2O, HHCl, aMHCl, akHCl, alClONO2, fClONO2,  &
     &     GClONO2rxn, GbHCl, Gs, FHCl, Gsp, gHOCl, c_HOCl, SHOCl,  &
     &     HHOCl, DHOCl, akHOCl_HCl, GHOClrxn, alHOCl, fHOCl, avgvel
!... set aerosol radius
          real*8 :: rad
!... number of particles/cm**3
          real*8, save :: sparts
          DATA sparts /10./


!=======================================================================
!     HOCl + HCl on stratospheric sulfate aerosol = Cl2 + H2O
!=======================================================================
!
        pi = acos(-1.0d0)
        do n=1,size(tk)
          if( sad(n) > 0.0 ) then
!... calculate mean radius of aerosol particles
            rad = sqrt(sad(n)/(4.*pi*sparts))
!... partial pressure of h2o, limit
            ph2o = pr(n) * h2o(n)/ad(n)
            if( ph2o < 5.d-11 ) ph2o = 5.d-11
!... Table 1: H2SO4 wt% from T and pH2O
            p0H2O = EXP( 18.452406985 - 3505.1578807/tk(n)  &
     &        - 330918.55082/tk(n)**2 + 12725068.262/tk(n)**3)
            aw = pH2O / p0H2O
!... limit aw to 1
            aw = min(aw,1.0d0)

            if (aw < 0.05) then
              a1 = 12.37208932
              b1 = -0.16125516114
              c1 = -30.490657554
              d1 = -2.1133114241
              a2 = 13.455394705
              b2 = -0.1921312255
              c2 = -34.285174607
              d2 = -1.7620073078
            elseif (aw >= 0.85) then
              a1 = -180.06541028
              b1 = -0.38601102592
              c1 = -93.317846778
              d1 = 273.88132245
              a2 = -176.95814097
              b2 = -0.3625704815410
              c2 = -90.469744201
              d2 = 267.45509988
            else  !where (aw >= 0.05)
              a1 = 11.820654354
              b1 = -0.20786404244
              c1 = -4.807306373
              d1 = -5.1727540348
              a2 = 12.891938068
              b2 = -0.23233847708
              c2 = -6.4261237757
              d2 = -4.9005471319
            endif

            y1 = a1*aw**b1+c1*aw+d1
            y2 = a2*aw**b2+c2*aw+d2
            am = y1+(tk(n)-190.)*(y2-y1)/70.
            wt = 9800.*am/(98.*am+1000.)
!... keep bounds wt pct H2SO4 in check
            wt = min(wt,1.0d+2)
            wt = max(wt,1.0d-7)

            pHCl = pr(n)*hcl(n)/ad(n)/1013.25
            pClONO2 = pr(n)*cln(n)/ad(n)/1013.25

!... Table 2: Parameters for H2SO4 Solution
            z1 = 0.12364 -5.6d-7*tk(n)**2
            z2 = -0.02954 +1.814d-7*tk(n)**2
            z3 = 2.343d-3 -1.487d-6*tk(n) -1.324d-8*tk(n)**2
            rho = 1. +z1*am +z2*am**1.5 +z3*am**2
            aMSO4 = rho*wt/9.8
            X = wt/(wt+(100.-wt)*98./18.)
            A = 169.5 +5.18*wt-0.0825*wt**2.+3.27d-3*wt**3.
            To = 144.11 +0.166*wt -0.015*wt**2. +2.18d-4*wt**3.
!... This has caused grief because T gets too close to To
            eval = tk(n)-To
            if(eval < 1.0)  eval = 1.0

            h = A*tk(n)**(-1.43)*exp(448./(eval))
            aH = exp(60.51-0.095*wt+0.0077*wt**2.-1.61d-5*wt**3.  &
     &       -(1.76+2.52d-4*wt**2)*sqrt(tk(n))  &
     &       +(-805.89 +253.05*wt**0.076)/sqrt(tk(n)))
!... Table 3. ClONO2 + H2O and ClONO2+HCl
            c_ClONO2 = 1474.*sqrt(tk(n))
            SClONO2 = 0.306+24.0/tk(n)
            HClONO2 = 1.6d-6*exp(4710./tk(n))*exp(-SClONO2*aMSO4)
            DClONO2 = 5.d-8*tk(n)/h
            akH2O = 1.95d10*exp(-2800./tk(n))
            akH = 1.22d12*exp(-6200./tk(n))
            akhydr = akH2O*aw + akH*aH*aw
            GbH2O = 4.*HClONO2*0.082*tk(n)*sqrt(DClONO2*akhydr)/c_ClONO2
            HHCl = (0.094-0.61*X+1.2*X**2.)  &
     &        *exp(-8.68+(8515-10718*X**0.7)/tk(n))
            aMHCl = HHCl *pHCl
            akHCl = 7.9d11*aH*DClONO2*amhcl
            alClONO2 = sqrt(DClONO2/(akhydr+akHCl))
            fClONO2 = 1./tanh(rad/alClONO2)- alClONO2/rad
            GClONO2rxn = fClONO2*GbH2O *sqrt(1.+akHCl/akhydr)
            GbHCl = GClONO2rxn* akHCl/( akHCl+ akhydr)
            Gs = 66.12*exp(-1374./tk(n))*HClONO2*aMHCl
            FHCl = 1.
            gHOCl = 0.0
            if(phcl > 0.0) then
              FHCl = 1./(1.+0.612*(Gs+GbHCl)* pClONO2/pHCl)

!... Table 4. HOCl + HCl
              c_HOCl = 2009.*sqrt(tk(n))
              SHOCl = 0.0776+59.18/tk(n)
              HHOCl = 1.91d-6*exp(5862.4/tk(n))*exp(-SHOCl*aMSO4)
              DHOCl = 6.4d-8*tk(n)/h
              akHOCl_HCl = 1.25d9*aH*DHOCl*amhcl
!              GHOClrxn = 4*HHOCl*0.082*tk(n)*(DHOCl*akHOCl_HCl)**0.5/c_HOCl
              GHOClrxn = 4.*HHOCl*0.082*tk(n)
              GHOClrxn = GHOClrxn*sqrt(DHOCl*akHOCl_HCl)
              GHOClrxn = GHOClrxn/c_HOCl
              alHOCl = sqrt(DHOCl/akHOCl_HCl)
              fHOCl = 1./tanh(rad/alHOCl)- alHOCl/rad
              gHOCl = 0.0
              eval = fHOCl*GHOClrxn*FHCl
              if (eval > 0.0) gHOCl = 1./( 1.+1./eval )
             endif
!... limit gamma to 1.0
            gHOCl = min(gHOCl,1.0d0)

            avgvel = 100.0*sqrt(8.0*8314.48*tk(n)/(pi*mw(IHOCL)))

            hcle = max(hcl(n),1.0d6*sad(n))
            sksts_hocl_hcl(n) = 0.25 * gHOCl * avgvel * sad(n) / hcle
!... make sure rate is (+)ve
            sksts_hocl_hcl(n) = max(sksts_hocl_hcl(n),0.0d0)
!... no aerosol
           else
            sksts_hocl_hcl(n) = 0.0
           endif

          if ( present(ptrop) ) then
            if( pr(n) > ptrop )  sksts_hocl_hcl(n) = 0.0
           endif
         enddo
        END FUNCTION sksts_hocl_hcl
!
!.... sksts_hobr_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOBr,:) ,
!
!_6_
!
!.... JPL 02-25
!
        FUNCTION sksts_hobr_hcl (tk,ad,pr,sad,hobr,hcl,h2o,ptrop)
          implicit none
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,ad(:) ,pr(:) ,sad(:) ,hobr(:) ,h2o(:) ,hcl(:)
          real*8, DIMENSION (size(tk)) :: sksts_hobr_hcl
!.          real*8
!.     &      adrop ,alpha ,d1 ,hsqrtd ,kii ,minconc ,pi
!.          real*8
!.     &      adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))
!.     &     ,fterm(size(tk))
!.     &     ,gcalc(size(tk))
!.     &     ,gprob_tot(size(tk))
!.     &     ,hstar(size(tk))
!.     &     ,k(size(tk))
!.     &     ,ph2o(size(tk)) ,phcl(size(tk))
!.     &     ,tk_150(size(tk))

!=======================================================================
!     HOBr + HCl on stratospheric sulfate aerosol = BrCl + H2O
!=======================================================================

!.sds taking out below parameterization for JPL02 update
!.          pi            = acos(-1.0d0)
!.c
!.c.... First order reaction rate constant
!.c
!.c.... Hanson and Ravi, JPC, 98, 5728, 1994
!.c.... DEK, 1/10/97
!.c
!.c....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!.c....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!.c....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!.c
!.          minconc       = 1.0d0
!.          adrop         = 1.0D-05
!.          alpha         = 1.0d0
!.          d1            = 1.2D-08
!.c
!.c....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!.c
!.          ph2o(:)       = (h2o(:) / ad(:)) * pr(:)
!.
!.          ph2o(:)       = ph2o(:) / 1013.25d0
!.
!.          where( hcl(:) <= minconc )
!.            phcl        = ((1.0d0 / ad(:)) * pr(:)) *
!.     &                    (1.0d0 / 1013.25d0)
!.          elsewhere
!.            phcl        = ((hcl(:) / ad(:)) * pr(:)) *
!.     &                    (1.0d0 / 1013.25d0)
!.          end where
!.c
!.c....    NOTE: Activity of H2O not allowed to exceed 1.1
!.c....          ah2o and Hstar taken from Table 2, JPC,
!.c....          Hanson and Ravi, 98, 5734, 1994
!.c
!.          ah2o(:)       = 1013.25d0 * ph2o(:) /
!.     &                    10.0d0**
!.     &                    (9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!.
!.          where( ah2o(:) > 1.1d0)
!.            ah2o        = 1.1d0
!.          end where
!.
!.          tk_150(:)     = tk(:)
!.
!.          where( tk(:) < 150.0d0 )
!.            tk_150      = 150.0d0
!.          endwhere
!.
!.          hstar(:)      = exp((6250.0d0 / tk_150(:)) - 10.414d0) *
!.     &                    (ah2o(:)**3.49d0)
!.
!.          kii           = 1.0D+05
!.
!.          k(:)          = kii * hstar(:) * phcl(:)
!.
!.          hsqrtd        = 110.0d0
!.
!.          gcalc(:)      = 2.2548D-05 * hsqrtd *
!.     &                    sqrt(tk(:) * mw(IHOBR) * k(:))
!.
!.          adivl(:)      = adrop / sqrt(d1 / k(:))
!.
!.          fterm(:)      = ((exp(adivl(:)) + exp(-adivl(:))) /
!.     &                     (exp(adivl(:)) - exp(-adivl(:)))) -
!.     &                    (1.0d0 / adivl(:))
!.c
!.c....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!.c
!.          where( fterm > 0.0d0 )
!.            gprob_tot   = 1.0d0 / (1.0d0 /
!.     &                    (fterm(:) * gcalc(:)) +
!.     &                    1.0d0 / alpha)
!.          elsewhere
!.            gprob_tot   = 0.0d0
!.          end where
!.
!.          avgvel(:)     = 100.0d0 *
!.     &                    (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /
!.     &                    (pi * mw(IHOBR)))**0.5d0
!.
!.          where( hcl > minconc )
!.            sksts_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad / hcl
!.          elsewhere
!.            sksts_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad
!.          end where
!.
!.          where( sad < 0.0d0 )
!.     &      sksts_hobr_hcl = 0.0d0
!.
!.
!.          if ( present(ptrop) ) then
!.            where( pr > ptrop )
!.     &        sksts_hobr_hcl = 0.0d0
!.          endif
!.sds... end of stuff taken out

!... JPL02 has lots of caveats and uncertainity - ignore for now
          sksts_hobr_hcl = 0.0d0

        END FUNCTION sksts_hobr_hcl
!
!.... sknat_clono2 (temperature ,sad_nat,ptrop)
!
!.... (1) JPL 00-003
!
        FUNCTION sknat_clono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: sknat_clono2
          real*8  &
     &      gprob ,pi
          real*8  &
     &      avgvel(size(tk))
!
          pi              = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + PSC Type I NAT particles = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob           = 0.004d0
          avgvel(:)       = 100.0d0 *  &
     &                      (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(ICLONO2)))**0.5d0
!
          sknat_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
          if ( present(ptrop) ) then
            where( pr > ptrop )  &
     &        sknat_clono2 = 0.0d0
          end if
        END FUNCTION sknat_clono2
!
!.... sknat_brono2 (temperature ,sad_nat)
!
!.... (2) JPL 00-003
!
        FUNCTION sknat_brono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: sknat_brono2
          real*8  &
     &      gprob ,pi
          real*8  &
     &      avgvel(size(tk))
!
          pi              = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + PSC Type I NAT particles = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob           = 0.004d0
          avgvel(:)       = 100.0d0 *  &
     &                      (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(IBRONO2)))**0.5d0
!
          sknat_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
          if ( present(ptrop) ) then
            where( pr > ptrop )  &
     &        sknat_brono2 = 0.0d0
          end if
        END FUNCTION sknat_brono2
!
!.... sknat_hcl_clono2 (temperature ,sad_nat ,specarr( HCl,:))
!
!.... (3) JPL 00-003
!
        FUNCTION sknat_hcl_clono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_clono2
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))
!
          pi                  = acos(-1.0d0)
!
!=======================================================================
!     HCl + ClONO2 on PSC Type I NAT particles = Cl2 + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc             = 1.0d0
          gprob               = 0.20d0
          avgvel(:)           = 100.0d0 *  &
     &                          (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                          (pi * mw(ICLONO2)))**0.5d0
!
          sknat_hcl_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_clono2  = sknat_hcl_clono2 / minconc
          elsewhere
            sknat_hcl_clono2  = sknat_hcl_clono2 / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop )  &
     &        sknat_hcl_clono2 = 0.0d0
          end if
        END FUNCTION sknat_hcl_clono2
!
!.... sknat_hcl_hocl (temperature ,sad_nat ,specarr( HCl,:))
!
!.... (4) JPL 00-003
!
        FUNCTION sknat_hcl_hocl (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_hocl
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))
!
          pi                  = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOCl on PSC Type I NAT particles = Cl2 + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc           = 1.0d0
          gprob             = 0.10d0
          avgvel(:)         = 100.0d0 *  &
     &                        (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                        (pi * mw(IHOCL)))**0.5d0
!
          sknat_hcl_hocl(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_hocl  = sknat_hcl_hocl / minconc
          elsewhere
            sknat_hcl_hocl  = sknat_hcl_hocl / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop )  &
     &        sknat_hcl_hocl = 0.0d0
          end if
        END FUNCTION sknat_hcl_hocl
!
!.... sknat_hcl_brono2 (temperature ,sad_nat ,specarr( HCl,:))
!
!.... (5) JPL 00-003
!
        FUNCTION sknat_hcl_brono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_brono2
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))
!
          pi                  = acos(-1.0d0)
!
!=======================================================================
!     HCl + BrONO2 on PSC Type I NAT particles = BrCl + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc             = 1.0d0
          gprob               = 0.20d0
          avgvel(:)           = 100.0d0 *  &
     &                          (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                          (pi * mw(IBRONO2)))**0.5d0
!
          sknat_hcl_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_brono2  = sknat_hcl_brono2 / minconc
          elsewhere
            sknat_hcl_brono2  = sknat_hcl_brono2 / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop )  &
     &        sknat_hcl_brono2 = 0.0d0
          end if
        END FUNCTION sknat_hcl_brono2
!
!.... sknat_hcl_hobr (temperature ,sad_nat ,specarr( HCl,:))
!
!.... (6) JPL 00-003
!
        FUNCTION sknat_hcl_hobr (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_hobr
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))
!
          pi                  = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOBr on PSC Type I NAT particles = BrCl + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc           = 1.0d0
          gprob             = 0.10d0
          avgvel(:)         = 100.0d0 *  &
     &                        (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                        (pi * mw(IHOBR)))**0.5d0
!
          sknat_hcl_hobr(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_hobr  = sknat_hcl_hobr / minconc
          elsewhere
            sknat_hcl_hobr  = sknat_hcl_hobr / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop )  &
     &        sknat_hcl_hobr = 0.0d0
          end if
        END FUNCTION sknat_hcl_hobr
!
!.... skice_clono2 (temperature ,sad_ice)
!
!.... (1) JPL 02-25
!
        FUNCTION skice_clono2 (tk,sad)
          real*8  &
     &      tk(:) ,sad(:)
          real*8, dimension(size(tk)) :: skice_clono2
          real*8  &
     &      gprob ,pi
          real*8  &
     &      avgvel(size(tk))

          pi              = acos(-1.0d0)

!=======================================================================
!     ClONO2 + PSC Type II ice particles = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob           = 0.30d0
          avgvel(:)       = 100.0d0 *  &
     &                      (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(ICLONO2)))**0.5d0

          skice_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
        END FUNCTION skice_clono2
!
!.... skice_brono2 (temperature ,sad_ice)
!
!.... (2) JPL 00-003
!
        FUNCTION skice_brono2 (tk,sad)
          real*8  &
     &      tk(:) ,sad(:)
          real*8, dimension(size(tk)) :: skice_brono2
          real*8  &
     &      gprob ,pi
          real*8  &
     &      avgvel(size(tk))

          pi              = acos(-1.0d0)

!=======================================================================
!     BrONO2 + PSC Type II ice particles = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob           = 0.30d0
          avgvel(:)       = 100.0d0 *  &
     &                      (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(IBRONO2)))**0.5d0

          skice_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
        END FUNCTION skice_brono2
!
!.... skice_hcl_clono2 (temperature ,sad_ice ,specarr( HCl,:))
!
!.... (3) JPL 00-003
!
        FUNCTION skice_hcl_clono2 (tk,sad,hcl)
          real*8  &
     &      tk(:) ,sad(:) ,hcl(:)
          real*8, dimension(size(tk)) :: skice_hcl_clono2
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))

          pi                  = acos(-1.0d0)

!=======================================================================
!     HCl + ClONO2 on PSC Type II ice particles = Cl2 + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc             = 1.0d0
          gprob               = 0.30d0
          avgvel(:)           = 100.0d0 *  &
     &                          (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                          (pi * mw(ICLONO2)))**0.5d0

          skice_hcl_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)

          where (hcl < minconc)
            skice_hcl_clono2  = skice_hcl_clono2 / minconc
          elsewhere
            skice_hcl_clono2  = skice_hcl_clono2 / hcl
          end where
        END FUNCTION skice_hcl_clono2
!
!.... skice_hcl_hocl (temperature ,sad_ice ,specarr( HCl,:))
!
!.... (4) JPL 00-003
!
        FUNCTION skice_hcl_hocl (tk,sad,hcl)
          real*8  &
     &      tk(:) ,sad(:) ,hcl(:)
          real*8, dimension(size(tk)) :: skice_hcl_hocl
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))

          pi                  = acos(-1.0d0)

!=======================================================================
!     HCl + HOCl on PSC Type II ice particles = Cl2 + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc           = 1.0d0
          gprob             = 0.20d0
          avgvel(:)         = 100.0d0 *  &
     &                        (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                        (pi * mw(IHOCL)))**0.5d0

          skice_hcl_hocl(:) = 0.25d0 * gprob * avgvel(:) * sad(:)

          where (hcl < minconc)
            skice_hcl_hocl  = skice_hcl_hocl / minconc
          elsewhere
            skice_hcl_hocl  = skice_hcl_hocl / hcl
          end where
        END FUNCTION skice_hcl_hocl
!
!.... skice_hcl_brono2 (temperature ,sad_ice ,specarr( HCl,:))
!
!.... (5) JPL 00-003
!
        FUNCTION skice_hcl_brono2 (tk,sad,hcl)
          real*8  &
     &      tk(:) ,sad(:) ,hcl(:)
          real*8, dimension(size(tk)) :: skice_hcl_brono2
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))

          pi                  = acos(-1.0d0)

!=======================================================================
!     HCl + BrONO2 on PSC Type II ice particles = BrCl + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc             = 1.0d0
          gprob               = 0.30d0
          avgvel(:)           = 100.0d0 *  &
     &                          (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                          (pi * mw(IBRONO2)))**0.5d0

          skice_hcl_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)

          where (hcl < minconc)
            skice_hcl_brono2  = skice_hcl_brono2 / minconc
          elsewhere
            skice_hcl_brono2  = skice_hcl_brono2 / hcl
          end where
        END FUNCTION skice_hcl_brono2
!
!.... skice_hcl_hobr (temperature ,sad_ice ,specarr( HCl,:))
!
!.... (6) JPL 00-003
!
        FUNCTION skice_hcl_hobr (tk,sad,hcl)
          real*8  &
     &      tk(:) ,sad(:) ,hcl(:)
          real*8, dimension(size(tk)) :: skice_hcl_hobr
          real*8  &
     &      gprob ,minconc ,pi
          real*8  &
     &      avgvel(size(tk))

          pi                  = acos(-1.0d0)

!=======================================================================
!     HCl + HOBr on PSC Type II ice particles = BrCl + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc           = 1.0d0
          gprob             = 0.30d0
          avgvel(:)         = 100.0d0 *  &
     &                        (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                        (pi * mw(IHOBR)))**0.5d0

          skice_hcl_hobr(:) = 0.25d0 * gprob * avgvel(:) * sad(:)

          where (hcl < minconc)
            skice_hcl_hobr  = skice_hcl_hobr / minconc
          elsewhere
            skice_hcl_hobr  = skice_hcl_hobr / hcl
          end where
        END FUNCTION skice_hcl_hobr
!
!.... sksoot_hno3 (temperature ,sad_soot)
!
!.... (1)
!
        FUNCTION sksoot_hno3 (tk,sad)
          real*8  &
     &      tk(:) ,sad(:)
          real*8, dimension(size(tk)) :: sksoot_hno3
          real*8  &
     &      gprob ,pi
          real*8  &
     &      avgvel(size(tk))

          pi              = acos(-1.0d0)

!=======================================================================
!     HNO3 + soot particles = OH + NO2
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 1/16/2002
!
          gprob           = 0.0d0
          avgvel(:)       = 100.0d0 *  &
     &                      (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(IHNO3)))**0.5d0

          sksoot_hno3(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
        END FUNCTION sksoot_hno3
      END
