!=======================================================================
!
! $Id: setkin_kcalc.F90,v 1.4 2011-08-09 22:12:57 mrdamon Exp $
!   kcalc - IMPACT 3D model (setkin_kcalc.F)
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
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
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

!cc   INTEGER               :: naltmax

      real*8  :: mw(NSP)

      real*8  &
!cc  &  nitrogen (npres0)
!cc  & ,oxygen   (npres0)
     &  water    (npres0)

      real*8  &
!cc  &  sad_ice  (npres0)
     &  sad_lbs  (npres0)
!cc  & ,sad_nat  (npres0)
!cc  & ,sad_soot (npres0)
!cc  & ,sad_sts  (npres0)

!cc   naltmax     = npres0

      mw(:) = mw_data(:)

!....          Define molecular nitrogen and oxygen number densities
!
!cc   nitrogen(:) = adcol(:) * MXRN2
!cc   oxygen(:)   = adcol(:) * MXRO2
      water(:)    = specarr(8 ,:)

      sad_lbs(:)  = sadcol(1, :)
!cc   sad_sts(:)  = sadcol(2, :)
!cc   sad_nat(:)  = sadcol(3, :)
!cc   sad_ice(:)  = sadcol(4, :)
!cc   sad_soot(:) = sadcol(5, :)
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
      rcarr(3,:) = skarr( 1.800D-11 ,-110.0D+00 ,temperature)
!
!....           O1D + O2 = O + O2
!
      rcarr(4,:) = skarr( 3.200D-11 ,-70.0D+00 ,temperature)
!
!....           O1D + O3 = 2 O2
!
      rcarr(5,:) = skarr( 1.200D-10 ,0.0D+00 ,temperature)
!
!....           H2 + O1D = H + OH
!
      rcarr(6,:) = skarr( 1.100D-10 ,0.0D+00 ,temperature)
!
!....           H2O + O1D = 2 OH
!
      rcarr(7,:) = skarr( 2.200D-10 ,0.0D+00 ,temperature)
!
!....           H + O2 = HO2
!
      rcarr(8,:) = sktroe( 5.700D-32 ,1.60D0  &
     &                     ,7.500D-11 ,0.00D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           H + O3 = O2 + OH
!
      rcarr(9,:) = skarr( 1.400D-10 ,470.0D+00 ,temperature)
!
!....           H + HO2 = 2 OH
!
      rcarr(10,:) = skarr( 7.000D-11 ,0.0D+00 ,temperature)
!
!....           O + OH = H + O2
!
      rcarr(11,:) = skarr( 2.200D-11 ,-120.0D+00 ,temperature)
!
!....           O3 + OH = HO2 + O2
!
      rcarr(12,:) = skarr( 1.500D-12 ,880.0D+00 ,temperature)
!
!....           H2 + OH = H + H2O
!
      rcarr(13,:) = skarr( 5.500D-12 ,2000.0D+00 ,temperature)
!
!....           OH + OH = H2O + O
!
      rcarr(14,:) = skarr( 4.200D-12 ,240.0D+00 ,temperature)
!
!....           HO2 + O = O2 + OH
!
      rcarr(15,:) = skarr( 3.000D-11 ,-200.0D+00 ,temperature)
!
!....           HO2 + O3 = 2 O2 + OH
!
      rcarr(16,:) = skarr( 2.000D-14 ,680.0D+00 ,temperature)
!
!....           HO2 + OH = H2O + O2
!
      rcarr(17,:) = skarr( 4.800D-11 ,-250.0D+00 ,temperature)
!
!....           HO2 + HO2 = H2O2 + O2
!
      rcarr(18,:) = skho2dis (temperature ,adcol)
!
!....           H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      rcarr(19,:) = skho2h2o (temperature ,adcol)
!
!....           H2O2 + OH = H2O + HO2
!
      rcarr(20,:) = skarr( 2.900D-12 ,160.0D+00 ,temperature)
!
!....           N2O + O1D = N2 + O2
!
      rcarr(21,:) = skarr( 4.900D-11 ,0.0D+00 ,temperature)
!
!....           N2O + O1D = 2 NO
!
      rcarr(22,:) = skarr( 6.700D-11 ,0.0D+00 ,temperature)
!
!....           N + O2 = NO + O
!
      rcarr(23,:) = skarr( 1.500D-11 ,3600.0D+00 ,temperature)
!
!....           NO + O3 = NO2 + O2
!
      rcarr(24,:) = skarr( 3.000D-12 ,1500.0D+00 ,temperature)
!
!....           NO2 + O = NO + O2
!
      rcarr(25,:) = skarr( 5.600D-12 ,-180.0D+00 ,temperature)
!
!....           NO2 + O3 = NO3 + O2
!
      rcarr(26,:) = skarr( 1.200D-13 ,2450.0D+00 ,temperature)
!
!....           NO + OH = HONO
!
      rcarr(27,:) = sktroe( 7.000D-31 ,2.60D0  &
     &                     ,3.600D-11 ,0.10D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           HO2 + NO = NO2 + OH
!
      rcarr(28,:) = skarr( 3.500D-12 ,-250.0D+00 ,temperature)
!
!....           NO2 + OH = HNO3
!
      rcarr(29,:) = sktroe( 2.400D-30 ,3.10D0  &
     &                     ,1.700D-11 ,2.10D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           HO2 + NO2 = HO2NO2
!
      rcarr(30,:) = sktroe( 1.800D-31 ,3.20D0  &
     &                     ,4.700D-12 ,1.40D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           H2O + N2O5 = 2 HNO3
!
      rcarr(31,:) = skarr( 5.000D-22 ,0.0D+00 ,temperature)
!
!....           HO2NO2 = HO2 + NO2
!
      rcarr(32,:) = sktroe( 8.570D-05 ,3.20D0  &
     &                     ,2.240D+15 ,1.40D0 ,10900.0D+00  &
     &                     ,temperature ,adcol)
!
!....           HONO + OH = H2O + NO2
!
      rcarr(33,:) = skarr( 1.800D-11 ,390.0D+00 ,temperature)
!
!....           HNO3 + OH = H2O + NO3
!
      rcarr(34,:) = skohhno3 (temperature ,adcol)
!
!....           HO2NO2 + OH = H2O + NO2 + O2
!
      rcarr(35,:) = skarr( 1.300D-12 ,-380.0D+00 ,temperature)
!
!....           N + NO = N2 + O
!
      rcarr(36,:) = skarr( 2.100D-11 ,-100.0D+00 ,temperature)
!
!....           NO + NO3 = 2 NO2
!
      rcarr(37,:) = skarr( 1.500D-11 ,-170.0D+00 ,temperature)
!
!....           NO2 + NO3 = N2O5
!
      rcarr(38,:) = sktroe( 2.000D-30 ,4.40D0  &
     &                     ,1.400D-12 ,0.70D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           N2O5 = NO2 + NO3
!
      rcarr(39,:) = sktroe( 6.670D-04 ,4.40D0  &
     &                     ,4.670D+14 ,0.70D0 ,10991.0D+00  &
     &                     ,temperature ,adcol)
!
!....           Cl + O3 = ClO + O2
!
      rcarr(40,:) = skarr( 2.300D-11 ,200.0D+00 ,temperature)
!
!....           ClO + O = Cl + O2
!
      rcarr(41,:) = skarr( 3.000D-11 ,-70.0D+00 ,temperature)
!
!....           ClONO2 + O = ClO + NO3
!
      rcarr(42,:) = skarr( 4.500D-12 ,900.0D+00 ,temperature)
!
!....           Cl + H2 = H + HCl
!
      rcarr(43,:) = skarr( 3.700D-11 ,2300.0D+00 ,temperature)
!
!....           Cl + HO2 = HCl + O2
!
      rcarr(44,:) = skarr( 1.800D-11 ,-170.0D+00 ,temperature)
!
!....           Cl + HO2 = ClO + OH
!
      rcarr(45,:) = skarr( 4.100D-11 ,450.0D+00 ,temperature)
!
!....           Cl + H2O2 = HCl + HO2
!
      rcarr(46,:) = skarr( 1.100D-11 ,980.0D+00 ,temperature)
!
!....           ClO + OH = Cl + HO2
!
      rcarr(47,:) = skarr( 7.400D-12 ,-270.0D+00 ,temperature)
!
!....           ClO + OH = HCl + O2
!
      rcarr(48,:) = skarr( 3.200D-13 ,-320.0D+00 ,temperature)
!
!....           ClO + HO2 = HOCl + O2
!
      rcarr(49,:) = skarr( 4.800D-13 ,-700.0D+00 ,temperature)
!
!....           ClO + HO2 = HCl + O3
!
      rcarr(50,:) = skarr( 0.000D+00 ,0.0D+00 ,temperature)
!
!....           HCl + OH = Cl + H2O
!
      rcarr(51,:) = skarr( 2.600D-12 ,350.0D+00 ,temperature)
!
!....           HOCl + OH = ClO + H2O
!
      rcarr(52,:) = skarr( 3.000D-12 ,500.0D+00 ,temperature)
!
!....           ClONO2 + OH = HOCl + NO3
!
      rcarr(53,:) = skarr( 1.200D-12 ,330.0D+00 ,temperature)
!
!....           ClO + NO = Cl + NO2
!
      rcarr(54,:) = skarr( 6.400D-12 ,-290.0D+00 ,temperature)
!
!....           ClO + NO2 = ClONO2
!
      rcarr(55,:) = sktroe( 1.800D-31 ,3.40D0  &
     &                     ,1.500D-11 ,1.90D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           CH4 + Cl = CH3O2 + HCl
!
      rcarr(56,:) = skarr( 1.100D-11 ,1400.0D+00 ,temperature)
!
!....           CH2O + Cl = CO + HCl + HO2
!
      rcarr(57,:) = skarr( 8.100D-11 ,30.0D+00 ,temperature)
!
!....           Cl + HOCl = ClO + HCl
!
      rcarr(58,:) = skarr( 2.500D-12 ,130.0D+00 ,temperature)
!
!....           ClO + ClO = 2 Cl + O2
!
      rcarr(59,:) = skarr( 3.000D-11 ,2450.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2 + O2
!
      rcarr(60,:) = skarr( 1.000D-12 ,1590.0D+00 ,temperature)
!
!....           ClO + ClO = Cl + OClO
!
      rcarr(61,:) = skarr( 3.500D-13 ,1370.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2O2
!
      rcarr(62,:) = sktroe( 2.200D-32 ,3.10D0  &
     &                     ,3.400D-12 ,1.00D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           Cl2O2 = 2 ClO
!
      rcarr(63,:) = sktroe( 1.690D-05 ,3.10D0  &
     &                     ,2.690D+15 ,1.00D0 ,8744.0D+00  &
     &                     ,temperature ,adcol)
!
!....           Cl + ClONO2 = Cl2 + NO3
!
      rcarr(64,:) = skarr( 6.500D-12 ,-135.0D+00 ,temperature)
!
!....           Br + O3 = BrO + O2
!
      rcarr(65,:) = skarr( 1.700D-11 ,800.0D+00 ,temperature)
!
!....           BrO + O = Br + O2
!
      rcarr(66,:) = skarr( 1.900D-11 ,-230.0D+00 ,temperature)
!
!....           Br + HO2 = HBr + O2
!
      rcarr(67,:) = skarr( 1.500D-11 ,600.0D+00 ,temperature)
!
!....           BrO + OH = Br + HO2
!
      rcarr(68,:) = skarr( 7.500D-11 ,0.0D+00 ,temperature)
!
!....           BrO + HO2 = HOBr + O2
!
      rcarr(69,:) = skarr( 3.400D-12 ,-540.0D+00 ,temperature)
!
!....           HBr + OH = Br + H2O
!
      rcarr(70,:) = skarr( 1.100D-11 ,0.0D+00 ,temperature)
!
!....           BrO + NO = Br + NO2
!
      rcarr(71,:) = skarr( 8.800D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + NO2 = BrONO2
!
      rcarr(72,:) = sktroe( 5.200D-31 ,3.20D0  &
     &                     ,6.900D-12 ,2.90D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           Br + CH2O = CO + HBr + HO2
!
      rcarr(73,:) = skarr( 1.700D-11 ,800.0D+00 ,temperature)
!
!....           BrO + ClO = Br + OClO
!
      rcarr(74,:) = skarr( 9.500D-13 ,-550.0D+00 ,temperature)
!
!....           BrO + ClO = Br + Cl + O2
!
      rcarr(75,:) = skarr( 2.900D-12 ,-220.0D+00 ,temperature)
!
!....           BrO + ClO = BrCl + O2
!
      rcarr(76,:) = skarr( 4.100D-13 ,-290.0D+00 ,temperature)
!
!....           BrO + BrO = 2 Br + O2
!
      rcarr(77,:) = skbrodis (temperature)
!
!....           CF2Cl2 + O1D = 2 Cl
!
      rcarr(78,:) = skarr( 1.200D-10 ,0.0D+00 ,temperature)
!
!....           CH3Cl + OH = Cl + H2O + HO2
!
      rcarr(79,:) = skarr( 4.000D-12 ,1400.0D+00 ,temperature)
!
!....           CH3Br + OH = Br + H2O
!
      rcarr(80,:) = skarr( 4.000D-12 ,1470.0D+00 ,temperature)
!
!....           CH3CCl3 + OH = 3 Cl + H2O
!
      rcarr(81,:) = skarr( 1.800D-12 ,1550.0D+00 ,temperature)
!
!....           CH3Cl + Cl = CO + 2 HCl + HO2
!
      rcarr(82,:) = skarr( 3.200D-11 ,1250.0D+00 ,temperature)
!
!....           N2O5 = 2 HNO3
!
      rcarr(83,:) = sksts_n2o5 ( temperature  &
     &                          ,pressure ,sad_lbs ,tropp)
!
!....           N2O5 = 2 HNO3
!
      rcarr(84,:) = sktrs_n2o5 ( temperature  &
     &                          ,pressure ,sad_lbs ,tropp)
!
!....           NO3 = HNO3
!
      rcarr(85,:) = sktrs_no3 ( temperature  &
     &                         ,pressure ,sad_lbs ,tropp)
!
!....           H2O + NO2 = HNO3 + HONO
!
      rcarr(86,:) = skarr( 4.000D-24 ,0.0D+00 ,temperature)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(87,:) = sksts_clono2 ( temperature  &
     &           ,adcol ,pressure ,sad_lbs ,specarr(  24,:) ,water  &
     &           ,tropp)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(88,:) = sksts_brono2 ( temperature  &
     &                            ,pressure ,sad_lbs ,tropp)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(89,:) = sksts_clono2_hcl (temperature  &
     &           ,adcol ,pressure ,sad_lbs ,specarr(    26,:)  &
     &           ,specarr(  24,:) ,water ,tropp)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(90,:) = sksts_hocl_hcl (temperature  &
     &           ,adcol ,pressure ,sad_lbs ,specarr(  25,:)  &
     &           ,specarr( 24,:) ,water ,tropp)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(91,:) = sksts_hobr_hcl (temperature  &
     &           ,adcol ,pressure ,sad_lbs ,specarr(  30,:)  &
     &           ,specarr( 24,:) ,water ,tropp)
!
!....           CH4 + O1D = CH3O2 + OH
!
      rcarr(92,:) = skarr( 1.125D-10 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H + HO2
!
      rcarr(93,:) = skarr( 3.000D-11 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H2
!
      rcarr(94,:) = skarr( 7.500D-12 ,0.0D+00 ,temperature)
!
!....           CH2O + O = CO + HO2 + OH
!
      rcarr(95,:) = skarr( 3.400D-11 ,1600.0D+00 ,temperature)
!
!....           CH4 + OH = CH3O2 + H2O
!
      rcarr(96,:) = skarr( 2.450D-12 ,1775.0D+00 ,temperature)
!
!....           CO + OH = H
!
      rcarr(97,:) = skcooh (temperature ,pressure)
!
!....           CH2O + OH = CO + H2O + HO2
!
      rcarr(98,:) = skarr( 1.000D-11 ,0.0D+00 ,temperature)
!
!....           CH3OH + OH = CH2O + HO2
!
      rcarr(99,:) = skarr( 6.700D-12 ,600.0D+00 ,temperature)
!
!....           CH3OOH + OH = CH3O2 + H2O
!
      rcarr(100,:) = skarr( 2.700D-12 ,-200.0D+00 ,temperature)
!
!....           CH3OOH + OH = CH2O + H2O + OH
!
      rcarr(101,:) = skarr( 1.100D-12 ,-200.0D+00 ,temperature)
!
!....           CH3O2 + HO2 = CH3OOH + O2
!
      rcarr(102,:) = skarr( 3.800D-13 ,-800.0D+00 ,temperature)
!
!....           CH3O2 + NO = CH2O + HO2 + NO2
!
      rcarr(103,:) = skarr( 3.000D-12 ,-280.0D+00 ,temperature)
!
!....           CH3O2 + NO2 = CH3O2NO2
!
      rcarr(104,:) = sktroe( 1.500D-30 ,4.00D0  &
     &                     ,6.500D-12 ,2.00D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           CH3O2NO2 = CH3O2 + NO2
!
      rcarr(105,:) = sktroe( 1.150D-02 ,4.00D0  &
     &                     ,5.000D+16 ,2.00D0 ,11200.0D+00  &
     &                     ,temperature ,adcol)
!
!....           CH3O2 + CH3O2 =  1.33 CH2O +  0.66 CH3OH +  0.80 HO2
!
      rcarr(106,:) = skarr( 2.500D-13 ,-190.0D+00 ,temperature)
!
!....           CH2O + HO2 = CH3O3
!
      rcarr(107,:) = skarr( 6.700D-15 ,-600.0D+00 ,temperature)
!
!....           CH3O3 = CH2O + HO2
!
      rcarr(108,:) = skarr( 2.400D+12 ,7000.0D+00 ,temperature)
!
!....           CH3O3 + NO = HCOOH + HO2 + NO2
!
      rcarr(109,:) = skarr( 3.000D-12 ,-280.0D+00 ,temperature)
!
!....           CH3O3 + HO2 = HCOOH
!
      rcarr(110,:) = skarr( 5.600D-15 ,-2300.0D+00 ,temperature)
!
!....           CH3O3 + CH3O3 = 2 HCOOH + 2 HO2
!
      rcarr(111,:) = skarr( 5.700D-14 ,-750.0D+00 ,temperature)
!
!....           HCOOH + OH = CO + HO2
!
      rcarr(112,:) = skarr( 4.000D-13 ,0.0D+00 ,temperature)
!
!....           C2H6 + OH = ETO2
!
      rcarr(113,:) = skarr( 8.700D-12 ,1070.0D+00 ,temperature)
!
!....           C2H6 + Cl = ETO2 + HCl
!
      rcarr(114,:) = skarr( 7.700D-11 ,90.0D+00 ,temperature)
!
!....           ETO2 + NO = ALD2 + HO2 + NO2
!
      rcarr(115,:) = skarr( 2.600D-12 ,-365.0D+00 ,temperature)
!
!....           ETO2 + ETO2 =  1.60 ALD2 +  1.60 HO2
!
      rcarr(116,:) = skarr( 6.800D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + HO2 = ETP
!
      rcarr(117,:) = skarr( 7.500D-13 ,-700.0D+00 ,temperature)
!
!....           ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH
!
      rcarr(118,:) = skarr( 3.000D-12 ,-190.0D+00 ,temperature)
!
!....           ALD2 + OH = MCO3
!
      rcarr(119,:) = skarr( 5.600D-12 ,-270.0D+00 ,temperature)
!
!....           ALD2 + NO3 = HNO3 + MCO3
!
      rcarr(120,:) = skarr( 1.400D-12 ,1900.0D+00 ,temperature)
!
!....           MCO3 + NO2 = PAN
!
      rcarr(121,:) = sktroe( 9.700D-29 ,5.60D0  &
     &                     ,9.300D-12 ,1.50D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           MCO3 + NO = CH3O2 + NO2
!
      rcarr(122,:) = skarr( 5.300D-12 ,-360.0D+00 ,temperature)
!
!....           HO2 + MCO3 =  0.75 MAP +  0.25 O3
!
      rcarr(123,:) = skarr( 4.500D-13 ,-1000.0D+00 ,temperature)
!
!....           MCO3 + MCO3 = 2 CH3O2
!
      rcarr(124,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....           CH3O2 + MCO3 = CH2O +  0.90 CH3O2 +  0.90 HO2
!
      rcarr(125,:) = skarr( 1.300D-12 ,-640.0D+00 ,temperature)
!
!....           PAN = MCO3 + NO2
!
      rcarr(126,:) = sktroe( 1.078D+00 ,5.60D0  &
     &                     ,1.022D+17 ,1.50D0 ,14000.0D+00  &
     &                     ,temperature ,adcol)
!
!....           GLYX + OH = 2 CO + HO2
!
      rcarr(127,:) = skarr( 1.100D-11 ,0.0D+00 ,temperature)
!
!....           C3H8 + OH = A3O2
!
      rcarr(128,:) = skarr( 1.000D-11 ,660.0D+00 ,temperature)
!
!....           C3H8 + Cl = A3O2 + HCl
!
      rcarr(129,:) = skarr( 1.200D-10 ,-40.0D+00 ,temperature)
!
!....           A3O2 + NO =  0.82 ACET +  0.97 HO2 +  0.97 NO2 +  0.16 RCHO
!
      rcarr(130,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           A3O2 + A3O2 =  0.40 ALD2 +  1.20 HO2 +  1.20 RCHO
!
      rcarr(131,:) = skarr( 3.000D-13 ,0.0D+00 ,temperature)
!
!....           A3O2 + HO2 = RA3P
!
      rcarr(132,:) = skarr( 1.660D-13 ,-1300.0D+00 ,temperature)
!
!....           OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      rcarr(133,:) = skarr( 3.000D-12 ,-190.0D+00 ,temperature)
!
!....           ACET + OH = ATO2
!
      rcarr(134,:) = skarr( 2.200D-12 ,685.0D+00 ,temperature)
!
!....           ATO2 + NO =  0.96 HO2 +  0.96 MGLY +  0.96 NO2
!
      rcarr(135,:) = skarr( 8.000D-12 ,0.0D+00 ,temperature)
!
!....           ATO2 + HO2 = CH3O2 + MCO3
!
      rcarr(136,:) = skarr( 1.150D-13 ,-1300.0D+00 ,temperature)
!
!....           MGLY + OH = CO + MCO3
!
      rcarr(137,:) = skarr( 1.500D-11 ,0.0D+00 ,temperature)
!
!....           OH + RCHO = ECO3
!
      rcarr(138,:) = skarr( 2.000D-11 ,0.0D+00 ,temperature)
!
!....           NO3 + RCHO = ECO3 + HNO3
!
      rcarr(139,:) = skarr( 5.700D-15 ,0.0D+00 ,temperature)
!
!....           ECO3 + NO2 = PPN
!
      rcarr(140,:) = sktroe( 9.700D-29 ,5.60D0  &
     &                     ,9.300D-12 ,1.50D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           ECO3 + NO = ETO2 + NO2
!
      rcarr(141,:) = skarr( 2.000D-11 ,0.0D+00 ,temperature)
!
!....           ECO3 + ECO3 = 2 ETO2
!
      rcarr(142,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....           CH3O2 + ECO3 =  0.15 ALD2 + CH2O +  0.85 ETO2 +  0.65 HO2
!
      rcarr(143,:) = skarr( 1.300D-12 ,-640.0D+00 ,temperature)
!
!....           ECO3 + MCO3 = CH3O2 + ETO2
!
      rcarr(144,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....           ECO3 + HO2 =  0.25 ALD2 +  0.25 HO2 +  0.25 O3 +  0.75 RP
!
      rcarr(145,:) = skarr( 4.500D-13 ,-1000.0D+00 ,temperature)
!
!....           PPN = ECO3 + NO2
!
      rcarr(146,:) = sktroe( 1.078D+00 ,5.60D0  &
     &                     ,1.022D+17 ,1.50D0 ,14000.0D+00  &
     &                     ,temperature ,adcol)
!
!....           OH + RP =  0.50 ALD2 +  0.50 ECO3 +  0.50 OH
!
      rcarr(147,:) = skarr( 3.000D-12 ,-190.0D+00 ,temperature)
!
!....           ISOP + OH = RIO2
!
      rcarr(148,:) = skarr( 2.550D-11 ,-409.0D+00 ,temperature)
!
!....           MVK + OH = VRO2
!
      rcarr(149,:) = skarr( 2.670D-12 ,-452.0D+00 ,temperature)
!
!....           MACR + OH = MAO3
!
      rcarr(150,:) = skarr( 3.900D-12 ,-379.0D+00 ,temperature)
!
!....           MACR + OH = MRO2
!
      rcarr(151,:) = skarr( 3.900D-12 ,-379.0D+00 ,temperature)
!
!....           HAC + OH = GLYX + HO2
!
      rcarr(152,:) = skarr( 6.700D-13 ,-270.0D+00 ,temperature)
!
!....           HAC + OH = HACO
!
      rcarr(153,:) = skarr( 2.380D-12 ,-270.0D+00 ,temperature)
!
!....           ISOP + O3 =  0.80 CH2O +  0.07 CH3OH +  0.05 CO +  0.06 HO2
!
      rcarr(154,:) = skarr( 5.590D-15 ,1814.0D+00 ,temperature)
!
!....           MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.16 CO +  0.06 HO2 +
!
      rcarr(155,:) = skarr( 6.910D-16 ,1519.0D+00 ,temperature)
!
!....           MACR + O3 =  0.70 CH2O +  0.09 CO +  0.21 HO2 +  0.80 MGLY +
!
      rcarr(156,:) = skarr( 1.300D-15 ,2112.0D+00 ,temperature)
!
!....           O3 + RIP =  0.70 CH2O
!
      rcarr(157,:) = skarr( 8.000D-18 ,0.0D+00 ,temperature)
!
!....           O3 + VRP =  0.70 CH2O
!
      rcarr(158,:) = skarr( 8.000D-18 ,0.0D+00 ,temperature)
!
!....           MRP + O3 =  0.70 CH2O
!
      rcarr(159,:) = skarr( 8.000D-18 ,0.0D+00 ,temperature)
!
!....           ISOP + NO3 = INO2
!
      rcarr(160,:) = skarr( 3.030D-12 ,450.0D+00 ,temperature)
!
!....           MACR + NO3 = HNO3 + MAO3
!
      rcarr(161,:) = skarr( 1.100D-15 ,0.0D+00 ,temperature)
!
!....           MACR + NO3 = MAN2
!
      rcarr(162,:) = skarr( 2.200D-15 ,0.0D+00 ,temperature)
!
!....           NO + RIO2 =  0.74 CH2O +  0.78 HO2 +  0.14 ISN1 +  0.32 MACR
!
      rcarr(163,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           NO + VRO2 =  0.27 CH2O +  0.68 HAC +  0.27 HO2 +  0.05 ISN1
!
      rcarr(164,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           MAO3 + NO = MAO2
!
      rcarr(165,:) = skarr( 8.690D-12 ,-290.0D+00 ,temperature)
!
!....           MAO2 + NO = CH2O + MCO3 + NO2
!
      rcarr(166,:) = skarr( 4.200D-12 ,-180.0D+00 ,temperature)
!
!....           MRO2 + NO = HACN + HO2 + NO2
!
      rcarr(167,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           HACO + NO = CH3O3 + NO2
!
      rcarr(168,:) = skarr( 5.300D-12 ,-360.0D+00 ,temperature)
!
!....           INO2 + NO =  0.15 CH2O +  0.80 HO2 +  0.75 ISN1 +  0.10 MACR
!
      rcarr(169,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           ISNR + NO =  0.95 ACET +  1.90 HAC +  0.05 HO2 +  0.05 ISN1
!
      rcarr(170,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           MAN2 + NO = HO2 + MGLY + 2 NO2
!
      rcarr(171,:) = skarr( 2.900D-12 ,-350.0D+00 ,temperature)
!
!....           HO2 + RIO2 = RIP
!
      rcarr(172,:) = skarr( 1.660D-13 ,-1300.0D+00 ,temperature)
!
!....           HO2 + VRO2 = VRP
!
      rcarr(173,:) = skarr( 1.660D-13 ,-1300.0D+00 ,temperature)
!
!....           HO2 + MAO3 = RP
!
      rcarr(174,:) = skarr( 1.150D-12 ,-550.0D+00 ,temperature)
!
!....           HO2 + MRO2 = MRP
!
      rcarr(175,:) = skarr( 1.660D-13 ,-1300.0D+00 ,temperature)
!
!....           HACO + HO2 = CH3OOH
!
      rcarr(176,:) = skarr( 1.150D-12 ,-550.0D+00 ,temperature)
!
!....           HO2 + ISNR = PRN2
!
      rcarr(177,:) = skarr( 1.400D-13 ,-1380.0D+00 ,temperature)
!
!....           HO2 + INO2 = PRN2
!
      rcarr(178,:) = skarr( 1.400D-13 ,-1380.0D+00 ,temperature)
!
!....           HO2 + MAN2 = PRN2
!
      rcarr(179,:) = skarr( 1.400D-13 ,-1380.0D+00 ,temperature)
!
!....           INO2 + NO2 = 2 PRN2
!
      rcarr(180,:) = skarr( 4.200D-13 ,-180.0D+00 ,temperature)
!
!....           HACN + OH = HO2 + MGLY
!
      rcarr(181,:) = skarr( 3.000D-12 ,0.0D+00 ,temperature)
!
!....           ISN1 + OH = ISNR
!
      rcarr(182,:) = skarr( 3.350D-11 ,0.0D+00 ,temperature)
!
!....           OH + RIP =  0.37 CH2O + HO2 +  0.16 MACR +  0.21 MVK +  0.50
!
      rcarr(183,:) = skarr( 3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + VRP =  0.50 HAC + HO2 + VRO2
!
      rcarr(184,:) = skarr( 3.800D-12 ,-200.0D+00 ,temperature)
!
!....           MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH
!
      rcarr(185,:) = skarr( 3.000D-12 ,-190.0D+00 ,temperature)
!
!....           MRP + OH =  0.50 HAC + HO2 +  0.50 MRO2
!
      rcarr(186,:) = skarr( 3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + PRN2 =  0.50 ISNR +  0.50 NO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(187,:) = skarr( 3.800D-12 ,-200.0D+00 ,temperature)
!
!....           HO2 + MAO3 = CH2O + ETO2 + O3
!
      rcarr(188,:) = skarr( 3.860D-16 ,-2640.0D+00 ,temperature)
!
!....           HO2 + MAO2 = ACET + O3
!
      rcarr(189,:) = skarr( 3.860D-16 ,-2640.0D+00 ,temperature)
!
!....           HACO + HO2 = ETO2 + O3
!
      rcarr(190,:) = skarr( 3.860D-16 ,-2640.0D+00 ,temperature)
!
!....           MAO3 + NO2 = MPAN
!
      rcarr(191,:) = sktroe( 9.700D-29 ,5.60D0  &
     &                     ,9.300D-12 ,1.50D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           MPAN = MAO3 + NO2
!
      rcarr(192,:) = sktroe( 1.078D+00 ,5.60D0  &
     &                     ,1.022D+17 ,1.50D0 ,14000.0D+00  &
     &                     ,temperature ,adcol)
!
!....           HACO + NO2 = IPAN
!
      rcarr(193,:) = sktroe( 9.700D-29 ,5.60D0  &
     &                     ,9.300D-12 ,1.50D0 ,0.0D0  &
     &                     ,temperature ,adcol)
!
!....           IPAN = HACO + NO2
!
      rcarr(194,:) = sktroe( 1.078D+00 ,5.60D0  &
     &                     ,1.022D+17 ,1.50D0 ,14000.0D+00  &
     &                     ,temperature ,adcol)
!
!....           CH3O2 + MAO3 = CH2O +  0.50 HO2 +  0.85 MAO2
!
      rcarr(195,:) = skarr( 1.300D-12 ,-640.0D+00 ,temperature)
!
!....           MAO3 + MCO3 = CH3O2 + MAO2
!
      rcarr(196,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....           MAO3 + MAO3 = 2 MAO2
!
      rcarr(197,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....           CH3O2 + HACO = CH2O +  0.85 CH3O2 +  0.50 HO2
!
      rcarr(198,:) = skarr( 1.300D-12 ,-640.0D+00 ,temperature)
!
!....           HACO + MCO3 = 2 CH3O2
!
      rcarr(199,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....           HACO + HACO = 2 CH3O2
!
      rcarr(200,:) = skarr( 2.900D-12 ,-500.0D+00 ,temperature)
!
!....          End thermal rate constants
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
     &                  1.4d-21 * exp(2200.d0 / tk(:))
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
          real*8, dimension(size(tk)) :: skohhno3
          real*8  &
     &      xxx1(size(tk)) ,xxx2(size(tk))
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
          skbrodis(:) = 1.5D-12 * exp(230.0d0 / tk(:))
        END FUNCTION skbrodis
!
!.... sksts_n2o5 (temperature ,pressure ,sad_lbs ,tropp)
!
!.... JPL 97-4
!
        FUNCTION sksts_n2o5 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:)
          real*8, dimension(size(tk)) :: sksts_n2o5
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk))

          pi           = acos(-1.0d0)

!=======================================================================
!     N2O5 + stratospheric sulfate aerosol = 2 HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 3/30/99
!
          gamma        = 0.10d0

          avgvel(:)    = 100.0d0 *  &
     &                   (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(IN2O5)))**0.5d0

          where( sad(:) > 0.0d0 )
            sksts_n2o5(:) = 0.25d0 * gamma * avgvel * sad(:)
          elsewhere
            sksts_n2o5(:) = 0.0d0
          end where

          if (present(ptrop)) then
            where( pr(:) > ptrop ) sksts_n2o5(:) = 0.0d0
          end if

        END FUNCTION sksts_n2o5
!
!.... sktrs_n2o5 (temperature ,pressure ,sad_lbs ,tropp)
! -1-
!
!.... JPL 97-4
!
        FUNCTION sktrs_n2o5 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:)
          real*8, dimension(size(tk)) :: sktrs_n2o5
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk))

          pi = acos(-1.0d0)

!=======================================================================
!     N2O5 + tropospheric sulfate aerosol = 2 HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 3/30/99
!
!.... This is the pure water gamma tabulated in JPL 97-4, pp. 231
!.... and Note 15.
!
          gamma        = 0.05d0

          avgvel(:)    = 100.0d0 *  &
     &                   (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(IN2O5)))**0.5d0

          where( sad(:) > 0.0d0 )
            sktrs_n2o5(:) = 0.25d0 * gamma * avgvel * sad(:)
          elsewhere
            sktrs_n2o5(:) = 0.0d0
          end where

          if (present(ptrop)) then
            where( pr(:) < ptrop ) sktrs_n2o5(:) = 0.0d0
          end if

        END FUNCTION sktrs_n2o5
!
!.... sktrs_no3 (temperature ,pressure ,sad_lbs ,tropp)
! -2-
!
!.... JPL 97-4
!
        FUNCTION sktrs_no3 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:)
          real*8, dimension(size(tk)) :: sktrs_no3
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk))

          pi = acos(-1.0d0)

!=======================================================================
!     NO3 + tropospheric sulfate aerosol = HNO3
!=======================================================================
!
!.... First order rate constant
!.... PSC 3/30/99
!
!.... This is the pure water gamma tabulated in JPL 97-4, pp. 230
!.... and Note 13.
!
          gamma       = 2.0d-04

          avgvel(:)   = 100.0d0 *  &
     &                  (8.0d0 * 8.31448d0 * tk(:) * 1000.d0 /  &
     &                  (pi * mw(INO3)))**0.5d0

          where( sad(:) > 0.0d0 )
            sktrs_no3(:) = 0.25d0 * gamma * avgvel * sad(:)
          elsewhere
            sktrs_no3(:) = 0.0d0
          end where

          if (present(ptrop)) then
            where( pr(:) < ptrop ) sktrs_no3(:) = 0.0d0
          end if

        END FUNCTION sktrs_no3
!
!.... sksts_clono2 (temperature ,adcol ,pressure ,sad_lbs ,specarr( HCl,:) ,wa
!
!.... JPL 97-4
!
        FUNCTION sksts_clono2 ( tk,ad,pr,sad,hcl,h2o  &
     &                                          ,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,ad(:) ,pr(:)  &
     &     ,sad(:) ,h2o(:) ,hcl(:)
          real*8, dimension(size(tk)) :: sksts_clono2
          real*8  &
     &      adrop ,alpha ,ksur ,pi ,ro
          real*8  &
     &      adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gamma(size(tk)) ,gam0(size(tk)) ,gcalc(size(tk))  &
     &     ,gprob_hcl(size(tk)) ,gprob_tot(size(tk)) ,gsurf(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) ,prate(size(tk))  &
     &     ,tk_150(size(tk))

          pi            = acos(-1.0d0)

!=======================================================================
!     ClONO2 + stratospheric sulfate aerosol = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 3/30/99
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          alpha         = 1.0d0
          ksur          = 576.0d0
          ro            = 2000.0d0
          adrop         = 1.0d-05
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:)       = ((h2o(:) / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)

          where( hcl(:) <= 1.0d0 )
            phcl(:)     = ((1.0d0 / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)
          elsewhere
            phcl(:)     = ((hcl(:) / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:)       = 1013.25d0 * ph2o(:) /  &
     &                    10.0d0**  &
     &                    (9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))

          where( ah2o(:) > 1.1d0)
            ah2o(:)     = 1.1d0
          end where

          tk_150(:)     = tk(:)

          where( tk(:) < 150.0d0 )
            tk_150(:)   = 150.0d0
          endwhere

          hstar(:)      = exp((6250.0d0 / tk_150(:)) - 10.414d0) *  &
     &                    (ah2o(:)**3.49d0)

          gsurf(:)      = ah2o(:) * ksur * hstar(:) * phcl(:)

          prate(:)      = ro * hstar(:) * phcl(:) / ah2o(:)

          gam0(:)       = 1.18d-04 + (9.1d-03 * ah2o(:)) +  &
     &                    (0.5d0 * ah2o(:)**2.0d0)

          gcalc(:)      = gam0(:) * sqrt(1.0d0 + prate(:))

          adivl(:)      = adrop / (1.4d-06 * sqrt(1.0d0 / ah2o(:)))

          fterm(:)      = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                     (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &                    (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for ClONO2
!
          gprob_tot(:) = 1.0d0 / (1.0d0 /  &
     &                   (gsurf(:) + fterm(:) * gcalc(:)) +  &
     &                   1.0d0 / alpha)

          gprob_hcl(:) = gprob_tot(:) *  &
     &                   (gsurf(:) +  &
     &                    fterm(:) * gcalc(:) * prate(:) /  &
     &                    (1.0d0 + prate(:))) /  &
     &                   (gsurf(:) + fterm(:) * gcalc(:))

          gamma(:)     = gprob_tot(:) - gprob_hcl(:)

          avgvel(:)    = 100.0d0 *  &
     &                   (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(ICLONO2)))**0.5d0

          where( sad(:) > 0.0d0 )
            sksts_clono2(:) = 0.25d0 * gamma * avgvel * sad(:)
          elsewhere
            sksts_clono2(:) = 0.0d0
          end where

          if (present(ptrop)) then
            where( pr(:) > ptrop ) sksts_clono2(:) = 0.0d0
          end if

        END FUNCTION sksts_clono2
!
!.... sksts_brono2 (temperature ,pressure ,sad_lbs ,tropp)
!
!.... JPL 97-4
!
        FUNCTION sksts_brono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,pr(:) ,sad(:)
          real*8, dimension(size(tk)) :: sksts_brono2
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk))

          pi     = acos(-1.0d0)

!=======================================================================
!     BrONO2 + stratospheric sulfate aerosol = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... David Hanson, personal communication, May 13, 1997
!
          avgvel = 100.0d0 *  &
     &             (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &              (pi * mw(IBRONO2)))**0.5d0

          gamma  = 0.8d0

          where( sad(:) > 0.0d0 )
            sksts_brono2(:) = 0.25d0 * gamma * avgvel * sad(:)
          elsewhere
            sksts_brono2(:) = 0.0d0
          end where

          if ( present(ptrop) ) then
            where( pr(:) > ptrop ) sksts_brono2(:) = 0.0d0
          end if

        END FUNCTION sksts_brono2
!
!.... sksts_clono2_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(ClONO2,
!
!.... JPL 97-4
!
        FUNCTION sksts_clono2_hcl ( tk,ad,pr,sad,clono2  &
     &                                              ,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,ad(:) ,pr(:)  &
     &     ,sad(:) ,h2o(:) ,hcl(:) ,clono2(:)
          real*8, dimension(size(tk)) :: sksts_clono2_hcl
          real*8  &
     &      adrop ,alpha ,ksur ,minconc ,pi ,ro
          real*8  &
     &      adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gam0(size(tk)) ,gcalc(size(tk))  &
     &     ,gprob_hcl(size(tk)) ,gprob_tot(size(tk)) ,gsurf(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) ,prate(size(tk))  &
     &     ,tk_150(size(tk))

          pi            = acos(-1.0d0)

!=======================================================================
!     ClONO2 + HCl on stratospheric sulfate aerosol = Cl2 + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc       = 1.0d0
          alpha         = 1.0d0
          ksur          = 576.0d0
          ro            = 2000.0d0
          adrop         = 1.0d-05
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:)       = ((h2o(:) / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)

          where( hcl(:) <= minconc )
            phcl(:)     = ((1.0d0 / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)
          elsewhere
            phcl(:)     = ((hcl(:) / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:)       = 1013.25d0 * ph2o(:) /  &
     &                    10.0d0**  &
     &                    (9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))

          where( ah2o(:) > 1.1d0)
            ah2o(:)     = 1.1d0
          end where

          tk_150(:)     = tk(:)

          where( tk(:) < 150.0d0 )
            tk_150(:)   = 150.0d0
          endwhere

          hstar(:)      = exp((6250.0d0 / tk_150(:)) - 10.414d0) *  &
     &                    (ah2o(:)**3.49d0)

          gsurf(:)      = ah2o(:) * ksur * hstar(:) * phcl(:)

          prate(:)      = ro * hstar(:) * phcl(:) / ah2o(:)

          gam0(:)       = 1.18d-04 + (9.1d-03 * ah2o(:)) +  &
     &                    (0.5d0 * ah2o(:)**2.0d0)

          gcalc(:)      = gam0(:) * sqrt(1.0d0 + prate(:))

          adivl(:)      = adrop / (1.4d-06 * sqrt(1.0d0 / ah2o(:)))

          fterm(:)      = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                     (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &                    (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for ClONO2
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          gprob_tot(:)  = 1.0d0 / (1.0d0 /  &
     &                    (gsurf(:) + fterm(:) * gcalc(:)) +  &
     &                    1.0d0 / alpha)

          gprob_hcl(:)  = gprob_tot(:) *  &
     &                    (gsurf(:) +  &
     &                     fterm(:) * gcalc(:) * prate(:) /  &
     &                     (1.0d0 + prate(:))) /  &
     &                    (gsurf(:) + fterm(:) * gcalc(:))

          avgvel(:)     = 100.0d0 *  &
     &                    (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                    (pi * mw(ICLONO2)))**0.5d0

          where( hcl(:) > minconc )
            sksts_clono2_hcl(:) =  &
     &        0.25d0 * gprob_hcl(:) * avgvel(:) * sad(:) / hcl(:)
          elsewhere
            sksts_clono2_hcl(:) =  &
     &        0.25d0 * gprob_hcl(:) * avgvel(:) * sad(:)
          end where

          where( sad(:) < 0.0d0 )
            sksts_clono2_hcl(:) = 0.0d0
          end where

          if ( present(ptrop) ) then
            where( pr(:) > ptrop ) sksts_clono2_hcl(:) = 0.0d0
          end if

        END FUNCTION sksts_clono2_hcl
!
!.... sksts_hocl_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOCl,:) ,
!
!.... JPL 97-4
!
        FUNCTION sksts_hocl_hcl ( tk,ad,pr,sad,hocl  &
     &                                            ,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,ad(:) ,pr(:) ,sad(:)  &
     &     ,h2o(:) ,hcl(:) ,hocl(:)
          real*8, dimension(size(tk)) :: sksts_hocl_hcl
          real*8  &
     &      adrop ,alpha ,d1 ,minconc ,pi
          real*8  &
     &      adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,c1(size(tk)) ,c2(size(tk)) ,c3(size(tk)) ,conv(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gcalc(size(tk))  &
     &     ,gprob_tot(size(tk))  &
     &     ,hhuth(size(tk)) ,hm(size(tk))  &
     &     ,hsqrtd(size(tk)) ,hstar(size(tk)) ,hstar_hocl(size(tk))  &
     &     ,k(size(tk)) ,kii(size(tk))  &
     &     ,mterm(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk))  &
     &     ,rho(size(tk))  &
     &     ,tk_150(size(tk))  &
     &     ,wtper(size(tk))  &
     &     ,z(size(tk))

          pi            = acos(-1.0d0)

!=======================================================================
!     HOCl + HCl on stratospheric sulfate aerosol = Cl2 + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc         = 1.0d0
          adrop           = 1.0D-05
          alpha           = 1.0d0
          d1              = 9.0D-09
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:)         = (h2o(:) / ad(:)) * pr(:)

          z(:)            = log(ph2o(:))

          ph2o(:)         = ph2o(:) / 1013.25d0


          where( hcl(:) <= minconc )
            phcl(:)       = ((1.0d0 / ad(:)) * pr(:)) *  &
     &                      (1.0d0 / 1013.25d0)
          elsewhere
            phcl(:)       = ((hcl(:) / ad(:)) * pr(:)) *  &
     &                      (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:)       = 1013.25d0 * ph2o(:) /  &
     &                    10.0d0**  &
     &                    (9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))

          where( ah2o(:) > 1.1d0)
            ah2o(:)     = 1.1d0
          end where

          wtper(:)      = ((-14.0508d0 + 0.708928d0 * z(:)) *  &
     &                     tk(:) + 3578.6d0) /  &
     &                    (45.5374d0 + 1.55981d0 * z(:) -  &
     &                     0.197298d0 * tk(:))

          where( wtper(:) < 40.0d0 )
            wtper(:)    = 40.0d0
          end where

          where( wtper(:) > 80.0d0 )
            wtper(:)    = 80.0d0
          end where

          kii(:)        = exp(2.303d0 *  &
     &                        (6.08d0 - 1050.0d0 / tk(:) +  &
     &                         0.0747d0 * wtper(:)))

          tk_150(:)     = tk(:)

          where( tk(:) < 150.0d0 )
            tk_150(:)   = 150.0d0
          endwhere

          hstar(:)      = exp((6250.0d0 / tk_150(:)) - 10.414d0) *  &
     &                    (ah2o(:)**3.49d0)

          k(:)          = kii(:) * hstar(:) * phcl(:)

          adivl(:)      = adrop / sqrt(d1 / k(:))

          fterm(:)      = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                     (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &                    (1.0d0 / adivl(:))

          mterm(:)      = 10.196d0 * wtper(:) / (100.0d0 - wtper(:))

          c1(:)         = 123.64d0 - 5.6D-04 * tk(:)**2.0d0

          c2(:)         = -29.54d0 + 1.814D-04 * tk(:)**2.0d0

          c3(:)         = 2.243d0 - 1.487D-03 * tk(:) +  &
     &                    1.324D-05 * tk(:)**2.0d0

          rho(:)        = 1000.0d0 + c1(:) * mterm(:) +  &
     &                               c2(:) * mterm(:)**1.5d0 +  &
     &                               c3(:) * mterm(:)**2.0d0

          conv(:)       = (rho(:) / 1000.0d0) /  &
     &                    (1.0d0 + mterm(:) * 0.09808d0)

          hhuth(:)      = exp(6.4946d0 -  &
     &                        mterm(:) *  &
     &                        (-0.04107d0 + 54.56d0 / tk(:)) -  &
     &                        5862.0d0 *  &
     &                        (1.0d0 / 298.15d0 - 1.0d0 / tk(:)))

          hm(:)         = hhuth(:) * conv(:)

          hstar_hocl(:) = hm(:) *  &
     &                    (1.0d0 + 1.052d0 *  &
     &                     exp(0.273d0 * (wtper(:) - 65.66d0)))

          hsqrtd(:)     = hstar_hocl(:) * sqrt(d1)

          gcalc(:)      = 2.2548D-05 * hsqrtd(:) *  &
     &                    sqrt(tk(:) * mw(IHOCL) * k(:))

!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm(:) > 0.0d0 )
            gprob_tot(:) = 1.0d0 / (1.0d0 /  &
     &                     (fterm(:) * gcalc(:)) +  &
     &                     1.0d0 / alpha)
          elsewhere
            gprob_tot(:) = 0.0d0
          end where

          avgvel(:)      = 100.0d0 *  &
     &                     (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                     (pi * mw(IHOCL)))**0.5d0

          where( hcl(:) > minconc )
            sksts_hocl_hcl(:) =  &
     &        0.25d0 * gprob_tot(:) * avgvel(:) * sad(:) / hcl(:)
          elsewhere
            sksts_hocl_hcl(:) =  &
     &        0.25d0 * gprob_tot(:) * avgvel(:) * sad(:)
          end where

          where( sad(:) < 0.0d0 )
            sksts_hocl_hcl(:) = 0.0d0
          end where

          if ( present(ptrop) ) then
            where( pr(:) > ptrop ) sksts_hocl_hcl(:) = 0.0d0
          end if

        END FUNCTION sksts_hocl_hcl
!
!.... sksts_hobr_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOBr,:) ,
!
!.... JPL 97-4
!
        FUNCTION sksts_hobr_hcl ( tk,ad,pr,sad,hobr  &
     &                                            ,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,ad(:) ,pr(:) ,sad(:)  &
     &     ,h2o(:) ,hcl(:) ,hobr(:)
          real*8, dimension(size(tk)) :: sksts_hobr_hcl
          real*8  &
     &      adrop ,alpha ,d1 ,hsqrtd ,kii ,minconc ,pi
          real*8  &
     &      adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gcalc(size(tk))  &
     &     ,gprob_tot(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,k(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk))  &
     &     ,tk_150(size(tk))

          pi            = acos(-1.0d0)

!=======================================================================
!     HOBr + HCl on stratospheric sulfate aerosol = BrCl + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc       = 1.0d0
          adrop         = 1.0D-05
          alpha         = 1.0d0
          d1            = 1.2D-08
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:)       = (h2o(:) / ad(:)) * pr(:)

          ph2o(:)       = ph2o(:) / 1013.25d0

          where( hcl(:) <= minconc )
            phcl(:)     = ((1.0d0 / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)
          elsewhere
            phcl(:)     = ((hcl(:) / ad(:)) * pr(:)) *  &
     &                    (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:)       = 1013.25d0 * ph2o(:) /  &
     &                    10.0d0**  &
     &                    (9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))

          where( ah2o(:) > 1.1d0)
            ah2o(:)     = 1.1d0
          end where

          tk_150(:)     = tk(:)

          where( tk(:) < 150.0d0 )
            tk_150(:)   = 150.0d0
          endwhere

          hstar(:)      = exp((6250.0d0 / tk_150(:)) - 10.414d0) *  &
     &                    (ah2o(:)**3.49d0)

          kii           = 1.0D+05

          k(:)          = kii * hstar(:) * phcl(:)

          hsqrtd        = 110.0d0

          gcalc(:)      = 2.2548D-05 * hsqrtd *  &
     &                    sqrt(tk(:) * mw(IHOBR) * k(:))

          adivl(:)      = adrop / sqrt(d1 / k(:))

          fterm(:)      = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                     (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &                    (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm(:) > 0.0d0 )
            gprob_tot(:) = 1.0d0 / (1.0d0 /  &
     &                     (fterm(:) * gcalc(:)) +  &
     &                     1.0d0 / alpha)
          elsewhere
            gprob_tot(:) = 0.0d0
          end where

          avgvel(:)      = 100.0d0 *  &
     &                     (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                     (pi * mw(IHOBR)))**0.5d0

          where( hcl(:) > minconc )
            sksts_hobr_hcl(:) =  &
     &        0.25d0 * gprob_tot(:) * avgvel(:) * sad(:) / hcl(:)
          elsewhere
            sksts_hobr_hcl(:) =  &
     &        0.25d0 * gprob_tot(:) * avgvel(:) * sad(:)
          end where

          where( sad(:) < 0.0d0 )
            sksts_hobr_hcl(:) = 0.0d0
          end where

          if ( present(ptrop) ) then
            where( pr(:) > ptrop ) sksts_hobr_hcl(:) = 0.0d0
          end if

        END FUNCTION sksts_hobr_hcl
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
      END
