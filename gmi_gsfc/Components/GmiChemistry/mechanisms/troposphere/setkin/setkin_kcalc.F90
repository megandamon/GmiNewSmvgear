!=======================================================================
!
! $Id: setkin_kcalc.F90,v 1.6 2013-07-31 16:35:49 ssteenro Exp $
!
! ROUTINE
!   kcalc - GMIMOD 3D model (setkin_kcalc.F)
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
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!=======================================================================
      subroutine kcalc( npres0,sadcol,sadcol2,pressure,tropp, &
     &  temperature,lwc,adcol,specarr,rcarr,radA,FRH)

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
      REAL*8,  INTENT (IN)  :: sadcol2     (NSADaer+NSADdust,npres0)
      REAL*8,  INTENT (IN)  :: specarr     (NMF   ,npres0)
      REAL*8,  INTENT (IN)  :: FRH         (       npres0)
      REAL*8,  INTENT (IN)  :: radA        (NSADaer+NSADdust,npres0)

      REAL*8,  INTENT (OUT) :: rcarr       (NUM_K ,npres0)

!.... Local variable declarations

      INTEGER               :: naltmax

      real*8 &
     &  nitrogen (npres0) &
     & ,oxygen   (npres0) &
     & ,water    (npres0)

! The sad_* variables are not used in the Tropospheric Mechanism.
      real*8 &
     &  sad_ice  (npres0) &
     & ,sad_lbs  (npres0) &
     & ,sad_nat  (npres0) &
     & ,sad_soot (npres0) &
     & ,sad_sts  (npres0)

      real*8 mw(NSP)

      mw(:) = mw_data(:)

      naltmax     = npres0
!===================================================
! Variables needed for gas/heterogeneous chemistry.
! adcol = air density (molec/cm3)
! FRH = relative humidity fraction (0-1)
! radA = effective radius of aerosol (cm)
! sadcol = surface area of aerosol/volume of air (cm2/cm3)
!===================================================

!....          Define molecular nitrogen and oxygen number densities
!
      nitrogen(:) = adcol(:) * MXRN2
      oxygen(:)   = adcol(:) * MXRO2
      water(:)    = specarr(83 ,:)

      sad_lbs(:)  = sadcol(1, :)
      sad_sts(:)  = sadcol(2, :)
      sad_nat(:)  = sadcol(3, :)
      sad_ice(:)  = sadcol(4, :)
      sad_soot(:) = sadcol(5, :)
!
!....          Start thermal rate constants
!
!....           NO + O3 = NO2
!
      rcarr(1,:) = skarr(  3.000D-12 ,1500.0D+00 ,temperature)
!
!....           O3 + OH = HO2
!
      rcarr(2,:) = skarr(  1.700D-12 ,940.0D+00 ,temperature)
!
!....           HO2 + O3 = OH
!
      rcarr(3,:) = skarr(  1.000D-14 ,490.0D+00 ,temperature)
!
!....           NO2 + O3 = NO3
!
      rcarr(4,:) = skarr(  1.200D-13 ,2450.0D+00 ,temperature)
!
!....           MO2 + O3 = CH2O + HO2
!
      rcarr(5,:) = skarr(  2.900D-16 ,1000.0D+00 ,temperature)
!
!....           OH + OH = H2O + O3
!
      rcarr(6,:) = skarr(  1.800D-12 ,0.0D+00 ,temperature)
!
!....           OH + OH = H2O2
!
      rcarr(7,:) = sktroe(  6.900D-31 ,1.00D0 & 
     &                     , 2.600D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HO2 + OH = H2O
!
      rcarr(8,:) = skarr(  4.800D-11 ,-250.0D+00 ,temperature)
!
!....           H2O2 + OH = H2O + HO2
!
      rcarr(9,:) = skarr(  1.800D-12 ,0.0D+00 ,temperature)
!
!....           HO2 + NO = NO2 + OH
!
      rcarr(10,:) = skarr(  3.500D-12 ,-250.0D+00 ,temperature)
!
!....           HO2 + HO2 = H2O2
!
      rcarr(11,:) = skho2dis (temperature ,adcol ,water)
!
!....           H2 + OH = H2O + HO2
!
      rcarr(12,:) = skarr(  2.800D-12 ,1800.0D+00 ,temperature)
!
!....           CO + OH = HO2
!
      rcarr(13,:) = skcooh (temperature ,adcol)
!
!....           CH4 + OH = H2O + MO2
!
      rcarr(14,:) = skohch4 (temperature)
!
!....           MO2 + NO = CH2O + HO2 + NO2
!
      rcarr(15,:) = skarr(  2.800D-12 ,-300.0D+00 ,temperature)
!
!....           HO2 + MO2 = MP
!
      rcarr(16,:) = skarr(  4.100D-13 ,-750.0D+00 ,temperature)
!
!....           MO2 + MO2 = CH2O + MOH
!
      rcarr(17,:) = skmo2dis_1 (temperature)
!
!....           MO2 + MO2 = 2 CH2O + 2 HO2
!
      rcarr(18,:) = skmo2dis_2 (temperature)
!
!....           MP + OH = H2O + MO2
!
      rcarr(19,:) = skarr(  2.700D-12 ,-200.0D+00 ,temperature)
!
!....           MP + OH = CH2O + H2O + OH
!
      rcarr(20,:) = skarr(  1.100D-12 ,-200.0D+00 ,temperature)
!
!....           CH2O + OH = CO + H2O + HO2
!
      rcarr(21,:) = skarr(  5.500D-12 ,-125.0D+00 ,temperature)
!
!....           NO2 + OH = HNO3
!
      rcarr(22,:) = sktroe(  1.800D-30 ,3.00D0 & 
     &                     , 2.800D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO3 + OH = H2O + NO3
!
      rcarr(23,:) = skohhno3 (temperature ,adcol)
!
!....           NO + OH = HNO2
!
      rcarr(24,:) = sktroe(  7.000D-31 ,2.60D0 & 
     &                     , 3.600D-11 ,0.10D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO2 + OH = H2O + NO2
!
      rcarr(25,:) = skarr(  1.800D-11 ,390.0D+00 ,temperature)
!
!....           HO2 + NO2 = HNO4
!
      rcarr(26,:) = sktroe(  2.000D-31 ,3.40D0 & 
     &                     , 2.900D-12 ,1.10D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO4 = HO2 + NO2
!
      rcarr(27,:) = sktroe(  9.520D-05 ,3.4D0 &
     &                     , 1.380D+15 ,1.10D0 ,10900.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HNO4 + OH = H2O + NO2
!
      rcarr(28,:) = skarr(  1.300D-12 ,-380.0D+00 ,temperature)
!
!....           HO2 + NO3 = NO2 + OH
!
      rcarr(29,:) = skarr(  3.500D-12 ,0.0D+00 ,temperature)
!
!....           NO + NO3 = 2 NO2
!
      rcarr(30,:) = skarr(  1.500D-11 ,-170.0D+00 ,temperature)
!
!....           NO3 + OH = HO2 + NO2
!
      rcarr(31,:) = skarr(  2.200D-11 ,0.0D+00 ,temperature)
!
!....           NO2 + NO3 = N2O5
!
      rcarr(32,:) = sktroe(  2.000D-30 ,4.40D0 & 
     &                     , 1.400D-12 ,0.70D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           N2O5 = NO2 + NO3
!
      rcarr(33,:) = sktroe(  7.400D-04 ,4.4D0 &
     &                     , 5.190D+14 ,0.70D0 ,11000.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HCOOH + OH = H2O + HO2
!
      rcarr(34,:) = skarr(  4.000D-13 ,0.0D+00 ,temperature)
!
!....           MOH + OH = CH2O + HO2
!
      rcarr(35,:) = skarr(  2.900D-12 ,345.0D+00 ,temperature)
!
!....           NO2 + NO3 = NO + NO2
!
      rcarr(36,:) = skarr(  4.500D-14 ,1260.0D+00 ,temperature)
!
!....           CH2O + NO3 = CO + HNO3 + HO2
!
      rcarr(37,:) = skarr(  5.800D-16 ,0.0D+00 ,temperature)
!
!....           ALD2 + OH = H2O + MCO3
!
      rcarr(38,:) = skarr(  5.600D-12 ,-270.0D+00 ,temperature)
!
!....           ALD2 + NO3 = HNO3 + MCO3
!
      rcarr(39,:) = skarr(  1.400D-12 ,1900.0D+00 ,temperature)
!
!....           MCO3 + NO2 = PAN
!
      rcarr(40,:) = sktroe(  9.700D-29 ,5.60D0 & 
     &                     , 9.300D-12 ,1.50D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           PAN = MCO3 + NO2
!
      rcarr(41,:) = skpanan (temperature,adcol)
!
!....           MCO3 + NO = MO2 + NO2
!
      rcarr(42,:) = skarr(  8.100D-12 ,-270.0D+00 ,temperature)
!
!....           C2H6 + OH = ETO2 + H2O
!
      rcarr(43,:) = skarr(  8.700D-12 ,1070.0D+00 ,temperature)
!
!....           ETO2 + NO = ALD2 + HO2 + NO2
!
      rcarr(44,:) = skarr(  2.600D-12 ,-365.0D+00 ,temperature)
!
!....           C3H8 + OH = B3O2
!
      rcarr(45,:) = skarr(  6.300D-12 ,580.0D+00 ,temperature)
!
!....           C3H8 + OH = A3O2
!
      rcarr(46,:) = skarr(  6.300D-12 ,1050.0D+00 ,temperature)
!
!....           A3O2 + NO = HO2 + NO2 + RCHO
!
      rcarr(47,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           NO + PO2 = ALD2 + CH2O + HO2 + NO2
!
      rcarr(48,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           ALK4 + OH = R4O2
!
      rcarr(49,:) = skarr(  9.100D-12 ,405.0D+00 ,temperature)
!
!....           NO + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O
!
      rcarr(50,:) = skro2noabs_1 (temperature ,adcol)
!
!....           NO + R4O2 = R4N2
!
      rcarr(51,:) = skro2noadd_1 (temperature ,adcol)
!
!....           NO + R4N1 =  0.75 ALD2 +  0.39 CH2O + 2 NO2 +  0.30 R4O2 +  
!
      rcarr(52,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           ATO2 + NO =  0.19 CH2O +  0.77 HO2 +  0.19 MCO3 +  0.77 MGLY
!
      rcarr(53,:) = skarr(  2.900D-12 ,-300.0D+00 ,temperature)
!
!....           KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2 +  0.07 R4N2
!
      rcarr(54,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           NO + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR
!
      rcarr(55,:) = skro2noabs_2 (temperature ,adcol)
!
!....           NO + RIO2 = HNO3
!
      rcarr(56,:) = skro2noadd_2 (temperature ,adcol)
!
!....           NO + RIO1 =  0.75 CH2O + HO2 + IALD + NO2
!
      rcarr(57,:) = skro2noabs_2 (temperature ,adcol)
!
!....           NO + RIO1 = HNO3
!
      rcarr(58,:) = skro2noadd_2 (temperature ,adcol)
!
!....           IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.24 GLYC +  0.17 GLYX 
!
      rcarr(59,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           NO + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3
!
      rcarr(60,:) = skro2noabs_3 (temperature ,adcol)
!
!....           NO + VRO2 = HNO3
!
      rcarr(61,:) = skro2noadd_3 (temperature ,adcol)
!
!....           MRO2 + NO =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 
!
      rcarr(62,:) = skro2noabs_3 (temperature ,adcol)
!
!....           MRO2 + NO = HNO3
!
      rcarr(63,:) = skro2noadd_3 (temperature ,adcol)
!
!....           MVN2 + NO =  0.30 CH2O +  0.60 GLYC +  0.10 HNO3 +  0.30 HO2
!
      rcarr(64,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           MAN2 + NO = CH2O + MGLY + 2 NO2
!
      rcarr(65,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           B3O2 + NO = ACET + HO2 + NO2
!
      rcarr(66,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR
!
      rcarr(67,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           NO + PRN1 = ALD2 + CH2O + 2 NO2
!
      rcarr(68,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           ALK4 + NO3 = HNO3 + R4O2
!
      rcarr(69,:) = skarr(  2.800D-12 ,3280.0D+00 ,temperature)
!
!....           OH + R4N2 = H2O + R4N1
!
      rcarr(70,:) = skarr(  1.300D-12 ,0.0D+00 ,temperature)
!
!....           ACTA + OH = H2O + MO2
!
      rcarr(71,:) = skarr(  4.000D-13 ,-200.0D+00 ,temperature)
!
!....           OH + RCHO = H2O + RCO3
!
      rcarr(72,:) = skarr(  4.900D-12 ,-405.0D+00 ,temperature)
!
!....           NO2 + RCO3 = PPN
!
      rcarr(73,:) = sktroe(  9.000D-28 ,8.90D0 & 
     &                     , 7.700D-12 ,0.20D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           GCO3 + NO2 = GPAN
!
      rcarr(74,:) = sktroe(  9.000D-28 ,8.90D0 & 
     &                     , 7.700D-12 ,0.20D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           MAO3 + NO2 = PMN
!
      rcarr(75,:) = sktroe(  9.000D-28 ,8.90D0 & 
     &                     , 7.700D-12 ,0.20D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           PPN = NO2 + RCO3
!
      rcarr(76,:) = skppndecomp (temperature,adcol)
!
!....           GPAN = GCO3 + NO2
!
      rcarr(77,:) = skppndecomp (temperature,adcol)
!
!....           NO + RCO3 = ETO2 + NO2
!
      rcarr(78,:) = skarr(  6.700D-12 ,-340.0D+00 ,temperature)
!
!....           GCO3 + NO = CH2O + HO2 + NO2
!
      rcarr(79,:) = skarr(  6.700D-12 ,-340.0D+00 ,temperature)
!
!....           MAO3 + NO = 4 CH2O + HO2 + NO2
!
      rcarr(80,:) = skarr(  6.700D-12 ,-340.0D+00 ,temperature)
!
!....           NO3 + RCHO = HNO3 + RCO3
!
      rcarr(81,:) = skarr(  1.820D-12 ,1680.0D+00 ,temperature)
!
!....           ACET + OH = ATO2 + H2O
!
      rcarr(82,:) = skacetoh (temperature)
!
!....           A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.
!
      rcarr(83,:) = skarr(  5.920D-13 ,0.0D+00 ,temperature)
!
!....           MO2 + PO2 =  0.50 ALD2 +  1.25 CH2O +  0.16 HAC + HO2 +  0.2
!
      rcarr(84,:) = skarr(  5.920D-13 ,0.0D+00 ,temperature)
!
!....           HO2 + R4O2 = R4P
!
      rcarr(85,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + R4N1 = R4N2
!
      rcarr(86,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           ATO2 + HO2 = MCO3 + MO2
!
      rcarr(87,:) = skarr(  8.600D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + KO2 = MGLY + MO2
!
      rcarr(88,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + RIO2 = RIP
!
      rcarr(89,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + RIO1 = RIP
!
      rcarr(90,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + IAO2 = IAP
!
      rcarr(91,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + ISN1 = ISNP
!
      rcarr(92,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + VRO2 = VRP
!
      rcarr(93,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + MRO2 = MRP
!
      rcarr(94,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + MVN2 = ISNP
!
      rcarr(95,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + MAN2 = ISNP
!
      rcarr(96,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           B3O2 + HO2 = RB3P
!
      rcarr(97,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + INO2 = INPN
!
      rcarr(98,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + PRN1 = PRPN
!
      rcarr(99,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           MEK + OH = H2O + KO2
!
      rcarr(100,:) = skohmek (temperature)
!
!....           ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O +  0.25 EOH + HO2 +  0.
!
      rcarr(101,:) = skarr(  3.000D-13 ,0.0D+00 ,temperature)
!
!....           MEK + NO3 = HNO3 + KO2
!
      rcarr(102,:) = skarr(  8.000D-16 ,0.0D+00 ,temperature)
!
!....           MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3
!
      rcarr(103,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH
!
      rcarr(104,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           ATO2 + MO2 =  0.85 CH2O +  0.90 HO2 +  0.10 MCO3 +  0.25 MEK
!
      rcarr(105,:) = skarr(  7.500D-13 ,-500.0D+00 ,temperature)
!
!....           KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3
!
      rcarr(106,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.06 IALD +  0.14 MAC
!
      rcarr(107,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + RIO1 =  1.13 CH2O + HO2 +  0.50 IALD +  0.25 MEK +  0.
!
      rcarr(108,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           IAO2 + MO2 =  0.95 CH2O +  0.15 CO +  0.13 GLYC +  0.09 GLYX
!
      rcarr(109,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           ISN1 + MO2 =  0.75 CH2O +  0.50 GLYC +  0.50 HAC +  0.50 HO2
!
      rcarr(110,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + VRO2 =  0.89 CH2O +  0.36 GLYC +  0.64 HO2 +  0.36 MCO
!
      rcarr(111,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + MRO2 =  0.84 CH2O +  0.42 CO +  0.42 HAC + HO2 +  0.25
!
      rcarr(112,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + MVN2 =  1.25 CH2O +  0.75 HO2 +  0.25 MCO3 +  0.25 MGL
!
      rcarr(113,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MAN2 + MO2 =  1.25 CH2O +  0.50 HO2 +  0.50 MGLY +  0.25 MOH
!
      rcarr(114,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.
!
      rcarr(115,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MAC
!
      rcarr(116,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH
!
      rcarr(117,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           EOH + OH = ALD2 + HO2
!
      rcarr(118,:) = skarr(  6.900D-12 ,230.0D+00 ,temperature)
!
!....           OH + ROH = HO2 + RCHO
!
      rcarr(119,:) = skarr(  3.300D-12 ,-180.0D+00 ,temperature)
!
!....           ETO2 + ETO2 = 2 ALD2 + 2 HO2
!
      rcarr(120,:) = skarr(  4.100D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + ETO2 = ALD2 + EOH
!
      rcarr(121,:) = skarr(  2.720D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + HO2 = ETP
!
      rcarr(122,:) = skarr(  7.500D-13 ,-700.0D+00 ,temperature)
!
!....           A3O2 + HO2 = RA3P
!
      rcarr(123,:) = skarr(  7.500D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + PO2 = PP
!
      rcarr(124,:) = skarr(  7.500D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + MCO3 = ACTA + O3
!
      rcarr(125,:) = skho2mco3_1 (temperature)
!
!....           HO2 + MCO3 = MAP
!
      rcarr(126,:) = skho2mco3_2 (temperature)
!
!....           HO2 + RCO3 =  0.30 O3 +  0.30 RCOOH +  0.70 RP
!
      rcarr(127,:) = skarr(  4.300D-13 ,-1040.0D+00 ,temperature)
!
!....           GCO3 + HO2 =  0.70 GP +  0.30 O3 +  0.30 RCOOH
!
      rcarr(128,:) = skarr(  4.300D-13 ,-1040.0D+00 ,temperature)
!
!....           HO2 + MAO3 =  0.70 MAOP +  0.30 O3 +  0.30 RCOOH
!
      rcarr(129,:) = skarr(  4.300D-13 ,-1040.0D+00 ,temperature)
!
!....           OH + PRPE = PO2
!
      rcarr(130,:) = sktroe(  8.000D-27 ,3.50D0 & 
     &                     , 3.000D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           O3 + PRPE =  0.50 ALD2 +  0.54 CH2O +  0.42 CO +  0.06 H2 + 
!
      rcarr(131,:) = skarr(  6.500D-15 ,1880.0D+00 ,temperature)
!
!....           PMN = MAO3 + NO2
!
      rcarr(132,:) = skppndecomp (temperature,adcol)
!
!....           OH + PMN =  2.23 CH2O +  0.59 HAC + 2 HO2 + NO2
!
      rcarr(133,:) = skarr(  3.200D-11 ,0.0D+00 ,temperature)
!
!....           O3 + PMN =  0.60 CH2O + HO2 + NO2
!
      rcarr(134,:) = skarr(  8.200D-18 ,0.0D+00 ,temperature)
!
!....           GLYC + OH =  0.80 GCO3 +  0.20 GLYX +  0.20 HO2
!
      rcarr(135,:) = skarr(  1.000D-11 ,0.0D+00 ,temperature)
!
!....           NO3 + PRPE = PRN1
!
      rcarr(136,:) = skarr(  4.590D-13 ,1156.0D+00 ,temperature)
!
!....           GLYX + OH = 2 CO + HO2
!
      rcarr(137,:) = skohglyx (temperature ,oxygen)
!
!....           MGLY + OH = CO + MCO3
!
      rcarr(138,:) = skarr(  1.700D-11 ,0.0D+00 ,temperature)
!
!....           GLYX + NO3 = 2 CO + HNO3 + HO2
!
      rcarr(139,:) = skno3glyx (temperature ,oxygen)
!
!....           MGLY + NO3 = CO + HNO3 + MCO3
!
      rcarr(140,:) = skarr(  1.400D-12 ,1860.0D+00 ,temperature)
!
!....           ISOP + OH = RIO2
!
      rcarr(141,:) = skarr(  2.700D-11 ,-390.0D+00 ,temperature)
!
!....           MVK + OH = VRO2
!
      rcarr(142,:) = skarr(  4.130D-12 ,-452.0D+00 ,temperature)
!
!....           MACR + OH =  0.50 MAO3 +  0.50 MRO2
!
      rcarr(143,:) = skarr(  1.860D-11 ,-175.0D+00 ,temperature)
!
!....           HAC + OH = HO2 + MGLY
!
      rcarr(144,:) = skarr(  3.000D-12 ,0.0D+00 ,temperature)
!
!....           A3O2 + MCO3 = HO2 + MO2 + RCHO
!
      rcarr(145,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + PO2 = ALD2 + CH2O + HO2 + MO2
!
      rcarr(146,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           A3O2 + MCO3 = ACTA + RCHO
!
      rcarr(147,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + PO2 = ACTA +  0.65 HAC +  0.35 RCHO
!
      rcarr(148,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.39 MACR +
!
      rcarr(149,:) = skarr(  1.050D-14 ,2000.0D+00 ,temperature)
!
!....           MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 + 
!
      rcarr(150,:) = skarr(  4.000D-15 ,2000.0D+00 ,temperature)
!
!....           MACR + O3 =  0.70 CH2O +  0.20 CO +  0.28 HO2 +  0.80 MGLY +
!
      rcarr(151,:) = skarr(  4.400D-15 ,2500.0D+00 ,temperature)
!
!....           ISOP + NO3 = INO2
!
      rcarr(152,:) = skarr(  3.030D-12 ,446.0D+00 ,temperature)
!
!....           MVK + NO3 = MVN2
!
      rcarr(153,:) = skarr(  2.000D-14 ,0.0D+00 ,temperature)
!
!....           MACR + NO3 = MAN2
!
      rcarr(154,:) = skarr(  6.700D-15 ,0.0D+00 ,temperature)
!
!....           MACR + NO3 = HNO3 + MAO3
!
      rcarr(155,:) = skarr(  3.300D-15 ,0.0D+00 ,temperature)
!
!....           MO2 + RCO3 = CH2O + ETO2 + HO2
!
      rcarr(156,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           GCO3 + MO2 = 2 CH2O + 2 HO2
!
      rcarr(157,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MAO3 + MO2 = 2 CH2O + HO2 + MCO3
!
      rcarr(158,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + RCO3 = CH2O + RCOOH
!
      rcarr(159,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           GCO3 + MO2 = CH2O + RCOOH
!
      rcarr(160,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MAO3 + MO2 = CH2O + RCOOH
!
      rcarr(161,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           INPN + OH = INO2
!
      rcarr(162,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + PRPN = PRN1
!
      rcarr(163,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH
!
      rcarr(164,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      rcarr(165,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RB3P =  0.50 B3O2 +  0.50 OH +  0.50 RCHO
!
      rcarr(166,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + R4P =  0.50 OH +  0.50 R4O2 +  0.50 RCHO
!
      rcarr(167,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3
!
      rcarr(168,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + PP =  0.50 OH +  0.50 PO2 +  0.50 RCHO
!
      rcarr(169,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           GP + OH =  0.50 CH2O +  0.50 GCO3 +  0.50 OH
!
      rcarr(170,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2
!
      rcarr(171,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(172,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           ISNP + OH =  0.50 ISN1 +  0.50 NO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(173,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           OH + VRP =  0.50 OH +  0.50 RCHO +  0.50 VRO2
!
      rcarr(174,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           MRP + OH =  0.50 MRO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(175,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           MAOP + OH =  0.50 MAO3 +  0.50 OH +  0.50 RCHO
!
      rcarr(176,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH
!
      rcarr(177,:) = skarr(  3.800D-12 ,-200.0D+00 ,temperature)
!
!....           C2H6 + NO3 = ETO2 + HNO3
!
      rcarr(178,:) = skarr(  1.400D-18 ,0.0D+00 ,temperature)
!
!....           IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3
!
      rcarr(179,:) = skarr(  3.700D-11 ,0.0D+00 ,temperature)
!
!....           IALD + O3 =  0.12 CH2O +  0.28 GLYC +  0.20 GLYX +  0.20 HAC
!
      rcarr(180,:) = skarr(  6.160D-15 ,1814.0D+00 ,temperature)
!
!....           MCO3 + MCO3 = 2 MO2
!
      rcarr(181,:) = skarr(  2.900D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MO2 = CH2O + HO2 + MO2
!
      rcarr(182,:) = skarr(  1.800D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MO2 = ACTA + CH2O
!
      rcarr(183,:) = skarr(  2.000D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B
!
      rcarr(184,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 +  0.80 MG
!
      rcarr(185,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           KO2 + MCO3 = ALD2 + MCO3 + MO2
!
      rcarr(186,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MA
!
      rcarr(187,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO1 =  0.75 CH2O + HO2 + IALD + MO2
!
      rcarr(188,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           IAO2 + MCO3 =  0.40 CH2O +  0.29 CO +  0.26 GLYC +  0.18 GLY
!
      rcarr(189,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ISN1 + MCO3 = GLYC + HAC + MO2 + NO2
!
      rcarr(190,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MC
!
      rcarr(191,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MRO2 =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.1
!
      rcarr(192,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           B3O2 + MCO3 = ACET + HO2 + MO2
!
      rcarr(193,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + R4N1 =  0.75 ALD2 +  0.39 CH2O + MO2 + NO2 +  0.30 R4
!
      rcarr(194,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MVN2 = CH2O +  0.50 HO2 +  0.50 MCO3 +  0.50 MGLY + M
!
      rcarr(195,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MAN2 + MCO3 = CH2O + MGLY + MO2 + NO2
!
      rcarr(196,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MA
!
      rcarr(197,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + PRN1 = ALD2 + CH2O + MO2 + NO2
!
      rcarr(198,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + R4O2 = ACTA + MEK
!
      rcarr(199,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ATO2 + MCO3 = ACTA + MEK
!
      rcarr(200,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           KO2 + MCO3 = ACTA + MEK
!
      rcarr(201,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO2 = ACTA + MEK
!
      rcarr(202,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO1 = ACTA + MEK
!
      rcarr(203,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           IAO2 + MCO3 = ACTA + MEK
!
      rcarr(204,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + VRO2 = ACTA + MEK
!
      rcarr(205,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MRO2 = ACTA + MEK
!
      rcarr(206,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + R4N1 = ACTA + NO2 + RCHO
!
      rcarr(207,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ISN1 + MCO3 = ACTA + NO2 + RCHO
!
      rcarr(208,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MVN2 = ACTA + NO2 + RCHO
!
      rcarr(209,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MAN2 + MCO3 = ACTA + NO2 + RCHO
!
      rcarr(210,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           INO2 + MCO3 = ACTA + NO2 + RCHO
!
      rcarr(211,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + PRN1 = ACTA + NO2 + RCHO
!
      rcarr(212,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           B3O2 + MCO3 = ACET + ACTA
!
      rcarr(213,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ETO2 + MCO3 = ALD2 + HO2 + MO2
!
      rcarr(214,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ETO2 + MCO3 = ACTA + ALD2
!
      rcarr(215,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RCO3 = ETO2 + MO2
!
      rcarr(216,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           GCO3 + MCO3 = CH2O + HO2 + MO2
!
      rcarr(217,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           MAO3 + MCO3 = CH2O + MCO3 + MO2
!
      rcarr(218,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           NO3 + NO3 = 2 NO2
!
      rcarr(219,:) = skarr(  8.500D-13 ,2450.0D+00 ,temperature)
!
!....           HO2 =  0.50 H2O2
!
      rcarr(220,:) = sktrs_ho2 (temperature, & 
     &           sadcol2,adcol,radA,NSADaer,NSADdust,tropp,pressure)
!
!....           NO2 =  0.50 HNO2 +  0.50 HNO3
!
      rcarr(221,:) = sktrs_no2 (temperature, & 
     &            sadcol2, adcol, radA, NSADaer,NSADdust,tropp & 
     &           , pressure)
!
!....           NO3 = HNO3
!
      rcarr(222,:) = sktrs_no3 (temperature, & 
     &            sadcol2, adcol, radA,NSADaer,NSADdust,tropp & 
     &           , pressure)
!
!....           N2O5 = 2 HNO3
!
      rcarr(223,:) = sktrs_n2o5 (temperature, & 
     &           sadcol2,adcol,radA,FRH,NSADaer,NSADdust, tropp & 
     &           , pressure)
!
!....          End thermal rate constants
!
      CONTAINS
        FUNCTION skarr (af,ae,tk)
          real*8 &
     &      af ,ae ,tk(:)
          real*8, dimension(size(tk)) :: skarr
          skarr(:) = af * exp(-ae / tk(:))
        END FUNCTION skarr
        FUNCTION sklp (af,npwr,tk,ad)
          real*8 &
     &      af ,npwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: sklp
          sklp(:) = ad(:) * af * (300.0d0/tk(:))**npwr
        END FUNCTION sklp
        FUNCTION skhp (ai,mpwr,tk)
          real*8 &
     &      ai ,mpwr ,tk(:)
          real*8, dimension(size(tk)) :: skhp
          skhp(:) = ai * (300.0d0/tk(:))**mpwr
        END FUNCTION skhp
        FUNCTION skfo (af,npwr,ai,mpwr,tk,ad)
          real*8 &
     &      af ,npwr ,ai ,mpwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skfo
          skfo(:) = sklp(af,npwr,tk,ad) / skhp(ai,mpwr,tk)
        END FUNCTION skfo
        FUNCTION skterlp (af,npwr,ae,tk,ad)
          real*8 &
     &      af ,npwr ,ae ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skterlp
          skterlp(:) = sklp(af,npwr,tk,ad) * exp(-ae / tk(:))
        END FUNCTION skterlp
        FUNCTION sktroe (af,npwr,ai,mpwr,ae,tk,ad,fc)
          real*8 &
     &      af ,npwr ,ai ,mpwr ,ae ,tk(:) ,ad(:)
          real*8, OPTIONAL :: fc
          real*8, dimension(size(tk)) :: sktroe
          real*8 fsubc
          real*8 skfo_local(size(tk))
          if (present (fc)) then; fsubc=fc; else; fsubc=0.6d0; end if
          skfo_local(:) = skfo(af,npwr,ai,mpwr,tk,ad)
          sktroe(:) = skterlp(af,npwr,ae,tk,ad) * fsubc** &
     &                (1.0d0/(1.0d0+log10(skfo_local(:))**2)) / &
     &                (1.0d0+skfo(af,npwr,ai,mpwr,tk,ad))
        END FUNCTION sktroe
!
!.... skho2dis (temperature ,adcol ,water)
!
!_3_
!
!.... JPL 06-2
!
        FUNCTION skho2dis (tk,ad,h2o)
!
!....       HO2 + HO2 (+ H2O) = H2O2 + O2
!....
!.... This routine returns a bimolecular rate constant that accounts
!.... for the pressure and the H2O dependence of the reaction.
!.... Water is assumed to be passed in as a number density.
!
          real*8  &
     &      tk(:) ,ad(:) ,h2o(:)
          real*8, DIMENSION (size(tk)) :: skho2dis
          skho2dis(:) = (3.5D-13 * exp(430.0D0/tk(:)) +  &
     &                   1.7D-33 * ad(:) * exp(1000.0D0/tk(:))) *  &
     &                  (1.0D0+1.4D-21*h2o(:)*exp(2200.0D0/tk(:)))
        END FUNCTION skho2dis
!
!.... skcooh (temperature ,adcol)
!
!_1_
!
!.... JPL 06-2 ; now CO + OH is composed of two separate reactions!
!.... Now it's density and temperature dependent.
!         M
! OH + CO -> HOCO , but HOCO + O2 -> HO2 + CO2 quickly ; termolecular
!         M
! OH + CO -> H + CO2, but H + O2 -> HO2 quickly ; chemical activation reaction
!
        FUNCTION skcooh (tk,ad)
!
!                  M
!....      1) OH + CO = H + CO2; assume H+O2->HO2 quick
!....      2) OH + CO = HO2 + CO2
!....
!.... Pressure in hPa
!
          real*8  &
     &      tk(:), ad(:), af ,npwr ,ai ,mpwr ,ae
          real*8, DIMENSION (size(tk)) :: skcooh
          real*8, DIMENSION (size(tk)) :: skcoohlp
          real*8, DIMENSION (size(tk)) :: skcoohhp
          real*8, DIMENSION (size(tk)) :: skcooh1
          real*8, DIMENSION (size(tk)) :: skcooh2
          real*8 fsubc

! 1)
          fsubc=0.6d0
          af=1.5d-13
          npwr=-0.6
          ai=2.1d9
          mpwr=-6.1

          skcoohlp(:) = af * (300.0d0/tk(:))**npwr
          skcoohhp(:) = skhp(ai,mpwr,tk)
          skcooh2(:)  =  skcoohlp(:)/( skcoohhp(:)/ad(:) )
          skcooh(:)   = ( skcoohlp(:)/(1.0d0+skcooh2(:)) ) * fsubc**  &
     & ( 1.0d0/( 1.0d0+(log10(skcooh2(:)))**2 ) )

          fsubc=0.6d0
          af=5.9d-33
          npwr=1.4
          ai=1.1d-12
          mpwr=-1.3

! 2)
          skcooh(:) = skcooh(:) +  &
     &   sktroe(af,npwr,ai,mpwr,0.0D0,tk(:),ad(:))

        END FUNCTION skcooh
!
!.... skohch4 (temperature)
! _1_
!
!.... Harvard/GMI JPL 06-2
!
        FUNCTION skohch4 (tk)
!
!....      OH + CH4 = MO2 + H2O
!
          real*8  &
     &      tk(:)
          real*8, DIMENSION(size(tk)) :: skohch4
          skohch4(:) = 2.80D-14 * tk(:)**0.667d0 *  &
     &                 exp(-1575.0d0 / tk(:))
        END FUNCTION skohch4
!
!.... skmo2dis_1 (temperature)
! _2_
!
!.... Harvard/GMI
!
        FUNCTION skmo2dis_1 (tk)
!
!....      MO2 + MO2 = MOH + CH2O + O2
!....
!======================================================================

          real*8  &
     &      tk(:)
          real*8, DIMENSION(size(tk)) :: skmo2dis_1
          skmo2dis_1(:) = 9.50D-14 * exp(390.0d0 / tk(:)) * (1.0d0 /  &
     &                    (1.0d0 + 26.2d0 * exp(-1130.0d0 / tk(:))))
        END FUNCTION skmo2dis_1
!
!.... skmo2dis_2 (temperature)
! _3_
!
!.... Harvard/GMI
!
        FUNCTION skmo2dis_2 (tk)
!
!....      MO2 + MO2 = 2 CH2O + 2 HO2
!....
!======================================================================
!
          real*8  &
     &      tk(:)
          real*8, DIMENSION(size(tk)) :: skmo2dis_2
          skmo2dis_2(:) = 9.50D-14 * exp(390.0d0 / tk(:)) * (1.0d0 /  &
     &                    (1.0d0 + 0.04d0 * exp(1130.0d0 / tk(:))))
        END FUNCTION skmo2dis_2
!
!.... skohhno3 (temperature ,adcol)
!
!_1_
!
!.... JPL 06-2
!
        FUNCTION skohhno3 (tk,ad)
!
!                    M
!....      OH + HNO3 = NO3 + H2O
!....
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION (size(tk)) :: skohhno3
          real*8  &
     &      xxx1(size(tk)) ,xxx2(size(tk))
          xxx1(:)     = ad(:) * 6.5d-34 * exp(1335.0d0 / tk(:))
          xxx2(:)     = 2.7d-17 * exp(2199.0d0 / tk(:))
          skohhno3(:) = 2.4d-14 * exp(460.0d0 / tk(:)) +  &
     &                  xxx1(:) / (1.0d0 + xxx1(:) / xxx2(:))
        END FUNCTION skohhno3
!
!.... skpanan (temperature,adcol)
! _19_
!
!.... Harvard/GMI JPL 06-2
!
        FUNCTION skpanan (tk,ad)
!
!.... PAN = MCO3 + NO2
!....
          real*8  &
     &      tk(:),ad(:)
          real*8, DIMENSION(size(tk)) :: skpanan
       skpanan(:) =  &
     &  sktroe(9.700D-29,5.60D0,9.300D-12,1.50D0,0.0D0,tk(:),ad(:))/  &
     &  skarr(9.000D-29, -14000.0D0,tk(:))
        END FUNCTION skpanan
!
!.... skro2noabs_1 (temperature ,adcol)
! _13_
!
!.... Harvard/GMI
!
        FUNCTION fyrno3(xcarbn,tk,ad)
          real*8  &
     &      xcarbn
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: fyrno3
          real*8  &
     &      aaa(size(tk)) ,rarb(size(tk))  &
     &     ,xxyn(size(tk)) ,yyyn(size(tk)) ,zzyn(size(tk))
          xxyn(:)   = 1.94D-22 * exp(0.97d0 * xcarbn) * ad(:) *  &
     &                (300.0d0 / tk(:))**0.0d0
          yyyn(:)   = 0.826d0 * (300.0d0 / tk(:))**8.1d0
          aaa(:)    = log10(xxyn(:) / yyyn(:))
          zzyn(:)   = 1.0d0 / (1.0d0 + aaa(:) * aaa(:))
          rarb(:)   = (xxyn(:) / (1.0d0 + (xxyn(:) / yyyn(:)))) *  &
     &                0.411d0**zzyn(:)
          fyrno3(:) = rarb(:) / (1.0d0 + rarb(:))
        END FUNCTION fyrno3
        FUNCTION skro2noabs_1 (tk,ad)
!
!....      R4O2 + NO = NO2 + 0.320 ACET + 0.190 MEK +
!....                        0.180 MO2 + 0.270 HO2 +
!....                        0.320 ALD2 + 0.130 RCHO +
!....                        0.050 A3O2 + 0.180 B3O2 +
!....                        0.320 ETO2
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_1
          skro2noabs_1(:) = 2.70D-12 * exp(360.0d0 / tk(:)) *  &
     &                      (1.0d0 - fyrno3(4.5d0,tk,ad))
        END FUNCTION skro2noabs_1
!
!.... skro2noadd_1 (temperature ,adcol)
! _16_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_1 (tk,ad)
!
!.... R4O2 + NO = R4N2
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_1
          skro2noadd_1(:) = 2.70D-12 * exp(360.0d0 / tk(:)) *  &
     &                      fyrno3(4.5d0,tk,ad)
        END FUNCTION skro2noadd_1
!
!.... skro2noabs_2 (temperature ,adcol)
! _14_
!
!.... Harvard/GMI
!
        FUNCTION skro2noabs_2 (tk,ad)
!
!.... RIO2 + NO = NO2 + 0.864 HO2 + 0.690 CH2O +
!....                   0.402 MVK + 0.288 MACR +
!....                   0.136 RIO1 + 0.127 IALD
!.... RIO1 + NO = NO2 + IALD + HO2 + 0.750 CH2O
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_2
          skro2noabs_2(:) = 2.70D-12 * exp(360.0d0 / tk(:)) *  &
     &                      (1.0d0 - fyrno3(5.0d0,tk,ad))
        END FUNCTION skro2noabs_2
!
!.... skro2noadd_2 (temperature ,adcol)
! _17_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_2 (tk,ad)
!
!.... RIO2 + NO = HNO3
!.... RIO1 + NO = HNO3
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_2
          skro2noadd_2(:) = 2.70D-12 * exp(360.0d0 / tk(:)) *  &
     &                      fyrno3(5.0d0,tk,ad)
        END FUNCTION skro2noadd_2
!
!.... skro2noabs_3 (temperature ,adcol)
! _15_
!
!.... Harvard/GMI
!
        FUNCTION skro2noabs_3 (tk,ad)
!
!.... VRO2 + NO = NO2 + 0.280 HO2 + 0.280 CH2O +
!....                   0.720 MCO3 + 0.720 GLYC +
!....                   0.280 MGLY
!.... MRO2 + NO = NO2 + HO2 + 0.170 MGLY + 0.830 HAC +
!....                         0.830 CO + 0.170 CH2O
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_3
          skro2noabs_3(:) = 2.70D-12 * exp(360.0d0 / tk(:)) *  &
     &                      (1.0d0 - fyrno3(4.0d0,tk,ad))
        END FUNCTION skro2noabs_3
!
!.... skro2noadd_3 (temperature ,adcol)
! _18_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_3 (tk,ad)
!
!.... VRO2 + NO = HNO3
!.... MRO2 + NO = HNO3
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_3
          skro2noadd_3(:) = 2.70D-12 * exp(360.0d0 / tk(:)) *  &
     &                      fyrno3(4.0d0,tk,ad)
        END FUNCTION skro2noadd_3
!
!.... skppndecomp (temperature,adcol)
! _25_
!
!.... Harvard/GMI JPL 06-2
!
        FUNCTION skppndecomp (tk,ad)
!
!....     PPN   = NO2 + RCO3
!....     GPAN  = NO2 + GCO3
!....     PMN   = NO2 + MAO3
!....
!======================================================================

          real*8  &
     &      tk(:), ad(:)
          real*8, dimension(size(tk)) :: skppndecomp
      skppndecomp(:)=sktroe(9.000D-28, 8.9d0,7.700D-12, 0.2d0, 0.0d0,  &
     &      tk(:), ad(:))/skarr(9.00D-29 ,-14000.0D+00, tk(:))

        END FUNCTION skppndecomp
!
!.... skacetoh (temperature)
! _24_
!
!.... Harvard/GMI JPL 06-2
!
        FUNCTION skacetoh (tk)
!
!....      ACET + OH = ATO2 + H2O
!....
!======================================================================

          real*8  &
     &      tk(:)
          real*8, dimension(size(tk)) :: skacetoh
      skacetoh(:)=skarr(3.8200D-11 ,2000.0D+00 ,tk(:)) +1.330D-13
        END FUNCTION skacetoh
!
!.... skohmek (temperature)
! _6_
!
!.... Harvard/GMI
!
        FUNCTION skohmek (tk)
!
!....      OH + MEK = KO2 + H2O
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:)
          real*8, DIMENSION(size(tk)) :: skohmek
          skohmek(:) = 2.92D-12 * (300.0d0 / tk(:))**(-2.0d0) *  &
     &                 exp(414.0d0 / tk(:))
        END FUNCTION skohmek
!
!.... skho2mco3_1 (temperature)
! _7_
!
!.... Harvard/GMI JPL 06-2
!
        FUNCTION skho2mco3_1 (tk)
!
!....      HO2 + MCO3 = ACTA + O3
!
          real*8  &
     &      tk(:)
          real*8, DIMENSION(size(tk)) :: skho2mco3_1
          skho2mco3_1(:) = 4.30D-13 * exp(1040.0d0 / tk(:)) * (1.0d0 /  &
     &                     (1.0d0 + 37.0d0 * exp(-660.0d0 / tk(:))))
        END FUNCTION skho2mco3_1
!
!.... skho2mco3_2 (temperature)
! _8_
!
!.... Harvard/GMI JPS 06-2
!
        FUNCTION skho2mco3_2 (tk)
!
!....      HO2 + MCO3 = MAP
!
          real*8  &
     &      tk(:)
          real*8, DIMENSION(size(tk)) :: skho2mco3_2
          skho2mco3_2(:) = 4.30D-13 * exp(1040.0d0 / tk(:)) * (1.0d0 /  &
     &                     (1.0d0 + 2.70D-02 * exp(660.0d0 / tk(:))))
        END FUNCTION skho2mco3_2
!
!.... skohglyx (temperature ,oxygen)
! _4_
!
!.... Harvard/GMI
!
        FUNCTION skohglyx (tk,o2)
!
!                    O2
!....      OH + GLYX = HO2 + 2 CO
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,o2(:)
          real*8, DIMENSION(size(tk)) :: skohglyx
          skohglyx(:) = 1.10D-11 * (o2(:) + 3.5D+18) /  &
     &                  (2.0d0 * o2(:) + 3.5D+18)
        END FUNCTION skohglyx
!
!.... skno3glyx (temperature ,oxygen)
! _5_
!
!.... Harvard/GMI
!
        FUNCTION skno3glyx (tk,o2)
!
!                     O2
!....      NO3 + GLYX = HO2 + 2 CO
!....
!.... PSC - 8/8/2002
!....
!
          real*8  &
     &      tk(:) ,o2(:)
          real*8, DIMENSION(size(tk)) :: skno3glyx
          skno3glyx(:) = 1.40D-12 * exp(-1860.0d0 / tk(:)) *  &
     &                   (o2(:) + 3.5D+18) / (2.0d0 * o2(:) + 3.5D+18)
        END FUNCTION skno3glyx
!
!.... sktrs_ho2 (temperature,sadcol2,adcol,radA,NSADaer,NSADdust,tropp,pressur
!
!_1_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_ho2 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_ho2
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk)), dfkg(size(tk))
          integer ntotA,jj,NSADaer,NSADdust

          ntotA=NSADaer+NSADdust

!=======================================================================
!     HO2 + tropospheric aerosol = 0.5 H2O2
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)

! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_ho2(:)= 0.d0

          gamma       = 0.2d0
          pi          = acos(-1.0d0)

! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &               + 1.D0/mw(IHO2) )**0.5d0

! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &                (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IHO2)))**0.5d0

! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA

          where( sad(jj,:) > 0.0d0 )
            sktrs_ho2(:) = sktrs_ho2(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where

       enddo

          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_ho2 = 0.0d0
          end if

        END FUNCTION sktrs_ho2
!
!.... sktrs_no2 (temperature, sadcol2, adcol, radA, NSADaer,NSADdust,tropp, pr
!
!_2_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_no2 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_no2
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk)), dfkg(size(tk))
          integer ntotA,jj,NSADaer,NSADdust

           ntotA=NSADaer+NSADdust

!=======================================================================
!     NO2 + tropospheric aerosol = 0.5 HNO3 + 0.5 HONO
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)

! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_no2(:)= 0.d0

          gamma       = 1.0d-04
          pi          = acos(-1.0d0)

! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &              + 1.D0/mw(INO2) )**0.5d0

! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(INO2)))**0.5d0

! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA

          where( sad(jj,:) > 0.0d0 )
            sktrs_no2(:) = sktrs_no2(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where

       enddo

          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_no2 = 0.0d0
          end if

        END FUNCTION sktrs_no2
!
!.... sktrs_no3 (temperature, sadcol2, adcol, radA,NSADaer,NSADdust,tropp, pre
!
!_3_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_no3 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_no3
          real*8  &
     &      gamma ,pi
          real*8  &
     &      avgvel(size(tk)), dfkg(size(tk))
          integer NSADaer,NSADdust,ntotA,jj

          ntotA=NSADaer+NSADdust

!=======================================================================
!     NO3 + tropospheric aerosol = HNO3
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)

! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_no3(:)= 0.d0

          gamma       = 1.0d-03
          pi          = acos(-1.0d0)

! calculate gas phase diffusion coefficient (cm2/s)
          dfkg(:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &              + 1.D0/mw(INO3) )**0.5d0

! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(INO3)))**0.5d0

! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA

          where( sad(jj,:) > 0.0d0 )
            sktrs_no3(:) =  sktrs_no3(:) +  &
     &           sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &           radA(jj,:) / dfkg(:) )**(-1.0d0)
          endwhere

       enddo

          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_no3 = 0.0d0
          end if

        END FUNCTION sktrs_no3
!
!.... sktrs_n2o5 (temperature,sadcol2,adcol,radA,FRH,NSADaer,NSADdust, tropp, 
!
!_4_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_n2o5 (tk, sad, ad, radA,FRH, NSADaer,NSADdust,ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  &
     &      tk(:) ,sad(:,:), ad(:), radA(:,:), FRH(:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_n2o5
          real*8  &
     &      pi, FRH_P(size(tk)),ttk(size(tk)),fact(size(tk))
          real*8  &
     &      avgvel(size(tk)),dfkg(size(tk)),gamma(size(tk))
          integer NSADaer,NSADdust,ntotA,jj

          ntotA=NSADaer+NSADdust
!=======================================================================
! N2O5 + tropospheric aerosol = 2 HNO3                  Branch 4
!=======================================================================
! FRH = relative humidity fraction (0-1)
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)

! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_n2o5(:)= 0.d0

          pi           = acos(-1.0d0)

! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) =  9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &                + 1.D0/mw(IN2O5) )**0.5d0

! calculate mean molecular speed (cm/s)
          avgvel(:)    = 100.0d0 *  &
     &                   (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(IN2O5)))**0.5d0

! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA

!***********************************************************
! calculate gamma which is function of aerosol type, T, & RH
! following Evans and Jacob, "Impact of new laboratory studies of N2O5
! hydrolysis on global model budgets of tropospheric nitrogen oxides,
! ozone, and OH"
!***********************************************************
      ! Convert RH to % (max = 100%)
      FRH_P(:)  = FRH(:) * 100.d0
      where( FRH_P(:) > 100.d0) FRH_P(:) = 100.d0

        gamma(:) = 0d0

! DUST
!    Based on unpublished Crowley work
      if(jj.le.7) gamma(:) = 0.01d0
! SULFATE
      if(jj.eq.8) then
!===========================================================
! RH dependence from Kane et al., Heterogenous uptake of
! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
! J. Phys. Chem. A , 2001, 105, 6465-6470
!===========================================================
            gamma(:) = 2.79d-4 + FRH_P(:)*(  1.30d-4 +  &
     &                        FRH_P(:)*( -3.43d-6 +  &
     &                        FRH_P(:)*(  7.52d-8 ) ) )

!===========================================================
! Temperature dependence factor (Cox et al, Cambridge UK)
! is of the form:
!
!          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
! FACT = -------------------------------------------------
!                     10^( LOG10( G294 ) )
!
! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
!
! For computational speed, replace LOG10( 1e-2 ) with -2
! and replace 10^( LOG10( G294 ) ) with G294
!===========================================================
            ttk(:) = tk(:)
            where( ttk(:) < 282d0) ttk(:) = 282d0
            fact(:) = 10d0**( -2d0 - 4d-2*(ttk(:) - 294.d0))/1d-2

            ! Apply temperature dependence
            gamma(:) = gamma(:) * fact(:)
       endif
! BLACK CARBON
!     from IUPAC
      if(jj.eq.9) gamma(:) = 0.005d0
! ORGANIC CARBON
      if(jj.eq.10) then
!===========================================================
! Based on Thornton, Braban and Abbatt, 2003
! N2O5 hydrolysis on sub-micron organic aerosol: the effect
! of relative humidity, particle phase and particle size
!===========================================================
            where ( FRH_P(:) >= 57d0 ) gamma(:) = 0.03d0
            where ( FRH_P(:) <  57d0 ) gamma(:) = FRH_P(:) * 5.2d-4
!
! Bryan & Jules 01/14/05
! Set gamma to very samll number when gamma=0, which occurs
! when relative humidity is 0 (e.g., UT/LS polar night).
!
            where ( gamma(:) == 0.0d0) gamma(:) = 0.03d-2

       endif
! SEA SALT
!     Based on IUPAC recomendation
      if(jj.ge.11) then
            where ( FRH_P(:) >= 62 ) gamma(:) = 0.03d0
            where ( FRH_P(:) <  62 ) gamma(:) = 0.005d0
      endif

! end calculation of gamma
!***********************************************************
! calculate loss rate
!***********************************************************
          where( sad(jj,:) > 0.0d0 )
            sktrs_n2o5(:) = sktrs_n2o5(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma(:) * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where

       enddo

          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_n2o5 = 0.0d0
          end if

        END FUNCTION sktrs_n2o5
      END
