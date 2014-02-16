!=======================================================================
!
! $Id: setkin_lchem.h,v 1.5 2011-08-09 22:12:58 mrdamon Exp $
!
! FILE
!   setkin_lchem.h - character labels for species and reactions
!             (setkin_lchem.h)
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   Include file that provides ascii strings identifying reactions
!   and species
!
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!=======================================================================


      integer i

      character*16 lchemvar(NSP), ldynvar(NDYN)
      character*80 lqkchem(NUM_K),  lqjchem(NUM_J)
!
!.... All species labels
!
      data lchemvar(1) /"A3O2"/
      data lchemvar(2) /"ACTA"/
      data lchemvar(3) /"ALD2"/
      data lchemvar(4) /"ALK4"/
      data lchemvar(5) /"ATO2"/
      data lchemvar(6) /"B3O2"/
      data lchemvar(7) /"C2H6"/
      data lchemvar(8) /"C3H8"/
      data lchemvar(9) /"CH2O"/
      data lchemvar(10) /"CO"/
      data lchemvar(11) /"EOH"/
      data lchemvar(12) /"ETO2"/
      data lchemvar(13) /"ETP"/
      data lchemvar(14) /"GCO3"/
      data lchemvar(15) /"GLYC"/
      data lchemvar(16) /"GLYX"/
      data lchemvar(17) /"GP"/
      data lchemvar(18) /"GPAN"/
      data lchemvar(19) /"H2O2"/
      data lchemvar(20) /"HAC"/
      data lchemvar(21) /"HCOOH"/
      data lchemvar(22) /"HNO2"/
      data lchemvar(23) /"HNO3"/
      data lchemvar(24) /"HNO4"/
      data lchemvar(25) /"HO2"/
      data lchemvar(26) /"IALD"/
      data lchemvar(27) /"IAO2"/
      data lchemvar(28) /"IAP"/
      data lchemvar(29) /"INO2"/
      data lchemvar(30) /"INPN"/
      data lchemvar(31) /"ISN1"/
      data lchemvar(32) /"ISNP"/
      data lchemvar(33) /"ISOP"/
      data lchemvar(34) /"KO2"/
      data lchemvar(35) /"MACR"/
      data lchemvar(36) /"MAN2"/
      data lchemvar(37) /"MAO3"/
      data lchemvar(38) /"MAOP"/
      data lchemvar(39) /"MAP"/
      data lchemvar(40) /"MCO3"/
      data lchemvar(41) /"MEK"/
      data lchemvar(42) /"MGLY"/
      data lchemvar(43) /"MO2"/
      data lchemvar(44) /"MOH"/
      data lchemvar(45) /"MP"/
      data lchemvar(46) /"MRO2"/
      data lchemvar(47) /"MRP"/
      data lchemvar(48) /"MVK"/
      data lchemvar(49) /"MVN2"/
      data lchemvar(50) /"N2O5"/
      data lchemvar(51) /"NO"/
      data lchemvar(52) /"NO2"/
      data lchemvar(53) /"NO3"/
      data lchemvar(54) /"O3"/
      data lchemvar(55) /"OH"/
      data lchemvar(56) /"PAN"/
      data lchemvar(57) /"PMN"/
      data lchemvar(58) /"PO2"/
      data lchemvar(59) /"PP"/
      data lchemvar(60) /"PPN"/
      data lchemvar(61) /"PRN1"/
      data lchemvar(62) /"PRPE"/
      data lchemvar(63) /"PRPN"/
      data lchemvar(64) /"R4N1"/
      data lchemvar(65) /"R4N2"/
      data lchemvar(66) /"R4O2"/
      data lchemvar(67) /"R4P"/
      data lchemvar(68) /"RA3P"/
      data lchemvar(69) /"RB3P"/
      data lchemvar(70) /"RCHO"/
      data lchemvar(71) /"RCO3"/
      data lchemvar(72) /"RCOOH"/
      data lchemvar(73) /"RIO1"/
      data lchemvar(74) /"RIO2"/
      data lchemvar(75) /"RIP"/
      data lchemvar(76) /"ROH"/
      data lchemvar(77) /"RP"/
      data lchemvar(78) /"VRO2"/
      data lchemvar(79) /"VRP"/
      data lchemvar(80) /"ACET"/
      data lchemvar(81) /"CH4"/
      data lchemvar(82) /"H2"/
      data lchemvar(83) /"H2O"/
      data lchemvar(84) /"Total density"/
      data lchemvar(85) /"SYNOZ"/
!
!.... Dynamic (transported) species labels
!
      data ldynvar(1) /"A3O2"/
      data ldynvar(2) /"ACTA"/
      data ldynvar(3) /"ALD2"/
      data ldynvar(4) /"ALK4"/
      data ldynvar(5) /"ATO2"/
      data ldynvar(6) /"B3O2"/
      data ldynvar(7) /"C2H6"/
      data ldynvar(8) /"C3H8"/
      data ldynvar(9) /"CH2O"/
      data ldynvar(10) /"CO"/
      data ldynvar(11) /"EOH"/
      data ldynvar(12) /"ETO2"/
      data ldynvar(13) /"ETP"/
      data ldynvar(14) /"GCO3"/
      data ldynvar(15) /"GLYC"/
      data ldynvar(16) /"GLYX"/
      data ldynvar(17) /"GP"/
      data ldynvar(18) /"GPAN"/
      data ldynvar(19) /"H2O2"/
      data ldynvar(20) /"HAC"/
      data ldynvar(21) /"HCOOH"/
      data ldynvar(22) /"HNO2"/
      data ldynvar(23) /"HNO3"/
      data ldynvar(24) /"HNO4"/
      data ldynvar(25) /"HO2"/
      data ldynvar(26) /"IALD"/
      data ldynvar(27) /"IAO2"/
      data ldynvar(28) /"IAP"/
      data ldynvar(29) /"INO2"/
      data ldynvar(30) /"INPN"/
      data ldynvar(31) /"ISN1"/
      data ldynvar(32) /"ISNP"/
      data ldynvar(33) /"ISOP"/
      data ldynvar(34) /"KO2"/
      data ldynvar(35) /"MACR"/
      data ldynvar(36) /"MAN2"/
      data ldynvar(37) /"MAO3"/
      data ldynvar(38) /"MAOP"/
      data ldynvar(39) /"MAP"/
      data ldynvar(40) /"MCO3"/
      data ldynvar(41) /"MEK"/
      data ldynvar(42) /"MGLY"/
      data ldynvar(43) /"MO2"/
      data ldynvar(44) /"MOH"/
      data ldynvar(45) /"MP"/
      data ldynvar(46) /"MRO2"/
      data ldynvar(47) /"MRP"/
      data ldynvar(48) /"MVK"/
      data ldynvar(49) /"MVN2"/
      data ldynvar(50) /"N2O5"/
      data ldynvar(51) /"NO"/
      data ldynvar(52) /"NO2"/
      data ldynvar(53) /"NO3"/
      data ldynvar(54) /"O3"/
      data ldynvar(55) /"OH"/
      data ldynvar(56) /"PAN"/
      data ldynvar(57) /"PMN"/
      data ldynvar(58) /"PO2"/
      data ldynvar(59) /"PP"/
      data ldynvar(60) /"PPN"/
      data ldynvar(61) /"PRN1"/
      data ldynvar(62) /"PRPE"/
      data ldynvar(63) /"PRPN"/
      data ldynvar(64) /"R4N1"/
      data ldynvar(65) /"R4N2"/
      data ldynvar(66) /"R4O2"/
      data ldynvar(67) /"R4P"/
      data ldynvar(68) /"RA3P"/
      data ldynvar(69) /"RB3P"/
      data ldynvar(70) /"RCHO"/
      data ldynvar(71) /"RCO3"/
      data ldynvar(72) /"RCOOH"/
      data ldynvar(73) /"RIO1"/
      data ldynvar(74) /"RIO2"/
      data ldynvar(75) /"RIP"/
      data ldynvar(76) /"ROH"/
      data ldynvar(77) /"RP"/
      data ldynvar(78) /"VRO2"/
      data ldynvar(79) /"VRP"/
      data ldynvar(80) /"SYNOZ"/
!
!.... Thermal reaction labels
!
      data (lqkchem(i), i=1,10) / &
     & 'NO + O3 = NO2', &
     & 'O3 + OH = HO2', &
     & 'HO2 + O3 = OH', &
     & 'NO2 + O3 = NO3', &
     & 'MO2 + O3 = CH2O + HO2', &
     & 'OH + OH = H2O + O3', &
     & 'OH + OH = H2O2', &
     & 'HO2 + OH = H2O', &
     & 'H2O2 + OH = H2O + HO2', &
     & 'HO2 + NO = NO2 + OH' /

      data (lqkchem(i), i=11,20) / &
     & 'HO2 + HO2 = H2O2', &
     & 'H2 + OH = H2O + HO2', &
     & 'CO + OH = HO2', &
     & 'CH4 + OH = H2O + MO2', &
     & 'MO2 + NO = CH2O + HO2 + NO2', &
     & 'HO2 + MO2 = MP', &
     & 'MO2 + MO2 = CH2O + MOH', &
     & 'MO2 + MO2 = 2 CH2O + 2 HO2', &
     & 'MP + OH = H2O + MO2', &
     & 'MP + OH = CH2O + H2O + OH' /

      data (lqkchem(i), i=21,30) / &
     & 'CH2O + OH = CO + H2O + HO2', &
     & 'NO2 + OH = HNO3', &
     & 'HNO3 + OH = H2O + NO3', &
     & 'NO + OH = HNO2', &
     & 'HNO2 + OH = H2O + NO2', &
     & 'HO2 + NO2 = HNO4', &
     & 'HNO4 = HO2 + NO2', &
     & 'HNO4 + OH = H2O + NO2', &
     & 'HO2 + NO3 = NO2 + OH', &
     & 'NO + NO3 = 2 NO2' /

      data (lqkchem(i), i=31,40) / &
     & 'NO3 + OH = HO2 + NO2', &
     & 'NO2 + NO3 = N2O5', &
     & 'N2O5 = NO2 + NO3', &
     & 'HCOOH + OH = H2O + HO2', &
     & 'MOH + OH = CH2O + HO2', &
     & 'NO2 + NO3 = NO + NO2', &
     & 'CH2O + NO3 = CO + HNO3 + HO2', &
     & 'ALD2 + OH = H2O + MCO3', &
     & 'ALD2 + NO3 = HNO3 + MCO3', &
     & 'MCO3 + NO2 = PAN' /

      data (lqkchem(i), i=41,50) / &
     & 'PAN = MCO3 + NO2', &
     & 'MCO3 + NO = MO2 + NO2', &
     & 'C2H6 + OH = ETO2 + H2O', &
     & 'ETO2 + NO = ALD2 + HO2 + NO2', &
     & 'C3H8 + OH = B3O2', &
     & 'C3H8 + OH = A3O2', &
     & 'A3O2 + NO = HO2 + NO2 + RCHO', &
     & 'NO + PO2 = ALD2 + CH2O + HO2 + NO2', &
     & 'ALK4 + OH = R4O2', &
     & 'NO + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O' /

      data (lqkchem(i), i=51,60) / &
     & 'NO + R4O2 = R4N2', &
     & 'NO + R4N1 =  0.75 ALD2 +  0.39 CH2O + 2 NO2 +  0.30 R4O2 +  ', &
     & 'ATO2 + NO =  0.19 CH2O +  0.77 HO2 +  0.19 MCO3 +  0.77 MGLY', &
     & 'KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2 +  0.07 R4N2', &
     & 'NO + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR', &
     & 'NO + RIO2 = HNO3', &
     & 'NO + RIO1 =  0.75 CH2O + HO2 + IALD + NO2', &
     & 'NO + RIO1 = HNO3', &
     & 'IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.24 GLYC +  0.17 GLYX ', &
     & 'NO + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3' /

      data (lqkchem(i), i=61,70) / &
     & 'NO + VRO2 = HNO3', &
     & 'MRO2 + NO =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 ', &
     & 'MRO2 + NO = HNO3', &
     & 'MVN2 + NO =  0.30 CH2O +  0.60 GLYC +  0.10 HNO3 +  0.30 HO2', &
     & 'MAN2 + NO = CH2O + MGLY + 2 NO2', &
     & 'B3O2 + NO = ACET + HO2 + NO2', &
     & 'INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR', &
     & 'NO + PRN1 = ALD2 + CH2O + 2 NO2', &
     & 'ALK4 + NO3 = HNO3 + R4O2', &
     & 'OH + R4N2 = H2O + R4N1' /

      data (lqkchem(i), i=71,80) / &
     & 'ACTA + OH = H2O + MO2', &
     & 'OH + RCHO = H2O + RCO3', &
     & 'NO2 + RCO3 = PPN', &
     & 'GCO3 + NO2 = GPAN', &
     & 'MAO3 + NO2 = PMN', &
     & 'PPN = NO2 + RCO3', &
     & 'GPAN = GCO3 + NO2', &
     & 'NO + RCO3 = ETO2 + NO2', &
     & 'GCO3 + NO = CH2O + HO2 + NO2', &
     & 'MAO3 + NO = 4 CH2O + HO2 + NO2' /

      data (lqkchem(i), i=81,90) / &
     & 'NO3 + RCHO = HNO3 + RCO3', &
     & 'ACET + OH = ATO2 + H2O', &
     & 'A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.', &
     & 'MO2 + PO2 =  0.50 ALD2 +  1.25 CH2O +  0.16 HAC + HO2 +  0.2', &
     & 'HO2 + R4O2 = R4P', &
     & 'HO2 + R4N1 = R4N2', &
     & 'ATO2 + HO2 = MCO3 + MO2', &
     & 'HO2 + KO2 = MGLY + MO2', &
     & 'HO2 + RIO2 = RIP', &
     & 'HO2 + RIO1 = RIP' /

      data (lqkchem(i), i=91,100) / &
     & 'HO2 + IAO2 = IAP', &
     & 'HO2 + ISN1 = ISNP', &
     & 'HO2 + VRO2 = VRP', &
     & 'HO2 + MRO2 = MRP', &
     & 'HO2 + MVN2 = ISNP', &
     & 'HO2 + MAN2 = ISNP', &
     & 'B3O2 + HO2 = RB3P', &
     & 'HO2 + INO2 = INPN', &
     & 'HO2 + PRN1 = PRPN', &
     & 'MEK + OH = H2O + KO2' /

      data (lqkchem(i), i=101,110) / &
     & 'ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O +  0.25 EOH + HO2 +  0.', &
     & 'MEK + NO3 = HNO3 + KO2', &
     & 'MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3', &
     & 'MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH', &
     & 'ATO2 + MO2 =  0.85 CH2O +  0.90 HO2 +  0.10 MCO3 +  0.25 MEK', &
     & 'KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3', &
     & 'MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.06 IALD +  0.14 MAC', &
     & 'MO2 + RIO1 =  1.13 CH2O + HO2 +  0.50 IALD +  0.25 MEK +  0.', &
     & 'IAO2 + MO2 =  0.95 CH2O +  0.15 CO +  0.13 GLYC +  0.09 GLYX', &
     & 'ISN1 + MO2 =  0.75 CH2O +  0.50 GLYC +  0.50 HAC +  0.50 HO2' /

      data (lqkchem(i), i=111,120) / &
     & 'MO2 + VRO2 =  0.89 CH2O +  0.36 GLYC +  0.64 HO2 +  0.36 MCO', &
     & 'MO2 + MRO2 =  0.84 CH2O +  0.42 CO +  0.42 HAC + HO2 +  0.25', &
     & 'MO2 + MVN2 =  1.25 CH2O +  0.75 HO2 +  0.25 MCO3 +  0.25 MGL', &
     & 'MAN2 + MO2 =  1.25 CH2O +  0.50 HO2 +  0.50 MGLY +  0.25 MOH', &
     & 'B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.', &
     & 'INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MAC', &
     & 'MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH', &
     & 'EOH + OH = ALD2 + HO2', &
     & 'OH + ROH = HO2 + RCHO', &
     & 'ETO2 + ETO2 = 2 ALD2 + 2 HO2' /

      data (lqkchem(i), i=121,130) / &
     & 'ETO2 + ETO2 = ALD2 + EOH', &
     & 'ETO2 + HO2 = ETP', &
     & 'A3O2 + HO2 = RA3P', &
     & 'HO2 + PO2 = PP', &
     & 'HO2 + MCO3 = ACTA + O3', &
     & 'HO2 + MCO3 = MAP', &
     & 'HO2 + RCO3 =  0.30 O3 +  0.30 RCOOH +  0.70 RP', &
     & 'GCO3 + HO2 =  0.70 GP +  0.30 O3 +  0.30 RCOOH', &
     & 'HO2 + MAO3 =  0.70 MAOP +  0.30 O3 +  0.30 RCOOH', &
     & 'OH + PRPE = PO2' /

      data (lqkchem(i), i=131,140) / &
     & 'O3 + PRPE =  0.50 ALD2 +  0.54 CH2O +  0.42 CO +  0.06 H2 + ', &
     & 'PMN = MAO3 + NO2', &
     & 'OH + PMN =  2.23 CH2O +  0.59 HAC + 2 HO2 + NO2', &
     & 'O3 + PMN =  0.60 CH2O + HO2 + NO2', &
     & 'GLYC + OH =  0.80 GCO3 +  0.20 GLYX +  0.20 HO2', &
     & 'NO3 + PRPE = PRN1', &
     & 'GLYX + OH = 2 CO + HO2', &
     & 'MGLY + OH = CO + MCO3', &
     & 'GLYX + NO3 = 2 CO + HNO3 + HO2', &
     & 'MGLY + NO3 = CO + HNO3 + MCO3' /

      data (lqkchem(i), i=141,150) / &
     & 'ISOP + OH = RIO2', &
     & 'MVK + OH = VRO2', &
     & 'MACR + OH =  0.50 MAO3 +  0.50 MRO2', &
     & 'HAC + OH = HO2 + MGLY', &
     & 'A3O2 + MCO3 = HO2 + MO2 + RCHO', &
     & 'MCO3 + PO2 = ALD2 + CH2O + HO2 + MO2', &
     & 'A3O2 + MCO3 = ACTA + RCHO', &
     & 'MCO3 + PO2 = ACTA +  0.65 HAC +  0.35 RCHO', &
     & 'ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.39 MACR +', &
     & 'MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 + ' /

      data (lqkchem(i), i=151,160) / &
     & 'MACR + O3 =  0.70 CH2O +  0.20 CO +  0.28 HO2 +  0.80 MGLY +', &
     & 'ISOP + NO3 = INO2', &
     & 'MVK + NO3 = MVN2', &
     & 'MACR + NO3 = MAN2', &
     & 'MACR + NO3 = HNO3 + MAO3', &
     & 'MO2 + RCO3 = CH2O + ETO2 + HO2', &
     & 'GCO3 + MO2 = 2 CH2O + 2 HO2', &
     & 'MAO3 + MO2 = 2 CH2O + HO2 + MCO3', &
     & 'MO2 + RCO3 = CH2O + RCOOH', &
     & 'GCO3 + MO2 = CH2O + RCOOH' /

      data (lqkchem(i), i=161,170) / &
     & 'MAO3 + MO2 = CH2O + RCOOH', &
     & 'INPN + OH = INO2', &
     & 'OH + PRPN = PRN1', &
     & 'ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH', &
     & 'OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO', &
     & 'OH + RB3P =  0.50 B3O2 +  0.50 OH +  0.50 RCHO', &
     & 'OH + R4P =  0.50 OH +  0.50 R4O2 +  0.50 RCHO', &
     & 'OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3', &
     & 'OH + PP =  0.50 OH +  0.50 PO2 +  0.50 RCHO', &
     & 'GP + OH =  0.50 CH2O +  0.50 GCO3 +  0.50 OH' /

      data (lqkchem(i), i=171,180) / &
     & 'OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2', &
     & 'IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO', &
     & 'ISNP + OH =  0.50 ISN1 +  0.50 NO2 +  0.50 OH +  0.50 RCHO', &
     & 'OH + VRP =  0.50 OH +  0.50 RCHO +  0.50 VRO2', &
     & 'MRP + OH =  0.50 MRO2 +  0.50 OH +  0.50 RCHO', &
     & 'MAOP + OH =  0.50 MAO3 +  0.50 OH +  0.50 RCHO', &
     & 'MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH', &
     & 'C2H6 + NO3 = ETO2 + HNO3', &
     & 'IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3', &
     & 'IALD + O3 =  0.12 CH2O +  0.28 GLYC +  0.20 GLYX +  0.20 HAC' /

      data (lqkchem(i), i=181,190) / &
     & 'MCO3 + MCO3 = 2 MO2', &
     & 'MCO3 + MO2 = CH2O + HO2 + MO2', &
     & 'MCO3 + MO2 = ACTA + CH2O', &
     & 'MCO3 + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B', &
     & 'ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 +  0.80 MG', &
     & 'KO2 + MCO3 = ALD2 + MCO3 + MO2', &
     & 'MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MA', &
     & 'MCO3 + RIO1 =  0.75 CH2O + HO2 + IALD + MO2', &
     & 'IAO2 + MCO3 =  0.40 CH2O +  0.29 CO +  0.26 GLYC +  0.18 GLY', &
     & 'ISN1 + MCO3 = GLYC + HAC + MO2 + NO2' /

      data (lqkchem(i), i=191,200) / &
     & 'MCO3 + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MC', &
     & 'MCO3 + MRO2 =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.1', &
     & 'B3O2 + MCO3 = ACET + HO2 + MO2', &
     & 'MCO3 + R4N1 =  0.75 ALD2 +  0.39 CH2O + MO2 + NO2 +  0.30 R4', &
     & 'MCO3 + MVN2 = CH2O +  0.50 HO2 +  0.50 MCO3 +  0.50 MGLY + M', &
     & 'MAN2 + MCO3 = CH2O + MGLY + MO2 + NO2', &
     & 'INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MA', &
     & 'MCO3 + PRN1 = ALD2 + CH2O + MO2 + NO2', &
     & 'MCO3 + R4O2 = ACTA + MEK', &
     & 'ATO2 + MCO3 = ACTA + MEK' /

      data (lqkchem(i), i=201,210) / &
     & 'KO2 + MCO3 = ACTA + MEK', &
     & 'MCO3 + RIO2 = ACTA + MEK', &
     & 'MCO3 + RIO1 = ACTA + MEK', &
     & 'IAO2 + MCO3 = ACTA + MEK', &
     & 'MCO3 + VRO2 = ACTA + MEK', &
     & 'MCO3 + MRO2 = ACTA + MEK', &
     & 'MCO3 + R4N1 = ACTA + NO2 + RCHO', &
     & 'ISN1 + MCO3 = ACTA + NO2 + RCHO', &
     & 'MCO3 + MVN2 = ACTA + NO2 + RCHO', &
     & 'MAN2 + MCO3 = ACTA + NO2 + RCHO' /

      data (lqkchem(i), i=211,220) / &
     & 'INO2 + MCO3 = ACTA + NO2 + RCHO', &
     & 'MCO3 + PRN1 = ACTA + NO2 + RCHO', &
     & 'B3O2 + MCO3 = ACET + ACTA', &
     & 'ETO2 + MCO3 = ALD2 + HO2 + MO2', &
     & 'ETO2 + MCO3 = ACTA + ALD2', &
     & 'MCO3 + RCO3 = ETO2 + MO2', &
     & 'GCO3 + MCO3 = CH2O + HO2 + MO2', &
     & 'MAO3 + MCO3 = CH2O + MCO3 + MO2', &
     & 'NO3 + NO3 = 2 NO2', &
     & 'HO2 =  0.50 H2O2' /

      data (lqkchem(i), i=221,223) / &
     & 'NO2 =  0.50 HNO2 +  0.50 HNO3', &
     & 'NO3 = HNO3', &
     & 'N2O5 = 2 HNO3' /
!
!.... Photolytic reaction labels
!
      data (lqjchem(i), i=1,10) / &
     & 'O3 + hv = 2 OH', &
     & 'NO2 + hv = NO + O3', &
     & 'H2O2 + hv = 2 OH', &
     & 'MP + hv = CH2O + HO2 + OH', &
     & 'CH2O + hv = CO + 2 HO2', &
     & 'CH2O + hv = CO + H2', &
     & 'HNO3 + hv = NO2 + OH', &
     & 'HNO2 + hv = NO + OH', &
     & 'HNO4 + hv = NO3 + OH', &
     & 'NO3 + hv = NO2 + O3' /

      data (lqjchem(i), i=11,20) / &
     & 'NO3 + hv = NO', &
     & 'N2O5 + hv = NO2 + NO3', &
     & 'N2O5 + hv = NO + NO3 + O3', &
     & 'HNO4 + hv = HO2 + NO2', &
     & 'ALD2 + hv = CO + HO2 + MO2', &
     & 'ALD2 + hv = CH4 + CO', &
     & 'PAN + hv = MCO3 + NO2', &
     & 'RCHO + hv = CO + ETO2 + HO2', &
     & 'ACET + hv = MCO3 + MO2', &
     & 'MEK + hv = ETO2 + MCO3' /

      data (lqjchem(i), i=21,30) / &
     & 'GLYC + hv = CH2O + CO + 2 HO2', &
     & 'GLYX + hv = 2 CO + H2', &
     & 'GLYX + hv = 2 CO + 2 HO2', &
     & 'GLYX + hv = CH2O + CO', &
     & 'MGLY + hv = CO + HO2 + MCO3', &
     & 'MGLY + hv = ALD2 + CO', &
     & 'MVK + hv = CO + PRPE', &
     & 'MVK + hv = CH2O + CO + HO2 + MCO3', &
     & 'MVK + hv = MAO3 + MO2', &
     & 'MACR + hv = HO2 + MAO3' /

      data (lqjchem(i), i=31,40) / &
     & 'MACR + hv =  0.20 CH2O + CO +  1.80 HO2 +  0.20 MCO3 +  0.80', &
     & 'HAC + hv = CH2O + HO2 + MCO3', &
     & 'INPN + hv = HO2 + NO2 + OH + RCHO', &
     & 'PRPN + hv = HO2 + NO2 + OH + RCHO', &
     & 'ETP + hv = ALD2 + HO2 + OH', &
     & 'RA3P + hv = HO2 + OH + RCHO', &
     & 'RB3P + hv = HO2 + OH + RCHO', &
     & 'R4P + hv = HO2 + OH + RCHO', &
     & 'PP + hv = HO2 + OH + RCHO', &
     & 'RP + hv = ALD2 + HO2 + OH' /

      data (lqjchem(i), i=41,49) / &
     & 'GP + hv = CH2O + HO2 + OH', &
     & 'RIP + hv =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR ', &
     & 'IAP + hv =  0.67 CO +  0.26 GLYC +  0.19 H2 +  0.36 HAC + HO', &
     & 'ISNP + hv = HO2 + NO2 + OH + RCHO', &
     & 'VRP + hv =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 ', &
     & 'MRP + hv =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 M', &
     & 'MAOP + hv = HO2 + OH + RCHO', &
     & 'R4N2 + hv =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O', &
     & 'MAP + hv = MO2 + OH' /

!                                  --^--

