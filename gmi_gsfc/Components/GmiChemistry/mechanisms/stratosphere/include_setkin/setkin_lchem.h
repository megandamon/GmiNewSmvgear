!=======================================================================
!
! $Id: setkin_lchem.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
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
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
!
!=======================================================================


      integer i

      character*16 lchemvar(NSP), ldynvar(NDYN)
      character*80 lqkchem(NUM_K),  lqjchem(NUM_J)
!
!.... All species labels
!
      data lchemvar(1) /"O"/
      data lchemvar(2) /"O1D"/
      data lchemvar(3) /"O3"/
      data lchemvar(4) /"N2O"/
      data lchemvar(5) /"N"/
      data lchemvar(6) /"NO"/
      data lchemvar(7) /"NO2"/
      data lchemvar(8) /"NO3"/
      data lchemvar(9) /"N2O5"/
      data lchemvar(10) /"HNO3"/
      data lchemvar(11) /"HO2NO2"/
      data lchemvar(12) /"H2O"/
      data lchemvar(13) /"H"/
      data lchemvar(14) /"OH"/
      data lchemvar(15) /"HO2"/
      data lchemvar(16) /"H2O2"/
      data lchemvar(17) /"H2"/
      data lchemvar(18) /"CH4"/
      data lchemvar(19) /"CH3O2"/
      data lchemvar(20) /"CH3OOH"/
      data lchemvar(21) /"CH2O"/
      data lchemvar(22) /"CO"/
      data lchemvar(23) /"Cl"/
      data lchemvar(24) /"Cl2"/
      data lchemvar(25) /"ClO"/
      data lchemvar(26) /"OClO"/
      data lchemvar(27) /"Cl2O2"/
      data lchemvar(28) /"HCl"/
      data lchemvar(29) /"HOCl"/
      data lchemvar(30) /"ClONO2"/
      data lchemvar(31) /"Br"/
      data lchemvar(32) /"BrCl"/
      data lchemvar(33) /"BrO"/
      data lchemvar(34) /"HBr"/
      data lchemvar(35) /"HOBr"/
      data lchemvar(36) /"BrONO2"/
      data lchemvar(37) /"CH3Cl"/
      data lchemvar(38) /"CH3Br"/
      data lchemvar(39) /"CFCl3"/
      data lchemvar(40) /"CF2Cl2"/
      data lchemvar(41) /"CFC113"/
      data lchemvar(42) /"CFC114"/
      data lchemvar(43) /"CFC115"/
      data lchemvar(44) /"HCFC22"/
      data lchemvar(45) /"CCl4"/
      data lchemvar(46) /"CH3CCl3"/
      data lchemvar(47) /"HCFC141b"/
      data lchemvar(48) /"HCFC142b"/
      data lchemvar(49) /"CF3Br"/
      data lchemvar(50) /"CF2ClBr"/
      data lchemvar(51) /"CF2Br2"/
      data lchemvar(52) /"H2402"/
      data lchemvar(53) /"N2"/
      data lchemvar(54) /"O2"/
      data lchemvar(55) /"Total density"/
      data lchemvar(56) /"DEHYD"/
      data lchemvar(57) /"H2OAIR"/
!
!.... Dynamic (transported) species labels
!
      data ldynvar(1) /"O"/
      data ldynvar(2) /"O1D"/
      data ldynvar(3) /"O3"/
      data ldynvar(4) /"N2O"/
      data ldynvar(5) /"N"/
      data ldynvar(6) /"NO"/
      data ldynvar(7) /"NO2"/
      data ldynvar(8) /"NO3"/
      data ldynvar(9) /"N2O5"/
      data ldynvar(10) /"HNO3"/
      data ldynvar(11) /"HO2NO2"/
      data ldynvar(12) /"H2O"/
      data ldynvar(13) /"H"/
      data ldynvar(14) /"OH"/
      data ldynvar(15) /"HO2"/
      data ldynvar(16) /"H2O2"/
      data ldynvar(17) /"H2"/
      data ldynvar(18) /"CH4"/
      data ldynvar(19) /"CH3O2"/
      data ldynvar(20) /"CH3OOH"/
      data ldynvar(21) /"CH2O"/
      data ldynvar(22) /"CO"/
      data ldynvar(23) /"Cl"/
      data ldynvar(24) /"Cl2"/
      data ldynvar(25) /"ClO"/
      data ldynvar(26) /"OClO"/
      data ldynvar(27) /"Cl2O2"/
      data ldynvar(28) /"HCl"/
      data ldynvar(29) /"HOCl"/
      data ldynvar(30) /"ClONO2"/
      data ldynvar(31) /"Br"/
      data ldynvar(32) /"BrCl"/
      data ldynvar(33) /"BrO"/
      data ldynvar(34) /"HBr"/
      data ldynvar(35) /"HOBr"/
      data ldynvar(36) /"BrONO2"/
      data ldynvar(37) /"CH3Cl"/
      data ldynvar(38) /"CH3Br"/
      data ldynvar(39) /"CFCl3"/
      data ldynvar(40) /"CF2Cl2"/
      data ldynvar(41) /"CFC113"/
      data ldynvar(42) /"CFC114"/
      data ldynvar(43) /"CFC115"/
      data ldynvar(44) /"HCFC22"/
      data ldynvar(45) /"CCl4"/
      data ldynvar(46) /"CH3CCl3"/
      data ldynvar(47) /"HCFC141b"/
      data ldynvar(48) /"HCFC142b"/
      data ldynvar(49) /"CF3Br"/
      data ldynvar(50) /"CF2ClBr"/
      data ldynvar(51) /"CF2Br2"/
      data ldynvar(52) /"H2402"/
      data ldynvar(53) /"DEHYD"/
      data ldynvar(54) /"H2OAIR"/
!
!.... Thermal reaction labels
!
      data (lqkchem(i), i=1,10) /  &
     & 'O + O2 = O3',  &
     & 'O + O3 = 2 O2',  &
     & 'N2 + O1D = N2 + O',  &
     & 'O1D + O2 = O + O2',  &
     & 'O1D + O3 = 2 O2',  &
     & 'H2O + O1D = 2 OH',  &
     & 'H2 + O1D = H + OH',  &
     & 'N2O + O1D = N2 + O2',  &
     & 'N2O + O1D = 2 NO',  &
     & 'CH4 + O1D = CH3O2 + OH' /

      data (lqkchem(i), i=11,20) /  &
     & 'CH4 + O1D = CH2O + H + HO2',  &
     & 'CH4 + O1D = CH2O + H2',  &
     & 'CF2Cl2 + O1D = 2 Cl',  &
     & 'CFC113 + O1D = 3 Cl',  &
     & 'CFC114 + O1D = 2 Cl',  &
     & 'CFC115 + O1D = Cl',  &
     & 'HCFC22 + O1D = Cl',  &
     & 'HCFC141b + O1D = 2 Cl',  &
     & 'HCFC142b + O1D = Cl',  &
     & 'H + O2 = HO2' /

      data (lqkchem(i), i=21,30) /  &
     & 'H + O3 = O2 + OH',  &
     & 'H2 + OH = H + H2O',  &
     & 'O3 + OH = HO2 + O2',  &
     & 'O + OH = H + O2',  &
     & 'OH + OH = H2O + O',  &
     & 'HO2 + O = O2 + OH',  &
     & 'HO2 + O3 = 2 O2 + OH',  &
     & 'H + HO2 = 2 OH',  &
     & 'HO2 + OH = H2O + O2',  &
     & 'HO2 + HO2 = H2O2 + O2' /

      data (lqkchem(i), i=31,40) /  &
     & 'H2O + HO2 + HO2 = H2O + H2O2 + O2',  &
     & 'H2O2 + OH = H2O + HO2',  &
     & 'N + O2 = NO + O',  &
     & 'N + NO = N2 + O',  &
     & 'NO + O3 = NO2 + O2',  &
     & 'NO2 + OH = HNO3',  &
     & 'HO2 + NO = NO2 + OH',  &
     & 'NO2 + O = NO + O2',  &
     & 'NO2 + O3 = NO3 + O2',  &
     & 'HO2 + NO2 = HO2NO2' /

      data (lqkchem(i), i=41,50) /  &
     & 'NO3 + O = NO2 + O2',  &
     & 'NO + NO3 = 2 NO2',  &
     & 'NO2 + NO3 = N2O5',  &
     & 'N2O5 = NO2 + NO3',  &
     & 'HNO3 + OH = H2O + NO3',  &
     & 'HO2NO2 = HO2 + NO2',  &
     & 'HO2NO2 + OH = H2O + NO2 + O2',  &
     & 'Cl + O3 = ClO + O2',  &
     & 'Cl + H2 = H + HCl',  &
     & 'Cl + H2O2 = HCl + HO2' /

      data (lqkchem(i), i=51,60) /  &
     & 'Cl + HO2 = HCl + O2',  &
     & 'Cl + HO2 = ClO + OH',  &
     & 'ClO + O = Cl + O2',  &
     & 'ClO + OH = Cl + HO2',  &
     & 'ClO + OH = HCl + O2',  &
     & 'ClO + HO2 = HOCl + O2',  &
     & 'ClO + HO2 = HCl + O3',  &
     & 'ClO + NO = Cl + NO2',  &
     & 'ClO + NO2 = ClONO2',  &
     & 'ClO + ClO = 2 Cl + O2' /

      data (lqkchem(i), i=61,70) /  &
     & 'ClO + ClO = Cl2 + O2',  &
     & 'ClO + ClO = Cl + OClO',  &
     & 'ClO + ClO = Cl2O2',  &
     & 'Cl2O2 = 2 ClO',  &
     & 'HCl + OH = Cl + H2O',  &
     & 'HOCl + OH = ClO + H2O',  &
     & 'ClONO2 + O = ClO + NO3',  &
     & 'ClONO2 + OH = HOCl + NO3',  &
     & 'Cl + ClONO2 = Cl2 + NO3',  &
     & 'Br + O3 = BrO + O2' /

      data (lqkchem(i), i=71,80) /  &
     & 'Br + HO2 = HBr + O2',  &
     & 'Br + CH2O = CO + HBr + HO2',  &
     & 'BrO + O = Br + O2',  &
     & 'BrO + HO2 = HOBr + O2',  &
     & 'BrO + NO = Br + NO2',  &
     & 'BrO + NO2 = BrONO2',  &
     & 'BrO + ClO = Br + OClO',  &
     & 'BrO + ClO = Br + Cl + O2',  &
     & 'BrO + ClO = BrCl + O2',  &
     & 'BrO + BrO = 2 Br + O2' /

      data (lqkchem(i), i=81,90) /  &
     & 'HBr + OH = Br + H2O',  &
     & 'CO + OH = H',  &
     & 'CH4 + OH = CH3O2 + H2O',  &
     & 'CH2O + OH = CO + H2O + HO2',  &
     & 'CH2O + O = CO + HO2 + OH',  &
     & 'CH4 + Cl = CH3O2 + HCl',  &
     & 'CH2O + Cl = CO + HCl + HO2',  &
     & 'CH3O2 + NO = CH2O + HO2 + NO2',  &
     & 'CH3O2 + HO2 = CH3OOH + O2',  &
     & 'CH3OOH + OH = CH3O2 + H2O' /

      data (lqkchem(i), i=91,100) /  &
     & 'CH3Cl + OH = Cl + H2O + HO2',  &
     & 'CH3CCl3 + OH = 3 Cl + H2O',  &
     & 'HCFC22 + OH = Cl + H2O',  &
     & 'HCFC141b + OH = 2 Cl + H2O',  &
     & 'HCFC142b + OH = Cl + H2O',  &
     & 'CH3Cl + Cl = CO + 2 HCl + HO2',  &
     & 'CH3Br + OH = Br + H2O + HO2',  &
     & 'N2O5 = 2 HNO3',  &
     & 'ClONO2 = HNO3 + HOCl',  &
     & 'BrONO2 = HNO3 + HOBr' /

      data (lqkchem(i), i=101,110) /  &
     & 'ClONO2 + HCl = Cl2 + HNO3',  &
     & 'HCl + HOCl = Cl2 + H2O',  &
     & 'HCl + HOBr = BrCl + H2O',  &
     & 'N2O5 = 2 HNO3',  &
     & 'ClONO2 = HNO3 + HOCl',  &
     & 'BrONO2 = HNO3 + HOBr',  &
     & 'ClONO2 + HCl = Cl2 + HNO3',  &
     & 'HCl + HOCl = Cl2 + H2O',  &
     & 'HCl + HOBr = BrCl + H2O',  &
     & 'ClONO2 = HNO3 + HOCl' /

      data (lqkchem(i), i=111,120) /  &
     & 'BrONO2 = HNO3 + HOBr',  &
     & 'ClONO2 + HCl = Cl2 + HNO3',  &
     & 'HCl + HOCl = Cl2 + H2O',  &
     & 'BrONO2 + HCl = BrCl + HNO3',  &
     & 'HCl + HOBr = BrCl + H2O',  &
     & 'ClONO2 = HNO3 + HOCl',  &
     & 'BrONO2 = HNO3 + HOBr',  &
     & 'ClONO2 + HCl = Cl2 + HNO3',  &
     & 'HCl + HOCl = Cl2 + H2O',  &
     & 'BrONO2 + HCl = BrCl + HNO3' /

      data (lqkchem(i), i=121,122) /  &
     & 'HCl + HOBr = BrCl + H2O',  &
     & 'HNO3 = NO2 + OH' /
!
!.... Photolytic reaction labels
!
      data (lqjchem(i), i=1,10) /  &
     & 'O2 + hv = 2 O',  &
     & 'O3 + hv = O + O2',  &
     & 'O3 + hv = O1D + O2',  &
     & 'HO2 + hv = O + OH',  &
     & 'H2O2 + hv = 2 OH',  &
     & 'H2O + hv = H + OH',  &
     & 'NO + hv = N + O',  &
     & 'NO2 + hv = NO + O',  &
     & 'N2O + hv = N2 + O1D',  &
     & 'NO3 + hv = NO2 + O' /

      data (lqjchem(i), i=11,20) /  &
     & 'NO3 + hv = NO + O2',  &
     & 'N2O5 + hv = NO2 + NO3',  &
     & 'HNO3 + hv = NO2 + OH',  &
     & 'HO2NO2 + hv = NO3 + OH',  &
     & 'HO2NO2 + hv = HO2 + NO2',  &
     & 'Cl2 + hv = 2 Cl',  &
     & 'OClO + hv = ClO + O',  &
     & 'Cl2O2 + hv = 2 Cl + O2',  &
     & 'HOCl + hv = Cl + OH',  &
     & 'ClONO2 + hv = Cl + NO3' /

      data (lqjchem(i), i=21,30) /  &
     & 'ClONO2 + hv = ClO + NO2',  &
     & 'BrCl + hv = Br + Cl',  &
     & 'BrO + hv = Br + O',  &
     & 'HOBr + hv = Br + OH',  &
     & 'BrONO2 + hv = Br + NO3',  &
     & 'BrONO2 + hv = BrO + NO2',  &
     & 'CH3OOH + hv = CH2O + HO2 + OH',  &
     & 'CH2O + hv = CO + H2',  &
     & 'CH2O + hv = CO + H + HO2',  &
     & 'CH3Cl + hv = CH3O2 + Cl' /

      data (lqjchem(i), i=31,40) /  &
     & 'CCl4 + hv = 4 Cl',  &
     & 'CH3CCl3 + hv = 3 Cl',  &
     & 'CFCl3 + hv = 3 Cl',  &
     & 'CF2Cl2 + hv = 2 Cl',  &
     & 'CFC113 + hv = 3 Cl',  &
     & 'CFC114 + hv = 2 Cl',  &
     & 'CFC115 + hv = Cl',  &
     & 'HCFC141b + hv = 2 Cl',  &
     & 'HCFC142b + hv = Cl',  &
     & 'CH3Br + hv = Br + CH3O2' /

      data (lqjchem(i), i=41,44) /  &
     & 'CF3Br + hv = Br',  &
     & 'CF2Br2 + hv = 2 Br',  &
     & 'H2402 + hv = 2 Br',  &
     & 'CF2ClBr + hv = Br + Cl' /

!                                  --^--

