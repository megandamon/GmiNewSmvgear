!=======================================================================
!
! $Id: setkin_ratecalc.F90,v 1.5 2011-08-09 22:12:58 mrdamon Exp $
!
! ROUTINE
!   Calc_rate_Setkin - GMI version (setkin_ratecalc.F)
!   13 JUN 02 - PSC
!
! DESCRIPTION
!   Calculates and returns rates of kinetic processes for block
!   of boxes
!
! ARGUMENTS
!   INPUT
!    nblock  : number of blocks
!    numj    : number of photolysis rate constants
!    numk    : number of thermal rate constants
!    numspc  ; number of species
!    qk     : thermal reaction rate constants (cm^3 molecule^-1 s^-1)
!    qj     : photolysis rate constants (s^-1)
!    y      : concentrations of transported species (molecules cm^-3)
!   OUTPUT
!    qqk    : rates of thermal processes (molecules cm^-3 s^-1)
!    qqj    : rates of photolytic processes (molecules cm^-3 s^-1)
!
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!=======================================================================
      subroutine Calc_rate_Setkin &
     &  (nblock, numj, numk, numspc, &
     &   qk, qj, y, qqk, qqj)

      implicit none

      integer &
     &  nblock, &
     &  numj, numk, numspc

      real*8 &
     &  qj    (nblock, numj), &
     &  qk    (nblock, numk), &
     &  qqj   (nblock, numj), &
     &  qqk   (nblock, numk), &
     &  y     (nblock, numspc)

      integer &
     &  kloop

      do kloop = 1, nblock
!
!                  start thermal reactions
!
!....         NO + O3 = NO2
!
      qqk(kloop,1)=qk(kloop,1)*y(kloop,51)*y(kloop,54)
!
!....         O3 + OH = HO2
!
      qqk(kloop,2)=qk(kloop,2)*y(kloop,54)*y(kloop,55)
!
!....         HO2 + O3 = OH
!
      qqk(kloop,3)=qk(kloop,3)*y(kloop,25)*y(kloop,54)
!
!....         NO2 + O3 = NO3
!
      qqk(kloop,4)=qk(kloop,4)*y(kloop,52)*y(kloop,54)
!
!....         MO2 + O3 = CH2O + HO2
!
      qqk(kloop,5)=qk(kloop,5)*y(kloop,43)*y(kloop,54)
!
!....         OH + OH = H2O + O3
!
      qqk(kloop,6)=qk(kloop,6)*y(kloop,55)*y(kloop,55)
!
!....         OH + OH = H2O2
!
      qqk(kloop,7)=qk(kloop,7)*y(kloop,55)*y(kloop,55)
!
!....         HO2 + OH = H2O
!
      qqk(kloop,8)=qk(kloop,8)*y(kloop,25)*y(kloop,55)
!
!....         H2O2 + OH = H2O + HO2
!
      qqk(kloop,9)=qk(kloop,9)*y(kloop,19)*y(kloop,55)
!
!....         HO2 + NO = NO2 + OH
!
      qqk(kloop,10)=qk(kloop,10)*y(kloop,25)*y(kloop,51)
!
!....         HO2 + HO2 = H2O2
!
      qqk(kloop,11)=qk(kloop,11)*y(kloop,25)*y(kloop,25)
!
!....         H2 + OH = H2O + HO2
!
      qqk(kloop,12)=qk(kloop,12)*y(kloop,82)*y(kloop,55)
!
!....         CO + OH = HO2
!
      qqk(kloop,13)=qk(kloop,13)*y(kloop,10)*y(kloop,55)
!
!....         CH4 + OH = H2O + MO2
!
      qqk(kloop,14)=qk(kloop,14)*y(kloop,81)*y(kloop,55)
!
!....         MO2 + NO = CH2O + HO2 + NO2
!
      qqk(kloop,15)=qk(kloop,15)*y(kloop,43)*y(kloop,51)
!
!....         HO2 + MO2 = MP
!
      qqk(kloop,16)=qk(kloop,16)*y(kloop,25)*y(kloop,43)
!
!....         MO2 + MO2 = CH2O + MOH
!
      qqk(kloop,17)=qk(kloop,17)*y(kloop,43)*y(kloop,43)
!
!....         MO2 + MO2 = 2 CH2O + 2 HO2
!
      qqk(kloop,18)=qk(kloop,18)*y(kloop,43)*y(kloop,43)
!
!....         MP + OH = H2O + MO2
!
      qqk(kloop,19)=qk(kloop,19)*y(kloop,45)*y(kloop,55)
!
!....         MP + OH = CH2O + H2O + OH
!
      qqk(kloop,20)=qk(kloop,20)*y(kloop,45)*y(kloop,55)
!
!....         CH2O + OH = CO + H2O + HO2
!
      qqk(kloop,21)=qk(kloop,21)*y(kloop,9)*y(kloop,55)
!
!....         NO2 + OH = HNO3
!
      qqk(kloop,22)=qk(kloop,22)*y(kloop,52)*y(kloop,55)
!
!....         HNO3 + OH = H2O + NO3
!
      qqk(kloop,23)=qk(kloop,23)*y(kloop,23)*y(kloop,55)
!
!....         NO + OH = HNO2
!
      qqk(kloop,24)=qk(kloop,24)*y(kloop,51)*y(kloop,55)
!
!....         HNO2 + OH = H2O + NO2
!
      qqk(kloop,25)=qk(kloop,25)*y(kloop,22)*y(kloop,55)
!
!....         HO2 + NO2 = HNO4
!
      qqk(kloop,26)=qk(kloop,26)*y(kloop,25)*y(kloop,52)
!
!....         HNO4 = HO2 + NO2
!
      qqk(kloop,27)=qk(kloop,27)*y(kloop,24)
!
!....         HNO4 + OH = H2O + NO2
!
      qqk(kloop,28)=qk(kloop,28)*y(kloop,24)*y(kloop,55)
!
!....         HO2 + NO3 = NO2 + OH
!
      qqk(kloop,29)=qk(kloop,29)*y(kloop,25)*y(kloop,53)
!
!....         NO + NO3 = 2 NO2
!
      qqk(kloop,30)=qk(kloop,30)*y(kloop,51)*y(kloop,53)
!
!....         NO3 + OH = HO2 + NO2
!
      qqk(kloop,31)=qk(kloop,31)*y(kloop,53)*y(kloop,55)
!
!....         NO2 + NO3 = N2O5
!
      qqk(kloop,32)=qk(kloop,32)*y(kloop,52)*y(kloop,53)
!
!....         N2O5 = NO2 + NO3
!
      qqk(kloop,33)=qk(kloop,33)*y(kloop,50)
!
!....         HCOOH + OH = H2O + HO2
!
      qqk(kloop,34)=qk(kloop,34)*y(kloop,21)*y(kloop,55)
!
!....         MOH + OH = CH2O + HO2
!
      qqk(kloop,35)=qk(kloop,35)*y(kloop,44)*y(kloop,55)
!
!....         NO2 + NO3 = NO + NO2
!
      qqk(kloop,36)=qk(kloop,36)*y(kloop,52)*y(kloop,53)
!
!....         CH2O + NO3 = CO + HNO3 + HO2
!
      qqk(kloop,37)=qk(kloop,37)*y(kloop,9)*y(kloop,53)
!
!....         ALD2 + OH = H2O + MCO3
!
      qqk(kloop,38)=qk(kloop,38)*y(kloop,3)*y(kloop,55)
!
!....         ALD2 + NO3 = HNO3 + MCO3
!
      qqk(kloop,39)=qk(kloop,39)*y(kloop,3)*y(kloop,53)
!
!....         MCO3 + NO2 = PAN
!
      qqk(kloop,40)=qk(kloop,40)*y(kloop,40)*y(kloop,52)
!
!....         PAN = MCO3 + NO2
!
      qqk(kloop,41)=qk(kloop,41)*y(kloop,56)
!
!....         MCO3 + NO = MO2 + NO2
!
      qqk(kloop,42)=qk(kloop,42)*y(kloop,40)*y(kloop,51)
!
!....         C2H6 + OH = ETO2 + H2O
!
      qqk(kloop,43)=qk(kloop,43)*y(kloop,7)*y(kloop,55)
!
!....         ETO2 + NO = ALD2 + HO2 + NO2
!
      qqk(kloop,44)=qk(kloop,44)*y(kloop,12)*y(kloop,51)
!
!....         C3H8 + OH = B3O2
!
      qqk(kloop,45)=qk(kloop,45)*y(kloop,8)*y(kloop,55)
!
!....         C3H8 + OH = A3O2
!
      qqk(kloop,46)=qk(kloop,46)*y(kloop,8)*y(kloop,55)
!
!....         A3O2 + NO = HO2 + NO2 + RCHO
!
      qqk(kloop,47)=qk(kloop,47)*y(kloop,1)*y(kloop,51)
!
!....         NO + PO2 = ALD2 + CH2O + HO2 + NO2
!
      qqk(kloop,48)=qk(kloop,48)*y(kloop,51)*y(kloop,58)
!
!....         ALK4 + OH = R4O2
!
      qqk(kloop,49)=qk(kloop,49)*y(kloop,4)*y(kloop,55)
!
!....         NO + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O
!
      qqk(kloop,50)=qk(kloop,50)*y(kloop,51)*y(kloop,66)
!
!....         NO + R4O2 = R4N2
!
      qqk(kloop,51)=qk(kloop,51)*y(kloop,51)*y(kloop,66)
!
!....         NO + R4N1 =  0.75 ALD2 +  0.39 CH2O + 2 NO2 +  0.30 R4O2 +  
!
      qqk(kloop,52)=qk(kloop,52)*y(kloop,51)*y(kloop,64)
!
!....         ATO2 + NO =  0.19 CH2O +  0.77 HO2 +  0.19 MCO3 +  0.77 MGLY
!
      qqk(kloop,53)=qk(kloop,53)*y(kloop,5)*y(kloop,51)
!
!....         KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2 +  0.07 R4N2
!
      qqk(kloop,54)=qk(kloop,54)*y(kloop,34)*y(kloop,51)
!
!....         NO + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR
!
      qqk(kloop,55)=qk(kloop,55)*y(kloop,51)*y(kloop,74)
!
!....         NO + RIO2 = HNO3
!
      qqk(kloop,56)=qk(kloop,56)*y(kloop,51)*y(kloop,74)
!
!....         NO + RIO1 =  0.75 CH2O + HO2 + IALD + NO2
!
      qqk(kloop,57)=qk(kloop,57)*y(kloop,51)*y(kloop,73)
!
!....         NO + RIO1 = HNO3
!
      qqk(kloop,58)=qk(kloop,58)*y(kloop,51)*y(kloop,73)
!
!....         IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.24 GLYC +  0.17 GLYX 
!
      qqk(kloop,59)=qk(kloop,59)*y(kloop,27)*y(kloop,51)
!
!....         NO + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3
!
      qqk(kloop,60)=qk(kloop,60)*y(kloop,51)*y(kloop,78)
!
!....         NO + VRO2 = HNO3
!
      qqk(kloop,61)=qk(kloop,61)*y(kloop,51)*y(kloop,78)
!
!....         MRO2 + NO =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 
!
      qqk(kloop,62)=qk(kloop,62)*y(kloop,46)*y(kloop,51)
!
!....         MRO2 + NO = HNO3
!
      qqk(kloop,63)=qk(kloop,63)*y(kloop,46)*y(kloop,51)
!
!....         MVN2 + NO =  0.30 CH2O +  0.60 GLYC +  0.10 HNO3 +  0.30 HO2
!
      qqk(kloop,64)=qk(kloop,64)*y(kloop,49)*y(kloop,51)
!
!....         MAN2 + NO = CH2O + MGLY + 2 NO2
!
      qqk(kloop,65)=qk(kloop,65)*y(kloop,36)*y(kloop,51)
!
!....         B3O2 + NO = ACET + HO2 + NO2
!
      qqk(kloop,66)=qk(kloop,66)*y(kloop,6)*y(kloop,51)
!
!....         INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR
!
      qqk(kloop,67)=qk(kloop,67)*y(kloop,29)*y(kloop,51)
!
!....         NO + PRN1 = ALD2 + CH2O + 2 NO2
!
      qqk(kloop,68)=qk(kloop,68)*y(kloop,51)*y(kloop,61)
!
!....         ALK4 + NO3 = HNO3 + R4O2
!
      qqk(kloop,69)=qk(kloop,69)*y(kloop,4)*y(kloop,53)
!
!....         OH + R4N2 = H2O + R4N1
!
      qqk(kloop,70)=qk(kloop,70)*y(kloop,55)*y(kloop,65)
!
!....         ACTA + OH = H2O + MO2
!
      qqk(kloop,71)=qk(kloop,71)*y(kloop,2)*y(kloop,55)
!
!....         OH + RCHO = H2O + RCO3
!
      qqk(kloop,72)=qk(kloop,72)*y(kloop,55)*y(kloop,70)
!
!....         NO2 + RCO3 = PPN
!
      qqk(kloop,73)=qk(kloop,73)*y(kloop,52)*y(kloop,71)
!
!....         GCO3 + NO2 = GPAN
!
      qqk(kloop,74)=qk(kloop,74)*y(kloop,14)*y(kloop,52)
!
!....         MAO3 + NO2 = PMN
!
      qqk(kloop,75)=qk(kloop,75)*y(kloop,37)*y(kloop,52)
!
!....         PPN = NO2 + RCO3
!
      qqk(kloop,76)=qk(kloop,76)*y(kloop,60)
!
!....         GPAN = GCO3 + NO2
!
      qqk(kloop,77)=qk(kloop,77)*y(kloop,18)
!
!....         NO + RCO3 = ETO2 + NO2
!
      qqk(kloop,78)=qk(kloop,78)*y(kloop,51)*y(kloop,71)
!
!....         GCO3 + NO = CH2O + HO2 + NO2
!
      qqk(kloop,79)=qk(kloop,79)*y(kloop,14)*y(kloop,51)
!
!....         MAO3 + NO = 4 CH2O + HO2 + NO2
!
      qqk(kloop,80)=qk(kloop,80)*y(kloop,37)*y(kloop,51)
!
!....         NO3 + RCHO = HNO3 + RCO3
!
      qqk(kloop,81)=qk(kloop,81)*y(kloop,53)*y(kloop,70)
!
!....         ACET + OH = ATO2 + H2O
!
      qqk(kloop,82)=qk(kloop,82)*y(kloop,80)*y(kloop,55)
!
!....         A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.
!
      qqk(kloop,83)=qk(kloop,83)*y(kloop,1)*y(kloop,43)
!
!....         MO2 + PO2 =  0.50 ALD2 +  1.25 CH2O +  0.16 HAC + HO2 +  0.2
!
      qqk(kloop,84)=qk(kloop,84)*y(kloop,43)*y(kloop,58)
!
!....         HO2 + R4O2 = R4P
!
      qqk(kloop,85)=qk(kloop,85)*y(kloop,25)*y(kloop,66)
!
!....         HO2 + R4N1 = R4N2
!
      qqk(kloop,86)=qk(kloop,86)*y(kloop,25)*y(kloop,64)
!
!....         ATO2 + HO2 = MCO3 + MO2
!
      qqk(kloop,87)=qk(kloop,87)*y(kloop,5)*y(kloop,25)
!
!....         HO2 + KO2 = MGLY + MO2
!
      qqk(kloop,88)=qk(kloop,88)*y(kloop,25)*y(kloop,34)
!
!....         HO2 + RIO2 = RIP
!
      qqk(kloop,89)=qk(kloop,89)*y(kloop,25)*y(kloop,74)
!
!....         HO2 + RIO1 = RIP
!
      qqk(kloop,90)=qk(kloop,90)*y(kloop,25)*y(kloop,73)
!
!....         HO2 + IAO2 = IAP
!
      qqk(kloop,91)=qk(kloop,91)*y(kloop,25)*y(kloop,27)
!
!....         HO2 + ISN1 = ISNP
!
      qqk(kloop,92)=qk(kloop,92)*y(kloop,25)*y(kloop,31)
!
!....         HO2 + VRO2 = VRP
!
      qqk(kloop,93)=qk(kloop,93)*y(kloop,25)*y(kloop,78)
!
!....         HO2 + MRO2 = MRP
!
      qqk(kloop,94)=qk(kloop,94)*y(kloop,25)*y(kloop,46)
!
!....         HO2 + MVN2 = ISNP
!
      qqk(kloop,95)=qk(kloop,95)*y(kloop,25)*y(kloop,49)
!
!....         HO2 + MAN2 = ISNP
!
      qqk(kloop,96)=qk(kloop,96)*y(kloop,25)*y(kloop,36)
!
!....         B3O2 + HO2 = RB3P
!
      qqk(kloop,97)=qk(kloop,97)*y(kloop,6)*y(kloop,25)
!
!....         HO2 + INO2 = INPN
!
      qqk(kloop,98)=qk(kloop,98)*y(kloop,25)*y(kloop,29)
!
!....         HO2 + PRN1 = PRPN
!
      qqk(kloop,99)=qk(kloop,99)*y(kloop,25)*y(kloop,61)
!
!....         MEK + OH = H2O + KO2
!
      qqk(kloop,100)=qk(kloop,100)*y(kloop,41)*y(kloop,55)
!
!....         ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O +  0.25 EOH + HO2 +  0.
!
      qqk(kloop,101)=qk(kloop,101)*y(kloop,12)*y(kloop,43)
!
!....         MEK + NO3 = HNO3 + KO2
!
      qqk(kloop,102)=qk(kloop,102)*y(kloop,41)*y(kloop,53)
!
!....         MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3
!
      qqk(kloop,103)=qk(kloop,103)*y(kloop,43)*y(kloop,66)
!
!....         MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH
!
      qqk(kloop,104)=qk(kloop,104)*y(kloop,43)*y(kloop,64)
!
!....         ATO2 + MO2 =  0.85 CH2O +  0.90 HO2 +  0.10 MCO3 +  0.25 MEK
!
      qqk(kloop,105)=qk(kloop,105)*y(kloop,5)*y(kloop,43)
!
!....         KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3
!
      qqk(kloop,106)=qk(kloop,106)*y(kloop,34)*y(kloop,43)
!
!....         MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.06 IALD +  0.14 MAC
!
      qqk(kloop,107)=qk(kloop,107)*y(kloop,43)*y(kloop,74)
!
!....         MO2 + RIO1 =  1.13 CH2O + HO2 +  0.50 IALD +  0.25 MEK +  0.
!
      qqk(kloop,108)=qk(kloop,108)*y(kloop,43)*y(kloop,73)
!
!....         IAO2 + MO2 =  0.95 CH2O +  0.15 CO +  0.13 GLYC +  0.09 GLYX
!
      qqk(kloop,109)=qk(kloop,109)*y(kloop,27)*y(kloop,43)
!
!....         ISN1 + MO2 =  0.75 CH2O +  0.50 GLYC +  0.50 HAC +  0.50 HO2
!
      qqk(kloop,110)=qk(kloop,110)*y(kloop,31)*y(kloop,43)
!
!....         MO2 + VRO2 =  0.89 CH2O +  0.36 GLYC +  0.64 HO2 +  0.36 MCO
!
      qqk(kloop,111)=qk(kloop,111)*y(kloop,43)*y(kloop,78)
!
!....         MO2 + MRO2 =  0.84 CH2O +  0.42 CO +  0.42 HAC + HO2 +  0.25
!
      qqk(kloop,112)=qk(kloop,112)*y(kloop,43)*y(kloop,46)
!
!....         MO2 + MVN2 =  1.25 CH2O +  0.75 HO2 +  0.25 MCO3 +  0.25 MGL
!
      qqk(kloop,113)=qk(kloop,113)*y(kloop,43)*y(kloop,49)
!
!....         MAN2 + MO2 =  1.25 CH2O +  0.50 HO2 +  0.50 MGLY +  0.25 MOH
!
      qqk(kloop,114)=qk(kloop,114)*y(kloop,36)*y(kloop,43)
!
!....         B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.
!
      qqk(kloop,115)=qk(kloop,115)*y(kloop,6)*y(kloop,43)
!
!....         INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MAC
!
      qqk(kloop,116)=qk(kloop,116)*y(kloop,29)*y(kloop,43)
!
!....         MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH
!
      qqk(kloop,117)=qk(kloop,117)*y(kloop,43)*y(kloop,61)
!
!....         EOH + OH = ALD2 + HO2
!
      qqk(kloop,118)=qk(kloop,118)*y(kloop,11)*y(kloop,55)
!
!....         OH + ROH = HO2 + RCHO
!
      qqk(kloop,119)=qk(kloop,119)*y(kloop,55)*y(kloop,76)
!
!....         ETO2 + ETO2 = 2 ALD2 + 2 HO2
!
      qqk(kloop,120)=qk(kloop,120)*y(kloop,12)*y(kloop,12)
!
!....         ETO2 + ETO2 = ALD2 + EOH
!
      qqk(kloop,121)=qk(kloop,121)*y(kloop,12)*y(kloop,12)
!
!....         ETO2 + HO2 = ETP
!
      qqk(kloop,122)=qk(kloop,122)*y(kloop,12)*y(kloop,25)
!
!....         A3O2 + HO2 = RA3P
!
      qqk(kloop,123)=qk(kloop,123)*y(kloop,1)*y(kloop,25)
!
!....         HO2 + PO2 = PP
!
      qqk(kloop,124)=qk(kloop,124)*y(kloop,25)*y(kloop,58)
!
!....         HO2 + MCO3 = ACTA + O3
!
      qqk(kloop,125)=qk(kloop,125)*y(kloop,25)*y(kloop,40)
!
!....         HO2 + MCO3 = MAP
!
      qqk(kloop,126)=qk(kloop,126)*y(kloop,25)*y(kloop,40)
!
!....         HO2 + RCO3 =  0.30 O3 +  0.30 RCOOH +  0.70 RP
!
      qqk(kloop,127)=qk(kloop,127)*y(kloop,25)*y(kloop,71)
!
!....         GCO3 + HO2 =  0.70 GP +  0.30 O3 +  0.30 RCOOH
!
      qqk(kloop,128)=qk(kloop,128)*y(kloop,14)*y(kloop,25)
!
!....         HO2 + MAO3 =  0.70 MAOP +  0.30 O3 +  0.30 RCOOH
!
      qqk(kloop,129)=qk(kloop,129)*y(kloop,25)*y(kloop,37)
!
!....         OH + PRPE = PO2
!
      qqk(kloop,130)=qk(kloop,130)*y(kloop,55)*y(kloop,62)
!
!....         O3 + PRPE =  0.50 ALD2 +  0.54 CH2O +  0.42 CO +  0.06 H2 + 
!
      qqk(kloop,131)=qk(kloop,131)*y(kloop,54)*y(kloop,62)
!
!....         PMN = MAO3 + NO2
!
      qqk(kloop,132)=qk(kloop,132)*y(kloop,57)
!
!....         OH + PMN =  2.23 CH2O +  0.59 HAC + 2 HO2 + NO2
!
      qqk(kloop,133)=qk(kloop,133)*y(kloop,55)*y(kloop,57)
!
!....         O3 + PMN =  0.60 CH2O + HO2 + NO2
!
      qqk(kloop,134)=qk(kloop,134)*y(kloop,54)*y(kloop,57)
!
!....         GLYC + OH =  0.80 GCO3 +  0.20 GLYX +  0.20 HO2
!
      qqk(kloop,135)=qk(kloop,135)*y(kloop,15)*y(kloop,55)
!
!....         NO3 + PRPE = PRN1
!
      qqk(kloop,136)=qk(kloop,136)*y(kloop,53)*y(kloop,62)
!
!....         GLYX + OH = 2 CO + HO2
!
      qqk(kloop,137)=qk(kloop,137)*y(kloop,16)*y(kloop,55)
!
!....         MGLY + OH = CO + MCO3
!
      qqk(kloop,138)=qk(kloop,138)*y(kloop,42)*y(kloop,55)
!
!....         GLYX + NO3 = 2 CO + HNO3 + HO2
!
      qqk(kloop,139)=qk(kloop,139)*y(kloop,16)*y(kloop,53)
!
!....         MGLY + NO3 = CO + HNO3 + MCO3
!
      qqk(kloop,140)=qk(kloop,140)*y(kloop,42)*y(kloop,53)
!
!....         ISOP + OH = RIO2
!
      qqk(kloop,141)=qk(kloop,141)*y(kloop,33)*y(kloop,55)
!
!....         MVK + OH = VRO2
!
      qqk(kloop,142)=qk(kloop,142)*y(kloop,48)*y(kloop,55)
!
!....         MACR + OH =  0.50 MAO3 +  0.50 MRO2
!
      qqk(kloop,143)=qk(kloop,143)*y(kloop,35)*y(kloop,55)
!
!....         HAC + OH = HO2 + MGLY
!
      qqk(kloop,144)=qk(kloop,144)*y(kloop,20)*y(kloop,55)
!
!....         A3O2 + MCO3 = HO2 + MO2 + RCHO
!
      qqk(kloop,145)=qk(kloop,145)*y(kloop,1)*y(kloop,40)
!
!....         MCO3 + PO2 = ALD2 + CH2O + HO2 + MO2
!
      qqk(kloop,146)=qk(kloop,146)*y(kloop,40)*y(kloop,58)
!
!....         A3O2 + MCO3 = ACTA + RCHO
!
      qqk(kloop,147)=qk(kloop,147)*y(kloop,1)*y(kloop,40)
!
!....         MCO3 + PO2 = ACTA +  0.65 HAC +  0.35 RCHO
!
      qqk(kloop,148)=qk(kloop,148)*y(kloop,40)*y(kloop,58)
!
!....         ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.39 MACR +
!
      qqk(kloop,149)=qk(kloop,149)*y(kloop,33)*y(kloop,54)
!
!....         MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 + 
!
      qqk(kloop,150)=qk(kloop,150)*y(kloop,48)*y(kloop,54)
!
!....         MACR + O3 =  0.70 CH2O +  0.20 CO +  0.28 HO2 +  0.80 MGLY +
!
      qqk(kloop,151)=qk(kloop,151)*y(kloop,35)*y(kloop,54)
!
!....         ISOP + NO3 = INO2
!
      qqk(kloop,152)=qk(kloop,152)*y(kloop,33)*y(kloop,53)
!
!....         MVK + NO3 = MVN2
!
      qqk(kloop,153)=qk(kloop,153)*y(kloop,48)*y(kloop,53)
!
!....         MACR + NO3 = MAN2
!
      qqk(kloop,154)=qk(kloop,154)*y(kloop,35)*y(kloop,53)
!
!....         MACR + NO3 = HNO3 + MAO3
!
      qqk(kloop,155)=qk(kloop,155)*y(kloop,35)*y(kloop,53)
!
!....         MO2 + RCO3 = CH2O + ETO2 + HO2
!
      qqk(kloop,156)=qk(kloop,156)*y(kloop,43)*y(kloop,71)
!
!....         GCO3 + MO2 = 2 CH2O + 2 HO2
!
      qqk(kloop,157)=qk(kloop,157)*y(kloop,14)*y(kloop,43)
!
!....         MAO3 + MO2 = 2 CH2O + HO2 + MCO3
!
      qqk(kloop,158)=qk(kloop,158)*y(kloop,37)*y(kloop,43)
!
!....         MO2 + RCO3 = CH2O + RCOOH
!
      qqk(kloop,159)=qk(kloop,159)*y(kloop,43)*y(kloop,71)
!
!....         GCO3 + MO2 = CH2O + RCOOH
!
      qqk(kloop,160)=qk(kloop,160)*y(kloop,14)*y(kloop,43)
!
!....         MAO3 + MO2 = CH2O + RCOOH
!
      qqk(kloop,161)=qk(kloop,161)*y(kloop,37)*y(kloop,43)
!
!....         INPN + OH = INO2
!
      qqk(kloop,162)=qk(kloop,162)*y(kloop,30)*y(kloop,55)
!
!....         OH + PRPN = PRN1
!
      qqk(kloop,163)=qk(kloop,163)*y(kloop,55)*y(kloop,63)
!
!....         ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH
!
      qqk(kloop,164)=qk(kloop,164)*y(kloop,13)*y(kloop,55)
!
!....         OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,165)=qk(kloop,165)*y(kloop,55)*y(kloop,68)
!
!....         OH + RB3P =  0.50 B3O2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,166)=qk(kloop,166)*y(kloop,55)*y(kloop,69)
!
!....         OH + R4P =  0.50 OH +  0.50 R4O2 +  0.50 RCHO
!
      qqk(kloop,167)=qk(kloop,167)*y(kloop,55)*y(kloop,67)
!
!....         OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3
!
      qqk(kloop,168)=qk(kloop,168)*y(kloop,55)*y(kloop,77)
!
!....         OH + PP =  0.50 OH +  0.50 PO2 +  0.50 RCHO
!
      qqk(kloop,169)=qk(kloop,169)*y(kloop,55)*y(kloop,59)
!
!....         GP + OH =  0.50 CH2O +  0.50 GCO3 +  0.50 OH
!
      qqk(kloop,170)=qk(kloop,170)*y(kloop,17)*y(kloop,55)
!
!....         OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2
!
      qqk(kloop,171)=qk(kloop,171)*y(kloop,55)*y(kloop,75)
!
!....         IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,172)=qk(kloop,172)*y(kloop,28)*y(kloop,55)
!
!....         ISNP + OH =  0.50 ISN1 +  0.50 NO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,173)=qk(kloop,173)*y(kloop,32)*y(kloop,55)
!
!....         OH + VRP =  0.50 OH +  0.50 RCHO +  0.50 VRO2
!
      qqk(kloop,174)=qk(kloop,174)*y(kloop,55)*y(kloop,79)
!
!....         MRP + OH =  0.50 MRO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,175)=qk(kloop,175)*y(kloop,47)*y(kloop,55)
!
!....         MAOP + OH =  0.50 MAO3 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,176)=qk(kloop,176)*y(kloop,38)*y(kloop,55)
!
!....         MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH
!
      qqk(kloop,177)=qk(kloop,177)*y(kloop,39)*y(kloop,55)
!
!....         C2H6 + NO3 = ETO2 + HNO3
!
      qqk(kloop,178)=qk(kloop,178)*y(kloop,7)*y(kloop,53)
!
!....         IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3
!
      qqk(kloop,179)=qk(kloop,179)*y(kloop,26)*y(kloop,55)
!
!....         IALD + O3 =  0.12 CH2O +  0.28 GLYC +  0.20 GLYX +  0.20 HAC
!
      qqk(kloop,180)=qk(kloop,180)*y(kloop,26)*y(kloop,54)
!
!....         MCO3 + MCO3 = 2 MO2
!
      qqk(kloop,181)=qk(kloop,181)*y(kloop,40)*y(kloop,40)
!
!....         MCO3 + MO2 = CH2O + HO2 + MO2
!
      qqk(kloop,182)=qk(kloop,182)*y(kloop,40)*y(kloop,43)
!
!....         MCO3 + MO2 = ACTA + CH2O
!
      qqk(kloop,183)=qk(kloop,183)*y(kloop,40)*y(kloop,43)
!
!....         MCO3 + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B
!
      qqk(kloop,184)=qk(kloop,184)*y(kloop,40)*y(kloop,66)
!
!....         ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 +  0.80 MG
!
      qqk(kloop,185)=qk(kloop,185)*y(kloop,5)*y(kloop,40)
!
!....         KO2 + MCO3 = ALD2 + MCO3 + MO2
!
      qqk(kloop,186)=qk(kloop,186)*y(kloop,34)*y(kloop,40)
!
!....         MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MA
!
      qqk(kloop,187)=qk(kloop,187)*y(kloop,40)*y(kloop,74)
!
!....         MCO3 + RIO1 =  0.75 CH2O + HO2 + IALD + MO2
!
      qqk(kloop,188)=qk(kloop,188)*y(kloop,40)*y(kloop,73)
!
!....         IAO2 + MCO3 =  0.40 CH2O +  0.29 CO +  0.26 GLYC +  0.18 GLY
!
      qqk(kloop,189)=qk(kloop,189)*y(kloop,27)*y(kloop,40)
!
!....         ISN1 + MCO3 = GLYC + HAC + MO2 + NO2
!
      qqk(kloop,190)=qk(kloop,190)*y(kloop,31)*y(kloop,40)
!
!....         MCO3 + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MC
!
      qqk(kloop,191)=qk(kloop,191)*y(kloop,40)*y(kloop,78)
!
!....         MCO3 + MRO2 =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.1
!
      qqk(kloop,192)=qk(kloop,192)*y(kloop,40)*y(kloop,46)
!
!....         B3O2 + MCO3 = ACET + HO2 + MO2
!
      qqk(kloop,193)=qk(kloop,193)*y(kloop,6)*y(kloop,40)
!
!....         MCO3 + R4N1 =  0.75 ALD2 +  0.39 CH2O + MO2 + NO2 +  0.30 R4
!
      qqk(kloop,194)=qk(kloop,194)*y(kloop,40)*y(kloop,64)
!
!....         MCO3 + MVN2 = CH2O +  0.50 HO2 +  0.50 MCO3 +  0.50 MGLY + M
!
      qqk(kloop,195)=qk(kloop,195)*y(kloop,40)*y(kloop,49)
!
!....         MAN2 + MCO3 = CH2O + MGLY + MO2 + NO2
!
      qqk(kloop,196)=qk(kloop,196)*y(kloop,36)*y(kloop,40)
!
!....         INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MA
!
      qqk(kloop,197)=qk(kloop,197)*y(kloop,29)*y(kloop,40)
!
!....         MCO3 + PRN1 = ALD2 + CH2O + MO2 + NO2
!
      qqk(kloop,198)=qk(kloop,198)*y(kloop,40)*y(kloop,61)
!
!....         MCO3 + R4O2 = ACTA + MEK
!
      qqk(kloop,199)=qk(kloop,199)*y(kloop,40)*y(kloop,66)
!
!....         ATO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,200)=qk(kloop,200)*y(kloop,5)*y(kloop,40)
!
!....         KO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,201)=qk(kloop,201)*y(kloop,34)*y(kloop,40)
!
!....         MCO3 + RIO2 = ACTA + MEK
!
      qqk(kloop,202)=qk(kloop,202)*y(kloop,40)*y(kloop,74)
!
!....         MCO3 + RIO1 = ACTA + MEK
!
      qqk(kloop,203)=qk(kloop,203)*y(kloop,40)*y(kloop,73)
!
!....         IAO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,204)=qk(kloop,204)*y(kloop,27)*y(kloop,40)
!
!....         MCO3 + VRO2 = ACTA + MEK
!
      qqk(kloop,205)=qk(kloop,205)*y(kloop,40)*y(kloop,78)
!
!....         MCO3 + MRO2 = ACTA + MEK
!
      qqk(kloop,206)=qk(kloop,206)*y(kloop,40)*y(kloop,46)
!
!....         MCO3 + R4N1 = ACTA + NO2 + RCHO
!
      qqk(kloop,207)=qk(kloop,207)*y(kloop,40)*y(kloop,64)
!
!....         ISN1 + MCO3 = ACTA + NO2 + RCHO
!
      qqk(kloop,208)=qk(kloop,208)*y(kloop,31)*y(kloop,40)
!
!....         MCO3 + MVN2 = ACTA + NO2 + RCHO
!
      qqk(kloop,209)=qk(kloop,209)*y(kloop,40)*y(kloop,49)
!
!....         MAN2 + MCO3 = ACTA + NO2 + RCHO
!
      qqk(kloop,210)=qk(kloop,210)*y(kloop,36)*y(kloop,40)
!
!....         INO2 + MCO3 = ACTA + NO2 + RCHO
!
      qqk(kloop,211)=qk(kloop,211)*y(kloop,29)*y(kloop,40)
!
!....         MCO3 + PRN1 = ACTA + NO2 + RCHO
!
      qqk(kloop,212)=qk(kloop,212)*y(kloop,40)*y(kloop,61)
!
!....         B3O2 + MCO3 = ACET + ACTA
!
      qqk(kloop,213)=qk(kloop,213)*y(kloop,6)*y(kloop,40)
!
!....         ETO2 + MCO3 = ALD2 + HO2 + MO2
!
      qqk(kloop,214)=qk(kloop,214)*y(kloop,12)*y(kloop,40)
!
!....         ETO2 + MCO3 = ACTA + ALD2
!
      qqk(kloop,215)=qk(kloop,215)*y(kloop,12)*y(kloop,40)
!
!....         MCO3 + RCO3 = ETO2 + MO2
!
      qqk(kloop,216)=qk(kloop,216)*y(kloop,40)*y(kloop,71)
!
!....         GCO3 + MCO3 = CH2O + HO2 + MO2
!
      qqk(kloop,217)=qk(kloop,217)*y(kloop,14)*y(kloop,40)
!
!....         MAO3 + MCO3 = CH2O + MCO3 + MO2
!
      qqk(kloop,218)=qk(kloop,218)*y(kloop,37)*y(kloop,40)
!
!....         NO3 + NO3 = 2 NO2
!
      qqk(kloop,219)=qk(kloop,219)*y(kloop,53)*y(kloop,53)
!
!....         HO2 =  0.50 H2O2
!
      qqk(kloop,220)=qk(kloop,220)*y(kloop,25)
!
!....         NO2 =  0.50 HNO2 +  0.50 HNO3
!
      qqk(kloop,221)=qk(kloop,221)*y(kloop,52)
!
!....         NO3 = HNO3
!
      qqk(kloop,222)=qk(kloop,222)*y(kloop,53)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,223)=qk(kloop,223)*y(kloop,50)
!
!                  start photolytic reactions
!
!....  O3 + hv = 2 OH
!
      qqj(kloop,1)=qj(kloop,1)*y(kloop,54)
!
!....  NO2 + hv = NO + O3
!
      qqj(kloop,2)=qj(kloop,2)*y(kloop,52)
!
!....  H2O2 + hv = 2 OH
!
      qqj(kloop,3)=qj(kloop,3)*y(kloop,19)
!
!....  MP + hv = CH2O + HO2 + OH
!
      qqj(kloop,4)=qj(kloop,4)*y(kloop,45)
!
!....  CH2O + hv = CO + 2 HO2
!
      qqj(kloop,5)=qj(kloop,5)*y(kloop,9)
!
!....  CH2O + hv = CO + H2
!
      qqj(kloop,6)=qj(kloop,6)*y(kloop,9)
!
!....  HNO3 + hv = NO2 + OH
!
      qqj(kloop,7)=qj(kloop,7)*y(kloop,23)
!
!....  HNO2 + hv = NO + OH
!
      qqj(kloop,8)=qj(kloop,8)*y(kloop,22)
!
!....  HNO4 + hv = NO3 + OH
!
      qqj(kloop,9)=qj(kloop,9)*y(kloop,24)
!
!....  NO3 + hv = NO2 + O3
!
      qqj(kloop,10)=qj(kloop,10)*y(kloop,53)
!
!....  NO3 + hv = NO
!
      qqj(kloop,11)=qj(kloop,11)*y(kloop,53)
!
!....  N2O5 + hv = NO2 + NO3
!
      qqj(kloop,12)=qj(kloop,12)*y(kloop,50)
!
!....  N2O5 + hv = NO + NO3 + O3
!
      qqj(kloop,13)=qj(kloop,13)*y(kloop,50)
!
!....  HNO4 + hv = HO2 + NO2
!
      qqj(kloop,14)=qj(kloop,14)*y(kloop,24)
!
!....  ALD2 + hv = CO + HO2 + MO2
!
      qqj(kloop,15)=qj(kloop,15)*y(kloop,3)
!
!....  ALD2 + hv = CH4 + CO
!
      qqj(kloop,16)=qj(kloop,16)*y(kloop,3)
!
!....  PAN + hv = MCO3 + NO2
!
      qqj(kloop,17)=qj(kloop,17)*y(kloop,56)
!
!....  RCHO + hv = CO + ETO2 + HO2
!
      qqj(kloop,18)=qj(kloop,18)*y(kloop,70)
!
!....  ACET + hv = MCO3 + MO2
!
      qqj(kloop,19)=qj(kloop,19)*y(kloop,80)
!
!....  MEK + hv = ETO2 + MCO3
!
      qqj(kloop,20)=qj(kloop,20)*y(kloop,41)
!
!....  GLYC + hv = CH2O + CO + 2 HO2
!
      qqj(kloop,21)=qj(kloop,21)*y(kloop,15)
!
!....  GLYX + hv = 2 CO + H2
!
      qqj(kloop,22)=qj(kloop,22)*y(kloop,16)
!
!....  GLYX + hv = 2 CO + 2 HO2
!
      qqj(kloop,23)=qj(kloop,23)*y(kloop,16)
!
!....  GLYX + hv = CH2O + CO
!
      qqj(kloop,24)=qj(kloop,24)*y(kloop,16)
!
!....  MGLY + hv = CO + HO2 + MCO3
!
      qqj(kloop,25)=qj(kloop,25)*y(kloop,42)
!
!....  MGLY + hv = ALD2 + CO
!
      qqj(kloop,26)=qj(kloop,26)*y(kloop,42)
!
!....  MVK + hv = CO + PRPE
!
      qqj(kloop,27)=qj(kloop,27)*y(kloop,48)
!
!....  MVK + hv = CH2O + CO + HO2 + MCO3
!
      qqj(kloop,28)=qj(kloop,28)*y(kloop,48)
!
!....  MVK + hv = MAO3 + MO2
!
      qqj(kloop,29)=qj(kloop,29)*y(kloop,48)
!
!....  MACR + hv = HO2 + MAO3
!
      qqj(kloop,30)=qj(kloop,30)*y(kloop,35)
!
!....  MACR + hv =  0.20 CH2O + CO +  1.80 HO2 +  0.20 MCO3 +  0.80
!
      qqj(kloop,31)=qj(kloop,31)*y(kloop,35)
!
!....  HAC + hv = CH2O + HO2 + MCO3
!
      qqj(kloop,32)=qj(kloop,32)*y(kloop,20)
!
!....  INPN + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,33)=qj(kloop,33)*y(kloop,30)
!
!....  PRPN + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,34)=qj(kloop,34)*y(kloop,63)
!
!....  ETP + hv = ALD2 + HO2 + OH
!
      qqj(kloop,35)=qj(kloop,35)*y(kloop,13)
!
!....  RA3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,36)=qj(kloop,36)*y(kloop,68)
!
!....  RB3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,37)=qj(kloop,37)*y(kloop,69)
!
!....  R4P + hv = HO2 + OH + RCHO
!
      qqj(kloop,38)=qj(kloop,38)*y(kloop,67)
!
!....  PP + hv = HO2 + OH + RCHO
!
      qqj(kloop,39)=qj(kloop,39)*y(kloop,59)
!
!....  RP + hv = ALD2 + HO2 + OH
!
      qqj(kloop,40)=qj(kloop,40)*y(kloop,77)
!
!....  GP + hv = CH2O + HO2 + OH
!
      qqj(kloop,41)=qj(kloop,41)*y(kloop,17)
!
!....  RIP + hv =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR 
!
      qqj(kloop,42)=qj(kloop,42)*y(kloop,75)
!
!....  IAP + hv =  0.67 CO +  0.26 GLYC +  0.19 H2 +  0.36 HAC + HO
!
      qqj(kloop,43)=qj(kloop,43)*y(kloop,28)
!
!....  ISNP + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,44)=qj(kloop,44)*y(kloop,32)
!
!....  VRP + hv =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 
!
      qqj(kloop,45)=qj(kloop,45)*y(kloop,79)
!
!....  MRP + hv =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 M
!
      qqj(kloop,46)=qj(kloop,46)*y(kloop,47)
!
!....  MAOP + hv = HO2 + OH + RCHO
!
      qqj(kloop,47)=qj(kloop,47)*y(kloop,38)
!
!....  R4N2 + hv =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O
!
      qqj(kloop,48)=qj(kloop,48)*y(kloop,65)
!
!....  MAP + hv = MO2 + OH
!
      qqj(kloop,49)=qj(kloop,49)*y(kloop,39)
      enddo
!
!.... End of kloop
!
      return
      end
