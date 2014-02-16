!=======================================================================
!
! $Id: setkin_ratecalc.F90,v 1.2 2011-08-09 22:12:57 mrdamon Exp $
!
! ROUTINE
!   Calc_rate_Setkin - IMPACT version (setkin_ratecalc.F)
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
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!
!=======================================================================
      subroutine Calc_rate_Setkin  &
     &  (nblock, numj, numk, numspc,  &
     &   qk, qj, y, qqk, qqj)

      implicit none

      integer  &
     &  nblock,  &
     &  numj, numk, numspc

      real*8  &
     &  qj    (nblock, numj),  &
     &  qk    (nblock, numk),  &
     &  qqj   (nblock, numj),  &
     &  qqk   (nblock, numk),  &
     &  y     (nblock, numspc)

      integer  &
     &  kloop

      do kloop = 1, nblock
!
!                  start thermal reactions
!
!....         O + O2 = O3
!
      qqk(kloop,1)=qk(kloop,1)*y(kloop,1)
!
!....         O + O3 = 2 O2
!
      qqk(kloop,2)=qk(kloop,2)*y(kloop,1)*y(kloop,3)
!
!....         N2 + O1D = N2 + O
!
      qqk(kloop,3)=qk(kloop,3)*y(kloop,2)
!
!....         O1D + O2 = O + O2
!
      qqk(kloop,4)=qk(kloop,4)*y(kloop,2)
!
!....         O1D + O3 = 2 O2
!
      qqk(kloop,5)=qk(kloop,5)*y(kloop,2)*y(kloop,3)
!
!....         H2 + O1D = H + OH
!
      qqk(kloop,6)=qk(kloop,6)*y(kloop,5)*y(kloop,2)
!
!....         H2O + O1D = 2 OH
!
      qqk(kloop,7)=qk(kloop,7)*y(kloop,8)*y(kloop,2)
!
!....         H + O2 = HO2
!
      qqk(kloop,8)=qk(kloop,8)*y(kloop,4)
!
!....         H + O3 = O2 + OH
!
      qqk(kloop,9)=qk(kloop,9)*y(kloop,4)*y(kloop,3)
!
!....         H + HO2 = 2 OH
!
      qqk(kloop,10)=qk(kloop,10)*y(kloop,4)*y(kloop,7)
!
!....         O + OH = H + O2
!
      qqk(kloop,11)=qk(kloop,11)*y(kloop,1)*y(kloop,6)
!
!....         O3 + OH = HO2 + O2
!
      qqk(kloop,12)=qk(kloop,12)*y(kloop,3)*y(kloop,6)
!
!....         H2 + OH = H + H2O
!
      qqk(kloop,13)=qk(kloop,13)*y(kloop,5)*y(kloop,6)
!
!....         OH + OH = H2O + O
!
      qqk(kloop,14)=qk(kloop,14)*y(kloop,6)*y(kloop,6)
!
!....         HO2 + O = O2 + OH
!
      qqk(kloop,15)=qk(kloop,15)*y(kloop,7)*y(kloop,1)
!
!....         HO2 + O3 = 2 O2 + OH
!
      qqk(kloop,16)=qk(kloop,16)*y(kloop,7)*y(kloop,3)
!
!....         HO2 + OH = H2O + O2
!
      qqk(kloop,17)=qk(kloop,17)*y(kloop,7)*y(kloop,6)
!
!....         HO2 + HO2 = H2O2 + O2
!
      qqk(kloop,18)=qk(kloop,18)*y(kloop,7)*y(kloop,7)
!
!....         H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      qqk(kloop,19)=qk(kloop,19)*y(kloop,8)*y(kloop,7)*y(kloop,7)
!
!....         H2O2 + OH = H2O + HO2
!
      qqk(kloop,20)=qk(kloop,20)*y(kloop,9)*y(kloop,6)
!
!....         N2O + O1D = N2 + O2
!
      qqk(kloop,21)=qk(kloop,21)*y(kloop,11)*y(kloop,2)
!
!....         N2O + O1D = 2 NO
!
      qqk(kloop,22)=qk(kloop,22)*y(kloop,11)*y(kloop,2)
!
!....         N + O2 = NO + O
!
      qqk(kloop,23)=qk(kloop,23)*y(kloop,10)
!
!....         NO + O3 = NO2 + O2
!
      qqk(kloop,24)=qk(kloop,24)*y(kloop,12)*y(kloop,3)
!
!....         NO2 + O = NO + O2
!
      qqk(kloop,25)=qk(kloop,25)*y(kloop,13)*y(kloop,1)
!
!....         NO2 + O3 = NO3 + O2
!
      qqk(kloop,26)=qk(kloop,26)*y(kloop,13)*y(kloop,3)
!
!....         NO + OH = HONO
!
      qqk(kloop,27)=qk(kloop,27)*y(kloop,12)*y(kloop,6)
!
!....         HO2 + NO = NO2 + OH
!
      qqk(kloop,28)=qk(kloop,28)*y(kloop,7)*y(kloop,12)
!
!....         NO2 + OH = HNO3
!
      qqk(kloop,29)=qk(kloop,29)*y(kloop,13)*y(kloop,6)
!
!....         HO2 + NO2 = HO2NO2
!
      qqk(kloop,30)=qk(kloop,30)*y(kloop,7)*y(kloop,13)
!
!....         H2O + N2O5 = 2 HNO3
!
      qqk(kloop,31)=qk(kloop,31)*y(kloop,8)*y(kloop,15)
!
!....         HO2NO2 = HO2 + NO2
!
      qqk(kloop,32)=qk(kloop,32)*y(kloop,18)
!
!....         HONO + OH = H2O + NO2
!
      qqk(kloop,33)=qk(kloop,33)*y(kloop,16)*y(kloop,6)
!
!....         HNO3 + OH = H2O + NO3
!
      qqk(kloop,34)=qk(kloop,34)*y(kloop,17)*y(kloop,6)
!
!....         HO2NO2 + OH = H2O + NO2 + O2
!
      qqk(kloop,35)=qk(kloop,35)*y(kloop,18)*y(kloop,6)
!
!....         N + NO = N2 + O
!
      qqk(kloop,36)=qk(kloop,36)*y(kloop,10)*y(kloop,12)
!
!....         NO + NO3 = 2 NO2
!
      qqk(kloop,37)=qk(kloop,37)*y(kloop,12)*y(kloop,14)
!
!....         NO2 + NO3 = N2O5
!
      qqk(kloop,38)=qk(kloop,38)*y(kloop,13)*y(kloop,14)
!
!....         N2O5 = NO2 + NO3
!
      qqk(kloop,39)=qk(kloop,39)*y(kloop,15)
!
!....         Cl + O3 = ClO + O2
!
      qqk(kloop,40)=qk(kloop,40)*y(kloop,19)*y(kloop,3)
!
!....         ClO + O = Cl + O2
!
      qqk(kloop,41)=qk(kloop,41)*y(kloop,21)*y(kloop,1)
!
!....         ClONO2 + O = ClO + NO3
!
      qqk(kloop,42)=qk(kloop,42)*y(kloop,26)*y(kloop,1)
!
!....         Cl + H2 = H + HCl
!
      qqk(kloop,43)=qk(kloop,43)*y(kloop,19)*y(kloop,5)
!
!....         Cl + HO2 = HCl + O2
!
      qqk(kloop,44)=qk(kloop,44)*y(kloop,19)*y(kloop,7)
!
!....         Cl + HO2 = ClO + OH
!
      qqk(kloop,45)=qk(kloop,45)*y(kloop,19)*y(kloop,7)
!
!....         Cl + H2O2 = HCl + HO2
!
      qqk(kloop,46)=qk(kloop,46)*y(kloop,19)*y(kloop,9)
!
!....         ClO + OH = Cl + HO2
!
      qqk(kloop,47)=qk(kloop,47)*y(kloop,21)*y(kloop,6)
!
!....         ClO + OH = HCl + O2
!
      qqk(kloop,48)=qk(kloop,48)*y(kloop,21)*y(kloop,6)
!
!....         ClO + HO2 = HOCl + O2
!
      qqk(kloop,49)=qk(kloop,49)*y(kloop,21)*y(kloop,7)
!
!....         ClO + HO2 = HCl + O3
!
      qqk(kloop,50)=qk(kloop,50)*y(kloop,21)*y(kloop,7)
!
!....         HCl + OH = Cl + H2O
!
      qqk(kloop,51)=qk(kloop,51)*y(kloop,24)*y(kloop,6)
!
!....         HOCl + OH = ClO + H2O
!
      qqk(kloop,52)=qk(kloop,52)*y(kloop,25)*y(kloop,6)
!
!....         ClONO2 + OH = HOCl + NO3
!
      qqk(kloop,53)=qk(kloop,53)*y(kloop,26)*y(kloop,6)
!
!....         ClO + NO = Cl + NO2
!
      qqk(kloop,54)=qk(kloop,54)*y(kloop,21)*y(kloop,12)
!
!....         ClO + NO2 = ClONO2
!
      qqk(kloop,55)=qk(kloop,55)*y(kloop,21)*y(kloop,13)
!
!....         CH4 + Cl = CH3O2 + HCl
!
      qqk(kloop,56)=qk(kloop,56)*y(kloop,42)*y(kloop,19)
!
!....         CH2O + Cl = CO + HCl + HO2
!
      qqk(kloop,57)=qk(kloop,57)*y(kloop,43)*y(kloop,19)
!
!....         Cl + HOCl = ClO + HCl
!
      qqk(kloop,58)=qk(kloop,58)*y(kloop,19)*y(kloop,25)
!
!....         ClO + ClO = 2 Cl + O2
!
      qqk(kloop,59)=qk(kloop,59)*y(kloop,21)*y(kloop,21)
!
!....         ClO + ClO = Cl2 + O2
!
      qqk(kloop,60)=qk(kloop,60)*y(kloop,21)*y(kloop,21)
!
!....         ClO + ClO = Cl + OClO
!
      qqk(kloop,61)=qk(kloop,61)*y(kloop,21)*y(kloop,21)
!
!....         ClO + ClO = Cl2O2
!
      qqk(kloop,62)=qk(kloop,62)*y(kloop,21)*y(kloop,21)
!
!....         Cl2O2 = 2 ClO
!
      qqk(kloop,63)=qk(kloop,63)*y(kloop,23)
!
!....         Cl + ClONO2 = Cl2 + NO3
!
      qqk(kloop,64)=qk(kloop,64)*y(kloop,19)*y(kloop,26)
!
!....         Br + O3 = BrO + O2
!
      qqk(kloop,65)=qk(kloop,65)*y(kloop,27)*y(kloop,3)
!
!....         BrO + O = Br + O2
!
      qqk(kloop,66)=qk(kloop,66)*y(kloop,28)*y(kloop,1)
!
!....         Br + HO2 = HBr + O2
!
      qqk(kloop,67)=qk(kloop,67)*y(kloop,27)*y(kloop,7)
!
!....         BrO + OH = Br + HO2
!
      qqk(kloop,68)=qk(kloop,68)*y(kloop,28)*y(kloop,6)
!
!....         BrO + HO2 = HOBr + O2
!
      qqk(kloop,69)=qk(kloop,69)*y(kloop,28)*y(kloop,7)
!
!....         HBr + OH = Br + H2O
!
      qqk(kloop,70)=qk(kloop,70)*y(kloop,29)*y(kloop,6)
!
!....         BrO + NO = Br + NO2
!
      qqk(kloop,71)=qk(kloop,71)*y(kloop,28)*y(kloop,12)
!
!....         BrO + NO2 = BrONO2
!
      qqk(kloop,72)=qk(kloop,72)*y(kloop,28)*y(kloop,13)
!
!....         Br + CH2O = CO + HBr + HO2
!
      qqk(kloop,73)=qk(kloop,73)*y(kloop,27)*y(kloop,43)
!
!....         BrO + ClO = Br + OClO
!
      qqk(kloop,74)=qk(kloop,74)*y(kloop,28)*y(kloop,21)
!
!....         BrO + ClO = Br + Cl + O2
!
      qqk(kloop,75)=qk(kloop,75)*y(kloop,28)*y(kloop,21)
!
!....         BrO + ClO = BrCl + O2
!
      qqk(kloop,76)=qk(kloop,76)*y(kloop,28)*y(kloop,21)
!
!....         BrO + BrO = 2 Br + O2
!
      qqk(kloop,77)=qk(kloop,77)*y(kloop,28)*y(kloop,28)
!
!....         CF2Cl2 + O1D = 2 Cl
!
      qqk(kloop,78)=qk(kloop,78)*y(kloop,37)*y(kloop,2)
!
!....         CH3Cl + OH = Cl + H2O + HO2
!
      qqk(kloop,79)=qk(kloop,79)*y(kloop,33)*y(kloop,6)
!
!....         CH3Br + OH = Br + H2O
!
      qqk(kloop,80)=qk(kloop,80)*y(kloop,38)*y(kloop,6)
!
!....         CH3CCl3 + OH = 3 Cl + H2O
!
      qqk(kloop,81)=qk(kloop,81)*y(kloop,34)*y(kloop,6)
!
!....         CH3Cl + Cl = CO + 2 HCl + HO2
!
      qqk(kloop,82)=qk(kloop,82)*y(kloop,33)*y(kloop,19)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,83)=qk(kloop,83)*y(kloop,15)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,84)=qk(kloop,84)*y(kloop,15)
!
!....         NO3 = HNO3
!
      qqk(kloop,85)=qk(kloop,85)*y(kloop,14)
!
!....         H2O + NO2 = HNO3 + HONO
!
      qqk(kloop,86)=qk(kloop,86)*y(kloop,8)*y(kloop,13)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,87)=qk(kloop,87)*y(kloop,26)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,88)=qk(kloop,88)*y(kloop,31)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,89)=qk(kloop,89)*y(kloop,26)*y(kloop,24)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,90)=qk(kloop,90)*y(kloop,24)*y(kloop,25)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,91)=qk(kloop,91)*y(kloop,24)*y(kloop,30)
!
!....         CH4 + O1D = CH3O2 + OH
!
      qqk(kloop,92)=qk(kloop,92)*y(kloop,42)*y(kloop,2)
!
!....         CH4 + O1D = CH2O + H + HO2
!
      qqk(kloop,93)=qk(kloop,93)*y(kloop,42)*y(kloop,2)
!
!....         CH4 + O1D = CH2O + H2
!
      qqk(kloop,94)=qk(kloop,94)*y(kloop,42)*y(kloop,2)
!
!....         CH2O + O = CO + HO2 + OH
!
      qqk(kloop,95)=qk(kloop,95)*y(kloop,43)*y(kloop,1)
!
!....         CH4 + OH = CH3O2 + H2O
!
      qqk(kloop,96)=qk(kloop,96)*y(kloop,42)*y(kloop,6)
!
!....         CO + OH = H
!
      qqk(kloop,97)=qk(kloop,97)*y(kloop,41)*y(kloop,6)
!
!....         CH2O + OH = CO + H2O + HO2
!
      qqk(kloop,98)=qk(kloop,98)*y(kloop,43)*y(kloop,6)
!
!....         CH3OH + OH = CH2O + HO2
!
      qqk(kloop,99)=qk(kloop,99)*y(kloop,47)*y(kloop,6)
!
!....         CH3OOH + OH = CH3O2 + H2O
!
      qqk(kloop,100)=qk(kloop,100)*y(kloop,48)*y(kloop,6)
!
!....         CH3OOH + OH = CH2O + H2O + OH
!
      qqk(kloop,101)=qk(kloop,101)*y(kloop,48)*y(kloop,6)
!
!....         CH3O2 + HO2 = CH3OOH + O2
!
      qqk(kloop,102)=qk(kloop,102)*y(kloop,45)*y(kloop,7)
!
!....         CH3O2 + NO = CH2O + HO2 + NO2
!
      qqk(kloop,103)=qk(kloop,103)*y(kloop,45)*y(kloop,12)
!
!....         CH3O2 + NO2 = CH3O2NO2
!
      qqk(kloop,104)=qk(kloop,104)*y(kloop,45)*y(kloop,13)
!
!....         CH3O2NO2 = CH3O2 + NO2
!
      qqk(kloop,105)=qk(kloop,105)*y(kloop,49)
!
!....         CH3O2 + CH3O2 =  1.33 CH2O +  0.66 CH3OH +  0.80 HO2
!
      qqk(kloop,106)=qk(kloop,106)*y(kloop,45)*y(kloop,45)
!
!....         CH2O + HO2 = CH3O3
!
      qqk(kloop,107)=qk(kloop,107)*y(kloop,43)*y(kloop,7)
!
!....         CH3O3 = CH2O + HO2
!
      qqk(kloop,108)=qk(kloop,108)*y(kloop,46)
!
!....         CH3O3 + NO = HCOOH + HO2 + NO2
!
      qqk(kloop,109)=qk(kloop,109)*y(kloop,46)*y(kloop,12)
!
!....         CH3O3 + HO2 = HCOOH
!
      qqk(kloop,110)=qk(kloop,110)*y(kloop,46)*y(kloop,7)
!
!....         CH3O3 + CH3O3 = 2 HCOOH + 2 HO2
!
      qqk(kloop,111)=qk(kloop,111)*y(kloop,46)*y(kloop,46)
!
!....         HCOOH + OH = CO + HO2
!
      qqk(kloop,112)=qk(kloop,112)*y(kloop,44)*y(kloop,6)
!
!....         C2H6 + OH = ETO2
!
      qqk(kloop,113)=qk(kloop,113)*y(kloop,50)*y(kloop,6)
!
!....         C2H6 + Cl = ETO2 + HCl
!
      qqk(kloop,114)=qk(kloop,114)*y(kloop,50)*y(kloop,19)
!
!....         ETO2 + NO = ALD2 + HO2 + NO2
!
      qqk(kloop,115)=qk(kloop,115)*y(kloop,57)*y(kloop,12)
!
!....         ETO2 + ETO2 =  1.60 ALD2 +  1.60 HO2
!
      qqk(kloop,116)=qk(kloop,116)*y(kloop,57)*y(kloop,57)
!
!....         ETO2 + HO2 = ETP
!
      qqk(kloop,117)=qk(kloop,117)*y(kloop,57)*y(kloop,7)
!
!....         ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH
!
      qqk(kloop,118)=qk(kloop,118)*y(kloop,58)*y(kloop,6)
!
!....         ALD2 + OH = MCO3
!
      qqk(kloop,119)=qk(kloop,119)*y(kloop,54)*y(kloop,6)
!
!....         ALD2 + NO3 = HNO3 + MCO3
!
      qqk(kloop,120)=qk(kloop,120)*y(kloop,54)*y(kloop,14)
!
!....         MCO3 + NO2 = PAN
!
      qqk(kloop,121)=qk(kloop,121)*y(kloop,52)*y(kloop,13)
!
!....         MCO3 + NO = CH3O2 + NO2
!
      qqk(kloop,122)=qk(kloop,122)*y(kloop,52)*y(kloop,12)
!
!....         HO2 + MCO3 =  0.75 MAP +  0.25 O3
!
      qqk(kloop,123)=qk(kloop,123)*y(kloop,7)*y(kloop,52)
!
!....         MCO3 + MCO3 = 2 CH3O2
!
      qqk(kloop,124)=qk(kloop,124)*y(kloop,52)*y(kloop,52)
!
!....         CH3O2 + MCO3 = CH2O +  0.90 CH3O2 +  0.90 HO2
!
      qqk(kloop,125)=qk(kloop,125)*y(kloop,45)*y(kloop,52)
!
!....         PAN = MCO3 + NO2
!
      qqk(kloop,126)=qk(kloop,126)*y(kloop,59)
!
!....         GLYX + OH = 2 CO + HO2
!
      qqk(kloop,127)=qk(kloop,127)*y(kloop,51)*y(kloop,6)
!
!....         C3H8 + OH = A3O2
!
      qqk(kloop,128)=qk(kloop,128)*y(kloop,61)*y(kloop,6)
!
!....         C3H8 + Cl = A3O2 + HCl
!
      qqk(kloop,129)=qk(kloop,129)*y(kloop,61)*y(kloop,19)
!
!....         A3O2 + NO =  0.82 ACET +  0.97 HO2 +  0.97 NO2 +  0.16 RCHO
!
      qqk(kloop,130)=qk(kloop,130)*y(kloop,68)*y(kloop,12)
!
!....         A3O2 + A3O2 =  0.40 ALD2 +  1.20 HO2 +  1.20 RCHO
!
      qqk(kloop,131)=qk(kloop,131)*y(kloop,68)*y(kloop,68)
!
!....         A3O2 + HO2 = RA3P
!
      qqk(kloop,132)=qk(kloop,132)*y(kloop,68)*y(kloop,7)
!
!....         OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,133)=qk(kloop,133)*y(kloop,6)*y(kloop,69)
!
!....         ACET + OH = ATO2
!
      qqk(kloop,134)=qk(kloop,134)*y(kloop,66)*y(kloop,6)
!
!....         ATO2 + NO =  0.96 HO2 +  0.96 MGLY +  0.96 NO2
!
      qqk(kloop,135)=qk(kloop,135)*y(kloop,65)*y(kloop,12)
!
!....         ATO2 + HO2 = CH3O2 + MCO3
!
      qqk(kloop,136)=qk(kloop,136)*y(kloop,65)*y(kloop,7)
!
!....         MGLY + OH = CO + MCO3
!
      qqk(kloop,137)=qk(kloop,137)*y(kloop,62)*y(kloop,6)
!
!....         OH + RCHO = ECO3
!
      qqk(kloop,138)=qk(kloop,138)*y(kloop,6)*y(kloop,81)
!
!....         NO3 + RCHO = ECO3 + HNO3
!
      qqk(kloop,139)=qk(kloop,139)*y(kloop,14)*y(kloop,81)
!
!....         ECO3 + NO2 = PPN
!
      qqk(kloop,140)=qk(kloop,140)*y(kloop,64)*y(kloop,13)
!
!....         ECO3 + NO = ETO2 + NO2
!
      qqk(kloop,141)=qk(kloop,141)*y(kloop,64)*y(kloop,12)
!
!....         ECO3 + ECO3 = 2 ETO2
!
      qqk(kloop,142)=qk(kloop,142)*y(kloop,64)*y(kloop,64)
!
!....         CH3O2 + ECO3 =  0.15 ALD2 + CH2O +  0.85 ETO2 +  0.65 HO2
!
      qqk(kloop,143)=qk(kloop,143)*y(kloop,45)*y(kloop,64)
!
!....         ECO3 + MCO3 = CH3O2 + ETO2
!
      qqk(kloop,144)=qk(kloop,144)*y(kloop,64)*y(kloop,52)
!
!....         ECO3 + HO2 =  0.25 ALD2 +  0.25 HO2 +  0.25 O3 +  0.75 RP
!
      qqk(kloop,145)=qk(kloop,145)*y(kloop,64)*y(kloop,7)
!
!....         PPN = ECO3 + NO2
!
      qqk(kloop,146)=qk(kloop,146)*y(kloop,70)
!
!....         OH + RP =  0.50 ALD2 +  0.50 ECO3 +  0.50 OH
!
      qqk(kloop,147)=qk(kloop,147)*y(kloop,6)*y(kloop,74)
!
!....         ISOP + OH = RIO2
!
      qqk(kloop,148)=qk(kloop,148)*y(kloop,83)*y(kloop,6)
!
!....         MVK + OH = VRO2
!
      qqk(kloop,149)=qk(kloop,149)*y(kloop,72)*y(kloop,6)
!
!....         MACR + OH = MAO3
!
      qqk(kloop,150)=qk(kloop,150)*y(kloop,73)*y(kloop,6)
!
!....         MACR + OH = MRO2
!
      qqk(kloop,151)=qk(kloop,151)*y(kloop,73)*y(kloop,6)
!
!....         HAC + OH = GLYX + HO2
!
      qqk(kloop,152)=qk(kloop,152)*y(kloop,55)*y(kloop,6)
!
!....         HAC + OH = HACO
!
      qqk(kloop,153)=qk(kloop,153)*y(kloop,55)*y(kloop,6)
!
!....         ISOP + O3 =  0.80 CH2O +  0.07 CH3OH +  0.05 CO +  0.06 HO2
!
      qqk(kloop,154)=qk(kloop,154)*y(kloop,83)*y(kloop,3)
!
!....         MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.16 CO +  0.06 HO2 +
!
      qqk(kloop,155)=qk(kloop,155)*y(kloop,72)*y(kloop,3)
!
!....         MACR + O3 =  0.70 CH2O +  0.09 CO +  0.21 HO2 +  0.80 MGLY +
!
      qqk(kloop,156)=qk(kloop,156)*y(kloop,73)*y(kloop,3)
!
!....         O3 + RIP =  0.70 CH2O
!
      qqk(kloop,157)=qk(kloop,157)*y(kloop,3)*y(kloop,85)
!
!....         O3 + VRP =  0.70 CH2O
!
      qqk(kloop,158)=qk(kloop,158)*y(kloop,3)*y(kloop,77)
!
!....         MRP + O3 =  0.70 CH2O
!
      qqk(kloop,159)=qk(kloop,159)*y(kloop,78)*y(kloop,3)
!
!....         ISOP + NO3 = INO2
!
      qqk(kloop,160)=qk(kloop,160)*y(kloop,83)*y(kloop,14)
!
!....         MACR + NO3 = HNO3 + MAO3
!
      qqk(kloop,161)=qk(kloop,161)*y(kloop,73)*y(kloop,14)
!
!....         MACR + NO3 = MAN2
!
      qqk(kloop,162)=qk(kloop,162)*y(kloop,73)*y(kloop,14)
!
!....         NO + RIO2 =  0.74 CH2O +  0.78 HO2 +  0.14 ISN1 +  0.32 MACR
!
      qqk(kloop,163)=qk(kloop,163)*y(kloop,12)*y(kloop,84)
!
!....         NO + VRO2 =  0.27 CH2O +  0.68 HAC +  0.27 HO2 +  0.05 ISN1
!
      qqk(kloop,164)=qk(kloop,164)*y(kloop,12)*y(kloop,75)
!
!....         MAO3 + NO = MAO2
!
      qqk(kloop,165)=qk(kloop,165)*y(kloop,71)*y(kloop,12)
!
!....         MAO2 + NO = CH2O + MCO3 + NO2
!
      qqk(kloop,166)=qk(kloop,166)*y(kloop,63)*y(kloop,12)
!
!....         MRO2 + NO = HACN + HO2 + NO2
!
      qqk(kloop,167)=qk(kloop,167)*y(kloop,76)*y(kloop,12)
!
!....         HACO + NO = CH3O3 + NO2
!
      qqk(kloop,168)=qk(kloop,168)*y(kloop,53)*y(kloop,12)
!
!....         INO2 + NO =  0.15 CH2O +  0.80 HO2 +  0.75 ISN1 +  0.10 MACR
!
      qqk(kloop,169)=qk(kloop,169)*y(kloop,87)*y(kloop,12)
!
!....         ISNR + NO =  0.95 ACET +  1.90 HAC +  0.05 HO2 +  0.05 ISN1
!
      qqk(kloop,170)=qk(kloop,170)*y(kloop,88)*y(kloop,12)
!
!....         MAN2 + NO = HO2 + MGLY + 2 NO2
!
      qqk(kloop,171)=qk(kloop,171)*y(kloop,80)*y(kloop,12)
!
!....         HO2 + RIO2 = RIP
!
      qqk(kloop,172)=qk(kloop,172)*y(kloop,7)*y(kloop,84)
!
!....         HO2 + VRO2 = VRP
!
      qqk(kloop,173)=qk(kloop,173)*y(kloop,7)*y(kloop,75)
!
!....         HO2 + MAO3 = RP
!
      qqk(kloop,174)=qk(kloop,174)*y(kloop,7)*y(kloop,71)
!
!....         HO2 + MRO2 = MRP
!
      qqk(kloop,175)=qk(kloop,175)*y(kloop,7)*y(kloop,76)
!
!....         HACO + HO2 = CH3OOH
!
      qqk(kloop,176)=qk(kloop,176)*y(kloop,53)*y(kloop,7)
!
!....         HO2 + ISNR = PRN2
!
      qqk(kloop,177)=qk(kloop,177)*y(kloop,7)*y(kloop,88)
!
!....         HO2 + INO2 = PRN2
!
      qqk(kloop,178)=qk(kloop,178)*y(kloop,7)*y(kloop,87)
!
!....         HO2 + MAN2 = PRN2
!
      qqk(kloop,179)=qk(kloop,179)*y(kloop,7)*y(kloop,80)
!
!....         INO2 + NO2 = 2 PRN2
!
      qqk(kloop,180)=qk(kloop,180)*y(kloop,87)*y(kloop,13)
!
!....         HACN + OH = HO2 + MGLY
!
      qqk(kloop,181)=qk(kloop,181)*y(kloop,67)*y(kloop,6)
!
!....         ISN1 + OH = ISNR
!
      qqk(kloop,182)=qk(kloop,182)*y(kloop,86)*y(kloop,6)
!
!....         OH + RIP =  0.37 CH2O + HO2 +  0.16 MACR +  0.21 MVK +  0.50
!
      qqk(kloop,183)=qk(kloop,183)*y(kloop,6)*y(kloop,85)
!
!....         OH + VRP =  0.50 HAC + HO2 + VRO2
!
      qqk(kloop,184)=qk(kloop,184)*y(kloop,6)*y(kloop,77)
!
!....         MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH
!
      qqk(kloop,185)=qk(kloop,185)*y(kloop,56)*y(kloop,6)
!
!....         MRP + OH =  0.50 HAC + HO2 +  0.50 MRO2
!
      qqk(kloop,186)=qk(kloop,186)*y(kloop,78)*y(kloop,6)
!
!....         OH + PRN2 =  0.50 ISNR +  0.50 NO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,187)=qk(kloop,187)*y(kloop,6)*y(kloop,82)
!
!....         HO2 + MAO3 = CH2O + ETO2 + O3
!
      qqk(kloop,188)=qk(kloop,188)*y(kloop,7)*y(kloop,71)
!
!....         HO2 + MAO2 = ACET + O3
!
      qqk(kloop,189)=qk(kloop,189)*y(kloop,7)*y(kloop,63)
!
!....         HACO + HO2 = ETO2 + O3
!
      qqk(kloop,190)=qk(kloop,190)*y(kloop,53)*y(kloop,7)
!
!....         MAO3 + NO2 = MPAN
!
      qqk(kloop,191)=qk(kloop,191)*y(kloop,71)*y(kloop,13)
!
!....         MPAN = MAO3 + NO2
!
      qqk(kloop,192)=qk(kloop,192)*y(kloop,79)
!
!....         HACO + NO2 = IPAN
!
      qqk(kloop,193)=qk(kloop,193)*y(kloop,53)*y(kloop,13)
!
!....         IPAN = HACO + NO2
!
      qqk(kloop,194)=qk(kloop,194)*y(kloop,60)
!
!....         CH3O2 + MAO3 = CH2O +  0.50 HO2 +  0.85 MAO2
!
      qqk(kloop,195)=qk(kloop,195)*y(kloop,45)*y(kloop,71)
!
!....         MAO3 + MCO3 = CH3O2 + MAO2
!
      qqk(kloop,196)=qk(kloop,196)*y(kloop,71)*y(kloop,52)
!
!....         MAO3 + MAO3 = 2 MAO2
!
      qqk(kloop,197)=qk(kloop,197)*y(kloop,71)*y(kloop,71)
!
!....         CH3O2 + HACO = CH2O +  0.85 CH3O2 +  0.50 HO2
!
      qqk(kloop,198)=qk(kloop,198)*y(kloop,45)*y(kloop,53)
!
!....         HACO + MCO3 = 2 CH3O2
!
      qqk(kloop,199)=qk(kloop,199)*y(kloop,53)*y(kloop,52)
!
!....         HACO + HACO = 2 CH3O2
!
      qqk(kloop,200)=qk(kloop,200)*y(kloop,53)*y(kloop,53)
!
!                  start photolytic reactions
!
!....  O2 + hv = 2 O
!
      qqj(kloop,1)=qj(kloop,1)*y(kloop,90)
!
!....  O3 + hv = O + O2
!
      qqj(kloop,2)=qj(kloop,2)*y(kloop,3)
!
!....  O3 + hv = O1D + O2
!
      qqj(kloop,3)=qj(kloop,3)*y(kloop,3)
!
!....  N2O + hv = N2 + O1D
!
      qqj(kloop,4)=qj(kloop,4)*y(kloop,11)
!
!....  NO + hv = N + O
!
      qqj(kloop,5)=qj(kloop,5)*y(kloop,12)
!
!....  NO2 + hv = NO + O
!
      qqj(kloop,6)=qj(kloop,6)*y(kloop,13)
!
!....  NO3 + hv = NO + O2
!
      qqj(kloop,7)=qj(kloop,7)*y(kloop,14)
!
!....  NO3 + hv = NO2 + O
!
      qqj(kloop,8)=qj(kloop,8)*y(kloop,14)
!
!....  N2O5 + hv = NO2 + NO3
!
      qqj(kloop,9)=qj(kloop,9)*y(kloop,15)
!
!....  HONO + hv = NO + OH
!
      qqj(kloop,10)=qj(kloop,10)*y(kloop,16)
!
!....  HNO3 + hv = NO2 + OH
!
      qqj(kloop,11)=qj(kloop,11)*y(kloop,17)
!
!....  HO2NO2 + hv = HO2 + NO2
!
      qqj(kloop,12)=qj(kloop,12)*y(kloop,18)
!
!....  HO2NO2 + hv = NO3 + OH
!
      qqj(kloop,13)=qj(kloop,13)*y(kloop,18)
!
!....  H2O + hv = H + OH
!
      qqj(kloop,14)=qj(kloop,14)*y(kloop,8)
!
!....  HO2 + hv = O + OH
!
      qqj(kloop,15)=qj(kloop,15)*y(kloop,7)
!
!....  H2O2 + hv = 2 OH
!
      qqj(kloop,16)=qj(kloop,16)*y(kloop,9)
!
!....  ClO + hv = Cl + O
!
      qqj(kloop,17)=qj(kloop,17)*y(kloop,21)
!
!....  ClO + hv = Cl + O1D
!
      qqj(kloop,18)=qj(kloop,18)*y(kloop,21)
!
!....  HCl + hv = Cl + H
!
      qqj(kloop,19)=qj(kloop,19)*y(kloop,24)
!
!....  HOCl + hv = Cl + OH
!
      qqj(kloop,20)=qj(kloop,20)*y(kloop,25)
!
!....  ClONO2 + hv = Cl + NO3
!
      qqj(kloop,21)=qj(kloop,21)*y(kloop,26)
!
!....  ClONO2 + hv = ClO + NO2
!
      qqj(kloop,22)=qj(kloop,22)*y(kloop,26)
!
!....  Cl2 + hv = 2 Cl
!
      qqj(kloop,23)=qj(kloop,23)*y(kloop,20)
!
!....  OClO + hv = ClO + O
!
      qqj(kloop,24)=qj(kloop,24)*y(kloop,22)
!
!....  Cl2O2 + hv = 2 Cl + O2
!
      qqj(kloop,25)=qj(kloop,25)*y(kloop,23)
!
!....  BrO + hv = Br + O
!
      qqj(kloop,26)=qj(kloop,26)*y(kloop,28)
!
!....  HBr + hv = Br + H
!
      qqj(kloop,27)=qj(kloop,27)*y(kloop,29)
!
!....  HOBr + hv = Br + OH
!
      qqj(kloop,28)=qj(kloop,28)*y(kloop,30)
!
!....  BrONO2 + hv = Br + NO3
!
      qqj(kloop,29)=qj(kloop,29)*y(kloop,31)
!
!....  BrONO2 + hv = BrO + NO2
!
      qqj(kloop,30)=qj(kloop,30)*y(kloop,31)
!
!....  BrCl + hv = Br + Cl
!
      qqj(kloop,31)=qj(kloop,31)*y(kloop,32)
!
!....  CH3Cl + hv = CH3O2 + Cl
!
      qqj(kloop,32)=qj(kloop,32)*y(kloop,33)
!
!....  CH3Br + hv = Br + CH3O2
!
      qqj(kloop,33)=qj(kloop,33)*y(kloop,38)
!
!....  CFCl3 + hv = 3 Cl
!
      qqj(kloop,34)=qj(kloop,34)*y(kloop,36)
!
!....  CF2Cl2 + hv = 2 Cl
!
      qqj(kloop,35)=qj(kloop,35)*y(kloop,37)
!
!....  CCl4 + hv = 4 Cl
!
      qqj(kloop,36)=qj(kloop,36)*y(kloop,35)
!
!....  CH3CCl3 + hv = 3 Cl
!
      qqj(kloop,37)=qj(kloop,37)*y(kloop,34)
!
!....  CF3Br + hv = Br
!
      qqj(kloop,38)=qj(kloop,38)*y(kloop,39)
!
!....  CF2ClBr + hv = Br + Cl
!
      qqj(kloop,39)=qj(kloop,39)*y(kloop,40)
!
!....  CH3OOH + hv = CH2O + HO2 + OH
!
      qqj(kloop,40)=qj(kloop,40)*y(kloop,48)
!
!....  CH2O + hv = CO + H + HO2
!
      qqj(kloop,41)=qj(kloop,41)*y(kloop,43)
!
!....  CH2O + hv = CO + H2
!
      qqj(kloop,42)=qj(kloop,42)*y(kloop,43)
!
!....  CH3O2NO2 + hv = CH2O + HO2 + NO3
!
      qqj(kloop,43)=qj(kloop,43)*y(kloop,49)
!
!....  ALD2 + hv = CH3O2 + CO + HO2
!
      qqj(kloop,44)=qj(kloop,44)*y(kloop,54)
!
!....  ALD2 + hv = CH4 + CO
!
      qqj(kloop,45)=qj(kloop,45)*y(kloop,54)
!
!....  GLYX + hv = CH2O + CO
!
      qqj(kloop,46)=qj(kloop,46)*y(kloop,51)
!
!....  RCHO + hv = CO + ETO2 + HO2
!
      qqj(kloop,47)=qj(kloop,47)*y(kloop,81)
!
!....  ACET + hv = CH3O2 + MCO3
!
      qqj(kloop,48)=qj(kloop,48)*y(kloop,66)
!
!....  MGLY + hv = CO + HO2 + MCO3
!
      qqj(kloop,49)=qj(kloop,49)*y(kloop,62)
!
!....  HAC + hv = CH2O + CO + 2 HO2
!
      qqj(kloop,50)=qj(kloop,50)*y(kloop,55)
!
!....  HACN + hv = CH2O + HO2 + MCO3
!
      qqj(kloop,51)=qj(kloop,51)*y(kloop,67)
!
!....  PAN + hv = MCO3 + NO2
!
      qqj(kloop,52)=qj(kloop,52)*y(kloop,59)
!
!....  PAN + hv = CH3O2 + NO3
!
      qqj(kloop,53)=qj(kloop,53)*y(kloop,59)
!
!....  ETP + hv = ALD2 + HO2 + OH
!
      qqj(kloop,54)=qj(kloop,54)*y(kloop,58)
!
!....  RA3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,55)=qj(kloop,55)*y(kloop,69)
!
!....  RP + hv = ALD2 + HO2 + OH
!
      qqj(kloop,56)=qj(kloop,56)*y(kloop,74)
!
!....  RIP + hv = HAC + HO2 + OH
!
      qqj(kloop,57)=qj(kloop,57)*y(kloop,85)
!
!....  VRP + hv = HO2 + OH + RCHO
!
      qqj(kloop,58)=qj(kloop,58)*y(kloop,77)
!
!....  MAP + hv = CH2O + HO2 + OH
!
      qqj(kloop,59)=qj(kloop,59)*y(kloop,56)
!
!....  MRP + hv = HO2 + OH + RCHO
!
      qqj(kloop,60)=qj(kloop,60)*y(kloop,78)
!
!....  PRN2 + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,61)=qj(kloop,61)*y(kloop,82)
      enddo
!
!.... End of kloop
!
      return
      end
