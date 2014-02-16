!=======================================================================
!
! $Id: setkin_ratecalc.F90,v 1.2 2011-08-08 19:37:04 mrdamon Exp $
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
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
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
!....         H2O + O1D = 2 OH
!
      qqk(kloop,6)=qk(kloop,6)*y(kloop,12)*y(kloop,2)
!
!....         H2 + O1D = H + OH
!
      qqk(kloop,7)=qk(kloop,7)*y(kloop,17)*y(kloop,2)
!
!....         N2O + O1D = N2 + O2
!
      qqk(kloop,8)=qk(kloop,8)*y(kloop,4)*y(kloop,2)
!
!....         N2O + O1D = 2 NO
!
      qqk(kloop,9)=qk(kloop,9)*y(kloop,4)*y(kloop,2)
!
!....         CH4 + O1D = CH3O2 + OH
!
      qqk(kloop,10)=qk(kloop,10)*y(kloop,18)*y(kloop,2)
!
!....         CH4 + O1D = CH2O + H + HO2
!
      qqk(kloop,11)=qk(kloop,11)*y(kloop,18)*y(kloop,2)
!
!....         CH4 + O1D = CH2O + H2
!
      qqk(kloop,12)=qk(kloop,12)*y(kloop,18)*y(kloop,2)
!
!....         CF2Cl2 + O1D = 2 Cl
!
      qqk(kloop,13)=qk(kloop,13)*y(kloop,40)*y(kloop,2)
!
!....         CFC113 + O1D = 3 Cl
!
      qqk(kloop,14)=qk(kloop,14)*y(kloop,41)*y(kloop,2)
!
!....         CFC114 + O1D = 2 Cl
!
      qqk(kloop,15)=qk(kloop,15)*y(kloop,42)*y(kloop,2)
!
!....         CFC115 + O1D = Cl
!
      qqk(kloop,16)=qk(kloop,16)*y(kloop,43)*y(kloop,2)
!
!....         HCFC22 + O1D = Cl
!
      qqk(kloop,17)=qk(kloop,17)*y(kloop,44)*y(kloop,2)
!
!....         HCFC141b + O1D = 2 Cl
!
      qqk(kloop,18)=qk(kloop,18)*y(kloop,47)*y(kloop,2)
!
!....         HCFC142b + O1D = Cl
!
      qqk(kloop,19)=qk(kloop,19)*y(kloop,48)*y(kloop,2)
!
!....         H + O2 = HO2
!
      qqk(kloop,20)=qk(kloop,20)*y(kloop,13)
!
!....         H + O3 = O2 + OH
!
      qqk(kloop,21)=qk(kloop,21)*y(kloop,13)*y(kloop,3)
!
!....         H2 + OH = H + H2O
!
      qqk(kloop,22)=qk(kloop,22)*y(kloop,17)*y(kloop,14)
!
!....         O3 + OH = HO2 + O2
!
      qqk(kloop,23)=qk(kloop,23)*y(kloop,3)*y(kloop,14)
!
!....         O + OH = H + O2
!
      qqk(kloop,24)=qk(kloop,24)*y(kloop,1)*y(kloop,14)
!
!....         OH + OH = H2O + O
!
      qqk(kloop,25)=qk(kloop,25)*y(kloop,14)*y(kloop,14)
!
!....         HO2 + O = O2 + OH
!
      qqk(kloop,26)=qk(kloop,26)*y(kloop,15)*y(kloop,1)
!
!....         HO2 + O3 = 2 O2 + OH
!
      qqk(kloop,27)=qk(kloop,27)*y(kloop,15)*y(kloop,3)
!
!....         H + HO2 = 2 OH
!
      qqk(kloop,28)=qk(kloop,28)*y(kloop,13)*y(kloop,15)
!
!....         HO2 + OH = H2O + O2
!
      qqk(kloop,29)=qk(kloop,29)*y(kloop,15)*y(kloop,14)
!
!....         HO2 + HO2 = H2O2 + O2
!
      qqk(kloop,30)=qk(kloop,30)*y(kloop,15)*y(kloop,15)
!
!....         H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      qqk(kloop,31)=qk(kloop,31)*y(kloop,12)*y(kloop,15)*y(kloop,15)
!
!....         H2O2 + OH = H2O + HO2
!
      qqk(kloop,32)=qk(kloop,32)*y(kloop,16)*y(kloop,14)
!
!....         N + O2 = NO + O
!
      qqk(kloop,33)=qk(kloop,33)*y(kloop,5)
!
!....         N + NO = N2 + O
!
      qqk(kloop,34)=qk(kloop,34)*y(kloop,5)*y(kloop,6)
!
!....         NO + O3 = NO2 + O2
!
      qqk(kloop,35)=qk(kloop,35)*y(kloop,6)*y(kloop,3)
!
!....         NO2 + OH = HNO3
!
      qqk(kloop,36)=qk(kloop,36)*y(kloop,7)*y(kloop,14)
!
!....         HO2 + NO = NO2 + OH
!
      qqk(kloop,37)=qk(kloop,37)*y(kloop,15)*y(kloop,6)
!
!....         NO2 + O = NO + O2
!
      qqk(kloop,38)=qk(kloop,38)*y(kloop,7)*y(kloop,1)
!
!....         NO2 + O3 = NO3 + O2
!
      qqk(kloop,39)=qk(kloop,39)*y(kloop,7)*y(kloop,3)
!
!....         HO2 + NO2 = HO2NO2
!
      qqk(kloop,40)=qk(kloop,40)*y(kloop,15)*y(kloop,7)
!
!....         NO3 + O = NO2 + O2
!
      qqk(kloop,41)=qk(kloop,41)*y(kloop,8)*y(kloop,1)
!
!....         NO + NO3 = 2 NO2
!
      qqk(kloop,42)=qk(kloop,42)*y(kloop,6)*y(kloop,8)
!
!....         NO2 + NO3 = N2O5
!
      qqk(kloop,43)=qk(kloop,43)*y(kloop,7)*y(kloop,8)
!
!....         N2O5 = NO2 + NO3
!
      qqk(kloop,44)=qk(kloop,44)*y(kloop,9)
!
!....         HNO3 + OH = H2O + NO3
!
      qqk(kloop,45)=qk(kloop,45)*y(kloop,10)*y(kloop,14)
!
!....         HO2NO2 = HO2 + NO2
!
      qqk(kloop,46)=qk(kloop,46)*y(kloop,11)
!
!....         HO2NO2 + OH = H2O + NO2 + O2
!
      qqk(kloop,47)=qk(kloop,47)*y(kloop,11)*y(kloop,14)
!
!....         Cl + O3 = ClO + O2
!
      qqk(kloop,48)=qk(kloop,48)*y(kloop,23)*y(kloop,3)
!
!....         Cl + H2 = H + HCl
!
      qqk(kloop,49)=qk(kloop,49)*y(kloop,23)*y(kloop,17)
!
!....         Cl + H2O2 = HCl + HO2
!
      qqk(kloop,50)=qk(kloop,50)*y(kloop,23)*y(kloop,16)
!
!....         Cl + HO2 = HCl + O2
!
      qqk(kloop,51)=qk(kloop,51)*y(kloop,23)*y(kloop,15)
!
!....         Cl + HO2 = ClO + OH
!
      qqk(kloop,52)=qk(kloop,52)*y(kloop,23)*y(kloop,15)
!
!....         ClO + O = Cl + O2
!
      qqk(kloop,53)=qk(kloop,53)*y(kloop,25)*y(kloop,1)
!
!....         ClO + OH = Cl + HO2
!
      qqk(kloop,54)=qk(kloop,54)*y(kloop,25)*y(kloop,14)
!
!....         ClO + OH = HCl + O2
!
      qqk(kloop,55)=qk(kloop,55)*y(kloop,25)*y(kloop,14)
!
!....         ClO + HO2 = HOCl + O2
!
      qqk(kloop,56)=qk(kloop,56)*y(kloop,25)*y(kloop,15)
!
!....         ClO + HO2 = HCl + O3
!
      qqk(kloop,57)=qk(kloop,57)*y(kloop,25)*y(kloop,15)
!
!....         ClO + NO = Cl + NO2
!
      qqk(kloop,58)=qk(kloop,58)*y(kloop,25)*y(kloop,6)
!
!....         ClO + NO2 = ClONO2
!
      qqk(kloop,59)=qk(kloop,59)*y(kloop,25)*y(kloop,7)
!
!....         ClO + ClO = 2 Cl + O2
!
      qqk(kloop,60)=qk(kloop,60)*y(kloop,25)*y(kloop,25)
!
!....         ClO + ClO = Cl2 + O2
!
      qqk(kloop,61)=qk(kloop,61)*y(kloop,25)*y(kloop,25)
!
!....         ClO + ClO = Cl + OClO
!
      qqk(kloop,62)=qk(kloop,62)*y(kloop,25)*y(kloop,25)
!
!....         ClO + ClO = Cl2O2
!
      qqk(kloop,63)=qk(kloop,63)*y(kloop,25)*y(kloop,25)
!
!....         Cl2O2 = 2 ClO
!
      qqk(kloop,64)=qk(kloop,64)*y(kloop,27)
!
!....         HCl + OH = Cl + H2O
!
      qqk(kloop,65)=qk(kloop,65)*y(kloop,28)*y(kloop,14)
!
!....         HOCl + OH = ClO + H2O
!
      qqk(kloop,66)=qk(kloop,66)*y(kloop,29)*y(kloop,14)
!
!....         ClONO2 + O = ClO + NO3
!
      qqk(kloop,67)=qk(kloop,67)*y(kloop,30)*y(kloop,1)
!
!....         ClONO2 + OH = HOCl + NO3
!
      qqk(kloop,68)=qk(kloop,68)*y(kloop,30)*y(kloop,14)
!
!....         Cl + ClONO2 = Cl2 + NO3
!
      qqk(kloop,69)=qk(kloop,69)*y(kloop,23)*y(kloop,30)
!
!....         Br + O3 = BrO + O2
!
      qqk(kloop,70)=qk(kloop,70)*y(kloop,31)*y(kloop,3)
!
!....         Br + HO2 = HBr + O2
!
      qqk(kloop,71)=qk(kloop,71)*y(kloop,31)*y(kloop,15)
!
!....         Br + CH2O = CO + HBr + HO2
!
      qqk(kloop,72)=qk(kloop,72)*y(kloop,31)*y(kloop,21)
!
!....         BrO + O = Br + O2
!
      qqk(kloop,73)=qk(kloop,73)*y(kloop,33)*y(kloop,1)
!
!....         BrO + HO2 = HOBr + O2
!
      qqk(kloop,74)=qk(kloop,74)*y(kloop,33)*y(kloop,15)
!
!....         BrO + NO = Br + NO2
!
      qqk(kloop,75)=qk(kloop,75)*y(kloop,33)*y(kloop,6)
!
!....         BrO + NO2 = BrONO2
!
      qqk(kloop,76)=qk(kloop,76)*y(kloop,33)*y(kloop,7)
!
!....         BrO + ClO = Br + OClO
!
      qqk(kloop,77)=qk(kloop,77)*y(kloop,33)*y(kloop,25)
!
!....         BrO + ClO = Br + Cl + O2
!
      qqk(kloop,78)=qk(kloop,78)*y(kloop,33)*y(kloop,25)
!
!....         BrO + ClO = BrCl + O2
!
      qqk(kloop,79)=qk(kloop,79)*y(kloop,33)*y(kloop,25)
!
!....         BrO + BrO = 2 Br + O2
!
      qqk(kloop,80)=qk(kloop,80)*y(kloop,33)*y(kloop,33)
!
!....         HBr + OH = Br + H2O
!
      qqk(kloop,81)=qk(kloop,81)*y(kloop,34)*y(kloop,14)
!
!....         CO + OH = H
!
      qqk(kloop,82)=qk(kloop,82)*y(kloop,22)*y(kloop,14)
!
!....         CH4 + OH = CH3O2 + H2O
!
      qqk(kloop,83)=qk(kloop,83)*y(kloop,18)*y(kloop,14)
!
!....         CH2O + OH = CO + H2O + HO2
!
      qqk(kloop,84)=qk(kloop,84)*y(kloop,21)*y(kloop,14)
!
!....         CH2O + O = CO + HO2 + OH
!
      qqk(kloop,85)=qk(kloop,85)*y(kloop,21)*y(kloop,1)
!
!....         CH4 + Cl = CH3O2 + HCl
!
      qqk(kloop,86)=qk(kloop,86)*y(kloop,18)*y(kloop,23)
!
!....         CH2O + Cl = CO + HCl + HO2
!
      qqk(kloop,87)=qk(kloop,87)*y(kloop,21)*y(kloop,23)
!
!....         CH3O2 + NO = CH2O + HO2 + NO2
!
      qqk(kloop,88)=qk(kloop,88)*y(kloop,19)*y(kloop,6)
!
!....         CH3O2 + HO2 = CH3OOH + O2
!
      qqk(kloop,89)=qk(kloop,89)*y(kloop,19)*y(kloop,15)
!
!....         CH3OOH + OH = CH3O2 + H2O
!
      qqk(kloop,90)=qk(kloop,90)*y(kloop,20)*y(kloop,14)
!
!....         CH3Cl + OH = Cl + H2O + HO2
!
      qqk(kloop,91)=qk(kloop,91)*y(kloop,37)*y(kloop,14)
!
!....         CH3CCl3 + OH = 3 Cl + H2O
!
      qqk(kloop,92)=qk(kloop,92)*y(kloop,46)*y(kloop,14)
!
!....         HCFC22 + OH = Cl + H2O
!
      qqk(kloop,93)=qk(kloop,93)*y(kloop,44)*y(kloop,14)
!
!....         HCFC141b + OH = 2 Cl + H2O
!
      qqk(kloop,94)=qk(kloop,94)*y(kloop,47)*y(kloop,14)
!
!....         HCFC142b + OH = Cl + H2O
!
      qqk(kloop,95)=qk(kloop,95)*y(kloop,48)*y(kloop,14)
!
!....         CH3Cl + Cl = CO + 2 HCl + HO2
!
      qqk(kloop,96)=qk(kloop,96)*y(kloop,37)*y(kloop,23)
!
!....         CH3Br + OH = Br + H2O + HO2
!
      qqk(kloop,97)=qk(kloop,97)*y(kloop,38)*y(kloop,14)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,98)=qk(kloop,98)*y(kloop,9)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,99)=qk(kloop,99)*y(kloop,30)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,100)=qk(kloop,100)*y(kloop,36)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,101)=qk(kloop,101)*y(kloop,30)*y(kloop,28)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,102)=qk(kloop,102)*y(kloop,28)*y(kloop,29)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,103)=qk(kloop,103)*y(kloop,28)*y(kloop,35)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,104)=qk(kloop,104)*y(kloop,9)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,105)=qk(kloop,105)*y(kloop,30)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,106)=qk(kloop,106)*y(kloop,36)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,107)=qk(kloop,107)*y(kloop,30)*y(kloop,28)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,108)=qk(kloop,108)*y(kloop,28)*y(kloop,29)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,109)=qk(kloop,109)*y(kloop,28)*y(kloop,35)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,110)=qk(kloop,110)*y(kloop,30)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,111)=qk(kloop,111)*y(kloop,36)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,112)=qk(kloop,112)*y(kloop,30)*y(kloop,28)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,113)=qk(kloop,113)*y(kloop,28)*y(kloop,29)
!
!....         BrONO2 + HCl = BrCl + HNO3
!
      qqk(kloop,114)=qk(kloop,114)*y(kloop,36)*y(kloop,28)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,115)=qk(kloop,115)*y(kloop,28)*y(kloop,35)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,116)=qk(kloop,116)*y(kloop,30)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,117)=qk(kloop,117)*y(kloop,36)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,118)=qk(kloop,118)*y(kloop,30)*y(kloop,28)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,119)=qk(kloop,119)*y(kloop,28)*y(kloop,29)
!
!....         BrONO2 + HCl = BrCl + HNO3
!
      qqk(kloop,120)=qk(kloop,120)*y(kloop,36)*y(kloop,28)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,121)=qk(kloop,121)*y(kloop,28)*y(kloop,35)
!
!....         HNO3 = NO2 + OH
!
      qqk(kloop,122)=qk(kloop,122)*y(kloop,10)
!
!                  start photolytic reactions
!
!....  O2 + hv = 2 O
!
      qqj(kloop,1)=qj(kloop,1)*y(kloop,54)
!
!....  O3 + hv = O + O2
!
      qqj(kloop,2)=qj(kloop,2)*y(kloop,3)
!
!....  O3 + hv = O1D + O2
!
      qqj(kloop,3)=qj(kloop,3)*y(kloop,3)
!
!....  HO2 + hv = O + OH
!
      qqj(kloop,4)=qj(kloop,4)*y(kloop,15)
!
!....  H2O2 + hv = 2 OH
!
      qqj(kloop,5)=qj(kloop,5)*y(kloop,16)
!
!....  H2O + hv = H + OH
!
      qqj(kloop,6)=qj(kloop,6)*y(kloop,12)
!
!....  NO + hv = N + O
!
      qqj(kloop,7)=qj(kloop,7)*y(kloop,6)
!
!....  NO2 + hv = NO + O
!
      qqj(kloop,8)=qj(kloop,8)*y(kloop,7)
!
!....  N2O + hv = N2 + O1D
!
      qqj(kloop,9)=qj(kloop,9)*y(kloop,4)
!
!....  NO3 + hv = NO2 + O
!
      qqj(kloop,10)=qj(kloop,10)*y(kloop,8)
!
!....  NO3 + hv = NO + O2
!
      qqj(kloop,11)=qj(kloop,11)*y(kloop,8)
!
!....  N2O5 + hv = NO2 + NO3
!
      qqj(kloop,12)=qj(kloop,12)*y(kloop,9)
!
!....  HNO3 + hv = NO2 + OH
!
      qqj(kloop,13)=qj(kloop,13)*y(kloop,10)
!
!....  HO2NO2 + hv = NO3 + OH
!
      qqj(kloop,14)=qj(kloop,14)*y(kloop,11)
!
!....  HO2NO2 + hv = HO2 + NO2
!
      qqj(kloop,15)=qj(kloop,15)*y(kloop,11)
!
!....  Cl2 + hv = 2 Cl
!
      qqj(kloop,16)=qj(kloop,16)*y(kloop,24)
!
!....  OClO + hv = ClO + O
!
      qqj(kloop,17)=qj(kloop,17)*y(kloop,26)
!
!....  Cl2O2 + hv = 2 Cl + O2
!
      qqj(kloop,18)=qj(kloop,18)*y(kloop,27)
!
!....  HOCl + hv = Cl + OH
!
      qqj(kloop,19)=qj(kloop,19)*y(kloop,29)
!
!....  ClONO2 + hv = Cl + NO3
!
      qqj(kloop,20)=qj(kloop,20)*y(kloop,30)
!
!....  ClONO2 + hv = ClO + NO2
!
      qqj(kloop,21)=qj(kloop,21)*y(kloop,30)
!
!....  BrCl + hv = Br + Cl
!
      qqj(kloop,22)=qj(kloop,22)*y(kloop,32)
!
!....  BrO + hv = Br + O
!
      qqj(kloop,23)=qj(kloop,23)*y(kloop,33)
!
!....  HOBr + hv = Br + OH
!
      qqj(kloop,24)=qj(kloop,24)*y(kloop,35)
!
!....  BrONO2 + hv = Br + NO3
!
      qqj(kloop,25)=qj(kloop,25)*y(kloop,36)
!
!....  BrONO2 + hv = BrO + NO2
!
      qqj(kloop,26)=qj(kloop,26)*y(kloop,36)
!
!....  CH3OOH + hv = CH2O + HO2 + OH
!
      qqj(kloop,27)=qj(kloop,27)*y(kloop,20)
!
!....  CH2O + hv = CO + H2
!
      qqj(kloop,28)=qj(kloop,28)*y(kloop,21)
!
!....  CH2O + hv = CO + H + HO2
!
      qqj(kloop,29)=qj(kloop,29)*y(kloop,21)
!
!....  CH3Cl + hv = CH3O2 + Cl
!
      qqj(kloop,30)=qj(kloop,30)*y(kloop,37)
!
!....  CCl4 + hv = 4 Cl
!
      qqj(kloop,31)=qj(kloop,31)*y(kloop,45)
!
!....  CH3CCl3 + hv = 3 Cl
!
      qqj(kloop,32)=qj(kloop,32)*y(kloop,46)
!
!....  CFCl3 + hv = 3 Cl
!
      qqj(kloop,33)=qj(kloop,33)*y(kloop,39)
!
!....  CF2Cl2 + hv = 2 Cl
!
      qqj(kloop,34)=qj(kloop,34)*y(kloop,40)
!
!....  CFC113 + hv = 3 Cl
!
      qqj(kloop,35)=qj(kloop,35)*y(kloop,41)
!
!....  CFC114 + hv = 2 Cl
!
      qqj(kloop,36)=qj(kloop,36)*y(kloop,42)
!
!....  CFC115 + hv = Cl
!
      qqj(kloop,37)=qj(kloop,37)*y(kloop,43)
!
!....  HCFC141b + hv = 2 Cl
!
      qqj(kloop,38)=qj(kloop,38)*y(kloop,47)
!
!....  HCFC142b + hv = Cl
!
      qqj(kloop,39)=qj(kloop,39)*y(kloop,48)
!
!....  CH3Br + hv = Br + CH3O2
!
      qqj(kloop,40)=qj(kloop,40)*y(kloop,38)
!
!....  CF3Br + hv = Br
!
      qqj(kloop,41)=qj(kloop,41)*y(kloop,49)
!
!....  CF2Br2 + hv = 2 Br
!
      qqj(kloop,42)=qj(kloop,42)*y(kloop,51)
!
!....  H2402 + hv = 2 Br
!
      qqj(kloop,43)=qj(kloop,43)*y(kloop,52)
!
!....  CF2ClBr + hv = Br + Cl
!
      qqj(kloop,44)=qj(kloop,44)*y(kloop,50)
      enddo
!
!.... End of kloop
!
      return
      end
