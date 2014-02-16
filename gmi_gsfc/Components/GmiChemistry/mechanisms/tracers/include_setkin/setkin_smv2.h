!=======================================================================
!
! FILE
!   setkin_smv2.h
!
! DESCRIPTION
!   This include file contains the primary description of the chemical
!   mechanism in the form of SMV2 data statements.
!
!
!=======================================================================
!
!.... Species indices for user within SMVGEAR
!
      IH2O      = -1
      INITROGEN = -1
      IOXYGEN   = -1
      IMGAS     = -1
!
!.... O2, N2, and M appear as reactants
!
      NMN2(1)   =  0
      NMO2(1)   =  0
      NMAIR(1)  =  0
      NM3BOD(1) =  0

      NREACN2(1, 1) = 0

      NREACO2(1, 1) = 0

      NREACAIR(:,:) = 0
      NREAC3B (:,:) = 0
      LGAS3BOD(:,:) = 0

!                                  --^--

