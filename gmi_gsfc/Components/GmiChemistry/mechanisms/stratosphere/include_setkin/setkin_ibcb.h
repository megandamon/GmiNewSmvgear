!=======================================================================
!
! $Id: setkin_ibcb.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_ibcb.h
!
! DESCRIPTION
!   This include file contains information about treatment of surface
!   boundary conditions.
!
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
!
!=======================================================================
!
!.... Set default boundary condition types
!
!.... Type 1 means fixed concentration
!.... surface boundary condition, Type 2 means
!.... fixed flux surface boundary condition
!
      ibcb(:)           = 2

      ibcb(NACT+1:IGAS) = 1
!
!.... Reset boundary condition type for special cases
!                N2O
      ibcb( 4) = 1
!                H2O
      ibcb(12) = 1
!                H2
      ibcb(17) = 1
!                CH4
      ibcb(18) = 1
!                CO
      ibcb(22) = 1
!                HCl
!cc      ibcb(28) = 1
!                CH3Cl
      ibcb(37) = 1
!                CH3Br
      ibcb(38) = 1
!                CFCl3
      ibcb(39) = 1
!                CF2Cl2
      ibcb(40) = 1
!                CFC113
      ibcb(41) = 1
!                CFC114
      ibcb(42) = 1
!                CFC115
      ibcb(43) = 1
!                HCFC22
      ibcb(44) = 1
!                CCl4
      ibcb(45) = 1
!                CH3CCl3
      ibcb(46) = 1
!                HCFC141b
      ibcb(47) = 1
!                HCFC142b
      ibcb(48) = 1
!                CF3Br
      ibcb(49) = 1
!                CF2ClBr
      ibcb(50) = 1
!                CF2Br2
      ibcb(51) = 1
!                H2402
      ibcb(52) = 1
!                                  --^--
