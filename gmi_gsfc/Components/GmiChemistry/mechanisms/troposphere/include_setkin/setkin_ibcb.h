!=======================================================================
!
! $Id: setkin_ibcb.h,v 1.5 2011-08-09 22:12:58 mrdamon Exp $
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
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
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
!                C2H6
!      ibcb(7) = 1
!                C3H8
!      ibcb(8) = 1
!                CO
!      ibcb(10) = 1
!                ISOP
!      ibcb(33) = 1
!                                  --^--
