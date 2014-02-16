!=======================================================================
!
! $Id: setkin_ibcb.h,v 1.2 2011-08-08 17:46:52 mrdamon Exp $
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
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
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
!                H2
      ibcb(5) = 1
!                HCl
      ibcb(24) = 1
!                CH3Cl
      ibcb(33) = 1
!                CH3CCl3
      ibcb(34) = 1
!                CCl4
      ibcb(35) = 1
!                CH3Br
      ibcb(38) = 1
!                CF3Br
      ibcb(39) = 1
!                CF2ClBr
      ibcb(40) = 1
!                                  --^--
