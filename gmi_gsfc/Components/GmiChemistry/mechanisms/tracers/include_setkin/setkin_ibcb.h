!=======================================================================
!
! FILE
!   setkin_ibcb.h
!
! DESCRIPTION
!   This include file contains information about treatment of surface
!   boundary conditions.
!
!
!=======================================================================
!
!.... Set default boundary condition types
!
!.... Type 1 means fixed concentration
!.... surface boundary condition, Type 2 means
!.... fixed flux surface boundary condition
!
      ibcb(:)           = 0

      ibcb(NACT+1:IGAS) = 0
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
