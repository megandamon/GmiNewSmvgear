!=======================================================================
!
! $Id: setkin_surf_bc.h,v 1.5 2011-08-09 22:12:58 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_surf_bc.h
!
! DESCRIPTION
!   This include file contains information about surface boundary
!   conditions.
!
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!=======================================================================

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -------------------------------------------------------
!     NUM_SBC : number of surface bounday conditions to reset
!     K_SBC   : max k to which species will be adjusted
!     -------------------------------------------------------

      integer, parameter :: NUM_SBC = 1

      integer, parameter :: K_SBC   = 2

!     ------------------
!     Integer variables.
!     ------------------

      integer, save :: sbc_map(NUM_SBC) = &
     &  (/ 0 /)

!                                  --^--

