!=======================================================================
!
! $Id: setkin_surf_bc.h,v 1.1 2013-07-31 15:35:22 ssteenro Exp $
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
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_SO2_JPL06.db
!  Setkin files generated:  Thu Feb 17 20:08:54 2011
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

