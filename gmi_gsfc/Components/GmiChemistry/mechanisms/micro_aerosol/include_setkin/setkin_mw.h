!=======================================================================
!
! $Id: setkin_mw.h,v 1.2 2011-08-09 22:02:42 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
!   "setkin_mw.h" was modified by Xiaohong Liu
!
! FILE
!   setkin_mw.h
!
! DESCRIPTION
!   This include file contains gram molecular weights
!   for the species in the model mechanism.
!
!  Chemistry input file generated: 2/23/00 10:12AM
!  Reaction dictionary:            3/29/99 2:55PM
!  Setkin files generated:         Wed Feb 23 10:30:41 2000
!
!=======================================================================

      real*8, save  :: mw_data(1:NSP)

      data mw_data(1:NSP) /  &
     &   34.01d0,   64.06d0,   64.06d0,   62.13d0,  &
     &   98.08d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   28.96d0,  &
     &   28.96d0,   28.96d0,   28.96d0,   48.00d0,  &
     &   17.01d0,   33.01d0,   62.00d0,   28.88d0 /

!   aerosol molecular weight is set to be MWTAIR,
!   so that mole_fraction = 1, and mass mixing ratio will be used.
!   unit of gas species: mol/mol,
!   unit of aerosol mass: kg/kg air, unit of aerosol number: #/kg air
