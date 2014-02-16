!=======================================================================
!
! $Id: setkin_smv2par.h,v 1.1 2013-07-31 15:35:22 ssteenro Exp $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_SO2_JPL06.db
!  Setkin files generated:  Thu Feb 17 20:08:54 2011
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   = 135)
      parameter (SK_IPHOT  =  81)
      parameter (SK_ITHERM = 336)
      parameter (SK_NACT   = 131)

!                                  --^--

