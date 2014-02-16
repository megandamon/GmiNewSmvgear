!=======================================================================
!
! $Id: setkin_smv2par.h,v 1.5 2011-08-09 22:12:58 mrdamon Exp $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files generated:  Wed Feb  6 12:28:39 2008
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   =  84)
      parameter (SK_IPHOT  =  49)
      parameter (SK_ITHERM = 223)
      parameter (SK_NACT   =  79)

!                                  --^--

