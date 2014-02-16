!=======================================================================
!
! $Id: setkin_par.h,v 1.5 2013-07-31 15:08:45 ssteenro Exp $
!
! FILE
!   setkin_par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!  "setkin_par.h" was modified for sulfur chemistry by Xiaohong Liu
!
!========1=========2=========3=========4=========5=========6=========7==


      integer  &
     &  NSP  &
     & ,NUM_J  &
     & ,NUM_K

      parameter (NSP   =  40)
      parameter (NUM_J =   1)
      parameter (NUM_K =   8)

!.... Species affiliations and assignments

      integer  &
     &  NACT  &
     & ,NCHEM  &
     & ,NCONST  &
     & ,NDYN  &
     & ,NFAM  &
     & ,NMF  &
     & ,NSS

      parameter (NACT   =  35)
      parameter (NCHEM  =  40)
      parameter (NCONST =   1)
      parameter (NDYN   =  39)
      parameter (NFAM   =   1)
      parameter (NMF    =  39)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer  &
     &  NSAD

      parameter (NSAD  =   0)


      integer IO3
      parameter (IO3 = 36)

      integer IOH
      parameter (IOH = 37)

      integer IHO2
      parameter (IHO2 = 38)

      integer INO3
      parameter (INO3 = 39)

      integer IMGAS
      parameter (IMGAS = 40)

      !mk
      integer IAN
      parameter (IAN       = 0)

!.... Species individual identifications

      integer IC3H6
      parameter (IC3H6     =  0)

      integer IH2O
      parameter (IH2O      =  0)

      integer IN2O
      parameter (IN2O      =  0)

      integer INO
      parameter (INO       =  0)

      integer INO2
      parameter (INO2      =  0)

      integer IN2O5
      parameter (IN2O5     =  0)

      integer IHNO3
      parameter (IHNO3     =  0)

      integer IHCL
      parameter (IHCL      =  0)

      integer IHOCL
      parameter (IHOCL     =  0)

      integer ICLONO2
      parameter (ICLONO2   =  0)

      integer IHOBR
      parameter (IHOBR     =  0)

      integer IBRONO2
      parameter (IBRONO2   =  0)

      integer ICH3CCL3
      parameter (ICH3CCL3  =  0)

      integer ICCL4
      parameter (ICCL4     =  0)

      integer ICFCL3
      parameter (ICFCL3    =  0)

      integer ICF2CL2
      parameter (ICF2CL2   =  0)

      integer ICH3BR
      parameter (ICH3BR    =  0)

      integer IH1301
      parameter (IH1301    =  0)

      integer IH1211
      parameter (IH1211    =  0)

      integer ICO
      parameter (ICO       =  0)

      integer ICH4
      parameter (ICH4      =  0)

      integer IC3H6O
      parameter (IC3H6O    =  0)

      integer IC5H8
      parameter (IC5H8     =  0)

      integer IN2
      parameter (IN2       =  0)

      integer INITROGEN
      parameter (INITROGEN =  0)

      integer IOXYGEN
      parameter (IOXYGEN   =  0)

      integer IO2
      parameter (IO2       =  0)

!.... Additional model-required parameters

      integer IDEHYD
      parameter (IDEHYD   =   0)

      integer IH2OAIR
      parameter (IH2OAIR  =   0)

      integer ISYNOZ
      parameter (ISYNOZ   =   0)

      integer JBRONO21
      parameter (JBRONO21 =   0)

      integer KHNO31
      parameter (KHNO31   =   0)

!.... Species placeholders for sulfchem solver

      integer IH2O2
      parameter (IH2O2    =   1)

      integer IFSO2
      parameter (IFSO2    =   2)

      integer INSO2
      parameter (INSO2    =   3)

      integer INDMS
      parameter (INDMS    =   4)

      integer ISO4G
      parameter (ISO4G    =   5)

      integer ISO4M1
      parameter (ISO4M1   =   6)

      integer ISO4N1
      parameter (ISO4N1   =   7)

      integer ISO4M2
      parameter (ISO4M2   =   8)

      integer ISO4N2
      parameter (ISO4N2   =   9)

      integer ISO4NOC
      parameter (ISO4NOC  =  10)

      integer ISO4FOC
      parameter (ISO4FOC  =  11)

      integer ISO4FBC
      parameter (ISO4FBC  =  12)

      integer ISO4BOC
      parameter (ISO4BOC  =  13)

      integer ISO4BBC
      parameter (ISO4BBC  =  14)

      integer ISO4D1
      parameter (ISO4D1   =  15)

      integer ISO4D2
      parameter (ISO4D2   =  16)

      integer ISO4D3
      parameter (ISO4D3   =  17)

      integer ISO4D4
      parameter (ISO4D4   =  18)

      integer ISO4S1
      parameter (ISO4S1   =  19)

      integer ISO4S2
      parameter (ISO4S2   =  20)

      integer ISO4S3
      parameter (ISO4S3   =  21)

      integer ISO4S4
      parameter (ISO4S4   =  22)

      integer INO3An1
      parameter (INO3An1     =   0)

      integer INO3An2
      parameter (INO3An2     =   0)

      integer INO3An3
      parameter (INO3An3     =   0)

      integer INOC
      parameter (INOC     =  23)

      integer IFOC
      parameter (IFOC     =  24)

      integer IFBC
      parameter (IFBC     =  25)

      integer IBOC
      parameter (IBOC     =  26)

      integer IBBC
      parameter (IBBC     =  27)

      integer IDUST1
      parameter (IDUST1   =  28)

      integer IDUST2
      parameter (IDUST2   =  29)

      integer IDUST3
      parameter (IDUST3   =  30)

      integer IDUST4
      parameter (IDUST4   =  31)

      integer ISSLT1
      parameter (ISSLT1   =  32)

      integer ISSLT2
      parameter (ISSLT2   =  33)

      integer ISSLT3
      parameter (ISSLT3   =  34)

      integer ISSLT4
      parameter (ISSLT4   =  35)

!                                  --^--
! for Jules' routine AerDust_Glob_To_Sub in gmi_glob_to_sub.F

      integer NSADaer
      parameter (NSADaer =   7)
      integer NSADdust
      parameter (NSADdust =   7)


      integer :: Iage     =  0

      integer :: Ie90     =  0

      integer :: Itm25    =  0

      integer :: IRn      =  0

      integer :: IPb      =  0

      integer :: IPbS     =  0

      integer :: IBe7     =  0

      integer :: IBe10    =  0

      integer :: IBe7S    =  0

      integer :: IBe10S   =  0

      integer :: ICH3I    =  0

      integer :: IfCO2    =  0

      integer :: ILINOZ   =  0

      integer :: ISF6     =  0

      integer :: ICLOCK   =  0

      integer :: Iuniform =  0 

