!=======================================================================
!
! $Id: setkin_par.h,v 1.10 2013-07-31 15:08:57 ssteenro Exp $
!
! FILE
!   setkin_par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_SO2_JPL06.db
!  Setkin files generated:  Thu Nov 19 21:21:30 2009
!
!========1=========2=========3=========4=========5=========6=========7==


      integer &
     &  NSP &
     & ,NUM_J &
     & ,NUM_K

      parameter (NSP   = 124)
      parameter (NUM_J =  81)
      parameter (NUM_K = 321)

!.... Species affiliations and assignments

      integer &
     &  NACT &
     & ,NCHEM &
     & ,NCONST &
     & ,NDYN &
     & ,NFAM &
     & ,NMF &
     & ,NSS

      parameter (NACT   = 117)
      parameter (NCHEM  = 121)
      parameter (NCONST =   4)
      parameter (NDYN   = 120)
      parameter (NFAM   =   2)
      parameter (NMF    = 118)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer &
     &  NSAD, NSADdust, NSADaer

      parameter (NSAD  =   5)
      parameter (NSADdust =   7)
      parameter (NSADaer =   7)


!.... Species individual identifications

      integer IBBC
      parameter (IBBC      =   0)

      integer IBOC
      parameter (IBOC      =   0)

      integer IDUST1
      parameter (IDUST1    =   0)

      integer IDUST2
      parameter (IDUST2    =   0)

      integer IDUST3
      parameter (IDUST3    =   0)

      integer IDUST4
      parameter (IDUST4    =   0)

      integer IDUST5
      parameter (IDUST5    =   0)

      integer IFBC
      parameter (IFBC      =   0)

      integer IFOC
      parameter (IFOC      =   0)

      integer IFSO2
      parameter (IFSO2     =   0)

      integer IFSO4A
      parameter (IFSO4A    =   0)

      integer IFSO4D1
      parameter (IFSO4D1   =   0)

      integer IFSO4D2
      parameter (IFSO4D2   =   0)

      integer IFSO4D3
      parameter (IFSO4D3   =   0)

      integer IFSO4D4
      parameter (IFSO4D4   =   0)

      integer IFSO4N1
      parameter (IFSO4N1   =   0)

      integer IFSO4N2
      parameter (IFSO4N2   =   0)

      integer IFSO4N3
      parameter (IFSO4N3   =   0)

      integer IFSO4N4
      parameter (IFSO4N4   =   0)

      integer IFSO4S1
      parameter (IFSO4S1   =   0)

      integer IFSO4S2
      parameter (IFSO4S2   =   0)

      integer IFSO4S3
      parameter (IFSO4S3   =   0)

      integer IFSO4S4
      parameter (IFSO4S4   =   0)

      integer INDMS
      parameter (INDMS     =   0)

      integer INOC
      parameter (INOC      =   0)

      integer INSO2
      parameter (INSO2     =   0)

      integer INSO4A
      parameter (INSO4A    =   0)

      integer INSO4D1
      parameter (INSO4D1   =   0)

      integer INSO4D2
      parameter (INSO4D2   =   0)

      integer INSO4D3
      parameter (INSO4D3   =   0)

      integer INSO4D4
      parameter (INSO4D4   =   0)

      integer INSO4N1
      parameter (INSO4N1   =   0)

      integer INSO4N2
      parameter (INSO4N2   =   0)

      integer INSO4N3
      parameter (INSO4N3   =   0)

      integer INSO4N4
      parameter (INSO4N4   =   0)

      integer INSO4S1
      parameter (INSO4S1   =   0)

      integer INSO4S2
      parameter (INSO4S2   =   0)

      integer INSO4S3
      parameter (INSO4S3   =   0)

      integer INSO4S4
      parameter (INSO4S4   =   0)

      integer ISSLT1
      parameter (ISSLT1    =   0)

      integer ISSLT2
      parameter (ISSLT2    =   0)

      integer ISSLT3
      parameter (ISSLT3    =   0)

      integer ISSLT4
      parameter (ISSLT4    =   0)

      integer INO3An1
      parameter (INO3An1   =   0)

      integer INO3An2
      parameter (INO3An2   =   0)

      integer INO3An3
      parameter (INO3An3   =   0)

      integer ICH2O
      parameter (ICH2O     =   1)

      integer ICH4
      parameter (ICH4      =   2)

      integer ICO
      parameter (ICO       =   3)

      integer IH2
      parameter (IH2       =   5)

      integer IHNO2
      parameter (IHNO2     =   7)

      integer IHNO3
      parameter (IHNO3     =   8)

      integer IHNO4
      parameter (IHNO4     =   9)

      integer IH2O
      parameter (IH2O      =  10)

      integer IHO2
      parameter (IHO2      =  11)

      integer IH2O2
      parameter (IH2O2     =  12)

      integer IMP
      parameter (IMP       =  15)

      integer IAN
      parameter (IAN       =  16)

      integer IN2O
      parameter (IN2O      =  17)

      integer INO
      parameter (INO       =  18)

      integer INO2
      parameter (INO2      =  19)

      integer INO3
      parameter (INO3      =  20)

      integer IN2O5
      parameter (IN2O5     =  21)

      integer IO3
      parameter (IO3       =  24)

      integer IOH
      parameter (IOH       =  25)

      integer IBRCL
      parameter (IBRCL     =  27)

      integer IBRONO2
      parameter (IBRONO2   =  29)

      integer IHOBR
      parameter (IHOBR     =  31)

      integer ICL2
      parameter (ICL2      =  33)

      integer ICLO
      parameter (ICLO      =  34)

      integer ICL2O2
      parameter (ICL2O2    =  35)

      integer ICLONO2
      parameter (ICLONO2   =  36)

      integer IHCL
      parameter (IHCL      =  37)

      integer IHOCL
      parameter (IHOCL     =  38)

      integer IOCLO
      parameter (IOCLO     =  39)

      integer IALD2
      parameter (IALD2     =  58)

      integer IALK4
      parameter (IALK4     =  59)

      integer IC2H6
      parameter (IC2H6     =  62)

      integer IC3H8
      parameter (IC3H8     =  63)

      integer IC5H8
      parameter (IC5H8     =  80)

      integer IMACR
      parameter (IMACR     =  82)

      integer IMCO3
      parameter (IMCO3     =  87)

      integer IMVK
      parameter (IMVK      =  92)

      integer IPMN
      parameter (IPMN      =  95)

      integer IC3H6
      parameter (IC3H6     = 100)

      integer IR4N2
      parameter (IR4N2     = 103)

      integer IRCHO
      parameter (IRCHO     = 108)

      integer IC3H6O
      parameter (IC3H6O    = 118)

      integer INITROGEN
      parameter (INITROGEN = 119)

      integer IOXYGEN
      parameter (IOXYGEN   = 120)

      integer IMGAS
      parameter (IMGAS     = 121)

!.... Additional model-required parameters

      integer IN2
      parameter (IN2      = 119)

      integer IO2
      parameter (IO2      = 120)

      integer IDEHYD
      parameter (IDEHYD   = 122)

      integer IH2OAIR
      parameter (IH2OAIR  = 123)

      integer ISYNOZ
      parameter (ISYNOZ   = 124)

      integer JBRONO21
      parameter (JBRONO21 =   0)

      integer KHNO31
      parameter (KHNO31   =   0)

      integer INH3
      parameter (INH3     =   0)

      integer INO3A
      parameter (INO3A    =   0)

      integer INH4A
      parameter (INH4A    =   0)

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

      integer :: Istrat_O3=  0

!                                  --^--

