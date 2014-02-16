!=======================================================================
!
! $Id: setkin_par.h,v 1.6 2013-07-31 15:09:33 ssteenro Exp $
!
! FILE
!   setkin_par.h
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
     &  NSP &
     & ,NUM_J &
     & ,NUM_K

      parameter (NSP   =  85)
      parameter (NUM_J =  49)
      parameter (NUM_K = 223)

!.... Species affiliations and assignments

      integer &
     &  NACT &
     & ,NCHEM &
     & ,NCONST &
     & ,NDYN &
     & ,NFAM &
     & ,NMF &
     & ,NSS

      parameter (NACT   =  79)
      parameter (NCHEM  =  84)
      parameter (NCONST =   5)
      parameter (NDYN   =  80)
      parameter (NFAM   =   0)
      parameter (NMF    =  83)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer &
     &  NSAD, NSADdust, NSADaer

      parameter (NSAD  =   5)
      parameter (NSADdust =   7)
      parameter (NSADaer =   7)


!.... Species individual identifications

      integer IBBC
      parameter (IBBC    =  0)

      integer IAN
      parameter (IAN     =  0)

      integer IBOC
      parameter (IBOC    =  0)

      integer IDUST1
      parameter (IDUST1  =  0)

      integer IDUST2
      parameter (IDUST2  =  0)

      integer IDUST3
      parameter (IDUST3  =  0)

      integer IDUST4
      parameter (IDUST4  =  0)

      integer IFBC
      parameter (IFBC    =  0)

      integer IFOC
      parameter (IFOC    =  0)

      integer IFSO2
      parameter (IFSO2   =  0)

      integer IFSO4A
      parameter (IFSO4A  =  0)

      integer IFSO4D1
      parameter (IFSO4D1 =  0)

      integer IFSO4D2
      parameter (IFSO4D2 =  0)

      integer IFSO4D3
      parameter (IFSO4D3 =  0)

      integer IFSO4D4
      parameter (IFSO4D4 =  0)

      integer IFSO4N1
      parameter (IFSO4N1 =  0)

      integer IFSO4N2
      parameter (IFSO4N2 =  0)

      integer IFSO4N3
      parameter (IFSO4N3 =  0)

      integer IFSO4N4
      parameter (IFSO4N4 =  0)

      integer IFSO4S1
      parameter (IFSO4S1 =  0)

      integer IFSO4S2
      parameter (IFSO4S2 =  0)

      integer IFSO4S3
      parameter (IFSO4S3 =  0)

      integer IFSO4S4
      parameter (IFSO4S4 =  0)

      integer INDMS
      parameter (INDMS   =  0)

      integer INOC
      parameter (INOC    =  0)

      integer INSO2
      parameter (INSO2   =  0)

      integer INSO4A
      parameter (INSO4A  =  0)

      integer INSO4D1
      parameter (INSO4D1 =  0)

      integer INSO4D2
      parameter (INSO4D2 =  0)

      integer INSO4D3
      parameter (INSO4D3 =  0)

      integer INSO4D4
      parameter (INSO4D4 =  0)

      integer INSO4N1
      parameter (INSO4N1 =  0)

      integer INSO4N2
      parameter (INSO4N2 =  0)

      integer INSO4N3
      parameter (INSO4N3 =  0)

      integer INSO4N4
      parameter (INSO4N4 =  0)

      integer INSO4S1
      parameter (INSO4S1 =  0)

      integer INSO4S2
      parameter (INSO4S2 =  0)

      integer INSO4S3
      parameter (INSO4S3 =  0)

      integer INSO4S4
      parameter (INSO4S4 =  0)

      integer ISSLT1
      parameter (ISSLT1  =  0)

      integer ISSLT2
      parameter (ISSLT2  =  0)

      integer ISSLT3
      parameter (ISSLT3  =  0)

      integer ISSLT4
      parameter (ISSLT4  =  0)

      integer IALD2
      parameter (IALD2   =  3)

      integer IALK4
      parameter (IALK4   =  4)

      integer IC2H6
      parameter (IC2H6   =  7)

      integer IC3H8
      parameter (IC3H8   =  8)

      integer ICH2O
      parameter (ICH2O   =  9)

      integer ICO
      parameter (ICO     = 10)

      integer IH2O2
      parameter (IH2O2   = 19)

      integer IHNO2
      parameter (IHNO2   = 22)

      integer IHNO3
      parameter (IHNO3   = 23)

      integer IHNO4
      parameter (IHNO4   = 24)

      integer IHO2
      parameter (IHO2    = 25)

      integer IC5H8
      parameter (IC5H8   = 33)

      integer IMACR
      parameter (IMACR   = 35)

      integer IMCO3
      parameter (IMCO3   = 40)

      integer IMP
      parameter (IMP     = 45)

      integer IMVK
      parameter (IMVK    = 48)

      integer IN2O5
      parameter (IN2O5   = 50)

      integer INO
      parameter (INO     = 51)

      integer INO2
      parameter (INO2    = 52)

      integer INO3
      parameter (INO3    = 53)

      integer IO3
      parameter (IO3     = 54)

      integer IOH
      parameter (IOH     = 55)

      integer IPMN
      parameter (IPMN    = 57)

      integer IC3H6
      parameter (IC3H6   = 62)

      integer IR4N2
      parameter (IR4N2   = 65)

      integer IRCHO
      parameter (IRCHO   = 70)

      integer IC3H6O
      parameter (IC3H6O  = 80)

      integer ICH4
      parameter (ICH4    = 81)

      integer IH2
      parameter (IH2     = 82)

      integer IH2O
      parameter (IH2O    = 83)

      integer IMGAS
      parameter (IMGAS   = 84)

!.... Additional model-required parameters

      integer ISYNOZ
      parameter (ISYNOZ    =  85)

      integer IBRONO2
      parameter (IBRONO2   =   0)

      integer ICLONO2
      parameter (ICLONO2   =   0)

      integer IDEHYD
      parameter (IDEHYD    =   0)

      integer IH2OAIR
      parameter (IH2OAIR   =   0)

      integer IHCL
      parameter (IHCL      =   0)

      integer IHOBR
      parameter (IHOBR     =   0)

      integer IHOCL
      parameter (IHOCL     =   0)

      integer IN2O
      parameter (IN2O      =   0)

      integer INITROGEN
      parameter (INITROGEN =   0)

      integer IOXYGEN
      parameter (IOXYGEN   =   0)

      integer IN2
      parameter (IN2       =   0)

      integer IO2
      parameter (IO2       =   0)

      integer JBRONO21
      parameter (JBRONO21  =   0)

      integer KHNO31
      parameter (KHNO31    =   0)

      integer INH3
      parameter (INH3      =   0)

      integer INH4A
      parameter (INH4A     =   0)

      integer INO3A
      parameter (INO3A     =   0)

      integer INO3An1
      parameter (INO3An1     =   0)

      integer INO3An2
      parameter (INO3An2     =   0)

      integer INO3An3
      parameter (INO3An3     =   0)


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

!                                  --^--

