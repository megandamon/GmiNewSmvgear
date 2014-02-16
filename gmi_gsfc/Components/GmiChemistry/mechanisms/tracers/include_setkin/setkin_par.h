!=======================================================================
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!
!========1=========2=========3=========4=========5=========6=========7==


      integer &
     &  NSP &
     & ,NUM_J &
     & ,NUM_K

      parameter (NSP   =  18)
      parameter (NUM_J =   0)
      parameter (NUM_K =   0)

!.... Species affiliations and assignments

      integer &
     &  NACT &
     & ,NCHEM &
     & ,NCONST &
     & ,NDYN &
     & ,NFAM &
     & ,NMF &
     & ,NSS

      parameter (NACT   =  0)
      parameter (NCHEM  =  0)
      parameter (NCONST =   0)
      parameter (NDYN   =  18)
      parameter (NFAM   =   0)
      parameter (NMF    =   0)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer &
     &  NSAD, NSADdust, NSADaer

      parameter (NSAD  =   0)
      parameter (NSADdust =   0)
      parameter (NSADaer =   0)


!.... Species individual identifications

      integer, parameter :: Iage     =  1

      integer, parameter :: Ie90     =  2

      integer, parameter :: Itm25    =  3

      integer :: IRn      =  4

      integer :: IPb      =  5

      integer, parameter :: IPbS     =  6

      integer :: IBe7     =  7

      integer :: IBe10    =  8

      integer, parameter :: IBe7S    =  9

      integer, parameter :: IBe10S   =  10

      integer, parameter :: ICH3I    =  11

      integer, parameter :: IfCO2    =  12

      integer, parameter :: ILINOZ   =  13

      integer, parameter :: ISYNOZ   =  14

      integer, parameter :: ISF6     =  15

      integer, parameter :: ICLOCK   =  16

      integer, parameter :: Iuniform =  17 

      integer, parameter :: IStrat_O3=  18

      integer, parameter :: IO3     =  0 

      integer, parameter :: IOH     =  0

      integer, parameter :: IHO2    =  0

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

      integer IDUST5
      parameter (IDUST5  =  0)

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
      parameter (IALD2   =  0)

      integer IALK4
      parameter (IALK4   =  0)

      integer IC2H6
      parameter (IC2H6   =  0)

      integer IC3H8
      parameter (IC3H8   =  0)

      integer ICH2O
      parameter (ICH2O   =  0)

      integer ICO
      parameter (ICO     =  0)

      integer IH2O2
      parameter (IH2O2   =  0)

      integer IHNO2
      parameter (IHNO2   =  0)

      integer IHNO3
      parameter (IHNO3   =  0)

      integer IHNO4
      parameter (IHNO4   =  0)

      integer IC5H8
      parameter (IC5H8   =  0)

      integer IMACR
      parameter (IMACR   =  0)

      integer IMCO3
      parameter (IMCO3   =  0)

      integer IMP
      parameter (IMP     =  0)

      integer IMVK
      parameter (IMVK    =  0)

      integer IN2O5
      parameter (IN2O5   =  0)

      integer INO
      parameter (INO     =  0)

      integer INO2
      parameter (INO2    =  0)

      integer INO3
      parameter (INO3    =  0)

      integer IPMN
      parameter (IPMN    =  0)

      integer IC3H6
      parameter (IC3H6   =  0)

      integer IR4N2
      parameter (IR4N2   =  0)

      integer IRCHO
      parameter (IRCHO   =  0)

      integer IC3H6O
      parameter (IC3H6O  =  0)

      integer ICH4
      parameter (ICH4    =  0)

      integer IH2
      parameter (IH2     =  0)

      integer IH2O
      parameter (IH2O    =  0)

      integer IMGAS
      parameter (IMGAS   =  0)

!.... Additional model-required parameters

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
      parameter (INO3An1   =   0)

      integer INO3An2
      parameter (INO3An2   =   0)

      integer INO3An3
      parameter (INO3An3   =   0)

!                                --^--

