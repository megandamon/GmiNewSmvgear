!=======================================================================
!
! $Id: setkin_par.h,v 1.3 2013-07-31 15:08:30 ssteenro Exp $
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
!
!========1=========2=========3=========4=========5=========6=========7==


      integer  &
     &  NSP  &
     & ,NUM_J  &
     & ,NUM_K

      parameter (NSP   =  31)
      parameter (NUM_J =  1)
      parameter (NUM_K = 8)

!.... Species affiliations and assignments

      integer  &
     &  NACT  &
     & ,NCHEM  &
     & ,NCONST  &
     & ,NDYN  &
     & ,NFAM  &
     & ,NMF  &
     & ,NSS

      parameter (NACT   =  25)
      parameter (NCHEM  =  31)
      parameter (NCONST =   1)
      parameter (NDYN   =  30)
      parameter (NFAM   =   1)
      parameter (NMF    =  30)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer  &
     &  NSAD, NSADdust, NSADaer

      parameter (NSAD  =   0)
      parameter (NSADdust =   7)
      parameter (NSADaer =   7)


!.... Species individual identifications

      integer IC3H6
      parameter (IC3H6     =  0)

      integer IO3
      parameter (IO3       =  27)

      integer IH2O
      parameter (IH2O      =  0)

      integer IN2O
      parameter (IN2O      = 0)
!... Atomic Nitrogen
      integer IAN
      parameter (IAN       = 0)

      integer INO
      parameter (INO       = 0)

      integer INO2
      parameter (INO2      = 0)

      integer INO3
      parameter (INO3      = 30)

      integer IN2O5
      parameter (IN2O5     = 0)

      integer IHNO3
      parameter (IHNO3     = 0)

      integer IHCL
      parameter (IHCL      = 0)

      integer IHOCL
      parameter (IHOCL     = 0)

      integer ICLONO2
      parameter (ICLONO2   = 0)

      integer IHOBR
      parameter (IHOBR     = 0)

      integer IBRONO2
      parameter (IBRONO2   = 0)

      integer ICH3CCL3
      parameter (ICH3CCL3  = 0)

      integer ICCL4
      parameter (ICCL4     = 0)

      integer ICFCL3
      parameter (ICFCL3    = 0)

      integer ICF2CL2
      parameter (ICF2CL2   = 0)

      integer ICH3BR
      parameter (ICH3BR    = 0)

      integer IH1301
      parameter (IH1301    = 0)

      integer IH1211
      parameter (IH1211    = 0)

      integer ICO
      parameter (ICO       = 0)

      integer ICH4
      parameter (ICH4      = 0)

      integer IC3H6O
      parameter (IC3H6O    = 0)

      integer IC5H8
      parameter (IC5H8     = 0)

      integer IN2
      parameter (IN2       = 0)

      integer INITROGEN
      parameter (INITROGEN = 0)

      integer IOXYGEN
      parameter (IOXYGEN   = 0)

      integer IO2
      parameter (IO2       = 0)

      integer IMGAS
      parameter (IMGAS     = 31)

!.... Additional model-required parameters

      integer IH2O2
      parameter (IH2O2    =   1)

      integer IDEHYD
      parameter (IDEHYD   =   0)

      integer IH2OAIR
      parameter (IH2OAIR  =   0)

      integer ISYNOZ
      parameter (ISYNOZ   =   0)

      integer JBRONO21
      parameter (JBRONO21 =  0)

      integer KHNO31
      parameter (KHNO31   =   0)

      integer INDMS
      parameter (INDMS    =   4)

      integer INSO2
      parameter (INSO2    =   3)

      integer IOH
      parameter (IOH      =   28)

      integer IHO2
      parameter (IHO2     =   29)

      integer IFSO2
      parameter (IFSO2    =   2)

      integer IFBC
      parameter (IFBC     =   15)

      integer INOC
      parameter (INOC     =   13)

      integer IBBC
      parameter (IBBC     =   17)

      integer IBOC
      parameter (IBOC     =   16)

      integer IFOC
      parameter (IFOC     =   14)

      integer IFSO4A
      parameter (IFSO4A   =   5)

      integer INSO4A
      parameter (INSO4A   =   9)

      integer IDUST1
      parameter (IDUST1   =   18)

      integer IDUST2
      parameter (IDUST2   =   19)

      integer IDUST3
      parameter (IDUST3   =   20)

      integer IDUST4
      parameter (IDUST4   =   21)

      integer IDUST5
      parameter (IDUST5   =   22)

      integer ISSLT1
      parameter (ISSLT1   =   23)

      integer ISSLT2
      parameter (ISSLT2   =   24)

      integer ISSLT3
      parameter (ISSLT3   =   25)

      integer ISSLT4
      parameter (ISSLT4   =   26)

      integer INSO4S1
      parameter (INSO4S1   =   0)

      integer INSO4S2
      parameter (INSO4S2   =   0)

      integer INSO4S3
      parameter (INSO4S3   =   0)

      integer INSO4S4
      parameter (INSO4S4   =   0)

      integer INSO4D1
      parameter (INSO4D1   =   0)

      integer INSO4D2
      parameter (INSO4D2   =   0)

      integer INSO4D3
      parameter (INSO4D3   =   0)

      integer INSO4D4
      parameter (INSO4D4   =   0)

      integer INSO4N1
      parameter (INSO4N1   =   10)

      integer INSO4N2
      parameter (INSO4N2   =   11)

      integer INSO4N3
      parameter (INSO4N3   =   12)

      integer INSO4N4
      parameter (INSO4N4   =   0)

      integer IFSO4S1
      parameter (IFSO4S1   =   0)

      integer IFSO4S2
      parameter (IFSO4S2   =   0)

      integer IFSO4S3
      parameter (IFSO4S3   =   0)

      integer IFSO4S4
      parameter (IFSO4S4   =   0)

      integer IFSO4D1
      parameter (IFSO4D1   =   0)

      integer IFSO4D2
      parameter (IFSO4D2   =   0)

      integer IFSO4D3
      parameter (IFSO4D3   =   0)

      integer IFSO4D4
      parameter (IFSO4D4   =   0)

      integer IFSO4N1
      parameter (IFSO4N1   =   6)

      integer IFSO4N2
      parameter (IFSO4N2   =   7)

      integer IFSO4N3
      parameter (IFSO4N3   =   8)

      integer INO3AN1
      parameter (INO3AN1   =   0)

      integer INO3AN2
      parameter (INO3AN2   =   0)

      integer INO3AN3
      parameter (INO3AN3   =   0)

      integer INO3A
      parameter (INO3A     =   0)

      integer INH3
      parameter (INH3      =   0)

      integer INH4A
      parameter (INH4A     =   0)
!                                  --^--

      integer ICH2O, IH2, IHNO4, IMP, IALD2, IALK4, IC2H6
      integer IC3H8, IMACR, IMCO3, IMVK, IPMN, IR4N2, IRCHO
      parameter (ICH2O = 0, IH2   = 0, IHNO4 = 0, IMP   = 0)
      parameter (IALD2 = 0, IALK4 = 0, IC2H6 = 0, IC3H8 = 0)
      parameter (IMACR = 0, IMCO3 = 0, IMVK  = 0, IPMN  = 0)
      parameter (IR4N2 = 0, IRCHO = 0                      )


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
