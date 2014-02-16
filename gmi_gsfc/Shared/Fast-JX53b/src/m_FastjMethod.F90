!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_FastjMethod
!
! !INTERFACE:
      module m_FastjMethod
!
! !USES:
      use m_FastjMetVars         , only : MONTH
      use m_FastjMetVars         , only : PMEAN
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializeFastJ
      public  :: RunFastJ
      public  :: FinalizeFastJ
!
! !DESCRIPTION:
!  Routines to Initialize, Run and Finalize the package fastJ.
!
! !AUTHOR:
! Jules Kouatchou
!
! !HISTORY:
!
!EOP
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializeFastJ
!
! !INTERFACE:
!
      subroutine InitializeFastJ &
                 ( LongDim, LatDim, VertDim, NumJval, VertDimJval, iMONTH, &
                   ILNG, JLAT, &
                   EtaaMet, EtabMet, GridBoxArea, LongGrid, LatGrid,       &
                   PressureIJ, &
                   Temp, CloudOptDepth, IndexCloud, SurfAlbedo, SurfPressure, &
                   CrossSection_InfileName  , ScatteringData_InfileName,   &
                   PhotRateLabels_InfileName, T_O3_Climatology_InfileName)
!
! !USES:
      use m_FastjCTMparameters   , only : I_, J_, L_, L1_, L2_, JVL_, JVN_
      use m_FastjCTMparameters   , only : SetFastjCTMdims
      use m_FastjJvaluesVars     , only : InitializationFastjJvaluesVars
      use m_FastjMetVars         , only : InitializationFastjMetVars
      use m_FastjCTMspecific     , only : INPHOT
      use m_FastjMetVars         , only : TCNAME
      use m_FastjMetVars         , only : STT
      use m_FastjJvaluesVars     , only : RAD, ZZHT
      use m_FastjMetVars         , only : ETAA, ETAB, P, PMEAN
      use m_FastjMetVars         , only : AREAXY
      use m_FastjMetVars         , only : YGRD, YDGRD
      use m_FastjMetVars         , only : XGRD, XDGRD
      use m_FastjCTMspecific     , only : SET_CLD0
      use m_FastjCTMspecific     , only : SET_AER0
      use m_FastjCTMspecific     , only : SET_ATM0
      use m_FastjJvaluesVars     , only : InitializationFastjJvaluesVars

      implicit none
!
! !INPUT PARAMETERS:
      integer            , intent(in) :: LongDim           ! Longitude dim of CTM grid
      integer            , intent(in) :: LatDim            ! latidude  dim of CTM grid
      integer            , intent(in) :: VertDim           ! altitude  dim of CTM grid
      integer            , intent(in) :: NumJval           ! num. photolysis reactions 
      integer            , intent(in) :: VertDimJval       ! vertical(levels) dim for J-values
      integer            , intent(in) :: iMONTH, ILNG, JLAT
      real*8             , intent(in) :: LongGrid(LongDim)
      real*8             , intent(in) :: LatGrid (LatDim)
      real*8             , intent(in) :: EtaaMet(VertDim+1) !  Eta(a) value for level boundaries
      real*8             , intent(in) :: EtabMet(VertDim+1) !  Eta(b) value for level boundaries
      real*8             , intent(in) :: GridBoxArea(LongDim, LatDim)
      real*8             , intent(in) :: SurfPressure (LongDim, LatDim) ! surface pressure
      integer            , intent(in) :: IndexCloud                     ! Cloud Index
      real*8             , intent(in) :: SurfAlbedo                     ! Surface Albedo
      real*8             , intent(in) :: Temp         (VertDim)         ! Temperature
      real*8             , intent(in) :: CloudOptDepth(VertDim)         ! Cloud optical depth
      real*8             , intent(in) :: PressureIJ
      character (len=128), intent(in) :: CrossSection_InfileName
!                             ! fast-J X-sections (spectral data) input file name
      character (len=128), intent(in) :: ScatteringData_InfileName
!                             ! Aerosol/cloud scattering data input file name
      character (len=128), intent(in) :: PhotRateLabels_InfileName
!                             ! Labels of photolysis rates required input file name
!                             ! keyed to chem code
      character (len=128), intent(in) :: T_O3_Climatology_InfileName
!                             ! Read in T & O3 climatology input file name
!                             ! general backup clim.
!
! !DESCRIPTION:
!  Carry out the initialization of the package FastJ.
!
! !LOCAL VARIABLES:
      integer :: I, J
      real*8  :: PI
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      PI   = 3.141592653589793d0

      call SetFastjCTMdims (LongDim, Latdim, VertDim, VertDimJval, NumJval)

      ! Initialize MetFields related variables.
      call InitializationFastjMetVars ()
 
      ! Initialize Jvalues related variables.
      call InitializationFastjJvaluesVars()

! Defaults & constants
      RAD  = 6375.d5
      ZZHT = 5.d5
      STT(:,:,:,:) = 1.d6
      TCNAME(:)    = 'CO2'

      MONTH       = iMONTH
      PMEAN  (:,:) = SurfPressure(:,:)
      AREAXY (:,:) = GridBoxArea (:,:)
      XDGRD  (:)   = LongGrid(:)
      YDGRD  (:)   = LatGrid(:)
      ETAA   (:)   = EtaaMet(:)
      ETAB   (:)   = EtabMet(:)
      P(ILNG,JLAT) = PressureIJ

      do I = 1,I_
         XGRD(I) = XDGRD(I) *PI/180.d0
      enddo
      do J = 1,J_
         YGRD(J) = YDGRD(J) *PI/180.d0
      enddo
 
      call INPHOT (CrossSection_InfileName,     &
     &             ScatteringData_InfileName,   &
     &             PhotRateLabels_InfileName,   &
     &             T_O3_Climatology_InfileName)

      call SET_CLD0(Temp, CloudOptDepth, IndexCloud, SurfAlbedo)

      call SET_ATM0

      PMEAN(ILNG,JLAT) = P(ILNG,JLAT)

      call SET_AER0

      return

      end subroutine InitializeFastJ
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RunFastJ
!
! !INTERFACE:
!
      subroutine RunFastJ(ILNG, JLAT, IDAY, GMTAU, iMONTH, &
                          SolarZenithAngle, PhotRate)
!
! !USE:
      use m_FastjMetVars         , only : PMEAN, MONTH
      use m_FastjCTM_FastjLinking, only : PHOTOJ
      use m_FastjCTMspecific     , only : SET_ATM
      use m_FastjCTMspecific     , only : SET_AER
      use m_FastjCTMspecific     , only : SET_CLD
      use m_FastjCTMparameters   , only : L_, JVL_, JVN_

      implicit none
      
! !INPUT PARAMETERS:
      integer, intent(in   ) :: ILNG
      integer, intent(in   ) :: JLAT
      real*8 , intent(in   ) :: GMTAU                ! Current time
      integer, intent(in   ) :: IDAY                 ! Current day
      integer, intent(in   ) :: imonth               ! Current Month
!
! !OUPUT PARAMETERS:
      real*8 , intent(  out) :: PhotRate(JVL_,JVN_)  ! array of J's indexed to CTM chemistry
      real*8 , intent(  out) :: SolarZenithAngle
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

     MONTH = imonth

!---reset the atmosphere, aerosols, clouds for the time step (now = dummy)

     call SET_ATM(GMTAU)
     call SET_AER(GMTAU)
     call SET_CLD(GMTAU)

     call PHOTOJ(GMTAU,IDAY,ILNG,JLAT, SolarZenithAngle, PhotRate)

      return

      end subroutine RunFastJ
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  FinalizeFastJ
!
! !INTERFACE:
!
      subroutine FinalizeFastJ()
!
      implicit none
!
! !DESCRIPTION:
!  Finalize the package FastJ by deallocating
!  (if necessary) variables.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

     return

      end subroutine FinalizeFastJ
!EOC
!-----------------------------------------------------------------------

 end module m_FastjMethod
