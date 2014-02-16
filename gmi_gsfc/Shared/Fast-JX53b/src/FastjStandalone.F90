!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
      program FastjStandalone
!
! !USES:
      use m_FastjMethod          , only : InitializeFastJ
      use m_FastjMethod          , only : RunFastJ
      use m_FastjMetVars         , only : XDGRD
      use m_FastjMetVars         , only : YDGRD
      use m_FastjJvaluesVars     , only : JLABEL
      use m_FastjReadNamelistFile, only : FastjReadNamelistFile
!
      implicit none
!
! !DESCRIPTION:
!  Routine to test Fast-J code by simulating model calls
!
!
! The following variables are required to initialize and run
! FastJ.
!
      integer              :: LongDim           ! Longitude dim of CTM grid
      integer              :: LatDim            ! latidude  dim of CTM grid
      integer              :: VertDim           ! altitude  dim of CTM grid
      integer              :: NumJval           ! num. photolysis reactions
      integer              :: VertDimJval       ! vertical(levels) dim for J-values
      real*8               :: GMTAU             ! Current time
      integer              :: iMONTH            ! Current month
      real*8               :: SurfAlbedo        ! Surface albedo
      real*8               :: SolarZenithAngle  ! Solar zenith angle
      integer              :: IndexCloud        ! Index Cloud 
      real*8 , allocatable :: PhotRate(:,:)     ! array of J's indexed to CTM chemistry
      real*8 , allocatable :: Temp(:)           ! Temperature
      real*8 , allocatable :: CloudOptDepth(:)  ! Cloud optical depth
      real*8 , allocatable :: LongGrid(:)       ! longitude grid
      real*8 , allocatable :: LatGrid (:)       ! latitude grid
      real*8 , allocatable :: EtaaMet(:)        ! Eta(a) value for level boundaries
      real*8 , allocatable :: EtabMet(:)        ! Eta(b) value for level boundaries
      real*8 , allocatable :: GridBoxArea(:,:)  ! Grid box area (m2)
      real*8 , allocatable :: SurfPressure(:,:) ! Surface pressure

      character (len=128)  :: CrossSection_InfileName     ! fast-J X-sections (spectral data) input file name
      character (len=128)  :: ScatteringData_InfileName   ! Aerosol/cloud scattering data input file name
      character (len=128)  :: PhotRateLabels_InfileName   ! Labels of photolysis rates required input file name
                                                          ! keyed to chem code
      character (len=128)  :: T_O3_Climatology_InfileName ! Read in T & O3 climatology input file name
                                                          ! general backup clim.
!
! The variables below are specific to the current test case.
!
      real*8               :: PressureIJ
      real*8               :: DELTAU
      integer              :: I, J, K, L
      integer              :: NSTEP
      integer              :: ILNG,JLAT,IDAY, IPH

      character (len=128)  :: namelist_file               ! Input namelist file
      character (len=128)  :: SampleRun_InfileName        ! Sample run input file name
      character (len=128)  :: CTMdata_InfileName          ! Input file setting up p, T, grid etc
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      write(6,*) ' fast-JX ver-5.3 standalone CTM code'

!!!!!
!!!!! Initialization Procedure
!!!!!
      ! Read in the namelist file.
      ! We read not only the name of the input files but also
      ! the variables setting the CTM dimensions.

      namelist_file = 'FastjNamelistFile'
      call FastjReadNamelistFile(namelist_file,          &
                             LongDim, Latdim, VertDim,   &
                             VertDimJval, NumJval,       &
                             CrossSection_InfileName,    &
                             ScatteringData_InfileName,  &
                             PhotRateLabels_InfileName,  &
                             T_O3_Climatology_InfileName,&
                             CTMdata_InfileName,         &
                             SampleRun_InfileName         )

     allocate( LongGrid     (LongDim))
     allocate( LatGrid      (LatDim))
     allocate( EtaaMet      (VertDim+1))
     allocate( EtabMet      (VertDim+1))
     allocate( Temp         (VertDim))
     allocate( CloudOptDepth(VertDim))
     allocate( GridBoxArea  (LongDim, LatDim))
     allocate( SurfPressure (LongDim, LatDim))
     allocate( PhotRate     (VertDimJval, NumJval))

!!!!!
!!!!! Read in CTM input files
!!!!!
      open(1, file=CTMdata_InfileName, status='old')
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) IDAY
          iMONTH=int(dble(IDAY)*12.d0/365.d0)+1    !  Approximately
      read(1,*) NSTEP
      read(1,*) GMTAU
      read(1,*) DELTAU
      read(1,*) JLAT
      read(1,*) ILNG
      read(1,*) PressureIJ
      read(1,*) SurfAlbedo
      read(1,*) IndexCloud
      read(1,*)
!--only for L_=37-layer EC
      do L=1,VertDim
        read(1,*) J,EtaaMet(L),EtabMet(L),Temp(L),CloudOptDepth(L)
      enddo
        read(1,*) J,EtaaMet(VertDim+1),EtabMet(VertDim+1)    !top layer
      close(1)

!! Read in annual mean pressure field
      IPH = 8
      open (IPH, file=SampleRun_InfileName, form='formatted',status='old')
      read (IPH,*)
      read (IPH,'(16F8.2)') (LongGrid(I),I=1,LongDim)
      read (IPH,*)
      read (IPH,'(16F8.2)') (LatGrid(J),J=1,LatDim)
      read (IPH,*)
      read (IPH,'(16F8.2)') (GridBoxArea (1,J),J=1,LatDim)
      read (IPH,*)
      read (IPH,'(16F8.1)') ((SurfPressure(I,J),I=1,LongDim),J=1,LatDim)
      close(IPH)

      do J=1,LatDim
        GridBoxArea (1,J) = GridBoxArea (1,J)*1.d9
        do I=2,LongDim
        GridBoxArea (I,J) = GridBoxArea (1,J)
        enddo
      enddo

      call InitializeFastJ &
                 ( LongDim, LatDim, VertDim, NumJval, VertDimJval, iMONTH, &
                   ILNG, JLAT, &
                   EtaaMet, EtabMet, GridBoxArea, LongGrid, LatGrid,       &
                   PressureIJ, &
                   Temp, CloudOptDepth, IndexCloud, SurfAlbedo, SurfPressure, &
                   CrossSection_InfileName  , ScatteringData_InfileName,   &
                   PhotRateLabels_InfileName, T_O3_Climatology_InfileName)

!!!!!
!!!!! Time stepping procedures
!!!!!

      GMTAU = GMTAU - DELTAU

      do I=1,NSTEP
         GMTAU = GMTAU + DELTAU

         if (GMTAU.ge.24.0) then
            IDAY=mod(IDAY,365)+1
            GMTAU=mod(GMTAU,24.d0)
            iMONTH=int(dble(IDAY)*12.d0/365.d0)+1   !do a better calendar 
         endif

         call  RunFastJ(ILNG, JLAT, IDAY, GMTAU, iMONTH, &
                          SolarZenithAngle, PhotRate)

         !---printout J's:
         write(6,1001) I,IDAY,GMTAU,SolarZenithAngle, YDGRD(JLAT),XDGRD(ILNG)
         write(6,1002) (JLABEL(K), K=1,NumJval)
         do L=VertDimJval,1,-1
            write(6,1003) L,(PhotRate(L,K),K=1,NumJval)
         enddo

 1001    format('Step=',i3,'  Day=',i3,' UT(hr)=',f4.1,'  SZA=',f7.2, &
     &                 '  LAT x LONG=',2f7.2)
 1002    format(1x,'L=  ',64(a7,2x))
 1003    format(i3,1p, 64e9.2)

      enddo
      stop

      end program FastjStandalone
!EOC
!-------------------------------------------------------------------------
