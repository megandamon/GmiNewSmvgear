
      module fastJX_Parameters_mod

      PUBLIC

      INTEGER, PARAMETER :: IPAR = 1       ! Number of lon points in CTM
      INTEGER, PARAMETER :: JPAR = 1       ! Number of lat points in CTM
      INTEGER, PARAMETER :: LPAR_MAX = 75  ! max number of vertical levels
      INTEGER, PARAMETER :: JPNL_MAX = LPAR_MAX

#ifdef GTmodule
      INTEGER, PARAMETER :: JPPJ_MAX = 1     !  Number of photolytic reactions supplied
      INTEGER, PARAMETER :: NS = 61 ! Max # species which require J-values calculating
      INTEGER, PARAMETER :: NW = 30 ! Max # wavelength bins that can be used
      INTEGER, PARAMETER :: NP = 56 ! Max # aerosol/cloud types that can be used
      INTEGER, PARAMETER :: NH = 7  ! Max # Herzberg X-sections that can be used
      INTEGER, PARAMETER :: MX = 35 ! aerosol/cloud types supplied from CTM
#else
      INTEGER, PARAMETER :: JPPJ_MAX = 82    !  Number of photolytic reactions supplied
      INTEGER, PARAMETER :: NS = 61 ! Max # species which require J-values calculating
      INTEGER, PARAMETER :: NW = 18 ! Max # wavelength bins that can be used
      INTEGER, PARAMETER :: NP = 56 ! Max # aerosol/cloud types that can be used
      INTEGER, PARAMETER :: NH = 7  ! Max # Herzberg X-sections that can be used
      INTEGER, PARAMETER :: MX = 45 ! aerosol/cloud types supplied from CTM
#endif
      LOGICAL, PARAMETER :: LDEG45 = .false. ! Logical flag for degraded CTM resolution

      INTEGER, parameter :: NB_MAX = LPAR_MAX+1   ! Added for GMI  {AOO, 8/04}
      INTEGER, parameter :: NC_MAX = 2*NB_MAX     ! Added for GMI  {AOO, 8/04}

!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: NL  = 1200 ! Max # levels after insertion of extra Mie levels
      INTEGER, PARAMETER :: N__ = 2*NL ! # levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
      INTEGER, PARAMETER :: M__ = 4    ! #Gauss points used


      end module fastJX_Parameters_mod
