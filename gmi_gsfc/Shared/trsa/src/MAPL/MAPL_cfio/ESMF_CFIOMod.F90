!==============================================================================
!BOP
! !MODULE: ESMF_CFIOMod.F90 - Source file for CFIO

       module ESMF_CFIOMod
!
! !DESCRIPTION:
!
! The code in this file provides data type definitions and interface 
! specifications
!
! This module provides all the necessary subroutines for users to write/read
! HDF format output using CF convention.
!
! !REVISION HISTORY:
!
!  Jan2004  Baoyu Yin  Initial design and prototyping.
!  Apr2004  Baoyu Yin  Implementation
!  Sep2004  Baoyu Yin  Modified return codes to make it more specific.
!  Sep2004  Baoyu Yin  Moved some utility routines to ESMF_CFIOUtil.F90.
!  Sep2004  Baoyu Yin  Modified station grid metadata.
!  Sep2004  Baoyu Yin  Added ptopUnit to ptop for eta and sigma coordinates.
!  Oct2004  Baoyu Yin  Migrated to Halem and fixed some bugs.
!  Oct2004  Baoyu Yin  Added timeString to ESMF_CFIOSet and ESMF_CFIOVarWrite.
!                      Rearranged the argument order in ESMF_CFIOVarWrite.
!  Jan2005  Baoyu Yin  Fixed some memory problems. Fixed scaleFactor and offset
!                      problem. Fixed standard_name problem in reading GFIO files.
!  Mar2005  Baoyu Yin  Moved some utility routines into ESMF_CFIOUtil.F90
!                      Modified error return codes.
!  Mar2005  Baoyu Yin  Added file name template
!  Apr2005  Baoyu Yin  Added time interpolation routine VarReadT
!  Apr2006  da Silva   Eliminated mpeu dependency.
!  Jun2006  Baoyu Yin  Added cyclic option for VarReadT
!  Jun2006  Baoyu Yin  Added reading 2D variable with VarReadT
!  Jul2006  da Silva   Eliminated read(str,fmt) to parse time; replaced
!                      with more robust mod() calculations.
!                      Made StrTemplate public.
!  Aug2006  da Silva   Added alternative interfaces to VarReadT.
!                      Included Baoyu patches in FileOpen() to handle double
!                      coordinate variables; previous merge of VarRead()
!                      and VarReadT() has been rolled back.
!  Dec2006  da Silva   Added ESMF_CFIODownBit() to downgrade precision for
!                      better gzipping.
!------------------------------------------------------------------------------
! !USES:
      use ESMF_CFIOUtilMod
      implicit none
!------------------------------------------------------------------------------
! !PRIVATE TYPES:
      private
!------------------------------------------------------------------------------
! !PUBLIC DATA TYPES:
!
      public :: ESMF_CFIO                ! A CFIO file object 
      public :: ESMF_CFIOVarInfo         ! A CFIO variable object
      public :: ESMF_CFIOGrid            ! A CFIO grid object

! !PUBLIC MEMBER FUNCTIONS:

      public :: ESMF_CFIOCreate          ! constructor for a CFIO object
      public :: ESMF_CFIOSet             ! set meta data for a CFIO object    
      public :: ESMF_CFIOGet             ! Get meta data
      public :: ESMF_CFIOFileCreate      ! Create a CFIO file for writing 
      public :: ESMF_CFIOFileOpen        ! Open a CFIO file 
      public :: ESMF_CFIOVarWrite        ! Write a variable to a file 
      public :: ESMF_CFIOVarRead         ! Read a variable from a file
      public :: ESMF_CFIOVarReadT        ! Read a variable from two files 
                                         ! with time interpolation
      public :: ESMF_CFIOFileClose       ! Close an existing CFIO file. 
      public :: ESMF_CFIODestroy         ! destructor for a CFIO object

      public :: ESMF_CFIOVarInfoCreate   ! constructor 
      public :: ESMF_CFIOVarInfoSet      ! set info for a CFIOVarInfo object
      public :: ESMF_CFIOVarInfoGet      ! Get info from a CFIOVarInfo object
      public :: ESMF_CFIOVarInfoDestroy  ! destructor 

      public :: ESMF_CFIOGridCreate      ! constructor 
      public :: ESMF_CFIOGridSet         ! set a CFIO grid
      public :: ESMF_CFIOGridGet         ! get a CFIO grid
      public :: ESMF_CFIOGridDestroy     ! destructor 

      public :: ESMF_CFIOstrTemplate     ! replacement for the one in mpeu

      public :: ESMF_CFIODownBit         ! Downgrade precision for better
                                         !  lossless compression

      interface ESMF_CFIOVarWrite; module procedure   &
        ESMF_CFIOVarWrite3D_,  &
        ESMF_CFIOVarWrite2D_,  &
        ESMF_CFIOVarWrite1D_
      end interface
                                                                                                 
      interface ESMF_CFIOVarRead; module procedure   &
        ESMF_CFIOVarRead3D_,  &
        ESMF_CFIOVarRead2D_,  &
        ESMF_CFIOVarRead1D_
      end interface

      interface ESMF_CFIOVarReadT; module procedure   &
        ESMF_CFIOVarReadT3D_,  &
        ESMF_CFIOVarReadT2D_,  &
        ESMF_CFIOVarReadT3D__, &
        ESMF_CFIOVarReadT2D__
      end interface

      interface ESMF_CFIOstrTemplate; module procedure   &
        strTemplate_
      end interface

      interface ESMF_CFIODownBit
        module procedure ESMF_CFIODownBit3D_ 
        module procedure ESMF_CFIODownBit2D_
      end interface

!
!EOP
!------------------------------------------------------------------------------

! Define a new data type "CFIO_Grid" -- contains grid information

      type ESMF_CFIOGrid
         private
         character(len=MVARLEN) :: gName ! name for this grid
         integer :: im              ! size of longitudinal dimension
         integer :: jm              ! size of latitudinal  dimension
         integer :: km              ! size of vertical dimension
         real, pointer :: lon(:)    ! longitude of center of gridbox in 
                                    ! degrees east of Greenwich (can be 
                                    ! -180 -> 180 or 0 -> 360)
         real, pointer :: lat(:)    ! latitude of center of gridbox in 
                                    ! degrees north of equator
         real, pointer :: lev(:)    ! Level (units given by levUnits) of
                                    ! center of gridbox
         character(len=MLEN) :: levUnits   ! units of level dimension, e.g.,
                                           ! "hPa", "sigma_level" 
         character(len=MLEN) :: coordinate ! string to indicate vertical coord
                                           ! (pressure, sigma, pressure_sigma)
         character(len=MLEN) :: standardName ! string for CF standard name
         character(len=MLEN) :: formulaTerm  ! string for CF formula terms
         real, pointer :: ak(:)     ! parameter for hybrid sigma prs coord.
         real, pointer :: bk(:)     ! parameter for hybrid sigma prs coord.
         real, pointer :: sigma(:)  ! parameter for sigma coordinate
         real :: ptop               ! parameter for sigma/eta coordinate
         character(len=MVARLEN) :: ptopUnit !  unit of ptop
         logical :: twoDimLat       ! support 2D lat/lon or not
         logical :: reduceGrid      ! support for reduced grid 
         logical :: stnGrid         ! support for station data
      end type ESMF_CFIOGrid 

! Define a new data type "CFIO_VarInfo" -- contains variable information

      type ESMF_CFIOVarInfo  
         private
         character(len=MVARLEN) :: vName       ! variable short name
         type(ESMF_CFIOGrid) :: grid           ! grid used for this var
         character(len=MLEN) :: vTitle         ! variable long name, e.g.,
                                               ! "Geopotential Height"
         character(len=MVARLEN):: vUnits       ! variable units, e.g., 
                                               ! "meter/second"
         real :: validRange(2)                 ! Variable valid range
         real :: packingRange(2) 
         real :: amiss                         ! Missing value such as 1.0E15
         real :: addOffSet                     ! optional
         real :: scaleFactor                   ! optional
         character(len=MLEN) :: standardName   ! optional, standard name 
                                               ! following CF convention
         logical :: twoDimVar     ! True for 2D; false for 3D
         logical :: timAve        ! True for time averaging file
         character :: aveMethod   ! 'c' for center averaging for time
                                  ! [-0.5*timeInc+time, 0.5*timeInc+time]
                                  ! Default: 'c'
                                  ! 'd' for downstream averaging 
                                  ! [time, time+timeInc]
                                  ! 'u' for upstream averaging 
                                  ! [time-timeInc, time]
         character(len=MVARLEN) :: cellMthd   ! Cell methmod 
         integer :: nVarAttInt    ! number of variable int attributes
         integer :: nVarAttReal   ! number of variable real attributes
         integer :: nVarAttChar   ! number of variable char attributes 
         integer, pointer :: attCharCnts(:)       ! length of char attributes
         integer, pointer :: attRealCnts(:)       ! length of real attributes
         integer, pointer :: attIntCnts(:)        ! length of int attributes
         character(len=MLEN), pointer :: attCharNames(:)! char attribute name
         character(len=MLEN), pointer :: attRealNames(:)! Real attribute name
         character(len=MLEN), pointer :: attIntNames(:) ! int attribute name
         character(len=MLEN), pointer :: varAttChars(:) ! char attributes
         real, pointer :: varAttReals(:,:)   ! User defined real attributes
         integer, pointer :: varAttInts(:,:) ! User defined integer attributes
         character(len=MVARLEN) :: ordering ! (time, lev, lat, lon) (default)
                                            ! can be any combination of xyzt
         type(iNode), pointer :: iList
         type(rNode), pointer :: rList
         type(cNode), pointer :: cList
      end type ESMF_CFIOVarInfo


! Define a new data type "ESMF_CFIO" -- a CFIO object(file) with file name, 
! CFIO variable objects, time, grid index and global attributes.

      type ESMF_CFIO
         private 
         character(len=MLEN) :: cfioObjName   ! name for this CFIO object
         character(len=MLEN) :: fName         ! file name in this CFIO obj. 
         character(len=MLEN) :: fNameTmplt    ! file name in this CFIO obj. 
         character(len=MLEN) :: expid         ! Experiment I
         integer :: mVars                     ! total number of variables
         type(ESMF_CFIOVarInfo), pointer :: varObjs(:) ! CFIO variable objects
         integer :: mGrids                    ! total number of grids 
         type(ESMF_CFIOGrid), pointer :: grids(:)     ! CFIO variable grid
         integer :: date                      ! yyyymmdd
         integer :: begTime                   ! hhmmss
         integer :: timeInc                   ! time step increment 
         integer :: tSteps                    ! total time steps
         character(len=MLEN) :: title    ! A title for the data set
         character(len=MLEN) :: source   ! Source of data, e.g. NASA/GMAO
         character(len=MLEN) :: contact  ! Who to contact about the data set
         character(len=MLEN) :: history  !
         character(len=MLEN) :: convention ! CFIO 
         character(len=MLEN) :: institution
         character(len=MLEN) :: references
         character(len=MLEN) :: comment

         integer :: nAttChar ! Number of char attributes
         integer :: nAttReal ! Number of Real attributes
         integer :: nAttInt  ! Number of int attributes
         integer, pointer :: attCharCnts(:)       ! length of char attributes
         integer, pointer :: attRealCnts(:)       ! length of real attributes
         integer, pointer :: attIntCnts(:)        ! length of int attributes
         character(len=MLEN), pointer :: attCharNames(:)! User defined char 
                                                       ! attribute name
         character(len=MLEN), pointer :: attRealNames(:)! Real attribute name
         character(len=MLEN), pointer :: attIntNames(:) ! int attribute name
         character(len=MLEN), pointer :: attChars(:) ! char attributes 
         real, pointer :: attReals(:,:)           ! global real attributes
         integer, pointer :: attInts(:,:)         ! global integer attributes
 
         integer :: prec                         ! Desired precision of data
         integer :: fid                          ! file ID for internal use
         type(iNode), pointer :: iList
         type(rNode), pointer :: rList
         type(cNode), pointer :: cList
         logical :: isOpen              ! flag to check fName is opened or not
         logical :: isCyclic            ! flag for cyclic for input files
      end type ESMF_CFIO

      contains

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOCreate -- ESMF_CFIO object constructor   

! !INTERFACE:
      type (ESMF_CFIO) function ESMF_CFIOCreate (cfioObjName, rc) 
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in), OPTIONAL :: cfioObjName  ! object name
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !DESCRIPTION:
!     Create a CFIO object and initialize vars . The required global metadata
!     title, institution, source, history, references, and comment are set to 
!     unknown.
!EOP
!------------------------------------------------------------------------------
      type(ESMF_CFIO) :: cfio                   ! a CFIO object
      integer :: rtcode

      if ( present(cfioObjName) ) then 
         cfio%cfioObjName = cfioObjName
      else
         cfio%cfioObjName = 'CFIO'
      end if

! Initializing variables

      cfio%nAttChar = 0
      cfio%nAttReal = 0
      cfio%nAttInt = 0

      cfio%fName = 'unknown'
      cfio%title = 'unknown'
      cfio%source = 'unknown'
      cfio%contact = 'unknown'
      cfio%history = 'unknown'
      cfio%convention = 'unknown'
      cfio%institution = 'unknown'
      cfio%references = 'unknown'
      cfio%comment = 'unknown'
      cfio%prec = 0
      cfio%date = -999
      cfio%begTime = 0
      cfio%timeInc = 60000
      cfio%mVars = 1
      cfio%mGrids = 1
      cfio%fNameTmplt = ''
      cfio%isOpen = .false.
      cfio%isCyclic = .false.
      cfio%expid = ''
      allocate(cfio%iList, cfio%rList, cfio%cList)
      nullify(cfio%iList)
      nullify(cfio%rList)
      nullify(cfio%cList)

      rtcode = 0

      if ( present(rc) ) rc = rtcode

      ESMF_CFIOCreate = cfio

      end function ESMF_CFIOCreate

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOSet  -- Set meta data for a CFIO object

! !INTERFACE:
      subroutine ESMF_CFIOSet(cfio, cfioObjName, varObjs, grids, grid,      &
                              fName, title, source, contact, history,       &
                              convention, institution, references, comment, &
                              date, begTime, timeInc, timeString, prec,     &
                              attCharNames, attCharCnts, attChars,          &
                              attRealNames, attRealCnts, attReals,          &
                              attIntNames, attIntCnts, attInts,             &
                              attCharName, attChar, attRealName, attReal,   &
                              attIntName, attInt, rc )
       implicit NONE

! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
       character(len=*), intent(in), OPTIONAL :: cfioObjName ! object name
       type(ESMF_CFIOVarInfo), OPTIONAL :: varObjs(:)! variable objects 
       type(ESMF_CFIOGrid), OPTIONAL :: grids(:)     ! grid array
       type(ESMF_CFIOGrid), OPTIONAL :: grid         

       character(len=*), intent(in), OPTIONAL :: fName      ! File name
       character(len=*), intent(in), OPTIONAL :: title      
       character(len=*), intent(in), OPTIONAL :: source     ! Source of data
       character(len=*), intent(in), OPTIONAL :: contact    ! Who to contact 
       character(len=*), intent(in), OPTIONAL :: history    !
       character(len=*), intent(in), OPTIONAL :: convention ! CFIO or COARDS
       character(len=*), intent(in), OPTIONAL :: institution! File name
       character(len=*), intent(in), OPTIONAL :: references
       character(len=*), intent(in), OPTIONAL :: comment

       integer, intent(in), OPTIONAL :: date          ! yyyymmdd
       integer, intent(in), OPTIONAL :: begTime       ! hhmmss
       integer, intent(in), OPTIONAL :: timeInc       ! time step increment
       character(len=*), intent(in), OPTIONAL :: timeString 
                                ! string expression of date and time  
       integer, intent(in), OPTIONAL :: prec      ! Desired precision of data:
                                                  ! 0 = 32 bit; 1 = 64 bit

       character(len=*), intent(in), OPTIONAL :: attCharNames(:) 
                                    ! User defined global char attribute names
       character(len=*), intent(in), OPTIONAL :: attRealNames(:)
                                    ! User defined global real attribute names 
       character(len=*), intent(in), OPTIONAL :: attIntNames(:)
                                    ! User defined global int attribute names
       integer, intent(in), OPTIONAL :: attCharCnts(:)! length of attributes
       integer, intent(in), OPTIONAL :: attRealCnts(:)! length of attributes
       integer, intent(in), OPTIONAL :: attIntCnts(:) ! length of attributes

       character(len=*), intent(in), OPTIONAL :: attChars(:) 
                                    ! User defined global char attribute 
       real,      intent(in), OPTIONAL :: attReals(:,:) 
                                    ! User defined global real attribute 
       integer,   intent(in), OPTIONAL :: attInts(:,:)
                                    ! User defined global int attribute 

       character(len=*), intent(in), OPTIONAL :: attCharName 
                                    ! User defined global char attribute name
       character(len=*), intent(in), OPTIONAL :: attRealName
                                    ! User defined global real attribute name
       character(len=*), intent(in), OPTIONAL :: attIntName
                                    ! User defined global int attribute name
       character(len=*), intent(in), OPTIONAL :: attChar 
                                    ! User defined global char attribute 
       real,    intent(in), OPTIONAL :: attReal(:)
                                    ! User defined global real attribute 
       integer, intent(in), OPTIONAL :: attInt(:)
                                    ! User defined global int attribute 
!
! !OUTPUT PARAMETERS:
!
       integer, intent(out), OPTIONAL :: rc 
                                    ! Error return code:
                                    ! 0   all is well
                                    ! -1  can't allocate memory for grid(s)
                                    ! -2  can't allocate memory: varObjs    
                                    ! -3  can't allocate mem: attIntCnts   
                                    ! -4  can't allocate mem: attIntNames  
                                    ! -5  can't allocate memory: attInts    
                                    ! -6  can't allocate mem: attRealCnts   
                                    ! -7  can't allocate mem: attRealNames  
                                    ! -8  can't allocate memory: attReals  
                                    ! -9  can't allocate mem: attCharCnts  
                                    ! -10  can't allocate mem: attCharNames  
                                    ! -11  can't allocate memory: attChars  
! !INPUT/OUTPUT PARAMETERS:
!
       type(ESMF_CFIO), intent(inout) :: cfio    ! a CFIO object
!
! !DESCRIPTION:
!     Set meta data for a CFIO object with detailed information. 
!EOP
!------------------------------------------------------------------------------
       integer :: iCnt, jCnt, count, rtcode

!      set required global meta data

       if ( present(cfioObjName) ) cfio%cfioObjName = cfioObjName
       if ( present(fName) ) cfio%fName = fName
       if ( present(title) ) cfio%title = title
       if ( present(source) ) cfio%source = source
       if ( present(contact) ) cfio%contact = contact
       if ( present(history) ) cfio%history = history
       if ( present(convention) ) cfio%convention = convention
       if ( present(institution) ) cfio%institution = institution
       if ( present(references) ) cfio%references = references
       if ( present(comment) ) cfio%comment = comment
       if ( present(date) ) cfio%date = date    
       if ( present(begTime) ) cfio%begTime = begTime 
       if ( present(timeInc) ) cfio%timeInc = timeInc 
       if ( present(timeString) ) then
          call strToInt(timeString, cfio%date, cfio%begTime)
       end  if
       if ( present(prec) ) cfio%prec = prec

!      set grid information
       if ( present(grids) ) then
          cfio%mGrids = size(grids) 
          allocate( cfio%grids(cfio%mGrids), stat = rtcode)
          if (err("can't allocate memory for grids",rtcode,-1) .lt. 0 ) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          cfio%grids = grids   
       end if
       if ( present(grid) ) then
          cfio%mGrids = 1
          allocate( cfio%grids(cfio%mGrids), stat = rtcode)
          if (err("can't allocate memory for grid",rtcode,-1) .lt. 0 ) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          cfio%grids = grid   
       end if

!      set variable
       if ( present(varObjs) ) then
          cfio%mVars = size(varObjs) 
          allocate( cfio%varObjs(cfio%mVars), stat = rtcode)
          if (err("can't allocate memory: varObjs",rtcode,-2) .lt. 0 ) then
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%varObjs = varObjs 
       end if

!      set integer names, counts and data
       if ( present(attIntCnts) )  then 
          allocate(cfio%attIntCnts(size(attIntCnts)), stat=rtcode)
          if (err("can't allocate mem: attIntCnts",rtcode,-3) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attIntCnts = attIntCnts 
          cfio%nAttInt = size(attIntCnts)
       end if
       if ( present(attIntNames) )  then 
          cfio%nAttInt = size(attIntNames)
          allocate(cfio%attIntNames(cfio%nAttInt), stat=rtcode)
          if (err("can't allocate mem: attIntNames",rtcode,-4) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attIntNames = attIntNames
       end if
       if ( present(attInts) )  then 
          iCnt = size(cfio%attIntCnts)
          jCnt = size(attInts)/size(cfio%attIntCnts)
          allocate(cfio%attInts(iCnt, jCnt), stat=rtcode)
          rtcode = err("can't allocate memory for attInts", rtcode, -1)
          if (err("can't allocate memory: attInts",rtcode,-5) .lt. 0) then   
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attInts= attInts
       end if

!      set real names, counts and data with array
       if ( present(attRealCnts) )  then 
          allocate(cfio%attRealCnts(size(attRealCnts)), stat=rtcode)
          rtcode = err("can't allocate memory for attRealCnts", rtcode, -1)
          if (err("can't allocate mem: attRealCnts",rtcode,-6) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attRealCnts = attRealCnts 
          cfio%nAttReal = size(attRealCnts)
       end if
       if ( present(attRealNames) )  then 
          cfio%nAttReal = size(attRealNames)
          allocate(cfio%attRealNames(cfio%nAttReal), stat=rtcode)
          if (err("can't allocate mem: attRealNames",rtcode,-7) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attRealNames = attRealNames
       end if
       if ( present(attReals) )  then 
          iCnt = size(cfio%attRealCnts)
          jCnt = size(attReals)/size(cfio%attRealCnts)
          allocate(cfio%attReals(iCnt, jCnt), stat=rtcode)
          if (err("can't allocate memory: attReals",rtcode,-8) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attReals= attReals
       end if

!      set character names, counts and data with array
       if ( present(attCharCnts) )  then 
          allocate(cfio%attCharCnts(size(attCharCnts)), stat=rtcode)
          if (err("can't allocate mem: attCharCnts",rtcode,-9) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attCharCnts = attCharCnts 
          cfio%nAttChar = size(attCharCnts)
       end if
       if ( present(attCharNames) )  then 
          cfio%nAttChar = size(attCharNames)
          allocate(cfio%attCharNames(cfio%nAttChar), stat=rtcode)
          if (err("can't allocate mem: attCharNames",rtcode,-10) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attCharNames = attCharNames
       end if
       if ( present(attChars) )  then 
          allocate(cfio%attChars(cfio%nAttChar), stat=rtcode)
          if (err("can't allocate memory: attChars",rtcode,-11) .lt. 0) then   
             if ( present(rc) ) rc = rtcode
             return
          end if

          cfio%attChars= attChars
       end if

!      set integer name, count and data into a list
       if ( present(attRealName) .and. present(attReal) ) then
          count = size(attReal)
          call addList(attRealName, count, attReal=attReal, &
                       rList=cfio%rList)
       end if

!      set real attribute name, count and data into a list
       if ( present(attIntName) .and. present(attInt) ) then
          count = size(attInt)
          call addList(attIntName, count, attInt=attInt, &
                       iList=cfio%iList)
       end if

!      set character attribute name, count and data into a list
       if ( present(attCharName) .and. present(attChar) ) then
          call addList(attCharName, len(attChar), attChar=attChar, &
                       cList=cfio%cList)
       end if

       rtcode = 0
       if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOSet


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOGet -- Get meta data from a CFIO object

! !INTERFACE:
      subroutine ESMF_CFIOGet (cfio, cfioObjName, nVars, varObjs, grid,     &
                              nGrids, grids, fName, title, source, contact, &
                              history, convention, institution, references, &
                              comment, date, begTime, timeInc, nSteps, prec,&
                              attCharNames, nAttChar, attCharCnts, attChars,&
                              attRealNames, nAttReal, attRealCnts, attReals,&
                              attIntNames, nAttInt, attIntCnts, attInts,    &
                              attCharName, attCharCnt, attChar, attRealName,&
                              attRealCnt, attReal, attIntName, attIntCnt,   &
                              attInt, rc )
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
       type(ESMF_CFIO), intent(in) :: cfio       ! a CFIO object
       character(len=*), intent(in), OPTIONAL :: attCharName
                                    ! User defined global char attribute name
       character(len=*), intent(in), OPTIONAL :: attRealName
                                    ! User defined global real attribute name
       character(len=*), intent(in), OPTIONAL :: attIntName
                                    ! User defined global int attribute name
!
! !OUTPUT PARAMETERS:
!
       character(len=*), intent(out), OPTIONAL :: cfioObjName ! CFIO Obj name
       integer, OPTIONAL :: nVars          ! number of variable objects
       type(ESMF_CFIOVarInfo), pointer, OPTIONAL :: varObjs(:)! var objects
       integer, OPTIONAL :: nGrids                   ! number of grids
       type(ESMF_CFIOGrid), pointer, OPTIONAL :: grids(:)    ! grid array
       type(ESMF_CFIOGrid), pointer, OPTIONAL :: grid
       character(len=*), intent(out), OPTIONAL :: fName      ! File name
       character(len=*), intent(out), OPTIONAL :: title
       character(len=*), intent(out), OPTIONAL :: source     ! Source of data
       character(len=*), intent(out), OPTIONAL :: contact    ! Who to contact
       character(len=*), intent(out), OPTIONAL :: history    !
       character(len=*), intent(out), OPTIONAL :: convention ! CFIO or COARDS
       character(len=*), intent(out), OPTIONAL :: institution! File name
       character(len=*), intent(out), OPTIONAL :: references
       character(len=*), intent(out), OPTIONAL :: comment
       integer, intent(out), OPTIONAL :: date          ! yyyymmdd
       integer, intent(out), OPTIONAL :: begTime       ! hhmmss
       integer, intent(out), OPTIONAL :: timeInc      ! time step increment
       integer, intent(out), OPTIONAL :: nSteps       ! number of time steps
       integer, intent(out), OPTIONAL :: prec     ! Desired precision of data:
                                                  ! 0 = 32 bit; 1 = 64 bit
       integer, intent(out), OPTIONAL :: nAttChar ! Number of char attributes
       integer, intent(out), OPTIONAL :: nAttReal ! Number of Real attributes
       integer, intent(out), OPTIONAL :: nAttInt  ! Number of int attributes
       character(len=*), pointer, OPTIONAL :: attCharNames(:)
                                    ! User defined global char attribute names
       character(len=*), pointer, OPTIONAL :: attRealNames(:)
                                    ! User defined global real attribute names
       character(len=*), pointer, OPTIONAL :: attIntNames(:)
                                    ! User defined global int attribute names
       integer, pointer, OPTIONAL :: attCharCnts(:)! length of attributes
       integer, pointer, OPTIONAL :: attRealCnts(:)! length of attributes
       integer, pointer, OPTIONAL :: attIntCnts(:) ! length of attributes
       character(len=*), pointer, OPTIONAL :: attChars(:)
                                    ! User defined global char attribute
       real,      pointer, OPTIONAL :: attReals(:,:)
                                    ! User defined global real attribute
       integer,   pointer, OPTIONAL :: attInts(:,:)
                                    ! User defined global int attribute

       integer, intent(out), OPTIONAL :: attIntCnt
       integer, intent(out), OPTIONAL :: attRealCnt
       integer, intent(out), OPTIONAL :: attCharCnt
       character(len=*), intent(out), OPTIONAL :: attChar
                                    ! User defined global char attribute
       real,    pointer, OPTIONAL :: attReal(:)
                                    ! User defined global real attribute
       integer, pointer, OPTIONAL :: attInt(:)
                                    ! User defined global int attribute

       integer, intent(out), OPTIONAL :: rc      ! Error return code:
                         !  0   all is well
                         ! -1  can't allocate memory for grid(s)
                         ! -2  can't allocate memory: varObjs
                         ! -3  can't allocate mem: attCharNames
                         ! -4  can't allocate mem: attRealNames
                         ! -5  can't allocate mem: attIntNames
                         ! -6  can't allocate mem: attCharCnts
                         ! -7  can't allocate mem: attRealCnts
                         ! -8  can't allocate mem: attIntCnts
                         ! -9  can't allocate mem: attChars
                         ! -10  can't allocate mem: attReals
                         ! -11  can't allocate mem: attInts
                         ! -12  can't allocate mem: attInt
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -47  error from ncdid or ncdinq (time)
                         !  rc = -48  error from ncinq
                         !  rc = -53  error from ncagtc/ncagt
!
! !DESCRIPTION:
!     Get meta data from a CFIO file
!EOP
!------------------------------------------------------------------------------
       integer :: rtcode
       integer :: i

       if ( present(cfioObjName) ) cfioObjName =cfio%cfioObjName

       if ( present(fName) ) fName =cfio%fName 
       if ( present(title) ) title = cfio%title
       if ( present(source) ) source = cfio%source 
       if ( present(contact) ) contact = cfio%contact
       if ( present(history) ) history = cfio%history 
       if ( present(convention) ) convention = cfio%convention
       if ( present(institution) ) institution = cfio%institution 
       if ( present(references) ) references = cfio%references 
       if ( present(comment) ) comment = cfio%comment 
       if ( present(date) ) date = cfio%date 
       if ( present(begTime) ) begTime = cfio%begTime
       if ( present(timeInc) ) timeInc = cfio%timeInc
       if ( present(prec) ) prec = cfio%prec
       if ( present(nSteps) ) nSteps = cfio%tSteps
       if ( present(nVars) ) nVars  = cfio%mVars   
       if ( present(nGrids) ) nGrids = cfio%mGrids  
       if ( present(grids) ) then
          allocate(grids(size(cfio%grids)), stat=rtcode)
          if (err("can't allocate memory for grids",rtcode,-1) .lt. 0 ) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          grids = cfio%grids
       end if
       if ( present(grid) ) then
          allocate(grid, stat=rtcode)
          if (err("can't allocate memory for grid",rtcode,-1) .lt. 0 ) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          grid = cfio%grids(1)
       end if
                                                                                     
       if ( present(varObjs) ) then
          allocate(varObjs(size(cfio%varObjs)), stat=rtcode)
          if (err("can't allocate memory: varObjs",rtcode,-2) .lt. 0 ) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          varObjs = cfio%varObjs
       end if

       if ( present(nAttChar) ) nAttChar = cfio%nAttChar
       if ( present(nAttReal) ) nAttReal = cfio%nAttReal
       if ( present(nAttInt) ) nAttInt = cfio%nAttInt

!      get global attribute names as an array.
       if ( present(attCharNames) ) then
          allocate(attCharNames(cfio%nAttChar), stat=rtcode)
          if (err("can't allocate mem: attCharNames",rtcode,-3) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attCharNames = cfio%attCharNames
       end if
       if ( present(attRealNames) ) then
          allocate(attRealNames(cfio%nAttReal), stat=rtcode)
          if (err("can't allocate mem: attRealNames",rtcode,-4) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attRealNames = cfio%attRealNames
       end if
       if ( present(attIntNames) ) then
          allocate(attIntNames(cfio%nAttInt), stat=rtcode)
          if (err("can't allocate mem: attIntNames",rtcode,-5) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attIntNames = cfio%attIntNames
       end if

!      get global attribute counts as an array.
       if ( present(attCharCnts) ) then
          allocate(attCharCnts(cfio%nAttChar), stat=rtcode)
          if (err("can't allocate mem: attCharCnts",rtcode,-6) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attCharCnts = cfio%attCharCnts
       end if
       if ( present(attRealCnts) ) then
          allocate(attRealCnts(cfio%nAttReal), stat=rtcode)
          if (err("can't allocate mem: attRealCnts",rtcode,-7) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attRealCnts = cfio%attRealCnts
       end if
       if ( present(attIntCnts) ) then
          allocate(attIntCnts(cfio%nAttInt), stat=rtcode)
          if (err("can't allocate mem: attIntCnts",rtcode,-8) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attIntCnts = cfio%attIntCnts
       end if

!      get global attributes as an array.
       if ( present(attChars) ) then
          allocate(attChars(cfio%nAttChar), stat=rtcode)
          if (err("can't allocate mem: attChars",rtcode,-9) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attChars= cfio%attChars
       end if
       if ( present(attReals) ) then
          allocate(attReals(cfio%nAttReal,size(cfio%attReals)/  &
                   cfio%nAttReal), stat=rtcode)
          if (err("can't allocate mem: attReals",rtcode,-10) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attReals= cfio%attReals
       end if
       if ( present(attInts) ) then
          allocate(attInts(cfio%nAttInt,size(cfio%attInts)/  &
                   cfio%nAttInt), stat=rtcode)
          if (err("can't allocate mem: attInts",rtcode,-11) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          attInts= cfio%attInts
       end if

!      provide attIntName and get its count and data
       if ( present(attIntName) ) then
          if ( present(attIntCnt) ) then
             do i = 1, cfio%nAttInt
                if (trim(attIntName) .eq. trim(cfio%attIntNames(i))) &
                    then
                   attIntCnt = cfio%attIntCnts(i)
                end if
             end do
          end if
          if ( present(attInt) ) then
             do i = 1, cfio%nAttInt
                if (trim(attIntName) .eq. trim(cfio%attIntNames(i)))&
                    then
                   allocate(attInt(cfio%attIntCnts(i)))
                   if (err("can't allocate mem: attInt",rtcode,-12) .lt. 0) &
                      then  
                      if ( present(rc) ) rc = rtcode
                      return
                   end if

                   attInt = cfio%attInts(i,1:cfio%attIntCnts(i))
                end if
             end do
          end if
       end if

!      provide attRealName and get its count and data
       if ( present(attRealName) ) then
          if ( present(attRealCnt) ) then
             do i = 1, cfio%nAttReal
                if (trim(attRealName) .eq. trim(cfio%attRealNames(i))) &
                    then
                   attRealCnt = cfio%attRealCnts(i)
                end if
             end do
          end if
          if ( present(attReal) ) then
             do i = 1, cfio%nAttReal
                if (trim(attRealName) .eq. trim(cfio%attRealNames(i)))&
                    then
                   allocate(attReal(cfio%attRealCnts(i)))
                   attReal = cfio%attReals(i,1:cfio%attRealCnts(i))
                end if
             end do
          end if
       end if

!      provide attCharName and get its count and data
       if ( present(attCharName) ) then
          if ( present(attCharCnt) ) then
             do i = 1, cfio%nAttChar
                if (trim(attCharName) .eq. trim(cfio%attCharNames(i))) &
                    then
                   attCharCnt = cfio%attCharCnts(i)
                end if
             end do
          end if
          if ( present(attChar) ) then
             do i = 1, cfio%nAttChar
                if (trim(attCharName) .eq. trim(cfio%attCharNames(i)))&
                    then
                   attChar = trim(cfio%attChars(i))
                end if
             end do
          end if
       end if

       rtcode = 0
       if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOGet


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOFileCreate -- Create a CFIO output file with meta data

! !INTERFACE:
      subroutine ESMF_CFIOFileCreate (cfio, rc, expid)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inout) :: cfio       ! a CFIO object
      character(len=*), intent(in), OPTIONAL  :: expid    ! Experiment ID
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                      ! 0   all is well
                      ! -1 Time increment is 0
                      ! -2  allocate memory error
                      ! -3  Num of int/char/real elements and Cnt don't match
                      ! -12  error determining default precision
                      ! -18 incorrect time increment
                      ! -30 can't open file
                      ! -31 error from ncddef
                      ! -32 error from ncvdef (dimension variable)
                      ! -33 error from ncapt(c) (dimension attribute)
                      ! -34 error from ncvdef (variable)
                      ! -35  error from ncapt(c) (variable attribute)
                      ! -36  error from ncaptc/ncapt (global attribute)
                      ! -37  error from ncendf
                      ! -38  error from ncvpt (dimension variable)
                      ! -39 Num of real var elements and Cnt differ
                      ! -55  error from ncredf (enter define mode)
                      ! -56  error from ncedf (exit define mode)
!
! !DESCRIPTION:
!     Create a CFIO output file with meta data
!EOP
!------------------------------------------------------------------------------
       integer :: i, n, rtcode
       integer :: maxLen
       character(len=MLEN) :: fNameTmp     ! file name 

!      checking file name template
       if (present(expid)) then 
          cfio%expid = expid
          call strTemplate_(fNameTmp,cfio%fName,xid=expid,nymd=cfio%date, &
                            nhms=cfio%begTime, stat=rtcode)
       else
          call strTemplate_(fNameTmp,cfio%fName,nymd=cfio%date, &
                            nhms=cfio%begTime, stat=rtcode)
       end if

       if (trim(fNameTmp) .ne. trim(cfio%fName)) then
          cfio%fNameTmplt = cfio%fName
          cfio%fName = fNameTmp
       end if

       call GFIO_Create_(cfio, rtcode)
       if (err("Error form GFIO_Create_",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

!      put global attributes
       call GFIO_PutCharAtt(cfio%fid, 'History', len(trim(cfio%history)),    &
                             cfio%history, rtcode )
       if (err("can't write History",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid, 'Source', len(trim(cfio%source)),      &
                             cfio%source, rtcode )
       if (err("can't write Source",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid, 'Title', len(trim(cfio%title)),        &
                             cfio%title, rtcode )
       if (err("can't write Title",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid, 'Contact', len(trim(cfio%contact)),    &
                             cfio%contact, rtcode )
       if (err("can't write Contact",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid,'Conventions',len(trim(cfio%convention))&
                             ,cfio%convention, rtcode )
       if (err("can't write Conventions",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid,'Institution',                          &
                            len(trim(cfio%institution)),                     &
                            cfio%institution, rtcode )
       if (err("can't write Institution",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid,'References',len(trim(cfio%references)),&
                             cfio%references, rtcode )
       if (err("can't write References",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if

       call GFIO_PutCharAtt(cfio%fid,'Comment',len(trim(cfio%comment)),      &
                             cfio%comment, rtcode )
       if (err("can't write Comment",rtcode,rtcode) .lt. 0) then  
          if ( present(rc) ) rc = rtcode
          return
       end if


!      get integer attributes from iList
       if ( associated(cfio%iList) ) then
          call getMaxLenCnt(maxLen, cfio%nAttInt, iList=cfio%iList)
          allocate(cfio%attIntNames(cfio%nAttInt),                           &
                   cfio%attIntCnts(cfio%nAttInt),                            &
                   cfio%attInts(cfio%nAttInt,maxLen), stat=rtcode)
          if (err("can't allocate mem: attIntCnts",rtcode,-2) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          call getList(iList=cfio%iList, intAttNames=cfio%attIntNames,       &
                       intAttCnts=cfio%attIntCnts, intAtts=cfio%attInts )
       end if

!      write user defined integer attributes
       if ( cfio%nAttInt .gt. 0 ) then
          do i = 1, cfio%nAttInt
             if ( cfio%attIntCnts(i) .gt. size(cfio%attInts(i,:)) )  then
                rtcode=err("FileCreate: Num of int elements and Cnt differ"  &
                            ,-3,-3)
                if ( present(rc) ) rc = rtcode
                return
             end if

             call GFIO_PutIntAtt(cfio%fid, cfio%attIntNames(i),              &
                                 cfio%attIntCnts(i), cfio%attInts(i,:),      &
                                 cfio%prec, rtcode )
             if (err("error in GFIO_PutIntAtt",rtcode,rtcode) .lt. 0) then
                if ( present(rc) ) rc = rtcode
                return
             end if

          end do
       end if

!      get real attributes from rList
       if ( associated(cfio%rList) ) then
          call getMaxLenCnt(maxLen, cfio%nAttReal, rList=cfio%rList)
          allocate(cfio%attRealNames(cfio%nAttReal),                       &
                   cfio%attRealCnts(cfio%nAttReal),                        &
                   cfio%attReals(cfio%nAttReal,maxLen), stat=rtcode)
          if (err("can't allocate mem: attRealNames",rtcode,-2) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if

          call getList(rList=cfio%rList, realAttNames=cfio%attRealNames,   &
                       realAttCnts=cfio%attRealCnts, realAtts=cfio%attReals )
          do i = 1, cfio%nAttReal
          end do
       end if

!      write user defined real attributes
       if ( cfio%nAttReal .gt. 0 ) then
          do i = 1, cfio%nAttReal
             if ( cfio%attRealCnts(i) .gt. size(cfio%attReals(i,:)) )  then
                rtcode=err("FileCreate: Num of real elements and Cnt differ" &
                            ,-3,-3)
                if ( present(rc) ) rc = rtcode
                return
             end if
             call GFIO_PutRealAtt(cfio%fid, cfio%attRealNames(i),            &
                                 cfio%attRealCnts(i),                        &
                                 cfio%attReals(i,1:cfio%attRealCnts(i)),     &
                                 cfio%prec, rtcode )
             if (err("error in GFIO_PutRealAtt",rtcode,rtcode) .lt. 0) then
                if ( present(rc) ) rc = rtcode
                return
             end if
          end do
       end if

!      get char attributes from cList
       if ( associated(cfio%cList) ) then
          call getMaxLenCnt(maxLen, cfio%nAttChar, cList=cfio%cList)
          allocate(cfio%attCharNames(cfio%nAttChar),                      &
                   cfio%attCharCnts(cfio%nAttChar),                       &
                   cfio%attChars(cfio%nAttChar), stat=rtcode)
          if (err("can't allocate mem: attCharNames",rtcode,-2) .lt. 0) then  
             if ( present(rc) ) rc = rtcode
             return
          end if
          call getList(cList=cfio%cList, charAttNames=cfio%attCharNames,  &
                       charAttCnts=cfio%attCharCnts, charAtts=cfio%attChars )
       end if

!      write user defined char attributes
       if ( cfio%nAttChar .gt. 0 ) then
          do i = 1, cfio%nAttChar
             call GFIO_PutCharAtt(cfio%fid, cfio%attCharNames(i),       &
                                 cfio%attCharCnts(i), cfio%attChars(i), &
                                 rtcode )
             if (err("error in GFIO_PutCharAtt",rtcode,rtcode) .lt. 0) then
                if ( present(rc) ) rc = rtcode
                return
             end if
          end do
       end if

       cfio%isOpen = .true.
 
       rtcode = 0
       if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOFileCreate

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOFileOpen -- open a CFIO file, and get CFIO meta data
!                                into a cfio Object.

! !INTERFACE:
      subroutine ESMF_CFIOFileOpen (cfio, fmode, rc, expid, cyclic)

!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      integer, intent(in) :: fmode              ! 0 for READ-WRITE
                                                ! non-zero for READ-ONLY
      character(len=*), intent(in), OPTIONAL :: expid   ! Experiment ID
      logical, intent(in), OPTIONAL :: cyclic           ! cyclic input file
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         ! -1   invalid count
                         ! -2   type mismatch
                         ! -12  error determining default precision
                         ! -10  ngatts is incompatible with file
                         ! -11  character string not long enough
                         ! -19  unable to identify coordinate variable
                         ! -36  error from ncaptc/ncapt (global attribute)
                         ! -39  error from ncopn (file open)
                         ! -40  error from ncvid
                         ! -41  error from ncdid or ncdinq (lat or lon)
                         ! -42  error from ncdid or ncdinq (lev)
                         ! -43  error from ncvid (time variable)
                         ! -47  error from ncdid or ncdinq (time)
                         ! -48  error from ncinq
                         ! -51  error from ncagtc/ncagt (global attribute)
                         ! -52  error from ncvinq
                         ! -53  error from ncagtc/ncagt
                         ! -57  error from ncanam
                         ! -58  error from ncainq

!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inout) :: cfio    ! a CFIO object
!
! !DESCRIPTION:
!     open a CFIO file, and get CFIO meta data into a cfio Object.
!EOP
!------------------------------------------------------------------------------
      integer :: ngatts, lm, i, ii, iv
      real*4 :: amiss
      real*4 :: vRange32(2)
      real*4, pointer :: lon(:), lat(:), lev(:)
      real*8, pointer :: lon_64(:), lat_64(:), lev_64(:)
      integer :: coXType = NCFLOAT
      integer :: coYType = NCFLOAT
      integer :: coZType = NCFLOAT
      character(len=MVARLEN) :: levunits
      character(len=MVARLEN) :: vAttName
      character(len=MVARLEN), pointer :: vname(:) 
      character(len=MLEN), pointer :: vtitle(:) 
      character(len=MVARLEN), pointer :: vunits(:) 
      integer, pointer :: kmvar(:)
      real, pointer :: valid_range(:,:), packing_range(:,:) 
      integer, pointer :: yyyymmdd(:), hhmmss(:)
      character(len=MLEN), pointer :: attNames(:)
      integer :: iCnt, rCnt, cCnt
      integer :: iMaxLen, rMaxLen, cMaxLen
      integer :: type, count, rtcode
      integer :: dimId
      integer :: varId
      integer :: datatype         ! variable type
      integer :: vtype            ! variable type
      integer :: nvDims           ! number of dimensions
      integer :: vDims(MAXVDIMS)  ! variable shape
      integer :: nvatts           ! number of attributes
      real*4, pointer :: rtmp(:)
      integer, pointer :: itmp(:)
      character(len=MVARLEN), pointer :: ctmp(:)
      logical :: esmf_file = .false.
      logical :: tmpLog
      logical :: new_grid
      integer :: nDims, allVars, recdim
      integer :: im, jm, km
      integer :: hour, min 
      integer :: fid, nVars, dimSize(4), myIndex
      character(len=MVARLEN) :: dimName(4), dimUnits(4), vnameTemp
      character(len=MVARLEN) :: nameAk, nameBk, namePtop
      integer :: loc1, loc2
      integer :: akid, bkid, ptopid
      integer :: icount
      real*4, pointer :: ak(:), bk(:)
      real*4 :: ptop
      real*4 :: scale, offset
      character, pointer ::  globalAtt(:)
      character(len=MLEN) :: fNameTmp     ! file name
  
      call ncpopt(0)

      fNameTmp = ''                                                                                   
!     checking file name template
      if (present(expid)) cfio%expid = expid
      if (present(cyclic)) cfio%isCyclic = cyclic
      if (present(expid) .and. cfio%date .gt. 0 .and. cfio%begTime .ge. 0) then
         call strTemplate_(fNameTmp,cfio%fName,xid=expid,nymd=cfio%date, &
                           nhms=cfio%begTime, stat=rtcode)
      else
         if (cfio%date .gt. 0 .and. cfio%begTime .ge. 0) then
            call strTemplate_(fNameTmp,cfio%fName,nymd=cfio%date, &
                              nhms=cfio%begTime, stat=rtcode)
         else   
            if (present(expid)) then
               call strTemplate_(fNameTmp,cfio%fName,xid=expid, stat=rtcode)
            end if
         end if
      end if
      if (trim(fNameTmp) .ne. trim(cfio%fName) .and. len(trim(fNameTmp)) .gt. 0) then
         cfio%fNameTmplt = cfio%fName
         cfio%fName = fNameTmp
      end if

!     open a cfio file
      call GFIO_Open ( cfio%fName, fmode, cfio%fid, rtcode )
      if (err("problem in GFIO_Open",rtcode,rtcode) .lt. 0 ) then
         if ( present(rc) ) rc = rtcode
         return
      end if
      cfio%isOpen = .true.
      fid =cfio%fid

!     get grid information and global meta data
                                                                                          
      call GFIO_DimInquire (cfio%fid, im, jm, km, lm, &
                            cfio%mVars, ngatts, rtcode)
      if (err("GFIO_DimInquire failed",rtcode,rtcode) .lt. 0) then  
         if ( present(rc) ) rc = rtcode
         return
      end if
      cfio%tSteps = lm

      call ncinq (cfio%fid,nDims,allVars,ngatts,recdim,rtcode)
      if (err("FileOpen: ncinq failed",rtcode,-48) .NE. 0) then  
         if ( present(rc) ) rc = rtcode
         return
      end if

      allocate(cfio%varObjs(cfio%mVars))
      nVars = 0
      cfio%mGrids = 0
      do i=1,allVars
        call ncvinq (fid,i,vnameTemp,vtype,nvDims,vDims,nvAtts,rtcode)
        if (err("Inquire: variable inquire error",rtcode,-52) .NE. 0) then  
           if ( present(rc) ) rc = rtcode
           return
        end if
        if (nvDims .EQ. 1 .and. (index(vnameTemp, 'lon') .gt. 0 .or.  &
            index(vnameTemp, 'XDim:EOSGRID') .gt. 0) ) then
           coXType = vtype
           cfio%mGrids = cfio%mGrids + 1
        end if
        if (nvDims .EQ. 1 .and. (index(vnameTemp, 'lat') .gt. 0 .or.  &
            index(vnameTemp, 'YDim:EOSGRID') .gt. 0) ) then
           coYType = vtype
        end if
        if (nvDims .EQ. 1 .and. (index(vnameTemp, 'lev') .gt. 0 .or.  &
            index(vnameTemp, 'Height:EOSGRID') .gt. 0) ) then
           coZType = vtype
        end if

        cfio%varObjs(nVars+1)%timAve = .false.
        if (trim(vnameTemp) .eq. 'time_bnds') then 
           cfio%varObjs(nVars)%timAve = .true.
           cycle
        end if
        if (nvDims .EQ. 1) cycle
        nVars = nVars + 1
        cfio%varObjs(nVars)%vName = trim(vnameTemp)
        cfio%varObjs(nVars)%grid%km = 0
        cfio%varObjs(nVars)%grid%stnGrid = .false.
        do iv = 1, nvDims
           call ncdinq(fid, vDims(iv), dimName(iv), dimSize(iv), rtcode)
           if (err("problem in ncdinq",rtcode,-41) .NE. 0) then  
              if ( present(rc) ) rc = rtcode
              return
           end if
           if (index(dimName(iv),'station') .gt. 0) then
              cfio%varObjs(nVars)%grid%im = dimSize(iv)
              cfio%varObjs(nVars)%grid%jm = dimSize(iv)
              cfio%varObjs(nVars)%grid%stnGrid = .true.
              cycle
           end if
           varId = ncvid (fid, dimName(iv), rtcode)
           dimUnits(iv) = ' '
           call ncagtc(fid,varId,'units',dimUnits(iv),MAXCHR,rtcode)
           if (err("problem in ncagtc",rtcode,-53) .NE. 0) then  
              if ( present(rc) ) rc = rtcode
              return
           end if
           myIndex = IdentifyDim (dimName(iv), dimUnits(iv))
           if (myIndex .EQ. 0) then
              cfio%varObjs(nVars)%grid%im = dimSize(iv)
              allocate(cfio%varObjs(nVars)%grid%lon(dimSize(iv)), &
                       lon(dimSize(iv)))
!              call ncvgt (fid, vDims(iv), 1, dimSize(iv), lon, rtcode)
              if ( coXType .eq. NCFLOAT ) then
                 call ncvgt (fid, varId, 1, dimSize(iv), lon, rtcode)
              else
                 allocate(lon_64(dimSize(iv)))
                 call ncvgt (fid, varId, 1, dimSize(iv), lon_64, rtcode)
                 lon =lon_64
                 deallocate(lon_64)
              end if
              if (err("problem in ncvgt",rtcode,-53) .NE. 0) then  
                 if ( present(rc) ) rc = rtcode
                 return
              end if
              cfio%varObjs(nVars)%grid%lon = lon
              deallocate(lon)
           end if
           if (myIndex .EQ. 1) then
              cfio%varObjs(nVars)%grid%jm = dimSize(iv)
              allocate(cfio%varObjs(nVars)%grid%lat(dimSize(iv)), &
                       lat(dimSize(iv)))
              if ( coYType .eq. NCFLOAT ) then
                 call ncvgt (fid, varId, 1, dimSize(iv), lat, rtcode)
              else
                 allocate(lat_64(dimSize(iv)))
                 call ncvgt (fid, varId, 1, dimSize(iv), lat_64, rtcode)
                 lat = lat_64
                 deallocate(lat_64)
              end if
!              call ncvgt (fid, vDims(iv), 1, dimSize(iv), lat, rtcode)
!print *, "vDims(iv) varId: ", vDims(iv), varId
!print *, "dimName dimUnits: ", trim(dimName(iv)), trim(dimUnits(iv))
              if (err("problem in ncvgt",rtcode,-51) .NE. 0) then  
                 if ( present(rc) ) rc = rtcode
                 return
              end if
              cfio%varObjs(nVars)%grid%lat = lat
              deallocate(lat)
           end if
           if (myIndex .EQ. 2) then
              cfio%varObjs(nVars)%grid%km = dimSize(iv)
              call ncpopt(0)
              call ncagtc(fid,varId,'standard_name',                   &
                          cfio%varObjs(nVars)%grid%standardName,           &
                          MAXCHR, rtcode)
              if (rtcode /= 0) cfio%varObjs(nVars)%grid%standardName="pressure"
              if ( index(cfio%varObjs(nVars)%grid%standardName,        &
                   'atmosphere_sigma_coordinate') .gt. 0  .or.         &
                   index(cfio%varObjs(nVars)%grid%standardName,        &
                   'atmosphere_hybrid_sigma_pressure_coordinate' )     &
                   .gt.  0 ) then

                 call ncagtc(fid,varId,'formula_term',                 &
                          cfio%varObjs(nVars)%grid%formulaTerm,            &
                          MAXCHR, rtcode)
                 if ( index(cfio%varObjs(nVars)%grid%standardName,     &
                   'atmosphere_sigma_coordinate') .gt. 0 ) then
                    loc1 = index(cfio%varObjs(nVars)%grid%formulaTerm,'ptop:')
                    icount = loc1  + 5
                    do icount = loc1+5, len(cfio%varObjs(nVars)%grid%formulaTerm)
                      if (cfio%varObjs(nVars)%grid%formulaTerm(icount:icount) &
                           .ne. ' ') exit
                    end do
                    namePtop=trim(cfio%varObjs(nVars)%grid%formulaTerm    &
                         (icount:len(cfio%varObjs(nVars)%grid%formulaTerm)))
                    ptopid = ncvid(cfio%fid, trim(namePtop), rtcode)
                    if (rtcode .ne. 0) print *, "problem in getting ptopid in ncvid"
                    if (rtcode .eq. 0) call ncvgt(cfio%fid,ptopid,1, 1, ptop, rtcode)
                    if (rtcode .eq. 0) cfio%varObjs(nVars)%grid%ptop = ptop
                 end if
              end if
              if (index(cfio%varObjs(nVars)%grid%standardName,             &
                        'atmosphere_hybrid_sigma_pressure_coordinate')     &
                        .gt. 0)  then
                 loc1 = index(cfio%varObjs(nVars)%grid%formulaTerm,'a:')
                 loc2 = index(cfio%varObjs(nVars)%grid%formulaTerm,'b:')
                 icount = 0
                 do icount = loc1+2, loc2
                   if (cfio%varObjs(nVars)%grid%formulaTerm(icount:icount) &
                        .ne. ' ') exit
                 end do
                 nameAk=trim(cfio%varObjs(nVars)%grid%formulaTerm          &
                             (icount:loc2-1))
                 loc1 = index(cfio%varObjs(nVars)%grid%formulaTerm,'b:')
                 loc2 = index(cfio%varObjs(nVars)%grid%formulaTerm,'ps:')
                 do icount = loc1+2, loc2
                   if (cfio%varObjs(nVars)%grid%formulaTerm(icount:icount) &
                        .ne. ' ') exit
                 end do
                 nameBk=trim(cfio%varObjs(nVars)%grid%formulaTerm          &
                             (icount:loc2-1))
                 loc1 = index(cfio%varObjs(nVars)%grid%formulaTerm,'p0:')
                 icount = loc1  + 4
                 namePtop=trim(cfio%varObjs(nVars)%grid%formulaTerm        &
                         (icount:len(cfio%varObjs(nVars)%grid%formulaTerm)))

                 akid = ncvid(cfio%fid, trim(nameAk), rtcode)
                 if (rtcode .ne. 0) print *, "problem in getting akid in ncvid"

                 allocate(cfio%varObjs(nVars)%grid%ak                      &
                          (cfio%varObjs(nVars)%grid%km+1),                 &
                          ak(cfio%varObjs(nVars)%grid%km+1))
                 call ncvgt(cfio%fid,akid,1,cfio%varObjs(nVars)%grid%km+1, &
                            ak, rtcode)
                 if (rtcode .ne. 0) print *, "problem in getting ak in ncvgt"
                 cfio%varObjs(nVars)%grid%ak = ak
                 deallocate(ak)
                 bkid = ncvid(cfio%fid, trim(nameBk), rtcode)
                 if (rtcode .ne. 0) print *, "problem in getting bkid in ncvid"
                 allocate(cfio%varObjs(nVars)%grid%bk                      &
                          (cfio%varObjs(nVars)%grid%km+1),                 &
                          bk(cfio%varObjs(nVars)%grid%km+1))
                 call ncvgt(cfio%fid,bkid,1,cfio%varObjs(nVars)%grid%km+1, &
                            bk, rtcode)
                 if (rtcode .ne. 0) print *, "problem in getting bk in ncvgt"
                 cfio%varObjs(nVars)%grid%bk = bk
                 deallocate(bk)

                 ptopid = ncvid(cfio%fid, trim(namePtop), rtcode)
                 if (rtcode .ne. 0) print *, "problem in getting ptopid in ncvid"
                 call ncvgt(cfio%fid,ptopid,1, 1, ptop, rtcode)
                 if (rtcode .ne. 0) print *, "problem in getting ptop in ncvgt"
                 cfio%varObjs(nVars)%grid%ptop = ptop
             end if
              call ncpopt(0)
              call ncagtc(fid,varId,'coordinate',                      &
                          cfio%varObjs(nVars)%grid%coordinate,             &
                          MAXCHR, rtcode)
              if (rtcode .ne. 0) cfio%varObjs(nVars)%grid%coordinate = "pressure"  
              cfio%varObjs(nVars)%grid%levUnits = trim(dimUnits(iv))

              allocate(cfio%varObjs(nVars)%grid%lev(dimSize(iv)), &
                       lev(dimSize(iv)))
              call ncpopt(0)
              if ( coZType .eq. NCFLOAT ) then
                 call ncvgt (fid, varId, 1, dimSize(iv), lev, rtcode) 
              else
                 allocate(lev_64(dimSize(iv)))
                 call ncvgt (fid, varId, 1, dimSize(iv), lev_64, rtcode) 
                 lev =lev_64
                 deallocate(lev_64)
              end if
              cfio%varObjs(nVars)%grid%lev = lev
              deallocate(lev)
           end if
        end do
        varId = ncvid (cfio%fid, cfio%varObjs(nVars)%vName, rtcode)
        if (rtcode .ne. 0) then 
           print *, "problem in getting varId in ncvid"
           if ( present(rc) ) rc = -40    
           return
        end if
        call ncagtc(fid,varId,'units',cfio%varObjs(nVars)%vunits,            &
                    MAXCHR,rtcode)
        if (rtcode .ne. 0) then
           print *, "ncagtc failed for units"
           if ( present(rc) ) rc = -53   
           return
        end if
        cfio%varObjs(nVars)%vtitle = ' '
        call ncpopt(0)
        call ncagtc(fid,varId,'long_name',cfio%varObjs(nVars)%vtitle,        &
                    MLEN,rtcode)
        if ( cfio%varObjs(nVars)%grid%km .gt. 0 ) then
            cfio%varObjs(nVars)%twoDimVar = .false.
        else
            cfio%varObjs(nVars)%twoDimVar = .true.
        end if
        call ncagt (fid, varId, '_FillValue', amiss, rtcode)
        if (rtcode .NE. 0) then
           call ncagt (fid, varId, 'missing_value', amiss, rtcode)
        end if
        cfio%varObjs(nVars)%amiss = amiss
        call ncpopt(0)
        call ncagt (fid, varId, 'scale_factor', scale, rtcode)
        if (rtcode .NE. 0) then
           cfio%varObjs(nVars)%scaleFactor = 1.0
        else
           cfio%varObjs(nVars)%scaleFactor = scale
        end if
        call ncpopt(0)
        call ncagt (fid, varId, 'add_offset', offset, rtcode)
        if (rtcode .NE. 0) then
           cfio%varObjs(nVars)%addOffset = 0.0
        else
           cfio%varObjs(nVars)%addOffset = offset
        end if
        call ncagt (fid, varId, 'vmin', vRange32(1), rtcode)
        if (rtcode .NE. 0) then
          cfio%varObjs(nVars)%validRange(1) = cfio%varObjs(nVars)%amiss
        else
          cfio%varObjs(nVars)%validRange(1) = vRange32(1)
        endif
        call ncagt (fid, varId, 'vmax', vRange32(2), rtcode)
        if (rtcode .NE. 0) then
          cfio%varObjs(nVars)%validRange(2) = cfio%varObjs(nVars)%amiss
        else
          cfio%varObjs(nVars)%validRange(2) = vRange32(2)
        endif
        
      end do
     
      call GetBegDateTime(fid,cfio%date,cfio%begTime,cfio%timeInc,rtcode)
      if (rtcode .ne. 0) then
         print *, "GetBegDateTime failed to get data/time/timeInc"
         if ( present(rc) ) rc = rtcode
         return
      end if

      hour = cfio%timeInc/3600
      min = mod(cfio%timeInc,3600*hour)/60
      cfio%timeInc = hour*10000 + min*100

      allocate(attNames(ngatts))
      call GFIO_GetAttNames ( cfio%fid, ngatts, attNames, rtcode )
      if (err("GFIO_GetAttNames failed",rtcode,rtcode) .lt. 0) then  
         if ( present(rc) ) rc = rtcode
         return
      end if
 
      iCnt = 0
      rCnt = 0
      cCnt = 0
      iMaxLen = 0
      rMaxLen = 0
      cMaxLen = 0

!     get how many int/real/char attributes in attNames
      do i =1, ngatts
         call GFIO_AttInquire (cfio%fid, attNames(i), type, count, rtcode)
         if (err("GFIO_AttInquire failed",rtcode,rtcode) .lt. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
         select case  (type)
            case ( 0 )
               iCnt = iCnt + 1
               if ( count .gt. iMaxLen ) iMaxLen = count
            case ( 1 )
               rCnt = rCnt + 1
               if ( count .gt. rMaxLen ) rMaxLen = count
            case ( 2 )
               cCnt = cCnt + 1
               if ( count .gt. cMaxLen ) cMaxLen = count
            case ( 3 )
               rCnt = rCnt + 1
               if ( count .gt. rMaxLen ) rMaxLen = count
            case ( 4 )
               iCnt = iCnt + 1
               if ( count .gt. iMaxLen ) iMaxLen = count
         end select
      end do

      cfio%nAttChar = cCnt
      cfio%nAttReal = rCnt
      cfio%nAttInt = iCnt

      allocate(cfio%attCharCnts(cCnt), cfio%attRealCnts(rCnt), &
               cfio%attIntCnts(iCnt))
      allocate(cfio%attCharNames(cCnt), cfio%attRealNames(rCnt), &
               cfio%attIntNames(iCnt))

      iCnt = 0
      rCnt = 0
      cCnt = 0
!     get attNames and count, then put them into a cfio obj
      do i =1, ngatts
         call GFIO_AttInquire (cfio%fid, attNames(i), type, count, rtcode)
         if (err("GFIO_AttInquire failed",rtcode,rtcode) .lt. 0) then  
            if ( present(rc) ) rc = rtcode
            return
         end if
         select case  (type)
            case ( 0 )
               iCnt = iCnt + 1
               cfio%attIntNames(iCnt) = attNames(i)         
               cfio%attIntCnts(iCnt) = count
            case ( 1 )
               rCnt = rCnt + 1
               cfio%attRealNames(rCnt) = attNames(i)
               cfio%attRealCnts(rCnt) = count
            case ( 2 )
               cCnt = cCnt + 1
               cfio%attCharNames(cCnt) = attNames(i)
               cfio%attCharCnts(cCnt) = count
            case ( 3 )
               rCnt = rCnt + 1
               cfio%attRealNames(rCnt) = attNames(i)
               cfio%attRealCnts(rCnt) = count
            case ( 4 )
               iCnt = iCnt + 1
               cfio%attIntNames(iCnt) = attNames(i)
               cfio%attIntCnts(iCnt) = count
         end select
      end do

      deallocate(attNames)

      allocate(cfio%attReals(rCnt, rMaxLen), cfio%attInts(iCnt, iMaxLen),    &
               cfio%attChars(cCnt))
!     get global integer attributes
      do i = 1, iCnt
         call GFIO_GetIntAtt(cfio%fid,cfio%attIntNames(i),cfio%attIntCnts(i) &
                            , cfio%attInts(i,:), rtcode)
         if (err("GFIO_GetIntAtt failed",rtcode,rtcode) .lt. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
      end do

!     get global real attributes
      do i = 1, rCnt
         call GFIO_GetRealAtt(cfio%fid,cfio%attRealNames(i),               &
                              cfio%attRealCnts(i),                         &
                              cfio%attReals(i,:), rtcode)
         if (err("GFIO_GetRealAtt",rtcode,rtcode) .lt. 0) then  
            if ( present(rc) ) rc = rtcode
            return
         end if
      end do
     
!     get global char attributes
      do i = 1, cCnt
         allocate(globalAtt(cfio%attCharCnts(i)))
         call GFIO_GetCharAtt(cfio%fid,cfio%attCharNames(i),   &
                              cfio%attCharCnts(i),             &
                              globalAtt, rtcode)
         if (err("GetCharAtt",rtcode,rtcode) .lt. 0) then  
            if ( present(rc) ) rc = rtcode
            return
         end if
!        cfio%attChars(i) can only hold MLEN characters.
         do ii = 1, cfio%attCharCnts(i)
            cfio%attChars(i)(ii:ii) = globalAtt(ii)
            if (ii .ge. MLEN) then
               print *,"global attribute ",trim(cfio%attCharNames(i)), &
                       " is longer than MLEN"
               exit
            end if
         end do
         cfio%attChars(i)(cfio%attCharCnts(i)+1:MLEN) = ' '
         if (index(cfio%attCharNames(i),'Conventions') .gt. 0 .and.  &
             index(cfio%attChars(i), 'ESMF') .gt. 0) esmf_file=.true.

         if (index(cfio%attCharNames(i),'History') .gt. 0)  &
            cfio%History=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'Source') .gt. 0)  &
            cfio%source=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'Title') .gt. 0)  &
            cfio%title=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'Contact') .gt. 0)  &
            cfio%contact=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'Conventions') .gt. 0)  &
            cfio%convention=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'Institution') .gt. 0)  &
            cfio%institution=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'References') .gt. 0)  &
            cfio%references=cfio%attChars(i)
         if (index(cfio%attCharNames(i),'Comment') .gt. 0)  &
            cfio%comment=cfio%attChars(i)
      end do


!     get variable meta data
      do i = 1, cfio%mVars
         varId = ncvid (cfio%fid, cfio%varObjs(i)%vName, rtcode)
         if (err("ncvid failed for vName",rtcode,rtcode) .lt. 0) then   
            if ( present(rc) ) rc = -40
            return
         end if
         call ncvinq(cfio%fid, varId, cfio%varObjs(i)%vName, datatype, &
                     nvdims, vdims, nvatts, rtcode)
         if (err("ncvinq failed for vName",rtcode,rtcode) .lt. 0) then  
            if ( present(rc) ) rc = -52
            return
         end if
         iCnt = 0
         rCnt = 0
         cCnt = 0
         iMaxLen = 0
         rMaxLen = 0
         cMaxLen = 0

!        get variable int/real/char attribute count
         do iv =1, nvatts
            call ncanam (cfio%fid, varId, iv, vAttName, rtcode)
            if (err("ncanam failed for vName",rtcode,rtcode) .lt. 0) then  
               if ( present(rc) ) rc = -57   
               return
            end if
            call ncainq (cfio%fid,varId,vAttName,vtype,count,rtcode)
            if (err("ncainq failed for vName",rtcode,rtcode) .lt. 0) then
               if ( present(rc) ) rc = -58   
               return
            end if
            select case  (vtype)
               case ( NCSHORT )
                  iCnt = iCnt + 1
                  if ( count .gt. iMaxLen ) iMaxLen = count
               case ( NCFLOAT )
                  rCnt = rCnt + 1
                  if ( count .gt. rMaxLen ) rMaxLen = count
               case ( NCCHAR )
                  cCnt = cCnt + 1
                  if ( count .gt. cMaxLen ) cMaxLen = count
               case ( NCDOUBLE )
                  rCnt = rCnt + 1
                  if ( count .gt. rMaxLen ) rMaxLen = count
               case ( NCLONG )
                  iCnt = iCnt + 1
                  if ( count .gt. iMaxLen ) iMaxLen = count
            end select
         end do
                                                                                            
         cfio%varObjs(i)%nVarAttChar = cCnt
         cfio%varObjs(i)%nVarAttReal = rCnt
         cfio%varObjs(i)%nVarAttInt = iCnt
                                                                                            
         allocate(cfio%varObjs(i)%attCharCnts(cCnt),  &
                  cfio%varObjs(i)%attRealCnts(rCnt),  &
                  cfio%varObjs(i)%attIntCnts(iCnt))     
         allocate(cfio%varObjs(i)%attCharNames(cCnt), &
                  cfio%varObjs(i)%attRealNames(rCnt),&
                  cfio%varObjs(i)%attIntNames(iCnt))

         iCnt = 0
         rCnt = 0
         cCnt = 0
!        get variable int/real/char attribute names and counts
         do iv =1, nvatts
            call ncanam (cfio%fid, varId, iv, vAttName, rtcode)
            if (err("ncanam failed for vName",rtcode,rtcode) .lt. 0) then  
               if ( present(rc) ) rc = -57
               return
            end if
            call ncainq (cfio%fid,varId,vAttName,vtype,count,rtcode)
            if (err("ncainq failed for vName",rtcode,rtcode) .lt. 0) then   
               if ( present(rc) ) rc = -58
               return
            end if
            select case  (vtype)
               case ( NCSHORT )
                  iCnt = iCnt + 1
                  cfio%varObjs(i)%attIntNames(iCnt) = vAttName
                  cfio%varObjs(i)%attIntCnts(iCnt) = count   
               case ( NCFLOAT )
                  rCnt = rCnt + 1
                  cfio%varObjs(i)%attRealNames(rCnt) = vAttName
                  cfio%varObjs(i)%attRealCnts(rCnt) = count   
               case ( NCCHAR )
                  cCnt = cCnt + 1
                  cfio%varObjs(i)%attCharNames(cCnt) = vAttName
                  cfio%varObjs(i)%attCharCnts(cCnt) = count   
               case ( NCDOUBLE )
                  rCnt = rCnt + 1
                  cfio%varObjs(i)%attRealNames(rCnt) = vAttName
                  cfio%varObjs(i)%attRealCnts(rCnt) = count   
               case ( NCLONG )
                  iCnt = iCnt + 1
                  cfio%varObjs(i)%attIntNames(iCnt) = vAttName
                  cfio%varObjs(i)%attIntCnts(iCnt) = count   
            end select
         end do
   
         allocate(cfio%varObjs(i)%varAttReals(rCnt, rMaxLen), &
                  cfio%varObjs(i)%varAttInts(iCnt, iMaxLen),  &
                  cfio%varObjs(i)%varAttChars(cCnt))

!        get int variable attributes
         do ii = 1, iCnt
            allocate(itmp(cfio%varObjs(i)%attIntCnts(ii)))
            call ncagt(cfio%fid,varId,cfio%varObjs(i)%attIntNames(ii),&
                       itmp, rtcode)
            if (err("ncagt failed for attIntNames",rtcode,rtcode) .lt. 0) then
               if ( present(rc) ) rc = -53   
               return
            end if
            cfio%varObjs(i)%varAttInts(ii,1:cfio%varObjs(i)%attIntCnts(ii))&
                       = itmp
            deallocate(itmp)
         end do

!        get real variable attributes
         do ii = 1, rCnt
            allocate(rtmp(cfio%varObjs(i)%attRealCnts(ii)))
            call ncagt(cfio%fid,varId,cfio%varObjs(i)%attRealNames(ii),      &
                       rtmp, rtcode)
            if (err("ncagt failed for attRealNames",rtcode,rtcode) .lt. 0) then
               if ( present(rc) ) rc = -53
               return
            end if
            cfio%varObjs(i)%varAttReals(ii,1:cfio%varObjs(i)%attRealCnts(ii))&
                       = rtmp
            deallocate(rtmp)
         end do

!        get char variable attributes
         do ii = 1, cCnt
            call ncagtc(cfio%fid,varId,cfio%varObjs(i)%attCharNames(ii),     &
                       cfio%varObjs(i)%varAttChars(ii),                      &
                       cfio%varObjs(i)%attCharCnts(ii), rtcode)
            if (err("ncagt failed for attCharNames",rtcode,rtcode) .lt. 0) then
               if ( present(rc) ) rc = -53   
               return
            end if
            cfio%varObjs(i)%varAttChars(ii)  &
                 (cfio%varObjs(i)%attCharCnts(ii)+1:MLEN) = ' '             
         end do

      end do

!     set grids objects in a CFIO object
      allocate( cfio%grids(cfio%mGrids), stat = rtcode)
      cfio%grids(1) = cfio%varObjs(1)%grid
      if ( cfio%mGrids .eq. 1 .and. cfio%varObjs(1)%grid%km .eq. 0) &
         cfio%grids(1)%km = km
      
      if ( cfio%mGrids .gt. 1 ) then
        do i = 2, cfio%mGrids
           iCnt = 1
           do iv = 2, cfio%mVars
              new_grid = .true.
              iCnt = iCnt + 1
              do ii = 2, i
                if (cfio%varObjs(iv)%grid%im .eq. cfio%grids(ii-1)%im .and.  &
                  cfio%varObjs(iv)%grid%jm .eq. cfio%grids(ii-1)%jm .and.  &
                  cfio%varObjs(iv)%grid%km .eq. cfio%grids(ii-1)%km ) then 
                  new_grid = .false.
                end if
              end do
              if ( new_grid ) exit
           end do
           cfio%grids(i) = cfio%varObjs(iCnt)%grid
        end do
      end if 

      rtcode = 0
      if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOFileOpen

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarWrite3D_ -- Write a variable to a output file

! !INTERFACE:
      subroutine ESMF_CFIOVarWrite3D_(cfio, vName, field, date, curTime, &
                                      kbeg, kount, timeString, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj  
      character(len=*), intent(in) :: vName       ! Variable name  
      real, intent(in) :: field(:,:,:)            ! array contains data
      integer, intent(in), OPTIONAL :: date       ! yyyymmdd
      integer, intent(in), OPTIONAL :: curTime    ! hhmmss
      integer, intent(in), OPTIONAL :: kbeg       ! first level to write
      integer, intent(in), OPTIONAL :: kount      ! number of levels to write
      character(len=*), intent(in), OPTIONAL :: timeString
                                  ! string expression for date and time


!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -15  data outside of valid range
                         !  rc = -16  data outside of packing range
                         !  rc = -17  data outside of pack and valid range
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -45  error from ncvpt
                         !  rc = -46  error from ncvgt
                         !  rc = -52  error from ncvinq
                         !  rc = -53  error from ncagtc/ncagt

!
! !DESCRIPTION:
!     Write a variable to file
!EOP
!------------------------------------------------------------------------------
      integer :: i, rtcode
      integer :: myKbeg, myKount
      integer :: myDate, myCurTime
      character(len=MLEN) :: fNameTmp     ! file name 
                                                                                         
      fNameTmp = ''
      if ( present(date) ) myDate = date
      if ( present(curTime) ) myCurTime = curTime
      if ( present(timeString) ) call strToInt(timeString,myDate,myCurTime)

      if (len(trim(cfio%fNameTmplt)) .gt. 1) then
         call strTemplate_(fNameTmp,cfio%fNameTmplt,xid=cfio%expid,nymd=myDate, &
                           nhms=myCurTime, stat=rtcode)
         if (trim(fNameTmp) .ne. trim(cfio%fName)) then
            call ESMF_CFIOFileClose(cfio)
            cfio%fName = fNameTmp
            call ESMF_CFIOSet(cfio, fName=cfio%fName)
            call ESMF_CFIOSet(cfio, date=myDate, begTime=myCurTime)
            if (len(trim(cfio%expid)) .gt. 0) then
               call ESMF_CFIOFileCreate(cfio, expid=cfio%expid)
            else
               call ESMF_CFIOFileCreate(cfio)
            end if
         end if
      end if
          
!
!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do

!     write 2D variable
      if ( cfio%varObjs(i)%twoDimVar ) then 
         call GFIO_PutVar (cfio%fid, vName, myDate, myCurTime,             &
                        cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,  &
                        0, 1, field, rtcode )
         if (err("GFIO_PutVar failed",rtcode,rtcode) .lt. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
!     write 3D variable
      else
         myKbeg = 1
         myKount = cfio%varObjs(i)%grid%km

         if ( present(kbeg) ) myKbeg = kbeg 
         if ( present(kount) ) myKount = kount

         call GFIO_PutVar (cfio%fid, vName, myDate, myCurTime,             &
                        cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,  &
                        myKbeg, myKount, field, rtcode )
         if (err("GFIO_PutVar failed",rtcode,rtcode) .lt. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
      end if

      if ( cfio%varObjs(i)%timAve ) then
         call writeBnds(cfio, vName, myDate, myCurTime, rtcode)
      end if

      if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOVarWrite3D_

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarWrite1D_ -- Write a variable to a output file
                                                                                
! !INTERFACE:
      subroutine ESMF_CFIOVarWrite1D_(cfio, vName, field, date, curTime,  &
                                      timeString, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj
      character(len=*), intent(in) :: vName       ! Variable name
      real, intent(in) :: field(:)            ! array contains data
      integer, intent(in), OPTIONAL :: date       ! yyyymmdd
      integer, intent(in), OPTIONAL :: curTime    ! hhmmss
      character(len=*), intent(in), OPTIONAL :: timeString
                                  ! string expression for date and time
                                                                                
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !DESCRIPTION:
!     Write a variable to file
!EOP
!------------------------------------------------------------------------------
      integer :: i, rtcode
      integer :: myDate, myCurTime
      character(len=MLEN) :: fNameTmp     ! file name
                                                                                         
      fNameTmp = ''
      if ( present(date) ) myDate = date
      if ( present(curTime) ) myCurTime = curTime
      if ( present(timeString) ) call strToInt(timeString,myDate,myCurTime)
                                                                                         
      if (len(trim(cfio%fNameTmplt)) .gt. 1) then
         call strTemplate_(fNameTmp,cfio%fNameTmplt,xid=cfio%expid,nymd=myDate, &
                           nhms=myCurTime, stat=rtcode)
         if (trim(fNameTmp) .ne. trim(cfio%fName)) then
            call ESMF_CFIOFileClose(cfio)
            cfio%fName = fNameTmp
            call ESMF_CFIOSet(cfio, fName=cfio%fName)
            call ESMF_CFIOSet(cfio, date=myDate, begTime=myCurTime)
            if (len(trim(cfio%expid)) .gt. 0) then
               call ESMF_CFIOFileCreate(cfio, expid=cfio%expid)
            else
               call ESMF_CFIOFileCreate(cfio)
            end if
         end if
      end if
!
!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
                                                                                
!     NEED WORK HERE
      if (index(cfio%varObjs(i)%grid%gName,'station') .gt. 0) then
         call GFIO_SPutVar (cfio%fid, vName, myDate, myCurTime,      &
                  cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,  &
                  0, 1, field, rtcode )
         if (err("GFIO_SPutVar failed",rtcode,rtcode) .lt. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
      else
         if (err("It isn't 1D station grid",rtcode,-1) .lt. 0 ) return
      end if

      if ( cfio%varObjs(i)%timAve ) then
         call writeBnds(cfio, vName, myDate, myCurTime, rtcode)
      end if

      if ( present(rc) ) rc = rtcode
                                                                                
      end subroutine ESMF_CFIOVarWrite1D_
                                                                                
                                                                                
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarWrite2D_ -- Write a variable to a output file
                                                                                
! !INTERFACE:
      subroutine ESMF_CFIOVarWrite2D_(cfio, vName, field, date, curTime, &
                                      kbeg, kount, timeString, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio    ! a CFIO obj
      character(len=*), intent(in) :: vName     ! Variable name
      real, intent(in) :: field(:,:)            ! array contains data
      integer, intent(in), OPTIONAL :: date     ! yyyymmdd
      integer, intent(in), OPTIONAL :: curTime  ! hhmmss
      integer, intent(in), OPTIONAL :: kbeg     ! first level to write
      integer, intent(in), OPTIONAL :: kount    ! number of levels to write
      character(len=*), intent(in), OPTIONAL :: timeString
                                  ! string expression for date and time
                                                                                
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !DESCRIPTION:
!     Write a variable to file
!EOP
!------------------------------------------------------------------------------
      integer :: i, rtcode
      integer :: myKbeg, myKount
      integer :: myDate, myCurTime
      character(len=MLEN) :: fNameTmp     ! file name
                                                                                         
      fNameTmp = ''
      if ( present(date) ) myDate = date
      if ( present(curTime) ) myCurTime = curTime
      if ( present(timeString) ) call strToInt(timeString,myDate,myCurTime)
                                                                                         
      if (len(trim(cfio%fNameTmplt)) .gt. 1) then
         call strTemplate_(fNameTmp,cfio%fNameTmplt,xid=cfio%expid,nymd=myDate, &
                           nhms=myCurTime, stat=rtcode)
         if (trim(fNameTmp) .ne. trim(cfio%fName)) then
            call ESMF_CFIOFileClose(cfio)
            cfio%fName = fNameTmp
            call ESMF_CFIOSet(cfio, fName=cfio%fName)
            call ESMF_CFIOSet(cfio, date=myDate, begTime=myCurTime)
            if (len(trim(cfio%expid)) .gt. 0) then
               call ESMF_CFIOFileCreate(cfio, expid=cfio%expid)
            else
               call ESMF_CFIOFileCreate(cfio)
            end if
         end if
      end if

!
!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
                                                                                
!     write 2D variable
      if (index(cfio%varObjs(i)%grid%gName,'station') .gt. 0) then
         if ( cfio%varObjs(i)%twoDimVar ) then
            call GFIO_SPutVar (cfio%fid, vName, myDate, myCurTime,      &
                     cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,  &
                     0, 1, field, rtcode )
            if (err("GFIO_SPutVar failed",rtcode,rtcode) .lt. 0) then
               if ( present(rc) ) rc = rtcode
               return
            end if
         else
            myKbeg = 1
            myKount = cfio%varObjs(i)%grid%km
            if ( present(kbeg) ) myKbeg = kbeg
            if ( present(kount) ) myKount = kount

            call GFIO_SPutVar (cfio%fid, vName, myDate, myCurTime,          &
                     cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,  &
                     myKbeg, myKount, field, rtcode )
            if (err("GFIO_SPutVar failed",rtcode,rtcode) .lt. 0) then
               if ( present(rc) ) rc = rtcode
               return
            end if
         end if
      else
         call GFIO_PutVar (cfio%fid, vName, myDate, myCurTime,              &
                     cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,  &
                     0, 1, field, rtcode )
         if (err("GFIO_PutVar failed",rtcode,rtcode) .lt. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if

      end if
                                                         
      if ( cfio%varObjs(i)%timAve ) then
         call writeBnds(cfio, vName, myDate, myCurTime, rtcode)
      end if

      if ( present(rc) ) rc = rtcode
                                                                                
      end subroutine ESMF_CFIOVarWrite2D_
                                                                                

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarRead3D_ -- Read a variable from an existing file

! !INTERFACE:
      subroutine ESMF_CFIOVarRead3D_(cfio, vName, field, date, curTime, &
                                     kBeg, kount, xBeg, xCount, yBeg,   &
                                     yCount, timeString, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      integer, intent(in), OPTIONAL :: date       ! yyyymmdd
      integer, intent(in), OPTIONAL :: curTime    ! hhmmss
      integer, intent(in), OPTIONAL :: kbeg       ! first level to write
      integer, intent(in), OPTIONAL :: kount      ! number of levels to write
      integer, intent(in), OPTIONAL :: xBeg       ! first point for lon 
      integer, intent(in), OPTIONAL :: xCount     ! number of points to read
      integer, intent(in), OPTIONAL :: yBeg       ! first point for lat 
      integer, intent(in), OPTIONAL :: yCount     ! number of points to read
      character(len=*), intent(in), OPTIONAL :: timeString
                                  ! string expression for date and time
                                                                                       
!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:,:,:)             ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -46  error from ncvgt
                         !  rc = -48  error from ncinq
                         !  rc = -52  error from ncvinq
!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------
      integer :: i, j, k, rtcode
      integer :: myKbeg, myKount
      integer :: myXbeg, myXount
      integer :: myYbeg, myYount
      integer :: myDate, myCurTime
      real, pointer :: tmp(:,:,:)          ! array contains data
      character(len=MLEN) :: fNameTmp     ! file name
                                                                                         
      fNameTmp = ''
                                                                                         
      if ( present(date) ) myDate = date
      if ( present(curTime) ) myCurTime = curTime
      if ( present(timeString) ) call strToInt(timeString,myDate,myCurTime)

      if (len(trim(cfio%fNameTmplt)) .gt. 1) then
         call strTemplate_(fNameTmp,cfio%fNameTmplt,xid=cfio%expid,nymd=MYdate, &
                           nhms=MYcurTime, stat=rtcode)
         if (trim(fNameTmp) .ne. trim(cfio%fName)) then
            call ESMF_CFIOFileClose(cfio)
            cfio%fName = fNameTmp
            if (len(trim(cfio%expid)) .gt. 0) then
               call ESMF_CFIOFileOpen(cfio, 1, expid=cfio%expid, cyclic=cfio%isCyclic)
            else
               call ESMF_CFIOFileOpen(cfio, 1, cyclic=cfio%isCyclic)
            end if
         end if
      end if

!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do

      myKbeg = 1
      myKount = 1

!     read 3D variable
      if ( cfio%varObjs(i)%grid%km .gt. 1 .and.                          &
           (.not. cfio%varObjs(i)%twoDimVar) ) then

         myKbeg = 1
         myKount = cfio%varObjs(i)%grid%km
         if ( present(kbeg) ) myKbeg = kbeg
         if ( present(kount) ) myKount = kount

         allocate(tmp(cfio%varObjs(i)%grid%im,cfio%varObjs(i)%grid%jm,   &
               myKount), stat=rtcode)
         if (rtcode /= 0) print *, "cannot allocate tmp in ESMF_CFIOVarRead3D"

         call GFIO_GetVar(cfio%fid,vName,mydate,mycurTime,                   &
                       cfio%varObjs(i)%grid%im,                          &
                       cfio%varObjs(i)%grid%jm,myKbeg,myKount,           &
                       cfio%tSteps, tmp, cfio%isCyclic, rtcode )
         if (rtcode .ne. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
!     read 2D variable
      else
         allocate(tmp(cfio%varObjs(i)%grid%im,cfio%varObjs(i)%grid%jm,1),&
                  stat=rtcode)
         if (rtcode /= 0) print *, "cannot allocate tmp in ESMF_CFIOVarRead3D"
         
         call GFIO_GetVar(cfio%fid,vName,mydate,MYcurTime,                   &
                       cfio%varObjs(i)%grid%im,                          &
                       cfio%varObjs(i)%grid%jm, 0, 1, cfio%tSteps, tmp,  &
                       cfio%isCyclic, rtcode )
         if (rtcode .ne. 0) then
            if ( present(rc) ) rc = rtcode
            return
         end if
      end if

      myXbeg = 1
      myXount = cfio%varObjs(i)%grid%im
      myYbeg = 1
      myYount = cfio%varObjs(i)%grid%jm
      if ( present(xBeg) ) myXbeg=xBeg
      if ( present(yBeg) ) myYbeg=yBeg
      if ( present(xCount) ) myXount = xCount
      if ( present(yCount) ) myYount = yCount

      if (.not. associated(field) ) then
         allocate(field(myXount,myYount,myKount),stat=rtcode)
      else
         deallocate(field,stat=rtcode)
         if (rtcode /= 0) print *, "Couldn't deallocate Field in VarRead3D"
         allocate(field(myXount,myYount,myKount),stat=rtcode)
      end if
!      allocate(field(myXount,myYount,myKount), stat=rtcode)
      if (rtcode /= 0) print *, "cannot allocate field in ESMF_CFIOVarRead3D_"
      do k = 1, myKount
         do j = 1, myYount
           do i = 1, myXount
              field(i,j,k) = tmp(myXbeg+i-1,myYbeg+j-1,k)
           end do
         end do
      end do

      deallocate(tmp)
      if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOVarRead3D_ 

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarRead2D_ -- Read a variable from an existing file
                                                                                
! !INTERFACE:
      subroutine ESMF_CFIOVarRead2D_(cfio, vName, field, date, curTime, &
                                     kbeg, kount, xBeg, xCount, yBeg,   &
                                     yCount, timeString, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inout) :: cfio         ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      integer, intent(in), OPTIONAL :: date       ! yyyymmdd
      integer, intent(in), OPTIONAL :: curTime    ! hhmmss
      integer, intent(in), OPTIONAL :: kbeg       ! first level to write
      integer, intent(in), OPTIONAL :: kount      ! number of levels to write
      integer, intent(in), OPTIONAL :: xBeg       ! first point for lon
      integer, intent(in), OPTIONAL :: xCount     ! number of points to read
      integer, intent(in), OPTIONAL :: yBeg       ! first point for lat
      integer, intent(in), OPTIONAL :: yCount     ! number of points to read
      character(len=*), intent(in), OPTIONAL :: timeString
                                  ! string expression for date and time

!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:,:)             ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -46  error from ncvgt
                         !  rc = -48  error from ncinq
                         !  rc = -52  error from ncvinq

!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------
      integer :: i, j, k, rtcode
      integer :: myKbeg, myKount
      integer :: myXbeg, myXount
      integer :: myYbeg, myYount
      integer :: myDate, myCurTime
      real, pointer :: tmp(:,:)          ! array contains data
      character(len=MLEN) :: fNameTmp     ! file name
                                                                                              
      fNameTmp = ''

      if ( present(date) ) myDate = date
      if ( present(curTime) ) myCurTime = curTime
      if ( present(timeString) ) call strToInt(timeString,myDate,myCurTime)
                                                                                              
      if (len(trim(cfio%fNameTmplt)) .gt. 1) then
         call strTemplate_(fNameTmp,cfio%fNameTmplt,xid=cfio%expid,nymd=MYdate, &
                           nhms=MYcurTime, stat=rtcode)
         if (trim(fNameTmp) .ne. trim(cfio%fName)) then
            call ESMF_CFIOFileClose(cfio)
            cfio%fName = fNameTmp
!            call ESMF_CFIOSet(cfio, fName=cfio%fName)
            if (len(trim(cfio%expid)) .gt. 0) then
               call ESMF_CFIOFileOpen(cfio, 1, expid=cfio%expid, cyclic=cfio%isCyclic)
            else
               call ESMF_CFIOFileOpen(cfio, 1, cyclic=cfio%isCyclic)
            end if
         end if
      end if
                                                                                
!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
                                                                                
      myXbeg = 1
      myXount = cfio%varObjs(i)%grid%im
      myYbeg = 1
      myYount = cfio%varObjs(i)%grid%jm
      myKbeg = 1
      myKount = cfio%varObjs(i)%grid%km
      if ( present(xBeg) ) myXbeg=xBeg
      if ( present(yBeg) ) myYbeg=yBeg
      if ( present(kbeg) ) myKbeg = kbeg
      if ( present(kount) ) myKount = kount
      if ( present(xCount) ) myXount = xCount
      if ( present(yCount) ) myYount = yCount

!     read 2D variable
      if ( cfio%varObjs(i)%twoDimVar .and.                              &
              .not. cfio%varObjs(i)%grid%stnGrid) then
        allocate(tmp(cfio%varObjs(i)%grid%im,cfio%varObjs(i)%grid%jm),  &
               stat=rtcode)
        call GFIO_GetVar(cfio%fid,vName,MYdate,MYcurTime,                   &
                    cfio%varObjs(i)%grid%im,                            &
                    cfio%varObjs(i)%grid%jm, 0, 1, cfio%tSteps, tmp,    &
                    cfio%isCyclic, rtcode )
        if (err("GFIO_GetVar failed",rtcode,rtcode) .lt. 0) then  
           if ( present(rc) ) rc = rtcode
           return
        end if
 
        allocate(field(myXount,myYount))
        do j = 1, myYount
           do i = 1, myXount
              field(i,j) = tmp(myXbeg+i-1,myYbeg+j-1)
           end do
        end do

      else
        if (cfio%varObjs(i)%twoDimVar ) then
           allocate(tmp(cfio%varObjs(i)%grid%im,1), stat=rtcode)
           call GFIO_SGetVar(cfio%fid,vName,MYdate,MYcurTime,               &
                    cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,   &
                    0,1, cfio%tSteps, tmp, cfio%isCyclic, rtcode )
           if (err("GFIO_SGetVar failed",rtcode,rtcode) .lt. 0) then  
              if ( present(rc) ) rc = rtcode
              return
           end if
           allocate(field(myXount,1))
           do i = 1, myXount
              field(i,1) = tmp(myXbeg+i-1,1)
           end do

        else
           allocate(tmp(cfio%varObjs(i)%grid%im,myKount),stat=rtcode)
           call GFIO_SGetVar(cfio%fid,vName,MYdate,MYcurTime,               &
                    cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,   &
                    myKbeg, myKount, cfio%tSteps, tmp, cfio%isCyclic, rtcode )
           if (err("GFIO_GetVar failed",rtcode,rtcode) .lt. 0) then  
              if ( present(rc) ) rc = rtcode
              return
           end if
           allocate(field(myXount,myKount))
           do k = 1, myKount
             do i = 1, myXount
                field(i,k) = tmp(myXbeg+i-1,k)
             end do
           end do

        end if
      end if
 
      deallocate(tmp)
                                                                                
      if ( present(rc) ) rc = rtcode
                                                                                
      end subroutine ESMF_CFIOVarRead2D_
                                                                                
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarRead1D_ -- Read a variable from an existing file
                                                                                
! !INTERFACE:
      subroutine ESMF_CFIOVarRead1D_(cfio, vName, field, date, curTime, &
                                     xBeg, xCount, timestring, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio         ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      integer, intent(in), OPTIONAL :: date       ! yyyymmdd
      integer, intent(in), OPTIONAL :: curTime    ! hhmmss
      integer, intent(in), OPTIONAL :: xBeg       ! first point for lon
      integer, intent(in), OPTIONAL :: xCount     ! number of points to read
      character(len=*), intent(in), OPTIONAL :: timeString
                                  ! string expression for date and time
!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:)             ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------

      integer :: i, j, rtcode
      integer :: myXbeg, myXount      
      integer :: myDate, myCurTime
      real, pointer :: tmp(:)          ! array contains data
      character(len=MLEN) :: fNameTmp     ! file name
                                                                                              
      fNameTmp = ''
      if ( present(date) ) myDate = date
      if ( present(curTime) ) myCurTime = curTime
      if ( present(timeString) ) call strToInt(timeString,myDate,myCurTime)
                                                                                              
      if (len(trim(cfio%fNameTmplt)) .gt. 1) then
         call strTemplate_(fNameTmp,cfio%fNameTmplt,xid=cfio%expid,nymd=MYdate, &
                           nhms=MYcurTime, stat=rtcode)
         if (trim(fNameTmp) .ne. trim(cfio%fName)) then
            call ESMF_CFIOFileClose(cfio)
            cfio%fName = fNameTmp
!            call ESMF_CFIOSet(cfio, fName=cfio%fName)
            if (len(trim(cfio%expid)) .gt. 0) then
               call ESMF_CFIOFileOpen(cfio, 1, expid=cfio%expid, cyclic=cfio%isCyclic)
            else
               call ESMF_CFIOFileOpen(cfio, 1, cyclic=cfio%isCyclic)
            end if
         end if
      end if

                                                                                
!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
                                                                                
      myXbeg = 1
      myXount = cfio%varObjs(i)%grid%im

      if (present(xBeg)) myXbeg = xBeg
      if (present(xCount)) myXount = xCount

!     read 1D variable
      allocate(tmp(cfio%varObjs(i)%grid%im), stat=rtcode)
      call GFIO_SGetVar(cfio%fid,vName,MYdate,MYcurTime,               &
               cfio%varObjs(i)%grid%im, cfio%varObjs(i)%grid%jm,   &
               0,1, cfio%tSteps, tmp, cfio%isCyclic, rtcode )

      do i = 1, myXount
         field(i) = tmp(myXbeg+i-1)
      end do

      deallocate(tmp)
                                                                                
      if ( present(rc) ) rc = rtcode
                                                                                
      end subroutine ESMF_CFIOVarRead1D_
                                                                                

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOFileClose -- close an open CFIO stream

! !INTERFACE:
      subroutine ESMF_CFIOFileClose (cfio, rc)
!
! !ARGUMENTS:
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                       ! 0   all is well
                                       ! -54  error from ncclos (file close)
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inout) :: cfio       ! CFIO object


!
! !DESCRIPTION:
!     close an open CFIO stream
!EOP
!------------------------------------------------------------------------------
       integer :: rtcode

       if ( cfio%isOpen ) then 
          call GFIO_Close(cfio%fid, rtcode)
          if (rtcode .ne. 0) then 
             print *, "GFIO_Close failed"
          else
             cfio%isOpen = .false.
          end if
       else
          rtcode = 0
       end if

       if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOFileClose


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIODestroy -- destructor for a CFIO object 

! !INTERFACE:
      subroutine ESMF_CFIODestroy (cfio, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inout) :: cfio       ! CFIO object

!
! !DESCRIPTION:
!     destructor for a CFIO object 
!EOP
!------------------------------------------------------------------------------
      integer :: rtcode

      if ( cfio%isOpen ) call GFIO_Close(cfio%fid, rtcode)
      if ( associated(cfio%varObjs) ) deallocate(cfio%varObjs, stat=rtcode)
      if ( associated(cfio%grids) ) deallocate(cfio%grids, stat=rtcode)

      if (associated(cfio%attCharCnts)) deallocate(cfio%attCharCnts,         &
         stat=rtcode)
      if (associated(cfio%attRealCnts)) deallocate(cfio%attRealCnts,         &
         stat=rtcode)
      if (associated(cfio%attIntCnts)) deallocate(cfio%attIntCnts,           &
         stat=rtcode)

      if (associated(cfio%attCharNames)) deallocate(cfio%attCharNames,       &
         stat=rtcode)
      if (associated(cfio%attRealNames)) deallocate(cfio%attRealNames,       &
         stat=rtcode)
      if (associated(cfio%attIntNames)) deallocate(cfio%attIntNames,         &
         stat=rtcode)

      if (associated(cfio%attChars)) deallocate(cfio%attChars, stat=rtcode)
      if (associated(cfio%attReals)) deallocate(cfio%attReals, stat=rtcode)
      if (associated(cfio%attInts)) deallocate(cfio%attInts, stat=rtcode)

      if (associated(cfio%iList)) deallocate(cfio%iList, stat=rtcode)
      if (associated(cfio%rList)) deallocate(cfio%rList, stat=rtcode)
      if (associated(cfio%cList)) deallocate(cfio%cList, stat=rtcode)

      rtcode = 0
      if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIODestroy


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarInfoCreate -- ESMF_CFIOVarInfo object constructor

! !INTERFACE:
      type(ESMF_CFIOVarInfo) function ESMF_CFIOVarInfoCreate (vName, rc)   
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in), OPTIONAL :: vName     ! variable name
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                                                ! -1  problem in creating Grid
!
! !DESCRIPTION:
!     Create a CFIO varInfo object and initialize variables
!EOP
!------------------------------------------------------------------------------
      type(ESMF_CFIOVarInfo) :: varObj ! a CFIO grid object
      integer :: rtcode = 0
      
      varObj%grid = ESMF_CFIOGridCreate(rc=rtcode)
      if (rtcode .ne. 0) then 
         print *, "problem in getting ESMF_CFIOGridCreate:lon"
         rtcode = -1
         if ( present(rc) ) rc = rtcode
         return
      end if

      varObj%nVarAttInt = 0
      varObj%nVarAttChar = 0
      varObj%nVarAttReal = 0

      varObj%twoDimVar = .false.
      varObj%timAve = .false.
      varObj%aveMethod = 'c'
      varObj%cellMthd = 'mean'
      varObj%amiss = 1.E15
      varObj%addOffSet = 0
      varObj%scaleFactor = 1
      varObj%validRange = 1.E15
      varObj%packingRange = 1.E15
      varObj%ordering = 'tzyx'

      varObj%vTitle = 'unknow'
      varObj%vUnits = 'unknow'
      varObj%standardName = 'unknow'

      allocate(varObj%iList, varObj%rList, varObj%cList)
      nullify(varObj%iList)
      nullify(varObj%rList)
      nullify(varObj%cList)

      if ( present(rc) ) rc = 0

      ESMF_CFIOVarInfoCreate = varObj

      end function ESMF_CFIOVarInfoCreate



!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarInfoSet  -- Set meta data for a CFIO variable

! !INTERFACE:
      subroutine ESMF_CFIOVarInfoSet (varObj, vName, grid, vTitle, vUnits,  &
                              twoDimVar, validRange, amiss, addOffSet,      &
                              scaleFactor, standardName, attCharNames,      &
                              vAttCharCnts, varAttChars, attRealNames,      &
                              vAttRealCnts, varAttReals, attIntNames,       &
                              vAttIntCnts, varAttInts, ordering,            &
                              attCharName, attChar, attRealName, attReal,   &
                              attIntName, attInt, packingRange, timAve,     &
                              aveMethod, cellMthd, rc )
       implicit NONE

! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
       character(len=*), intent(in), OPTIONAL :: vName  ! variable name 
       type(ESMF_CFIOGrid), intent(in), OPTIONAL :: grid  ! grid 
       character(len=*), intent(in), OPTIONAL :: vTitle ! variable long name 
       character(len=*), intent(in), OPTIONAL :: vUnits ! variable units      
       logical, OPTIONAL :: twoDimVar                 ! True for 2D
       logical, OPTIONAL :: timAve                    ! True for time average
       character, OPTIONAL :: aveMethod  ! 'c': center, 'd': downstream
                                         ! 'u': upstream
       character(len=*), OPTIONAL :: cellMthd    ! Cell methmod units
       real, intent(in), OPTIONAL :: validRange(2)    ! Variable valid range
       real, intent(in), OPTIONAL :: packingRange(2)  ! Variable valid range
       real, intent(in), OPTIONAL :: amiss            ! FILL value
       real, intent(in), OPTIONAL :: addOffSet    
       real, intent(in), OPTIONAL :: scaleFactor
       character(len=*), intent(in), OPTIONAL :: standardName 

       character(len=*), intent(in), OPTIONAL :: attCharNames(:)
       character(len=*), intent(in), OPTIONAL :: attRealNames(:)
       character(len=*), intent(in), OPTIONAL :: attIntNames(:)

       integer, intent(in), OPTIONAL :: vAttCharCnts(:) ! length of attributes
       integer, intent(in), OPTIONAL :: vAttRealCnts(:) ! length of attributes
       integer, intent(in), OPTIONAL :: vAttIntCnts(:)  ! length of attributes

       character(len=*), intent(in), OPTIONAL :: varAttChars(:)
       real, intent(in), OPTIONAL :: varAttReals(:,:)
       integer, intent(in), OPTIONAL :: varAttInts(:,:)
       
       character(len=*), OPTIONAL :: ordering 
                                 ! (time, lev, lat, lon) (default)
                                 ! can be any combination of xyzt

       character(len=*), intent(in), OPTIONAL :: attCharName
                                 ! User defined variable attribute name
       character(len=*), intent(in), OPTIONAL :: attRealName
                                 ! User defined variable real attribute name
       character(len=*), intent(in), OPTIONAL :: attIntName
                                 ! User defined variable int attribute name
       character(len=*), intent(in), OPTIONAL :: attChar
                                 ! User defined variable char attribute
       real,    intent(in), OPTIONAL :: attReal(:)
                                 ! User defined variable real attribute
       integer, intent(in), OPTIONAL :: attInt(:)
                                 ! User defined variable int attribute
!
! !OUTPUT PARAMETERS:
!
       integer, intent(out), OPTIONAL :: rc   
                                 ! 0  all is well
                                 ! -1 Allocation for attCharCnts failed
                                 ! -2 Allocation for attRealCnts failed
                                 ! -3 Allocation for attIntCnts failed
                                 ! -4 Allocation for varAttChars failed
                                 ! -5 Allocation for varAttReals failed
                                 ! -6 Allocation for varAttInts failed
                                 ! -7 Allocation for attCharNames failed
                                 ! -8 Allocation for attRealNames failed
                                 ! -9 Allocation for attIntNames failed
                                               

! !INPUT/OUTPUT PARAMETERS:
!
       type(ESMF_CFIOVarInfo), intent(inout) :: varObj        ! variable obj 
 
!
! !DESCRIPTION:
!     Set meta data for a CFIO variable         
!EOP
!------------------------------------------------------------------------------
       integer :: iCnt, jCnt, count, rtcode = 0
                                                                  
       if ( present(vName) ) varObj%vName = vName
       if ( present(grid) ) varObj%grid = grid 
       if ( present(vTitle) ) varObj%vTitle = vTitle
       if ( present(vUnits) ) varObj%vUnits = vUnits
       if ( present(twoDimVar) ) varObj%twoDimVar = twoDimVar
       if ( present(timAve) ) varObj%timAve = timAve
       if ( present(aveMethod) ) varObj%aveMethod = aveMethod
       if ( present(cellMthd) ) varObj%cellMthd = trim(cellMthd) 
       if ( present(validRange) ) varObj%validRange = validRange
       if ( present(packingRange) ) varObj%packingRange = packingRange
       if ( present(amiss) ) varObj%amiss = amiss
       if ( present(addOffSet) ) varObj%addOffSet = addOffSet
       if ( present(scaleFactor) ) varObj%scaleFactor = scaleFactor
       if ( present(standardName) )  varObj%standardName = standardName
       if ( present(ordering) ) varObj%ordering = ordering 
       
!      user provide int/real/char attribute counts as arrays 
       if ( present(vAttCharCnts) ) then
          allocate(varObj%attCharCnts(size(vAttCharCnts)), stat=rtcode)
          if (err("Allocation for attCharCnts failed",rtcode,-1) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if

          varObj%attCharCnts = vAttCharCnts
       end if
       if ( present(vAttRealCnts) ) then
          allocate(varObj%attRealCnts(size(vAttRealCnts)), stat=rtcode)
          if (err("Allocation for attRealCnts failed",rtcode,-2) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%attRealCnts = vAttRealCnts
          varObj%nVarAttReal = size(vAttRealCnts)
       end if
       if ( present(vAttIntCnts) ) then
          allocate(varObj%attIntCnts(size(vAttIntCnts)), stat=rtcode)
          if (err("Allocation for attIntCnts failed",rtcode,-3) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%attIntCnts = vAttIntCnts
          varObj%nVarAttInt = size(vAttIntCnts)
       end if

!      user provide int/real/char attribute data as arrays 
       if ( present(varAttChars) ) then
          allocate(varObj%varAttChars(size(varAttChars)), stat=rtcode)
          if (err("Allocation for varAttChars failed",rtcode,-4) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%varAttChars = varAttChars
          varObj%nVarAttChar = size(vAttCharCnts)
       end if
       if ( present(varAttReals) ) then
          iCnt = size(varObj%attRealCnts)
          jCnt = size(varAttReals)/size(varObj%attRealCnts)
          allocate(varObj%varAttReals(iCnt, jCnt), stat=rtcode)
          if (err("Allocation for varAttReals failed",rtcode,-5) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%varAttReals = varAttReals
       end if
       if ( present(varAttInts) ) then
          iCnt = size(varObj%attIntCnts)
          jCnt = size(varAttInts)/size(varObj%attIntCnts)
          allocate(varObj%varAttInts(iCnt, jCnt), stat=rtcode)
          if (err("Allocation for varAttInts failed",rtcode,-6) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%varAttInts = varAttInts  
       end if
    
!      user provide int/real/char attribute names as arrays 
       if ( present(attCharNames)) then
          allocate(varObj%attCharNames(size(attCharNames)), stat=rtcode)
          if (err("Allocation for attCharNames failed",rtcode,-7) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%attCharNames = attCharNames
       end if
       if ( present(attRealNames)) then
          allocate(varObj%attRealNames(size(attRealNames)), stat=rtcode)
          if (err("Allocation for attRealNames failed",rtcode,-8) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%attRealNames = attRealNames
       end if
       if ( present(attIntNames)) then
          allocate(varObj%attIntNames(size(attIntNames)), stat=rtcode)
          if (err("Allocation for attIntNames failed",rtcode,-9) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varObj%attIntNames = attIntNames
       end if

!      user provides real attribute name and data. Put them into rList
       if ( present(attRealName) .and. present(attReal) ) then
          count = size(attReal)
          call addList(attRealName, count, attReal=attReal, &
                       rList=varObj%rList)
       end if

!      user provides int attribute name and data. Put them into iList
       if ( present(attIntName) .and. present(attInt) ) then
          count = size(attInt)
          call addList(attIntName, count, attInt=attInt, &
                       iList=varObj%iList)
       end if

!      user provides char attribute name and data. Put them into cList
       if ( present(attCharName) .and. present(attChar) ) then
          call addList(attCharName, len(attChar), attChar=attChar, &
                       cList=varObj%cList)
       end if

       if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOVarInfoSet 


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarInfoGet -- get information from a CFIO variable object

! !INTERFACE:
      subroutine ESMF_CFIOVarInfoGet (varObj, vName, grid, vTitle, vUnits,  &
                              twoDimVar, validRange, amiss, addOffSet,      &
                              scaleFactor, standardName, nVarAttChar,       &
                              attCharNames,  vAttCharCnts, varAttChars,     &
                              nVarAttReal, attRealNames, vAttRealCnts,      &
                              varAttReals, nVarAttInt, attIntNames,         &
                              vAttIntCnts, varAttInts, ordering,            &
                              attCharName, attCharCnt, attChar, attRealName,&
                              attRealCnt, attReal, attIntName, attIntCnt,   &
                              attInt, packingRange, timAve, aveMethod,      &
                              cellMthd, rc )
       implicit NONE

! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
       type(ESMF_CFIOVarInfo), intent(in) :: varObj      ! variable obj 
       character(len=*), intent(in), OPTIONAL :: attCharName
                                    ! User defined  char attribute name
       character(len=*), intent(in), OPTIONAL :: attRealName
                                    ! User defined  real attribute name
       character(len=*), intent(in), OPTIONAL :: attIntName
                                    ! User defined  int attribute name


! !OUTPUT PARAMETERS:
!
       character(len=*), intent(out), OPTIONAL :: vName  ! variable short name
       type(ESMF_CFIOGrid), intent(out), OPTIONAL :: grid  ! grid 
       character(len=*), intent(out), OPTIONAL :: vTitle ! variable long name 
       character(len=*), intent(out), OPTIONAL :: vUnits ! variable units      
       logical, OPTIONAL :: twoDimVar                  ! True for 2D
       logical, OPTIONAL :: timAve                    ! True for time average
       character, OPTIONAL :: aveMethod  ! 'c': center, 'd': downstream
                                         ! 'u': upstream
       character(len=MVARLEN), OPTIONAL :: cellMthd    ! Cell methmod units
       real, intent(out), OPTIONAL :: validRange(2)    ! Variable valid range
       real, intent(out), OPTIONAL :: packingRange(2)
       real, intent(out), OPTIONAL :: amiss            ! FILL value
       real, intent(out), OPTIONAL :: addOffSet    
       real, intent(out), OPTIONAL :: scaleFactor
       character(len=*), intent(out), OPTIONAL :: standardName 
       character(len=*), OPTIONAL :: ordering
                                 ! (time, lev, lat, lon) (default)
                                 ! can be any combination of xyzt
       integer, intent(out), OPTIONAL :: nVarAttInt
       integer, intent(out), OPTIONAL :: nVarAttReal
       integer, intent(out), OPTIONAL :: nVarAttChar

       character(len=*), pointer, OPTIONAL :: attCharNames(:)
       character(len=*), pointer, OPTIONAL :: attRealNames(:)
       character(len=*), pointer, OPTIONAL :: attIntNames(:)

       integer, pointer,OPTIONAL::vAttCharCnts(:) ! length of attributes
       integer, pointer,OPTIONAL::vAttRealCnts(:) ! length of attributes
       integer, pointer,OPTIONAL::vAttIntCnts(:)  ! length of attributes

       character(len=*), pointer, OPTIONAL :: varAttChars(:)
       real, pointer, OPTIONAL :: varAttReals(:,:)
       integer, pointer, OPTIONAL :: varAttInts(:,:)
                                                                                       
       integer, intent(out), OPTIONAL :: attIntCnt
       integer, intent(out), OPTIONAL :: attRealCnt
       integer, intent(out), OPTIONAL :: attCharCnt
       character(len=*), intent(out), OPTIONAL :: attChar
                                    ! User defined  char attribute
       real,    pointer, OPTIONAL :: attReal(:)
                                    ! User defined  real attribute
       integer, pointer, OPTIONAL :: attInt(:)
                                    ! User defined  int attribute

       integer, intent(out), OPTIONAL :: rc 
                                    ! Error return code:
                                    ! 0   all is well
                                    ! -1  Allocation for attCharNames failed
                                    ! -2  Allocation for attRealNames failed
                                    ! -3  Allocation for attIntNames failed 
                                    ! -4  Allocation for vAttCharCnts failed
                                    ! -5  Allocation for vAttRealCnts failed 
                                    ! -6  Allocation for vAttIntCnts failed 
                                    ! -7  Allocation for varAttChars failed 
                                    ! -8  Allocation for varAttReals failed 
                                    ! -9  Allocation for varAttInts failed  

!
! !DESCRIPTION:
!     get information from a CFIO variable object
!EOP
!------------------------------------------------------------------------------
       integer :: rtcode = 0
       integer :: i
                                                                                     
       if ( present(vName) ) vName = varObj%vName
       if ( present(grid) ) grid = varObj%grid 
       if ( present(vTitle) ) vTitle = varObj%vTitle 
       if ( present(vUnits) ) vUnits = varObj%vUnits
       if ( present(twoDimVar) ) twoDimVar = varObj%twoDimVar 
       if ( present(timAve) ) timAve = varObj%timAve
       if ( present(aveMethod) ) aveMethod = varObj%aveMethod
       if ( present(cellMthd) ) cellMthd = varObj%cellMthd  
       if ( present(validRange) ) validRange = varObj%validRange 
       if ( present(packingRange) ) packingRange = varObj%packingRange
       if ( present(amiss) ) amiss = varObj%amiss 
       if ( present(addOffSet) ) addOffSet = varObj%addOffSet
       if ( present(scaleFactor) ) scaleFactor = varObj%scaleFactor 
       if ( present(standardName) )  standardName = varObj%standardName 
       if ( present(ordering) ) ordering = varObj%ordering 
                                                                                     
       if ( present(nVarAttInt) ) nVarAttInt = varObj%nVarAttInt
       if ( present(nVarAttReal) ) nVarAttReal = varObj%nVarAttReal
       if ( present(nVarAttChar) ) nVarAttChar = varObj%nVarAttChar

!      get all attribute names
       if ( present(attCharNames) ) then
          allocate(attCharNames(varObj%nVarAttChar), stat=rtcode)
          if (err("Allocation for attCharNames failed",rtcode,-1) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          attCharNames = varObj%attCharNames
       end if
       if ( present(attRealNames) ) then
          allocate(attRealNames(varObj%nVarAttReal), stat=rtcode)
          if (err("Allocation for attRealNames failed",rtcode,-2) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          attRealNames = varObj%attRealNames
       end if
       if ( present(attIntNames) ) then
          allocate(attIntNames(varObj%nVarAttInt), stat=rtcode)
          if (err("Allocation for attIntNames failed",rtcode,-3) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          attIntNames = varObj%attIntNames
       end if

!      get all attribute counts
       if ( present(vAttCharCnts) ) then
          allocate(vAttCharCnts(varObj%nVarAttChar), stat=rtcode)
          if (err("Allocation for vAttCharCnts failed",rtcode,-4) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          vAttCharCnts = varObj%attCharCnts
       end if
       if ( present(vAttRealCnts) ) then
          allocate(vAttRealCnts(varObj%nVarAttReal), stat=rtcode)
          if (err("Allocation for vAttRealCnts failed",rtcode,-5) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          vAttRealCnts = varObj%attRealCnts
       end if
       if ( present(vAttIntCnts) ) then
          allocate(vAttIntCnts(varObj%nVarAttInt), stat=rtcode)
          if (err("Allocation for vAttIntCnts failed",rtcode,-6) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          vAttIntCnts = varObj%attIntCnts
       end if

!      get all attribute data   
       if ( present(varAttChars) ) then 
          allocate(varAttChars(varObj%nVarAttChar), stat=rtcode)
          if (err("Allocation for varAttChars failed",rtcode,-7) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varAttChars = varObj%varAttChars 
       end if
       if ( present(varAttReals) ) then
          allocate(varAttReals(varObj%nVarAttReal,size(varObj%varAttReals) &
                   / varObj%nVarAttReal), stat=rtcode)
          if (err("Allocation for varAttReals failed",rtcode,-8) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varAttReals = varObj%varAttReals 
       end if
       if ( present(varAttInts) ) then 
          allocate(varAttInts(varObj%nVarAttInt,size(varObj%varAttInts) &
                   / varObj%nVarAttInt), stat=rtcode)
          if (err("Allocation for varAttInts failed",rtcode,-9) .lt. 0) then
             if ( present(rc) ) rc = rtcode
             return
          end if
          varAttInts = varObj%varAttInts  
       end if

!      user provides integer attribute name to get its count and data 
       if ( present(attIntName) ) then
          if ( present(attIntCnt) ) then
             do i = 1, varObj%nVarAttInt
                if (trim(attIntName) .eq. trim(varObj%attIntNames(i))) &
                    then
                   attIntCnt = varObj%attIntCnts(i)
                end if
             end do
          end if
          if ( present(attInt) ) then
             do i = 1, varObj%nVarAttInt
                if (trim(attIntName) .eq. trim(varObj%attIntNames(i)))&
                    then
                   allocate(attInt(varObj%attIntCnts(i)))
                   attInt = varObj%varAttInts(i,1:varObj%attIntCnts(i))
                end if
             end do
          end if
       end if

!      user provides real attribute name to get its count and data 
       if ( present(attRealName) ) then
          if ( present(attRealCnt) ) then
             do i = 1, varObj%nVarAttReal
                if (trim(attRealName) .eq. trim(varObj%attRealNames(i))) &
                    then
                   attRealCnt = varObj%attRealCnts(i)
                end if
             end do
          end if
          if ( present(attReal) ) then
             do i = 1, varObj%nVarAttReal
                if (trim(attRealName) .eq. trim(varObj%attRealNames(i)))&
                    then
                   allocate(attReal(varObj%attRealCnts(i)))
                   attReal = varObj%varAttReals(i,1:varObj%attRealCnts(i))
                end if
             end do
          end if
       end if

!      user provides char attribute name to get its count and data 
       if ( present(attCharName) ) then
          if ( present(attCharCnt) ) then
             do i = 1, varObj%nVarAttChar
                if (trim(attCharName) .eq. trim(varObj%attCharNames(i))) &
                    then
                   attCharCnt = varObj%attCharCnts(i)
                end if
             end do
          end if
          if ( present(attChar) ) then
             do i = 1, varObj%nVarAttChar
                if (trim(attCharName) .eq. trim(varObj%attCharNames(i)))&
                    then
                   attChar = trim(varObj%varAttChars(i))
                end if
             end do
          end if
       end if

       if ( present(rc) ) rc = rtcode


      end subroutine ESMF_CFIOVarInfoGet

!------------------------------------------------------------------------------!BOP
! !ROUTINE: ESMF_CFIOVarInfoDestroy -- destructor for a CFIO varInfo object

! !INTERFACE:
      subroutine ESMF_CFIOVarInfoDestroy (varObj, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_CFIOVarInfo), intent(inout) :: varObj ! CFIOVarInfo object

!
! !DESCRIPTION:
!     destructor for a CFIO varInfo object
!EOP
!------------------------------------------------------------------------------
      integer :: rtcode = 0

      if ( associated(varObj%attCharCnts) ) deallocate(varObj%attCharCnts, &
                                             stat=rtcode)
      if ( associated(varObj%attRealCnts) ) deallocate(varObj%attRealCnts, &
                                             stat=rtcode)
      if ( associated(varObj%attIntCnts) ) deallocate(varObj%attIntCnts,   &
                                             stat=rtcode)

      if ( associated(varObj%attCharNames) ) deallocate(varObj%attCharNames,&
                                             stat=rtcode)
      if ( associated(varObj%attRealNames) ) deallocate(varObj%attRealNames,&
                                             stat=rtcode)
      if ( associated(varObj%attIntNames) ) deallocate(varObj%attIntNames,  &
                                             stat=rtcode)

      if ( associated(varObj%varAttChars) ) deallocate(varObj%varAttChars,  &
                                             stat=rtcode)
      if ( associated(varObj%varAttReals) ) deallocate(varObj%varAttReals,  &
                                             stat=rtcode)
      if ( associated(varObj%varAttInts) ) deallocate(varObj%varAttInts,    &
                                             stat=rtcode)

      if ( associated(varObj%grid%lon) ) deallocate(varObj%grid%lon)
      if ( associated(varObj%grid%lat) ) deallocate(varObj%grid%lat)
      if ( associated(varObj%grid%lev) ) deallocate(varObj%grid%lev)
      if ( associated(varObj%grid%ak) ) deallocate(varObj%grid%ak)
      if ( associated(varObj%grid%bk) ) deallocate(varObj%grid%bk)
      if ( associated(varObj%grid%sigma) ) deallocate(varObj%grid%sigma)

      if ( associated(varObj%iList) ) deallocate(varObj%iList, stat=rtcode)
      if ( associated(varObj%rList) ) deallocate(varObj%rList, stat=rtcode)
      if ( associated(varObj%cList) ) deallocate(varObj%cList, stat=rtcode)

      if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOVarInfoDestroy

!------------------------------------------------------------------------------!BOP
! !ROUTINE: ESMF_CFIOGridCreate -- ESMF_Grid object constructor

! !INTERFACE:
      type(ESMF_CFIOGrid) function ESMF_CFIOGridCreate (gName, rc)   
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in), OPTIONAL :: gName  ! grid object name
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !DESCRIPTION:
!     Create a CFIO grid object and initialize vars
!EOP
!------------------------------------------------------------------------------
      type(ESMF_CFIOGrid) :: grid               ! a CFIO grid object

      grid%im = 0
      grid%jm = 0
      grid%km = 0

      grid%levUnits = 'unknown'
      grid%coordinate = 'unknown'
      grid%standardName = 'unknown'
      grid%formulaTerm = 'unknown'

      grid%ptop = 0
      grid%ptopUnit = 'Pa' 
      grid%twoDimLat = .false.
      grid%reduceGrid = .false.
      grid%stnGrid = .false.

      if ( present(gName) ) grid%gName = gName

      if ( present(rc) ) rc = 0

      ESMF_CFIOGridCreate = grid

      end function ESMF_CFIOGridCreate



!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOGridSet -- set up a grid

! !INTERFACE:
      subroutine ESMF_CFIOGridSet (grid, gName, im, jm, km, lat, lon, lev,  &
                                   coordinate, standardName, formulaTerm,   &
                                   levUnit, ak, bk, sigma, ptop, ptopUnit,  &
                                   twoDimLat, reduceGrid, stnGrid, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in), OPTIONAL :: gName  ! grid name
      integer, intent(in), OPTIONAL :: im  ! size of longitudinal dimension
      integer, intent(in), OPTIONAL :: jm  ! size of latitudinal  dimension
      integer, intent(in), OPTIONAL :: km  ! size of vertical dimension
      real, intent(in), OPTIONAL :: lon(:) ! longitude 
      real, intent(in), OPTIONAL :: lat(:) ! latitude 
      real, intent(in), OPTIONAL :: lev(:) ! Level   
      character(len=*), intent(in), OPTIONAL :: levUnit   
                                 ! units of level dimension, e.g., "hPa".
      character(len=*), intent(in), OPTIONAL :: coordinate
                                 ! string to indicate vertical coord
                                 ! (pressure, sigma, pressure_sigma)
      character(len=*), intent(in), OPTIONAL :: standardName 
                                 ! string for standard name
      character(len=*), intent(in), OPTIONAL :: formulaTerm  
                                 ! formula terms
      real, intent(in), OPTIONAL :: ak(:)     
                                 ! parameter for hybrid sigma prs coord.
      real, intent(in), OPTIONAL :: bk(:)     
                                 ! parameter for hybrid sigma prs coord.
      real, intent(in), OPTIONAL :: sigma(:)  
                                 ! parameter for sigma coordinate
      real, intent(in), OPTIONAL :: ptop              
                                 ! parameter for sigma coordinate
      character(len=*), intent(in), OPTIONAL :: ptopUnit   
                                 ! unit of ptop
      logical, intent(in), OPTIONAL :: twoDimLat       
                                 ! support 2D lat/lon or not
      logical, intent(in), OPTIONAL :: reduceGrid      
                                 ! support for reduced grid 
      logical, intent(in), OPTIONAL :: stnGrid      
                                 ! support for statio grid 

! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc ! Error return code:
                                           ! 0   all is well
                                           ! -1  Problem in setting lon
                                           ! -2  Problem in setting lat
                                           ! -3  Problem in setting lev
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_CFIOGrid), intent(inout) :: grid   ! CFIO grid   
!
! !DESCRIPTION:
!     Initializing a CFIO grid
!EOP
!------------------------------------------------------------------------------
       integer :: rtcode  = 0
       integer :: i, j

       if ( present(gName) ) grid%gName = gName
       if ( present(im) ) grid%im = im
       if ( present(jm) ) grid%jm = jm
       if ( present(km) ) grid%km = km

        if ( present(lon) ) then
           grid%im = size(lon)
           allocate(grid%lon(grid%im), stat = rtcode)
           grid%lon = lon
        else if ( present(im) ) then
                allocate(grid%lon(im), stat = rtcode)
                do i = 1, im
                   grid%lon(i) = 360./im * (i-1)
                end do
        end if    
        if (rtcode .ne. 0) then 
           print *, "problem in setting ESMF_CFIOGridSet:lon"
           rtcode = -1
           if ( present(rc) ) rc = rtcode
           return
        end if

        if ( present(lat) ) then
           grid%jm = size(lat)
           allocate(grid%lat(grid%jm), stat = rtcode)
           grid%lat = lat
        else if ( present(jm) ) then
                allocate(grid%lat(jm), stat = rtcode)
                do j = 1, jm
                   grid%lat(j) = 180./(jm-1) * (j-1) - 90
                end do
        end if    
        if (rtcode .ne. 0) then 
           print *, "problem in setting ESMF_CFIOGridSet:lat"
           rtcode = -2
           if ( present(rc) ) rc = rtcode
           return
        end if

        if ( present(lev) ) then
           grid%km = size(lev)
           allocate(grid%lev(grid%km), stat = rtcode)
           grid%lev = lev
           if (rtcode .ne. 0) then 
           print *, "problem in setting ESMF_CFIOGridSet:lev"
              rtcode = -3
              if ( present(rc) ) rc = rtcode
              return
           end if
        end if    

        if ( present(levUnit) ) grid%levUnits = levUnit

       if ( present(coordinate) ) grid%coordinate = coordinate
       if ( present(standardName) ) grid%standardName = standardName
       if ( present(formulaTerm) ) grid%formulaTerm = formulaTerm
       if ( present(ak) ) then
          allocate(grid%ak(size(ak)), stat=rtcode)
          grid%ak = ak 
       end if
       if ( present(bk) ) then
          allocate(grid%bk(size(bk)), stat=rtcode)
          grid%bk = bk 
       end if
       if ( present(ptop) ) then
          grid%ptop = ptop
       end if
       if ( present(ptopUnit) ) grid%ptopUnit = ptopUnit

       if ( present(stnGrid) ) then
          grid%stnGrid = stnGrid
       end if

       if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOGridSet


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOGridGet -- get grid info 

! !INTERFACE:
      subroutine ESMF_CFIOGridGet (grid, gName, im, jm, km, lat, lon, lev,  &
                                   coordinate, standardName, formulaTerm,   &
                                   levUnit, ak, bk, sigma, ptop, twoDimLat, &
                                   reduceGrid, stnGrid, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIOGrid), intent(in) :: grid   ! CFIO grid

! !OUTPUT PARAMETERS:
!
      character(len=*), intent(out), OPTIONAL :: gName  ! grid name
      integer, intent(out), OPTIONAL :: im  ! size of longitudinal dimension
      integer, intent(out), OPTIONAL :: jm  ! size of latitudinal  dimension
      integer, intent(out), OPTIONAL :: km  ! size of vertical dimension
      real, pointer, OPTIONAL :: lon(:) ! longitude
      real, pointer, OPTIONAL :: lat(:) ! latitude
      real, pointer, OPTIONAL :: lev(:) ! Level
      character(len=*), intent(out), OPTIONAL :: levUnit
                                 ! units of level dimension, e.g., "hPa".
      character(len=*), intent(out), OPTIONAL :: coordinate
                                 ! string to indicate vertical coord
                                 ! (pressure, sigma, pressure_sigma)
      character(len=*), intent(out), OPTIONAL :: standardName 
                                 ! string for standard name
      character(len=*), intent(out), OPTIONAL :: formulaTerm  
                                 ! formula terms
      real, intent(out), OPTIONAL :: ak(:)
                                 ! parameter for hybrid sigma prs coord.
      real, intent(out), OPTIONAL :: bk(:)
                                 ! parameter for hybrid sigma prs coord.
      real, intent(out), OPTIONAL :: sigma(:)
                                 ! parameter for sigma coordinate
      real, intent(out), OPTIONAL :: ptop
                                 ! parameter for sigma coordinate
      logical, intent(out), OPTIONAL :: twoDimLat
                                 ! support 2D lat/lon or not
      logical, intent(out), OPTIONAL :: reduceGrid
                                 ! support for reduced grid 
      logical, intent(out), OPTIONAL :: stnGrid
                                 ! support for station grid 
      integer, intent(out), OPTIONAL :: rc ! Error return code:
                                           ! 0   all is well
                                           ! -1  problem in getting lon
                                           ! -2  problem in getting lat
                                           ! -3  problem in getting lev
!
! !DESCRIPTION:
!     Get grid info
!EOP
!------------------------------------------------------------------------------
       integer :: rtcode = 0
                                                                                     
       if ( present(gName) ) gName = grid%gName
       if ( present(im) ) im = grid%im
       if ( present(jm) ) jm = grid%jm
       if ( present(km) ) km = grid%km 
                                                                                     
        if ( present(lon) ) then
           allocate(lon(size(grid%lon)), stat=rtcode)
           lon = grid%lon
           if (rtcode .ne. 0) then 
              print *, "problem in getting ESMF_CFIOGridGet:lon"
              rtcode = -1
              if ( present(rc) ) rc = rtcode
              return
           end if
        end if
        if ( present(lat) ) then
           allocate(lat(size(grid%lat)), stat=rtcode)
           lat = grid%lat
           if (rtcode .ne. 0) then 
              print *, "problem in getting ESMF_CFIOGridGet:lat"
              rtcode = -2
              if ( present(rc) ) rc = rtcode
              return
           end if
        end if
        if ( present(lev) ) then
           allocate(lev(size(grid%lev)), stat=rtcode)
           if (rtcode .ne. 0) then 
              print *, "problem in getting ESMF_CFIOGridGet:lev"
              rtcode = -3
              if ( present(rc) ) rc = rtcode
              return
           end if
           lev = grid%lev
        end if
         
        if ( present(ak) ) then
           if ( associated(grid%ak) ) then
              ak = grid%ak
           else
              print *, "ak was not defined in the input file"
           end if
        end if
        if ( present(bk) ) then
           if ( associated(grid%bk) ) then
              bk = grid%bk
           else
              print *, "bk was not defined in the input file"
           end if
        end if
        if ( present(sigma) ) then
           if ( associated(grid%sigma) ) then
              sigma = grid%sigma
           else
              print *, "sigam was not defined in the input file"
           end if
        end if
        if ( present(ptop) ) then
           ptop = grid%ptop
        end if

        if ( present(twoDimLat) ) then
           twoDimLat = grid%twoDimLat
        end if
        if ( present(reduceGrid) ) then
           reduceGrid = grid%reduceGrid
        end if
         
        if ( present(standardName) ) standardName = grid%standardName
        if ( present(coordinate) ) coordinate = grid%coordinate
        if ( present(formulaTerm) ) formulaTerm = grid%formulaTerm
        if ( present(levUnit) ) levUnit = grid%levUnits 
        if ( present(stnGrid) ) stnGrid = grid%stnGrid  

        if ( present(rc) ) rc = rtcode


      end subroutine ESMF_CFIOGridGet

!------------------------------------------------------------------------------!BOP
! !ROUTINE: ESMF_CFIOGridDestroy -- destructor for a CFIO grid object

! !INTERFACE:
      subroutine ESMF_CFIOGridDestroy (grid, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_CFIOGrid), intent(inout) :: grid ! CFIOGrid object


!
! !DESCRIPTION:
!     destructor for a CFIO grid object
!EOP
!------------------------------------------------------------------------------
      integer :: rtcode = 0

      if ( associated(grid%lon) ) deallocate(grid%lon, stat=rtcode)
      if ( associated(grid%lat) ) deallocate(grid%lat, stat=rtcode)
      if ( associated(grid%lev) ) deallocate(grid%lev, stat=rtcode)

      if ( associated(grid%ak) ) deallocate(grid%ak, stat=rtcode)
      if ( associated(grid%bk) ) deallocate(grid%bk, stat=rtcode)
      if ( associated(grid%sigma) ) deallocate(grid%sigma, stat=rtcode)

      if ( present(rc) ) rc = rtcode

      end subroutine ESMF_CFIOGridDestroy


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_Create_ -- Creates a DAO gridded file for writing
! 
! !DESCRIPTION: This routine is used to open a new file for a GFIO stream.
!
!, im, jm, km !INTERFACE:
!
      subroutine GFIO_Create_ ( cfio, rc )
!
! !USES:
!
      Implicit NONE  
!
! !INPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
      integer        fid     ! File handle
      integer        rc      ! Error return code:
                             ! 0  All is well
                             ! -1 Time increment is 0
                             ! -18 incorrect time increment
                             ! -30 can't open file
                             ! -31 error from ncddef
                             ! -32 error from ncvdef (dimension variable)
                             ! -33 error from ncapt(c) (dimension attribute)
                             ! -34 error from ncvdef (variable)
                             ! -35  error from ncapt(c) (variable attribute)
                             ! -36  error from ncaptc/ncapt (global attribute)
                             ! -37  error from ncendf
                             ! -38  error from ncvpt (dimension variable)
                             ! -39 Num of real var elements and Cnt differ

!
! !INPUT/OUTPUT PARAMETERS:
!
     type(ESMF_CFIO), intent(inout) :: cfio 
!
! !REVISION HISTORY: 
!
!EOP
!-------------------------------------------------------------------------

      ! REAL*4 variables for 32-bit output to netCDF file.

      integer :: im, jm, km, nst
      real*4, pointer :: lon_32(:), lat_32(:), levs_32(:)
      character(len=MVARLEN) :: levunits
      integer :: yyyymmdd_beg, hhmmss_beg, timinc
      real :: missing_val
      integer :: nvars
      character(len=MLEN), pointer :: vname(:)
      character(len=MVARLEN), pointer :: vtitle(:)
      character(len=MVARLEN), pointer :: vunits(:)
      integer, pointer :: kmvar(:), station(:)
      real, pointer :: valid_range(:,:), packing_range(:,:)
      integer, pointer :: akid(:), bkid(:), ptopid(:)
      integer :: prec
      integer, pointer ::  vid(:)

      real*4 amiss_32
      real*4 scale_32, offset_32
      real*4 high_32,low_32
      real*4, pointer :: ak_32(:), bk_32(:), layer(:)
      real*4 :: ptop_32(1)
      integer i, j
      integer timeid, timedim
      integer, pointer :: latid(:), lonid(:), stationid(:)
      integer, pointer :: levid(:), layerid(:)
      integer, pointer :: latdim(:), londim(:), stationdim(:)
      integer, pointer :: levdim(:), layerdim(:)
      integer, pointer :: gDims3D(:,:), gDims2D(:,:)
      integer dims3D(4), dims2D(3), dims1D(1), ptopdim
      integer corner(1), edges(1)
!      integer corner(4), edges(4)
      character*80 timeUnits 
      logical surfaceOnly
      character*8 strBuf
      character*14 dateString
      integer year,mon,day,hour,min,sec
      integer count
      integer maxLen
      integer rtcode
      logical :: aveFile = .false.
      character cellMthd
!      real*4 bndsdata(2)
      integer bndsid, dimsbnd(2), bndsdim
      integer ig
      integer ndim
      character cig

! Variables for packing

      integer*2 amiss_16
      real*4, pointer ::  pRange_32(:,:),vRange_32(:,:)
      logical packflag

! Set metadata strings.  These metadata values are specified in the 
! COARDS conventions

      character (len=50) :: lonName = "longitude"
      character (len=50) :: lonUnits = "degrees_east"
      character (len=50) :: latName = "latitude"
      character (len=50) :: latUnits = "degrees_north"
      character (len=50) :: levName = "vertical level"
!                           levUnits: specified by user in argument list
      character (len=50) :: layerName = "edges"
      character (len=50) :: layerUnits = "layer"
      character (len=50) :: timeName = "time"
!                           timeUnits: string is built below
      integer :: iCnt
      real*4, pointer :: realVarAtt(:)
      integer, pointer :: intVarAtt(:)
      real*4 :: scale_factor, add_offset
      character (len=50) :: nameLat, nameLon, nameLev, nameEdge
      character (len=50) :: nameAk, nameBk, namePtop, nameStation 
        
      nvars = cfio%mVars
      yyyymmdd_beg = cfio%date
      hhmmss_beg = cfio%begTime
      timinc = cfio%timeInc
      missing_val = cfio%varObjs(1)%amiss
      allocate(vname(nvars), vtitle(nvars), vunits(nvars), kmvar(nvars), &
            valid_range(2,nvars), packing_range(2,nvars), vid(nvars),    &
            vRange_32(2,nvars), pRange_32(2,nvars), stat = rtcode)

      allocate(latid(cfio%mGrids), lonid(cfio%mGrids),                   &
               levid(cfio%mGrids), layerid(cfio%mGrids),                 &
               latdim(cfio%mGrids), londim(cfio%mGrids),                 &
               levdim(cfio%mGrids), layerdim(cfio%mGrids),               &
               akid(cfio%mGrids),bkid(cfio%mGrids),ptopid(cfio%mGrids),  &
               gDims3D(4,cfio%mGrids), gDims2D(3,cfio%mGrids),           &
               stationdim(cfio%mGrids), stationid(cfio%mGrids) )

      do i=1,nvars
         vname(i) = cfio%varObjs(i)%vName
         vtitle(i) = cfio%varObjs(i)%vTitle
         vunits(i) = cfio%varObjs(i)%vUnits
         kmvar(i) = cfio%varObjs(i)%grid%km
         if ( cfio%varObjs(i)%twoDimVar ) kmvar(i) = 0
         valid_range(1, i) = cfio%varObjs(i)%validRange(1)
         valid_range(2, i) = cfio%varObjs(i)%validRange(2)
         packing_range(1, i) = cfio%varObjs(i)%packingRange(1)
         packing_range(2, i) = cfio%varObjs(i)%packingRange(2)
         if ( cfio%varObjs(i)%timAve ) then
            aveFile = .true.
            cellMthd = cfio%varObjs(i)%aveMethod
         end if
      enddo

      do j=1,nvars
        do i=1,2
           vRange_32(i,j) = valid_range(i,j)
           pRange_32(i,j) = packing_range(i,j)
        enddo
      enddo

      amiss_32 = cfio%varObjs(1)%amiss
      amiss_16 = PACK_FILL

! Variable initialization

      surfaceOnly = .TRUE.

! Basic error-checking.

      if (timinc .eq. 0) then
        rc=-1
        return
      endif

! Check to see if there is only surface data in this file definition

      do i=1,nvars
        if (kmvar(i) .NE. 0) then
          surfaceOnly = .FALSE.
          exit
        endif
      enddo

! Make NetCDF errors non-fatal, and do not warning messages.

      call ncpopt(0)

! Create the new NetCDF file. [ Enter define mode. ]

      fid = nccre (cfio%fName, NCCLOB, rc)
      if (err("Create: can't open file",rc,-30) .LT. 0) return

! Convert double-precision output variables to single-precision
   do ig = 1, cfio%mGrids
      im = cfio%grids(ig)%im
      jm = cfio%grids(ig)%jm
      km = cfio%grids(ig)%km
      if ( index(cfio%grids(ig)%gName, 'station') .gt. &
           0 ) then
         if (im .ne. jm) rtcode = err("It isn't station grid",-1,-1)
         nst = im
      end if

      levunits = trim(cfio%grids(ig)%levUnits)

      allocate(station(im))
      do i=1,im
         station(i) = i
      enddo

! Define dimensions.

      if ( ig .eq. 1 ) then
         if (cfio%mGrids .eq. 1) then
            nameLon = 'lon'
            nameLat = 'lat'
            nameLev = 'lev'
            nameEdge = 'edges'
            nameStation = 'station'
         else
            nameLon = 'lon0'
            nameLat = 'lat0'
            nameLev = 'lev0'
            nameEdge = 'edges0'
            nameStation = 'station0'
         end if
      else
        write (cig,"(I1)") ig-1
        nameLon = 'lon'//cig
        nameLat = 'lat'//cig
        nameLev = 'lev'//cig
        nameEdge = 'edges'//cig
        nameStation = 'station'//cig
      end if

      if (index(cfio%grids(ig)%gName,'station') .gt. 0) then
         stationdim(ig) = ncddef (fid, nameStation, im, rc)
         if (err("Create: error defining station",rc,-31) .LT. 0) return
!         londim(ig) = ncddef (fid, nameLon, im, rc)
!         if (err("Create: error defining lon",rc,-31) .LT. 0) return
!         latdim(ig) = ncddef (fid, nameLat, jm, rc)
!         if (err("Create: error defining lat",rc,-31) .LT. 0) return
      else
         londim(ig) = ncddef (fid, nameLon, im, rc)
         if (err("Create: error defining lon",rc,-31) .LT. 0) return
         latdim(ig) = ncddef (fid, nameLat, jm, rc)
         if (err("Create: error defining lat",rc,-31) .LT. 0) return
      end if

      if (.NOT. surfaceOnly) then
        levdim(ig) = ncddef (fid, nameLev, km, rc)
        if (err("Create: error defining lev",rc,-31) .LT. 0) return
      endif
      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_hybrid_sigma_pressure_coordinate' ) then
         layerdim(ig) = ncddef (fid, nameEdge, km+1, rc)
         if (err("Create: error defining edges",rc,-31) .LT. 0) return
      endif

! Define dimension variables.

      if (index(cfio%grids(ig)%gName,'station') .gt. 0) then
!         stationid(ig) = ncvdef (fid, nameStation, NCFLOAT, 1,       &
!                                 stationdim(ig), rc)
!         if (err("Create: error defining station",rc,-32) .LT. 0) return
         lonid(ig) = ncvdef (fid, nameLon, NCFLOAT, 1, stationdim(ig), rc)
         if (err("Create: error creating lon",rc,-32) .LT. 0) return
         latid(ig) = ncvdef (fid, nameLat, NCFLOAT, 1, stationdim(ig), rc)
         if (err("Create: error creating lat",rc,-32) .LT. 0) return
      else
         lonid(ig) = ncvdef (fid, nameLon, NCFLOAT, 1, londim(ig), rc)
         if (err("Create: error creating lon",rc,-32) .LT. 0) return
         latid(ig) = ncvdef (fid, nameLat, NCFLOAT, 1, latdim(ig), rc)
         if (err("Create: error creating lat",rc,-32) .LT. 0) return
      end if

      if (.NOT. surfaceOnly) then
        levid(ig) = ncvdef (fid, nameLev, NCFLOAT, 1, levdim(ig), rc)
        if (err("Create: error creating lev",rc,-32) .LT. 0) return
      endif
      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_hybrid_sigma_pressure_coordinate' ) then
         layerid(ig) = ncvdef (fid, nameEdge, NCFLOAT, 1, layerdim(ig), rc)
         if (err("Create: error creating edges",rc,-32) .LT. 0) return
      endif

! Set attributes for dimensions.

      call ncaptc (fid,lonid(ig),'long_name',NCCHAR,LEN_TRIM(lonName), &
                  lonName,rc)
      if (err("Create: error creating lon attribute",rc,-33) .LT. 0) &
        return
      call ncaptc (fid,lonid(ig),'units',NCCHAR,LEN_TRIM(lonUnits), &
                  lonUnits,rc)
      if (err("Create: error creating lon attribute",rc,-33) .LT. 0)  &
        return

      call ncaptc (fid,latid(ig),'long_name',NCCHAR,LEN_TRIM(latName),&
                  latName,rc)
      if (err("Create: error creating lat attribute",rc,-33) .LT. 0) &
        return
      call ncaptc (fid,latid(ig),'units',NCCHAR,LEN_TRIM(latUnits),&
                  latUnits,rc)
      if (err("Create: error creating lat attribute",rc,-33) .LT. 0) &
        return

      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_hybrid_sigma_pressure_coordinate' ) then
         call ncaptc (fid,layerid(ig),'long_name',NCCHAR,LEN_TRIM(layerName),&
                    layerName,rc)
         if (err("Create: error creating layer attribute",rc,-33) .LT. 0)&
            return
         call ncaptc (fid,layerid(ig),'units',NCCHAR,LEN_TRIM(layerUnits),&
                      layerUnits, rc)
         if (err("Create: error creating layer attribute",rc,-33) .LT. 0)&
           return
      endif
      if (.NOT. surfaceOnly) then
        call ncaptc (fid,levid(ig),'long_name',NCCHAR,LEN_TRIM(levName),&
                    levName,rc)
        if (err("Create: error creating lev attribute",rc,-33) .LT. 0)&
           return
        call ncaptc (fid,levid(ig),'units',NCCHAR,LEN_TRIM(levunits),&
                    levunits,rc)
        if (err("Create: error creating lev attribute",rc,-33) .LT. 0)&
           return
        call ncaptc (fid,levid(ig),'positive',NCCHAR,LEN_TRIM('down'),&
                    'down',rc)
        if (err("Create: error creating lev attribute",rc,-33) .LT. 0)&
           return
        call ncaptc (fid,levid(ig),'coordinate',NCCHAR,LEN_TRIM(  & 
                     cfio%grids(ig)%coordinate), cfio%grids(ig)%coordinate &
                     , rc)
        if (err("Create: error creating lev attribute",rc,-33) .LT. 0)&
           return
        call ncaptc (fid,levid(ig),'standard_name',NCCHAR,LEN_TRIM(  & 
                     cfio%grids(ig)%standardName),cfio%grids(ig)%standardName&
                     , rc)
        if (err("Create: error creating lev attribute",rc,-33) .LT. 0)&
           return
        if ( len(cfio%grids(ig)%formulaTerm) .gt. 0 .and. &
             trim(cfio%grids(ig)%formulaTerm) .ne. 'unknown') then
           call ncaptc (fid,levid(ig),'formula_term',NCCHAR,LEN_TRIM(  & 
                     cfio%grids(ig)%formulaTerm), cfio%grids(ig)%formulaTerm &
                     , rc)
           if (err("Create: error creating lev attribute",rc,-33) .LT. 0)&
           return
        end if
      endif
! end of mGrid loop
  end do

      timedim = ncddef(fid, 'time', NCUNLIM, rc)
      if (err("Create: error defining time",rc,-31) .LT. 0) return
      if ( aveFile ) then
         bndsdim = ncddef(fid, 'nv', 2, rc)
         if (err("Create: error defining time bounds",rc,-31) .LT. 0)&
              return
      end if
      do ig =1, cfio%mGrids
        if ( trim(cfio%grids(ig)%standardName) .eq.           &
           'atmosphere_hybrid_sigma_pressure_coordinate' .or. &
             trim(cfio%grids(ig)%standardName) .eq.           &
           'atmosphere_sigma_coordinate' ) then
            if (ig .eq. 1) then
              if (cfio%mGrids .eq. 1) then
                 ptopdim = ncddef (fid, "ptop", 1, rc)
              else
                 ptopdim = ncddef (fid, "ptop0", 1, rc)
              end if
            end if
        endif
      end do

      timeid = ncvdef (fid, 'time', NCLONG, 1, timedim, rc)
      if (err("Create: error creating time",rc,-32) .LT. 0) return
      call ncaptc (fid, timeid, 'long_name', NCCHAR, LEN_TRIM(timeName),&
                  timeName, rc)
      if (err("Create: error creating time attribute",rc,-33) .LT. 0)&
        return

!ams       write (dateString,200) yyyymmdd_beg, hhmmss_beg
!ams 200   format (I8,I6)
!ams       read (dateString,201) year,mon,day,hour,min,sec
!ams 201   format (I4,5I2)

      call GFIO_parseIntTime ( yyyymmdd_beg, year, mon, day )
      call GFIO_parseIntTime ( hhmmss_beg, hour,min,sec )

      write (timeUnits,202) year,mon,day,hour,min,sec
202   format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
              I2.2,':',I2.2)
      call ncaptc (fid, timeid, 'units', NCCHAR, LEN_TRIM(timeUnits),  &
                  timeUnits, rc)
      if (err("Create: error creating time attribute",rc,-33) .LT. 0) &
        return
      
!ams       write (strBuf,203) timinc
!ams 203   format (I6)
!ams       read (strBuf,204) hour, min, sec
!ams 204   format (3I2)

      call GFIO_parseIntTime ( timinc, hour, min, sec ) 

      if ( sec .NE. 0) then
        print *, 'GFIO_Create: Time increments not on minute', &
                ' boundaries are not currently allowed.'
        rc = -18
        return
      endif
      call ncapt (fid, timeid, 'time_increment', NCLONG, 1, timInc, rc)
      if (err("Create: error creating time attribute",rc,-33) .LT. 0) &
        return
      call ncapt (fid,timeid,'begin_date',NCLONG,1,yyyymmdd_beg,rc)
      if (err("Create: error creating time attribute",rc,-33) .LT. 0) &
        return
      call ncapt (fid,timeid,'begin_time',NCLONG,1,hhmmss_beg,rc)
      if (err("Create: error creating time attribute",rc,-33) .LT. 0) &
        return

      if ( aveFile ) then
         call ncapt (fid,timeid,'bounds',NCCHAR,9,'time_bnds',rc)
         if (err("Create: error creating time attribute",rc,-33) .LT. 0) &
             return
      end if

  do ig = 1, cfio%mGrids
      im = cfio%grids(ig)%im
      jm = cfio%grids(ig)%jm
      km = cfio%grids(ig)%km
      if ( index(cfio%grids(ig)%gName, 'station') .gt. &
           0 ) then
         if (im .ne. jm) rtcode = err("It isn't station grid",-1,-1)
         nst = im
      end if

      gDims3D(4,ig) = timedim
      gDims3D(3,ig) = levdim(ig)
      gDims3D(2,ig) = latdim(ig)
      gDims3D(1,ig) = londim(ig)
      
      gDims2D(3,ig) = timedim
      gDims2D(2,ig) = latdim(ig)
      gDims2D(1,ig) = londim(ig)

      if (index(cfio%grids(ig)%gName,'station') .gt. 0) then
         gDims3D(4,ig) = 0
         gDims3D(3,ig) = timedim
         gDims3D(2,ig) = levdim(ig)
         gDims3D(1,ig) = stationdim(ig)
      
         gDims2D(3,ig) = 0
         gDims2D(2,ig) = timedim
         gDims2D(1,ig) = stationdim(ig)
      end if

      if ( ig .eq. 1 ) then
         if (cfio%mGrids .eq. 1) then
           nameAk = 'ak'
           nameBk = 'bk'
           namePtop = 'ptop'
         else
           nameAk = 'ak0'
           nameBk = 'bk0'
           namePtop = 'ptop0'
         end if
      else
        write (cig,"(I1)") ig-1
        nameAk = 'ak'//cig
        nameBk = 'bk'//cig
        namePtop = 'ptop'//cig
      end if

      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_hybrid_sigma_pressure_coordinate' ) then

         dims1D = layerdim(ig)

         akid(ig) = ncvdef (fid, nameAk, NCFLOAT, 1, dims1D, rc)
         call ncaptc (fid,akid(ig),'long_name',NCCHAR,34,&
                     'ak component of hybrid coordinate',rc)
         if (err("Create: error creating ak attribute",rc,-33) .LT. 0)&
            return
         call ncaptc (fid,akid(ig),'units',NCCHAR,14,&
                     'dimensionless',rc)
         if (err("Create: error creating ak attribute",rc,-33) .LT. 0)&
            return

         bkid(ig) = ncvdef (fid, nameBk, NCFLOAT, 1, dims1D, rc)
         call ncaptc (fid,bkid(ig),'long_name',NCCHAR,34,&
                     'bk component of hybrid coordinate',rc)
         if (err("Create: error creating bk attribute",rc,-33) .LT. 0)&
            return
         call ncaptc (fid,bkid(ig),'units',NCCHAR,14,&
                     'dimensionless',rc)
         if (err("Create: error creating bk attribute",rc,-33) .LT. 0)&
            return

         ptopid(ig) = ncvdef (fid, namePtop, NCFLOAT, 1, ptopdim, rc)
         if (err("Create: error define ptopid",rc,-34) .LT. 0) return
         call ncaptc (fid,ptopid(ig),'long_name',NCCHAR,36,&
                     'ptop component of hybrid coordinate',rc)
         if (err("Create: error creating ptop attribute",rc,-33) .LT. 0)&
            return
         call ncaptc (fid,ptopid(ig),'units',NCCHAR,       &
                      len(trim(cfio%grids(ig)%ptopUnit)), &
                     trim(cfio%grids(ig)%ptopUnit),rc)
         if (err("Create: error creating ptop attribute",rc,-33) .LT. 0)&
            return
      end if

      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_sigma_coordinate' ) then
         ptopid(ig) = ncvdef (fid, namePtop, NCFLOAT, 1, ptopdim, rc)
         if (err("Create: error define ptopid",rc,-34) .LT. 0) return
         call ncaptc (fid,ptopid(ig),'long_name',NCCHAR,36,&
                     'ptop component of sigma coordinate',rc)
         if (err("Create: error creating ptop attribute",rc,-33) .LT. 0)&
            return
         call ncaptc (fid,ptopid(ig),'units',NCCHAR,      &
                     len(trim(cfio%grids(ig)%ptopUnit)), &
                     trim(cfio%grids(ig)%ptopUnit),rc)
         if (err("Create: error creating ptop attribute",rc,-33) .LT. 0)&
            return
      end if

! end of mGrids loop
  end do

      scale_32 = 1.0     ! No packing for now.
      offset_32 = 0.0    ! No packing for now.

! Set up packing attributes for each variable.  
! Define physical variables.  Set attributes for physical variables.

      do i=1,nvars
        scale_32 = 1.0                        ! default to no packing.
        offset_32 = 0.0
        if (pRange_32(1,i) .NE. amiss_32 .OR. pRange_32(2,i) .NE.  &
       amiss_32) then
          if (pRange_32(1,i) .GT. pRange_32(2,i)) then
            high_32 = pRange_32(1,i)
            low_32  = pRange_32(2,i)
          else
            high_32 = pRange_32(2,i)
            low_32  = pRange_32(1,i)
          endif
          scale_32 = (high_32 - low_32)/PACK_BITS*2
          offset_32 = high_32 - scale_32*PACK_BITS
          if (scale_32 .EQ. 0.0) then              ! If packing range is 0,
             scale_32 = 1.0                        ! no packing.
             offset_32 = 0.0
             packflag = .FALSE.
          else
             packflag = .TRUE.
          endif
        else
          packflag = .FALSE.
        endif
        do ig = 1, cfio%mGrids
           if (trim(cfio%varObjs(i)%grid%gName) .eq.              &
               trim(cfio%grids(ig)%gName)) then
              dims3D = gDims3D(:,ig)
              dims2D = gDims2D(:,ig)
           end if
        end do
        if ( kmvar(i) .eq. 0 ) then
          ndim = 3
          if (index(cfio%varObjs(i)%grid%gName,'station') .gt. 0) ndim = 2
          if (packflag) then
            vid(i) = ncvdef (fid, vname(i), NCSHORT, ndim, dims2D, rc)
          else if (cfio%prec .EQ. 1) then
            vid(i) = ncvdef (fid, vname(i), NCDOUBLE, ndim, dims2D, rc)
          else
            vid(i) = ncvdef (fid, vname(i), NCFLOAT, ndim, dims2D, rc)
          endif
        else
          ndim = 4
          if (index(cfio%varObjs(i)%grid%gName,'station') .gt. 0) ndim = 3
          if (packflag) then
            vid(i) = ncvdef (fid, vname(i), NCSHORT, ndim, dims3D, rc)
          else if (cfio%prec .EQ. 1) then
            vid(i) = ncvdef (fid, vname(i), NCDOUBLE, ndim, dims3D, rc)
          else
            vid(i) = ncvdef (fid, vname(i), NCFLOAT, ndim, dims3D, rc)
          endif
        endif
        if (err("Create: error defining variable",rc,-34) .LT. 0)  &
         return

        call ncaptc (fid, vid(i), 'long_name', NCCHAR,  &
                    LEN_TRIM(vtitle(i)),vtitle(i), rc)
        if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return
        call ncaptc (fid, vid(i), 'units', NCCHAR,  &
                    LEN_TRIM(vunits(i)),vunits(i), rc)
        if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return

        if (packflag) then
          call ncapt (fid,vid(i),'_FillValue',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return
          if ( scale_32 .ne. 1.0 .or. offset_32 .ne. 0.0 ) then
          call ncapt (fid,vid(i),'scale_factor',NCFLOAT,1,scale_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
              .LT. 0) return
          call ncapt (fid,vid(i),'add_offset',NCFLOAT,1,offset_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
              .LT. 0) return
          call ncapt (fid,vid(i),'packmin',NCFLOAT,1,low_32,rc)
          if (err("Create: error defining variable attribute",rc,-35)  &
             .LT. 0) return
          call ncapt (fid,vid(i),'packmax',NCFLOAT,1,high_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
             .LT. 0) return
          end if
          call ncapt (fid,vid(i),'missing_value',NCSHORT,1,amiss_16,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return
          call ncapt (fid,vid(i),'fmissing_value',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return
        else
          call ncapt (fid,vid(i),'_FillValue',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return
          if ( scale_32 .ne. 1.0 .or. offset_32 .ne. 0.0 ) then
          call ncapt (fid,vid(i),'scale_factor',NCFLOAT,1,scale_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
              .LT. 0) return
          call ncapt (fid,vid(i),'add_offset',NCFLOAT,1,offset_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
              .LT. 0) return
          end if
          call ncapt (fid,vid(i),'missing_value',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return
          call ncapt (fid,vid(i),'fmissing_value',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
          .LT. 0) return

! ADDED BY BYIN for more variable meta data
         cfio%fid = fid
!        get real variable attributes from rList
         do iCnt = 1, cfio%mVars
            if ( associated(cfio%varObjs(i)%rList) ) then
               call getMaxLenCnt(maxLen, cfio%varObjs(i)%nVarAttReal, &
                                 rList=cfio%varObjs(i)%rList)
               count = cfio%varObjs(i)%nVarAttReal
               allocate(cfio%varObjs(i)%attRealNames(count),     &
                     cfio%varObjs(i)%attRealCnts(count),      &
                     cfio%varObjs(i)%varAttReals(count,maxLen), stat=rtcode)
               call getList(rList=cfio%varObjs(i)%rList,    &
                       realAttNames=cfio%varObjs(i)%attRealNames, &
                       realAttCnts=cfio%varObjs(i)%attRealCnts,   &
                       realAtts=cfio%varObjs(i)%varAttReals)
            end if
         end do

!        write real variable attributes to output file
         do iCnt = 1, cfio%varObjs(i)%nVarAttReal
            allocate(realVarAtt(size(cfio%varObjs(i)%varAttReals)/ &
                     cfio%varObjs(i)%nVarAttReal), stat=rc)
            realVarAtt = cfio%varObjs(i)%varAttReals(iCnt,:)
            if (cfio%varObjs(i)%attRealCnts(iCnt) .ne. size(realVarAtt)) then
              rc=err("FileCreate: Num of real var elements and Cnt differ",-39,-39) 
              return
            end if
            call ncapt (cfio%fid,vid(i),cfio%varObjs(i)%attRealNames(iCnt),&
                        NCFLOAT, cfio%varObjs(i)%attRealCnts(iCnt),        &
                        realVarAtt, rc)
            if (err("FileCreate: error from ncapt for real att",rc,-35) &
                .LT. 0) return
            deallocate(realVarAtt)
         end do

!        get integer variable attributes from iList
         do iCnt = 1, cfio%mVars
            if ( associated(cfio%varObjs(i)%iList) ) then
               call getMaxLenCnt(maxLen, cfio%varObjs(i)%nVarAttInt, &
                                 iList=cfio%varObjs(i)%iList)
               count = cfio%varObjs(i)%nVarAttInt
               allocate(cfio%varObjs(i)%attIntNames(count),     &
                     cfio%varObjs(i)%attIntCnts(count),      &
                     cfio%varObjs(i)%varAttInts(count,maxLen), stat=rtcode)
               call getList(iList=cfio%varObjs(i)%iList,    &
                       intAttNames=cfio%varObjs(i)%attIntNames, &
                       intAttCnts=cfio%varObjs(i)%attIntCnts,   &
                       intAtts=cfio%varObjs(i)%varAttInts)
            end if
         end do

!        write int variable attributes to output file
         do iCnt = 1, cfio%varObjs(i)%nVarAttInt
            allocate(intVarAtt(size(cfio%varObjs(i)%varAttInts)/ &
                     cfio%varObjs(i)%nVarAttInt), stat=rc)
            intVarAtt = cfio%varObjs(i)%varAttInts(iCnt,:)
            if (cfio%varObjs(i)%attIntCnts(iCnt) .gt. size(intVarAtt)) then
              rc=err("FileCreate: Num of int var elements and Cnt differ",-39,-39) 
              return
            end if
            call ncapt (cfio%fid,vid(i),cfio%varObjs(i)%attIntNames(iCnt),&
                        NCLONG, cfio%varObjs(i)%attIntCnts(iCnt),         &
                        intVarAtt, rc)
            if (err("FileCreate: error from ncapt for int att",rc,-35) &
                .LT. 0) return
            deallocate(intVarAtt)
         end do

!        get char variable attributes from cList
         do iCnt = 1, cfio%mVars
            if ( associated(cfio%varObjs(i)%cList) ) then
               call getMaxLenCnt(maxLen, cfio%varObjs(i)%nVarAttChar, &
                                 cList=cfio%varObjs(i)%cList)
               count = cfio%varObjs(i)%nVarAttChar
               allocate(cfio%varObjs(i)%attCharNames(count),     &
                     cfio%varObjs(i)%attCharCnts(count),      &
                     cfio%varObjs(i)%varAttChars(count), stat=rtcode)
               call getList(cList=cfio%varObjs(i)%cList,    &
                       charAttNames=cfio%varObjs(i)%attCharNames, &
                       charAttCnts=cfio%varObjs(i)%attCharCnts,   &
                       charAtts=cfio%varObjs(i)%varAttChars)
            end if
         end do

!        write char variable attributes to output file
         do iCnt = 1, cfio%varObjs(i)%nVarAttChar
            call ncapt (cfio%fid,vid(i),cfio%varObjs(i)%attCharNames(iCnt),&
                        NCCHAR, cfio%varObjs(i)%attCharCnts(iCnt),         &
                        cfio%varObjs(i)%varAttChars(iCnt), rc)
            if (err("FileCreate: error from ncapt for char att",rc,-35) &
                .LT. 0) return
         end do

!         write scaleFactor, addOffSet, and standardName to output

!         if ( cfio%varObjs(i)%scaleFactor /= 0 ) then
            scale_factor = cfio%varObjs(i)%scaleFactor 
            call ncapt (cfio%fid, vid(i), 'scale_factor', NCFLOAT,   &
                        1, scale_factor, rc)
            if (err("FileCreate: error from ncapt for scale_factor",rc,-35) &
                .LT. 0) return
!         end if
!         if ( cfio%varObjs(i)%addOffSet /= 0 ) then
            add_offset = cfio%varObjs(i)%addOffSet   
            call ncapt (cfio%fid, vid(i), 'add_offset', NCFLOAT,   &
                        1, add_offset, rc)
            if (err("FileCreate: error from ncapt for add_offset",rc,-35) &
                .LT. 0) return
!           end if
              
         if ( LEN_TRIM(cfio%varObjs(i)%standardName) .gt. 0 ) then
            call ncaptc (cfio%fid, vid(i), 'standard_name', NCCHAR,   &
                        LEN_TRIM(cfio%varObjs(i)%standardName),       &
                        cfio%varObjs(i)%standardName, rc)
            if (err("FileCreate: error from ncapt for standard_name",rc,-35) &
                .LT. 0) return
         end if
        end if

        if (vRange_32(1,i) .NE. amiss_32 .OR. vRange_32(2,i) .NE.  &
           amiss_32) then
          if (vRange_32(1,i) .GT. vRange_32(2,i)) then
            high_32 = vRange_32(1,i)
            low_32  = vRange_32(2,i)
          else
            high_32 = vRange_32(2,i)
            low_32  = vRange_32(1,i)
          endif
          call ncapt (fid,vid(i),'vmin',NCFLOAT,1,low_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
             .LT. 0) return
          call ncapt (fid,vid(i),'vmax',NCFLOAT,1,high_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
             .LT. 0) return
        else
          call ncapt (fid,vid(i),'vmin',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
             .LT. 0) return
          call ncapt (fid,vid(i),'vmax',NCFLOAT,1,amiss_32,rc)
          if (err("Create: error defining variable attribute",rc,-35) &
             .LT. 0) return

        endif

        call ncapt (fid,vid(i),'valid_range',NCFLOAT,2,vRange_32(:,i),rc)
        if (err("Create: error defining variable attribute",rc,-35) &
           .LT. 0) return

        if ( cfio%varObjs(i)%timAve ) then
           call ncaptc (fid, vid(i), 'cell_methods', NCCHAR,  &
                    len(trim(cfio%varObjs(i)%cellMthd))+6,     &
                    'time: '//trim(cfio%varObjs(i)%cellMthd), rc)
           if (err("Create: error defining variable attribute",rc,-35) &
                   .LT. 0) return
        end if
      enddo
 
      if ( aveFile ) then
         dimsbnd(1) = bndsdim
         dimsbnd(2) = timedim
         bndsid = ncvdef (fid, 'time_bnds', NCFLOAT, 2, dimsbnd, rc)
      end if

! Exit define mode.

      call ncendf (fid, rc)
      if (err("Create: error exiting define mode",rc,-37) .LT. 0)  &
       return

! Write out dimension variables.

  do ig = 1, cfio%mGrids
      im = cfio%grids(ig)%im
      jm = cfio%grids(ig)%jm
      km = cfio%grids(ig)%km

      allocate(lon_32(im), lat_32(jm), levs_32(km), ak_32(km+1),         &
            bk_32(km+1), layer(km+1), stat = rtcode) 

      ptop_32(1) = cfio%grids(ig)%ptop
      do i=1,im
         lon_32(i) = cfio%grids(ig)%lon(i)
      enddo
      do i=1,jm
         lat_32(i) = cfio%grids(ig)%lat(i)
      enddo
      do i=1,km
         levs_32(i) = cfio%grids(ig)%lev(i)
      enddo
      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_hybrid_sigma_pressure_coordinate' ) then
         if (associated(cfio%grids(ig)%ak) .and. &
             associated(cfio%grids(ig)%bk) ) then
            do i=1,km+1
               layer(i) = i
               ak_32(i) = cfio%grids(ig)%ak(i)
               bk_32(i) = cfio%grids(ig)%bk(i)
            enddo
         else
             if (err(": ak or bk is not set",-1,-1) .lt. 0 ) return
         end if
      end if

      corner(1) = 1
      edges(1) = im
      call ncvpt (fid, lonid(ig), corner, edges, lon_32, rc)
      if (err("Create: error writing lons",rc,-38) .LT. 0) return

      corner(1) = 1
      edges(1) = jm
      call ncvpt (fid, latid(ig), corner, edges, lat_32, rc)
      if (err("Create: error writing lats",rc,-38) .LT. 0) return

      if (.NOT. surfaceOnly) then
        corner(1) = 1
        edges(1) = km
        call ncvpt (fid, levid(ig), corner, edges, levs_32, rc)
        if (err("Create: error writing levs",rc,-38) .LT. 0) return
      endif

      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_hybrid_sigma_pressure_coordinate' ) then
        corner(1) = 1
        edges(1) = 1
        call ncvpt (fid, ptopid(ig), corner, edges, ptop_32, rc)
        corner(1) = 1
        edges(1) = km+1
        call ncvpt (fid, layerid(ig), corner, edges, layer, rc)
        if (err("Create: error writing layers",rc,-38) .LT. 0) return
        call ncvpt (fid, akid(ig), corner, edges, ak_32, rc)
        call ncvpt (fid, bkid(ig), corner, edges, bk_32, rc)
      endif

      if ( trim(cfio%grids(ig)%standardName) .eq. &
           'atmosphere_sigma_coordinate' ) then
        corner(1) = 1
        edges(1) = 1
        call ncvpt (fid, ptopid(ig), corner, edges, ptop_32, rc)
      endif

! end of mGrids loop
 end do
      corner(1) = 1
      edges(1) = 1
      call ncvpt (fid, timeid, corner, edges, 0, rc)
      if (err("Create: error writing times",rc,-38) .LT. 0) return

      rc=0
      return
      end subroutine GFIO_Create_


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: writeBnds -- write time bounds

! !INTERFACE:
      subroutine writeBnds(cfio, vName, date, curTime, rc)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type (ESMF_CFIO), intent(in) :: cfio
      character(len=*), intent(in) :: vName
      integer, intent(in) :: date 
      integer, intent(in) :: curTime
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                                                ! 1   ...
!
!
! !DESCRIPTION:
!     write time bounds for time averaging variable
!EOP
!------------------------------------------------------------------------------

      integer :: vid, corner(4), edges(4)
      integer :: hour, min, sec, incSecs, timeIndex
      integer :: seconds, timeinc, curSecs
      real*4 :: bndsdata(2)
      character*8 :: strBuf
      integer :: i, rtcode=0

!     make sure user provides the right variable name
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
      if ( cfio%varObjs(i)%timAve ) then
         seconds = DiffDate (cfio%date, cfio%begTime, date, curTime)
         timeinc = cfio%timeInc

!ams          write (strBuf,203) timeinc
!ams 203      format (I6)
!ams          read (strBuf,204) hour, min, sec
!ams 204      format (3I2)

         call GFIO_parseIntTime ( timeinc, hour, min, sec ) 

         incSecs = hour*3600 + min*60 + sec

!ams         write (strBuf,203) curTime
!ams         read (strBuf,204) hour, min, sec

         call GFIO_parseIntTime ( curTime, hour, min, sec ) 

         curSecs = hour*3600 + min*60 + sec
                                                                     
         timeIndex = seconds/incSecs + 1
         corner(1) = 1
         corner(2) = timeIndex
         edges(1) = 2
         edges(2) = 1
         bndsdata(1) = (-incSecs + curSecs)/60.
         bndsdata(2) = curSecs/60.
         if ( cfio%varObjs(i)%aveMethod .eq. 'c' ) then
            bndsdata(1) = (-incSecs/2. + curSecs)/60.
            bndsdata(2) = (incSecs/2. + curSecs)/60.
         end if
         if ( cfio%varObjs(i)%aveMethod .eq. 'd' ) then
            bndsdata(1) = curSecs/60.
            bndsdata(2) = (incSecs + curSecs)/60.
         end if
 
         vid = ncvid (cfio%fid, 'time_bnds', rtcode)
         if ( rtcode .ne. 0 ) then 
            print *, "ncvid failed in ncvid for time_bnds"
            if ( present(rc) ) rc = rtcode
            return
         end if
         call ncvpt (cfio%fid, vid, corner, edges, bndsdata, rtcode)
         if ( rtcode .ne. 0 ) then 
            print *, "ncvid failed in ncvpt for time_bnds"
            if ( present(rc) ) rc = rtcode
            return
         end if
      end if
 
      if ( present(rc) ) rc = rtcode

      end subroutine writeBnds

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarReadT3D_ -- Read a variable from an existing file
                                                                                                              
! !INTERFACE:
      subroutine ESMF_CFIOVarReadT3D_ ( cfio, vName, field, &
                                        timeString, cfio2, rc )
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      type(ESMF_CFIO), intent(inOut), OPTIONAL :: cfio2  ! second CFIO obj
      character(len=*), intent(in) :: timeString
                                  ! string expression for date and time
                                                                                                        
!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:,:,:)             ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -46  error from ncvgt
                         !  rc = -48  error from ncinq
                         !  rc = -52  error from ncvinq
                         !  rc = -99  must specify date/curTime of timeString
!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------

    integer :: date_, curTime_

!     Resolve date/time
!     -----------------
      date_ = -1
      curTime_ = -1
      call strToInt(timeString,date_,curTime_)

      if ( date_ < 0 .OR. curTime_ < 0 ) then
           if ( present(rc) ) rc = -99
           return
      end if

      call ESMF_CFIOVarReadT3D__ ( cfio, vName, date_, curTime_, field, & 
                                   cfio2=cfio2, rc=rc )

    end subroutine ESMF_CFIOVarReadT3D_

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarReadT3D_ -- Read a variable from an existing file
                                                                                                              
! !INTERFACE:

      subroutine ESMF_CFIOVarReadT3D__(cfio, vName, date, curTime, field, rc, cfio2)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      integer, intent(in) :: date                 ! yyyymmdd
      integer, intent(in) :: curTime              ! hhmmss
      type(ESMF_CFIO), intent(inOut), OPTIONAL :: cfio2  ! second CFIO obj
                                                                                                              
!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:,:,:)             ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -46  error from ncvgt
                         !  rc = -48  error from ncinq
                         !  rc = -52  error from ncvinq
!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------

      integer rtcode
      integer begDate, begTime, incSecs, timeIndex1, timeIndex2
      integer secs, secs1, secs2, nymd1, nymd2, nhms1, nhms2
      integer i, j, k
      integer im, jm, km
                                                                                         
      real    alpha, amiss
      real, pointer ::  field2(:,:,:) => null() ! workspace for interpolation

      rtcode = 0

!     find the right variable obj.
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
      im = cfio%varObjs(i)%grid%im
      jm = cfio%varObjs(i)%grid%jm
      km = cfio%varObjs(i)%grid%km
      if (km .lt. 1) km = 1

      if ( .not. associated(field) ) allocate(field(im,jm,km))

!     Get beginning time & date.  Calculate offset seconds from start.
!     ----------------------------------------------------------------
      call GetBegDateTime ( cfio%fid, begDate, begTime, incSecs, rtcode )
      if (err("GetVar: could not determine begin_date/begin_time",rtcode,-44)&
         .NE. 0) go to 999
                                                                                         
      secs = DiffDate (begDate, begTime, date, curTime)
                                                                                         
!      if (date .LT. begDate .OR. (begDate .EQ. date .AND.  &
!         curTime .LT. begTime) .or. secs .LT. 0) then
!         rc = -7
!         return
!      endif
 
!     Determine brackting times
!     -------------------------
      if ( secs >= 0 ) then
         timeIndex1 = secs/incSecs + 1
      else
         timeIndex1 = secs/incSecs
      end if
      timeIndex2 = timeIndex1 + 1
      secs1 = (timeIndex1-1) * incSecs
      secs2 = (timeIndex2-1) * incSecs
      call GetDate ( begDate, begTime, secs1, nymd1, nhms1, rtcode )
      call GetDate ( begDate, begTime, secs2, nymd2, nhms2, rtcode )
 
!     Read grids at first time with GetVar()
!     --------------------------------------
      call ESMF_CFIOVarRead(cfio, vName, field, date=nymd1, curtime=nhms1,  rc=rtcode)
      if ( rtcode .ne. 0 ) goto 999
                                                                    
      if ( secs1 .eq. secs ) goto 999   ! no interpolation needed

      allocate(field2(im,jm,km))
                                                                                     
!     Read grids at second time with GetVar()
!     ---------------------------------------
      call ESMF_CFIOVarRead(cfio, vName, field2, date=nymd2, curtime=nhms2, rc=rtcode)
      if ( rtcode .ne. 0 ) then
         if ( present(cfio2) )     &
            call ESMF_CFIOVarRead(cfio2, vName, field2, &
                                  date=nymd2, curtime=nhms2, rc=rtcode)
         if ( rtcode .ne. 0 ) return
      end if
                                                                                         
!     Get missing value
!     -----------------
      amiss = GFIO_GetMissing ( cfio%fid, rtcode )
      if ( rtcode .ne. 0 ) goto 999

!     Do interpolation
!     ----------------
      alpha = float(secs - secs1)/float(secs2 - secs1)
!ams  print *, ' nymd = ', nymd1, nymd2
!ams  print *, ' nhms = ', nhms1, nhms2
!ams  print *, 'alpha = ', alpha
      do k = 1, km
         do j = 1, jm
            do i = 1, im
               if ( abs(field(i,j,k)-amiss) .gt. 0.001 .and.   &
                    abs(field2(i,j,k)-amiss) .gt. 0.001 ) then
                  field(i,j,k) = field(i,j,k)        &
                             + alpha * (field2(i,j,k) - field(i,j,k))
               else
                  field(i,j,k) = amiss
               end if
            end do
         end do
      end do
                                                        
      rtcode = 0

!     All done
!     --------
999   continue
      if ( associated(field2) ) deallocate(field2)
      if ( present(rc) ) rc = rtcode
                                                                         
      end subroutine ESMF_CFIOVarReadT3D__


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarReadT2D_ -- Read a variable from an existing file
                                                                                                              
! !INTERFACE:
      subroutine ESMF_CFIOVarReadT2D_ ( cfio, vName, field, &
                                        timeString, cfio2, rc )
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      type(ESMF_CFIO), intent(inOut), OPTIONAL :: cfio2  ! second CFIO obj
      character(len=*), intent(in) :: timeString
                                  ! string expression for date and time
                                                                                                        
!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:,:)               ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -46  error from ncvgt
                         !  rc = -48  error from ncinq
                         !  rc = -52  error from ncvinq
                         !  rc = -99  must specify date/curTime of timeString
!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------

    integer :: date_, curTime_

!     Resolve date/time
!     -----------------
      date_ = -1
      curTime_ = -1
      call strToInt(timeString,date_,curTime_)
      if ( date_ < 0 .OR. curTime_ < 0 ) then
           if ( present(rc) ) rc = -99
           return
      end if

      call ESMF_CFIOVarReadT2D__ ( cfio, vName, date_, curTime_, field, &
                                   cfio2=cfio2, rc=rc )

    end subroutine ESMF_CFIOVarReadT2D_

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ESMF_CFIOVarReadT2D_ -- Read a variable from an existing file
                                                                                                              
! !INTERFACE:

      subroutine ESMF_CFIOVarReadT2D__(cfio, vName, date, curTime, field, rc, cfio2)
!
! !ARGUMENTS:
!
! !INPUT PARAMETERS:
!
      type(ESMF_CFIO), intent(inOut) :: cfio      ! a CFIO obj
      character(len=*), intent(in) :: vName       ! variable name
      integer, intent(in) :: date                 ! yyyymmdd
      integer, intent(in) :: curTime              ! hhmmss
      type(ESMF_CFIO), intent(inOut), OPTIONAL :: cfio2  ! second CFIO obj
                                                                                                              
!
! !OUTPUT PARAMETERS:
!
      real, pointer :: field(:,:)             ! array contains data
      integer, intent(out), OPTIONAL :: rc      ! Error return code:
                                                ! 0   all is well
                         !  rc = -2  time is inconsistent with increment
                         !  rc = -3  number of levels is incompatible with file
                         !  rc = -4  im is incompatible with file
                         !  rc = -5  jm is incompatible with file
                         !  rc = -6  time must fall on a minute boundary
                         !  rc = -7  error in diffdate
                         !  rc = -12  error determining default precision
                         !  rc = -13  error determining variable type
                         !  rc = -19  unable to identify coordinate variable
                         !  rc = -38  error from ncvpt (dimension variable)
                         !  rc = -40  error from ncvid
                         !  rc = -41  error from ncdid or ncdinq (lat or lon)
                         !  rc = -42  error from ncdid or ncdinq (lev)
                         !  rc = -43  error from ncvid (time variable)
                         !  rc = -44  error from ncagt (time attribute)
                         !  rc = -46  error from ncvgt
                         !  rc = -48  error from ncinq
                         !  rc = -52  error from ncvinq
!
! !DESCRIPTION:
!     Read a variable from an existing file
!EOP
!------------------------------------------------------------------------------

      integer rtcode
      integer begDate, begTime, incSecs, timeIndex1, timeIndex2
      integer secs, secs1, secs2, nymd1, nymd2, nhms1, nhms2
      integer i, j, k
      integer im, jm, km
                                                                                         
      real    alpha, amiss
      real, pointer ::  field2(:,:) => null() ! workspace for interpolation

      rtcode = 0

!     find the right variable obj.
      do i = 1, cfio%mVars
         if ( trim(vName) .eq. trim(cfio%varObjs(i)%vName) ) exit
      end do
      im = cfio%varObjs(i)%grid%im
      jm = cfio%varObjs(i)%grid%jm
      km = cfio%varObjs(i)%grid%km
      if (km .lt. 1) km = 1

      if ( .not. associated(field) ) allocate(field(im,jm))

!     Get beginning time & date.  Calculate offset seconds from start.
!     ----------------------------------------------------------------
      call GetBegDateTime ( cfio%fid, begDate, begTime, incSecs, rtcode )
      if (err("GetVar: could not determine begin_date/begin_time",rtcode,-44)&
         .NE. 0) go to 999
                                                                                         
      secs = DiffDate (begDate, begTime, date, curTime)
                                                                                         
!      if (date .LT. begDate .OR. (begDate .EQ. date .AND.  &
!         curTime .LT. begTime) .or. secs .LT. 0) then
!         rc = -7
!         return
!      endif
 
!     Determine brackting times
!     -------------------------
      if ( secs >= 0 ) then
         timeIndex1 = secs/incSecs + 1
      else
         timeIndex1 = secs/incSecs
      end if
      timeIndex2 = timeIndex1 + 1
      secs1 = (timeIndex1-1) * incSecs
      secs2 = (timeIndex2-1) * incSecs
      call GetDate ( begDate, begTime, secs1, nymd1, nhms1, rtcode )
      call GetDate ( begDate, begTime, secs2, nymd2, nhms2, rtcode )
 
!     Read grids at first time with GetVar()
!     --------------------------------------
      call ESMF_CFIOVarRead(cfio, vName, field, date=nymd1, curtime=nhms1,  rc=rtcode)
      if ( rtcode .ne. 0 ) goto 999
                                                                    
      if ( secs1 .eq. secs ) goto 999   ! no interpolation needed

      allocate(field2(im,jm))
                                                                                     
!     Read grids at second time with GetVar()
!     ---------------------------------------
      call ESMF_CFIOVarRead(cfio, vName, field2, date=nymd2, curtime=nhms2, rc=rtcode)
      if ( rtcode .ne. 0 ) then
         if ( present(cfio2) )     &
            call ESMF_CFIOVarRead(cfio2, vName, field2, &
                                  date=nymd2, curtime=nhms2, rc=rtcode)
         if ( rtcode .ne. 0 ) return
      end if
                                                                                         
!     Get missing value
!     -----------------
      amiss = GFIO_GetMissing ( cfio%fid, rtcode )
      if ( rtcode .ne. 0 ) goto 999

!     Do interpolation
!     ----------------
      alpha = float(secs - secs1)/float(secs2 - secs1)
      do j = 1, jm
         do i = 1, im
            if ( abs(field(i,j)-amiss) .gt. 0.001 .and.   &
                 abs(field2(i,j)-amiss) .gt. 0.001 ) then
               field(i,j) = field(i,j) + alpha * (field2(i,j) - field(i,j))
            else
               field(i,j) = amiss
            end if
         end do
      end do
                                                        
      rtcode = 0

!     All done
!     --------
999   continue
      if ( associated(field2) ) deallocate(field2)
      if ( present(rc) ) rc = rtcode
                                                                         
      end subroutine ESMF_CFIOVarReadT2D__

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ESMF_CFIODownBit - GRIB-based compression pre-conditioner 
!
! !INTERFACE:

   subroutine ESMF_CFIODownBit3D_ ( x, xr, nbits, undef, flops, rc )

     implicit NONE

!
! !INPUT PARAMETERS:
!
     real, intent(in)    ::  x(:,:,:)       ! input array 
     integer, intent(in) :: nbits           ! number of bits per word to retain
                                            ! - no action if nbits<1
     real, OPTIONAL, intent(in) :: undef    ! missing value
     logical, OPTIONAL, intent(in) :: flops ! if true, uses slower float point
                                            !  based algorithm
!
! !OUTPUT PARAMETERS:
!
     real*4, intent(out)   :: xr(:,:,:) ! precision reduced array; can
!                                       ! share storage with input array
                                        ! if it has same kind
     integer, intent(out)  :: rc        ! error code
                                        !  = 0 - all is well
                                        ! /= 0 - something went wrong 
!
! !DESCRIPTION:  
!
!  This routine returns a lower precision version of the input array
!  {\tt x} which retains {\tt nbits} of precision. See routine
!  {\tt ESMF\_CFIODownBit2D} for additional details. This version for
!  rank 3 arrays, calls {\tt ESMF\_CFIODownBit2D()} for each vertical
!  level.
!
! !REVISION HISTORY:
!
!  06Dec2006  da Silva  Initial version.
!
!EOP
!------------------------------------------------------------------------------

   integer :: k

   do k = lbound(x,3), ubound(x,3)
      call ESMF_CFIODownBit2D_ ( x(:,:,k), xr(:,:,k), nbits, &
                                 undef=undef, flops=flops, rc=rc )
   end do

   end subroutine ESMF_CFIODownBit3D_


!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ESMF_CFIODownBit - GRIB-based compression pre-conditioner 
!
! !INTERFACE:

   subroutine ESMF_CFIODownBit2D_ ( x, xr, nbits, undef, flops, rc )

     implicit NONE

!
! !INPUT PARAMETERS:
!
     real, intent(in)    ::  x(:,:)         ! input array 
     integer, intent(in) :: nbits           ! number of bits per word to retain
     real, OPTIONAL, intent(in) :: undef    ! missing value
     logical, OPTIONAL, intent(in) :: flops ! if true, uses slower float point
                                            !  based algorithm
!
! !OUTPUT PARAMETERS:
!
     real*4, intent(out)   :: xr(:,:)   ! precision reduced array; can
!                                       !  share storage with input array
!                                       !  if it has same kind
     integer, intent(out)  :: rc        ! error code
                                        !  = 0 - all is well
                                        ! /= 0 - something went wrong 
!
! !DESCRIPTION:  
!
!  This routine returns a lower precision version of the input array
!  {\tt x} which retains {\tt nbits} of precision. Two algorithms are
!  implemented: 1) a fast one writen in C which downgrades precision
!  by shifting {\tt xbits = 24 - nbits} bits of the mantissa, and 2) a slower
!  float point based algorithm which is the same algorithm as GRIB 
!  with fixed number of bits packing. Notice that as in GRIB the scaling 
!  factor is forced to be a power of 2 rather than a generic float.  
!  Using this power of 2 binary scaling has the advantage of improving 
!  the GZIP compression rates.
!
!  This routine returns an array of the same type and kind as the input array, 
!  so no data compression has taken place. The goal here is to reduce the
!  entropy in the input array, thereby improving compression rates 
!  by the lossless algorithms implemented internally by HDF-4/5 when writing 
!  these data to a file. In fact, these GZIP'ed and pre-conditioned files 
!  have sizes comparable to the equivalent GRIB file, while being a bonafide 
!  self-describing HDF/NetCDF file.
!
! !TO DO:
!
!  Perhaps implement GRIB decimal scaling (variable number of bits).
!
! !REVISION HISTORY:
!
!  06Dec2006  da Silva  Initial version.
!
!EOP
!------------------------------------------------------------------------------

    integer   :: E, xbits, has_undef
    real*4    :: scale, xmin, xmax, tol, undef_
    logical   :: shave_mantissa
    integer, external :: ShaveMantissa32

    rc = 0

!   Defaults for optinal arguments
!   ------------------------------
    if ( present(undef) ) then
         undef_ = undef
         has_undef = 1
    else
         undef_ = huge(1.0)   ! why not?
         has_undef = 0
    endif
    if ( present(flops) ) then
         shave_mantissa = .not. flops
    else
         shave_mantissa = .true.
    endif

!   Fast, bit shifting in C
!   -----------------------
    if ( shave_mantissa ) then

       xr = x   ! compiled r8 this will convert to r4.
       xbits = 24 - nbits
       rc = ShaveMantissa32 ( xr, xr, size(x), xbits, has_undef, undef_ )
       return

!   Slow, flops in FORTRAN (GRIB inspired)
!   --------------------------------------
    else 

       if ( nbits < 1 ) then
          xr = x
          rc = 1
          return
       end if

       tol = 0.0001 * undef_
       xmin = minval(x,mask=(abs(undef_-x)>tol))
       xr = x - xmin     ! As in GRIB, force non-negative values 
       xmax = maxval(xr,mask=(abs(undef_-x)>tol)) ! max of positive

       if ( xmax <= 0.0 ) then
            xr = x
            rc = 0
            return  ! this means field is constant 
       end if

       E = nint(log(xmax)/log(2.)) - nbits ! GRIB binary scale factor
       scale = 2.**E                       ! GRIB requires power of 2

       if ( present(undef) ) then
          where ( abs(x - undef_) > tol )
             xr  = xmin + nint(xr/scale) * scale
          endwhere
       else
          xr  = xmin + nint(xr/scale) * scale
       end if

    end if

   end subroutine ESMF_CFIODownBit2D_

!..........................................................................

      end module ESMF_CFIOMod

