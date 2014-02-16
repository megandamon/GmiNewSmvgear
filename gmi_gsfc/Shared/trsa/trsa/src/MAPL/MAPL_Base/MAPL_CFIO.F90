!  $Id: MAPL_CFIO.F90,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

#include "MAPL_Generic.h"

module MAPL_CFIOMod

!BOP

! !MODULE: MAPL_CFIO

! !DESCRIPTION:  {\tt MAPL\_CFIO}  interfaces CFIO to the ESMF data types.
!      It currently includes read-write support for ESMF Fields and States,
!      and read support for ESMF Fields and Fortran arrays. It has only four methods:
!      MAPL\_CFIOCreate, MAPL\_CFIOWrite, MAPL\_CFIORead, MAPL\_CFIODestroy. Except
!      for MAPL\_Read, all work on the MAPL\_CFIO object. Reading is done
!      directly from a file to the appropriate ESMF object.
!
!      {\tt MAPL\_CFIO} is designed for two modes of I/O: self-describing formats (SDF),
!      of which it supports HDF-4 and netcdf-3, and flat files, which includes support
!      for GrADS readable files. The GrADS support is still und4er construction. There
!      are also plans to add GRIB support. 
!
!      In SDF mode, capability of {\tt MAPL\_CFIO} depends on which library (HDF or netcdf)
!      is linked with the application, since both cannot be used because
!      name conflicts between these two libraries. If netcdf is linked, only netcdf 
!      files may be read or written. If HDF is linked, both HDF and netcdf files may
!      be read, but only HDF files may be written.
!      

! !USES:

  use ESMF_Mod

  use MAPL_BaseMod
  use MAPL_CommsMod
  use MAPL_ConstantsMod
  use ESMF_CFIOMod  
  use ESMF_CFIOUtilMod
  use MAPL_IOMod
  use MAPL_HorzTransformMod

  implicit none
  private

! !PUBLIC TYPES:

  public MAPL_CFIO

! !PUBLIC MEMBER FUNCTIONS:

! MAPL-style names
! ----------------

  public MAPL_CFIOCreate
  public MAPL_CFIOWrite
  public MAPL_CFIORead
  public MAPL_CFIODestroy

! ESMF-style names
! ----------------

  public ESMF_ioRead
  public ESMF_ioCreate
  public ESMF_ioWrite  
  public ESMF_ioDestroy

! !INTERFACE:

  interface MAPL_CFIOCreate
     module procedure MAPL_CFIOCreateFromBundle
     module procedure MAPL_CFIOCreateFromState
  end interface

  interface MAPL_CFIOWrite
     module procedure MAPL_CFIOWriteState
     module procedure MAPL_CFIOWriteBundle
  end interface

  interface MAPL_CFIORead
     module procedure MAPL_CFIOReadState
     module procedure MAPL_CFIOReadBundle
     module procedure MAPL_CFIOReadField
     module procedure MAPL_CFIOReadArray3D
     module procedure MAPL_CFIOReadArray2D
  end interface

! ESMF consistent naming convention
! ---------------------------------

  interface ESMF_ioCreate
     module procedure MAPL_CFIOCreateFromBundle
     module procedure MAPL_CFIOCreateFromState
  end interface

  interface ESMF_ioRead
     module procedure MAPL_CFIOReadState
     module procedure MAPL_CFIOReadBundle
     module procedure MAPL_CFIOReadField
     module procedure MAPL_CFIOReadArray3D
     module procedure MAPL_CFIOReadArray2D
  end interface

  interface ESMF_ioWrite
     module procedure MAPL_CFIOWriteState
     module procedure MAPL_CFIOWriteBundle
  end interface

  interface ESMF_ioDestroy
     module procedure MAPL_CFIODestroy
  end interface

!EOP


  type Ptr3Arr
     real, pointer              :: Ptr(:,:,:)
  end type Ptr3Arr

  type MAPL_CFIO
     private
     logical                    :: CREATED=.false.
     character(len=ESMF_MAXSTR) :: NAME
     type(ESMF_CFIO)            :: CFIO
     integer                    :: XYOFFSET
     real                       :: VSCALE
     character(len=ESMF_MAXSTR) :: VVAR
     integer, pointer           :: RESOLUTION(:)
     type(ESMF_TIMEINTERVAL)    :: OFFSET
     type(ESMF_CLOCK)           :: CLOCK
     type(ESMF_BUNDLE)          :: BUNDLE
  end type MAPL_CFIO

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: MAPL_CFIOCreateFromBundle

! !DESCRIPTION: Creates a MAPL\_CFIO object from a Bundle. The MAPL\_CFIO
!           objects is opaque and its properties can only be set by this
!           method at creation. Its properties cannot be queried. The object
!           is used only as a handle in write operations. It is not needed for 
!           reading.
!
!           Its non-optional arguments associate a {\tt NAME}, an ESMF  {\tt BUNDLE},
!           and a {\tt CLOCK} with the object. An ESMF TimeInterval {\tt OFFSET} is
!           an optional argument that sets an offset between
!           the time on the clock when eriting and the time stamp used for the data 
!           (defaults to no offset).
!           The {\tt format} optional argument determines whether the write
!           will use the linked self-describing format (SDF) library (HDF or netcdf) or write
!           GrADS readable flat files. Currently only the SDF library option
!           is supported.
!
!           The remaining (optional) arguments are
!           especialized and used primarily to support MAPL\_History, or to provide
!           documentation in the form of character strings that will be placed
!           in corresponding attributes in the SDF file.
!
! !INTERFACE:

  subroutine MAPL_CFIOCreateFromBundle ( MCFIO, NAME, CLOCK, BUNDLE, OFFSET,    &
                                         RESOLUTION, LEVELS, DESCR,             &
                                         XYOFFSET, VCOORD, VUNIT, VSCALE,       &
                                         SOURCE, INSTITUTION, COMMENT, CONTACT, &
                                         FORMAT,                                &
                                         RC )

! !ARGUMENTS:

    type(MAPL_CFIO),             intent(OUT) :: MCFIO
    character(LEN=*),            intent(IN)  :: NAME
    type(ESMF_BUNDLE),           intent(IN)  :: BUNDLE
    type(ESMF_CLOCK),            intent(IN)  :: CLOCK
    type(ESMF_TIMEINTERVAL), &
                     optional,   intent(IN)  :: OFFSET
    integer,         optional,   pointer     :: RESOLUTION(:)
    real,            optional,   pointer     :: LEVELS(:)
    character(LEN=*),optional,   intent(IN)  :: DESCR
    integer,         optional,   intent(IN)  :: XYOFFSET
    real,            optional,   intent(IN)  :: VSCALE
    character(len=*),optional,   intent(IN)  :: VUNIT     
    character(len=*),optional,   intent(IN)  :: VCOORD     
    character(len=*),optional,   intent(IN)  :: source
    character(len=*),optional,   intent(IN)  :: institution     
    character(len=*),optional,   intent(IN)  :: comment
    character(len=*),optional,   intent(IN)  :: contact     
    character(len=*),optional,   intent(IN)  :: format

    integer,         optional,   intent(OUT) :: RC


!EOP

    character(len=ESMF_MAXSTR)     :: Iam="MAPL_CFIOCreateFromBundle"
    integer          :: STATUS

    type(ESMF_FIELD ) :: FIELD
    type(ESMF_GRID  ) :: ESMFGRID
    type(ESMF_ARRAY ), pointer :: EARRAY(:)
    type(ESMF_ARRAY ) :: ARRAY
    type(ESMF_TIME  ) :: TIME
!!ALT    type(ESMF_CFIO  ), pointer :: CFIO

    type(ESMF_CFIOVarInfo), pointer :: vars(:)
    type(ESMF_CFIOGrid),    pointer :: grid

    real, pointer  :: lats(:,:)
    real, pointer  :: lons(:,:)
    real, pointer  :: lats1d(:)
    real, pointer  :: lons1d(:)
    real, pointer  :: latsLocal(:,:)
    real, pointer  :: lonsLocal(:,:)
    real, pointer  :: lev (:  )
    real, pointer  :: ulevels(:  )
    integer        :: L, WriteInterval, counts(5)
    integer        :: NumVars
    integer        :: IM,JM,LM
    integer        :: gridRank
    integer        :: arrayRank
    logical        :: twoDimVar
    real           :: RANGE(2)
    real           :: xoff, yoff
    integer        :: hours, mins, secs, timeInc
    integer        :: I, J
    integer        :: VLOCATION
    integer        :: LOCATION
    logical        :: HASVLEVS

    real(KIND=8), pointer  :: R8D2(:,:)

    character(len=ESMF_MAXSTR)     :: VarName, Vvar, Vunits
    character(len=ESMF_MAXSTR)     :: LongName
    character(len=ESMF_MAXSTR)     :: Units
    character(len=ESMF_MAXSTR)     :: StartTime
    logical                        :: change_resolution
    real(KIND=8)                   :: dlam, dphi

    type (ESMF_Array), target :: tarray(2)


    character(len=esmf_maxstr)  :: fName
    character(len=esmf_maxstr)  :: Usource
    character(len=esmf_maxstr)  :: Uinstitution     
    character(len=esmf_maxstr)  :: Ucomment
    character(len=esmf_maxstr)  :: Ucontact     
    character(len=esmf_maxstr)  :: Uformat
    character(len=esmf_maxstr)  :: Utitle

    if(present(source)) then
       Usource = source
    else
       Usource = ""
    endif

    if(present(institution)) then
       Uinstitution = institution
    else
       Uinstitution = ""
    endif

    if(present(comment)) then
       Ucomment = comment
    else
       Ucomment = ""
    endif

    if(present(contact)) then
       Ucontact = contact
    else
       Ucontact = ""
    endif

    if(present(format)) then
       Uformat = format
    else
       Uformat = "HDF"
    endif

    if(present(descr )) then
       Utitle  = descr 
    else
       Utitle  = ""
    endif

    if(present(LEVELS)) then
       ulevels => LEVELS
    else
       nullify(ulevels)
    endif

    earray=>tarray

    MCFIO%NAME   = NAME 
    MCFIO%CLOCK  = CLOCK
    MCFIO%BUNDLE = BUNDLE
    

    Vvar = ""
    if(present(Vcoord)) then
       Vvar = adjustl(vcoord)
       if    (Vvar(1:3)=='log') then
          Vvar  = adjustl(Vvar(index(vvar,'(')+1:index(vvar,')')-1))
       elseif(Vvar(1:3)=='pow') then 
          Vvar  = adjustl(Vvar(index(vvar,'(')+1:index(vvar,',')-1))
       endif
    end if

    MCFIO%VVAR   = VVAR

    if(present(vscale)) then
       MCFIO%Vscale = Vscale
    else
       MCFIO%Vscale = 1.0
    endif

    call ESMF_BundleGet(BUNDLE, FieldCount=NumVars, RC=STATUS)
    VERIFY_(STATUS)

    ASSERT_(NumVars>0)

    LOCATION = MAPL_VLocationNone
    DO I = 1, NumVars

       call ESMF_BundleGetField(BUNDLE, I, FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGetAttribute(FIELD, NAME="REFRESH_INTERVAL", VALUE=WriteInterval, RC=STATUS)
       if ( status /= 0 ) WriteInterval = 0

       call ESMF_FieldGetAttribute(FIELD, NAME="VLOCATION", VALUE=VLOCATION, RC=STATUS)
       if ( status /= 0 ) VLOCATION = MAPL_VLocationCenter
       
       if (.not.associated(ULEVELS) ) then
          if ( VLOCATION /= MAPL_VLocationNone ) then
               if (Location == MAPL_VLocationNone) then
                  LOCATION = VLOCATION ! first time initialization
               elseif (LOCATION /= VLOCATION) then
                  print *, 'ERROR: Mixed Vlocation in CFIO mode not allowed unless LEVELS are specified'
                  RETURN_(ESMF_FAILURE)
               end if
          end if
       else
          if ( VLOCATION /= MAPL_VLocationNone ) LOCATION = VLOCATION
       end if

    end DO

! ALT: Now we are using the last field from the stream to get the grid           

    call ESMF_FieldGet(FIELD, grid=ESMFGRID, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_GridGet(ESMFGRID, dimCount=gridRank, rc=STATUS)
    VERIFY_(STATUS)

    if (gridRank == 3) then
       call ESMF_GridGet(ESMFGRID, horzRelLoc=ESMF_CELL_CENTER, &
            vertRelLoc=ESMF_CELL_CELL, &
            globalCellCountPerDim=COUNTS,  RC=STATUS)
       VERIFY_(STATUS)
       IM = COUNTS(1)
       JM = COUNTS(2)
       call ESMF_GridGetDELocalInfo(ESMFGRID, &
            horzRelLoc=ESMF_CELL_CENTER, &
            vertRelLoc=ESMF_CELL_CELL, &
            localCellCountPerDim=COUNTS,RC=STATUS)
       VERIFY_(STATUS)
       LM = COUNTS(3)

    else
       call ESMF_GridGet(ESMFGRID, horzRelLoc=ESMF_CELL_CENTER, &
            globalCellCountPerDim=COUNTS, RC=STATUS)
       VERIFY_(STATUS)
       IM = COUNTS(1)
       JM = COUNTS(2)
       LM = 1
    endif

    call ESMF_GridGetCoord(ESMFGRID,          &
            horzRelLoc =ESMF_CELL_CENTER,     &
            centerCoord=EARRAY,               &
                                     RC=STATUS) 
    VERIFY_(STATUS) 

    change_resolution = .FALSE.

    if (present(RESOLUTION)) then
       if (associated(RESOLUTION)) then
          if (IM /= resolution(1) .or. JM /= resolution(2)) then
             change_resolution = .TRUE.
          end if
          IM = resolution(1)
          JM = resolution(2)
       endif
       MCFIO%RESOLUTION=>RESOLUTION
    else
       nullify(MCFIO%RESOLUTION)
    endif

    allocate(LONS1D(IM), STAT=status)
    VERIFY_(status)
    allocate(LATS1D(JM), STAT=status)
    VERIFY_(status)

    if (change_resolution) then

       if(present(xyoffset)) then
          select case(xyoffset)
          case(0)
             xoff = 0.0
             yoff = 0.0
          case(1)
             xoff = 0.5
             yoff = 0.0
          case(2)
             xoff = 0.0
             yoff = 0.5
          case(3)
             xoff = 0.5
             yoff = 0.5
          case default
             ASSERT_(.false.)
          end select
          mcfio%xyoffset = xyoffset
       else
          xoff = 0.0
          yoff = 0.0
          mcfio%xyoffset = 0
       endif

       dlam = 2*MAPL_PI/IM

       if(yoff>0) then
          dphi = MAPL_PI/(JM  )
       else
          dphi = MAPL_PI/(JM-1)
       end if
  
       do i=1,IM
          lons1d(i) = -MAPL_PI     + (i-1+xoff)*dlam
       enddo
       
       do j=1,JM
          lats1d(j) = -MAPL_PI/2.0 + (j-1+yoff)*dphi
       enddo

       lats1d = lats1d*(180/MAPL_PI)
       lons1d = lons1d*(180/MAPL_PI)
    else

       allocate(LONS(IM,JM), STAT=status)
       VERIFY_(status)
       allocate(LATS(IM,JM), STAT=status)
       VERIFY_(status)

       call ESMF_ArrayGetData(EARRAY(1), R8D2, RC=STATUS)
       VERIFY_(STATUS)

       allocate(LONSLocal(size(R8D2,1),size(R8D2,2)), STAT=status)
       VERIFY_(status)
       allocate(LATSLocal(size(R8D2,1),size(R8D2,2)), STAT=status)
       VERIFY_(status)

       LONSLOCAL = R8D2*(180/MAPL_PI)
       call ArrayGather(LONSLOCAL, LONS, ESMFGRID, RC=STATUS)

       call ESMF_ArrayGetData(EARRAY(2), R8D2, RC=STATUS)
       VERIFY_(STATUS)

       LATSLOCAL = R8D2*(180/MAPL_PI)
       call ArrayGather(LATSLOCAL, LATS, ESMFGRID, RC=STATUS)
       VERIFY_(STATUS)

       LONS1D = LONS(:,1)
       LATS1D = LATS(1,:)

       DEALLOCATE(LONS)
       DEALLOCATE(LATS)
       DEALLOCATE(LONSlocal)
       DEALLOCATE(LATSlocal)

    endif

    allocate(GRID)

    GRID = ESMF_CFIOGridCreate(gName=trim(NAME)//"Grid",                   RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_CFIOGridSet(GRID, IM=IM, JM=JM, LON=LONS1D,LAT=LATS1D, RC=STATUS)
    VERIFY_(STATUS)

    deallocate(LONS1D)
    deallocate(LATS1D)

    if (LOCATION == MAPL_VLocationNone) then
       LM = 1 
       gridRank = 2
    end if

    if (gridRank == 3) then
       if (associated(ULEVELS)) then
          LM = size(ULEVELS)
       elseif(LOCATION == MAPL_VLocationEdge) then
          LM = LM+1
       end if

       allocate(LEV(LM))

       if (associated(ULEVELS)) then
          LEV = ULEVELS
       else
          LEV = (/(L, L=1,LM)/)
       end if

       if(present(VUNIT)) then
          vunits = trim(vunit)
       else
          vunits = ""
       endif

       if(Vvar/="") then
          call ESMF_BundleGetField(BUNDLE, name=Vvar, field=FIELD, RC=STATUS)
          VERIFY_(STATUS)
          if( trim(vunits).eq."" ) then
             call ESMF_FieldGetAttribute(FIELD, NAME="UNITS", VALUE=units, RC=STATUS)
             VERIFY_(STATUS)
             call ESMF_CFIOGridSet(grid, lev=lev, levUnit=trim(units), RC=STATUS)
             VERIFY_(STATUS)
          else
             call ESMF_CFIOGridSet(grid, lev=lev, levUnit=trim(vunits), RC=STATUS)
             VERIFY_(STATUS)
          endif
          call ESMF_CFIOGridSet(grid, standardName =trim(VVAR)//'_level', RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_CFIOGridSet(grid, coordinate   =trim(VVAR),   RC=STATUS)
          VERIFY_(STATUS)
       else
          call ESMF_CFIOGridSet(grid, lev=lev, levUnit='layer',     RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_CFIOGridSet(grid, standardName ='model_layers', RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_CFIOGridSet(grid, coordinate   ='eta',          RC=STATUS)
          VERIFY_(STATUS)
       end if

       deallocate(LEV)
    else
       call ESMF_CFIOGridSet(grid, standardName ='2d_fields',    RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_CFIOGridSet(grid, coordinate   ='N/A',          RC=STATUS)
       VERIFY_(STATUS)
    end if

! Create variable objects
!------------------------

    allocate(vars(NumVars))

    do L=1,NumVars

       call ESMF_BundleGetField(BUNDLE, L, FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet (FIELD, NAME=VarName,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGetAttribute(FIELD, NAME="LONG_NAME",VALUE=LongName, RC=STATUS)
!ams       VERIFY_(STATUS)
       if ( status /= 0 ) LongName = VarName
       call ESMF_FieldGetAttribute(FIELD, NAME="UNITS"    ,VALUE=Units,    RC=STATUS)
!ams       VERIFY_(STATUS)
       if ( status /= 0 ) Units = 'unknown'

       call ESMF_FieldGet (FIELD, ARRAY=array,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ArrayGet(array, rank=arrayRank, rc=status)
       VERIFY_(STATUS)
       if (arrayRank == 2) then
          twoDimVar = .true.
       else if (arrayRank == 3) then
          twoDimVar = .false.
       else
          RETURN_(ESMF_FAILURE)
       end if

       VARS(L) = ESMF_CFIOVarInfoCreate(vName=trim(VARNAME),               RC=STATUS)
       VERIFY_(STATUS)

! ALT: initialize range (NEEDs ATTENTION)
       RANGE(2) = 1.0e30
       RANGE(1) = -RANGE(2)

       call ESMF_CFIOVarInfoSet(VARS(L),           &
            vName        = VarName,                &
            vTitle       = LongName,               &
            grid         = GRID,                   &
            amiss        = MAPL_Undef,             &
            scaleFactor  = 1.,                     &
            addOffSet    = 0.,                     &
            standardName = LongName,               &
            twoDimVar    = twoDimVar,              &
            validRange   = RANGE,                  &
            vUnits       = UNITS,                  &
                                         RC=STATUS )
       VERIFY_(STATUS)
    end do

    call ESMF_ClockGet(CLOCK, CurrTime =TIME, RC=STATUS)
    VERIFY_(STATUS)

    if(present(OFFSET)) then
       TIME = TIME - OFFSET
       MCFIO%OFFSET = OFFSET
    else
      call ESMF_TimeIntervalSet( MCFIO%OFFSET, S=0, rc=status )
      VERIFY_(STATUS)
    endif
    call ESMF_TimeGet (TIME,  timeString=StartTime, RC=STATUS)
    VERIFY_(STATUS)

! Create CFIO object
!-------------------

    MCFIO%cfio =  ESMF_CFIOCreate(cfioObjName=trim(Name))

!ALT: currently CFIO and GFIO expect timeIncrement to be in HHMMSS format
!     this imposes severe limitations to the frequency of the output:
!     no writes should be done less frequently than once every 4 days (99 hours)

    if (writeInterval == 0     ) writeInterval = 21600  ! Default of 6 hours Output Frequency
    if (writeInterval > 4*86400) then
       RETURN_(ESMF_FAILURE)
    end if
    
    hours = writeInterval/3600
    writeInterval = writeInterval - 3600*hours
    mins = writeInterval/60
    secs = writeInterval - 60*mins

    timeinc = 10000*hours + 100*mins + secs

    if(format=="HDF") then
       fName = trim(Name)//'.hdf'
    else
       fName = trim(Name)
    endif

! Set global attributes
!----------------------

    call ESMF_CFIOSet(MCFIO%CFIO,                                 &
         fName       = fName,                                     &
         varObjs     = VARS,                                      &
         grid        = GRID,                                      &
         TimeString  = trim(StartTime),                           &
         timeInc     = timeInc,                                   &
         title       = trim(Utitle),                              &
         source      = Usource,                                   &
         history     = 'File written by MAPL_CFIO',               &
         institution = Uinstitution,                              &
         convention  = "COARDS",                                  &
         contact     = Ucontact,                                  &
         references  = "see MAPL documentation",                  &
         comment     = Ucomment, prec=0,                          & 
         RC=STATUS )
    VERIFY_(STATUS)

! Create FILE
!------------

    if (MAPL_AM_I_ROOT()) then
       call ESMF_CFIOFileCreate(MCFIO%CFIO, RC=STATUS)
       VERIFY_(STATUS)
    end if

! All Done
!---------

    deallocate(vars)
    deallocate(grid)

    MCFIO%CREATED = .true.

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_CFIOCreateFromBundle


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: MAPL_CFIOCreateFromState

! !DESCRIPTION: Creates a MAPL\_CFIO object from a State. States are
!    written by ``serializing'' all Fields in them, whether they are
!    directly in the State or are contained within a hierarchy of embedded
!    Bundles and States, into a single Bundle.
!
!     The Method optionally
!    returns a pointer to the serialized ESMF Bundle, but this is not
!    needed for MAPL\_Write operations. Otherwise arguments are the same as for 
!    CreateFromBundle.
!
!


! !INTERFACE:

  subroutine MAPL_CFIOCreateFromState ( MCFIO, NAME, CLOCK, STATE, OFFSET,  &
                                        RESOLUTION, LEVELS, DESCR, BUNDLE, &
                                        XYOFFSET, VCOORD, VUNIT, VSCALE,   &
                                        SOURCE, INSTITUTION, COMMENT, CONTACT, &
                                        FORMAT,                                &
                                        RC )

! !ARGUMENTS:

    type(MAPL_CFIO),             intent(OUT) :: MCFIO
    character(LEN=*),            intent(IN)  :: NAME
    type(ESMF_State),            intent(IN)  :: STATE
    type(ESMF_Clock),            intent(IN)  :: CLOCK
    type(ESMF_TimeInterval), &
                     optional,   intent(IN)  :: OFFSET
    integer,         optional,   pointer     :: RESOLUTION(:)
    real,            optional,   pointer     :: LEVELS(:)
    character(LEN=*),optional,   intent(IN)  :: DESCR
    real,            optional,   intent(IN)  :: VSCALE
    character(len=*),optional,   intent(IN)  :: VUNIT     
    character(len=*),optional,   intent(IN)  :: VCOORD     
    integer,         optional,   intent(IN)  :: XYOFFSET
    character(len=*),optional,   intent(IN)  :: source
    character(len=*),optional,   intent(IN)  :: institution     
    character(len=*),optional,   intent(IN)  :: comment
    character(len=*),optional,   intent(IN)  :: contact     
    character(len=*),optional,   intent(IN)  :: format
    type(ESMF_Bundle), optional,  pointer    :: BUNDLE
    integer, optional,           intent(OUT) :: RC

!
!EOP

    character(len=ESMF_MAXSTR)     :: Iam="MAPL_CFIOCreate"
    integer          :: STATUS

! Locals

    type(ESMF_Bundle), target :: tBUNDLE


!   Create an empty bundle
!   ----------------------

    tBUNDLE = ESMF_BundleCreate ( name=Iam, rc=STATUS )
    VERIFY_(STATUS)
    
!   Serialize the state
!   -------------------

    call ESMFL_BundleAddState_ ( tBUNDLE, STATE, rc=STATUS, VALIDATE=.true. )
    VERIFY_(STATUS)

!   Create the mapl_CFIO object
!   ----------------------

    call MAPL_CFIOCreateFromBundle ( MCFIO, NAME, CLOCK, tBUNDLE,         &
                                     OFFSET = OFFSET,                    & 
                                     RESOLUTION=RESOLUTION,              &
                                     LEVELS=LEVELS,                      &
                                     DESCR=DESCR,                        &
                                     XYOFFSET= XYOFFSET,                 &
                                     VCOORD  = VCOORD,                   &
                                     VUNIT   = VUNIT,                    &
                                     VSCALE  = VSCALE,                   &
                                     SOURCE  = SOURCE, &
                                     INSTITUTION = INSTITUTION, &
                                     COMMENT = COMMENT, &
                                     CONTACT = CONTACT, &
                                     FORMAT = FORMAT,                    &
                                                               RC=STATUS )
    VERIFY_(STATUS)

    if ( present(BUNDLE) ) then
         BUNDLE => tBUNDLE
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_CFIOCreateFromState
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: MAPL_CFIOWriteBundle

! !DESCRIPTION:  Writes an ESMF Bundle to a File. Only the MAPL\_CFIO
!       object is a required argument, and it must have been created.
!       {\tt CLOCK, BUNDLE} can be used to override the choice
!       made at creation, but this is of dubious value, particularly
!       for {\tt BUNDLE} since it must be excatly conformant with the
!       creation {\tt BUNDLE}. {\tt NBITS} if the number of bits of 
!       the mantissa to retain. This is used to write files with degraded
!       precision, which can then be compressed with standard utilities.
!       The default is no degradation of precision.
!
! !INTERFACE:

  subroutine MAPL_CFIOWriteBundle ( MCFIO, CLOCK, Bundle, &
                                    VERBOSE, NBITS, RC    )

! !ARGUMENTS:

    type(MAPL_CFIO  ),           intent(INOUT) :: MCFIO
    type(ESMF_CLOCK), optional,  intent(IN   ) :: CLOCK
    type(ESMF_BUNDLE),optional,  intent(IN   ) :: BUNDLE

    logical,          optional,  intent(  IN)  :: VERBOSE
    integer,          optional,  intent(  IN)  :: NBITS
    integer,          optional,  intent(  OUT) :: RC

!
!EOP

    character(len=*), parameter:: Iam="MAPL_CFIOWriteBundle"
    integer                    :: status

    type(ESMF_FIELD)           :: FIELD
    type(ESMF_ARRAY )          :: ARRAY
    type(ESMF_TIME )           :: TIME
    type(ESMF_GRID )           :: ESMFGRID
    type(ESMF_deLAYOUT)        :: LAYOUT
    type(ESMF_CFIOGrid), pointer :: GRID
    type(MAPL_HorzTransform)   :: Trans
    type(ESMF_BUNDLE)          :: UBUNDLE
    type(ESMF_CLOCK)           :: uCLOCK

    integer                    :: L, K, NumVars
    integer                    :: IM, JM, LM
    integer                    :: IM0, JM0, LM0
    integer                    :: K1, KK
    integer                    :: arrayRank
    integer                    :: gridRank
    character(len=ESMF_MAXSTR) :: VarName, DATE, vvar
    real, pointer              :: PTR2(:,:), PTR3(:,:,:), lev(:)
    real, pointer              :: GPTR2(:,:), GPTR2Out(:,:), GPTR3(:,:,:)
    integer                    :: counts(5)
    real, allocatable          :: ple3D(:,:,:), pl3D(:,:,:), levs(:)
    real                       :: pow
    real                       :: scale
    logical :: VERB = .false.

!                              ---

    ASSERT_(MCFIO%CREATED)

    scale = MCFIO%vscale
    vvar  = mCFIO%VVAR
 
    if(present(BUNDLE)) then
       UBUNDLE = BUNDLE
    else
       UBUNDLE = mcfio%BUNDLE
    endif

    if(present(CLOCK)) then
       UCLOCK = CLOCK
    else
       UCLOCK = mcfio%CLOCK
    endif

    if ( present(VERBOSE) ) VERB = VERBOSE

    call ESMF_BundleGet      (Ubundle, FieldCount=NumVars,  RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_BundleGetField (UBUNDLE,  1, FIELD,            RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldGet       (FIELD,    GRID=ESMFGRID,       RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_GridGet(ESMFGRID, dimCount=gridRank, rc=STATUS)
    VERIFY_(STATUS)

    if (gridRank == 3) then
       call ESMF_GridGet(ESMFGRID, horzRelLoc=ESMF_CELL_CENTER, &
            vertRelLoc=ESMF_CELL_CELL, &
            globalCellCountPerDim=COUNTS,  RC=STATUS)
       VERIFY_(STATUS)
       IM0 = COUNTS(1)
       JM0 = COUNTS(2)
       call ESMF_GridGetDELocalInfo(ESMFGRID, &
            horzRelLoc=ESMF_CELL_CENTER, &
            vertRelLoc=ESMF_CELL_CELL, &
            localCellCountPerDim=COUNTS,RC=STATUS)
       VERIFY_(STATUS)
       LM0 = COUNTS(3)
    else
       call ESMF_GridGet(ESMFGRID, horzRelLoc=ESMF_CELL_CENTER, &
            globalCellCountPerDim=COUNTS, RC=STATUS)
       VERIFY_(STATUS)
       IM0 = COUNTS(1)
       JM0 = COUNTS(2)
       LM0 = 1
    endif

    call ESMF_GridGet        (ESMFGRID, DELAYOUT=LAYOUT,     RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ClockGet       (UCLOCK,    CurrTime =TIME,      RC=STATUS)
    VERIFY_(STATUS)
    TIME = TIME - MCFIO%OFFSET
    call ESMF_TimeGet        (TIME,     timeString=DATE,     RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_CFIOGet        (MCFIO%CFIO,     grid=GRID,           RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_CFIOGridGet    (GRID,     IM=IM, JM=JM, KM=LM, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_CFIOGridGet    (GRID,     lev=lev,             RC=STATUS)
    VERIFY_(STATUS)

    allocate(GPtr2Out(IM,JM),stat=STATUS)
    VERIFY_(STATUS)
    allocate(GPtr2(IM0,JM0),stat=STATUS)
    VERIFY_(STATUS)
    allocate(GPtr3(IM,JM,LM),stat=STATUS)
    VERIFY_(STATUS)

    if(Vvar/="") then

       nullify (ptr3)
       call ESMF_BundleGetDataPointer(Ubundle, FieldName=Vvar, DataPointer=Ptr3,  RC=STATUS)
       VERIFY_(STATUS)

       allocate(ple3D(size(Ptr3,1),size(Ptr3,2),size(Ptr3,3)  ),stat=status)
       allocate( pl3D(size(Ptr3,1),size(Ptr3,2),size(Ptr3,3)-1),stat=status)
       VERIFY_(STATUS)
       allocate(Levs(LM),stat=STATUS)
       VERIFY_(STATUS)

       if(vvar(1:3)=='log') then
          ple3D = log(Ptr3)
           pl3D = log( 0.5*(Ptr3(:,:,1:)+Ptr3(:,:,0:ubound(Ptr3,3)-1)) )
           levs = log(lev*scale)
       elseif(vvar(1:3)=='pow') then
          ple3D = Ptr3**pow
           pl3D =    ( 0.5*(Ptr3(:,:,1:)+Ptr3(:,:,0:ubound(Ptr3,3)-1)) )**pow
           levs = (lev*scale)**pow
       else
          ple3D = Ptr3
           pl3D =    ( 0.5*(Ptr3(:,:,1:)+Ptr3(:,:,0:ubound(Ptr3,3)-1)) )
           levs = lev*scale
       end if

    endif

    if (MAPL_AM_I_ROOT(LAYOUT)) then
       if (IM /= IM0 .or. JM /= JM0) then
          if (IM <  IM0 .or. JM < JM0) then
             call MAPL_HorzTransformCreate (Trans, im0, jm0, im, jm, &
                  XYOFFSET=MCFIO%XYOFFSET, rc=STATUS)
             VERIFY_(STATUS)
          else
             call MAPL_HorzTransformCreate (Trans, im0, jm0, im, jm, &
                  XYOFFSET=MCFIO%XYOFFSET, order=1, rc=STATUS)
             VERIFY_(STATUS)
          end if
       end if
    end if

    do L=1,NumVars

       call ESMF_BundleGetField     (UBUNDLE, L, FIELD,           RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet           (FIELD, NAME=VarName,         RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet (FIELD, ARRAY=array,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ArrayGet(array, rank=arrayRank, rc=status)
       VERIFY_(STATUS)

       if (arrayRank == 2) then
          call ESMF_ArrayGetData(Array, PTR2, RC=STATUS)
          VERIFY_(STATUS)
          

          call ArrayGather(PTR2, GPTR2, ESMFGRID, RC=STATUS)
          VERIFY_(STATUS)


          if (MAPL_AM_I_ROOT(LAYOUT)) then
             if (IM /= IM0 .or. JM /= JM0) then
!!$                if (IM <  IM0 .or. JM < JM0) then
                   call MAPL_HorzTransformRun(Trans, Gptr2, Gptr2out, MAPL_undef, rc=STATUS)
                   VERIFY_(STATUS)
!!$                else
!!$                   call hinterp ( Gptr2,im0,jm0, Gptr2out,im,jm,1,MAPL_Undef )
!!$                end if
             else
                Gptr2Out = Gptr2
             end if
             if(present(NBITS)) then
               call ESMF_CFIODownBit ( GPTR2Out, GPTR2Out, NBITS, undef=MAPL_undef, rc=STATUS )
               VERIFY_(STATUS)
               if ( VERB ) print *, Iam // ': Writing '//trim(VARNAME)// &
                           ' at ' // trim(date) // ' with ', NBITS, ' bits'
             else
               if ( VERB ) print *, Iam // ': Writing '//trim(VARNAME)// &
                           ' at ' // trim(date)
             end if
             call ESMF_CFIOVarWrite(MCFIO%CFIO, trim(VARNAME), GPTR2Out, timeString=DATE, RC=STATUS)
             VERIFY_(STATUS)
          end if
       else
          call ESMF_ArrayGetData(array, PTR3, RC=STATUS)
          VERIFY_(STATUS)

          K1 = LBOUND(PTR3,3)-1

          ALLOCATE( ptr2( size(ptr3,1),size(ptr3,2) ),stat=status )
          VERIFY_(STATUS)
          do K=1,LM          
             if(Vvar/="") then
                if( size(ptr3,3).eq.size(ple3d,3) ) then
                    call VertInterp(PTR2,PTR3,ple3d,LEVS(K),rc=status)
                    VERIFY_(STATUS)
                else
                    call VertInterp(PTR2,PTR3, pl3d,LEVS(K),ple3d(:,:,ubound(ple3d,3)),rc=status)
                    VERIFY_(STATUS)
                endif
                call ArrayGather(PTR2, Gptr2, ESMFGRID, RC=STATUS)
                VERIFY_(STATUS)
             else
                call ArrayGather(PTR3(:,:,nint(lev(k))+K1), Gptr2, ESMFGRID, RC=STATUS)
                VERIFY_(STATUS)
             endif
             if (MAPL_AM_I_ROOT(LAYOUT)) then
                if (IM /= IM0 .or. JM /= JM0) then
!!$                   if (IM <  IM0 .or. JM < JM0) then
                      call MAPL_HorzTransformRun(Trans, Gptr2, Gptr2out, MAPL_undef, rc=STATUS)
                      VERIFY_(STATUS)
!!$                   else
!!$                      call hinterp ( Gptr2,im0,jm0, Gptr2out,im,jm,1,MAPL_Undef )
!!$                   end if
                   Gptr3(:,:,K) = Gptr2Out
                else
                   Gptr3(:,:,K) = Gptr2
                end if
             endif
          end do
          DEALLOCATE( ptr2 )
          
          if (MAPL_AM_I_ROOT(LAYOUT)) then
            if(present(NBITS)) then
               call ESMF_CFIODownBit ( GPTR3, GPTR3, NBITS, undef=MAPL_undef, rc=STATUS )
               VERIFY_(STATUS)
               if ( VERB ) print *, Iam // ': Writing '//trim(VARNAME)// &
                           ' at ' // trim(date) // ' with ', NBITS, ' bits'
            else
               if ( VERB ) print *, Iam // ': Writing '//trim(VARNAME)// &
                           ' at ' // trim(date)
            end if
            call ESMF_CFIOVarWrite(MCFIO%CFIO, trim(VARNAME), GPTR3, timeString=DATE, RC=STATUS)
             VERIFY_(STATUS)
          end if
       end if

    end do

    if (MAPL_AM_I_ROOT(LAYOUT)) then
       if (IM /= IM0 .or. JM /= JM0) then
!!$          if (IM <  IM0 .or. JM < JM0) then
             call MAPL_HorzTransformDestroy(Trans,rc=STATUS)
             VERIFY_(STATUS)
!!$          end if
       end if
    end if

! Clean up
!---------

    if(allocated(ple3D))deallocate(ple3D)
    if(allocated( pl3D))deallocate( pl3D)
    if(allocated(Levs ))deallocate(Levs )

    deallocate(GPtr3)
    deallocate(GPtr2)
    deallocate(GPtr2Out)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_CFIOWriteBundle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: MAPL_CFIOWriteState

! !INTERFACE:

  subroutine MAPL_CFIOWriteState ( MCFIO, CLOCK, State, &
                                   VERBOSE, NBITS, RC   )

! !ARGUMENTS:

    type(MAPL_CFIO  ),           intent(INOUT) :: MCFIO
    type(ESMF_State),  intent(IN   ) :: STATE
    type(ESMF_CLOCK), optional,  intent(IN   ) :: CLOCK
    integer, optional,           intent(  OUT) :: RC
    logical, optional,           intent(  IN)  :: VERBOSE
    integer, optional,           intent(  IN)  :: NBITS

! !DESCRIPTION: Serializes an ESMF state into a Bundle and writes it
!               to a file.
!
!EOP

    character(len=*), parameter  :: Iam="MAPL_CFIOWriteState"
    integer                      :: STATUS

! Locals

    type(ESMF_Bundle) :: tBUNDLE

! Get the appropriate bundle
!---------------------------

!!ALT    if(present(STATE)) then
       tBUNDLE = ESMF_BundleCreate ( name=Iam, rc=STATUS )
       VERIFY_(STATUS)
       call ESMFL_BundleAddState_ ( tBUNDLE, STATE, rc=STATUS, VALIDATE=.true. )
       VERIFY_(STATUS)
!!ALT    else
!!ALT       tBUNDLE = MCFIO%BUNDLE
!!ALT    end if

!   Write the Bundle
!   ----------------

    call MAPL_CFIOWriteBundle ( MCFIO, CLOCK=CLOCK, BUNDLE=tBUNDLE, &
                                VERBOSE=VERBOSE, NBITS=NBITS, RC=STATUS   )
    VERIFY_(STATUS)

!!ALT    if(present(STATE)) then
       call ESMF_BundleDestroy ( tBUNDLE, rc=STATUS )
       VERIFY_(STATUS)
!!ALT    endif

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)

 end subroutine MAPL_CFIOWriteState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: MAPL_CFIOReadBundle

! !DESCRIPTION: Reads a bundle

! !INTERFACE:

  subroutine MAPL_CFIOReadBundle( FILETMPL, TIME, BUNDLE, NOREAD, RC, &
                                  VERBOSE, FORCE_REGRID, ONLY_VARS,   &
                                  TIME_IS_CYCLIC, TIME_INTERP )

! !ARGUMENTS:

    character(len=*),            intent(IN   ) :: FILETMPL
    type(ESMF_TIME),             intent(IN   ) :: TIME
    type(ESMF_BUNDLE),           intent(INOUT) :: BUNDLE
    logical, optional,           intent(IN   ) :: NOREAD
    integer, optional,           intent(  OUT) :: RC
    logical, optional,           intent(  IN)  :: VERBOSE
    logical, optional,           intent(IN)    :: FORCE_REGRID ! obsolete
    logical, optional,           intent(IN)    :: TIME_IS_CYCLIC
    logical, optional,           intent(IN)    :: TIME_INTERP
    character(len=*), optional,  intent(IN   ) :: ONLY_VARS ! comma separated,
                                                            ! no spaces


!EOP

    character(len=*), parameter  :: Iam="MAPL_CFIOReadBundle"
    integer                      :: STATUS

! Locals


    type(ESMF_CFIO)              :: CFIO
    type(ESMF_CFIOGrid), pointer :: CFIOGRID
    type(ESMF_GRID)              :: ESMFGRID
    type(ESMF_FIELD)             :: FIELD
    type(ESMF_ARRAY)             :: ARRAY
    type(ESMF_DELAYOUT)          :: LAYOUT
    type(ESMF_FieldDataMap)      :: DATAMAP2
    type(ESMF_FieldDataMap)      :: DATAMAP3

    type(ESMF_CFIOVarInfo), pointer :: VARS(:)

    type(MAPL_HORZTRANSFORM)          :: Trans
    integer                      :: IM,  JM,  LM
    integer                      :: IM0, JM0
    integer                      :: L1, L, K
    integer                      :: NumVars, nVars
    integer                      :: counts(5)
    integer                      :: dims(3)
    integer                      :: arrayRank

    logical                      :: IamRoot, twoD

    real, pointer                ::  PTR2      (:,:),  PTR3      (:,:,:)
    real, pointer                :: GPTR2bundle(:,:), GPTR3bundle(:,:,:)
    real, pointer                :: GPTR2file  (:,:), GPTR3file  (:,:,:)

    character(len=ESMF_MAXSTR)   :: NAME
    character(len=ESMF_MAXSTR)   :: DATE
    character(len=ESMF_MAXSTR)   :: BundleVARNAME
    character(len=ESMF_MAXSTR)   :: CFIOVARNAME

    real, pointer :: LONSfile(:),   LATSfile(:)
    real, pointer :: LONSbundle(:), LATSbundle(:)

    character(len=ESMF_MAXSTR) :: FILENAME
    integer :: nymd, nhms
    logical :: timeInterp=.false., VERB = .false., change_resolution

!                              ---
    
    if ( present(VERBOSE) )     VERB = VERBOSE
    if ( present(TIME_INTERP) ) timeInterp = TIME_INTERP

! Create a CFIO object names after the bundle
!--------------------------------------------
    call ESMF_BundleGet     (Bundle,   name=NAME,                         RC=STATUS)
    VERIFY_(STATUS)
    cfio =  ESMF_CFIOCreate (          cfioObjName=trim(Name),            RC=STATUS)
    VERIFY_(STATUS)

! Transform ESMF time to string for use in CFIO
!----------------------------------------------
    call ESMF_TimeGet       (TIME,     timeString=DATE,                   RC=STATUS)
    VERIFY_(STATUS)

    call strToInt(DATE, nymd, nhms)
    call ESMF_CFIOstrTemplate ( filename, filetmpl, 'GRADS', &
                                nymd=nymd, nhms=nhms, stat=status )
    VERIFY_(STATUS)
    call WRITE_PARALLEL("CFIO: Reading " // trim(filename))
                                                                                              
! Set its filename and open it for reading
!-----------------------------------------
    call ESMF_CFIOSet       (CFIO, fName=trim(fileName), RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_CFIOFileOpen  (CFIO, FMODE=1, cyclic=TIME_IS_CYCLIC, RC=STATUS)
    VERIFY_(STATUS)

! Get info from the bundle
!-------------------------

    call ESMF_BundleGet     (Bundle,   Grid=ESMFGRID, FieldCount=NUMVARS, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_GridGet       (ESMFGRID, horzRelLoc=ESMF_CELL_CENTER, &
                             vertRelloc = ESMF_CELL_CELL,           &
                             globalCellCountPerDim=COUNTS,                RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_GridGet       (ESMFGRID, delayout=LAYOUT,                   RC=STATUS)
    VERIFY_(STATUS)

    IamRoot = MAPL_AM_I_ROOT(LAYOUT)

! Get info from the CFIO object
!------------------------------
    call ESMF_CFIOGet       (CFIO,     grid=CFIOGRID,                     RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_CFIOGridGet   (CFIOGRID, IM=IM, JM=JM, KM=LM,               RC=STATUS)
    VERIFY_(STATUS)

    allocate(LONSfile(IM), LATSfile(JM), stat=status )
    VERIFY_(STATUS)
    call ESMF_CFIOGridGet    (CFIOGRID, LON=LONSFILE, LAT=LATSFILE, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_GridGetDELocalInfo(ESMFGRID, &
         horzRelLoc=ESMF_CELL_CENTER, &
         vertRelLoc=ESMF_CELL_CELL, &
         localCellCountPerDim=DIMS, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_CFIOGet (CFIO,varObjs=VARS, nVars=nVars, RC=STATUS)
    VERIFY_(STATUS)

! Assert compatibility of file and bundle
!----------------------------------------
    ASSERT_( LM==0 .or. counts(3) == 0 .or. LM==counts(3) )

! Get lat/lons of input bundle
! ----------------------------
    call GridGetLatLons_ ( ESMFGRID, LONSbundle, LATSbundle )

! If the bundle is empty, read entire varlist from file
!------------------------------------------------------

    if(NUMVARS==0) then

       NUMVARS = nVars

       call ESMF_FieldDataMapSetDefault(datamap2, 2, rc=status)
       VERIFY_(STATUS)
       call ESMF_FieldDataMapSetDefault(datamap3, 3, rc=status)
       VERIFY_(STATUS)

       L1 = 0
       do L=1,NUMVARS

          call ESMF_CFIOVarInfoGet(VARS(L),vname=CFIOVARNAME, twoDimVar=twoD, RC=STATUS)   
          VERIFY_(STATUS)

          if ( present(ONLY_VARS) ) then
               if ( index(','//trim(ONLY_VARS)  //',', &
                          ','//trim(CFIOVARNAME)//',') < 1 ) cycle 
          endif

          L1 = L1 + 1

          BundleVarName = CFIOVARNAME
          if(twoD) then
            allocate(PTR2(DIMS(1),DIMS(2)),stat=STATUS)
            VERIFY_(STATUS)
            PTR2  = 0.0
            FIELD = ESMF_FieldCreate(grid=ESMFGRID, copyflag=ESMF_DATA_REF,   &
                            fptr=PTR2, name=BundleVARNAME, datamap=datamap2, RC=STATUS)
            VERIFY_(STATUS) 
!ALT: for now we add only HorzOnly (no tiles)
            call ESMF_FieldSetAttribute(FIELD, NAME='DIMS', VALUE=MAPL_DimsHorzOnly, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldSetAttribute(FIELD, NAME='VLOCATION', &
                                        VALUE=MAPL_VLocationNone, RC=STATUS)
            VERIFY_(STATUS) 
          else 
            allocate(PTR3(DIMS(1),DIMS(2),LM),stat=STATUS)
            VERIFY_(STATUS)
            PTR3  = 0.0
            FIELD = ESMF_FieldCreate(grid=ESMFGRID, copyflag=ESMF_DATA_REF,   &
                            fptr=PTR3, name=BundleVARNAME, datamap=datamap3, RC=STATUS)
            VERIFY_(STATUS)
!ALT: for now we add only HorzVert (no tiles)
            call ESMF_FieldSetAttribute(FIELD, NAME='DIMS', VALUE=MAPL_DimsHorzVert, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldSetAttribute(FIELD, NAME='VLOCATION', &
                                        VALUE=MAPL_VLocationCenter, RC=STATUS)
            VERIFY_(STATUS)
          end if
          call ESMF_BundleAddField(BUNDLE,FIELD,                          RC=STATUS)
          VERIFY_(STATUS)

       end do
       NUMVARS = L1  ! could be less than on file if user chooses to

    else
       
       do L=1,NumVars
          call ESMF_BundleGetField (BUNDLE, L, FIELD,                     RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_FieldGet       (FIELD,NAME=BundleVarName,array=ARRAY, RC=STATUS)
          VERIFY_(STATUS)
          do K=1,size(VARS)
             call ESMF_CFIOVarInfoGet(VARS(K),vname=CFIOVARNAME,          RC=STATUS)
             VERIFY_(STATUS)
             if(trim(BUNDLEVARNAME)==trim(CFIOVARNAME)) exit
          end do
!ams      ASSERT_(K<=size(VARS)) ! K is generally not defined at this point!
       end do

    end if

    if(present(NOREAD)) then
       if(NOREAD) then
          RETURN_(ESMF_SUCCESS)
       end if
    end if


! Allocate space for global arrays. Only root will use these
!-----------------------------------------------------------

    IM0 = counts(1)
    JM0 = counts(2)

    allocate(Gptr2bundle(IM0,JM0   ), stat=STATUS)
    VERIFY_(STATUS)
    allocate(Gptr3bundle(IM0,JM0,LM), stat=STATUS)
    VERIFY_(STATUS)
    allocate(Gptr2file  (IM ,JM    ), stat=STATUS)
    VERIFY_(STATUS)
    allocate(Gptr3file  (IM ,JM ,LM), stat=STATUS)
    VERIFY_(STATUS)

    if (IM /= IM0 .or. JM /= JM0)  then
        change_resolution = .true.
    else                              
        change_resolution = .false.
    end if

!    if ( present(FORCE_REGRID) ) then
!         change_resolution = change_resolution .OR. FORCE_REGRID
!    endif

    if (MAPL_AM_I_ROOT(LAYOUT)) then
       if ( change_resolution ) then
          if (IM0 <  IM .or. JM0 < JM) then
             call MAPL_HorzTransformCreate (Trans, im, jm, im0, jm0, rc=STATUS)
             VERIFY_(STATUS)
          end if
       end if
    end if

! Read each variable
!-------------------

    do L=1,NumVars

       call ESMF_BundleGetField (BUNDLE, L, FIELD,                       RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet       (FIELD, NAME=BundleVarName, array=ARRAY, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ArrayGet       (array, rank=arrayRank,                  RC=STATUS)
 
       VERIFY_(STATUS)

       if ( VERB .and. IamRoot ) &
            print *, Iam // ': Reading '//trim(BundleVARNAME)// &
                            ' at ' // trim(date)
       select case (arrayRank)

       case (2)

          call ESMF_ArrayGetData(Array, PTR2, RC=STATUS)
          VERIFY_(STATUS)

          if (IamRoot) then
             if ( timeInterp ) then
                call ESMF_CFIOVarReadT(CFIO, trim(BundleVARNAME), GPTR2file, &
                                       timeString=DATE, RC=STATUS)
             else
                call ESMF_CFIOVarRead (CFIO, trim(BundleVARNAME), GPTR2file, &
                                       timeString=DATE, RC=STATUS)
             endif
             VERIFY_(STATUS)
             if (change_resolution) then  
                if (IM0 <  IM .OR. JM0 < JM) then
                   if ( VERB ) print *, Iam // ': Binning... '
                   call MAPL_HorzTransformRun(Trans, Gptr2file, Gptr2bundle, MAPL_undef, rc=STATUS )
                   VERIFY_(STATUS)
                else
                   if ( VERB ) print *, Iam // ': Interpolating... '
                   call hinterp ( Gptr2file,im,jm, Gptr2bundle,im0,jm0,1,MAPL_Undef )
                end if ! coarsening
             else
                Gptr2bundle = Gptr2file
             end if ! change resolution
             if ( abs(LONSbundle(1)-LONSfile(1)+180.) <= 10.*tiny(1.) ) then
                 if ( VERB ) print *, Iam // &
                      ': shifting input longitudes by 180 degrees'
                 call shift180Lon2D_ ( Gptr2bundle, im0, jm0 )
             end if
          end if

          call ArrayScatter(PTR2, GPTR2bundle, ESMFGRID, RC=STATUS)
          VERIFY_(STATUS)

       case(3)

          call ESMF_ArrayGetData(array, PTR3, RC=STATUS)
          VERIFY_(STATUS)

!         TO DO: Read one level at a time and scatter
!         -------------------------------------------
          if (IamRoot) then
             if ( timeInterp ) then
                call ESMF_CFIOVarReadT(CFIO, trim(BundleVARNAME), GPTR3file, &
                                             timeString=DATE, RC=STATUS)
             else
                call ESMF_CFIOVarRead (CFIO, trim(BundleVARNAME), GPTR3file, &
                                             timeString=DATE, RC=STATUS)
             end if
             VERIFY_(STATUS)
             if (change_resolution) then
                if (IM0 <  IM .or. JM0 < JM) then
                    if ( VERB ) print *, Iam // ': Binning... '
                    do k = 1, LM
                       call MAPL_HorzTransformRun(Trans, Gptr3file(:,:,k), &
                                                Gptr3bundle(:,:,k), MAPL_undef, rc=STATUS)
                       VERIFY_(STATUS)
                    end do
                else
                   if ( VERB ) print *, Iam // ': Interpolating... '
                   call hinterp ( Gptr3file,im,jm, Gptr3bundle,im0,jm0,LM,MAPL_Undef )
                end if ! coarsening
             else
                Gptr3bundle = Gptr3file
             end if
             if ( abs(LONSbundle(1)-LONSfile(1)+180.) <= 10.*tiny(1.) ) then
                 if ( VERB ) print *, Iam // &
                      ': shifting input longitudes by 180 degrees'
                  do k = 1, LM
                     call shift180Lon2D_ ( Gptr3bundle(:,:,K), im0, jm0 )
                  end do 
             end if

          end if ! I am root

          L1 = LBOUND(PTR3,3)-1

          do K=1,LM
             call ArrayScatter(PTR3(:,:,K+L1), Gptr3bundle(:,:,K), ESMFGRID, RC=STATUS)
             VERIFY_(STATUS)
          end do

       case default

       end select

    end do

    if (MAPL_AM_I_ROOT(LAYOUT)) then
       if ( change_resolution ) then
          if (IM <  IM0 .or. JM < JM0) then
             call MAPL_HorzTransformDestroy(Trans,rc=STATUS)
             VERIFY_(STATUS)
          end if
       end if
    end if

    deallocate(GPtr2bundle)
    deallocate(GPtr3bundle)
    deallocate(GPtr2file  )
    deallocate(GPtr3file  )
    deallocate(LONSfile,LATSfile)
    deallocate(LONSbundle,LATSbundle)

    RETURN_(ESMF_SUCCESS)

CONTAINS

    subroutine shift180Lon2D_ ( c, im, jm )
    integer, intent(in) :: im, jm
    real, intent(inout) :: c(im,jm)
    real :: cj(im)
    integer :: m(4), n(4), imh, j
    imh = nint(im/2.)
    m = (/ 1,      imh, 1+imh,    im   /)
    n = (/ 1,   im-imh, 1+im-imh, im   /)
    do j = 1, jm
       cj(n(1):n(2)) = c(m(3):m(4),j)
       cj(n(3):n(4)) = c(m(1):m(2),j)
       c(:,j) = cj
    end do
    return
    end subroutine shift180Lon2D_

  end subroutine MAPL_CFIOReadBundle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: MAPL_CFIOReadState

! !DESCRIPTION: Serializes an ESMF state into a Bundle and reads its content
!               from a file.

! !INTERFACE:

  subroutine MAPL_CFIOReadState ( FILETMPL, TIME, STATE, NOREAD, RC, &
                                  VERBOSE, FORCE_REGRID, ONLY_VARS,   &
                                  TIME_IS_CYCLIC, TIME_INTERP )

! !ARGUMENTS:

    character(len=*),            intent(IN   ) :: FILETMPL
    type(ESMF_TIME),             intent(IN   ) :: TIME
    type(ESMF_STATE),            intent(INOUT) :: STATE
    logical, optional,           intent(IN   ) :: NOREAD
    integer, optional,           intent(  OUT) :: RC
    logical, optional,           intent(  IN)  :: VERBOSE
    logical, optional,           intent(IN)    :: FORCE_REGRID ! obsolete
    logical, optional,           intent(IN)    :: TIME_IS_CYCLIC
    logical, optional,           intent(IN)    :: TIME_INTERP
    character(len=*), optional,  intent(IN   ) :: ONLY_VARS ! comma separated,
                                                            ! no spaces

!
!EOP

    character(len=*), parameter  :: Iam="MAPL_CFIOWriteState"
    integer                      :: STATUS

! Locals

    type(ESMF_Bundle) :: tBUNDLE

!                          ----

!   Create an empty bundle
!   ----------------------
    tBUNDLE = ESMF_BundleCreate ( name=Iam, rc=STATUS )
    VERIFY_(STATUS)
    
!   Serialize the state
!   -------------------
    call ESMFL_BundleAddState_ ( tBUNDLE, STATE, rc=STATUS, VALIDATE=.true. )
    VERIFY_(STATUS)

!   Read the Bundle
!   ---------------
    call MAPL_CFIOReadBundle( FILETMPL, TIME, tBUNDLE,         &
                              NOREAD = NOREAD,                 &
                              VERBOSE = VERBOSE,               &
                              FORCE_REGRID=FORCE_REGRID,       &
                              ONLY_VARS = ONLY_VARS,           &
                              TIME_IS_CYCLIC = TIME_IS_CYCLIC, &
                              TIME_INTERP = TIME_INTERP,       &
                              RC = STATUS )

    VERIFY_(STATUS)

!   All done
!   --------
    call ESMF_BundleDestroy ( tBUNDLE, rc=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

 end subroutine MAPL_CFIOReadState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: MAPL_CFIOReadField

! !DESCRIPTION: Reads a variable into

! !INTERFACE:

  subroutine MAPL_CFIOReadField     ( VARN, FILETMPL, TIME, GRID, FIELD, RC, &
                                      VERBOSE, FORCE_REGRID, TIME_IS_CYCLIC, &
                                      TIME_INTERP )

! !ARGUMENTS:

    character(len=*),            intent(IN   ) :: VARN       ! Variable name
    character(len=*),            intent(IN   ) :: FILETMPL   ! File name
    type(ESMF_TIME),             intent(IN   ) :: TIME
    type(ESMF_GRID),             intent(IN   ) :: GRID
    type(ESMF_FIELD),            intent(IN)    :: FIELD
    integer, optional,           intent(  OUT) :: RC
    logical, optional,           intent(  IN)  :: VERBOSE
    logical, optional,           intent(IN)    :: FORCE_REGRID
    logical, optional,           intent(IN)    :: TIME_IS_CYCLIC
    logical, optional,           intent(IN)    :: TIME_INTERP

!EOP

    character(len=*), parameter  :: Iam="MAPL_CFIOReadField"
    integer                      :: STATUS

! Locals

    type(ESMF_BUNDLE)  :: BUNDLE
 
!   Create a temporary empty bundle
!   -------------------------------
    BUNDLE =  ESMF_BundleCreate ( grid=GRID, name=Iam, rc=STATUS )
    VERIFY_(STATUS)

!   Add the input field to the bundle
!   ---------------------------------
    call ESMF_BundleAddField ( BUNDLE, FIELD, rc=STATUS )
    VERIFY_(STATUS)

!   Now, we read the variable into the bundle, which in turn will put
!    the data inside the input array
!   -----------------------------------------------------------------
    call MAPL_CFIOReadBundle( FILETMPL, TIME, BUNDLE,                    &
                              VERBOSE=VERBOSE,                           &
                              FORCE_REGRID=FORCE_REGRID,                 &
                              ONLY_VARS = trim(varn),                    &
                              TIME_IS_CYCLIC=TIME_IS_CYCLIC,             &
                              TIME_INTERP=TIME_INTERP, RC=STATUS         )
    VERIFY_(STATUS)    


!   Destroy temporary bundle; field data will be preserved
!   ------------------------------------------------------
    call ESMF_BundleDestroy ( BUNDLE, rc=STATUS )
    VERIFY_(STATUS)    

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_CFIOReadField

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: MAPL_CFIOReadArray3D

! !DESCRIPTION: Reads a variable into FORTRAN 3D array

! !INTERFACE:

  subroutine MAPL_CFIOReadArray3D ( VARN, FILETMPL, TIME, GRID, fARRAY, RC, &
                                    VERBOSE, FORCE_REGRID, TIME_IS_CYCLIC,  &
                                    TIME_INTERP )

! !ARGUMENTS:

    character(len=*),            intent(IN   ) :: VARN       ! Variable name
    character(len=*),            intent(IN   ) :: FILETMPL   ! File name
    type(ESMF_TIME),             intent(IN   ) :: TIME
    type(ESMF_GRID),             intent(IN   ) :: GRID
    real, pointer                              :: fARRAY(:,:,:)
    integer, optional,           intent(  OUT) :: RC
    logical, optional,           intent(  IN)  :: VERBOSE
    logical, optional,           intent(IN)    :: FORCE_REGRID
    logical, optional,           intent(IN)    :: TIME_IS_CYCLIC
    logical, optional,           intent(IN)    :: TIME_INTERP

!EOP

    character(len=*), parameter  :: Iam="MAPL_CFIOReadArray3D"
    integer                      :: STATUS

    type(ESMF_Field)             :: FIELD
    type(ESMF_FieldDataMap)      :: DATAMAP

    real    :: const = 0.0
    integer :: ios, k

!                            ----

!   Special case: when filename is "/dev/null" it is assumed the user 
!   wants to set the variable to a constant
!   -----------------------------------------------------------------
    if ( FILETMPL(1:9) == '/dev/null' ) then    
         ios = -1
         k = index(FILETMPL,':')
         if ( k > 9 ) read(FILETMPL(k+1:),*,iostat=ios) const
         if ( ios /= 0 ) const = 0.0
         if ( MAPL_am_I_root() ) &
            print *, Iam // ': setting variable ' // trim(varn) // &
                            ' to constant = ', const
         RETURN_(ESMF_SUCCESS)
    end if

!   Create Field with input array
!   -----------------------------
    call ESMF_FieldDataMapSetDefault(datamap, 3, rc=status)
    VERIFY_(STATUS)
    FIELD = ESMF_FieldCreate(grid=GRID, copyflag=ESMF_DATA_REF,   &
            fptr=fARRAY, name=trim(varn), datamap=datamap, RC=STATUS)
    VERIFY_(STATUS)

   
!   Read array data from file
!   -------------------------
    call MAPL_CFIOReadField ( VARN, FILETMPL, TIME, GRID, FIELD,          &
                              VERBOSE=VERBOSE, FORCE_REGRID=FORCE_REGRID, &
                              TIME_IS_CYCLIC=TIME_IS_CYCLIC,              &
                              TIME_INTERP=TIME_INTERP, RC=STATUS         )
    VERIFY_(STATUS)

!   Destroy the ESMF array (data will be preserved since we own it)
!   --------------------------------------------------------------
    call ESMF_FieldDestroy ( FIELD, STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_CFIOReadArray3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: MAPL_CFIOReadArray2D

! !DESCRIPTION: Reads a variable into FORTRAN 2D array using an auxilliary ESMF grid.

! !INTERFACE:

  subroutine MAPL_CFIOReadArray2D ( VARN, FILETMPL, TIME, GRID, fARRAY, RC, &
                                    VERBOSE, FORCE_REGRID, TIME_IS_CYCLIC,  &
                                    TIME_INTERP )

! !ARGUMENTS:

    character(len=*),            intent(IN)  :: VARN       ! Variable name
    character(len=*),            intent(IN)  :: FILETMPL   ! File name
    type(ESMF_TIME),             intent(IN)  :: TIME
    type(ESMF_GRID),             intent(IN)  :: GRID
    real, pointer                            :: fARRAY(:,:)
    integer, optional,           intent(OUT) :: RC
    logical, optional,           intent(IN)  :: VERBOSE
    logical, optional,           intent(IN)  :: FORCE_REGRID
    logical, optional,           intent(IN)  :: TIME_IS_CYCLIC
    logical, optional,           intent(IN)  :: TIME_INTERP



!EOP

    character(len=*), parameter  :: Iam="MAPL_CFIOReadArray2D"
    integer                      :: STATUS

    type(ESMF_Field)             :: FIELD
    type(ESMF_FieldDataMap)      :: DATAMAP

    real    :: const = 0.0
    integer :: ios, k

!                            ----


!   Special case: when filename is "/dev/null" it is assumed the user 
!   wants to set the variable to a constant
!   -----------------------------------------------------------------
    if ( FILETMPL(1:9) == '/dev/null' ) then    
         ios = -1
         k = index(FILETMPL,':')
         if ( k > 9 ) read(FILETMPL(k+1:),*,iostat=ios) const
         if ( ios /= 0 ) const = 0.0
         if ( MAPL_am_I_root() ) &
            print *, Iam // ': setting variable ' // trim(varn) // &
                            ' to constant = ', const
         RETURN_(ESMF_SUCCESS)
    end if

!   Create Field with input array
!   -----------------------------
    call ESMF_FieldDataMapSetDefault(datamap, 2, rc=status)
    VERIFY_(STATUS)
    FIELD = ESMF_FieldCreate(grid=GRID, copyflag=ESMF_DATA_REF,   &
            fptr=fARRAY, name=trim(varn), datamap=datamap, RC=STATUS)
    VERIFY_(STATUS)
   
!   Read array data from file
!   -------------------------
    call MAPL_CFIOReadField ( VARN, FILETMPL, TIME, GRID, FIELD,          &
                              VERBOSE=VERBOSE, FORCE_REGRID=FORCE_REGRID, &
                              TIME_INTERP=TIME_INTERP,                    &
                              TIME_IS_CYCLIC=TIME_IS_CYCLIC, RC=STATUS    )
    VERIFY_(STATUS)

!   Destroy the ESMF array (data will be preserved since we own it)
!   --------------------------------------------------------------
    call ESMF_FieldDestroy ( FIELD, STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

   end subroutine MAPL_CFIOReadArray2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MAPL_CFIODestroy( MCFIO, RC )
  type(MAPL_CFIO),             intent(INOUT) :: MCFIO
  integer, optional,           intent(  OUT) :: RC
  integer :: status
  character(len=*), parameter  :: Iam="MAPL_CFIODestroy"


  if (MAPL_AM_I_ROOT()) then
     call ESMF_CFIOFileClose(MCFIO%CFIO,rc=status)
     if ( present(rc) ) then
        rc = status
     else
        VERIFY_(STATUS)
     endif
  end if
  call ESMF_CFIODestroy(MCFIO%CFIO,rc=status)
  if ( present(rc) ) then
       rc = status
  else
       VERIFY_(STATUS)
  endif
  end subroutine MAPL_CFIODestroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is a candidate for ESMFL, here for dependency reasons
!  

  subroutine GridGetLatLons_ ( grid, lons, lats )

    implicit NONE
    type(ESMF_Grid) :: grid

    real, pointer   :: lons(:), lats(:)

!                     ---

    character(len=*), parameter :: Iam = 'GridGetLatLons'

    type(ESMF_Array)       :: eARRAY(2)

    real(KIND=8), pointer  :: R8D2(:,:)
    real, pointer          :: lons2d(:,:), lats2d(:,:)
    real, pointer          :: LONSLocal(:,:), LATSlocal(:,:)
    integer                :: IM_WORLD, JM_WORLD, dims(3), STATUS, RC

!                          ----

!      Get world dimensions
!      --------------------
       call ESMF_GridGet ( grid, horzRelloc=ESMF_CELL_CENTER, &
                           vertRelLoc=ESMF_CELL_CENTER, &
                           globalCellCountPerDim=DIMS, RC=STATUS)
       VERIFY_(STATUS)

       IM_WORLD = dims(1)
       JM_WORLD = dims(2)

!      Allocate memory for output if necessary
!      ---------------------------------------
       if ( .not. associated(lons) ) then
            allocate(lons(IM_WORLD), stat=STATUS)
       else
            if(size(LONS,1) /= IM_WORLD) STATUS = 1
       end if
       VERIFY_(status)
       if ( .not. associated(lats) ) then
            allocate(lats(JM_WORLD), stat=STATUS)
       else
            if(size(LATS,1) /= JM_WORLD) STATUS = 1
       end if
       VERIFY_(status)

!      Retrieve the ESMF array with coordinates 
!      ----------------------------------------
       call ESMF_GridGetCoord ( grid, horzRelLoc =ESMF_CELL_CENTER, &
                                centerCoord=eARRAY, RC=STATUS ) 
       VERIFY_(STATUS) 

!      Local work space
!      ----------------
       allocate(LONS2d(IM_WORLD,JM_WORLD), LATS2d(IM_WORLD,JM_WORLD), &
                STAT=status)             
       VERIFY_(status)

!      Get the local longitudes and gather them into a global array
!      ------------------------------------------------------------
       call ESMF_ArrayGetData(EARRAY(1), R8D2, RC=STATUS)
       VERIFY_(STATUS)

       allocate(LONSLOCAL(size(R8D2,1),size(R8D2,2)), STAT=status)             
       VERIFY_(status)

       LONSLOCAL = R8D2*(180/MAPL_PI)

       call ArrayGather(LONSLOCAL, LONS2D, GRID, RC=STATUS)

!      Get the local longitudes and gather them into a global array
!      ------------------------------------------------------------
       call ESMF_ArrayGetData(eARRAY(2), R8D2, RC=STATUS)
       VERIFY_(STATUS)

       allocate(LATSLOCAL(size(R8D2,1),size(R8D2,2)), STAT=status)             
       VERIFY_(status)

       LATSlocal = R8D2*(180/MAPL_PI)

       call ArrayGather(LATSLOCAL, LATS2D, GRID, RC=STATUS)
       VERIFY_(STATUS)

!      Return 1D arrays
!      ----------------
       LONS = LONS2D(:,1)
       LATS = LATS2D(1,:)

       DEALLOCATE(LONSLOCAL, LATSLOCAL, LONS2d, LATS2d )
        
     end subroutine GridGetLatLons_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is another candidate for ESMFL, here for dependency reasons
!  
!
! !IROUTINE: ESMFL_BundleAddState - Addes contents of State to Bundle

! !INTERFACE:

    RECURSIVE subroutine ESMFL_BundleAddState_ ( BUNDLE, STATE, rc, &
                                                 GRID, VALIDATE )

! !ARGUMENTS:

    implicit NONE
    type(ESMF_Bundle), intent(inout)         :: BUNDLE
    type(ESMF_State),  intent(in)            :: STATE
    integer, optional                        :: rc
    type(ESMF_State),  optional, intent(in)  :: GRID
    logical, optional, intent(in)            :: VALIDATE


! !DESCRIPTION: Extracts fields from an ESMF State and adds them to a
!               ESMF Bundle. In essesence, it serializes an ESMF state
!  in a flat Bundle. The BUNDLE must have been created prior to calling 
!  this routine.
!
!EOP

    character(len=*), parameter          :: Iam="ESMFL_StateSerialize"
    integer                              :: STATUS

    type(ESMF_State)                     :: tSTATE
    type(ESMF_Bundle)                    :: tBUNDLE
    type(ESMF_Field)                     :: tFIELD
    type(ESMF_Grid)                      :: tGRID

    integer                              :: I, J
    integer                              :: ItemCount, FieldCount
    type (ESMF_StateItemType), pointer   :: ItemTypes(:)
    character(len=ESMF_MAXSTR ), pointer :: ItemNames(:), FieldNames(:)
    logical                              :: needGrid = .true.
    logical                              :: validate_ = .false.

!                           ---

    
    if ( present(validate) ) validate_ = validate

!   Query state for number of items and allocate space for them
!   -----------------------------------------------------------
    call ESMF_StateGet(STATE,ItemCount=ItemCount,RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(ItemCount>0)
    allocate ( ItemNames(ItemCount), stat=STATUS)
    VERIFY_(STATUS)
    allocate ( ItemTypes(ItemCount), stat=STATUS)
    VERIFY_(STATUS)

!   Next, retrieve the names and types of each item in the state
!   ------------------------------------------------------------
    call ESMF_StateGet ( STATE,      ItemNameList = ItemNames, &
                                StateItemtypeList = ItemTypes, &
                                rc=STATUS)
    VERIFY_(STATUS)

!   Loop over each item on STATE
!   ----------------------------
    do I = 1, ItemCount

!         Got a field
!         -----------
          if (ItemTypes(I) == ESMF_StateItem_Field) THEN

             call ESMF_StateGetField ( STATE, ItemNames(i), tFIELD, rc=status)
             VERIFY_(STATUS)
             call AddThisField_()

!         Got a Bundle
!         ------------
          else if (ItemTypes(I) == ESMF_StateItem_Bundle) then

             call ESMF_StateGetBundle(STATE, ItemNames(i), tBUNDLE, rc=STATUS)
             VERIFY_(STATUS)
             call ESMF_BundleGet ( tBUNDLE, FieldCount = FieldCount, rc=STATUS)
             VERIFY_(STATUS)
             ASSERT_(FieldCount>0)
             do j = 1, FieldCount
                call ESMF_BundleGetField ( tBUNDLE, j, tFIELD, rc=STATUS)
                VERIFY_(STATUS)
                call AddThisField_()
             end do

!         Got another state
!         -----------------
          else if (ItemTypes(I) == ESMF_StateItem_State) then

             call ESMF_StateGetState(STATE, ItemNames(i), tSTATE, rc=STATUS)
             VERIFY_(STATUS)
             call ESMFL_BundleAddState_ ( BUNDLE, tSTATE, rc=STATUS )
             VERIFY_(STATUS)

!         What is this?
!         ------------
          else

             STATUS = -1
             VERIFY_(STATUS)
           
          end IF

    end do

!   Make sure the Bundle is not empty
!   ---------------------------------
    call ESMF_BundleGet ( BUNDLE, FieldCount = FieldCount, rc=STATUS)
    VERIFY_(STATUS)
    ASSERT_(FieldCount>0)

!   Set the grid
!   ------------
    if ( present(GRID) ) then
       call ESMF_BundleSetGrid ( BUNDLE, tGRID, rc=STATUS )
       VERIFY_(STATUS)
    else if ( needGrid ) then
       STATUS = -1              ! could not find a Grid
       VERIFY_(STATUS)
    else
       call ESMF_BundleSetGrid ( BUNDLE, tGRID, rc=STATUS )
       VERIFY_(STATUS)
    end if

!   Make sure field names are unique
!   --------------------------------
    allocate ( FieldNames(FieldCount), stat=STATUS )
    VERIFY_(STATUS) 
    call ESMF_BundleGetFieldNames ( BUNDLE, FieldNames, rc=STATUS )
    VERIFY_(STATUS) 

!   Make sure field names are unique
!   --------------------------------
    if ( validate_ ) then
       do j = 1, FieldCount
          do i = j+1, FieldCount
             if ( trim(FieldNames(i)) == trim(FieldNames(j)) ) then
                STATUS = -1              ! same name
                VERIFY_(STATUS)
             end if
          end do
       end do
    end if

!   All done
!   --------
    deallocate(ItemNames)
    deallocate(ItemTypes)
    deallocate(FieldNames)

    RETURN_(ESMF_SUCCESS)

CONTAINS

    subroutine AddThisField_()
      call ESMF_BundleAddField ( BUNDLE, tFIELD, rc=STATUS )
      VERIFY_(STATUS)
      if ( needGrid ) then
         call ESMF_FieldGet ( tFIELD, grid=tGRID, rc=STATUS )
         VERIFY_(STATUS)
         needGrid = .false.
      end if
    end subroutine AddThisField_

  end subroutine ESMFL_BundleAddState_

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine VertInterp(v2,v3,pl,pp,ps,rc)

    real,              intent(OUT) :: v2(:,:)
    real,              intent(IN ) :: v3(:,:,:)
    real,     target,  intent(IN ) :: pl(:,:,:)
    real,              intent(IN ) :: pp
    real,    optional, intent(IN ) :: ps(:,:)
    integer, optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer km, K
    logical flip
    real    ppx
    real, pointer   :: plx(:,:,:)

    integer        :: status
    character*(10) :: Iam='VertInterp'

    km   = size(pl,3)

    flip = pl(1,1,km) < pl(1,1,km-1)

    if(flip) then
       allocate(plx(size(pl,1),size(pl,2),size(pl,3)),stat=status)
       VERIFY_(STATUS)
       plx = -pl
       ppx = -pp
    else
       plx => pl
       ppx = pp
    end if

    v2   = MAPL_UNDEF

       pb   = plx(:,:,km)
       do k=km-1,1,-1
          pt = plx(:,:,k)
          if(all(pb<ppx)) exit
          where(ppx>pt .and. ppx<=pb)
             al = (pb-ppx)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do

! Extend Lowest Level Value to the Surface
! ----------------------------------------
    if( present(ps) ) then
        where( (ppx>plx(:,:,km).and.ppx<=ps) )
                v2 = v3(:,:,km)
        end where
    end if

    if(flip) then
       deallocate(plx,stat=status)
       VERIFY_(STATUS)
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp


end module MAPL_CFIOMod
