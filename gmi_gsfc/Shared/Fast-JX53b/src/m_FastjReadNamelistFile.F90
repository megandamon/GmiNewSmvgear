!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: FastjReadNamelistFile
!
! !INTERFACE:
!
   module m_FastjReadNamelistFile
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  FastjReadNamelistFile
!
! !DESCRIPTION:
!  Provides a routine to read the input namelist file.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
 CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FastjReadNamelistFile
!
! !INTERFACE: 
!
      subroutine FastjReadNamelistFile(namelist_file               &
                                      ,LongDim, Latdim, VertDim    &
                                      ,VertDimJval, NumJval        &
                                      ,CrossSection_InfileName     &
                                      ,ScatteringData_InfileName   &
                                      ,PhotRateLabels_InfileName   &
                                      ,T_O3_Climatology_InfileName &
                                      ,CTMdata_InfileName          &
                                      ,SampleRun_InfileName        ) 
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=128), intent (in ) :: namelist_file
!
! !OUTPUT PARAMETERS:
      integer            , intent (out) :: LongDim, Latdim, VertDim
      integer            , intent (out) :: VertDimJval, NumJval
      character (len=128), intent (out) :: CrossSection_InfileName
      character (len=128), intent (out) :: ScatteringData_InfileName
      character (len=128), intent (out) :: PhotRateLabels_InfileName
      character (len=128), intent (out) :: T_O3_Climatology_InfileName
      character (len=128), intent (out) :: CTMdata_InfileName
      character (len=128), intent (out) :: SampleRun_InfileName
!
! !DESCRIPTION:
!  Routines to read in the namelist file.
!
! !LOCAL VARIABLES:
      integer             :: nllun, ierr
      character (len=128) :: err_msg
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
      namelist / NL_InFile / CrossSection_InfileName      &
                            ,ScatteringData_InfileName    &  
                            ,PhotRateLabels_InfileName    &
                            ,T_O3_Climatology_InfileName  &
                            ,CTMdata_InfileName           &
                            ,SampleRun_InfileName

      namelist / NL_Dims  /  LongDim                      &
                            ,LatDim                       &
                            ,VertDim                      &
                            ,VertDimJval                  &
                            ,NumJval

      CrossSection_InfileName     = ' '
      ScatteringData_InfileName   = ' '
      PhotRateLabels_InfileName   = ' '
      T_O3_Climatology_InfileName = ' '
      CTMdata_InfileName          = ' '
      SampleRun_InfileName        = ' '

      nllun = 15
      Open (unit=nllun, file=namelist_file, iostat=ierr)
      if (ierr /= 0) then
         write(6,*) 'Problem opening the file', namelist_file
         stop
      endif

      Read (nllun, nml=NL_InFile, iostat=ierr)
      if (ierr /= 0) then
         err_msg = "Namelist error in section NL_InFile"
         write(6,*) err_msg
      endif

      Read (nllun, nml=NL_Dims, iostat=ierr)
      if (ierr /= 0) then
         err_msg = "Namelist error in section NL_Dims"
         write(6,*) err_msg
      endif

      close(nllun)

      end subroutine FastjReadNamelistFile
!EOC
!-------------------------------------------------------------------------
      end  module m_FastjReadNamelistFile
