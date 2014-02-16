!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSurfaceTemperature_mod
!
! !INTERFACE:
!
      module GmiSurfaceTemperature_mod
!
! !USES:
      use ESMF_Mod
      use m_netcdf_io_open      , only : Ncop_Rd
      use m_netcdf_io_read      , only : Ncrd_3d
      use m_netcdf_io_close     , only : Nccl
      use m_netcdf_io_get_dimlen, only : Ncget_Unlim_Dimlen
      use GmiMetFieldsAttribute_mod, only : getMetFieldsAttribute
      use GmiTimeControl_mod        , only : t_GmiClock, Get_curGmiDate
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use GmiGrid_mod, only : t_gmiGrid
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_j2
      use GmiGrid_mod, only : Get_i1_gl, Get_ju1_gl, Get_i2_gl, Get_j2_gl
      use GmiDiagnosticsMethod_mod     , only : t_Diagnostics, Get_pr_diag
      use GmiMetFieldsControl_mod, only : t_metFields, Get_metNumMEGAN, &
     &       Get_mdt, Get_met_num_infiles, Get_met_infile_names, &
     &       Get_surfTemp15DayAvg, Set_surfTemp15DayAvg
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: updateSurfaceTempMEGAN

#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains routines to read in the surface temperatures, compute
! and update (at the beginning of each day) the 15-day average of surface
! temperature for the MEGAN emissions.
!
! !DEFINED PARAMETERS:
      integer, parameter  :: NUM_DAYS       = 15
!
      character(len=MAX_LENGTH_VAR_NAME), save :: surfTempName         ! name of the surface temperature
                                                      ! in the metFields file
      real*8 , pointer, save :: surfTemp15Days(:,:,:) ! Average daily surface
                                                      ! temperature over the
                                                      ! past 15 days.

      integer, save :: numRecsPerDay                  ! number of records per day
      integer, save :: numRecsPerFile                 ! number of records per file
      integer, save :: numFilesPerDay                 ! number of files   per day
      integer, save :: currentFileNumber              ! current file number to be read
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateSurfaceTempMEGAN
!
! !INTERFACE
!
      subroutine updateSurfaceTempMEGAN (metFields, gmiGrid, gmiDomain, &
     &                 gmiClock, Diagnostics)
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_GmiClock   ), intent(in) :: gmiClock
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_metFields), intent(inOut) :: metFields
!
! !DESCRIPTION:
!  Updates the 15-day average surface temperature for the MEGAN emission.
!     
! !LOCAL VARIABLES:
      logical, save        :: first = .true.
      logical              :: pr_diag
      integer              :: procID, nymd, met_num_infiles
      integer              :: i1, i2, ju1, j2, i1_gl, ju1_gl
      integer, save        :: oldDay = -999
      integer              :: curDay, metNumMEGAN
      real*8               :: mdt
      real*8 , allocatable :: surfTemp(:,:)
      real*8 , allocatable :: surfTemp15DayAvg(:,:)
      character (len=MAX_LENGTH_FILE_NAME), pointer :: met_infile_names(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID      (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'updateSurfaceTempMEGAN called by ', procID

      call Get_i1     (gmiGrid, i1    )   
      call Get_i2     (gmiGrid, i2    )
      call Get_ju1    (gmiGrid, ju1   )
      call Get_j2     (gmiGrid, j2    )
      call Get_i1_gl  (gmiGrid, i1_gl )
      call Get_ju1_gl (gmiGrid, ju1_gl)

      allocate(surfTemp15DayAvg(i1:i2,ju1:j2))
      call Get_surfTemp15DayAvg(metFields, surfTemp15DayAvg)

      call Get_curGmiDate (gmiClock, nymd)

      call Get_met_num_infiles  (metFields, met_num_infiles)

      allocate(met_infile_names(met_num_infiles))
      call Get_met_infile_names (metFields, met_infile_names)

      if (first) then
         call Get_mdt         (metFields, mdt)
         call Get_metNumMEGAN (metFields, metNumMEGAN)

         call initSurfTemp15DayAvg(surfTemp15DayAvg, met_infile_names, &
     &            met_num_infiles, mdt, metNumMEGAN, pr_diag, procID, &
     &            i1, i2, ju1, j2, i1_gl, ju1_gl)

         oldDay = Mod (nymd, 100)
         first = .false.
      else
         ! Determine the current day
         curDay = Mod (nymd, 100)

         ! Update the surface temperature array if a new day
         if (oldDay /= curDay) then
            allocate(surfTemp(i1:i2, ju1:j2))
            oldDay = curDay
            call readSurfaceTemperature(surfTemp, met_infile_names, &
     &               met_num_infiles, pr_diag, procID, &
     &               i1, i2, ju1, j2, i1_gl, ju1_gl)
            call updateSurfTemp15DayAvg(surfTemp15DayAvg, surfTemp, &
     &                    i1, i2, ju1, j2)
            deallocate(surfTemp)
         end if
      end if

      call Set_surfTemp15DayAvg(metFields, surfTemp15DayAvg)
      deallocate(surfTemp15DayAvg)

      end subroutine updateSurfaceTempMEGAN
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeSurfTemp15Days
!
! !INTERFACE:
!
      subroutine initSurfTemp15DayAvg(surfTemp15DayAvg, metFieldsList, numFiles, &
     &                                mdt, metNumMEGAN, pr_diag, procID, &
     &                                i1, i2, ju1, j2, i1_gl, ju1_gl)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: i1_gl, ju1_gl
      integer, intent(in) :: numFiles, procID
      integer, intent(in) :: metNumMEGAN       ! file number in the list where 
                                               ! GMI model will start from
      real*8 , intent(in) :: mdt               ! Time increment for reading new 
                                               ! metFields data
      character (len=*), intent(in) :: metFieldsList(numFiles) ! avail. metFields files.
!
! !OUTPUT PARAMETERS:
      ! 24h average surface temperature over past 15 days
      real*8 , intent(out) :: surfTemp15DayAvg(i1:i2, ju1:j2)
!
! !DESCRIPTION:
! Computes the daily averages of surface temperature for the 15
! days preceding the start of the model integration. And it then computes
! the average surface tempratures over the 15 day period.
!
! !LOCAL VARIABLES:
      integer :: it, iday, i, j, ifreq, fid, numFiles15Days
      real*8 , allocatable :: surfTemp(:,:)
      character(len=50)    :: metdata_name
      character(len=MAX_LENGTH_MET_NAME)    :: metdata_name_org
      character(len=MAX_LENGTH_MET_NAME)    :: metdata_name_model
      character(len=MAX_LENGTH_MET_NAME)    :: metdata_name_dims
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'initSurfTemp15DayAvg called by ', procID

      allocate(surfTemp(i1:i2, ju1:j2))

      allocate(surfTemp15Days(i1:i2, ju1:j2,NUM_DAYS))

      surfTemp15Days(:,:,:) = 0.0d0

      iday = 0

      ! Determine the number of records per day
      ifreq         = mdt/3600
      numRecsPerDay = 24/ifreq

      ! Open the first metFields file
      call Ncop_Rd (fid, metFieldsList(metNumMEGAN)) 

      ! Get the number of records in the file
      call Ncget_Unlim_Dimlen (fid, numRecsPerFile)

      ! Get the name of surface temperature from the metField file
      call getMetFieldsAttribute (fid, metdata_name, metdata_name_org,  &
     &               metdata_name_model, metdata_name_dims, pr_diag, procID)

      if (metdata_name_model(1:5) == 'GEOS3') then
         surfTempName = 't10m'
      elseif (metdata_name_model(1:5) == 'GEOS4') then
         surfTempName = 't2m'
      elseif (metdata_name_model(1:5) == 'GEOS5') then
         surfTempName = 'T2M'
      else
         surfTempName = 'ts'
      end if

      call Nccl(fid)

      ! Determine the number of files per day
      numFilesPerDay = numRecsPerDay / numRecsPerFile

      ! Determine the number of files over 15 days
      numFiles15Days = NUM_DAYS*numFilesPerDay

      ! File to be read in
      currentFileNumber = metNumMEGAN - numFiles15Days + 1

      if (currentFileNumber < 1) then
         write(6,*) 'Insufficient number of metFields files to initialize'
         write(6,*) 'MEGAN emissions. You need to provide at least' 
         write(6,*)  numFiles15Days, ' metFields files'
         STOP
      end if

      ! Read in the 15 days worth of surface temperature.
      do it = 1, NUM_DAYS
         call readSurfaceTemperature(surfTemp, metFieldsList, numFiles, &
     &                       pr_diag, procID, i1, i2, ju1, j2, i1_gl, ju1_gl)

         surfTemp15Days(:,:,it) = surfTemp(:,:)
      end do

      ! Compute the average over 15 days
      do j = ju1, j2
         do i = i1, i2
            surfTemp15DayAvg(i,j) = sum( surfTemp15Days(i,j,:) ) / dble(NUM_DAYS)
         end do
      end do

      deallocate(surfTemp)

      return

      end subroutine initSurfTemp15DayAvg
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readSurfaceTemperature
!
! !INTERFACE:
!
      subroutine readSurfaceTemperature(surfTemp, metFieldsList, numFiles, &
     &                      pr_diag, procID, i1, i2, ju1, j2, i1_gl, ju1_gl)
!
      implicit none

#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
      logical         , intent(in) :: pr_diag
      integer         , intent(in) :: procID, numFiles
      integer         , intent(in) :: i1, i2, ju1, j2, i1_gl, ju1_gl
      character (len=*), intent(in) :: metFieldsList(numFiles)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: surfTemp(i1:i2, ju1:j2)
! 
! !DESCRIPTION:
! Reads in all the records of the surface temperature and computes
! the daily average.
!
! !LOCAL VARIABLES:
      integer              :: fid, it, ifile
      integer              :: il, ij, inb, jnb
      integer              :: count3d(3), start3d(3)
      real*8               :: counter
      real*8, allocatable  :: temp2d(:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) "readSurfaceTemperature called by", procID

      allocate(temp2d(i2-i1+1, j2-ju1+1))
      temp2d  (:,:) = 0.0d0

      surfTemp(:,:) = 0.0d0

      start3d(:) = (/ i1 - i1_gl + 1, ju1 - ju1_gl + 1, 1/)
      count3d(:) = (/ i2-i1+1, j2-ju1+1, 1 /)

      do ifile = 1, numFilesPerDay
         ! Open the file
         call Ncop_Rd (fid, metFieldsList(currentFileNumber))

         ! Loop over the number of records
         do it = 1, numRecsPerFile
            start3d(3) = it
            call Ncrd_3d (temp2d, fid, surfTempName, start3d, count3d)
         
            do ij = ju1, j2
               jnb = ij - ju1 + 1
               do il = i1, i2
                  inb = il - i1 + 1
                  surfTemp(il,ij) = surfTemp(il,ij) + temp2d(inb,jnb)
               end do
            end do
         end do

         ! Close the file
         call Nccl(fid)

         ! Next file number to be read in
         currentFileNumber = currentFileNumber + 1
      end do

      counter = numRecsPerFile * numFilesPerDay
      
      ! Compute the daily average surface temperature
      surfTemp(:,:) = surfTemp(:,:) / counter

      deallocate(temp2d)

      return

      end subroutine readSurfaceTemperature      
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateSurfTemp15DayAvg
!
! !INTERFACE:
!
      subroutine updateSurfTemp15DayAvg(surfTemp15DayAvg, surfTemp, &
     &                                  i1, i2, ju1, j2)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2
      real*8 , intent(in) :: surfTemp(i1:i2, ju1:j2)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: surfTemp15DayAvg(i1:i2, ju1:j2)
!
! !DESCRIPTION:
!  MEGAN requires the surface temperature averaged over the past 15 days,
!  which is stored in the array {\em surfTemp15DayAvg}. Temperature is in Celsius.
!  This subroutine should be called at the beginning of each day.
!  It does the following:
!  \begin{enumerate}
!   \item Push the daily average surface temperature values through {\em surfTemp15Days}, 
!         throwing out the oldest and putting the newest (the {\em surfTemp15DayAvg} 
!         average) in the last spot 
!   \item Get {\em surfTemp15DayAvg} by averaging {\em surfTemp15Days} over the 15 day period.
!  \end{enumerate}
!
! !LOCAL VARIABLES:
      integer   :: it, i, j
!EOP
!------------------------------------------------------------------------------
!BOC

      do it = NUM_DAYS, 2, -1
         surfTemp15Days(:,:,it) = surfTemp15Days(:,:,it-1)
      end do

      surfTemp15Days(:,:,1) = surfTemp(:,:)

      do j = ju1, j2
         do i = i1, i2
            surfTemp15DayAvg(i,j) = sum( surfTemp15Days(i,j,:) ) / dble(NUM_DAYS)
         end do
      end do

      return

      end subroutine updateSurfTemp15DayAvg
!EOC
!------------------------------------------------------------------------------
      end module GmiSurfaceTemperature_mod
