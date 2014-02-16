!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
      module GmiBuffer_mod
!
! !USES:
      use GmiSub2Glob_mod      , only : subDomain2Global
      use GmiNcOutputSlice_mod , only : prepSlice4d    , prepSlice4d_Nomean
      use GmiNcOutputSlice_mod , only : prepSlice4d_Map, prepSlice4d_Map_Nomean
      use GmiNcOutputSlice_mod , only : prepSlice4d_Emiss
      use GmiNcOutputSlice_mod , only : writeSlice4d, writeSlice5d
      use GmiMessagePassing_mod, only : synchronizeGroup
      use m_netcdf_io_create   , only : Ncdo_Sync
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: bufferOutput4d, bufferOutput4d_Nomean, bufferOutput4d_Emiss
      public  :: bufferOutput4d_Map, bufferOutput4d_Map_Nom
      public  :: bufferOutput5d, bufferOutput5d_Tend
!
! !DESCRIPTION:
! Routines to buffer and write variables with at least 4 dimensions.
! The routines have the following arguments:
! \begin{verbatim}
!   nc_var_name : NetCDF variable name
!   do_mean     : should means values be put in periodic output files?
!   kx1, kx2    : altitude  dimensions
!   rootProc    : root processor id
!   procID      : processor id
!   commuWorld  : world MPI communicator
!   msg_num     : message number
!   ncid_out    : NetCDF output file id
!   num_out     : 4th dimension of norm_array and mean_array
!   rnum_out    : current NetCDF record number to write to
!   num_outrecs : number of species that will be output
!   outrec_map  : mapping of species output number to species input number
!   norm_array  : array of snapshot values
!   mean_array  : array of mean values
!   nc_array    : 3d slice of 4d norm_array or mean_array
! \end{verbatim}
!EOP
!------------------------------------------------------------------------------
       contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput4d
!
! !INTERFACE:
!
      subroutine bufferOutput4d (nc_var_name, do_mean, kx1, kx2, commuWorld, &
     &       msg_num, ncid_out, num_out, rnum_out, norm_array, mean_array,   &
     &       nc_array, map1_u, numDomains, rootProc, procID, &
     &       ix1_gl, ix2_gl, jx1_gl, jx2_gl, ix1, ix2, jx1, jx2)
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=*), intent(in) :: nc_var_name ! netCDF variable name
      integer, intent(in) :: commuWorld            ! message passing communicator
      integer, intent(in) :: ix1, ix2, jx1, jx2    ! local  horizontal dimensions
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl ! global horizontal dimensions
      integer, intent(in) :: kx1, kx2 ! altitude dimension
      integer, intent(in) :: rootProc, procID      ! processor ids
      integer, intent(in) :: numDomains            ! number of domains
      integer, intent(in) :: map1_u(2, 2, numDomains)
      logical, intent(in) :: do_mean  ! should means values be output?
      integer, intent(in) :: msg_num  ! message number
      integer, intent(in) :: ncid_out ! netCDF file id
      integer, intent(in) :: num_out  ! 4th dimension of norm_array and mean_array
      integer, intent(in) :: rnum_out ! current NetCDF record number to write to
      real*8 , intent(in) :: norm_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out) ! array of snapshot values
      real*8 , intent(in) :: mean_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out) ! array of mean values
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(ix1:ix2,jx1:jx2,kx1:kx2)
!
! !DESCRIPTION:
! Buffers and writes out all the 3D slices of a 4d array.
!
! !LOCAL VARIABLES:
      integer :: ic 
      real*8, allocatable :: globArray(:,:,:)
!EOC
!-----------------------------------------------------------------------------
!BOP
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do ic = 1, num_out

         if (procID == rootProc) globArray(:,:,:) = 0.0d0

         call prepSlice4d (do_mean, ix1, ix2, jx1, jx2, kx1, kx2, ic, num_out, &
     &            norm_array, mean_array, nc_array)

         call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

         if (procID == rootProc) then
            call writeSlice4d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2,        &
     &               globArray, ncid_out, nc_var_name, ic, rnum_out)
         end if

         call synchronizeGroup (commuWorld)
      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput4d
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput4d_Nomean
!
! !INTERFACE:
!
      subroutine bufferOutput4d_Nomean  &
     &  (nc_var_name, kx1, kx2, commuWorld, msg_num, ncid_out, num_out,        &
     &   rnum_out, norm_array, nc_array, map1_u, numDomains, rootProc, procID, &
     &   ix1_gl, ix2_gl, jx1_gl, jx2_gl, ix1, ix2, jx1, jx2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ix1, ix2, jx1, jx2
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl
      integer, intent(in) :: rootProc, procID ! processor ids
      character (len=*), intent(in) :: nc_var_name
      integer, intent(in) :: kx1, kx2 ! alitude dimensions
      integer, intent(in) :: commuWorld ! message passing communicator
      integer, intent(in) :: msg_num  ! message number
      integer, intent(in) :: ncid_out ! netCDF file id
      integer, intent(in) :: num_out  ! 4th dimension of norm_array
      integer, intent(in) :: rnum_out ! current NetCDF record number to write to
      real*8 , intent(in) :: norm_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out) ! array of snapshot values
      integer, intent(in) :: numDomains
      integer, intent(in) :: map1_u(2, 2, numDomains)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(ix1:ix2,jx1:jx2,kx1:kx2)
!
! !DESCRIPTION:
!  Buffers and writes out all the 3D slices of a 4d array (no mean values).
!
! !LOCAL VARIABLES:
      integer :: ic
      real*8, allocatable :: globArray(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do ic = 1, num_out

         if (procID == rootProc) globArray(:,:,:) = 0.0d0

         nc_array(:,:,:) = norm_array(:,:,:,ic)

!          call prepSlice4d_Nomean (ix1, ix2, jx1, jx2, kx1, kx2, ic, num_out, &
!     &             norm_array, nc_array)

         call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

        if (procID == rootProc) then
             call writeSlice4d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2, &
     &                 globArray, ncid_out, nc_var_name, ic, rnum_out)
        end if

          call synchronizeGroup (commuWorld)
      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput4d_Nomean
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput4d_Map
!
! !INTERFACE:
!
      subroutine bufferOutput4d_Map  &
     &  (nc_var_name, do_mean, kx1, kx2, commuWorld, msg_num,  &
     &   ncid_out, num_out, rnum_out, num_outrecs, outrec_map,  &
     &   norm_array, mean_array, nc_array, &
     &   map1_u, numDomains, rootProc, procID, &
     &   ix1_gl, ix2_gl, jx1_gl, jx2_gl, &
     &   ix1, ix2, jx1, jx2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ix1, ix2, jx1, jx2
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl
      integer, intent(in) :: rootProc, procID
      character (len=*) :: nc_var_name
      logical, intent(in) :: do_mean
      integer, intent(in) :: kx1, kx2
      integer, intent(in) :: commuWorld
      integer, intent(in) :: msg_num
      integer, intent(in) :: ncid_out
      integer, intent(in) :: num_out
      integer, intent(in) :: rnum_out
      integer, intent(in) :: num_outrecs
      integer, intent(in) :: outrec_map(num_outrecs)
      real*8 , intent(in) :: norm_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out)
      real*8 , intent(in) :: mean_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out)
      integer, intent(in) :: numDomains
      integer, intent(in) :: map1_u(2, 2, numDomains)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(ix1:ix2,jx1:jx2,kx1:kx2)
!
! !DESCRIPTION:
! Buffers and writes out all the 3D slices of a 4d mapped array.
!
! !LOCAL VARIABLES:
      integer :: ic
      real*8, allocatable :: globArray(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do ic = 1, num_outrecs

         if (procID == rootProc) globArray(:,:,:) = 0.0d0

         call prepSlice4d_Map (do_mean, ix1, ix2, jx1, jx2, kx1, kx2, ic, &
     &            num_out, num_outrecs, outrec_map, norm_array, mean_array, &
     &            nc_array)

         call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

         if (procID == rootProc) then
            call writeSlice4d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2, &
     &              globArray, ncid_out, nc_var_name, ic, rnum_out)
         end if

         call synchronizeGroup (commuWorld)
      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput4d_Map
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   bufferOutput4d_Map_Nom
!
! DESCRIPTION
!   This routine buffers and writes out all the 3D slices of a 4d mapped
!   array (no mean values).
!
      subroutine bufferOutput4d_Map_Nom  &
     &  (nc_var_name, kx1, kx2, commuWorld, msg_num, ncid_out,  &
     &   num_out, rnum_out, num_outrecs, outrec_map, norm_array,  &
     &   nc_array, map1_u, numDomains, &
     &   rootProc, procID, ix1_gl, ix2_gl, jx1_gl, jx2_gl, &
     &   ix1, ix2, jx1, jx2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ix1, ix2, jx1, jx2
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl
      integer, intent(in) :: rootProc, procID
      character (len=*), intent(in) :: nc_var_name
      integer, intent(in) :: kx1, kx2
      integer, intent(in) :: commuWorld
      integer, intent(in) :: msg_num
      integer, intent(in) :: ncid_out
      integer, intent(in) :: num_out
      integer, intent(in) :: rnum_out
      integer, intent(in) :: num_outrecs
      integer, intent(in) :: outrec_map(num_outrecs)
      real*8 , intent(in) :: norm_array(ix1:ix2,jx1:jx2,kx1:kx2, num_out)
      integer, intent(in) :: numDomains
      integer, intent(in) :: map1_u(2, 2, numDomains)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(:,:,:)
!
! !LOCAL VARIABLES:
      integer :: ic
      real*8, allocatable :: globArray(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do ic = 1, num_outrecs

         if (procID == rootProc) globArray(:,:,:) = 0.0d0

         call prepSlice4d_Map_Nomean (ix1, ix2, jx1, jx2, kx1, kx2, ic, &
     &             num_out, num_outrecs, outrec_map, norm_array, nc_array)

         call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

         if (procID == rootProc) then
            call writeSlice4d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2,        &
     &               globArray, ncid_out, nc_var_name, ic, rnum_out)
         end if

!         ===============
          call synchronizeGroup (commuWorld)
!         ===============

      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput4d_Map_Nom
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput5d
!
! !INTERFACE:
!
      subroutine bufferOutput5d  &
     &  (nc_var_name, do_mean, kx1, kx2, commuWorld, msg_num,  &
     &   ncid_out, num_out, num_out5, rnum_out, norm_array,  &
     &   mean_array, nc_array, map1_u, numDomains, &
     &   rootProc, procID, ix1_gl, ix2_gl, jx1_gl, jx2_gl, &
     &   ix1, ix2, jx1, jx2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ix1, ix2, jx1, jx2
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl
      integer, intent(in) :: rootProc, procID
      character (len=*), intent(in) :: nc_var_name
      integer, intent(in) :: kx1, kx2
      logical, intent(in) :: do_mean
      integer, intent(in) :: commuWorld
      integer, intent(in) :: msg_num
      integer, intent(in) :: ncid_out
      integer, intent(in) :: num_out
      integer, intent(in) :: num_out5
      integer, intent(in) :: rnum_out
      real*8 , intent(in) :: norm_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out)
      real*8 , intent(in) :: mean_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out)
      integer, intent(in) :: numDomains
      integer, intent(in) :: map1_u(2, 2, numDomains)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(ix1:ix2,jx1:jx2,kx1:kx2)
!
! !DESCRIPTION:
! Buffers and writes out all the 3D slices of a 5d array.
!
! !LOCAL VARIABLES:
      integer :: ic
      real*8, allocatable :: globArray(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do ic = 1, num_out
         if (procID == rootProc) globArray(:,:,:) = 0.0d0

         call prepSlice4d (do_mean, ix1, ix2, jx1, jx2, kx1, kx2, ic, num_out, &
     &            norm_array, mean_array, nc_array)

         call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

         if (procID == rootProc) then
            call writeSlice5d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2,        &
     &               globArray, ncid_out, nc_var_name, ic, num_out5, rnum_out)
         end if

         call synchronizeGroup (commuWorld)

      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput5d
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput5d_Tend
!
! !INTERFACE:
!
      subroutine bufferOutput5d_Tend  &
     &  (nc_var_name, kx1, kx2, commuWorld, msg_num, ncid_out,  &
     &   num_out, num_out5, rnum_out, norm_array, nc_array , &
     &   map1_u, numDomains, rootProc, procID, ix1_gl, ix2_gl, jx1_gl, jx2_gl, &
     &   ix1, ix2, jx1, jx2)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ix1, ix2, jx1, jx2
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl
      integer, intent(in) :: rootProc, procID
      character (len=*) :: nc_var_name
      integer :: kx1, kx2
      integer :: commuWorld
      integer :: msg_num
      integer :: ncid_out
      integer :: num_out
      integer :: num_out5
      integer :: rnum_out
      type (t_GmiArrayBundle) :: norm_array(num_out,num_out5)
      integer, intent(in) :: numDomains
      integer, intent(in) :: map1_u(2, 2, numDomains)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(ix1:ix2,jx1:jx2,kx1:kx2)
!
! !DESCRIPTION:
! Buffers and writes out all the 3D slices of the 5d tend 
! (net cumulative tendencies) array.
!
! !LOCAL VARIABLES:
      integer :: ic, io
      real*8, allocatable :: globArray(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do io = 1, num_out5
         do ic = 1, num_out
            if (procID == rootProc) globArray(:,:,:) = 0.0d0

            nc_array(:,:,:) = norm_array(ic,io)%pArray3D(ix1:ix2,jx1:jx2,kx1:kx2)

            call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &              rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

           if (procID == rootProc) then
              call writeSlice5d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2,     &
     &                  globArray, ncid_out, nc_var_name, ic, io, rnum_out)
           end if

           call synchronizeGroup (commuWorld)
 
         end do
      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput5d_Tend
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput4d_Emiss
!
! !INTERFACE:
!
      subroutine bufferOutput4d_Emiss  &
     &  (nc_var_name, kx1, kx2, commuWorld, msg_num,  &
     &   ncid_out, num_out, rnum_out, num_outrecs, outrec_map,  &
     &   norm_array, nc_array, map1_u, numDomains, &
     &   rootProc, procID, ix1_gl, ix2_gl, jx1_gl, jx2_gl, &
     &   ix1, ix2, jx1, jx2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ix1, ix2, jx1, jx2
      integer, intent(in) :: ix1_gl, ix2_gl, jx1_gl, jx2_gl
      integer, intent(in) :: rootProc, procID
      character (len=*) :: nc_var_name
      integer :: kx1, kx2
      integer :: commuWorld
      integer :: msg_num
      integer :: ncid_out
      integer :: num_out
      integer :: rnum_out
      integer :: num_outrecs
      integer :: outrec_map(num_outrecs)
      real*8  :: norm_array(ix1:ix2,jx1:jx2,kx1:kx2,num_out)
      integer, intent(in) :: numDomains
      integer, intent(in) :: map1_u(2, 2, numDomains)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: nc_array(ix1:ix2,jx1:jx2,kx1:kx2)
!
! !DESCRIPTION:
!  Buffers and writes out all the 3D slices of a 4d mapped array.
!
! !LOCAL VARIABLES:
      integer :: ic
      real*8, allocatable :: globArray(:,:,:)
!   
!EOP
!------------------------------------------------------------------------------
!BOC
      if (procID == rootProc) &
     &   allocate(globArray(ix1_gl:ix2_gl,jx1_gl:jx2_gl,kx1:kx2))

      do ic = 1, num_outrecs
         if (procID == rootProc) globArray(:,:,:) = 0.0d0

         call prepSlice4d_Emiss (ix1, ix2, jx1, jx2, kx1, kx2, ic, num_out, &
     &             num_outrecs, outrec_map, norm_array, nc_array)

         call subDomain2Global (globArray, nc_array, ix1_gl, ix2_gl, jx1_gl, &
     &           jx2_gl, ix1, ix2, jx1, jx2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, msg_num, commuWorld)

         if (procID == rootProc) then
            call writeSlice4d (ix1_gl, ix2_gl, jx1_gl, jx2_gl, kx1, kx2,        &
     &               globArray, ncid_out, nc_var_name, ic, rnum_out)
         end if

         call synchronizeGroup (commuWorld)
      end do

      if (procID == rootProc) then
         deallocate(globArray)
         call Ncdo_Sync (ncid_out)
      end if

      return

      end subroutine bufferOutput4d_Emiss
!EOC
!-----------------------------------------------------------------------------
      end module GmiBuffer_mod
