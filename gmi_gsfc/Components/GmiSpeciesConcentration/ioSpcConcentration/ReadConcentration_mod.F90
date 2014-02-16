     module ReadConcentration_mod

      use GmiPrintError_mod     , only : GmiPrintError
      use m_netcdf_io_close     , only : Nccl
      use m_netcdf_io_open      , only : Ncop_Rd
      use m_netcdf_io_read      , only : Ncrd_4d, Ncrd_5d
      use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
      use GmiCheckRange_mod     , only : checkRange3d
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

     implicit none

     private
     public  :: setAllConcentration
     public  :: setSliceConcentration, setSliceConcentrationRst
     public  :: getConcentrationFileInfo


     CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   setAllConcentration
!
! DESCRIPTION
!   This routine sets the concentrations for all species.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------

      subroutine setAllConcentration &
                    (const_infile_name, const_var_name, const_init_val, concentration, &
                     const_opt, num_const_inrecs, coscen, num_molefrac, &
                     pr_diag, procID, do_full_chem, &
                     i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
                     numSpecies)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: const_opt, procID
      integer, intent(in) :: numSpecies, num_const_inrecs, num_molefrac
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl
      logical, intent(in) :: do_full_chem, pr_diag
      character (len=*)   :: const_var_name
      character (len=*)   :: const_infile_name
      real*8              :: const_init_val(numSpecies)
      real*8              :: coscen(ju1_gl:j2_gl)
      type (t_GmiArrayBundle) :: concentration(numSpecies)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) then
        Write (6,*) 'setAllConcentration called by ', procID
      end if

      ! Set the concentration one species at the time
      do ic = 1, numSpecies
         call setSliceConcentration &
                 (ic, const_infile_name, const_var_name, const_init_val, &
                  concentration, const_opt, num_const_inrecs, coscen, &
                  num_molefrac,  pr_diag, procID, do_full_chem, &
                  i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
                  numSpecies)
      end do

      return

      end subroutine setAllConcentration
!EOC
!-----------------------------------------------------------------------------
!BOP
      subroutine setSliceConcentration &
                    (ic, const_infile_name, const_var_name, const_init_val, concentration, &
                     const_opt, num_const_inrecs, coscen, num_molefrac, &
                     pr_diag, procID, do_full_chem, &
                     i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
                     numSpecies)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: ic
      integer, intent(in) :: const_opt, procID
      integer, intent(in) :: numSpecies, num_const_inrecs, num_molefrac
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl
      logical, intent(in) :: do_full_chem, pr_diag
      character (len=*)   :: const_var_name
      character (len=*)   :: const_infile_name
      real*8              :: const_init_val(numSpecies)
      real*8              :: coscen(ju1_gl:j2_gl)
      type (t_GmiArrayBundle) :: concentration(numSpecies)
!
! !LOCAL VARIABLES:
      integer :: ivert, ilat, ilong

#ifdef MICRO_AEROSOL
      real*8, parameter :: threshold = 1.0d20
#elif GOCARTaerosol
      real*8, parameter :: threshold = 1.0d20
#else
      real*8, parameter :: threshold = 1.0d20
#endif
!EOP
!-----------------------------------------------------------------------------
!BOC

      if (pr_diag) then
        Write (6,*) 'setSliceConcentration called by ', procID
      end if
      
      ilong = i2 - i1 + 1
      ilat  = j2 - ju1 + 1
      ivert = k2 - k1 + 1

      if (const_opt == 1) then
         concentration(ic)%pArray3D(:,:,:) = const_init_val(ic)
      else if (const_opt == 2) then
         if (ic <= num_const_inrecs) then
            call  readSliceConcentration (const_infile_name, const_var_name, &
     &                    ic, concentration(ic)%pArray3D(:,:,:), &
     &                    i1, i2, ju1, j2, k1, k2, &
     &                    i1_gl, ju1_gl, ilong, ilat, ivert)

            if (.not. (do_full_chem .and. (ic > num_molefrac))) then
               call checkRange3d ('const_infile', procID, i1, i2, ju1, j2, k1, k2,  &
                              concentration(ic)%pArray3D(:,:,:), 0.0d0, threshold)
            end if
         else
            Write (6,*) 'WARNING:  No data for species #', ic, ' ; will set to 0.'
            concentration(ic)%pArray3D(:,:,:) = 1.0d-30
         end if

      else if (const_opt >= 3) then
         call calcSliceConcentration (const_opt, ic, coscen, &
                  concentration(ic)%pArray3D(:,:,:), &
                  i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)
      end if

      return

      end subroutine setSliceConcentration
!EOC
!------------------------------------------------------------------------------------------
!BOP
      subroutine setSliceConcentrationRst &
                    (ic, restart_infile_name, const_var_name, restart_inrec, &
                     concentration, tracer_opt, num_molefrac, &
                     pr_diag, procID, do_full_chem, &
                     i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
                     numSpecies)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: ic
      integer, intent(in) :: tracer_opt, procID
      integer, intent(in) :: numSpecies, num_molefrac, restart_inrec
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl
      logical, intent(in) :: do_full_chem, pr_diag
      character (len=*)   :: const_var_name
      character (len=*)   :: restart_infile_name
      real*8              :: const_init_val(numSpecies)
      type (t_GmiArrayBundle) :: concentration(numSpecies)
!
! !LOCAL VARIABLES:
      integer :: ivert, ilat, ilong
      real*8 :: min_range, max_range

!EOP
!-----------------------------------------------------------------------------
!BOC

      if (pr_diag) Write (6,*) 'setSliceConcentrationRst called by ', procID

      ilong = i2 - i1 + 1
      ilat  = j2 - ju1 + 1
      ivert = k2 - k1 + 1

#ifdef MICRO_AEROSOL
       min_range = 0.0d0
       max_range = 1.0d20
#else
       if(tracer_opt .gt. 0) then
         min_range = -1.0d-10
         max_range = 1.0d20
       else
         min_range = -1.0d-10
         max_range = 1.0d0
       endif
#endif

       call readSliceConcentration (restart_infile_name, const_var_name, &
     &                    ic, concentration(ic)%pArray3D(:,:,:), &
     &                    i1, i2, ju1, j2, k1, k2, &
     &                    i1_gl, ju1_gl, ilong, ilat, ivert, restart_inrec)


!
       if (.not. (do_full_chem .and. (ic > num_molefrac))) then
          call checkRange3d ('const_rst', procID, i1, i2, ju1, j2, k1, k2,  &
                    concentration(ic)%pArray3D(:,:,:), min_range, max_range)
       end if


      return

      end subroutine setSliceConcentrationRst
!EOC
!-----------------------------------------------------------------------------
!
! ROUTINE
!   readSliceConcentration
!
! DESCRIPTION
!   This routine reads in a single set of constituent data (i.e., one
!   species' worth).  Recall that all NetCDF arrays always start at
!   subscript 1, and this is not always the case for Gmimod dimensions.
!
! ARGUMENTS
!   const_infile_name : const input file name
!   const_var_name    : NetCDF const variable name
!   ic     : const species to read
!   const1 : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine readSliceConcentration  &
     &  (infile_name, const_var_name, ic, const1, &
     &   i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, ilong, ilat, ivert, restart_inrec)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl
      integer, intent(in) :: ilong, ilat, ivert
      character (len=*) :: infile_name
      character (len=*) :: const_var_name
      integer :: ic
      integer, optional, intent(in) :: restart_inrec
      real*8  :: const1(i1:i2, ju1:j2, k1:k2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik
      integer :: inc, jnc
      integer :: ncid

      integer :: cnt4d (4), cnt5d (5)
      integer :: strt4d(4), strt5d(5)

      real*8, allocatable  :: const4d_nc(:,:,:,:)
      real*8, allocatable  :: const5d_nc(:,:,:,:,:)

!     ----------------
!     Begin execution.
!     ----------------

      if (.not. present(restart_inrec)) then
         ! Reading from a regular species concentration file
         allocate(const4d_nc(ilong, ilat, ivert, 1))

         strt4d(1) = i1  -  i1_gl + 1
         strt4d(2) = ju1 - ju1_gl + 1
         strt4d(3) = k1
         strt4d(4) = ic

         cnt4d (:) = (/ ilong, ilat, ivert, 1 /)

         call Ncop_Rd (ncid, infile_name)

         call Ncrd_4d (const4d_nc, ncid, const_var_name, strt4d, cnt4d)

         call Nccl (ncid)

         do ik = k1, k2
            do ij = ju1, j2
               jnc = ij - ju1 + 1
               do il = i1, i2
                  inc = il - i1 + 1
                  const1(il,ij,ik) = const4d_nc(inc,jnc,ik,1)
               end do
            end do
         end do

         deallocate(const4d_nc)
      else


         ! Reading from the restart file

         allocate(const5d_nc(ilong, ilat, ivert, 1, 1))

         strt5d(1) = i1  -  i1_gl + 1
         strt5d(2) = ju1 - ju1_gl + 1
         strt5d(3) = k1
         strt5d(4) = ic
         strt5d(5) = restart_inrec

         cnt5d (:) = (/ ilong, ilat, ivert, 1, 1 /)

         ncid = 0

         call Ncop_Rd (ncid, infile_name)


         call Ncrd_5d (const5d_nc, ncid, const_var_name, strt5d, cnt5d)
         

         call Nccl (ncid)


         do ik = k1, k2
           do ij = ju1, j2
             jnc = ij - ju1 + 1
             do il = i1, i2
               inc = il - i1 + 1

               const1(il,ij,ik) = const5d_nc(inc,jnc,ik,1,1)

             end do
           end do
         end do

         deallocate(const5d_nc)

      end if


      return

      end subroutine readSliceConcentration

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcSliceConcentration
!
! DESCRIPTION
!   This routine sets a simple constituent pattern.  It covers options
!   3-8 for const_opt =>
!
!     3:  solid body rotation test field
!     4:  dummy test pattern with linear slope in each dimension
!     5:  exponential in vertical (decays with height)
!     6:  sin in latitude (largest at equator)
!     7:  linear vertical gradient
!     8:  (sin in latitude (largest at equator) + vert. gradient) * species
!
! ARGUMENTS
!   const_opt : constituent input option
!   ic        : const_species to set
!   coscen    : cosine of latitude of zone centers = cos(dlatr)
!   const1    : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine calcSliceConcentration  &
                    (const_opt, ic, coscen, const1, &
                     i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer :: const_opt
      integer :: ic
      real*8  :: coscen(ju1_gl:j2_gl)
      real*8  :: const1(i1:i2, ju1:j2, k1:k2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg
      integer            :: i2d2, il, ij, ik
      real*8             :: d2, rcap, rij, rik, rj2_gl, rr
      real*8             :: dlong(i1_gl:i2_gl)
      real*8             :: fd(i1_gl:i2_gl, ju1_gl:j2_gl)

!     ----------------
!     Begin execution.
!     ----------------

!     =======================
      select case (const_opt)
!     =======================

!       ========
        case (3)  ! solid body rotation
!       ========

          fd(:,:) = 0.0d0


          d2 = 2.0d0 * GMI_PI / i2_gl

          do il = i1_gl, i2_gl
            dlong(il) = (il - 1) * d2
          end do


          i2d2 = i2_gl / 2
          rcap = RADEAR / 3.0d0


          do ij = ju1_gl, j2_gl
            do il = i2d2 + 1, i2_gl

              rr =  &
     &          RADEAR *  &
     &          (Acos (Cos (dlong(il) - 1.5d0 * GMI_PI) * coscen(ij)))

              if (rr < rcap) then

                fd(il,ij) =  &
     &            500.0d0 * (1.0d0 + Cos (GMI_PI * rr / rcap))

              end if

            end do
          end do


          do ik = k1, k2

            const1(:,:,ik) = fd(:,:)

          end do

!       ========
        case (4)  ! dummy test pattern with linear slope in each dimension
!       ========

          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2

                const1(il,ij,ik) = il + ij + ik + ic

              end do
            end do
          end do

!       ========
        case (5)  ! exponential in vertical (decays with height)
!       ========

          do ik = k1, k2

            rik = ik

            const1(:,:,ik) = Exp (k2 - rik) * 1.0d-20

          end do

!       ========
        case (6)  ! sin in latitude (largest at equator)
!       ========

          rj2_gl = j2_gl

          do ij = ju1, j2

            rij = ij

            const1(:,ij,:) =  &
     &        Max (Sin ((rij / rj2_gl) * GMI_PI), 0.0d0)

          end do

!       ========
        case (7)  ! linear vertical gradient
!       ========

          do ik = k1, k2

            const1(:,:,ik) = ik

          end do

!       ========
        case (8)  ! sin in latitude (largest at equator) + vertical gradient
!       ========

          rj2_gl = j2_gl

          do ik = k1, k2
            do ij = ju1, j2

              rij = ij

              const1(:,ij,ik) =  &
     &          k2 * Max (Sin ((rij / rj2_gl) * GMI_PI), 0.0d0) + ik

            end do
          end do

          const1(:,:,:) = const1(:,:,:) * ic

!       ============
        case default
!       ============

          err_msg = 'Problems in Calc_Slice_Const.'
          call GmiPrintError (err_msg, .true., 1, const_opt, 0, 0, 0.0d0, 0.0d0)

!     ==========
      end select
!     ==========

      return

      end subroutine calcSliceConcentration

!-----------------------------------------------------------------------------
!
! ROUTINE
!   getConcentrationFileInfo
!
! DESCRIPTION
!   This routine also determines the number of const input sets
!   (i.e., the 4th dimension).
!
! ARGUMENTS
!   file_name : name of input file to get const data from
!
!-----------------------------------------------------------------------------

      subroutine getConcentrationFileInfo &
            (file_name, spc_dim_name, numSpecies, num_const_inrecs, pr_diag, procID)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      character (len=*), intent(in) :: file_name
      character (len=*), intent(in) :: spc_dim_name
      integer          , intent(in) :: numSpecies, procID
      logical          , intent(in) :: pr_diag
      integer          , intent(out) :: num_const_inrecs

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ncid

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'getConcentrationFileInfo called by ', procID
      end if

      call Ncop_Rd (ncid, file_name)

      call Ncget_Dimlen (ncid, spc_dim_name, num_const_inrecs)

      call Nccl (ncid)

      if (num_const_inrecs > numSpecies) then

         Write (6,*) 'WARNING:  # const input recs is > # species.'
         Write (6,*) 'Will ignore the extra data.'

         num_const_inrecs = numSpecies

      end if

      return

      end subroutine getConcentrationFileInfo

     end module ReadConcentration_mod
