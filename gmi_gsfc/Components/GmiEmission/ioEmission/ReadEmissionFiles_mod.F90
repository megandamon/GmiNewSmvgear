module ReadEmissionFiles_mod
!
  use m_netcdf_io_open , only : Ncop_Rd
  use m_netcdf_io_close, only : Nccl
  use m_netcdf_io_read , only : Ncrd_1d, Ncrd_1d_Int, Ncrd_2d, Ncrd_2d_Char, Ncrd_4d, Ncrd_5d
  use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
  use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
!
    implicit none
!
    private
    public  :: ReadEmiss
    public  :: ReadEmissDust
    public  :: ReadEmissAero
    public  :: readGalacticCosmisRayPar
!
#   include "GmiParameters.h"
#   include "setkin_par.h"
!
    contains
!
!-------------------------------------------------------------------------
! This routine reads in a set of emission data.
!-------------------------------------------------------------------------
!
  subroutine ReadEmiss (emissionArray, emiss_infile_name, emiss_var_name, emiss_map,  &
                        emittedSpeciesNames, num_emiss, num_emiss3d, lightning_opt,   &
                        i_no_lgt, curRecord, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl,  &
                        ilong, ilat, ivert, numSpecies, firstReadEmiss)
!
  use GmiStringManipulation_mod, only : convertCharacters2String
!
  implicit none
!
  integer, parameter  :: speciesCharLen = 10
!
  integer, intent(in) :: curRecord   ! current record (month or day) for reading
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: i1_gl, ju1_gl, ilong, ilat, ivert, numSpecies
  integer, intent(in) :: emiss_map(num_emiss)
  integer, intent(in) :: num_emiss, num_emiss3d, lightning_opt, i_no_lgt
  character (len=MAX_LENGTH_VAR_NAME) , intent(in) :: emiss_var_name
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: emiss_infile_name
  character (len=speciesCharLen) , intent(inout) :: emittedSpeciesNames(1:num_emiss)
  logical, intent(inOut) :: firstReadEmiss
  type (t_GmiArrayBundle), intent(inOut) :: emissionArray(num_emiss)
!
  integer :: ic, it, il, ij, ik, inb, jnb
  integer :: ncid_em
!
  integer :: count5d(5), start5d(5), start5d_3d(5)
!
  integer :: count2d(2), start2d(2)
!
  real*8  :: emiss_tmp(ilong, ilat, ivert)
!.sds.. for 2d emissions file
  integer :: count4d(4), start4d(4)
  real*8  :: emiss_tmp_2d(ilong, ilat)
  real*8  :: emiss_tmp_lgt(ilong, ilat, ivert)
  integer :: dimid, len3d, cnt3d2d(2), strt3d2d(2), did_no_air, did_no_lgt, did_so2_3d, idx
  character :: temp3dVarName(speciesCharLen)
  character (len=speciesCharLen) :: tempVar3dNames
!
  character (len=MAX_LENGTH_SPECIES_NAME), parameter :: speciesName = 'species'
  character (len=MAX_LENGTH_SPECIES_NAME), parameter :: speciesName3D = 'species3d'
  character (len=MAX_LENGTH_SPECIES_NAME), parameter :: speciesDim3D = 'species3d_dim'
  character  :: tempVarName(speciesCharLen)
  character (len=speciesCharLen) :: speciesNameVar
!
! ----------------
! Begin execution.
! ----------------
!
!
  start5d(1) = i1  -  i1_gl + 1
  start5d(2) = ju1 - ju1_gl + 1
  start5d(3) = k1
  start5d_3d(1) = i1  -  i1_gl + 1
  start5d_3d(2) = ju1 - ju1_gl + 1
  start5d_3d(3) = k1
!
  count5d (:) = (/ ilong, ilat, ivert, 1, 1 /)
!.sds.. for 2d emissions file
  start4d(1) = i1  -  i1_gl + 1
  start4d(2) = ju1 - ju1_gl + 1
  count4d (:) = (/ ilong, ilat, 1, 1 /)
!
  call Ncop_Rd (ncid_em, emiss_infile_name)
!
  it = curRecord
  start5d(5) = it
!.sds.. start point for 3d array when reading emiss_2d fields - only have nox from aircraft
!...    when using 2d emissions (surface only)
  start5d_3d(5) = it
  start4d(4) = it
!
  count2d(:) = (/speciesCharLen, 1/)
  start2d(:) = (/1, 1/)
!
!... read NO_airand NO_lgt only once
  did_no_air = 0
  did_no_lgt = 0
  did_so2_3d = 0
!... zero out emission array before starting
  do ic = 1, num_emiss
    emissionArray(ic)%pArray3D(:,:,:) = 0.0d0
  enddo
!
  emiss_tmp_2d = 0.0d0
  emiss_tmp = 0.0d0
! ===============================
  SPCLOOP: do ic = 1, num_emiss
! ===============================
!
!   -------------------------------------------------------------
!   Check emiss_map to see if this is a chunk of data that we are
!   interested in.  If not, move on to next one.
!   -------------------------------------------------------------
!
    if (emiss_map(ic) > 0) then
!
!     ---------------------------------------------------------
!     Complete start5d array and read netCDF emissions data into
!     emiss_tmp.
!     ---------------------------------------------------------
!
      start5d(4) = ic
      start4d(3) = ic
!
!     -------------------------------------------------------------
!     Increment index into emiss array; recall that the emiss array
!     is only allocated for num_emiss species to readuce space.
!     -------------------------------------------------------------
!
!... do 2D emissions - all at surface layer
      if(emiss_var_name .eq. 'emiss_2d') then
         if(ic <= num_emiss-num_emiss3d) then
            call Ncrd_4d (emiss_tmp_2d, ncid_em, emiss_var_name, start4d, count4d)
      ! Get the name of the species read in
!
            if (firstReadEmiss) then
               start2d(2) = ic
!               call Ncrd_2d_Char (tempVarName, ncid_em, speciesName, start2d, count2d)
!               write(speciesNameVar, '(32A)') convertCharacters2String(tempVarName)
               call Ncrd_2d_Char (speciesNameVar, ncid_em, speciesName, start2d, count2d)
               emittedSpeciesNames(ic) = TRIM(speciesNameVar)
            endif
         else
            start5d_3d(4) = ic-(num_emiss-num_emiss3d)
            call Ncrd_5d (emiss_tmp, ncid_em, 'emiss_3d', start5d_3d, count5d)
       ! Get the name of the species read in
!
            if (firstReadEmiss) then
               call nf_inq_dimid (ncid_em, speciesDim3D, dimid)
               call nf_inq_dimlen (ncid_em, dimid, len3d)
               strt3d2d(:) = (/1, start5d_3d(4)/)
               cnt3d2d(:) = (/speciesCharLen, 1/)
               call Ncrd_2d_Char (tempVar3dNames, ncid_em, speciesName3D, strt3d2d, cnt3d2d)
               emittedSpeciesNames(ic) = TRIM(tempVar3dNames)
            endif
         endif
!.. do the original way - all are 3D emissions
      else
         call Ncrd_5d (emiss_tmp, ncid_em, emiss_var_name, start5d, count5d)
!... Get the name of the species read in
         if (firstReadEmiss) then
            start2d(2) = ic
!            call Ncrd_2d_Char (tempVarName, ncid_em, speciesName, start2d, count2d)
!            write(speciesNameVar, '(32A)') convertCharacters2String(tempVarName)
            call Ncrd_2d_Char (speciesNameVar, ncid_em, speciesName, start2d, count2d)
            emittedSpeciesNames(ic) = TRIM(speciesNameVar)
         endif
      endif
!
!     zero out the lightning contribution if parameterized
!     lightning is on
!
!     -------------------------------------------------------
!     Force no NO from lightning emissions unless it is
!     specifically requested (lightning_opt=0)
!     -------------------------------------------------------
      if ((lightning_opt .ne. 0) .and. ((emittedSpeciesNames(ic)(1:6)) .eq. "NO_lgt")) then
         ! print '(''Zeroing NO_lgt emission array read in in ReadEmiss'')'
         emiss_tmp = 0.
      endif
!
!     -------------------------------------------------------
!     Load emiss_tmp (i.e., emiss data that was just read in)
!     into proper location in emiss array.
!     -------------------------------------------------------
       if(emiss_var_name .eq. 'emiss_2d') then
         if(ic <= (num_emiss - num_emiss3d) ) then
!... do 2D emissions - all at surface layer
            do ij = ju1, j2
               jnb = ij - ju1 + 1
               do il = i1, i2
                  inb = il - i1 + 1
                  emissionArray(ic)%pArray3D(il,ij,k1) = emiss_tmp_2d(inb,jnb)  &
                       + emissionArray(ic)%pArray3D(il,ij,k1)
               end do
            end do
!... do 3D emissions
         else
            do ik = k1, k2
               do ij = ju1, j2
                  jnb = ij - ju1 + 1
                  do il = i1, i2
                     inb = il - i1 + 1
                     emissionArray(ic)%pArray3D(il,ij,ik) = emiss_tmp(inb,jnb,ik)  &
                       + emissionArray(ic)%pArray3D(il,ij,ik)
                  enddo
               enddo
            enddo
         endif
!.. do the original way - all are 3D emissions
      else
        do ik = k1, k2
            do ij = ju1, j2
               jnb = ij - ju1 + 1
               do il = i1, i2
                  inb = il - i1 + 1
                  emissionArray(ic)%pArray3D(il,ij,ik) = emiss_tmp(inb,jnb,ik)  &
                       + emissionArray(ic)%pArray3D(il,ij,ik)
               end do
            end do
         end do
      endif
!
!     ---------------------------------------------------------
!     If all of the emissions data has been read in, no need to
!     check the rest of the species, so exit species loop.
!     ---------------------------------------------------------
!
!
    end if
!       ==============
  end do SPCLOOP
!       ==============
!
  call Nccl (ncid_em)
!
  if (firstReadEmiss) firstReadEmiss = .false.
!
  return
!
  end subroutine ReadEmiss
!
!-------------------------------------------------------------------------
! This routine reads in the aero (carbon & sslt) emission data.
!-------------------------------------------------------------------------
!
      subroutine ReadEmissAero (emissAero_t, emiss_aero_infile_name, naero, emiss_timpyr, &
                         pr_diag, loc_proc, i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat
      integer, intent(in) :: loc_proc, naero, emiss_timpyr
      logical, intent(in) :: pr_diag
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: emiss_aero_infile_name
      real*8, intent(inout) :: emissAero_t(i1:i2, ju1:j2, 1:naero, emiss_timpyr)
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      character*30, parameter :: EMISS_AERO_NAME = 'emiss_carbon_sslt'
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ic, it, inb, jnb, inum, ncid_em_aero
!
      integer :: count4d (4), start4d(4)
!
      real*8  :: emissAero_tmp(ilong, ilat, 1, 1)
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'ReadEmissAero called by ', loc_proc
      end if
!
      start4d(1) = i1  -  i1_gl + 1
      start4d(2) = ju1 - ju1_gl + 1
!
      count4d (:) = (/ ilong, ilat, 1, 1 /)
!
      call Ncop_Rd (ncid_em_aero, emiss_aero_infile_name)
!
      do it = 1, emiss_timpyr
        inum = 0
        start4d(4) = it
!
        do ic = 1, naero
!
          start4d(3) = ic
!
          call Ncrd_4d (emissAero_tmp, ncid_em_aero, EMISS_AERO_NAME,  &
     &       start4d, count4d)
!
          inum = inum+1
!
          do ij = ju1, j2
            jnb = ij - ju1 + 1
            do il = i1, i2
              inb = il - i1 + 1
!
              emissAero_t(il,ij,inum,it) = emissAero_tmp(inb,jnb,1,1)
!
            end do
          end do
!
          if (inum == naero) exit
!
        end do
      end do
!
      call Nccl (ncid_em_aero)
!
      return
!
      end subroutine ReadEmissAero
!
!-------------------------------------------------------------------------
! This routine reads in the dust emission data.
!-------------------------------------------------------------------------
!
      subroutine ReadEmissDust (emissDust_t, emiss_dust_infile_name, &
                                ndust, nt_dust, nst_dust, &
                                do_aerocom, pr_diag, loc_proc, &
                                i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat
      integer, intent(in) :: loc_proc, ndust, nt_dust, nst_dust
      logical, intent(in) :: pr_diag, do_aerocom
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: emiss_dust_infile_name
      real*8, intent(inout) :: emissDust_t(i1:i2, ju1:j2, 1:ndust, nst_dust:nst_dust+nt_dust-1)
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      character(len=30) :: EMISS_DUST_NAME
      integer :: il, ij, ic, inb, jnb, ist, it, ncid_em_dust
      integer :: count4d (4),  start4d(4)
      real*8  :: emissDust_tmp(ilong, ilat, ndust, 1)
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'ReadEmissDust called by ', loc_proc
      end if
!
      if (do_aerocom) then
           EMISS_DUST_NAME = 'emiss_dms_dust_ss'
      else
           EMISS_DUST_NAME = 'emiss_dust'
      end if
!
      start4d(1) = i1  -  i1_gl + 1
      start4d(2) = ju1 - ju1_gl + 1
      start4d(3) = 1
!
      count4d (:) = (/ ilong, ilat, ndust, 1 /)
!
      call Ncop_Rd (ncid_em_dust, emiss_dust_infile_name)
!
      do it = 1, nt_dust
!
        ist = nst_dust + it - 1
!
        start4d(4) = it
!
        call Ncrd_4d(emissDust_tmp, ncid_em_dust, EMISS_DUST_NAME, start4d, count4d)
!
        do ic = 1, ndust
          do ij = ju1, j2
            jnb = ij - ju1 + 1
            do il = i1, i2
              inb = il - i1 + 1
              emissDust_t(il,ij,ic,ist) = emissDust_tmp(inb,jnb,ic,1)
            end do
          end do
        end do
!
      end do
!
      call Nccl (ncid_em_dust)
!
      return
!
      end subroutine ReadEmissDust
!-----------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readGalacticCosmisRayPar
!
! CODE DEVELOPER
!   Stephen Steenrod ; GSFC
!   stephen.d.steenrod@nasa.gov
!
! !INTERFACE:
!
  subroutine readGalacticCosmisRayPar(gcr_infile_name, gcr_sunspot, gcr_slope, gcr_aintcp, &
       nymd, am, bm, pt, latdeg, ju1, j2, k1, k2, ju1_gl, pr_diag, procID)
!
!
! !USES
  use GmiInterpolation_mod, only : Interp_Bilinear
!
!
  implicit none
!
! !INPUT PARAMETERS:
  character (len=MAX_LENGTH_FILE_NAME), intent(in) :: gcr_infile_name  ! Galactic Cosmic Ray input filename
  real*8, intent(out) :: gcr_sunspot                 ! Galactic Cosmic Ray parameter
  real*8, intent(out) :: gcr_slope (ju1:j2, k1:k2)   ! Galactic Cosmic Ray parameter
  real*8, intent(out) :: gcr_aintcp(ju1:j2, k1:k2)   ! Galactic Cosmic Ray parameter
  integer, intent(in) :: ju1, j2, k1, k2, ju1_gl, nymd, procID
  real*8 , intent(in) :: pt, am(k1:), bm(k1:), latdeg(ju1_gl:)
  logical, intent(in) :: pr_diag
!
! !DESCRIPTION:
! Reads a 2D parameters of Galactic Cosmic Ray NOx emissions.
!   Data from Charley Jackman via the Goddard CTM.
!
! !DEFINED PARAMETERS:
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_DNAM      = 'latitude_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRES_DNAM     = 'pressure_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: TIME_DNAM     = 'time_dim'
!
  character (len=MAX_LENGTH_VAR_NAME), parameter :: GCS_spt_VNAM  = 'GCR_Sunspot_coeff'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: GCS_slp_VNAM  = 'GCR_Slope_coeff'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: GCS_aint_VNAM = 'GCR_Aintcp_coeff'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: GCS_time_VNAM = 'time_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: GCS_pres_VNAM = 'pressure_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: GCS_lat_VNAM  = 'latitude_dim'
!
! !LOCAL VARIABLES:
!
  character (len=75) :: err_msg
  integer :: ii, ij, ik, icnt
  integer :: ii_num, jj_num, kk_num
  integer :: idxyr, inyr
  integer :: ncid_gcr
  integer :: cnt1d (1)
  integer :: strt1d(1)
  integer :: cnt2d (2)
  integer :: strt2d(2)
  integer, allocatable :: itemp (:)
  real*8, allocatable :: temp (:)
  real*8, allocatable :: temp2d (:,:)
  real*8  :: rtemp
!... variables read in from file
  real*8              :: logplev
  real*8, allocatable :: in_gcr_lat (:)
  real*8, allocatable :: in_gcr_logplev (:)
  integer,allocatable :: in_gcr_date (:)
  real*8, allocatable :: in_gcr_sunspot (:)
  real*8, allocatable :: in_gcr_slope (:,:)
  real*8, allocatable :: in_gcr_aintcp (:,:)
!
!EOP
!-------------------------------------------------------------------------
!BOC
!     ---------------------------------------------------------
!     Do some error checking to make sure the dimensions in the
!     table are consistant with the dimensions in the include
!     file that describes the table and with lookup table file.
!     ---------------------------------------------------------
!
  if (pr_diag)  Write (6,*) 'readGalacticCosmisRayPar called by ', procID
!
  call Ncop_Rd (ncid_gcr, gcr_infile_name)
!
!... set dimensions of netCDF file for GCR parameterization
  call Ncget_Dimlen (ncid_gcr, LAT_DNAM , ii_num)
  call Ncget_Dimlen (ncid_gcr, PRES_DNAM, jj_num)
  call Ncget_Dimlen (ncid_gcr, TIME_DNAM, kk_num)
!
!... read in data
  strt1d(:) = (/ 1 /)
  cnt1d(1) = kk_num
!
!... get years of input data
  allocate(itemp(kk_num))
  call Ncrd_1d_Int (itemp, ncid_gcr, GCS_time_VNAM, strt1d, cnt1d)
!... calculate index for date we want
  inyr = int(nymd/10000)
  idxyr = inyr-itemp(1)+1
!... make sure indice is in range, if not adjust wrt 11 year solar cycle
  do while (idxyr <= 0)
    idxyr = idxyr+11
  enddo
  do while (idxyr >= kk_num)
    idxyr = idxyr-11
  enddo
!
!... get sunspot coeff of input data
  allocate(temp(kk_num))
  call Ncrd_1d (temp, ncid_gcr, GCS_spt_VNAM, strt1d, cnt1d)
  gcr_sunspot = temp(idxyr)
!
!... get pressures of input data
  allocate(in_gcr_logplev(jj_num))
  cnt1d(1) = jj_num
  call Ncrd_1d (in_gcr_logplev, ncid_gcr, GCS_pres_VNAM, strt1d, cnt1d)
!... take log for linear-log(p) interpolation
  do ik = 1, jj_num
    in_gcr_logplev(ik) = log(in_gcr_logplev(ik))
  enddo
!
!... get latitudes of input data
  allocate(in_gcr_lat(ii_num))
  cnt1d(1) = ii_num
  call Ncrd_1d (in_gcr_lat, ncid_gcr, GCS_lat_VNAM, strt1d, cnt1d)
!
!... read in the GCR parameters
  strt2d(:) = (/ 1, 1 /)
  cnt2d(:) = (/ ii_num, jj_num /)
!
  allocate(in_gcr_slope(ii_num,jj_num))
  call Ncrd_2d (in_gcr_slope, ncid_gcr, GCS_slp_VNAM, strt2d, cnt2d)
!
  allocate(in_gcr_aintcp(ii_num,jj_num))
  call Ncrd_2d (in_gcr_aintcp, ncid_gcr, GCS_aint_VNAM, strt2d, cnt2d)
!
  call Nccl (ncid_gcr)
!
!... interpolate input data to current model grid
  do ik = k1, k2
!... approximate the grid pt pressure by using sfc pressure of 1000 hPa
    logplev = log(am(ik)*pt + bm(ik)*1000.)
!
!... Need to interpolate gcr_slope and gcr_aintcp to current grid
    do ij = ju1, j2
!... interpolate gcr_slope to current grid
      call Interp_Bilinear(latdeg(ij), logplev, rtemp,  &
                 in_gcr_lat, in_gcr_logplev, ii_num, jj_num, in_gcr_slope)
      gcr_slope(ij,ik) = rtemp
!
!... interpolate gcr_aintcp to current grid
      call Interp_Bilinear(latdeg(ij), logplev, rtemp,  &
                 in_gcr_lat, in_gcr_logplev, ii_num, jj_num, in_gcr_aintcp)
      gcr_aintcp(ij,ik) = rtemp
    enddo
!
  enddo
!
!
  return
!
  end subroutine readGalacticCosmisRayPar
!
!-----------------------------------------------------------------------------------
!
end module ReadEmissionFiles_mod
!
