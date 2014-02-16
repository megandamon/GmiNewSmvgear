module GmiEmissionLightning_mod

implicit none

private

public :: emiss_lightning

# include "GmiParameters.h"

contains

!=============================================================================
!
!  Dale Allen's Lightning option, added Sep 29, 2005
!  Latest version: June 16, 2009 
!
! FILE
!   emiss_lightning.F
!
! ROUTINES
!   emiss_lightning
!   noxlightning_gmi
!   partitionnox  (Latest version 6/16/2009)
!                  (Variable cloud top added 2/3/2011) 
!   flashfit        (Latest version: 2/3/2011; cldmas 0.025 constraint added)
!             
!   NOTE: Code assumes metdata_name_org for GEOS-4 FDAS simulation is G4FDAS.  Change if necessary.
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   emiss_lightning
!
! DESCRIPTION
!   This routine adds lightning emissions to const.
!
!   input:
!     1) cldmas0 : normalized cldmas at desired level (2d array),
!                    for fvgcm: zmmu (gmi var)
!                    for others: cmf (gmi var)
!   output:
!     1) flashrate : flash rate (flashes per grid box per s)
!     2) lightning_no : 3d array of pnox production in kg./sec
!     3) Adds lightning emissions to existing mixing ratio of const
!        gmi array const is modified.
!        const(no) = const(no) + pnox3d * tdt
!        important: in the following, it is assumed that NO is
!                   species 51
!
! ARGUMENTS
!     pnox = pnox * desired_g_N_prod_rate /3.42
!   desired_g_N_prod_rate: namelist variable, results scaled to produce
!                          desired_g_N_prod_rate in Tg form assumed 3.42 Tg
!   cldmas0      : normalized cldmas at desired level
!   press3c      : atmospheric pressure at the center of each grid box (mb)
!   press3e      : atmospheric pressure at the edge of each grid box (mb)
!   cmi_flags    : array of flags that indicate continental, marine, or ice
!   mass         : mass of air in each grid box (kg)
!   dtrn         : detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
!   const        : species concentration, known at zone centers (mixing ratio)
!   cldmas       : convective mass flux in     updraft (kg/m^2*s)
!   flashrate    : flash rate (flashes per 4x5 box per s)
!   lightning_no : 3d array of pnox production in kg./sec
!
!  !REVISION HISTORY:
!   - July 28, 2008 - Dale Allen (first version)  

!-----------------------------------------------------------------------------

      subroutine emiss_lightning(metdata_name_org, metdata_name_model, nymd, &
     &  i1,i2 ,ju1,j2, k1, k2, num_species, ilo, ihi, julo, jhi, tdt, &
     &  procID, pr_diag, desired_g_N_prod_rate, cldmas0, press3c, press3e, &
     &  cmi_flags,  cmi_flags1, i2_gl, j2_gl,  &
     &  mass, dtrn, pnox3d, flashrate, lightning_no, &
     &  threshold,ratio_global,ratio_local,midLatAdj)

!      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod,        only : GmiSplitDateTime

      implicit none
!     ----------------------
!     Argument declarations.
!     ----------------------

      character (len=MAX_LENGTH_MET_NAME), intent(in)  :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME), intent(in)  :: metdata_name_model
      logical, intent(in) :: pr_diag
      integer,intent(in)   :: nymd
      integer,intent(in)   :: i2_gl, j2_gl, i1, i2, ju1, j2, k1, k2, num_species
!      integer,intent(in)   :: ino_num
      integer,intent(in)   :: procID
      integer,intent(in)   :: ilo, ihi, julo, jhi
      real*8,intent(in)    :: tdt
      real*8,  intent(in)  :: desired_g_N_prod_rate
      real*8,intent(in)    :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8,intent(in)    :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      integer, intent(in)  :: cmi_flags(i1:i2, ju1:j2)
      integer, intent(in)  :: cmi_flags1(i1:i2, ju1:j2)
      real*8,intent(in)    :: mass (i1:i2, ju1:j2, k1:k2)
      real*8 ,intent(in)   :: dtrn(i1:i2, ju1:j2, k1:k2)
      !real*8 ,intent(in)   :: pctm2(ilo:ihi, julo:jhi)
      
      real*8, intent(in) :: threshold                          
      real*8, intent(in) :: ratio_global             
      real*8, intent(in) :: ratio_local(i1:i2,ju1:j2) 
      real*8, intent(in) :: midLatAdj(i1:i2,ju1:j2) 

      real*8,intent(inout)    :: cldmas0(i1:i2, ju1:j2)
!      type (t_GmiArrayBundle), intent (inout) :: concentration(num_species)

      real*8,intent(out)   :: flashrate (i1:i2, ju1:j2)
      real*8,intent(out)   :: lightning_no (i1:i2, ju1:j2, k1:k2)
!     ----------------------
!     Local Variable declarations.
!     ----------------------

      real*8  :: pnox2d (i1:i2, ju1:j2)         ! Lightning NO production (g N/s)
      real*8  :: pnox3d ( i1:i2, ju1:j2, k1:k2) ! Lightning NO production (ppv/s)

      real*8  :: dlatcen2d(i1:i2, ju1:j2)
      real*8  :: dloncen2d(i1:i2, ju1:j2)

      real*8  :: pi, dtr, factor, d_ptop, ppmtoppv
      integer :: i, j, ispec, year,imonth,day, k
      real*8  :: dilon, djlat
      integer :: local_rank
      integer :: ierr
      
      real*8  :: cldmas1(i1:i2, ju1:j2)
      real*8 :: thresholdTmp
      
!EOP
!-------------------------------------------------------------------------------
!BOC

      if (pr_diag) then
         Write(6,*) " emiss_lightning called by ", procID
      end if

      ! adjust threshold to the time units of cldmas s-1
      thresholdTmp = threshold / 60.

      pi = 4.0 * atan(1.); dtr = pi / 180.;
      !ppmtoppv --> converts from mass of air to mass of N.
      ppmtoppv = 28.97 / 14.

      if (metdata_name_org(1:3) == 'DAO') then
         factor = 60.;
      else if(metdata_name_org(1:4) == 'NCAR') then
         factor = 60.;
      else if(metdata_name_org(1:4) == 'GISS') then
         factor = 60. * 9.8;
      else if(metdata_name_org(1:4) == 'GMAO') then
         factor = 60. 
      else if(metdata_name_org(1:6) == 'G4FDAS') then
         factor=60. 
      else
          write (6,*) ' Lightning parameterization does',  &
     &                '  not exist for this model'
      endif

      !get month
      call GmiSplitDateTime (nymd, year, imonth, day)

! Define values for longitude and latitude. note: applies to 4 x 5 grid
      if (i2_gl .eq. 144) then
        dilon= 2.5
      else
        dilon = 5.0
      endif

      if (j2_gl .eq. 91) then
        djlat= 2.0
      else
        djlat = 4.0
      endif

      !longitude

      do i = i1,i2
         dloncen2d(i,:) = (i-1) * dilon
      end do
      !latitude
      do j=ju1, j2
         dlatcen2d(:,j) = -90. + (j-1) * djlat
      enddo

!      if(ju1.eq.1)  dlatcen2d(:,1)  = -89.0
!      if(j2.eq.j2_gl)  dlatcen2d(:,46) =  89.0

      
      call noxlightning_gmi(i2_gl, j2_gl, i1,i2, ju1,j2, ilo,ihi, julo,jhi,  &
     & k1,k2, imonth, desired_g_N_prod_rate,cldmas0,  &
     & dlatcen2d, dloncen2d,  &
     & press3c, press3e, cmi_flags, cmi_flags1, mass, dtrn,  &
     & flashrate, pnox2d, pnox3d,thresholdTmp,ratio_global,ratio_local,midLatAdj)

! in gem/actm/gmimod/chem/troposphere/include_setkin/setkin_lchem.h
! data lchemvar(51) /"NO"/ ; data lchemvar(52) /"NO2"/   new

!      ispec = ino_num
!      concentration(ispec)%pArray3D(:,:,:) = &
!     &          concentration(ispec)%pArray3D(:,:,:) + pnox3d(:,:,:) * tdt

      !diagnostic for lighting emissions
!make sure that lightning_nc = 0.  at the start of every writing interval


!     Nitrogen production in  kg. /sec.
      lightning_no(:,:,:) =  (pnox3d(:,:,:)*mass(:,:,:)) / ppmtoppv

      return
      end subroutine emiss_lightning
!=============================================================================

!-----------------------------------------------------------------------------
      subroutine noxlightning_gmi(i2_gl, j2_gl, i1,i2, ju1,j2, ilo,ihi, julo,jhi,  &
     & k1,k2, imonth, desired_g_N_prod_rate, cldmas0,  &
     &        dlatcen, dloncen,  &
     &        pcen, lprslay, cmi_flags, cmi_flags1, mass, dtrn,  &
     &        flashrate, pnox2d, pnox3d,threshold,ratio_global,ratio_local, midLatAdj)

      implicit none


      integer, intent(in) :: i2_gl, j2_gl, i1,i2, ju1,j2, ilo,ihi, julo,jhi,  &
     &                       k1,k2, imonth
      real*8,  intent(in) :: desired_g_N_prod_rate
      real*8,  intent(inout) :: cldmas0(i1:i2, ju1:j2)
      real*8,  intent(in) :: dlatcen(i1:i2, ju1:j2)
      real*8,  intent(in) :: dloncen(i1:i2, ju1:j2)
      real*8,  intent(in) :: pcen(ilo:ihi, julo:jhi, k1:k2)
      real*8,  intent(in) :: lprslay(ilo:ihi, julo:jhi, k1-1:k2)
      integer, intent(in) :: cmi_flags(i1:i2, ju1:j2)
      integer, intent(in) :: cmi_flags1(i1:i2, ju1:j2)
      real*8,  intent(in) :: mass(i1:i2, ju1:j2, k1:k2)
      real*8,  intent(in) :: dtrn(i1:i2, ju1:j2, k1:k2)
      !real*8,  intent(in) :: pctm2   (ilo:ihi, julo:jhi)


!     Mass flux threshold (kg m-2 min-1) below which flash rate is assumed to be zero (Read in from file) 
      real*8, intent(in) :: threshold                          

!     Adjustment factor local flash rates must be multiplied so that globally averaged 
!     flash rate matches v2.2 OTD/LIS climatological globally averaged flash rate (Read in from file) 
!     Note: ratio_global varies monthly. 
      real*8, intent(in) :: ratio_global             

!     Adjustment factor local flash rates must be multiplied so that monthly average local flash rates (after
!     "ratio_global" adjustment) match monthly average local v2.2 climatological OTD/LIS flash rates (Read in from file).  
!     Note: ratio_local varies monthly. 
      real*8, intent(in) :: ratio_local(i1:i2,ju1:j2) 

      real*8, intent(in) :: midLatAdj(i1:i2,ju1:j2) 

!     OUTPUT FIELDS
!     flashrate: 2-d: lightning frequency (flashes per 4x5 box per s)
!     pnox2d:    2-d: lightning NOx production rate in column (g N/s)
!     pnox3d:    3-d: lightning NOx production rate (ppv/s)

      real*8, intent(inout) :: flashrate(i1:i2, ju1:j2)
      real*8, intent(out) :: pnox2d(i1:i2, ju1:j2)
      real*8, intent(out) :: pnox3d(i1:i2, ju1:j2, k1:k2)


!     INPUT FIELDS:
!     i1,i2:         0-d: number of east-west grid points
!     ju1,j2:        0-d: number of north-south grid points
!     ilo,ihi, julo,jhi:
!     k1,k2:         0-d: number of vertical layers
!     imonth
!     dlatcen        2-d: Latitude at center of grid boxes (degrees)
!     dloncen        2-d: Longitude at center of grid boxes (degrees)

!     pcen:          3-d: Approx. pressure at centers of layers: (hPa)
!     lprslay:       3-d: Approx. pressure at top edges of layers: (hPa)

!------------------------------------------------------------------------
!     mass:      3-d: Mass of grid volumes (kg) [100*dx*dy*dp]/g]
!     cmf or cldmas0:        3-d: Convective mass flux (kg m-2 min-1)
!     zmmu or cldmas0:       3-d: Convective mass flux (kg m-2 s-1) updraft
!     dtrn:        3-d: Detrainment (kg m-2 s-1) (units not crucial)
!     lwi_flags:      2-d: Marine(1), Continental (2) ice(3) flag
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate flash rate and resulting lightning NOx prod. rate in column.      
      call flashfit(i2_gl,j2_gl,i1,i2, ju1,j2,&
     & desired_g_N_prod_rate, threshold,ratio_global,ratio_local,&
     & midLatAdj, cldmas0, flashrate, pnox2d)
     
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Partition lightning NOx in the vertical within each column.

      call partitionnox(i1,i2, ju1,j2, ilo,ihi, julo,jhi, k1,k2,  &
     & dlatcen, lprslay, mass, dtrn, cmi_flags,  &
     & pnox2d, pnox3d,imonth,flashrate, i2_gl, j2_gl)
!-----------------------------------------------------------------------
      return
      end subroutine noxlightning_gmi

 
      subroutine partitionnox(i1,i2, ju1,j2, ilo,ihi, julo,jhi, k1,k2, &
     & dlatcen, lprslay, mass, dtrn, cmi_flags,  &
     & pnox2d, pnox3d,imonth,flashrate, i2_gl, j2_gl)

!     Modeling studies of Pickering (2007) are used as a basis
!      to partition lightning NO in the vertical.

!     Update 4/4/08 (cont. tropical partitioning changed) 
!     Update 6/16/09 (Midlat model cloud profile shifted upward) 

      implicit none

!     INPUTS
!     dlatcen:  2-d: Latitude at grid box centers (degrees N)
!     lprslay:  3-d: Pressure at grid box edges (hPa)
!     mass:     3-d: Background (atmospheric) mass in each grid volume (kg)
!     dtrn:     3-d: Detrainment rate (kg m-2 s-1) (units not important)
!     cmi_flags 2-d: continental(2), marine(1), or ice flag.
!     pnox2d:   2-d: lightning NOx production rate in column (g N/s)

      integer, intent(in)  :: i1,i2, ju1,j2, ilo,ihi,julo,jhi,k1,k2
      integer, intent(in)  :: i2_gl, j2_gl
      integer, intent(in)  :: imonth
      real*8,  intent(in)  :: dlatcen(i1:i2, ju1:j2)
      real*8,  intent(in)  :: lprslay(ilo:ihi, julo:jhi, k1-1:k2)
      real*8,  intent(in)  :: mass(i1:i2, ju1:j2, k1:k2)
      real*8,  intent(in)  :: dtrn(i1:i2, ju1:j2, k1:k2)
      integer, intent(in)  :: cmi_flags(i1:i2, ju1:j2)
      !real*8,  intent(in)  :: pctm2   (ilo:ihi, julo:jhi)
      real*8,  intent(in)  :: pnox2d(i1:i2, ju1:j2)

!     OUTPUTS
!     pnox3d:  lightning NOx production rate in each grid volume (ppv/s)

      real*8,  intent(out) :: pnox3d(i1:i2, ju1:j2, k1:k2)
      real*8, intent(inout) :: flashrate(i1:i2, ju1:j2)

!     LOCAL VARIABLES:
      integer :: il,ij,ik,ii,ikk,iktop

      real*8 :: ppmtoppv                  ! conversion from parts per mass to parts per vol.

      real*8 :: ztop(i1:i2,ju1:j2,k1:k2 ) ! Heights at top edges of model layers.
      real*8 :: htedge2(k1:k2)            ! Top edge heights after adjusting to match
                                          ! "model-cloud".
      real*8 :: fd0(k1-1:k2)              ! Used to calc cont of each "model-cloud"
                                          ! layer to each CTM cloud layer
      real*8 :: yout(k1:k2)               ! % of lightning NO mass deposited into
                                          ! each CTM cloud layer
      integer :: inox1(i1:i2, ju1:j2)     ! "Model-cloud" type indices
      integer :: ntop(i1:i2, ju1:j2)      ! Top cloud layer indices
      integer :: ikmm                     !Number of vertical layers in model

      integer,parameter :: pmax(4) = (/ 17, 17, 16, 15 /)
      real*8 :: r0(17,4)


      ikmm = 40

!     r0(*,i): % of lightning NOy mass dep. into each layer for cld type i
!     "Model cloud" has 17 1-km layers. Percent of total NOx mass
!     deposited into each layer is given for 1) tropical continental,  
!     2) tropical marine, 3) subtropical, and 4) midlatitude "model clouds". 

      data r0/  0.23, 0.47, 0.56, 1.40, 2.70, 4.00, 5.03, 6.24, &
     &          8.60,10.28,11.62,12.34,12.70,12.34, 7.63, 3.02, 0.84, &        
     &          0.60, 1.50, 2.90, 4.30, 5.40, 6.70, 7.70, 8.50, &
     &          9.60,10.20,10.50,10.20, 8.20, 6.50, 4.50, 2.20, 0.50, & 
     &          1.00, 2.10, 3.90, 5.80, 7.70, 9.30,10.50,11.00, &
     &         11.00,10.40, 9.20, 7.50, 5.50, 3.40, 1.50, 0.20, 0.00, &
     &          1.00, 2.37, 4.95, 7.32, 9.21,10.49,11.29,11.39, &
     &         10.89, 9.80, 8.22, 6.24, 4.16, 2.18, 0.49, 0.00, 0.00/ 

!     Assuming a scale height of 8 km, estimate heights at layer tops.
      do ij = ju1,j2
         do il = i1,i2
            do ik= k1,k2
            ztop(il,ij,ik) = -8.* dlog(lprslay(il,ij,ik)/lprslay(il,ij,0))
            !ztop(il,ij,ik) = -8.* dlog(lprslay(il,ij,ik)/pctm2(il,ij))
            end do
         end do
      end do


!     Determine which "model cloud" to use at each grid point. Decision
!     is based on whether pts are tropical continental, tropical marine, 
!     subtropical, or midlatitude. 
      select case (imonth)
      case (1,2,3,12)
!     Southern Hemisphere Summer      
         where     ((abs(dlatcen) <= 15.).and.(cmi_flags == 2))
            inox1 = 1      ! Tropical continental
         elsewhere     ((abs(dlatcen) <= 15.).and.(cmi_flags /= 2))
            inox1 = 2      ! Tropical marine
         elsewhere ((dlatcen > 15.).and.(dlatcen <= 30.))
            inox1 = 3      ! Subtropical 
         elsewhere ((dlatcen >= -40.).and.(dlatcen < -15.))
            inox1 = 3      ! Subtropical 
         elsewhere 
            inox1 = 4      ! Midlatitude 
         end where
!     In-between months     
      case (4,5,10,11)
         where     ((abs(dlatcen) <= 15.).and.(cmi_flags == 2))
            inox1 = 1      ! Tropical continental
         elsewhere     ((abs(dlatcen) <= 15.).and.(cmi_flags /= 2))
            inox1 = 2      ! Tropical marine
         elsewhere ((dlatcen > 15.).and.(dlatcen <= 30.))
            inox1 = 3      ! Subtropical 
         elsewhere ((dlatcen >= -30.).and.(dlatcen < -15.))
            inox1 = 3      ! Subtropical 
         elsewhere 
            inox1 = 4      ! Midlatitude 
         end where
!     Northern Hemisphere Summer      
      case (6,7,8,9)
         where     ((abs(dlatcen) <= 15.).and.(cmi_flags == 2))
            inox1 = 1      ! Tropical continental
         elsewhere     ((abs(dlatcen) <= 15.).and.(cmi_flags /= 2))
            inox1 = 2      ! Tropical marine
         elsewhere ((dlatcen > 15.).and.(dlatcen <= 40.))
            inox1 = 3      ! Subtropical 
         elsewhere ((dlatcen >= -30.).and.(dlatcen < -15.))
            inox1 = 3      ! Subtropical 
         elsewhere 
            inox1 = 4      ! Midlatitude 
         end where
      case default
         print*, "CASE default in partitionnox"
      end select 	 

!     Define cloud top to be highest layer with dtrn > 0.
!     Caution: Code assumes cloud does not extend above model layer 30.
!     Warning.  This statement should be examined as vertical res is increased.
      iktop = min(40,k2)

!     Determine model layer of cloud top.
      do ij=ju1,j2
      do il=i1,i2
         do ik=iktop,2,-1
            if (dtrn(il,ij,ik).gt.0) then
               ntop(il,ij) = ik
               goto 213
            end if
            ntop(il,ij) = 2
         end do
 213     continue
      end do
      end do


!     note: ppmtoppv converts from ppm to ppv.  N(MW) = 14.
      pnox3d = 0.

      ppmtoppv = 28.97 / 14.
      do ij=ju1,j2
      do il=i1,i2
         if (pnox2d(il,ij).gt.0.0) then
            iktop = ntop(il,ij)                 !Index of top cloud layer
            ii = inox1(il,ij)                   !"Model-cloud index"

!           CTM "cloud" is expanded or contracted so that its depth = 17 km, which is the
!           assumed depth of the "model-cloud".
            htedge2 =  ztop(il,ij,:)*real(pmax(ii))/ztop(il,ij,iktop)

            yout = 0.
            fd0 = 0.
	    do ik=1,pmax(ii) 
               do ikk=1,ikmm
                  fd0(ikk) = max(min(htedge2(ikk)-(ik-1),1.),0.)
               end do
               do ikk=1,ikmm-1
                  yout(ikk) = yout(ikk)+r0(ik,ii)*(fd0(ikk)-fd0(ikk-1))
               end do
            enddo
            yout = yout * 1. / sum(yout)

            do ikk=1,iktop
               pnox3d(il,ij,ikk)=yout(ikk)*pnox2d(il,ij)
            end do

         end if
      end do
      end do


!     Convert from g N s-1 to ppv s-1.
      pnox3d = 0.001 * pnox3d * ppmtoppv / mass

 
 

      return
      end subroutine partitionnox 


!=============================================================================

!-----------------------------------------------------------------------------
      subroutine flashfit(i2_gl,j2_gl,i1,i2, ju1,j2,&
     & desired_g_N_prod_rate, threshold,ratio_global,ratio_local,&
     & midLatAdj, cldmas, flashrate, pnox)

      implicit none

      integer, intent(in) :: i2_gl, j2_gl,i1,i2,ju1,j2
      
!     Desired mean lightning NO production rate (Tg yr-1) (specified in namelist)       	
      real*8, intent(in) :: desired_g_N_prod_rate
                  
!     Mass flux threshold (kg m-2 min-1) below which flash rate is assumed to be zero (Read in from file) 
      real*8, intent(in) :: threshold                          

!     Adjustment factor local flash rates must be multiplied so that globally averaged 
!     flash rate matches v2.2 OTD/LIS climatological globally averaged flash rate (Read in from file) 
!     Note: ratio_global varies monthly. 
      real*8, intent(in) :: ratio_global             

!     Adjustment factor local flash rates must be multiplied so that monthly average local flash rates (after
!     "ratio_global" adjustment) match monthly average local v2.2 climatological OTD/LIS flash rates (Read in from file).  
!     Note: ratio_local varies monthly. 
      real*8, intent(in) :: ratio_local(i1:i2,ju1:j2) 

      real*8, intent(in) :: midLatAdj(i1:i2,ju1:j2) 

!     Convective mass flux at desired layer (kg m-2 min-1)
!     Warning: Make sure units are correct and that you've accessed the proper vertical layer. 
      real*8, intent(in) :: cldmas(i1:i2, ju1:j2)    

!     Total (CG+IC) flash rate (flashes per grid box per s)
      real*8, intent(out) :: flashrate(i1:i2, ju1:j2)	

!     Lightning NO production rate (g N per grid box per s)
      real*8, intent(out) :: pnox(i1:i2, ju1:j2)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 :: prodfac, convfac
      real*8 :: cldmas_local(i1:i2, ju1:j2)
      integer :: i,j

      prodfac = 1.0e26        ! (joules per flash) * molecules of NO per joule
      convfac = 2.33e-23      ! 14g of N per mole of NO/ 6.02e23 molecules of NO/mole of NO.

!     If deep convection exists

      cldmas_local(:,:) = cldmas(:,:)
      where(cldmas_local > 0.025) cldmas_local = 0.025 
      where(cldmas_local > threshold) ! comparision is made in s-1 
         ! calculation is made in min-1
         flashrate = ratio_global*(cldmas_local(:,:)*60.-threshold*60.)*ratio_local(:,:)/60.
      elsewhere
         flashrate = 0.  
      end where

!Calculate N production rate (g s-1).
!For specified flashrate and prodfac, global production rate = 3.41 Tg N yr-1.
!Assume IC and CG flash rates produce the same amount of N / flash.
!      pnox = flashrate(:,:)*prodfac*convfac
      pnox = midLatAdj(:,:)*flashrate(:,:)*prodfac*convfac 
      
!adjust for desired global Nitrogen production rate
      pnox = pnox * desired_g_N_prod_rate /3.41

      return
      end subroutine flashfit

end module GmiEmissionLightning_mod
