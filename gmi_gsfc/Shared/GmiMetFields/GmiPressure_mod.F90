  module GmiPressure_mod

      implicit none

      private
      public  :: CalcPress3dCenter
      public  :: CalcPress3dEdge
      public  :: calcTropopausePress_Stobie
      public  :: calcTropopausePress_Ertel
      public  :: calcTropopausePress_WMO
      public  :: CalcAveragePressEdge

  contains

!-----------------------------------------------------------------------------
      subroutine CalcAveragePressEdge (ai, bi, pt, averagePressEdge, k1, k2)

      implicit none

#     include "gmi_phys_constants.h"

      integer, intent(in ) :: k1, k2
      real*8 , intent(in ) :: pt
      real*8 , intent(in ) :: ai(k1-1:k2)
      real*8 , intent(in ) :: bi(k1-1:k2)
      real*8 , intent(out) :: averagePressEdge(k1-1:k2)

      averagePressEdge(:) = (ai(:) * pt) + (bi(:) * AVG_SRFPRS)

      return

      end subroutine CalcAveragePressEdge
!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcPress3dCenter
!
! DESCRIPTION
!   This routine calculates the atmospheric pressure at the center of each
!   grid box assuming a hybrid sigma/pressure coordinate system.
!
!   pressure = (a * pt) + (b * psurf)
!
! ARGUMENTS
!   pt      : pressure = (am * pt) + (bm * psf) (mb)
!   am      : pressure = (am * pt) + (bm * psf), am at zone midpoint
!   bm      : pressure = (am * pt) + (bm * psf), bm at zone midpoint
!   psf     : surface pressure field at t1, known at zone centers (mb)
!   press3c : atmospheric pressure at the center of each grid box (mb)
!
!-----------------------------------------------------------------------------

      subroutine CalcPress3dCenter  &
     &  (pt, am, bm, psf, press3c, pr_diag, procID, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      real*8  :: pt
      real*8  :: am(k1:k2)
      real*8  :: bm(k1:k2)
      real*8  :: psf    (ilo:ihi, julo:jhi)
      real*8  :: press3c(ilo:ihi, julo:jhi, k1:k2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'CalcPress3dCenter called by ', procID
      end if


      do ik = k1, k2
        do ij = ju1, j2
          do il = i1, i2
              press3c(il,ij,ik) = (am(ik) * pt) + (bm(ik) * psf(il,ij))
          end do
        end do
      end do


      return

      end  subroutine CalcPress3dCenter

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcPress3dEdge
!
! DESCRIPTION
!   This routine calculates the atmospheric pressure at the edge of each
!   grid box assuming a hybrid sigma/pressure coordinate system.
!
!   pressure = (a * pt) + (b * psurf)
!
! ARGUMENTS
!   pt      : pressure = (ai * pt) + (bi * psf) (mb)
!   ai      : pressure = (ai * pt) + (bi * psf), ai at zone interface
!   bi      : pressure = (ai * pt) + (bi * psf), bi at zone interface
!   psf     : surface pressure field at t1, known at zone centers (mb)
!   press3e : atmospheric pressure at the edge of each grid box   (mb)
!
!-----------------------------------------------------------------------------

      subroutine CalcPress3dEdge  &
     &  (pt, ai, bi, psf, press3e, pr_diag, procID, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      real*8  :: pt
      real*8  :: ai(k1-1:k2)
      real*8  :: bi(k1-1:k2)
      real*8  :: psf    (ilo:ihi, julo:jhi)
      real*8  :: press3e(ilo:ihi, julo:jhi, k1-1:k2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'CalcPress3dEdge called by ', procID
      end if

      do ik = k1-1, k2
        do ij = ju1, j2
          do il = i1, i2
              press3e(il,ij,ik) = (ai(ik) * pt) + (bi(ik) * psf(il,ij))
          end do
        end do
      end do

      return

      end subroutine CalcPress3dEdge

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcTropopausePress_Stobie
!
! DESCRIPTION
!   This routine calculates the pressure at the tropopause using the Jim
!   Stobie algorithm.  That is, find, between 40 and 550 mb, the minimum
!   value of ALPHA (0.03) times the temperature minus log base 10 of the
!   pressure.
!
! ARGUMENTS
!   INPUT:
!     press3c : atmospheric pressure at the center of each grid box (mb)
!     kel     : temperature (degK)
!   OUTPUT:
!     tropp2d : 2D field of tropopause pressure in each column (mb)
!
!-----------------------------------------------------------------------------

      subroutine calcTropopausePress_Stobie  &
     &  (press3c, kel, tropp2d, pr_diag, procID, &
     &       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      real*8 :: kel    (ilo:ihi, julo:jhi, k1:k2)
      real*8 :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8 :: tropp2d(i1:i2, ju1:j2)


!     ----------------------
!     Parameter declarations.
!     ----------------------

      real*8, parameter :: ALPHA = 0.03d0


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical :: tropp_mask(i1:i2, ju1:j2, k1:k2)

      real*8  :: eqn    (i1:i2, ju1:j2, k1:k2)
      real*8  :: tropp3d(i1:i2, ju1:j2, k1:k2)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'calcTropopausePress_Stobie called by ', procID
      end if

      tropp3d(:,:,:) = 0.0d0

      eqn(:,:,:) = (ALPHA * kel  (i1:i2,ju1:j2,:)) -  &
     &             Log10 (press3c(i1:i2,ju1:j2,:))

      tropp_mask(:,:,:) =  &
     &  ((press3c(i1:i2,ju1:j2,:) >  40.0d0) .and.  &
     &   (press3c(i1:i2,ju1:j2,:) < 550.0d0))

      where (Spread  &
     &         (Minval (eqn, 3, tropp_mask(:,:,:)),  &
     &         3,  &
     &         k2) ==  &
     &       eqn(:,:,:))
        tropp3d(:,:,:) = press3c(i1:i2,ju1:j2,:)
      end where

      tropp2d(:,:) = Maxval (tropp3d(:,:,:), 3)

      return

      end subroutine calcTropopausePress_Stobie
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcTropopausePress_Ertel
!
! DESCRIPTION
!   This routine calculates the pressure at the tropopause using Ertel's 
!     Potential Vorticity (EPV) and Potential Temperature (theta). The tropopause is
!     determined by which has the greater pressure, the EPV surface equal to
!     EPV_LIMIT or the theta surface equal to THETA_LIMIT. The tropopause can
!     not be lower that the MINLEV level of the model and a lower pressure than
!     MINPRS.
!
! ARGUMENTS
!   INPUT:
!     press3c : atmospheric pressure at the center of each grid box (mb)
!     u       : zonal wind (m/s)
!     v       : meridional wind (m/s)
!     kel     : temperature (degK)
!     i's and j's : various model domain dimentions
!   OUTPUT:
!     trop_prs : 2D field of tropopause pressure in each column (hPa)
!     epv      : 3D field of Ertel's Potential Vorticity
!     theta    : 3D field of Potential Temperature (degK)
!
!-----------------------------------------------------------------------------

      subroutine calcTropopausePress_Ertel  &
     &  ( press3c, u, v, kel, trop_prs, epv, theta, &
     &    dlatr, coscen, cosp, londeg, latdeg, &
     &    pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &    i1_gl, i2_gl, ju1_gl, j2_gl )

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"

!
!... variables:
!... INPUT
!...  press3c - 3-D grid point pressure (hPa) consistant with the "kel" grid
!...  u - u-wind (m/s)
!...  v - v-wind (m/s)
!...  kel - temperature (K) - **ASSUMED UPSIDE DOWN WRT U AND V**
!... OUTPUT
!...  epv      - Ertel Potential Vorticity (3D)
!...  trop_prs - tropopause pressure (2D)
!
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID

      real*8, intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: kel    (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: u      (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: v      (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: londeg (i1_gl:i2_gl)
      real*8, intent(in) :: latdeg (ju1_gl:j2_gl)
      real*8, intent(in) :: coscen (ju1_gl:j2_gl)
      real*8, intent(in) :: cosp   (ju1_gl:j2_gl)
      real*8, intent(in) :: dlatr  (ju1_gl:j2_gl)

      real*8, intent(out) :: trop_prs(i1:i2, ju1:j2)
      real*8, intent(out) :: epv(i1:i2, ju1:j2, k1:k2)
      real*8, intent(out) :: theta(i1:i2, ju1:j2, k1:k2)

      integer :: i, j, k
      real*8  :: dudy, dvdx, thkp1, thkm1, dtpdp
      real*8  :: del_lat, del_lon

!... PARAMETERS
!... lowest level that can be stratosphere
      integer :: MINLEV=5
!... pressure limits for tropopause
      real*8 :: MINPRS=70.
!... threshold EPV value, higher than abs() => stratosphere
      real*8 :: EPV_LIMIT=3.6d-6
!... threshold potential temperature value, higher than value => stratosphere
      real*8 :: THETA_LIMIT=385.0d0

!... Angular velocity of rotation of Earth (rad/s)
      real*8 :: omega
!... Gas constant/specific heat for dry air (unitless)
      real*8 :: rocp = 297.0d0/1004

      real*8 :: tropp_epv(i1:i2, ju1:j2)
      real*8 :: tropp_pt(i1:i2, ju1:j2)
      real*8 :: f(ju1:j2)
      real*8 :: odacos(ju1:j2)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'calcTropopausePress_Ertel called by ', procID
      end if
  
!
!..................start of executable statements.......................
!... find lat in radians, coriolus force, cos(lat) and 1/(RADEAR*cos(lat))
      omega = 2.d0 * GMI_PI / SECPDY
      f(ju1:j2) = 2.d0 * omega * sin(DLATR(ju1:j2))

!      odacos(ju1:j2) = 1.0 / (RADEAR * COSP(ju1:j2))
!!... fix infinite problem at poles of 1/(ae*cos(lat))
!      if (ju1 .eq. julo)  &
!     &  odacos(julo) = odacos(julo+1) - (odacos(julo+2)-odacos(julo+1)) * 0.5
!      if (ju1 .eq. jhi)  &
!     &  odacos(jhi)  = odacos(jhi-1)  - (odacos(jhi-2) -odacos(jhi-1))  * 0.5

      odacos(ju1+1:j2-1) = 1.0 / (RADEAR * COSP(ju1+1:j2-1))
!... do end points with special cases for poles (extrapolation)
      if (cosp(ju1) == 0.0) then
         odacos(ju1) = odacos(ju1+1) - (odacos(ju1+2)-odacos(ju1+1)) * 0.5
      else
         odacos(ju1) = 1.0 / (RADEAR * COSP(ju1))
      endif
      if (cosp(j2) == 0.0) then
         odacos(j2)  = odacos(j2-1)  - (odacos(j2-2) -odacos(j2-1))  * 0.5
      else
         odacos(j2) = 1.0 / (RADEAR * COSP(j2))
      endif
!... or possibly use this?
!      if (ju1 .eq. ju1_gl) then
!         odacos(julo) = odacos(julo+1) - (odacos(julo+2)-odacos(julo+1)) * 0.5
!       else
!         odacos(ju1) = 1.0 / (RADEAR * COSP(ju1))
!       endif
!      if (j2 .eq. j2_gl) then
!         odacos(jhi)  = odacos(jhi-1)  - (odacos(jhi-2) -odacos(jhi-1))  * 0.5
!       else
!         odacos(j2) = 1.0 / (RADEAR * COSP(j2))
!       endif

!... calculate potential temperature (3D)
      theta(:,:,:) = kel(i1:i2,ju1:j2,k1:k2) &
     &       * ( (1000./press3c(i1:i2,ju1:j2,k1:k2))**rocp )

!... assuming regular grid in domain, we need 2*dlon in radians
      del_lon = 2.*(londeg(i1+1)-londeg(i1))*RADPDEG
!... assuming regular grid in domain, we need 2*dlat in radians
      del_lat = 2.*(latdeg(ju1+1)-latdeg(ju1))*RADPDEG
!... calculate EPV and find tropopause pressure
      do k=k1,k2,1
        do j=ju1,j2,1
          do i=i1,i2,1
!... calc du/d(lat) (centered differences)
            if (j .eq. ju1_gl) then
                dudy = ( (u(i,j+1,k) - u(i,j-1,k))*cos(latdeg(j+1)*RADPDEG) ) &
     &            / del_lat
              elseif( j .eq. j2_gl) then
                dudy = ( (u(i,j+1,k) - u(i,j-1,k) )*cos(latdeg(j-1)*RADPDEG) ) &
     &            / del_lat
              else
                dudy = (u(i,j+1,k)*coscen(j+1) - u(i,j-1,k)*coscen(j-1)) &
     &            / del_lat
            endif

!... calc dv/d(long) (centered differences)
            if(i .eq. i1_gl) then
                dvdx = (v(i+1,j,k) - v(i-1,j,k)) / del_lon
              elseif(i .eq. i2_gl) then
                dvdx = (v(i+1,j,k) - v(i-1,j,k)) / del_lon
              else
                dvdx = (v(i+1,j,k) - v(i-1,j,k)) / del_lon
            endif

!... calc d(theta)/dp (centered differences)
            if(k .eq. k1) then
                thkp1 = theta(i,j,k+1)
                thkm1 = theta(i,j,k)
                dtpdp = 0.01*((thkp1-thkm1)/(press3c(i,j,k+1)-press3c(i,j,k)))
              elseif(k .eq. k2) then
                thkp1 = theta(i,j,k)
                thkm1 = theta(i,j,k-1)
                dtpdp = 0.01*((thkp1-thkm1)/(press3c(i,j,k)-press3c(i,j,k-1)))
              else
                thkp1 = theta(i,j,k+1)
                thkm1 = theta(i,j,k-1)
                dtpdp = 0.01*((thkp1-thkm1)/(press3c(i,j,k+1)-press3c(i,j,k-1)))
            endif

!... calc EPV
            epv(i,j,k) = -GMI_G*((dvdx-dudy)*odacos(j)+f(j))*dtpdp

          enddo
        enddo
      enddo

!... get tropopause pressure
      tropp_epv = 0.0
      tropp_pt = 0.0
      do k=MINLEV,k2,1
        do j=ju1,j2,1
          do i=i1,i2,1

!... make sure we don't go too high in atm
            if (MINPRS .le. press3c(i,j,k)) then

!... calculate the EPV based tropopause pressure
              if(abs(epv(i,j,k)) .gt. EPV_LIMIT &
     &           .and. abs(epv(i,j,k-1)) .le. EPV_LIMIT) then

                tropp_epv(i,j) = press3c(i,j,k) +          &
     &            ((press3c(i,j,k-1)  - press3c(i,j,k))  *  &
     &            (EPV_LIMIT         - abs(epv(i,j,k)))) /  &
     &            (abs(epv(i,j,k-1)) - abs(epv(i,j,k)))
              endif

!... calculate the Potential Temp based tropopause pressure
              if(theta(i,j,k) .gt. THETA_LIMIT &
     &           .and. theta(i,j,k-1) .le. THETA_LIMIT) then

                tropp_pt(i,j) = press3c(i,j,k) +          &
     &            ((press3c(i,j,k-1)  - press3c(i,j,k)) *  &
     &            (THETA_LIMIT    - theta(i,j,k))) /     &
     &            (theta(i,j,k-1) - theta(i,j,k))
              endif

            endif

          enddo
        enddo
      enddo

!... choose higher pressure for tropopause pressure
      do j=ju1,j2,1
        do i=i1,i2,1
          trop_prs(i,j) = max(tropp_pt(i,j),tropp_epv(i,j))
        enddo
      enddo

      return

      end subroutine calcTropopausePress_Ertel


!-------------------------------------------------------------------------
      SUBROUTINE calcTropopausePress_WMO  &
     &  (metdata_name_org, metdata_name_model, &
     &   kel, press3c, grid_height, latdeg, &
!    &   ltpause, htpause, &
     &   ptpause, &
     &   pr_diag, procID, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   ilo, ihi, julo, jhi, ju1_gl, j2_gl)
!
!******************************************************************************
!  Subroutine calcTropopausePress_WMO defines the tropopause layer in terms of temperature
!  lapse rates. (hyl, bmy, 11/30/99, 6/18/03; hyl, 11/30/07)
!******************************************************************************

      IMPLICIT NONE

      character (len=16) ,intent(in) :: metdata_name_org
      character (len=16) ,intent(in) :: metdata_name_model
      integer, intent(in) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2

      real*8 , intent(in   ) :: kel    (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: grid_height(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: procID
!
!     integer, intent(out  ) :: ltpause(i1:i2, ju1:j2)
!     real*8,  intent(out  ) :: htpause(i1:i2, ju1:j2)
      real*8,  intent(out  ) :: ptpause(i1:i2, ju1:j2)
      integer :: ltpause(i1:i2, ju1:j2)
      real*8  :: htpause(i1:i2, ju1:j2)

      ! Local variables
      INTEGER :: I, J, L, K
      REAL*8  :: HTPLIMIT, Y
      REAL*8  :: H(i1:i2, ju1:j2, k1:k2)
      REAL*8  :: LAPSE_R(k1:k2)
      REAL*8  :: LAPSE_T(k1:k2)

      if (pr_diag) then
        Write (6,*) 'calcTropopausePress_WMO called by ', procID
      end if

      !=================================================================
      ! calcTropopausePress_WMO begins here!
      !
      ! H (in m) is the height of the midpoint of layer L (hyl, 03/28/99)
      !=================================================================
      ! Find height of the midpoint of the first level
      DO J = ju1, j2
      DO I = i1, i2
         H(I,J,1) = grid_height(I,J,1) / 2.d0
      ENDDO
      ENDDO

      ! Add to H 1/2 of the sum of the two adjacent boxheights
      DO L = k1, k2-1
      DO J = ju1, j2
      DO I = i1, i2
         H(I,J,L+1) = H(I,J,L) +  &
     &              ( grid_height(I,J,L) + grid_height(I,J,L+1) ) / 2.d0
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! Initialize ltpause and htpause arrays first.
      !
      ! ltpause  = the tropopause layer #
      ! htpause  = the tropopause height [ km ]
      ! HTPLIMIT = maximum tropopause height [ m ]
      !
      ! We need the factor of 1d3 to convert HTPLIMIT from km --> m.
      !=================================================================

      DO J = ju1, j2

         ! Latitude [degrees]
           Y = latdeg( J )

      DO I = i1,i2 
         ltpause(I,J) = 0
         htpause(I,J) = 0.0d0

         ! In the tropics, the tropopause maxes out at 7 km,
         ! elsewhere the tropopause maxes out at 5 km
         IF ( ABS( Y ) <= 30.d0 ) THEN
            HTPLIMIT = 7.0d3
         ELSE
            HTPLIMIT = 5.0d3
         ENDIF

         !==============================================================
         ! Tropopause: 15-20km (equator); 8-12km (temperate latitudes
         ! and poles).  According to WMO, the "1st tropopause" is
         ! defined as the lowest level at which the lapse rate decreases
         ! to 2 C/km or less, provided also the average lapse rate
         ! between this level and all higher levels within 2 km doesn't
         ! exceed 2 C/km. It is noted that average lapse rate is the
         ! difference between the temperatures at the respective end
         ! points divided by the height interval irrespective of lapse
         ! rate variations in the layer between the end points
         ! (hyl, 03/28/99).
         !
         ! NOTE: 2 C/km is equivalent to 2.0d-3 C/m.  We have to keep
         !       this in C/m since H is in meters, and therefore, the
         !       lapse rate will have units of K/m (= C/m).
         !       (hyl, bmy, 11/30/99)
         !==============================================================
          DO L = 2, k2-1
            IF ( H(I,J,L) >= HTPLIMIT ) THEN

               ! Lapse rate in the current level L
               LAPSE_R(L) = ( kel(I,J,L  ) - kel(I,J,L+1) ) /  &
     &                      ( H(I,J,L+1) - H(I,J,L  ) )

               ! Test for lapse rate in the Lth level < 2 C/km
               IF ( LAPSE_R(L) <= 2.0d-3 ) THEN
                  K = 1

                  ! Compute the lapse rate at level L+K
 10               LAPSE_T(K) = ( kel(I,J,L)   - kel(I,J,L+K) ) /  &
     &                         ( H(I,J,L+K) - H(I,J,L)   )

                  ! If the lapse rate at L+K is still less than 2 C/km,
                  ! and level L+K is within 2 km of level K, we have found
                  ! the tropopause!  Go to the next surface box.
                  IF ( ( H(I,J,L+K) - H(I,J,L) ) > 2.0d3 .and.  &
     &                 LAPSE_T(K) < 2.0d-3 ) THEN
                     ltpause(I,J) = L
                     htpause(I,J) = H(I,J,L) / 1.0d3 ! m --> km
                     GOTO 30
                  ENDIF

                  ! If the lapse rate at L+K is greater than 2 C/km
                  ! then we are not high enough.  Go to the next level L.
                  IF ( ( H(I,J,L+K) - H(I,J,L) ) > 2.0d3 .and.  &
     &                 LAPSE_T(K) > 2.0d-3 ) GOTO 20

                  ! If level L+K is not within 2km of level K, then
                  ! we are not high enough.  Go to the next level L.
                  IF ( ( H(I,J,L+K) - H(I,J,L) ) < 2.0d3 .and.  &
     &                 LAPSE_T(K) > 2.0d-3 ) GOTO 20

                  ! Increment level K until the lapse rate at L+K < 2 C/km
                  K = K + 1

!-----------------------------------------------------------------------------
!          IF ( K>LM )
!             LAPSE_T(K) = (T(I,J,L)-T(I,J,L+K)) / (H(I,J,L+K)-H(I,J,L))
!          IF((H(I,J,L+K)-H(I,J,L)<2000 .and. LAPSE_T(K)>0.002) GOTO 20
!-----------------------------------------------------------------------------
                  ! If none of the above conditions were met, then
                  ! test here to make sure we don't violate array bounds
                  !---------------------------------------
                  ! Prior to 7/20/04:
                  !IF ( ( L + K ) <= LM ) THEN
                  !---------------------------------------
                   IF ( ( L + K ) <= k2 ) THEN

                     ! If Level L+K is within 2 km of level K,
                     ! go back and compute the lapse rate
                     IF ( ( H(I,J,L+K) - H(I,J,L) ) < 2.0d3 ) THEN
                        GOTO 10

                     ! Otherwise we have found the tropopause level
                     ELSE
                        ltpause (I,J) = L
                        htpause (I,J) = H(I,J,L) / 1.0d3 ! m --> km
                        GOTO 30
                     ENDIF
                  ELSE

                     ! If level L+K is higher than LM, then we are
                     ! going out of array bounds, so call L=LM
                     ! the tropopause level
                     ltpause (I,J) = L
                     htpause (I,J) = H(I,J,L) / 1.0d3 ! m --> km
                     GOTO 30
                  ENDIF
               ENDIF
            ENDIF
 20      ENDDO      ! L -- go to next level
 30   ENDDO         ! I -- go to next grid box
      ENDDO         ! J

      !=================================================================
      ! Sometimes a tropopause cannot be located in terms of the
      ! above definition. For the time being, set it to that of
      ! the nearest grid.
      !=================================================================
      DO J = ju1, j2
      DO I = i1, i2 
         IF ( ltpause(I,J) == 0 ) THEN

            IF ( I /= 1 ) THEN
               ltpause(I,J) = ltpause(I-1,J)
               IF ( ltpause(I,J) == 0 ) THEN
                  WRITE(6,*) 'ltpause = 0 in calcTropopausePress_WMO!' 
                  WRITE(6,*) 'I,J= ', I, J
                  STOP
               ENDIF
               htpause(I,J) = H(I,J,ltpause(I,J)) / 1.0d3

            ELSE IF ( J /= 1 ) THEN
               ltpause(I,J) = ltpause(I,J-1)
               IF ( ltpause(I,J) == 0 ) THEN
                  WRITE(6,*) 'ltpause = 0 in calcTropopausePress_WMO!' 
                  WRITE(6,*) 'I,J= ', I, J
                  STOP
               ENDIF
               htpause(I,J) = H(I,J,ltpause(I,J)) / 1.0d3

            ! South polar boxes
            ELSE IF ( J == 1 ) THEN

             ! Select the proper polar tropopause level
             ! for DAO, GISS II', fvGCM or GEOS4-DAS (hyl, 11/30/07)
             if (metdata_name_org(1:3) .eq. 'DAO') then
               ltpause(1,1) = 11
             endif
             if (metdata_name_org(1:4) .eq. 'GISS') then
               ltpause(1,1) = 9
             endif
             if ((metdata_name_org(1:4) == 'GMAO') ) then
                if (metdata_name_model(1:5) == 'GEOS4') then
                  ltpause(1,1) = 11    !fvGCM or GEOS4-DAS
                elseif (metdata_name_model(1:5) == 'GEOS5') then
                  ltpause(1,1) = 28    !GEOS5-DAS
                else
                  write(6,*) 'stop in calcTropopausePress_WMO!'
                endif
             endif
             if ((metdata_name_model(1:3) == 'AM3') ) then
               ltpause(1,1) = 19       !277 hPa
             endif

               htpause(1,1) = H(I,J,ltpause(1,1)) / 1.0d3

            ENDIF
         ENDIF

         !--------------------------------------------------------------------
         ! Debug output...check if ltpause is 0 (hyl, 11/30/99)
         IF ( ltpause (I,J) <= 0 ) THEN
            write(6,*) 'ltpause(I,J)= ', ltpause (I,J),I,J
            stop 'ltpause <= 0 in calcTropopausePress_WMO !'
         ENDIF
         !--------------------------------------------------------------------

         ptpause (I, J) = press3c ( I, J, ltpause (I,J) )

      ENDDO  ! I
      ENDDO  ! J

      ! debug (hyl,12/3/07)
      !write(6,*) 'ltpause(:,:)= ', ltpause(:,:)
      !write(6,*) 'MinVal(ltpause)= ', MinVal(ltpause)
      !write(6,*) 'MaxVal(ltpause)= ', MaxVal(ltpause)

      ! Return to calling program
      END SUBROUTINE calcTropopausePress_WMO
!-------------------------------------------------------------------------

  end module GmiPressure_mod
