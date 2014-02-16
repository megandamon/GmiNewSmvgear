module GmiSolar_mod

      implicit none

      private
      public  :: CalcSolarDeclination
!.      public  :: CalcSolarZenithAngle
      public  :: CalcCosSolarZenithAngle
      public  :: computeSolarZenithAngle_Photolysis
      public  :: computeSolarZenithAngle_Photolysis_2

      contains

!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Keith Grant, LLNL)
!   jrt@llnl.gov
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcSolarDeclination
!
! DESCRIPTION
!   Given the Julian day, this routine calculates the solar declination and
!   the square of the ratio of the mean Earth-Sun distance to the current
!   Earth-Sun distance.   Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, pp. 57-63.
!
! ARGUMENTS
!   julday  : Julian day counting ranging from 0 on Jan. 1st to 364 on
!             Dec. 31st
!   decl    : solar declination (deg)
!   rdistsq : the square of the ratio of the mean Earth-Sun distance
!             to the current Earth-Sun distance
!
!-----------------------------------------------------------------------------

      subroutine CalcSolarDeclination  &
     &  (julday, decl, rdistsq)

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8, intent(in )  :: julday
      real*8, intent(out)  :: decl
      real*8, intent(out)  :: rdistsq

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: TWO_PI = 2.0d0 * GMI_PI

!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: theta           ! year angle on specified day
      real*8  :: theta2, theta3

!     ----------------
!     Begin execution.
!     ----------------

!c?   Leap year?
      theta  = (TWO_PI * julday) / 365.0d0

      theta2 = 2.0d0 * theta
      theta3 = 3.0d0 * theta

      decl =  &
     &  (360.0d0 / TWO_PI) *  &
     &  (0.006918d0 -  &
     &   0.399912d0 * Cos (theta)  + 0.070257d0 * Sin (theta)  -  &
     &   0.006758d0 * Cos (theta2) + 0.000907d0 * Sin (theta2) -  &
     &   0.002697d0 * Cos (theta3) + 0.001480d0 * Sin (theta3))

      rdistsq =  &
     &  1.000110d0  +  &
     &  0.034221d0  * Cos (theta)  + 0.001280d0 * Sin (theta)  +  &
     &  0.000719d0  * Cos (theta2) + 0.000077d0 * Sin (theta2)

      return

      end subroutine CalcSolarDeclination

!.!-----------------------------------------------------------------------------
!.!
!.! ROUTINE
!.!   Solrza
!.!
!.! DESCRIPTION
!.!   Given a Greenwich time, a solar declination and reference latitude and
!.!   longitudes,  this routine returns a corresponding list of cosines of solar
!.!   zenith angles.  Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!.!   "Radiative Processes in Meteorology and Climatology", Elsevier, p. 62.
!.!
!.! ARGUMENTS
!.!   time   : Greenwich time since Jan 1, counting zero from midnight (days)
!.!   decl   : solar declination (deg)
!.!   lat    : latitude   (deg)
!.!   lon    : longitudes (deg)
!.!   nn     : number of longitudes for which to calculate cosines of the
!.!            solar zenith angle
!.!   cossza : cosines of the solar zenith angle (output)
!.!
!.!-----------------------------------------------------------------------------
!.
!.      subroutine CalcSolarZenithAngle  &
!.     &  (time, decl, lat, lon, nn, cossza)
!.
!.      implicit none
!.
!.#     include "gmi_phys_constants.h"
!.
!.!     ----------------------
!.!     Argument declarations.
!.!     ----------------------
!.
!.      integer, intent(in ) :: nn
!.      real*8 , intent(in ) :: time
!.      real*8 , intent(in ) :: decl
!.      real*8 , intent(in ) :: lat
!.      real*8 , intent(in ) :: lon   (nn)
!.      real*8 , intent(out) :: cossza(nn)
!.
!.!     -----------------------
!.!     Parameter declarations.
!.!     -----------------------
!.
!.      real*8,  parameter :: TWO_PI = 2.0d0 * GMI_PI
!.      real*8,  parameter :: PI_180 = TWO_PI / 360.0d0
!.
!.!     ----------------------
!.!     Variable declarations.
!.!     ----------------------
!.
!.      integer :: ii
!.
!.      real*8  :: cosha
!.
!.!     ----------------
!.!     Begin execution.
!.!     ----------------
!.
!.      do ii = 1, nn
!.
!.!       -------------------------------------------------------------
!.!       Calculate the cosine of the hour angle which is referenced to
!.!       noon.
!.!       -------------------------------------------------------------
!.
!.        cosha =  &
!.     &    Cos (TWO_PI *  &
!.     &         (Mod (time + (lon(ii) / 360.d0), 1.0d0) - 0.5d0))
!.
!.
!.!       -----------------------------------------------
!.!       Calculate the cosine of the solar zenith angle.
!.!       -----------------------------------------------
!.
!.        cossza(ii) =  &
!.     &    Sin (PI_180 * lat) * Sin (PI_180 * decl) +  &
!.     &    Cos (PI_180 * lat) * Cos (PI_180 * decl) * cosha
!.
!.      end do
!.
!.
!.      return
!.
!.      end subroutine CalcSolarZenithAngle
!.
!-----------------------------------------------------------------------------
!
! ROUTINE
!  CalcCosSolarZenithAngle 
!
! DESCRIPTION
!   Given a Greenwich time, a solar declination and reference latitude and
!   longitudes,  this routine returns a corresponding list of cosines of solar
!   zenith angles.  Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, p. 62.
!
! ARGUMENTS
!   time   : Greenwich time since Jan 1, counting zero from midnight (days)
!   lat    : latitude   (deg)
!   lon    : longitudes (deg)
!   nn     : number of longitudes for which to calculate cosines of the
!            solar zenith angle
!   cossza : cosines of the solar zenith angle (output)
!
!-----------------------------------------------------------------------------

      subroutine CalcCosSolarZenithAngle  &
     &  (time, latDeg, lonDeg, cosSolarZenithAngle, &
     &   i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in ) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
      real*8 , intent(in ) :: time
      real*8 , intent(in ) :: latDeg(ju1_gl:j2_gl)
      real*8 , intent(in ) :: lonDeg(i1_gl:i2_gl) 
      real*8 , intent(out) :: cosSolarZenithAngle(i1:i2,ju1:j2)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8,  parameter :: TWO_PI = 2.0d0 * GMI_PI
      real*8,  parameter :: PI_180 = TWO_PI / 360.0d0

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ij
      real*8  :: decl, theta, theta2, theta3
      real*8  :: cosha(i1:i2), coslat, sinlat

!     ----------------
!     Begin execution.
!     ----------------

!c?   Leap year?
      theta  = (TWO_PI * time) / 365.0d0

      theta2 = 2.0d0 * theta
      theta3 = 3.0d0 * theta

!     Calculate the solar declination (deg)

      decl =  &
     &  (360.0d0 / TWO_PI) *  &
     &  (0.006918d0 -  &
     &   0.399912d0 * Cos (theta)  + 0.070257d0 * Sin (theta)  -  &
     &   0.006758d0 * Cos (theta2) + 0.000907d0 * Sin (theta2) -  &
     &   0.002697d0 * Cos (theta3) + 0.001480d0 * Sin (theta3))

!     -------------------------------------------------------------------
!     Calculate the cosine of the hour angle which is referenced to noon.
!     -------------------------------------------------------------------

      cosha (:) =  Cos(TWO_PI * (Mod (time + (lonDeg(i1:i2) / 360.d0), 1.0d0) - 0.5d0))

!     -----------------------------------------------
!     Calculate the cosine of the solar zenith angle.
!     -----------------------------------------------

      do ij = ju1, j2
        sinlat =  Sin (PI_180 * latDeg(ij)) * Sin (PI_180 * decl)
        coslat =  Cos (PI_180 * latDeg(ij)) * Cos (PI_180 * decl)

        cosSolarZenithAngle(:,ij) =  sinlat + coslat*cosha(:)
      end do

      return

      end subroutine CalcCosSolarZenithAngle
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: computeSolarZenithAngle_Photolysis
!
! !INTERFACE:
!
      function computeSolarZenithAngle_Photolysis (iday, time_sec, &
                   fastj_offset_sec, latDeg, lonDeg, i1, i2, j1, j2) &
                   result(this_)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, j1, j2
      integer, intent(in) :: iday             ! day of year (1-365)
      real*8 , intent(in) :: time_sec         ! current time in seconds
      real*8 , intent(in) :: fastj_offset_sec ! offset from tau at which to do photolysis (s)
      real*8 , intent(in) :: latDeg(j1:j2)    ! latitude (deg)
      real*8 , intent(in) :: lonDeg(i1:i2)    ! longitude (deg)
!
! !RETURNED VALUE
      real*8  :: this_(i1:i2,j1:j2)
!
! !DESCRIPTION:
!  Computes the solar zenith angle used in the photolysis package.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: locPI      = 3.141592653589793D0
      real*8, parameter :: Deg2Rad    = locPI / 180.0d0
      real*8, parameter :: secPerHour = 3600.0d0
!
! !LOCAL VARIABLES:
    REAL*8  :: sindec, soldek, cosdec
    REAL*8  :: sinlat, sollat, coslat
    REAL*8  :: cosz, latRad, lonRad
    real*8  :: tau, timej, loct
    integer :: ii, jj
!EOP
!------------------------------------------------------------------------------
!BOC
      tau   = time_sec         / secPerHour
      timej = fastj_offset_sec / secPerHour

      sindec = 0.3978d0*sin(0.9863d0*(dble(iday)-80.d0)*Deg2Rad)
      soldek = asin(sindec)
      cosdec = cos(soldek)

      do jj = j1, j2
         sinlat = sin(latDeg(jj)*Deg2Rad)
         sollat = asin(sinlat)
         coslat = cos(sollat)

         do ii = i1, i2
            loct = (((tau+timej)*15.d0)-180.d0)*Deg2Rad + lonDeg(ii)*Deg2Rad
            cosz = cosdec*coslat*cos(loct) + sindec*sinlat
            this_(ii,jj) = acos(cosz)/Deg2Rad
         enddo
      enddo


      end function computeSolarZenithAngle_Photolysis
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: computeSolarZenithAngle_Photolysis_2
!
! !INTERFACE:
!
      function computeSolarZenithAngle_Photolysis_2 (iday, time_sec, &
                   fastj_offset_sec, lat, lon) 
!
! !INPUT PARAMETERS:
      integer, intent(in) :: iday             ! day of year (1-365)
      real*8 , intent(in) :: time_sec         ! current time in seconds
      real*8 , intent(in) :: fastj_offset_sec ! offset from tau at which to do photolysis (s)
      real*8 , intent(in) :: lat    ! latitude (deg)
      real*8 , intent(in) :: lon    ! longitude (deg)
!
! !RETURNED VALUE
      real*8  :: computeSolarZenithAngle_Photolysis_2
!
! !DESCRIPTION:
!  Computes the solar zenith angle used in the photolysis package.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: locPI      = 3.141592653589793D0
      real*8, parameter :: pi180    = locPI / 180.0d0
      real*8, parameter :: secPerHour = 3600.0d0
!
! !LOCAL VARIABLES:
    REAL*8  :: sindec, soldek, cosdec
    REAL*8  :: sinlat, sollat, coslat
    REAL*8  :: cosz
    real*8  :: tau, timej, loct
    integer :: ii, jj
!EOP
!------------------------------------------------------------------------------
!BOC
      tau   = time_sec         / secPerHour
      timej = fastj_offset_sec / secPerHour

      sindec=0.3978d0*sin(0.9863d0*(dble(iday)-80.d0)*pi180)
      soldek=asin(sindec)
      cosdec=cos(soldek)

      sinlat = sin(lat*pi180)
      sollat=asin(sinlat)
      coslat=cos(sollat)
!
      loct = (((tau+timej)*15.d0)-180.d0)*pi180 + lon*pi180
      cosz = cosdec*coslat*cos(loct) + sindec*sinlat
      computeSolarZenithAngle_Photolysis_2 = acos(cosz)/pi180

      end function computeSolarZenithAngle_Photolysis_2
!EOC
!------------------------------------------------------------------------------

end module GmiSolar_mod
