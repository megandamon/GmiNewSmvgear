!=======================================================================
!
! $Id: setkin_group_specifics.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_group_specifics.h
!
! DESCRIPTION
!   This include file contains information about grouping species for
!   transport purposes, that is specific to the mechanism. In general,
!   this file can not be generated automatically, but must be manually
!   produced.
!
!  Chemistry input file:    GMIS2 4:59 PM 7/30/2002
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Mon Jan 13 22:23:01 2003
!
!=======================================================================

      qq2(:,:,:)                   = 0.0d0

      max_bry_adjust               = 0.0d0
      max_cly_adjust               = 0.0d0
      max_nox_adjust               = 0.0d0

      select case (ig)
        case (1)

!.... Bry

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        qq2(:,:,:)           = qq2(:,:,:) +  &
     &                               const(:,:,:,imsgrp) *  &
     &                               sgrp_fac(im,ig)
          end do

          group_factor(:,:,:)      = qq1(i1:i2,ju1:j2,:) / qq2(:,:,:)

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        const(:,:,:,imsgrp)  = const(:,:,:,imsgrp) *  &
     &                              group_factor(:,:,:)
          end do

        case (2)

!.... Cly

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(:,:,:)             = qq2(:,:,:) +  &
     &                               const(:,:,:,imsgrp) *  &
     &                               sgrp_fac(im,ig)
        end do

        group_factor(:,:,:)        = qq1(i1:i2,ju1:j2,:) / qq2(:,:,:)

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)

!.... Exclude BrCl (#32) from adjustment

          if ((imsgrp > 0) .and. (imsgrp /= 32))  &
     &      const(:,:,:,imsgrp)    = const(:,:,:,imsgrp) *  &
     &                               group_factor(:,:,:)
        end do

!.... Adjust Cl2 (#24) in place of BrCl (#32)

        group_adjust(:,:,:)        = (group_factor(:,:,:) - 1.0d0) *  &
     &                               0.5d0 * const(:,:,:,32)

        const(:,:,:,24)            = const(:,:,:,24) +  &
     &                               group_adjust(:,:,:)

        if (Any(const(:,:,:,24) < 0.0d0)) then

!.... In locations where Cl2 (#24) is too small to absorb the entire BrCl
!.... decrement, extend it to OClO (#26)

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,24) < 0.0d0)
            group_adjust(:,:,:)    = 2.0d0 * const(:,:,:,24)
            const(:,:,:,24)        = 0.0d0
            const(:,:,:,26)        = const(:,:,:,26) +  &
     &                               group_adjust(:,:,:)
          end where

        end if

        if (Any(const(:,:,:,26) < 0.0d0)) then

!.... In locations where OClO (#26) is too small to absorb the BrCl
!.... residual decrement,extend it to Cl2O2 (#27)

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,26) < 0.0d0)
            group_adjust(:,:,:)    = 0.5d0 * const(:,:,:,26)
            const(:,:,:,26)        = 0.0d0
            const(:,:,:,27)        = const(:,:,:,27) +  &
     &                               group_adjust(:,:,:)
          end where

        end if

        if (Any(const(:,:,:,27) < 0.0d0)) then

!.... In locations where Cl2O2 (#27) is too small to absorb the BrCl
!.... residual decrement,extend it to ClO (#25)

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,27) < 0.0d0)
            group_adjust(:,:,:)    = 2.0d0 * const(:,:,:,27)
            const(:,:,:,27)        = 0.0d0
            const(:,:,:,25)        = const(:,:,:,25) +  &
     &                               group_adjust(:,:,:)
          end where

        end if

        if (Any(const(:,:,:,25) < 0.0d0)) then

!.... In locations where ClO (#25) is too small to absorb the BrCl
!.... residual decrement,extend it to HCl (#28)

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,25) < 0.0d0)
            group_adjust(:,:,:)    = const(:,:,:,25)
            const(:,:,:,25)        = 0.0d0
            const(:,:,:,28)        = const(:,:,:,28) +  &
     &                               group_adjust(:,:,:)
          end where

          if (Any(const(:,:,:,28) < 0.0d0)) then
            max_cly_adjust         = min(max_cly_adjust,  &
     &                                   Minval(group_factor - 1.0d0,  &
     &                                   MASK=((group_factor - 1.0d0)  &
     &                                         < 0.0d0)))
            Write (6,*) 'BrCl adjustment drove HCl negative',  &
     &        max_cly_adjust, loc_proc
            Write (6,*) Minloc(const(:,:,:,28))
          end if

        end if

!.... Account for changes in NOy reservoir ClONO2 (#30)

        group_adjust(:,:,:)       = const(:,:,:,30) *  &
     &                              (1.0d0 - group_factor(:,:,:)) /  &
     &                              group_factor(:,:,:)

!.... Park difference in N2O5 (#9)

        const(:,:,:,9)            = const(:,:,:,9) +  &
     &                              0.5d0 * group_adjust(:,:,:)

        if (Any(const(:,:,:,9) < 0.0d0)) then

!.... Park remainder in NOx and then HNO3 if necessary

          group_adjust(:,:,:)     = 0.0d0
          Where (const(:,:,:,9) < 0.0d0)
            group_adjust(:,:,:)   = 2.0d0 * const(:,:,:,9)
            const(:,:,:,9)        = 0.0d0
            const(:,:,:,6)        = const(:,:,:,6) +  &
     &                              group_adjust(:,:,:) *  &
     &                              (const(:,:,:,6) /  &
     &                              (const(:,:,:,6) + const(:,:,:,7)))
            const(:,:,:,7)        = const(:,:,:,7) +  &
     &                              group_adjust(:,:,:) *  &
     &                              (const(:,:,:,7) /  &
     &                              (const(:,:,:,6) + const(:,:,:,7)))
          end where

          if (Any(const(:,:,:,6) < 0.0d0)) then

              group_adjust(:,:,:)   = 0.0d0
            Where (const(:,:,:,6) < 0.0d0)
              group_adjust(:,:,:) = const(:,:,:,6)
              const(:,:,:,6)      = 0.0d0
              const(:,:,:,10)     = const(:,:,:,10) +  &
     &                              group_adjust(:,:,:)
            end where

          endif

          if (Any(const(:,:,:,7) < 0.0d0)) then

              group_adjust(:,:,:)   = 0.0d0
            Where (const(:,:,:,7) < 0.0d0)
              group_adjust(:,:,:) = const(:,:,:,7)
              const(:,:,:,7)      = 0.0d0
              const(:,:,:,10)     = const(:,:,:,10) +  &
     &                              group_adjust(:,:,:)
            end where

          endif

          if (Any(const(:,:,:,10) < 0.0d0)) then
            max_nox_adjust        = min(max_nox_adjust,  &
     &                                  minval(group_adjust(:,:,:)))
            Write (6,*) 'ClONO2 adjustment drove NOy negative',  &
     &        max_nox_adjust, loc_proc
          end if

          end if

      end select

