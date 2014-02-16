!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlASCII_mod
!
      module GmiControlASCII_mod
!
! !USES:
      use GmiReduce_mod, only : Gmi_Maxpij_Reduce, Gmi_Minpij_Reduce
      use GmiReduce_mod, only : Gmi_Sum_Reduce
      use  GmiSub2Glob_mod, only : subDomain2Global
      use GmiASCIIoperations_mod, only : AsciiOpenWrite
      use GmiFileUnit_mod, only : ReleaseFileUnitNumber
      use GmiFileOperations_mod, only : makeOutfileName
      use GmiFlush_mod          , only : GmiFlush
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiDiagnosticsMethod_mod, only :  Get_pr_ascii, Get_pr_ascii1,       &
     &       Get_pr_ascii2, Get_pr_ascii3, Get_pr_ascii4, Get_pr_ascii5,       &
     &       Get_ascii_out_i, Get_ascii_out_n, Get_asclun, Get_problem_name,   &
     &       Get_pr_ascii_step_interval, Get_pr_diag, t_Diagnostics
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID, Get_londeg
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies,  Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl,           &
     &       Get_ju1_gl, Get_j2_gl
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID, Get_iAmRootProc,          &
     &       Get_numDomains, Get_communicatorWorld, Get_map1_u
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       Get_concentration, Get_prod, Get_loss, Get_fixed_const_map,       &
     &       Get_num_fixed_const
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_const_labels,        &
     &      Get_num_molefrac, Get_num_chem, Get_mw
      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime,             &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep,  &
     &       Get_numTimeSteps
      use GmiMetFieldsControl_mod, only: t_metFields, Get_mass
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: initializeOutputASCII, controlOutputASCII, finalizeOutputASCII
!
#     include "gmi_diag_constants_llnl.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
!
      logical, save :: pr_ascii, pr_ascii1, pr_ascii2, pr_ascii3 
      logical, save :: pr_ascii4, pr_ascii5 , pr_diag
      integer, save :: ascii_out_i, ascii_out_n, asclun
      integer, save :: pr_ascii_step_interval  
      integer, save :: numSpecies
      integer, save :: i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl
      logical,          save :: iAmRootProc
      integer,          save :: procID
      integer,          save :: rootProc       
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null() 

                          ! longitude index where max_const_asc (see below) 
                          ! occurs, by layer, for output to ".asc" file
      integer, pointer, save :: maxi_asc(:) => null()
                          ! latitude  index where max_const_asc (see below) 
                          ! occurs, by layer, for output to ".asc" file
      integer, pointer, save :: maxj_asc(:) => null()
                          ! longitude index where max_const_asc (see below) 
                          ! occurs, by layer, for output to ".asc" file
      integer, pointer, save :: mini_asc(:) => null()
                          ! latitude  index where min_const_asc (see below) 
                          ! occurs, by layer, for output to ".asc" file
      integer, pointer, save :: minj_asc(:) => null()
                          ! one selectable specie's concentration for output
                          ! known at zone centers (mixing ratio)
      real*8 , pointer, save :: const_asc     (:,:,:) => null()
                          ! mass of selected species, by layer, for output (kg)
      real*8 , pointer, save :: mass_const_asc    (:) => null()
                          ! total mass of each species for output (kg)
      real*8 , pointer, save :: mass_const_all_asc(:) => null()
                          ! maximum of const_asc, by layer, for output (mixing ratio)
      real*8 , pointer, save :: max_const_asc     (:) => null()
                          ! minimum of const_asc, by layer, for output (mixing ratio)
      real*8 , pointer, save :: min_const_asc     (:) => null()
                          ! total production of each species for output (kg)
      real*8 , pointer, save :: prod_const_all_asc(:,:,:) => null()
                          ! total loss of each species for output (kg)
      real*8 , pointer, save :: loss_const_all_asc(:,:,:) => null()
                          ! net loss of each species for output (kg)
      real*8 , pointer, save :: net_const_all_asc (:,:,:) => null()
      real*8 , pointer, save :: londeg (:) => null()
      real*8 , pointer, save :: mw (:) => null()
      integer, pointer, save :: fixed_const_map(:) => null()
      integer         , save :: num_fixed_const
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOutputASCII
!
! !INTERFACE:
!
      subroutine initializeOutputASCII(gmiGrid, gmiDomain, Diagnostics, &
     &                     Chemistry, speciesConcentration)
!
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid    
      type(t_gmiDomain  ), intent(in) :: gmiDomain  
      type(t_Chemistry  ), intent(in) :: Chemistry  
      type(t_Diagnostics), intent(in) :: Diagnostics
      type(t_SpeciesConcentration), intent(in) :: speciesConcentration

! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) 'initializeOutputASCII called by ', procID

      call Get_pr_ascii(Diagnostics, pr_ascii)

      call Get_i1(gmiGrid, i1)
      call Get_i2(gmiGrid, i2)
      call Get_j2(gmiGrid, j2)
      call Get_k1(gmiGrid, k1)
      call Get_k2(gmiGrid, k2)
      call Get_ju1(gmiGrid, ju1)
      call Get_i1_gl(gmiGrid, i1_gl)
      call Get_i2_gl(gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl )
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_numSpecies(gmiGrid, numSpecies)

      call Get_rootProc         (gmiDomain, rootProc  )
      call Get_numDomains       (gmiDomain, numDomains )
      call Get_iAmRootProc      (gmiDomain, iAmRootProc )
      call Get_communicatorWorld(gmiDomain, commuWorld)

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      call Get_pr_ascii1(Diagnostics, pr_ascii1)
      call Get_pr_ascii2(Diagnostics, pr_ascii2)
      call Get_pr_ascii3(Diagnostics, pr_ascii3)
      call Get_pr_ascii4(Diagnostics, pr_ascii4)
      call Get_pr_ascii5(Diagnostics, pr_ascii5)
      call Get_ascii_out_i(Diagnostics, ascii_out_i)
      call Get_ascii_out_n(Diagnostics, ascii_out_n)
      call Get_pr_ascii_step_interval(Diagnostics, pr_ascii_step_interval)

      allocate(const_labels(numSpecies))

      allocate(mw(numSpecies))
      call Get_mw(Chemistry, mw)

      if (pr_ascii4) then
         call Get_num_fixed_const(speciesConcentration, num_fixed_const)

         allocate(fixed_const_map(num_fixed_const))
         call Get_fixed_const_map (speciesConcentration, fixed_const_map, num_fixed_const)
      end if

      if (iAmRootProc) then
         allocate(londeg(i1_gl:i2_gl))
         call Get_londeg(gmiDomain, londeg)

         call Get_problem_name(Diagnostics, problem_name)
         call makeOutfileName (fname, '.asc', problem_name)

         call AsciiOpenWrite (asclun, fname)

         call Get_const_labels(Chemistry, const_labels)

         Write (asclun, 800) ascii_out_n, const_labels(ascii_out_n)
         Write (asclun, *)

 800     format (1x, 'Species chosen to look at:  ', i3, ' => ', a)
      end if

      call allocateOutputASCII()

      return

      end subroutine initializeOutputASCII
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutputASCII
!
! !INTERFACE
!
      subroutine controlOutputASCII (SpeciesConcentration, Chemistry, &
     &                  gmiClock, metFields)
!
! !USES:
      use GmiTimeControl_mod  , only : GmiSplitDateTime
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiClock            ), intent(in) :: gmiClock            
      type(t_metFields           ), intent(in) :: metFields           
      type(t_Chemistry           ), intent(in) :: Chemistry           
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! Controls the ascii file output.
!
! !LOCAL VARIABLES:
      logical :: time_for_asc
      integer :: day, month, num_time_steps, nymd
      real*8  :: gmi_sec
      integer :: idumday, idumyear
      integer, save :: month_save = -999
      real*8, allocatable :: prod(:,:,:)
      real*8, allocatable :: loss(:,:,:)
      real*8, allocatable :: const_ascGlob(:,:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
      character (len=MAX_LENGTH_SPECIES_NAME) :: const_labels(numSpecies)
      integer :: num_chem, num_molefrac
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputASCII called by ', procID

      time_for_asc = .false.

      call Get_curGmiDate(gmiClock, nymd)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)

      if (pr_ascii_step_interval < 0) then  ! monthly output
         call GmiSplitDateTime (nymd, idumyear, month, idumday)

         if (month_save == -999) month_save = month

         if (month /= month_save) time_for_asc = .true.

          month_save = month
      else if (Mod (num_time_steps, pr_ascii_step_interval) == 0) then
         time_for_asc = .true.
      end if

      if (time_for_asc) then
         call Get_concentration(SpeciesConcentration, concentration)

         if (pr_ascii5) then
            Allocate (loss(2, NUM_OPERATORS, 1:numSpecies))
            call Get_loss(SpeciesConcentration, loss)

            Allocate (prod(2, NUM_OPERATORS, 1:numSpecies))
            call Get_prod(SpeciesConcentration, prod)
         end if

         call prepASCIIoutput (metFields, prod, loss, concentration)

         if (pr_ascii5) then
            deallocate (loss)
            deallocate (prod)
         end if

         if (iAmRootProc) then
            if (pr_ascii2) allocate(const_ascGlob(i1_gl:i2_gl,ju1_gl:j2_gl,k1:k2))
         end if

         call Sub_To_Glob_Asc (const_ascGlob)

         if (iAmRootProc) then
            call Get_num_chem    (Chemistry, num_chem    )
            call Get_num_molefrac(Chemistry, num_molefrac)

            call Get_const_labels(Chemistry, const_labels)

            call writeASCIIoutput (const_ascGlob, const_labels, num_chem, &
     &                num_molefrac, gmi_sec)

            if (pr_ascii2) deallocate(const_ascGlob)
         end if

      end if

      return

      end subroutine controlOutputASCII
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeOutputASCII
! 
! !INTERFACE:
!
      subroutine finalizeOutputASCII ()
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce qj outputs.
!
! !LOCAL VARIABLES:
      integer :: ierr
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) ' finalizeOutputASCII called by ', procID

      deallocate (maxi_asc)
      deallocate (maxj_asc)
      deallocate (mini_asc)
      deallocate (minj_asc)

      deallocate (mass_const_all_asc)

      deallocate (mass_const_asc)
      deallocate (max_const_asc )
      deallocate (min_const_asc )

      if (pr_ascii2) then
         deallocate (const_asc)
      end if

      if (pr_ascii5) then
         deallocate (loss_const_all_asc)
         deallocate (prod_const_all_asc)
         deallocate (net_const_all_asc )
      end if

      if (iAmRootProc) call releaseFileUnitNumber(asclun, ierr)

      return

      end subroutine finalizeOutputASCII
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateOutputASCII
!
! !INTERFACE:
!
      subroutine allocateOutputASCII ()

!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateOutputASCII called by ', procID

      Allocate (maxi_asc(k1:k2))
      Allocate (maxj_asc(k1:k2))
      Allocate (mini_asc(k1:k2))
      Allocate (minj_asc(k1:k2))
      maxi_asc = 0
      maxj_asc = 0
      mini_asc = 0
      minj_asc = 0

      Allocate (mass_const_all_asc(1:numSpecies))
      mass_const_all_asc = 0.0d0

      Allocate (mass_const_asc(k1:k2))
      Allocate (max_const_asc (k1:k2))
      Allocate (min_const_asc (k1:k2))
      mass_const_asc = 0.0d0
      max_const_asc  = 0.0d0
      min_const_asc  = 0.0d0

      if (pr_ascii2) then
        Allocate (const_asc(i1:i2, ju1:j2, k1:k2))
        const_asc = 0.0d0
      end if

      if (pr_ascii5) then
         Allocate (loss_const_all_asc(2, NUM_OPERATORS, 1:numSpecies))
         Allocate (prod_const_all_asc(2, NUM_OPERATORS, 1:numSpecies))
         Allocate (net_const_all_asc (2, NUM_OPERATORS, 1:numSpecies))
         loss_const_all_asc = 0.0d0
         prod_const_all_asc = 0.0d0
         net_const_all_asc  = 0.0d0
      end if

      return

      end subroutine allocateOutputASCII
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: prepASCIIoutput
!
! !INTERFACE:
!
      subroutine prepASCIIoutput (metFields, prod, loss, concentration)

      implicit none
!
! !INPUT PARAMETERS:
      real*8                 , intent(in) :: loss(i1:,ju1:,1:)
      real*8                 , intent(in) :: prod(i1:,ju1:,1:)
      type (t_metFields     ), intent(in) :: metFields
      type (t_GmiArrayBundle), intent(in) :: concentration(*)
!
! !DESCRIPTION:
!  This routine calls Do_Prep_Ascii, passing mass as an argument rather
!  than as a pointer array in a common block.  This facilitates
!  vectorization, etc.
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik, ic
      real*8  :: mw_fac
      real*8, allocatable :: mass(:,:,:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'prepASCIIoutput called by ', procID

      allocate(mass(i1:i2,ju1:j2,k1:k2))
      call Get_mass(metFields, mass)

!     ==============
      if (pr_ascii1) then
!     ==============

        mw_fac = mw(ascii_out_n) / MWTAIR

        do ik = k1, k2

          mass_const_asc(ik) = 0.0d0

          do ij = ju1, j2
            do il = i1, i2

              mass_const_asc(ik) =  &
     &          mass_const_asc(ik) +  &
     &          ((mass(il,ij,ik) * mw_fac) *  &
     &           concentration(ascii_out_n)%pArray3D(il,ij,ik))

            end do
          end do

        end do

      end if


!     ==============
      if (pr_ascii2) then
!     ==============

        const_asc(:,:,:) = concentration(ascii_out_n)%pArray3D(:,:,:)

      end if


!     ==============
      if (pr_ascii3) then
!     ==============

!       -----------------------------
!       Update const min's and max's.
!       -----------------------------

        do ik = k1, k2

          min_const_asc(ik) = concentration(ascii_out_n)%pArray3D(i1,ju1,ik)
          mini_asc(ik)      = i1
          minj_asc(ik)      = ju1

          max_const_asc(ik) = concentration(ascii_out_n)%pArray3D(i1,ju1,ik)
          maxi_asc(ik)      = i1
          maxj_asc(ik)      = ju1

          do ij = ju1, j2
            do il = i1, i2

              if (concentration(ascii_out_n)%pArray3D(il,ij,ik) <  &
     &            min_const_asc(ik)) then

                min_const_asc(ik) = concentration(ascii_out_n)%pArray3D(il,ij,ik)
                mini_asc(ik)      = il
                minj_asc(ik)      = ij

              else if (concentration(ascii_out_n)%pArray3D(il,ij,ik) >  &
     &                 max_const_asc(ik)) then

                max_const_asc(ik) = concentration(ascii_out_n)%pArray3D(il,ij,ik)
                maxi_asc(ik)      = il
                maxj_asc(ik)      = ij

              end if

            end do
          end do

        end do

      end if


!     ==============
      if (pr_ascii4) then
!     ==============

        do ic = 1, numSpecies

           mw_fac = mw(ic) / MWTAIR

          mass_const_all_asc(ic) = 0.0d0

          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2

                mass_const_all_asc(ic) =  &
     &            mass_const_all_asc(ic) +  &
     &            ((mass(il,ij,ik) * mw_fac) * concentration(ic)%pArray3D(il,ij,ik))

              end do
            end do
          end do

        end do

      end if


!     ==============
      if (pr_ascii5) then
!     ==============

        do ic = 1, numSpecies

           mw_fac = mw(ic) / MWTAIR

          prod_const_all_asc(:,:,ic) =  &
     &      (mw_fac * prod(:,:,ic))

          loss_const_all_asc(:,:,ic) =  &
     &      (mw_fac * loss(:,:,ic))

        end do

        net_const_all_asc(:,:,:) =  &
     &    prod_const_all_asc(:,:,:) - loss_const_all_asc(:,:,:)

      end if

      deallocate(mass)

      return

      end subroutine prepASCIIoutput
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeASCIIoutput
!
! !INTERFACE
!
      subroutine writeASCIIoutput (const_ascGlob, const_labels, &
     &                num_chem, num_molefrac, gmi_sec)
!
! !USES:

      implicit none
!
! !INPUT PARAMETERS:
      character (len=MAX_LENGTH_SPECIES_NAME) , intent(in) :: const_labels(:)
      integer            , intent(in) :: num_chem, num_molefrac
      real*8             , intent(in) :: const_ascGlob(i1_gl:,ju1_gl:,k1:)
      real*8             , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! Outputs specified ASCII data.
!
! !DEFINED PARAMETERS:
      logical, parameter :: MORE_PRECISION = .false.
!
! !LOCAL VARIABLES:
      character (len=38)  :: fmt
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      integer :: ij, ik, ic, io
      integer :: nfixed
      real*8  :: dmtotl
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'writeASCIIoutput called by ', procID

!     ==============
      if (pr_ascii1) then
!     ==============

        dmtotl = Sum (mass_const_asc(:))

        Write (asclun, 802)
        Write (asclun, 803) gmi_sec / SECPDY
        Write (asclun, 802)
        Write (asclun, *)
        Write (asclun, 805) ascii_out_n, dmtotl
        Write (asclun, *)
        Write (asclun, 810) (mass_const_asc(ik), ik=k1,k2)
        Write (asclun, *)

 802    format (1x, '=============================')
 803    format (1x, 'Mass of consts at time (days):  ', f11.5)
 805    format (1x, 'Total sum of const ', i3, ' in atm:  ', es18.8)
 810    format (1x, 'Sum by level:'/, (3es20.8))

      end if


!     ==============
      if (pr_ascii2) then
!     ==============

        Write (asclun, 830) gmi_sec / SECPDY
        Write (asclun, *)
        Write (asclun, 832) ascii_out_i, londeg(ascii_out_i)
        Write (asclun, *)

        do ik = k1, k2
          Write (asclun, 835) ik
          Write (asclun, 840)  &
     &      (const_ascGlob(ascii_out_i,ij,ik), ij=ju1_gl,j2_gl)
        end do

        Write (asclun, *)

 830    format (1x, 'Final mixing ratios at time (days):  ', f11.5)
 832    format (1x, 'Longitude chosen to look at:  ', i3, ' => ', f7.1)
 835    format (1x, 'k = ', i3)
 840    format (3es20.8)

      end if


!     ==============
      if (pr_ascii3) then
!     ==============

        Write (asclun, 860)
        Write (asclun, 865)

        do ik = k1, k2
          Write (asclun, 870)  &
     &      ik, maxi_asc(ik), maxj_asc(ik), max_const_asc(ik),  &
     &          mini_asc(ik), minj_asc(ik), min_const_asc(ik)
        end do

        Write (asclun, *)

 860    format ('   k ', '  imax ', ' jmax ', t28, 'amax', 6x,  &
     &                   '  imin ', ' jmin ', t60, 'amin')
 865    format ('  --- ', ' ---- ', ' ---- ', t28, '----', 7x,  &
     &                    ' ---- ', ' ---- ', t60, '----')
 870    format (i5, 2(2i6, es20.8))

      end if


!     ==============
      if (pr_ascii4) then
!     ==============

        Write (asclun, 880) gmi_sec / SECPDY
        Write (asclun, *)

        Write (asclun, 885)
        Write (asclun, 890)

        if (MORE_PRECISION) then
          fmt(1:38) = '(2x, i3, 2x, a16, 2x, i4, 6x, es23.16)'
        else
          fmt(1:38) = '(2x, i3, 2x, a16, 2x, i4, 6x, es20.8)'
        end if

        do ic = 1, numSpecies

          nfixed = 0
          if (Any (fixed_const_map(:) == ic)) nfixed = 1

          if ((ic <= num_molefrac) .or. (ic > num_chem)) then

!           -----------------------------------------------
!           Do not write out total mass of O2, N2, ad, etc.
!           -----------------------------------------------

            Write (asclun, fmt)  &
     &        ic, const_labels(ic), nfixed, mass_const_all_asc(ic)

          end if

        end do

        Write (asclun, *)

 880    format(1x, 'Total mass of each const in atm. at time (days):  ', f11.5)
 885    format(2x, ' ic', 2x, '      name      ', 2x, ' fixed? ',  &
     &     2x, '      total mass')
 890    format(2x, '---', 2x, '----------------', 2x, '--------',  &
     &     2x, '-----------------------')

      end if

!     ==============
      if (pr_ascii5) then
!     ==============

        Write (asclun, *) 'TOTAL cumulative prod/loss:'

        Write (asclun, 900)
        Write (asclun, 905)

        ! We do not want to account for ForcedBC (that is part of Chemistry)
        ! twice. We need to compute the cumulative prod/loss for the first
        ! NUM_OPERATORS-1 operators only.

        Write (asclun, 910)  &
     &    (ic, const_labels(ic),  &
     &     Sum (prod_const_all_asc(CUMULATIVE,1:NUM_OPERATORS-1,ic)),  &
     &     Sum (loss_const_all_asc(CUMULATIVE,1:NUM_OPERATORS-1,ic)),  &
     &     Sum (net_const_all_asc (CUMULATIVE,1:NUM_OPERATORS-1,ic)),  &
     &     ic=1,numSpecies)

        Write (asclun, *)

        do io = 1, NUM_OPERATORS

          Write (asclun, *)  &
     &      'Cumulative prod/loss for operator:  ', OPERATOR_NAME(io)

          Write (asclun, 900)
          Write (asclun, 905)

          Write (asclun, 910)  &
     &      (ic, const_labels(ic),  &
     &       prod_const_all_asc(CUMULATIVE,io,ic),  &
     &       loss_const_all_asc(CUMULATIVE,io,ic),  &
     &       net_const_all_asc (CUMULATIVE,io,ic),  &
     &       ic=1,numSpecies)

          Write (asclun, *)

        end do

      end if

 900    format  &
     &    (2x, ' ic',  &
     &     2x, '      name      ',  &
     &     2x, '   production   ',  &
     &     2x, '      loss      ',  &
     &     2x, '      net       ')
 905    format  &
     &    (2x, '---',  &
     &     2x, '----------------',  &
     &     2x, '----------------',  &
     &     2x, '----------------',  &
     &     2x, '----------------')
 910    format  &
     &    (2x, i3, 2x, a16, 2x, es16.8, 2x, es16.8, 2x, es16.8)

      call GmiFlush (asclun)

      return

      end subroutine writeASCIIoutput
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Sub_To_Glob_Asc
!
! !INTERFACE:
!
      subroutine Sub_To_Glob_Asc (const_ascGlob )
!
      implicit none
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: const_ascGlob(i1_gl:,ju1_gl:,k1:)
!
! !DESCRIPTION:
! Collects ASCII data from the worker PEs and sends it to the root PE for output.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_CONST_ASC = 5001
      integer, parameter :: SG_MAX_ASC(3)  = (/ 5002, 50003, 5004 /)
      integer, parameter :: SG_MIN_ASC(3)  = (/ 5005, 50006, 5007 /)
!
! !LOCAL VARIABLES:
      integer :: s2g_size, ivert
      real*8  :: local_redarray(2*NUM_OPERATORS*numSpecies)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) Write (6,*) 'Asc_Sub_To_Glob called by ', procID

      ivert = k2 - k1 + 1

      if (pr_ascii1) then

!       ===================
        call Gmi_Sum_Reduce (k1, k2, mass_const_asc, commuWorld, rootProc, procID)
!       ===================

      end if

      if (pr_ascii2) then

!       ========================
        call subDomain2Global  &
!       ========================
     &    (const_ascGlob, const_asc, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, &
     &     j2, k1, k2, rootProc, procID, &
     &     map1_u, numDomains, SG_CONST_ASC, commuWorld)

      end if

      if (pr_ascii3) then

!       ======================
        call Gmi_Maxpij_Reduce  &
!       ======================
     &    (ivert, SG_MAX_ASC, maxi_asc, maxj_asc, max_const_asc, &
     &     commuWorld, rootProc, procID, numDomains, ivert)

!       ======================
        call Gmi_Minpij_Reduce  &
!       ======================
     &    (ivert, SG_MIN_ASC, mini_asc, minj_asc, min_const_asc, &
     &     commuWorld, rootProc, procID, numDomains, ivert)

      end if

      if (pr_ascii4) then

!       ===================
        call Gmi_Sum_Reduce(1, numSpecies, mass_const_all_asc, commuWorld, rootProc, procID)
!       ===================

      end if


      if (pr_ascii5) then

        s2g_size = 2 * NUM_OPERATORS * numSpecies

        local_redarray(:) =  Reshape (prod_const_all_asc(:,:,:), (/ s2g_size /))

!       ===================
        call Gmi_Sum_Reduce  &
!       ===================
     &    (1, s2g_size, local_redarray(:), commuWorld, rootProc, procID)

        prod_const_all_asc(:,:,:) =  Reshape (local_redarray(:),  &
     &             (/ 2, NUM_OPERATORS, numSpecies /))

        local_redarray(:) =   Reshape (loss_const_all_asc(:,:,:), (/ s2g_size /))

!       ===================
        call Gmi_Sum_Reduce  &
!       ===================
     &    (1, s2g_size, local_redarray(:), commuWorld, rootProc, procID)

        loss_const_all_asc(:,:,:) =  Reshape (local_redarray(:),  &
     &             (/ 2, NUM_OPERATORS, numSpecies /))


        local_redarray(:) = Reshape (net_const_all_asc(:,:,:), (/ s2g_size /))

!       ===================
        call Gmi_Sum_Reduce  &
!       ===================
     &    (1, s2g_size, local_redarray(:), commuWorld, rootProc, procID)

        net_const_all_asc(:,:,:) =   Reshape (local_redarray(:),  &
     &             (/ 2, NUM_OPERATORS, numSpecies /))

      end if

      return

      end subroutine Sub_To_Glob_Asc
!EOC
!------------------------------------------------------------------------------
      end module GmiControlASCII_mod
