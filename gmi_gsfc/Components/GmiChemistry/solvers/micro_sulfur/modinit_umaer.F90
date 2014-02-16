      module modinit_umaer

      contains

      subroutine umaerinit(nso4_i,nmomso4_i,nnon_i,nmomnon_i,ng_i,  &
     &                     sigmod_i,radmerg_i,  &
     &                     radgsrc_i,siggsrc_i,  &
     &                     radgso4_i,siggso4_i,fracso4_i,  &
     &                     radgnon_i,siggnon_i,fracnon_i,  &
     &                     rhonon_i,xmolnon_i,nrh_flg_i,  &
     &                     crh1_i,crh2_i,crh3_i,crh4_i)
!
!---->initialize umaerosol: input for number of modes and momemts
!                           info about background aerosol size distributions
!
!     nso4             number of sulfate aerosol modes
!     nmomso4          number of sulfate aerosol moments
!     nnon             number of non-sulfate aerosol modes
!     nmomnon          number of non-sulfate aerosol moments
!
!     ng               number of lognormal distr. in background aerosol
!
!---->assumed width and position of each mode
!     order: from smallest to largest mode
!
!     sigmod(nso4)     assumed geometric standard deviation  [1]
!     radmerg(nso4)    merging radius                        [m]
!
!---->sulfate aerosol in near SO2 source regions:
!
!     treated as direct particle emissions
!     radgsrc(ng)      geometric radius                      [m]
!     siggsrc(ng)      geometric standard deviation          [1]
!
!---->background sulfate aerosol:
!
!     radgso4(ng)      geometric radius                      [m]
!     siggso4(ng)      geometric standard deviation          [1]
!     fracso4(ng)      relative fraction                     [1]
!
!---->background non-sulfate aerosol:
!
!     radgnon(ng,nnon) geometric radius                      [m]
!     siggnon(ng,nnon) geometric standard deviation          [1]
!     fracnon(ng,nnon) relative fraction                     [1]
!     rhonon(nnon)     aerosol density                       [kg/m3]
!     xmolnon(nnon)    molecular weight                      [kg/mol]
!     nrh_flg(nnon)    flag for humidity growth (1->sea salt)
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      integer ::  &
     &     nso4_i,nmomso4_i,nnon_i,nmomnon_i,ng_i
      real(kdef) ::  &
     &     sigmod_i(nso4_i),radmerg_i(nso4_i),  &
     &     radgsrc_i(ng_i),siggsrc_i(ng_i),  &
     &     radgso4_i(ng_i),siggso4_i(ng_i),fracso4_i(ng_i),  &
     &     radgnon_i(ng_i,nnon_i),siggnon_i(ng_i,nnon_i),  &
     &     fracnon_i(ng_i,nnon_i),rhonon_i(nnon_i),xmolnon_i(nnon_i)
      integer :: nrh_flg_i(nnon_i)
      real(kdef) :: crh1_i(nnon_i),crh2_i(nnon_i),  &
     &              crh3_i(nnon_i),crh4_i(nnon_i)
!
!---->define local arrays and parameter
!
      parameter (r10=10.d0,r100=100.d0,  &
     &           r1h=0.5d0,r3h=1.5d0,  &
     &           r4td=4.d0/3.d0,r1td=1.d0/3.d0)
      parameter (r5=5.d0,r25=25.d0)
      parameter (r6=6.d0,r36=36.d0)
      dimension radvolm(ng_i)
!$$$c
!$$$c---->density of dry sulfate aerosol (pure sulfuric acid) [kg/m3]
!$$$c
!$$$      parameter (rhoso4=1861.d0)
!
!---->so4 density of dry ammonium sulfate aerosol [kg/m3]
!
      parameter (rhoso4=1313.d0)
!
!---->allocate and set data for module modpar_umaer
!
      nso4   =nso4_i
      nmomso4=nmomso4_i
      nnon   =nnon_i
      nmomnon=nmomnon_i
      ng     =ng_i

      naer=nmomso4*nso4
      ntot=nnon+naer+1

      allocate(radgsrc(1:ng),siggsrc(1:ng),  &
     &         radgso4(1:ng),siggso4(1:ng),fracso4(1:ng))
      radgsrc(1:ng)=radgsrc_i(1:ng)
      siggsrc(1:ng)=siggsrc_i(1:ng)
      radgso4(1:ng)=radgso4_i(1:ng)
      siggso4(1:ng)=siggso4_i(1:ng)
      fracso4(1:ng)=fracso4_i(1:ng)

      allocate (radgnon(1:ng,1:nnon),siggnon(1:ng,1:nnon),  &
     &          fracnon(1:ng,1:nnon),  &
     &          rhonon(1:nnon),xmolnon(1:nnon))
      radgnon(1:ng,1:nnon)=radgnon_i(1:ng,1:nnon)
      siggnon(1:ng,1:nnon)=siggnon_i(1:ng,1:nnon)
      fracnon(1:ng,1:nnon)=fracnon_i(1:ng,1:nnon)
      rhonon(1:nnon)=rhonon_i(1:nnon)
      xmolnon(1:nnon)=xmolnon_i(1:nnon)

      allocate (nrh_flg(1:nnon))
      nrh_flg(1:nnon)=nrh_flg_i(1:nnon)

      allocate (crh1(1:nnon),crh2(1:nnon),  &
     &          crh3(1:nnon),crh4(1:nnon))
      crh1(1:nnon)=crh1_i(1:nnon)
      crh2(1:nnon)=crh2_i(1:nnon)
      crh3(1:nnon)=crh3_i(1:nnon)
      crh4(1:nnon)=crh4_i(1:nnon)

!
!---->allocate memory for module modpar_umaer
!
      allocate(xscalemin(ntot))
      allocate(fracmerg(nmerg+1),radgmerg(nmerg+1,nso4),  &
     &         radmerg(nso4),pmsfac(nso4))
      allocate(fracaccu(0:naccu+1,nso4),radgaccu(0:naccu+1,nso4,2))
      allocate(fraccnon(0:naccu+1),radgaccnon(0:naccu+1,nnon))
      allocate(pmssrc(ng))
      allocate(pmsnon(nnon),radvnon(nnon))
      allocate(etaso4(nso4),  &
     &         etanon(nnon),alphanon(nnon))
!
!---->allocate memory for module modbox_umaer
!
      allocate(pmsmerg(nvec,nso4+1))
      allocate(condnon(nvec,nnon),xnonum(nvec,nnon),  &
     &         rwetnon(nvec,nnon),xnonkelvin(nvec,nnon),  &
     &         diffnon(nvec,nnon))

      allocate(so4sig(nvec,nso4),xnonsig(nvec,nnon))
      allocate(fracsrc(nvec,nso4),fracsrcnon(nvec,nnon))
      allocate(so4alpha(nvec,nso4))
!
!---->set maximum allowed error
!
      xscalemin(1:nnon)=1000.d0
      xscalemin(nnon+1:nnon+naer)=10.d0
!$$$      xscalemin(nnon+1:nnon+naer:nmomso4)=1000.d0
      do inon=nnon+1,nnon+naer,nmomso4
       xscalemin(inon)=1000.d0
      enddo
      xscalemin(nnon+(nso4-1)*nmomso4+2:ntot-1)=200.d0
      xscalemin(ntot)=10000.d0
!
!---->set accommodation coefficient
!
#ifdef LOWACCOEF
      etaso4(1:nso4)=0.04d0
      etanon(1:nnon)=0.04d0
#else
      etaso4(1:nso4)=1.d0
      etanon(1:nnon)=1.d0
#endif
!
!---->calculate particle mass of directly emitted particles
!     in near SO2 source regions
!
      radvolm(:)=radgsrc(:)*exp(r3h*log(siggsrc(:))**2)
      pmssrc(:)=r4td*pi*radvolm(:)**3*rhoso4*avno/so4mol
!
!---->calculate volume mean radii for background aerosol in [m]
!     and particle mass in [kg]
!
      radvolm(:)=radgso4(:)*exp(r3h*log(siggso4(:))**2)
      radvso4=sum(fracso4(:)*radvolm(:)**3)
      radvso4=radvso4**r1td

      do inon=1,nnon
       radvolm(:)=radgnon(:,inon)*exp(r3h*log(siggnon(:,inon))**2)
       radvnon(inon)=sum(fracnon(:,inon)*radvolm(:)**3)
       pmsnon(inon)=r4td*pi*rhonon(inon)*radvnon(inon)
       radvnon(inon)=radvnon(inon)**r1td
       xnonsig(1,inon)=sum(fracnon(:,inon)*log(siggnon(:,inon))**2)  &
     &                /sum(fracnon(:,inon))
       xnonsig(:,inon)=exp(sqrt(xnonsig(1,inon)))
      enddo
!
!---->correction factor for condensation for non-sulfate aerosol
!     for using the volume mean radius instead of the mean radius
!
      do inon=1,nnon
       alphanon(inon)=exp(-sum(fracnon(:,inon)*log(siggnon(:,inon))**2))
      enddo
!
!---->so4sig  - assumed geometric standard deviation [1]
!     radmerg - radius of merged particles [m]
!     pmsfac  - factor for mass of a merged particle [#molec m3/kg]
!
      do iso4=1,nso4
       iinv=nso4-iso4+1
       so4sig(:,iinv)=sigmod_i(iso4)
       radmerg(iinv)=radmerg_i(iso4)
       pmsfac(iinv)=r4td*pi*radmerg(iinv)**3/so4mol*avno
      enddo
!
!---->correction factor for condensation for sulfate aerosol
!     for using the volume mean radius instead of the mean radius
!
      do iso4=1,nso4
       so4alpha(:,iso4)=exp(-log(so4sig(:,iso4))**2)
      enddo
!
!---->calculate integrals over size distribution for merging
!
      call mergeinit
!
!---->calculate integrals over size distribution for accumulation mode
!
      call accuminit

      end subroutine umaerinit
!-----------------------------------------------------------------------
      subroutine mergeinit
!
!---->calculate integrals over size distribution for merging
!     determine geometric mean radii for each bin
!     for a given fractional change due to merging
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->define local arrays and parameter
!
      parameter (r5=5.d0,r10=10.d0,r1h=0.5d0,r3h=1.5d0)

      dimension func(3),radg(3)
!
!---->search radg for given value of integral over size distribution
!
      fracint=r1/real(nmerg,kdef)
      do iso4=2,nso4
       sigma=so4sig(1,iso4)
       do imerg=1,nmerg
        fracmerg(imerg)=fracint*(real(imerg,kdef)-r1h)

        radg(1)=radmerg(iso4)/(r5*sigma)
        radg(3)=radmerg(iso4)* r5*sigma
        func(1)=fracmerg_f(radg(1),sigma,radmerg(iso4))-fracmerg(imerg)
        func(3)=fracmerg_f(radg(3),sigma,radmerg(iso4))-fracmerg(imerg)

        do iter=1,20
         radg(2)=r1h*(radg(1)+radg(3))
         func(2)=fracmerg_f(radg(2),sigma,radmerg(iso4))-fracmerg(imerg)

         flag=r1h-sign(r1h,func(2)*func(3))
         flaginv=r1-flag

         radg(1)=radg(1)*flaginv+radg(2)*flag
         radg(3)=radg(2)*flaginv+radg(3)*flag
         func(1)=func(1)*flaginv+func(2)*flag
         func(3)=func(2)*flaginv+func(3)*flag
        enddo
!
!---->do linear interpolation
!
        radg(2)=radg(1)-(radg(3)-radg(1))  &
     &         *func(1)/(func(3)-func(1))
        radg(2)=max(radg(1),min(radg(3),radg(2)))

        radgmerg(imerg,iso4)=radg(2)
       enddo
      enddo
!
!---->set upper boundary
!
      fracmerg(nmerg+1)=r1
      radgmerg(nmerg+1,2:nso4)=radgmerg(nmerg  ,2:nso4)  &
     &                        *radgmerg(nmerg  ,2:nso4)  &
     &                        /radgmerg(nmerg-1,2:nso4)

      end subroutine mergeinit
!-----------------------------------------------------------------------
      function fracmerg_f(radg,sigma,radmrg)
!
!---->calculate integral over part of lognormal distribution
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->set local parameter
!
      parameter (r10=10.d0,r1h=0.5d0)
      parameter (dfrac=0.05d0)

      if (sigma<=r1) then
       if (radg>radmrg) then
        fracmerg_f=r1
       else
        fracmerg_f=r0
       endif
      else
       xlnsg=log(sigma)
       xln2sg=xlnsg*xlnsg
       xln2sginv=-r1h/xln2sg
       rscale=max(radg,radmrg)
       drad=dfrac*rscale*(r1-r1/sigma)
       sqr2pi=drad/(xlnsg*sqrt(r2*pi))

       fracmerg_f=r0
       radloc=radmrg+r1h*drad
       do while (radloc<r10*rscale)
        fracmerg_f=fracmerg_f  &
     &            +exp(xln2sginv*(log(radloc/radg))**2)/radloc
        radloc=radloc+drad
       enddo
       fracmerg_f=fracmerg_f*sqr2pi
      endif

      end function fracmerg_f
!-----------------------------------------------------------------------
      subroutine accuminit
!
!---->calculate integrals over size distribution for accumulation mode
!     determine geometric mean radii for each bin
!     for a given fraction in the accumulation mode
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->define local arrays and parameter
!
      parameter (r5=5.d0,r10=10.d0,r1h=0.5d0,r3h=1.5d0)

      dimension func(3),radg(3)
      dimension fracmax(nso4)
!
!---->define accumulation mode as radius interval [m]
!     calculate limit radius (maximum overlap)    [m]
!
      radaccu(1:2)=(/0.05d-6,1.25d-6/)
      racculim=sqrt(radaccu(1)*radaccu(2))
!
!---->search radg for given value of integral over size distribution
!
      fracint=r1/real(naccu,kdef)
      do iso4=1,nso4
       sigma=so4sig(1,iso4)
       fracmax(iso4)=fracaccu_f(racculim,sigma)
       do iaccu=1,naccu
        fracaccu(iaccu,iso4)=fracint*(real(iaccu,kdef)-r1h)
        if (fracaccu(iaccu,iso4)>=fracmax(iso4)) then
         radgaccu(iaccu,iso4,1:2)=racculim
        else
         do iloc=1,2
          if (iloc==1) then
           radg(1)=radaccu(iloc)/(r5*sigma)
           radg(3)=racculim
          else
           radg(1)=racculim
           radg(3)=radaccu(iloc)* r5*sigma
          endif

          func(1)=fracaccu_f(radg(1),sigma)-fracaccu(iaccu,iso4)
          func(3)=fracaccu_f(radg(3),sigma)-fracaccu(iaccu,iso4)
          if (func(1)*func(3)>r0) then
           print*, 'ACCUMINIT: no zero found'
           stop
          endif

          do iter=1,20
           radg(2)=r1h*(radg(1)+radg(3))
           func(2)=fracaccu_f(radg(2),sigma)-fracaccu(iaccu,iso4)

           flag=r1h-sign(r1h,func(2)*func(3))
           flaginv=r1-flag

           radg(1)=radg(1)*flaginv+radg(2)*flag
           radg(3)=radg(2)*flaginv+radg(3)*flag
           func(1)=func(1)*flaginv+func(2)*flag
           func(3)=func(2)*flaginv+func(3)*flag
          enddo
!
!---->do linear interpolation
!
          radg(2)=radg(1)-(radg(3)-radg(1))  &
     &           *func(1)/(func(3)-func(1))
          radg(2)=max(radg(1),min(radg(3),radg(2)))

          radgaccu(iaccu,iso4,iloc)=radg(2)
         enddo
        endif
       enddo
      enddo
!
!---->set lower/upper boundary
!
      fracaccu(0,1:nso4)=r0
      radgaccu(0,1:nso4,1:2)=r0
      fracaccu(naccu+1,1:nso4)=fracmax(1:nso4)
      radgaccu(naccu+1,1:nso4,1:2)=racculim
#ifndef NONON
!
!---->for non-sulfate aerosol:
!     all particles larger than radaccu(1) are in accumulation mode
!
!---->search radg for given value of integral over size distribution
!
      fracint=r1/real(naccu,kdef)
      do inon=1,nnon
       sigma=xnonsig(1,inon)
       do iaccu=1,naccu
        fraccnon(iaccu)=fracint*(real(iaccu,kdef)-r1h)

        radg(1)=radaccu(1)/(r5*sigma)
        radg(3)=radaccu(1)* r5*sigma
        func(1)=fracmerg_f(radg(1),sigma,radaccu(1))-fraccnon(iaccu)
        func(3)=fracmerg_f(radg(3),sigma,radaccu(1))-fraccnon(iaccu)

        do iter=1,20
         radg(2)=r1h*(radg(1)+radg(3))
         func(2)=fracmerg_f(radg(2),sigma,radaccu(1))-fraccnon(iaccu)

         flag=r1h-sign(r1h,func(2)*func(3))
         flaginv=r1-flag

         radg(1)=radg(1)*flaginv+radg(2)*flag
         radg(3)=radg(2)*flaginv+radg(3)*flag
         func(1)=func(1)*flaginv+func(2)*flag
         func(3)=func(2)*flaginv+func(3)*flag
        enddo
!
!---->do linear interpolation
!
        radg(2)=radg(1)-(radg(3)-radg(1))  &
     &         *func(1)/(func(3)-func(1))
        radg(2)=max(radg(1),min(radg(3),radg(2)))

        radgaccnon(iaccu,inon)=radg(2)
       enddo
      enddo
!
!---->set lower/upper boundary
!
      fraccnon(0)=r0
      radgaccnon(0,1:nnon)=r0
      fraccnon(naccu+1)=r1
      radgaccnon(naccu+1,1:nnon)=radgaccnon(naccu  ,1:nnon)  &
     &                          *radgaccnon(naccu  ,1:nnon)  &
     &                          /radgaccnon(naccu-1,1:nnon)
#endif
      end subroutine accuminit
!-----------------------------------------------------------------------
      function fracaccu_f(radg,sigma)
!
!---->calculate integral over part of lognormal distribution
!     between radaccu(1) and radaccu(2)
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->set local parameter
!
      parameter (r10=10.d0,r1h=0.5d0)
      parameter (dfrac=0.05d0)

      if (sigma<=r1) then
       if (radg>radaccu(1) .and. radg<radaccu(2)) then
        fracaccu_f=r1
       else
        fracaccu_f=r0
       endif
      else
       xlnsg=log(sigma)
       xln2sg=xlnsg*xlnsg
       xln2sginv=-r1h/xln2sg
       drad=dfrac*radg*(r1-r1/sigma)
       sqr2pi=drad/(xlnsg*sqrt(r2*pi))

       fracaccu_f=r0
       radloc=radaccu(1)+r1h*drad
       do while (radloc<radaccu(2))
        fracaccu_f=fracaccu_f  &
     &            +exp(xln2sginv*(log(radloc/radg))**2)/radloc
        radloc=radloc+drad
       enddo
       fracaccu_f=fracaccu_f*sqr2pi
      endif

      end function fracaccu_f
!-----------------------------------------------------------------------
      end module modinit_umaer
