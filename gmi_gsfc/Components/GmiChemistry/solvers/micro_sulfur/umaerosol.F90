      subroutine umaerosol(dt,nbox,temp,press,relhum,cldfrac,  &
     &                     so4gas,srcgas,srcmas,srcpar,  &
     &                     so4aer,so4non,aernon)
!
!---->do condensation, nucleation, and coagulation for sulfate aerosol,
!     calculate condensation of so4 on non-sulfate aerosol
!
!     dt                time step                             [s]
!     nbox              number of parcels
!     temp(nbox)        temperature                           [K]
!     press(nbox)       air pressure                          [Pa]
!     relhum(nbox)      relative humidity                     [0-1]
!     so4gas(nbox)      so4 gas concentration                 [#molec/m3]
!     srcgas(nbox)      so4 source rate (gas phase)           [#so4 molec/m3/s]
!     srcmas(nbox)      so4 source rate (aqueous phase)       [#so4 molec/m3/s]
!     srcpar            h2so4 source rate (particulate)       [#so4 molec/m3/s]
!
!     so4aer(nmomso4,nso4,nbox) sulfate aerosol               [#/m3]
!     so4non(nnon,nbox) so4 attached to non-sulfate aerosol   [#so4 molec/m3]
!     aernon(nmomnon,nnon,nbox) non-sulfate aerosol           [kg/m3]
!
!---->number of modes and moments set in umaerinit and stored in modpar_umaer
!
!     nso4              number of sulfate aerosol modes
!     nmomso4           number of sulfate aerosol moments
!     nnon              number of non-sulfate aerosol modes
!     nmomnon           number of non-sulfate aerosol moments
!
!     moments:   1 - mass
!                2 - number
!                3 - surface
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
#ifdef DIAG
      use modiag_umaer
#endif
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nbox),press(nbox),relhum(nbox),cldfrac(nbox),  &
     &          so4gas(nbox),srcgas(nbox),srcmas(nbox),srcpar(nbox),  &
     &          so4aer(nmomso4,nso4,nbox)
      dimension so4non(nnon,nbox),  &
     &          aernon(nmomnon,nnon,nbox)
!
!---->local arrays and parameter
!     xaer              temporary array for sulfate aerosol
!                       (all modes and moments)
!
      dimension xaer(nvec,nnon+nmomso4*nso4+1)

      dimension iterbox(nvec),locbox(nvec)
      dimension dtbox(nvec),time(nvec)
      dimension tempbox(nvec),pressbox(nvec),relhumbox(nvec),  &
     &          srcgasbox(nvec),srcmasbox(nvec),srcparbox(nvec)
#ifdef DIAG
      dimension iterstat(15)
      parameter (r1h=0.5d0)
#endif
      parameter (r4=4.d0,r6=6.d0,r10=10.d0,r3q=0.75d0,r1td=1.d0/3.d0)
      parameter (a1=1.333d0,a2=0.71d0,b1=0.572127d0)
      parameter (c1=1.249d0,c2=0.42d0,c3=-0.87d0)
      parameter (rhmin=1.d-2,rhmax=0.99d0)

#ifdef DIAG
      iterstat(:)=0
#endif
      maxstat=0
!
!---->loop over all parcels, do calculation for chunk of nvec parcels
!
      nv=nvec
      ibox=0
      idone=0
      time(1:nvec)=dt
      do while (idone==0)
!
!---->fill in new data as long as not all boxes are considered
!
       if (ibox<nbox) then
!
!---->fill local arrays
!
        do ivec=1,nvec
         if (time(ivec)>=dt) then
!
!---->update parcel counter
!
          ibox=ibox+1
          time(ivec)=r0
          iterbox(ivec)=0
!
!---->remember parcel index and preset time step
!
          locbox(ivec)=ibox
          dtbox(ivec)=-r1
!
!---->store sulfate stuff
!
          xaer(ivec,1:nnon)=so4non(1:nnon,ibox)
          do iso4=1,nso4
           do imom=1,nmomso4
            iaer=ntot-(iso4-1)*nmomso4-nmomso4+imom-1
            xaer(ivec,iaer)=so4aer(imom,iso4,ibox)
           enddo
          enddo
          xaer(ivec,ntot)=so4gas(ibox)
#ifdef DIAG
          call diagfill(ivec,ibox)
#endif
!
!---->store non-so4 aerosol number concentrations       [#part/m3]
!
          xnonum(ivec,1:nnon)=aernon(1,1:nnon,ibox)/pmsnon(1:nnon)
!
!---->store meteorological data and sources
!
          tempbox(ivec)=temp(ibox)
          pressbox(ivec)=press(ibox)
          relhumbox(ivec)=max(rhmin,relhum(ibox))
          srcgasbox(ivec)=srcgas(ibox)
          srcmasbox(ivec)=srcmas(ibox)
          srcparbox(ivec)=srcpar(ibox)
!
!---->the following quantities don't change during time integration
!     and are stored in common block /parcel/
!---->vapor pressures:
!     h2osat        h2o saturation pressure             [#molec/m3]
!     h2ogas        h2o gas concentration               [#molec/m3]
!     so4sat        h2so4 saturation pressure           [#molec/m3]
!
          h2osat(ivec)=h2osat_f(temp(ibox))
          h2ogas(ivec)=relhumbox(ivec)*h2osat(ivec)
          so4sat(ivec)=so4sat_f(temp(ibox))
!
!---->aerosol composition and Kelvin effect:
!     so4mfrac      so4 mass fraction in aerosol        [1]
!     so4frac       so4 mole fraction in aerosol        [1]
!     so4dens       aerosol density                     [kg/m3]
!     so4kelvin     factor for Kelvin effect            [m]
!
          so4mfrac(ivec)=so4mfrac_f(temp(ibox),h2osat(ivec),  &
     &                                         h2ogas(ivec))
          so4frac(ivec)=so4mfrac(ivec)*h2omol  &
     &             /(so4mfrac(ivec)*h2omol+(r1-so4mfrac(ivec))*so4mol)
          so4dens(ivec)=so4dens_f(temp(ibox),so4mfrac(ivec))

          so4kelvin(ivec)=  &
     &             r2*so4sig_f(temp(ibox),so4frac(ivec),so4mfrac(ivec))  &
     &                *(so4frac(ivec)*so4mol+(r1-so4frac(ivec))*h2omol)  &
     &                /(rgas*temp(ibox)*so4dens(ivec))
!
!---->rwetnon       humid volume mean radius            [m]
!
!     if nrh_flg=1 assume it is sea salt,
!     humidity growth for sea salt according to Gerber (1985)
!
          rhnew=max(rhmin,min(rhmax,relhum(ibox)))
          factor=r3q*so4min/(pi*so4mfrac(ivec)*so4dens(ivec))

          do inon=1,nnon
           if (nrh_flg(inon)==1) then
!           rwetnon(ivec,inon)=(radvnon(inon)**3
!    $                   +(crh1*radvnon(inon)**crh2)
!    $                  /((crh3*radvnon(inon)**crh4)
!    $                   -log10(rhnew))
!    $                   +factor*xaer(ivec,inon)
!    $                 /max(r1,xnonum(ivec,inon)))**r1td
            rwetnon(ivec,inon)=(radvnon(inon)**3  &
     &                   +(crh1(inon)*radvnon(inon)**crh2(inon))  &
     &                  /((crh3(inon)*radvnon(inon)**crh4(inon))  &
     &                   -log10(rhnew))  &
     &                   +factor*xaer(ivec,inon)  &
     &                 /max(r1,xnonum(ivec,inon)))**r1td
           else
            rwetnon(ivec,inon)=(radvnon(inon)**3  &
     &                +factor*xaer(ivec,inon)  &
     &              /max(r1,xnonum(ivec,inon)))**r1td
           endif
          enddo

!$$$          rwetnon(ivec,1:nnon-2)=(radvnon(1:nnon-2)**3
!$$$     $                  +factor*xaer(ivec,1:nnon-2)
!$$$     $                /max(r1,xnonum(ivec,1:nnon-2)))**r1td
!$$$          rwetnon(ivec,nnon-1:nnon)=(radvnon(nnon-1:nnon)**3
!$$$     $                   +(crh1*radvnon(nnon-1:nnon)**crh2)
!$$$     $                  /((crh3*radvnon(nnon-1:nnon)**crh4)
!$$$     $                   -log(rhnew))
!$$$     $                   +factor*xaer(ivec,nnon-1:nnon)
!$$$     $                 /max(r1,xnonum(ivec,nnon-1:nnon)))**r1td

          rwetnon(ivec,1:nnon)=min(rwetnon(ivec,1:nnon),  &
     &                              r10*radvnon(1:nnon))
!
!---->Kelvin effect for non-sulfate aerosol
!     assuming that non-sulfate aerosol is completely coated with so4
!
          xnonkelvin(ivec,1:nnon)=exp(so4kelvin(ivec)  &
     &                                 /rwetnon(ivec,1:nnon))
!
!---->thermodynamic quantities:
!     pathair       mean free path of air               [m]
!     pathso4       mean free path of h2so4 in air      [m]
!     diffso4       diffusion coeff. for h2so4 in air   [m2/s]
!     viscos        dynamic viscosity of air            [Ns/m2]
!
          pathair(ivec)=pathair_f(temp(ibox),press(ibox))
          pathso4(ivec)=b1*pathair(ivec)
          diffso4(ivec)=r4*pi*diffso4_f(temp(ibox),press(ibox))
          viscos (ivec)=viscos_f(temp(ibox))
#ifdef KULMALA
!
!---->quantities for nucleation:
!     so4crit       h2so4 conc with nuc' rate of 1/cm3  [#molec/m3]
!
          so4crit(ivec)=so4crit_f(temp(ibox),relhum(ibox))
#else
          so4crit(ivec)=r0
#endif
!
!---->quantities for condensation:
!     so4sol        h2so4 sat pressure over solution    [#molec/m3]
!     condnon       condensation coef's for non-so4 aer [m3/#molec]
!                   incl. correction for collision geometry
!                   following Fuchs and Sutugin (1971)
!
          so4sol(ivec)=so4sol_f(temp(ibox),so4sat(ivec),so4frac(ivec))
          so4crit(ivec)=max(r2*so4sol(ivec),so4crit(ivec))

          do inon=1,nnon
           xknu=pathso4(ivec)/rwetnon(ivec,inon)
           condnon(ivec,inon)=diffso4(ivec)*alphanon(inon)  &
     &                       /(r1+((a1*xknu+a2)/(r1+xknu)  &
     &                       +a1*(r1-etanon(inon))/etanon(inon))*xknu)  &
     &               *rwetnon(ivec,inon)
          enddo
!
!---->quantities for coagulation
!
          do inon=1,nnon
           xknu=pathair(ivec)/rwetnon(ivec,inon)
           diffnon(ivec,inon)=r4*boltz*temp(ivec)  &
     &                  /(r6*pi*viscos(ivec)*rwetnon(ivec,inon))  &
     &                  *(r1+xknu*(c1+c2*exp(c3/xknu)))
          enddo
!
!---->quantities for merging:
!     pmsmerg       limit particle mass for merging     [#molec]
!
          do iso4=1,nso4
           pmsmerg(ivec,iso4)=so4dens(ivec)*so4mfrac(ivec)*pmsfac(iso4)
          enddo
          pmsmerg(ivec,nso4+1)=r2
         endif
!
!---->exit if all boxes considered
!
         if (ibox>=nbox) then
          exit
         endif
        enddo
       endif
!
!----->shrink local array is all boxes are considered
!
       if (ibox>=nbox) then
        nvnew=nv
        do iv=1,nv
         if (iv>nvnew) exit
         if (time(iv)>=dt) then
          do jv=nvnew,iv+1,-1
           if (time(jv)<dt) then
!
!---->shift parcel data
!
            iterbox(iv)=iterbox(jv)
            locbox(iv)=locbox(jv)
            dtbox(iv)=dtbox(jv)
            time(iv)=time(jv)
            time(jv)=dt

            xaer(iv,1:ntot)=xaer(jv,1:ntot)
            xnonum(iv,1:nnon)=xnonum(jv,1:nnon)

            tempbox(iv)=tempbox(jv)
            pressbox(iv)=pressbox(jv)
            relhumbox(iv)=relhumbox(jv)
            srcgasbox(iv)=srcgasbox(jv)
            srcmasbox(iv)=srcmasbox(jv)
            srcparbox(iv)=srcparbox(jv)

            h2osat(iv)=h2osat(jv)
            h2ogas(iv)=h2ogas(jv)
            so4sat(iv)=so4sat(jv)

            so4mfrac(iv)=so4mfrac(jv)
            so4frac(iv)=so4frac(jv)
            so4dens(iv)=so4dens(jv)

            so4kelvin(iv)=so4kelvin(jv)

            rwetnon(iv,1:nnon)=rwetnon(jv,1:nnon)
            xnonkelvin(iv,1:nnon)=xnonkelvin(jv,1:nnon)
            diffnon(iv,1:nnon)=diffnon(jv,1:nnon)

            pathair(iv)=pathair(jv)
            pathso4(iv)=pathso4(jv)
            diffso4(iv)=diffso4(jv)
            viscos (iv)=viscos (jv)

            so4crit(iv)=so4crit(jv)
            so4sol(iv)=so4sol(jv)

            condnon(iv,1:nnon)=condnon(jv,1:nnon)

            pmsmerg(iv,1:nso4+1)=pmsmerg(jv,1:nso4+1)
#ifdef DIAG
            call diagshift(iv,jv)
#endif
            exit
           endif
          enddo
          nvnew=jv-1
         endif
        enddo
        nv=nvnew
       endif
!
!---->do one iteration of for time integration using Runge-Kutta
!
       call aertime(dt,nv,iterbox,dtbox,time,  &
     &              tempbox,pressbox,relhumbox,  &
     &              srcgasbox,srcmasbox,srcparbox,  &
     &              xaer)
!
!---->if iter>itermax set parcel's time to dt
!
       if (maxval(iterbox(1:nv))>itermax) then
        print*,'warning! iter>itermax!'
        do iv=1,nv
         if (iterbox(iv)>itermax) exit
        enddo
        write(*,200) locbox(iv),dt,time(iv),dtbox(iv)
        time(iv)=dt
       endif
!
!---->check overall success
!
       if (ibox>=nbox) then
        if (minval(time(1:nv))>=dt) idone=1
       endif
!
!---->check success for nv boxes
!
       do iv=1,nv
        if (time(iv)>=dt) then
#ifdef DIAG
!
!---->store number of iterations for statistics
!
         istat=(r2*iterbox(iv))**r1td
         istat=max(1,min(15,istat))
         iterstat(istat)=iterstat(istat)+1
         maxstat=max(maxstat,iterbox(iv))
#endif
!
!---->restore new concentrations
!
         ib=locbox(iv)
         so4non(1:nnon,ib)=xaer(iv,1:nnon)
         do iso4=1,nso4
          do imom=1,nmomso4
           iaer=ntot-(iso4-1)*nmomso4-nmomso4+imom-1
           so4aer(imom,iso4,ib)=xaer(iv,iaer)
          enddo
         enddo
         so4gas(ib)=xaer(iv,ntot)
#ifdef DIAG
         call diagrefill(iv,ib)
#endif
        endif
       enddo
!
!---->end over parcel loop
!
      enddo
#ifdef DIAG
!
!---->write statistics for iterations
!
      write(*,100) dt
      write(*,110) maxstat
      do i=1,15
       istat=r1h*(i+1)**3-mod(i,2)
       write(*,120) iterstat(i),istat
      enddo
#endif
      return
 100  format(/,'timestep: ',f10.0)
 110  format('maximum number of iteration: ',i7)
 120  format(i8,' less equal ',i5,' iterations')
 200  format(' parcel: ',i6,/,  &
     &       ' requested total time: ',f10.0,'s',/,  &
     &       ' reached total time:   ',f10.0,'s',/,  &
     &       ' last time step:       ',1pe10.2,'s',/)
      end
!-----------------------------------------------------------------------
      subroutine aertime(dt,nv,iterbox,dtbox,time,  &
     &                   temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                   xaer)
!
!---->time integration for nv parcels
!     Runge Kutta 4th order with adaptive time stepping:
!     W.H. Press et al., Numerical Recipes, Cambridge University Press (1988),
!     chapter 15.
!
!     xaer     array of nv parcels with initial values:
!              xaer(1:nv,1:nnon)            so4 on non-sulfate aerosol
!              xaer(1:nv,nnon+1:nnon+naer)  so4 aerosol (nso4*nmomso4)
!              xaer(1:nv,nnon+naer+1)       so4 gas
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension iterbox(nv)
      dimension dtbox(nv),time(nv),  &
     &          temp(nv),press(nv),relhum(nv),  &
     &          srcgas(nv),srcmas(nv),srcpar(nv)
      dimension xaer(nvec,ntot)
!
!---->define local arrays and parameter
!
      dimension xtmp(nvec,ntot),fxaer(nv,ntot),  &
     &          xsav(nvec,ntot),fxsav(nv,ntot),  &
     &          dthalf(nv),dt6inv(nv),errmax(nv),  &
     &          so4tot(nv),so4cor(nv)

      parameter (r4=4.d0,r1000=1000.d0,r1h=0.5d0,r1q=0.25d0)
      parameter (r6inv=r1/6.d0)
      parameter (predu=-0.25d0,pgrow=-0.2d0,fcor=r1/15.d0,  &
     &           safety=0.90d0,errcon=(r4/safety)**(-5))
      parameter (relacc=5.0d-2)
      parameter (so4max=1.d12)
!
!---->declare interfaces
!
      interface
      subroutine derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                      xaer,fxaer,dt,iflg)
      use precision
      use modpar_umaer
      use modbox_umaer
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nv),press(nv),relhum(nv),  &
     &          srcgas(nv),srcmas(nv),srcpar(nv)
      dimension xaer(nvec,ntot),fxaer(nv,ntot),dt(nv)
      target xaer
      end subroutine derivative

      subroutine rk4step(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                   xaer,fxaer,dtbox,xout)
      use precision
      use modpar_umaer
      use modbox_umaer
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nv),press(nv),relhum(nv),  &
     &          srcgas(nv),srcmas(nv),srcpar(nv)
      dimension xaer(nvec,ntot),fxaer(nv,ntot),dtbox(nv),  &
     &          xout(nvec,ntot)
      target xaer,xtmp
      end subroutine rk4step

      subroutine merge(ivec,xaer)
      use precision
      use modpar_umaer
      use modbox_umaer
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension xaer(nvec,ntot)
      target xaer
      end subroutine merge
      end interface
!
!---->calculate initial time derivates (fxaer)
!
#ifdef DIAG
      call diagpoint(nv,1,1)
#endif
      dt6inv(1:nv)=max(r0,r6inv*dtbox(1:nv))
      call derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                xaer,fxaer,dt6inv,1)
!
!---->estimate initial time step if time=0
!
      do iv=1,nv
       if (time(iv)==r0 .and. iterbox(iv)==0) then
        dtbox(iv)=min(dtmax,dt)
        do itot=1,ntot-1
         if (xaer(iv,itot)*fxaer(iv,itot)<r0)  &
     &    dtbox(iv)=min(dtbox(iv),max(dtmin,  &
     &             -safety*xaer(iv,itot)/fxaer(iv,itot)))
        enddo
        if (xaer(iv,ntot)*fxaer(iv,ntot)<r0)  &
     &    dtbox(iv)=min(dtbox(iv),max(dtmin,  &
     &             -safety*min(so4max,xaer(iv,ntot))/fxaer(iv,ntot)))
       endif
      enddo
!
!---->save old values and derivatives
!
      xsav(1:nv,1:ntot)=xaer(1:nv,1:ntot)
      fxsav(1:nv,1:ntot)=fxaer(1:nv,1:ntot)
!
!---->take first half Runge-Kutta step
!
      dthalf(1:nv)=r1h*dtbox(1:nv)
      call rk4step(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &             xaer,fxaer,dthalf,xtmp)
!
!---->ensure positive definiteness and mass conservation
!
      do iv=1,nv
       so4tot(iv)=xtmp(iv,ntot)+sum(xtmp(iv,1:nnon))
       so4cor(iv)=max(r0,xtmp(iv,ntot))+sum(max(r0,xtmp(iv,1:nnon)))
      enddo
      do iso4=1,nso4
       inum=nnon+iso4*nmomso4
       imas=inum-1
       do iv=1,nv
        so4tot(iv)=so4tot(iv)+xtmp(iv,imas)
        so4cor(iv)=so4cor(iv)+max(r0,xtmp(iv,imas))
        if (xtmp(iv,inum)<=r0) xtmp(iv,inum)=xtmp(iv,imas)  &
     &                                      /pmsmerg(iv,iso4+1)
        if (xtmp(iv,imas)<=r0) xtmp(iv,inum)=r0
       enddo
      enddo
      do iv=1,nv
       xtmp(iv,1:ntot)=max(r0,xtmp(iv,1:ntot))
      enddo
      do iv=1,nv
       if (so4cor(iv)>so4tot(iv)) then
        factor=so4tot(iv)/so4cor(iv)
        xtmp(iv,1:ntot)=factor*xtmp(iv,1:ntot)
       endif
      enddo
!
!---->take second half Runge-Kutta step
!
#ifdef DIAG
      call diagstore(nv,1)
      call diagpoint(nv,0,1)
#endif
      call derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                xtmp,fxaer,dt6inv,0)
      call rk4step(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &             xtmp,fxaer,dthalf,xaer)
!
!---->take one large Runge-Kutta step
!
#ifdef DIAG
      call diagrestore(nv,1)
      call diagpoint(nv,2,1)
#endif
      call rk4step(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &             xsav,fxsav,dtbox,xtmp)
!
!---->calculate error, new time step and concentrations,
!     advance time counter
!
      errmax(1:nv)=r0
      do itot=1,ntot-1
       do iv=1,nv
        dxaer=xaer(iv,itot)-xtmp(iv,itot)
        xscale=abs(xsav(iv,itot))+dtbox(iv)*abs(fxsav(iv,itot))
        xscale=max(xscalemin(itot),xscale)
        errmax(iv)=max(errmax(iv),abs(dxaer/xscale))
       enddo
      enddo
      do iv=1,nv
       dxaer=xaer(iv,ntot)-xtmp(iv,ntot)
       xscale=abs(xsav(iv,ntot))  &
     &       +max(dtbox(iv)*abs(fxsav(iv,ntot)),r1000*srcgas(iv))
       xscale=max(xscalemin(ntot),xscale)
       errmax(iv)=max(errmax(iv),abs(dxaer/xscale))
      enddo
      do iv=1,nv
       iterbox(iv)=iterbox(iv)+1
       errmax(iv)=errmax(iv)/relacc
       if (errmax(iv)>r1) then
        dtbox(iv)=max(r1q*dtbox(iv),  &
     &                safety*dtbox(iv)*(errmax(iv)**predu))
        xaer(iv,1:ntot)=xsav(iv,1:ntot)
        xtmp(iv,1:ntot)=r0
       else
        dtold=dtbox(iv)
        time(iv)=time(iv)+dtbox(iv)
        dtend=max(r0,dt-time(iv)+epsilon)
        if (errmax(iv)>errcon) then
         dtbox(iv)=min(dtend,safety*dtbox(iv)*(errmax(iv)**pgrow))
        else
         dtbox(iv)=min(dtend,r4*dtbox(iv))
        endif
!
!---->use xtmp to achieve 5th order accuracy
!     and ensure positive definiteness and mass conservation
!
        so4tot(iv)=xaer(iv,ntot)+sum(xaer(iv,1:nnon))
        do itot=1,nso4
         inum=nnon+itot*nmomso4
         imas=inum-1
         so4tot(iv)=so4tot(iv)+xaer(iv,imas)
        enddo
        xtmp(iv,1:ntot)=xaer(iv,1:ntot)  &
     &          +fcor*(xaer(iv,1:ntot)-xtmp(iv,1:ntot))

        if (minval((xaer(iv,1:ntot)-xsav(iv,1:ntot))  &
     &            *(xtmp(iv,1:ntot)-xsav(iv,1:ntot)))>=r0) then
         xaer(iv,1:ntot)=xtmp(iv,1:ntot)
#ifdef DIAG
         call diagupdate(time(iv),dtold,iv,1)
        else
         call diagupdate(time(iv),dtold,iv,0)
#endif
        endif

        so4cor(iv)=max(r0,xaer(iv,ntot))+sum(max(r0,xaer(iv,1:nnon)))
        do itot=1,nso4
         inum=nnon+itot*nmomso4
         imas=inum-1
         so4cor(iv)=so4cor(iv)+max(r0,xaer(iv,imas))
         if (xaer(iv,inum)<=r0) xaer(iv,inum)=xaer(iv,imas)  &
     &                                       /pmsmerg(iv,itot+1)
         if (xaer(iv,imas)<=r0) xaer(iv,inum)=r0
        enddo
        xaer(iv,1:ntot)=max(r0,xaer(iv,1:ntot))
        if (so4cor(iv)>so4tot(iv)) then
         factor=so4tot(iv)/so4cor(iv)
         xaer(iv,1:ntot)=factor*xaer(iv,1:ntot)
        endif
!
!---->shift particles to larger bins if overlap
!
        call merge(iv,xaer)
       endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine rk4step(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                   xaer,fxaer,dtbox,xout)
!
!---->do explicit time step using 4th order Runge-Kutta
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nv),press(nv),relhum(nv),  &
     &          srcgas(nv),srcmas(nv),srcpar(nv)
      dimension xaer(nvec,ntot),fxaer(nv,ntot),dtbox(nv),  &
     &          xout(nvec,ntot)
      target xaer,xout
!
!---->define local arrays and parameter
!
      dimension dthalf(nv),dt6inv(nv),  &
     &          fxtmp(nv,ntot),fxmid(nv,ntot)

      parameter (r1h=0.5d0,r6inv=r1/6.d0)
!
!---->declare interfaces
!
      interface
      subroutine derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                      xaer,fxaer,dt,iflg)
      use precision
      use modpar_umaer
      use modbox_umaer
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nv),press(nv),relhum(nv),  &
     &          srcgas(nv),srcmas(nv),srcpar(nv)
      dimension xaer(nvec,ntot),fxaer(nv,ntot),dt(nv)
      target xaer
      end subroutine derivative
      end interface
!
!---->start Runge-Kutta
!
      dthalf(1:nv)=r1h*dtbox(1:nv)
      dt6inv(1:nv)=r6inv*dtbox(1:nv)
!
!---->first step
!
      do itot=1,ntot
       do iv=1,nv
        xout(iv,itot)=max(r0,xaer(iv,itot)+dthalf(iv)*fxaer(iv,itot))
       enddo
      enddo
      do iso4=1,nso4
       inum=nnon+iso4*nmomso4
       imas=inum-1
       do iv=1,nv
        xout(iv,inum)=min(xout(iv,inum),xout(iv,imas)  &
     &                              /pmsmerg(iv,iso4+1))
        if (xout(iv,inum)<=r0) xout(iv,inum)=xout(iv,imas)  &
     &                              /pmsmerg(iv,iso4+1)
       enddo
      enddo
!
!---->second step
!
#ifdef DIAG
      call diagpoint(nv,0,3)
#endif
      call derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                xout,fxtmp,dt6inv,0)
      do itot=1,ntot
       do iv=1,nv
        xout(iv,itot)=max(r0,xaer(iv,itot)+dthalf(iv)*fxtmp(iv,itot))
       enddo
      enddo
      do iso4=1,nso4
       inum=nnon+iso4*nmomso4
       imas=inum-1
       do iv=1,nv
        xout(iv,inum)=min(xout(iv,inum),xout(iv,imas)  &
     &                              /pmsmerg(iv,iso4+1))
        if (xout(iv,inum)<=r0) xout(iv,inum)=xout(iv,imas)  &
     &                              /pmsmerg(iv,iso4+1)
       enddo
      enddo
!
!---->third step
!
#ifdef DIAG
      call diagpoint(nv,0,2)
#endif
      call derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                xout,fxmid,dt6inv,0)
      do itot=1,ntot
       do iv=1,nv
        xout(iv,itot)=max(r0,xaer(iv,itot)+dtbox(iv)*fxmid(iv,itot))
        fxmid(iv,itot)=fxmid(iv,itot)+fxtmp(iv,itot)
       enddo
      enddo
      do iso4=1,nso4
       inum=nnon+iso4*nmomso4
       imas=inum-1
       do iv=1,nv
        xout(iv,inum)=min(xout(iv,inum),xout(iv,imas)  &
     &                              /pmsmerg(iv,iso4+1))
        if (xout(iv,inum)<=r0) xout(iv,inum)=xout(iv,imas)  &
     &                              /pmsmerg(iv,iso4+1)
       enddo
      enddo
!
!---->fourth step
!
#ifdef DIAG
      call diagadd(nv,2,3)
      call diagpoint(nv,0,3)
#endif
      call derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                xout,fxtmp,dt6inv,0)
      do itot=1,ntot
       do iv=1,nv
        xout(iv,itot)=xaer(iv,itot)  &
     &               +dt6inv(iv)*(fxaer(iv,itot)+r2*fxmid(iv,itot)  &
     &                                             +fxtmp(iv,itot))
       enddo
      enddo
#ifdef DIAG
      call diagchange(nv,dt6inv)
#endif
      return
      end
!-----------------------------------------------------------------------
      subroutine derivative(nv,temp,press,relhum,srcgas,srcmas,srcpar,  &
     &                      xaer,fxaer,dt,iflg)
!
!---->calculate time derivates (fxaer)
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension xaer(nvec,ntot),fxaer(nv,ntot),dt(nv)
      dimension temp(nv),press(nv),relhum(nv),  &
     &          srcgas(nv),srcmas(nv),srcpar(nv)
      target xaer
!
!---->local arrays and parameter
!
      dimension so4radv(nv,nso4),  &
     &          so4knu(nv,nso4),so4vel(nv,nso4)

      parameter (r8=8.d0,r3q=0.75d0,r1td=1.d0/3.d0)
!
!---->define pointer
!
      dimension so4non(:,:),so4aer(:,:),so4gas(:)
      pointer so4non,so4aer,so4gas
!
!---->set pointer
!
      so4non=>xaer(1:nvec,1:nnon)
      so4aer=>xaer(1:nvec,nnon+1:nnon+naer)
      so4gas=>xaer(1:nvec,ntot)
!
!---->initialize derivatives
!
      fxaer(1:nv,1:ntot)=r0
!
!---->calculate wet aerosol size (Tang and Munkelwitz (1994)),
!               Knudsen number,
!               (square of) characteristic velocity
!
      do iso4=1,nso4
       inum=iso4*nmomso4
       imas=inum-1
       do iv=1,nv
        wetmas=max(r2,so4aer(iv,imas)  &
     &        /max(epsilon/pmsmerg(iv,iso4+1),so4aer(iv,inum)))  &
     &        *so4min/so4mfrac(iv)
        so4vel (iv,iso4)=r8*boltz*temp(iv)/(pi*wetmas)
        so4radv(iv,iso4)=(r3q*wetmas/(pi*so4dens(iv)))**r1td
        so4knu (iv,iso4)=pathair(iv)/so4radv(iv,iso4)
       enddo
      enddo
!
!---->calculate fraction in accumulation mode
!
      if (iflg==1) call accumode(nv,so4aer,so4radv,so4non)
!
!---->nucleation
!
      call nucleation(nv,temp,relhum,so4gas,fxaer)
!
!---->condensation on so4 aerosol
!
      call condensso4(nv,so4gas,so4aer,so4radv,so4knu,dt,fxaer)
!
!---->condensation on non-so4 aerosol
!
#ifndef NONON
      call condensnon(nv,so4gas,so4non,dt,fxaer)
#endif
!
!---->coagulation
!
      call coagulation(nv,temp,so4gas,so4aer,so4non,  &
     &                 so4radv,so4knu,so4vel,fxaer)
!
!---->external (chemical) sources
!
      call sources(nv,so4aer,so4radv,srcgas,srcmas,srcpar,fxaer)
#ifdef GRAV
!
!---->gravitational settling
!
      call settling(nv,so4aer,so4radv,fxaer)
#endif
      return
      end
!-----------------------------------------------------------------------
      subroutine nucleation(nv,temp,relhum,so4gas,fxaer)
!
!---->calculate binary nucleation rate
!     if -DKULMALA       according to Kulmala et al. (1998)
!     elif -DKULNEW      according to Vehkamaeki et al. (2002)
!     elif -DFITZGERALD  according to Fitzgerald et al. (1998)
!     else               according to classsical theory Zhao, Turco (1995)
!
!     nucleation rate   rnucnum [#particles/m3]
!       and mass rate   rnucmas [#molecules so4/m3]
!

!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nv),relhum(nv)
      dimension so4gas(nvec),fxaer(nv,ntot)
!
!---->define limit values for nucleation rate and mass of nucleus
!
      parameter (rnucmax=1.d12)
      parameter (rmasmin=2.d0,rmasmax=300.d0)
#ifdef KULMALA
!
!---->calculate number and mass of nucleated particles per time
!     following Kulmala et al. (1998)
!     JGR (103) pp 8301-8307
!

!
!---->define local parameter
!
      parameter (r3h=1.5d0,r3q=0.75d0,r4=4.d0)
      parameter (tmin=233.d0,tmax=299.d0,  &
     &           rhmin=0.1d0,rhmax=1.d0)
      parameter (h2osatmin=5.79d21,h2osatmax=8.056d23)
      parameter (ramax=1.d2)
      parameter (so4satmin=4.26d13,so4satmax=5.008d17)
      parameter (dmfracnew=0.005d0)
      parameter (smfracint=0.0005d0)
      parameter (fracmin=0.005d0)
!
!---->indices for so4gas and smallest particles
!
      inum=nso4*nmomso4
      ixnum=nnon+inum
      ixmas=ixnum-1
      igso4=ntot
!
!---->start loop over parcels
!
      do iv=1,nv
!
!---->limit temperature and relative humidity
!     calculate h2ogas that matches relhnew and tempnew
!
       tempnew=max(tmin,min(tmax,temp(iv)))
       relhnew=max(rhmin,min(rhmax,relhum(iv)))
       h2osatnew=max(h2osatmin,min(h2osatmax,h2osat(iv)))
       h2onew=relhnew*h2osatnew
!
!---->calculate relative acidity
!     and so4gas that matches relacd and tempnew
!
       relacd=min(ramax,so4gas(iv)/so4sat(iv))
       so4new=relacd*max(so4satmin,min(so4satmax,so4sat(iv)))
!
!---->calculate h2so4 mole/mass fraction in critical nucleus
!
       so4nuc=so4fracnuc_f(tempnew,relhnew,relacd,h2onew,so4new)
       if (so4nuc>fracmin) then
        so4mnuc=so4nuc*so4mol  &
     &        /(so4nuc*so4mol+(r1-so4nuc)*h2omol)
        so4mnuc=min(r1,max(r0,so4mnuc))
!
!---->calculate nucleation rate: rnucnum [#particles/m3]
!
        rnucnum=rnucnum_f(tempnew,relhnew,so4gas(iv),so4nuc,so4crit(iv))
!.deb   rnucnum=min(rnucmax,rnucnum)
!
!---->calculate density of nucleus
!     Jaecker-Voirol, 1988
!
        rhonuc=so4dens_f(tempnew,so4mnuc)
!
!---->calculate partial molecular volumes
!
        dmfrac=dmfracnew*sign(r1,so4mnuc-smfracint)
        drhonuc=(rhonuc-so4dens_f(tempnew,so4mnuc-dmfrac))/dmfrac
        h2ovol=h2omol/avno*(r1/rhonuc+drhonuc*so4mnuc/rhonuc**2)
!       so4vol=so4mol/avno*(r1/rhonuc-drhonuc*(r1-so4mnuc)/rhonuc**2)
!
!---->calculate change in chemical potential
!     using activity coefficients for water and sulfuric acid
!     Taleb et al. (1996)
!
        h2osol=h2osol_f(tempnew,h2osatnew,so4nuc)
        h2odmu=-boltz*tempnew*log(h2onew/h2osol)
!
!---->calculate surface tension of nucleus
!     Sabinina and Terbugow, 1935
!
        signuc=so4sig_f(tempnew,so4nuc,so4mnuc)
!
!---->calculate size and mass of critical nucleus
!     rnucmasnuc [#so4 molecules / nucleus]
!
        radnuc=-r2*signuc*h2ovol/h2odmu
        radnuc=max(radmin,radnuc)
        rnucmas=pi/r3q*radnuc**3*rhonuc*so4mnuc/so4mol*avno
        rnucmas=rnucnum*min(rmasmax,max(rmasmin,rnucmas))
!
!---->add to derivatives
!
        fxaer(iv,ixmas)=fxaer(iv,ixmas)+rnucmas ! mass   (small part)
        fxaer(iv,ixnum)=fxaer(iv,ixnum)+rnucnum ! number (small part)
        fxaer(iv,igso4)=fxaer(iv,igso4)-rnucmas ! so4gas
#ifdef DIAG
        fxnucn(iv)=rnucnum
        fxnucm(iv)=rnucmas
       else
        fxnucn(iv)=r0
        fxnucm(iv)=r0
#endif
       endif
      enddo
#elif KULNEW
!
!---->calculate number and mass of nucleated particles per time
!     following Vehkamaki (2002)
!
!
!---->define local parameter
!
      parameter (so4gmin=1.d7,so4gmax=1.d18)
      parameter (tmin=195.d0)
!
!---->indices for so4gas and smallest particles
!
      inum=nso4*nmomso4
      ixnum=nnon+inum
      ixmas=ixnum-1
      igso4=ntot
!
!---->start loop over parcels
!
      do iv=1,nv
!
!---->calculate h2so4 mole/mass fraction in critical nucleus
!
       if (so4gas(iv)>so4gmin) then
        tempnew=max(tmin,temp(iv))
        so4new=min(so4gmax,so4gas(iv))
        so4nuc=so4fracnuc_f(tempnew,relhum(iv),so4new)
!$$$        so4mnuc=so4nuc*so4mol
!$$$     $        /(so4nuc*so4mol+(r1-so4nuc)*h2omol)
!$$$        so4mnuc=min(r1,max(r0,so4mnuc))
!
!---->calculate nucleation rate: rnucnum [#particles/m3]
!
        rnucnum=rnucnum_f(tempnew,relhum(iv),so4new,so4nuc)
!.deb   rnucnum=min(rnucmax,rnucnum)
!
!---->calculate the number of so4 molecules in the critical cluster
!
        rnucmas=so4nuc*rnucmas_f(tempnew,relhum(iv),so4new,so4nuc)
        rnucmas=rnucnum*min(rmasmax,max(rmasmin,rnucmas))
!
!---->add to derivatives
!
        fxaer(iv,ixmas)=fxaer(iv,ixmas)+rnucmas ! mass   (small part)
        fxaer(iv,ixnum)=fxaer(iv,ixnum)+rnucnum ! number (small part)
        fxaer(iv,igso4)=fxaer(iv,igso4)-rnucmas ! so4gas
#ifdef DIAG
        fxnucn(iv)=rnucnum
        fxnucm(iv)=rnucmas
       else
        fxnucn(iv)=r0
        fxnucm(iv)=r0
#endif
       endif
      enddo
#elif FITZGERALD
!
!---->Jaeker-Voirol Mirabel nucleation rate for sulfuric acid-water as
!     parameterized in Fitzgerald et al., JGR Vol.103,p.16085-16102,1998.
!
!     input variables:
!
!       so4gas   = H2SO4(g) concentration   [molec/m3]
!       temp     = temperature              [K]
!       relhum   = relative humidity        [1]
!

!
!---->define local arrays and parameter
!
      dimension xmas(101)
      data xmas /40*207.6778d0,  &
     &           129.3836d0,127.3543d0,125.4111d0,123.5449d0,121.7469d0,  &
     &           120.0091d0,118.3238d0,116.6841d0,115.0834d0,113.5156d0,  &
     &           111.9748d0,110.4558d0,108.9535d0,107.4633d0,105.9807d0,  &
     &           104.5019d0,103.0231d0,101.5410d0,100.0523d0, 98.5545d0,  &
     &            97.0450d0, 95.5217d0, 93.9826d0, 92.4263d0, 90.8513d0,  &
     &            89.2569d0, 87.6421d0, 86.0068d0, 84.3507d0, 82.6740d0,  &
     &            80.9772d0, 79.2609d0, 77.5261d0, 75.7740d0, 74.0061d0,  &
     &            72.2262d0, 69.9497d0, 67.9782d0, 66.2147d0, 64.5760d0,  &
     &            62.9890d0, 61.3897d0, 59.7228d0, 57.9416d0, 56.0092d0,  &
     &            53.8991d0, 51.5964d0, 49.0984d0, 46.4150d0, 43.5677d0,  &
     &            40.1265d0, 35.4424d0, 30.6926d0, 26.2674d0, 22.3686d0,  &
     &            19.0614d0, 16.3291d0, 14.1144d0, 12.3469d0, 10.9579d0,  &
     &             9.8881d0/
      save xmas

      parameter (r3=3.d0,r4=4.d0,r7=7.d0,r10=10.d0,r100=100.d0)
      parameter (e1=0.8352d0,e2=1.09058d-2,e3=1.59865d-4,  &
     &           e4=2.7221d0,e5=0.16203d0 ,e6=36.01547d0)
      parameter (d1=1.3720d0,d2=1.50374d-2,d3=1.09676d-4,  &
     &           d4=3.3753d0,d5=0.137314d0,d6=31.12747d0)
      parameter (chi11=3.3212d0 ,chi12=8.0094d-4 ,chi13=6.4893d-5,  &
     &           chi14=11.9669d0,chi15=0.038134d0,chi16=17.0665d0)
      parameter (chi21=0.5368d0 ,chi22=4.13158d-3,chi23=5.0189d-5,  &
     &           chi24=0.6532d0 ,chi25=0.024719d0,chi26=4.8880d0)
      parameter (chi31=4.4133d0 ,chi32=8.8952d-4 ,chi33=3.5587d-5,  &
     &           chi34=17.1074d0,chi35=0.023770d0,chi36=11.7168d0)
      parameter (epsi=1.d-22)
      parameter (expmax=10.d0)
!
!---->indices for so4gas and smallest particles
!
      inum=nso4*nmomso4
      ixnum=nnon+inum
      ixmas=ixnum-1
      igso4=ntot

      do iv=1,nv
!
!---->calculate mass of critical nucleus using the magic numbers
!     from Mod6M (Brookhaven Model)
!
       rh=r100*relhum(iv)
       irh=min(100,int(rh)+1)
       weight=real(irh,kdef)-rh
       rnucmas=xmas(irh)*weight+xmas(irh+1)*(r1-weight)
!
!---->calculate number nucleation rate
!
       temp2=temp(iv)**2
       xlogrh=log10(rh)
       xlogrh2=xlogrh**2
       xlogso4=log10(max(epsi,so4gas(iv))*so4mol/avno)

       evalue=xlogso4-e1*xlogrh2 +e2*temp(iv)*xlogrh  &
     &               +e3*temp2   +e4*xlogrh  &
     &               -e5*temp(iv)+e6
       if (evalue<=r0) then
        chi1=chi11*xlogrh2-chi12*temp(iv)*xlogrh  &
     &      +chi13*temp2  -chi14*xlogrh  &
     &      -chi15*temp(iv)+chi16
        expon=min(expmax,r4+r3*evalue/chi1)
        rnucnum=r10**expon
       else
        dvalue=xlogso4-d1*xlogrh2 +d2*temp(iv)*xlogrh  &
     &                +d3*temp2   +d4*xlogrh  &
     &                -d5*temp(iv)+d6
        if (dvalue<=r0) then
         chi2=chi21*xlogrh2-chi22*temp(iv)*xlogrh  &
     &       +chi23*temp2  -chi24*xlogrh  &
     &       -chi25*temp(iv)+chi26
         expon=min(expmax,r7+r3*dvalue/chi2)
         rnucnum=r10**expon
        else
         chi3=-chi31*xlogrh2-chi32*temp(iv)*xlogrh  &
     &        +chi33*temp2  +chi34*xlogrh  &
     &        -chi35*temp(iv)-chi36
         expon=min(expmax,r7+r3*dvalue/chi3)
         rnucnum=r10**expon
        endif
       endif

       rnucmas=rnucnum*rnucmas
!
!---->add to derivatives
!
       fxaer(iv,ixmas)=fxaer(iv,ixmas)+rnucmas ! mass   (small part)
       fxaer(iv,ixnum)=fxaer(iv,ixnum)+rnucnum ! number (small part)
       fxaer(iv,igso4)=fxaer(iv,igso4)-rnucmas ! so4gas
#ifdef DIAG
       fxnucn(iv)=rnucnum
       fxnucm(iv)=rnucmas
#endif
      enddo
#else
!
!---->calculate number and mass of nucleated particles per time
!     following classical nucleation theory
!     Zhao and Turco, 1995
!     J. Aerosol Sci. (26) pp 779-795
!
!     method: find so4 mass fraction with func(so4mfrac)=0.
!

!
!---->define local arrays and parameter
!
      dimension func(nv,3),smfrac(nv,3)
      parameter (r4=4.d0,r1h=0.5d0,r1td=1.d0/3.d0,r4td=4.d0/3.d0)
      parameter (smfracmin=0.00005d0,smfracmax=0.99995d0,  &
     &           smfracint=0.005d0)
      parameter (racmin=1.d-6)
      parameter (epsi=1.d-22)
!
!---->indices for so4gas and smallest particles
!
      inum=nso4*nmomso4
      ixnum=nnon+inum
      ixmas=ixnum-1
      igso4=ntot
!
!---->calculate function for min and max so4mfrac
!
      smfrac(:,1)=smfracmin
      smfrac(:,3)=smfracmax
      do j=1,3,2
       do iv=1,nv
        sfrac=smfrac(iv,j)*h2omol  &
     &      /(smfrac(iv,j)*h2omol+(r1-smfrac(iv,j))*so4mol)

        srho1=so4dens_f(temp(iv),smfrac(iv,j))
        dsmfrac=smfracint*sign(r1,smfrac(iv,j)-smfracint)
        srho2=so4dens_f(temp(iv),smfrac(iv,j)-dsmfrac)
        dlogsrho=log(srho1/srho2)/dsmfrac

        ssol=so4sol_f(temp(iv),so4sat(iv),sfrac)
        wsol=h2osol_f(temp(iv),h2osat(iv),sfrac)

        sgas=max(racmin*ssol,so4gas(iv))
        func(iv,j)=log(  sgas   /ssol)*(r1+dlogsrho*smfrac(iv,j))  &
     &           -log(h2ogas(iv)/wsol)*(r1-dlogsrho*(r1-smfrac(iv,j)))  &
     &                               *so4mol/h2omol
       enddo
      enddo
!
!---->search so4mfrac for which function equals zero
!
      do iv=1,nv
       if (func(iv,1)*func(iv,3)<=r0) then
        do iter=1,10
         smfrac(iv,2)=r1h*(smfrac(iv,3)+smfrac(iv,1))
         sfrac=smfrac(iv,2)*h2omol  &
     &       /(smfrac(iv,2)*h2omol+(r1-smfrac(iv,2))*so4mol)

         srho1=so4dens_f(temp(iv),smfrac(iv,2))
         dsmfrac=smfracint*sign(r1,smfrac(iv,2)-smfracint)
         srho2=so4dens_f(temp(iv),smfrac(iv,2)-dsmfrac)
         dlogsrho=log(srho1/srho2)/dsmfrac

         ssol=so4sol_f(temp(iv),so4sat(iv),sfrac)
         wsol=h2osol_f(temp(iv),h2osat(iv),sfrac)

         sgas=max(racmin*ssol,so4gas(iv))
         func(iv,2)=log(  sgas   /ssol)*(r1+dlogsrho*smfrac(iv,2))  &
     &            -log(h2ogas(iv)/wsol)*(r1-dlogsrho*(r1-smfrac(iv,2)))  &
     &                                 *so4mol/h2omol

         flag=r1h-sign(r1h,func(iv,2)*func(iv,3))
         flaginv=r1-flag

         func(iv,1)=func(iv,1)*flaginv+func(iv,2)*flag
         func(iv,3)=func(iv,2)*flaginv+func(iv,3)*flag

         smfrac(iv,1)=smfrac(iv,1)*flaginv+smfrac(iv,2)*flag
         smfrac(iv,3)=smfrac(iv,2)*flaginv+smfrac(iv,3)*flag
        enddo
!
!---->do linear interpolation
!
        smfrac(iv,2)=smfrac(iv,1)-(smfrac(iv,3)-smfrac(iv,1))  &
     &               *func(iv,1)/(  func(iv,3)  -func(iv,1))
        smfrac(iv,2)=max(smfrac(iv,1),min(smfrac(iv,3),smfrac(iv,2)))
        sfrac=smfrac(iv,2)*h2omol  &
     &      /(smfrac(iv,2)*h2omol+(r1-smfrac(iv,2))*so4mol)

        srho=so4dens_f(temp(iv),smfrac(iv,2))
        ssig=so4sig_f(temp(iv),sfrac,smfrac(iv,2))
        ssol=so4sol_f(temp(iv),so4sat(iv),sfrac)
        wsol=h2osol_f(temp(iv),h2osat(iv),sfrac)
!
!---->calculate radius of nucleus
!     limit radius so that there are at least 2 so4 molecules in nucleus
!
        factor=r4td*pi*srho*smfrac(iv,2)/so4mol*avno
        rmin=(rmasmin/factor)**r1td

        sgas=max(racmin*ssol,so4gas(iv))
        xtmp=srho*rgas*temp(iv)*(smfrac(iv,2)/so4mol*log(sgas/ssol)  &
     &               +(r1-smfrac(iv,2))/h2omol*log(h2ogas(iv)/wsol))
        xtmp=max(epsi,xtmp)
        radnuc=max(rmin,r2*ssig/xtmp)
        radnuc2=radnuc*radnuc

        gamma=r4td*pi*radnuc2*ssig/(boltz*temp(iv))
        beta=so4gas(iv)*sqrt(rgas*temp(iv)/(r2*pi*so4mol))
        rnucnum=r4*pi*radnuc2*beta*h2ogas(iv)*exp(-gamma)
!.deb   rnucnum=min(rnucmax,rnucnum)

        rnucmas=factor*radnuc2*radnuc
        rnucmas=rnucnum*min(rmasmax,max(rmasmin,rnucmas))
!
!---->add to derivatives
!
        fxaer(iv,ixmas)=fxaer(iv,ixmas)+rnucmas ! mass   (small part)
        fxaer(iv,ixnum)=fxaer(iv,ixnum)+rnucnum ! number (small part)
        fxaer(iv,igso4)=fxaer(iv,igso4)-rnucmas ! so4gas
#ifdef DIAG
        fxnucn(iv)=rnucnum
        fxnucm(iv)=rnucmas
       else
        fxnucn(iv)=r0
        fxnucm(iv)=r0
#endif
       endif
      enddo
#endif
      return
      end
!-----------------------------------------------------------------------
      subroutine condensso4(nv,so4gas,so4aer,so4radv,so4knu,dt,fxaer)
!
!---->calculte condensation of h2so4
!     on pre-existing sulfat aerosol including Kelvin effect,
!     assuming a lognormal particle distribution with prescribed sigma,
!     correction for collision geometry follows Fuchs, Sutugin (1971),
!     see also Jacobson (1999), Fundamentals of Atmospheric Modeling, p. 457.
!     the accommodation coefficient eta:
!          if -DLOWACCOEF 0.04 (van Dingenen and Raes, 1991)
!          else           1.00
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy, input/output arrays and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension so4gas(nvec),so4aer(nvec,naer)
      dimension so4radv(nv,nso4),so4knu(nv,nso4)
      dimension fxaer(nv,ntot),dt(nv)

      parameter (r4=4.d0)
      parameter (a1=1.333d0,a2=0.71d0)
      parameter (b1=0.572127d0)
#ifdef DIAG
      fxcon(1:nv)=r0
#endif
!
!---->index for so4gas
!
      igso4=ntot
!
!---->start loop over bins
!
      do iso4=1,nso4
       inum=iso4*nmomso4
       imas=inum-1
       ixnum=nnon+inum
       ixmas=ixnum-1
!
!---->start loop over parcels
!     Euler backward used to damp overshooting of an explicit time scheme
!
       do iv=1,nv

        sulknu=b1*so4knu(iv,iso4)
!$$$        alpha=exp(-log(so4sig(iv,iso4))**2)
        beta=so4alpha(iv,iso4)/(r1+((a1*sulknu+a2)  &
     &      /(r1+sulknu)+a1*(r1-etaso4(iso4))/etaso4(iso4))*sulknu)
        so4supsat=so4gas(iv)  &
     &           -so4sol(iv)*exp(so4kelvin(iv)/so4radv(iv,iso4))
        coefso4=beta*diffso4(iv)*so4radv(iv,iso4)*so4aer(iv,inum)
        so4cond=coefso4/(r1+dt(iv)*coefso4)*so4supsat

        so4part=max(r1,so4aer(iv,imas)/max(epsilon,so4aer(iv,inum)))
!
!---->add to derivatives
!     reduce number concentration in case of evaporation
!
        fxaer(iv,ixmas)=fxaer(iv,ixmas)+so4cond               ! mass   (bin iso4)
        fxaer(iv,ixnum)=fxaer(iv,ixnum)+min(r0,so4cond)         & ! number (bin iso4)
     &                        /(so4part*max(r1,so4part-r4))
        fxaer(iv,igso4)=fxaer(iv,igso4)-so4cond               ! so4gas
#ifdef DIAG
        fxcon(iv)=fxcon(iv)+so4cond
#endif
       enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine condensnon(nv,so4gas,so4non,dt,fxaer)
!
!---->calculte condensation of h2so4
!     on non-sulfat aerosol including Kelvin effect,
!     assuming a lognormal particle distribution with prescribed sigma,
!     correction for collision geometry follows Fuchs, Sutugin (1971),
!     see also Jacobson (1999), Fundamentals of Atmospheric Modeling, p. 457.
!     the accommodation coefficient eta:
!          if -DLOWACCOEF 0.04 (van Dingenen and Raes, 1991)
!          else           1.00
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy, input/output arrays and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension so4gas(nvec),so4non(nvec,nnon)
      dimension fxaer(nv,ntot),dt(nv)

      parameter (r1h=0.5d0)
      parameter (ralim=10.d0)
#ifdef DIAG
        fxconn(1:nv)=r0
#endif
!
!---->index for so4gas
!
      igso4=ntot
!
!---->start loop over non-sulfate aerosol
!
      do inon=1,nnon
       do iv=1,nv
        so4supsat=so4gas(iv)-so4sol(iv)*xnonkelvin(iv,inon)
        if (so4supsat>r0) then
         eta=min(r1,ralim*so4supsat/max(r1,so4gas(iv)))
         coefnon=condnon(iv,inon)*min(max(r1h*so4non(iv,inon),  &
     &                                    eta*xnonum(iv,inon)),  &
     &                                        xnonum(iv,inon))
        else
         coefnon=condnon(iv,inon)*min(r1h*so4non(iv,inon),  &
     &                                    xnonum(iv,inon))
        endif
        cond=coefnon/(r1+dt(iv)*coefnon)*so4supsat
!
!---->add to derivatives
!
        fxaer(iv,inon )=fxaer(iv,inon )+cond                    ! non-so4
        fxaer(iv,igso4)=fxaer(iv,igso4)-cond                    ! so4gas
#ifdef DIAG
        fxconn(iv)=fxconn(iv)+cond
#endif
       enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine coagulation(nv,temp,so4gas,so4aer,so4non,  &
     &                       so4radv,so4knu,so4vel,fxaer)
!
!---->calculate coagulation of sulfate aerosol
!     following Seinfeld (1997) table 12.1
!     using Fuchs form of Brownian coagulation coefficient
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension temp(nv)
      dimension so4gas(nvec),so4aer(nvec,naer),so4non(nvec,nnon)
      dimension so4radv(nv,nso4),so4knu(nv,nso4),so4vel(nv,nso4)
      dimension fxaer(nv,ntot)
!
!---->define local arrays and parameter
!
      dimension diff(nv,nso4),glam(nv,nso4)
      dimension so4coag(nv,nso4,nso4)

      parameter (r3=3.d0,r4=4.d0,r6=6.d0,r8=8.d0,r1h=0.5d0,r3h=1.5d0)
      parameter (a1=1.249d0,a2=0.42d0,a3=-0.87d0)
#ifdef DIAG
      fxcoag(1:nv)=r0
#endif
!
!---->precalculate diffusion coefficient for particle,
!     characteristic length scale (glam)
!
      do iso4=1,nso4
       do iv=1,nv
        diff(iv,iso4)=r4*boltz*temp(iv)  &
     &            /(r6*pi*viscos(iv)*so4radv(iv,iso4))  &
     &            *(r1+so4knu(iv,iso4)*(a1+a2*exp(a3/so4knu(iv,iso4))))
        xlam=diff(iv,iso4)/(pi*sqrt(so4vel(iv,iso4)))
        glam(iv,iso4)=r2/(r3*so4radv(iv,iso4)*xlam)  &
     &                *((so4radv(iv,iso4)+xlam)**3  &
     &                 -(so4radv(iv,iso4)**2+xlam**2)**r3h)  &
     &               -r2*so4radv(iv,iso4)
        glam(iv,iso4)=glam(iv,iso4)*glam(iv,iso4)
       enddo
      enddo
!
!---->coagulation coefficients including kinetic effects
!
      do jso4=1,nso4
       do iso4=jso4,nso4
        do iv=1,nv
         radsum=so4radv(iv,iso4)+so4radv(iv,jso4)
         diffsum=diff(iv,iso4)+diff(iv,jso4)
         so4coag(iv,iso4,jso4)=pi*diffsum*radsum  &
     &        /(radsum/(radsum+sqrt(glam(iv,iso4)+glam(iv,jso4)))  &
     &         +diffsum/(sqrt(so4vel(iv,iso4)+so4vel(iv,jso4))*radsum))
        enddo
       enddo
      enddo
!
!---->autocoagulation
!
      do iso4=1,nso4
       inum=iso4*nmomso4
       imas=inum-1
       ixnum=nnon+inum
       ixmas=ixnum-1
       do iv=1,nv
        dnumcoag=r1h*so4coag(iv,iso4,iso4)*so4aer(iv,inum)  &
     &                                    *so4aer(iv,inum)
!
!---->add to derivatives
!
        fxaer(iv,ixnum)=fxaer(iv,ixnum)-dnumcoag         ! number (bin iso4)
#ifdef DIAG
        fxcoag(iv)=fxcoag(iv)-dnumcoag
#endif
       enddo
      enddo
!
!---->(non-auto) coagulation
!
      do jso4=1,nso4
       jnum=jso4*nmomso4
       jmas=jnum-1
       jxnum=nnon+jnum
       jxmas=jxnum-1
       do iso4=jso4+1,nso4
        inum=iso4*nmomso4
        imas=inum-1
        ixnum=nnon+iso4*nmomso4
        ixmas=ixnum-1
        do iv=1,nv
         dmascoag=so4coag(iv,iso4,jso4)*so4aer(iv,imas)*so4aer(iv,jnum)
         dnumcoag=so4coag(iv,iso4,jso4)*so4aer(iv,inum)*so4aer(iv,jnum)
!
!---->add to derivatives
!
         fxaer(iv,jxmas)=fxaer(iv,jxmas)+dmascoag        ! mas    (bin jso4)
         fxaer(iv,ixmas)=fxaer(iv,ixmas)-dmascoag        ! mass   (bin iso4)
         fxaer(iv,ixnum)=fxaer(iv,ixnum)-dnumcoag        ! number (bin iso4)
#ifdef DIAG
         fxcoag(iv)=fxcoag(iv)-dnumcoag
#endif
        enddo
       enddo
      enddo
!
!---->calculate coagulation coefficients
!     for coagulation between sulfate aerosol and non-sulfate aerosol
!
#ifndef NONON
#ifdef DIAG
      fxcoagn(1:nv)=r0
#endif
      igso4=ntot
      do inon=1,nnon
       xvelnon=r8*boltz/(pi*pmsnon(inon))
       do iso4=1,nso4
        inum=iso4*nmomso4
        imas=inum-1
        ixnum=nnon+inum
        ixmas=ixnum-1
        do iv=1,nv
         factor=radvnon(inon)/rwetnon(iv,inon)
         velnon=xvelnon*temp(iv)
         xlamnon=diffnon(iv,inon)/(pi*factor*sqrt(velnon))
         glamnon=r2/(r3*rwetnon(iv,inon)*xlam)  &
     &               *((rwetnon(iv,inon)+xlam)**3  &
     &                -(rwetnon(iv,inon)**2+xlam**2)**r3h)  &
     &              -r2*rwetnon(iv,inon)
         glamnon=glamnon*glamnon

         radsum=so4radv(iv,iso4)+rwetnon(iv,inon)
         diffsum=diff(iv,iso4)+diffnon(iv,inon)
         so4noncoag=pi*diffsum*radsum  &
     &        /(radsum/(radsum+sqrt(glam(iv,iso4)+glamnon))  &
     &         +diffsum/(radsum*(so4vel(iv,iso4)+velnon)))
!
!---->coagulation of sulfate aerosol and non-sulfate aerosol
!
!$$$         xnum=min(max(r1,r1h*so4non(iv,inon)),xnonum(iv,inon))
         xnum=xnonum(iv,inon)
         dmascoag=so4noncoag*so4aer(iv,imas)*xnum
         dnumcoag=so4noncoag*so4aer(iv,inum)*xnum
!
!---->add to derivatives
!     assume immediate evaporation if supsat<0
!
         so4supsat=so4gas(iv)-so4sol(iv)*xnonkelvin(iv,inon)
         if (so4supsat>r0) then
          fxaer(iv,inon) =fxaer(iv,inon) +dmascoag       ! mas    (bin inon)
         else
          fxaer(iv,igso4)=fxaer(iv,igso4)+dmascoag
         endif
         fxaer(iv,ixmas)=fxaer(iv,ixmas)-dmascoag        ! mass   (bin iso4)
         fxaer(iv,ixnum)=fxaer(iv,ixnum)-dnumcoag        ! number (bin iso4)
#ifdef DIAG
         fxcoagn(iv)=fxcoagn(iv)-dnumcoag
#endif
        enddo
       enddo
      enddo
#endif
      return
      end
!-----------------------------------------------------------------------
      subroutine sources(nv,so4aer,so4radv,srcgas,srcmas,srcpar,fxaer)
!
!---->external sources due to clear-sky and aqeous chemistry,
!     and direct particle emisssions due to fast near source conversion
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension so4aer(nvec,naer),so4radv(nv,nso4)
      dimension srcgas(nv),srcmas(nv),srcpar(nv)
      dimension fxaer(nv,ntot)
!
!---->define local arrays and parameter
!
      dimension srcm(nso4),srcn(nso4),fsrc(nso4),  &
     &          srcmnon(nnon),fsrcnon(nnon)
!
!---->partitioning between small and large particles
!     in near SO2 source regions
!
      parameter (rat1=0.85d0,rat2=r1-rat1)
!
!---->indices for so4gas etc.
!
      igso4=ntot
      do iv=1,nv
       srcn(:)=r0
!
!---->aqueous chemistry source feeds accumulation mode
!     if non-existent feeds so4gas
!
#ifdef NONON
       total=sum(fracsrc(iv,1:nso4))
       if (total>r1)then
        fsrc(1:nso4)=fracsrc(iv,1:nso4)/total
        srcm(1:nso4)=fsrc(1:nso4)*srcmas(iv)
        srcg=r0
       else
        srcg=srcmas(iv)
        srcm(:)=r0
       endif
#else
       total=sum(fracsrc(iv,1:nso4))+sum(fracsrcnon(iv,1:nnon))
       if (total>r1)then
        fsrc(1:nso4)=fracsrc(iv,1:nso4)/total
        fsrcnon(1:nnon)=fracsrcnon(iv,1:nnon)/total
        srcm(1:nso4)=fsrc(1:nso4)*srcmas(iv)
        srcmnon(1:nnon)=fsrcnon(1:nnon)*srcmas(iv)
        srcg=r0
       else
        srcg=srcmas(iv)
        srcm(:)=r0
        srcmnon(:)=r0
       endif
#endif
!
!---->gasphase chemistry source feeds so4gas
!
       srcg=srcg+srcgas(iv)
!
!---->particulate source in near SO2 source regions
!
       xmas=rat2*srcpar(iv)
       srcm(nso4)=srcm(nso4)+xmas
       srcn(nso4)=srcn(nso4)+xmas/pmssrc(1)

       xmas=rat1*srcpar(iv)
       srcm(nso4-1)=srcm(nso4-1)+xmas
       srcn(nso4-1)=srcn(nso4-1)+xmas/pmssrc(2)
!
!---->add to derivatives
!
#ifndef NONON
       do inon=1,nnon
        fxaer(iv,inon)=fxaer(iv,inon)+srcmnon(inon)        ! non-so4
       enddo
#endif
       do iso4=1,nso4
        ixnum=nnon+iso4*nmomso4
        ixmas=ixnum-1
        fxaer(iv,ixmas)=fxaer(iv,ixmas)+srcm(iso4)         ! mass
        fxaer(iv,ixnum)=fxaer(iv,ixnum)+srcn(iso4)         ! number
       enddo

       fxaer(iv,igso4)=fxaer(iv,igso4)+srcg                ! so4gas

      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine settling(nv,so4aer,so4radv,fxaer)
!
!---->gravitational settling
!     Stokes law with Cunningham slip correction factor
!     see Seinfeld (1997) eq. 8.34 and eq. 8.42
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy, input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension so4aer(nvec,naer),so4radv(nv,nso4)
      dimension fxaer(nv,ntot)
!
!---->define local parameter
!
      parameter (r1h=0.5d0,r5h=2.5d0)
      parameter (a1=2.d0/9.d0)
      parameter (b1=1.257d0,b2=0.4d0,b3=1.1d0)
      parameter (r3=3.d0)
      parameter (dzinv=r1/65.d0)
#ifdef DIAG
      fxgrvn(1:nv)=r0
      fxgrvm(1:nv)=r0
#endif
!
!---->start loop over bins
!
      do iso4=1,nso4
       inum=iso4*nmomso4
       imas=inum-1
       ixnum=nnon+inum
       ixmas=ixnum-1
!
!---->start loop over parcels
!
       do iv=1,nv
        xlnsg=log(so4sig(iv,iso4))
        xlnsg2=xlnsg*xlnsg
!
!---->calculate terminal fall velocities for mass and number
!
        radm=so4radv(iv,iso4)*exp(r5h*xlnsg2)
        so4knum=pathair(iv)/radm
        wmas=a1*grav/viscos(iv)*radm*radm*so4dens(iv)  &
     &      *(r1+so4knum*(b1+b2*exp(-b3/so4knum)))

        radn=so4radv(iv,iso4)*exp(-r1h*xlnsg2)
        so4knun=pathair(iv)/radn
        wnum=a1*grav/viscos(iv)*radn*radn*so4dens(iv)  &
     &      *(r1+so4knun*(b1+b2*exp(-b3/so4knun)))
!
!---->add to derivatives
!
        fxaer(iv,ixmas)=fxaer(iv,ixmas)-so4aer(iv,imas)        & ! mass  (bin iso4)
     &                                 *wmas*dzinv
        fxaer(iv,ixnum)=fxaer(iv,ixnum)-so4aer(iv,inum)        & ! nuber (bin iso4)
     &                                 *wnum*dzinv
#ifdef DIAG
        fxgrvm(iv)=fxgrvm(iv)-so4aer(iv,imas)*wmas*dzinv
        fxgrvn(iv)=fxgrvn(iv)-so4aer(iv,inum)*wnum*dzinv
#endif
       enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine merge(ivec,xaer)
!
!---->merging
!     move particles greater radmerg in diameter to next mode
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension xaer(nvec,ntot)
      target xaer
!
!---->define local arrays and parameter
!
      dimension so4radv(nso4)

      parameter (r8=8.d0,r1h=0.5d0,r3h=1.5d0,  &
     &           r3q=0.75d0,r1td=1.d0/3.d0)
!
!---->define pointer
!
      dimension so4aer(:)
      pointer so4aer
!
!---->set pointer
!
      so4aer=>xaer(ivec,nnon+1:nnon+naer)
!
!---->calculate wet aerosol size (Tang and Munkelwitz (1994)),
!
      do iso4=1,nso4
       inum=iso4*nmomso4
       imas=inum-1
       wetmas=max(r2,so4aer(imas)  &
     &       /max(epsilon/pmsmerg(ivec,iso4+1),so4aer(inum)))  &
     &       *so4min/so4mfrac(ivec)
       so4radv(iso4)=(r3q*wetmas/(pi*so4dens(ivec)))**r1td
      enddo
!
!---->start loop over aerosol bins
!
      do iso4=2,nso4
       xlnsg=log(so4sig(ivec,iso4))
       drgvdrv=exp(r3h*xlnsg*xlnsg)
       drgdrv=r1/drgvdrv
       radg=so4radv(iso4)*drgdrv
       if (radg>radgmerg(1,iso4)) then
        inum=iso4*nmomso4
        imas=inum-1
        inmerg=inum-nmomso4
        immerg=inmerg-1
!
!---->calculate change in mass and number due to merging
!
        radgv=so4radv(iso4)*drgvdrv
        kmergn=1
        kmergv=1
        do k=1,nmerg
         if (radgmerg(k,iso4)<radg)  kmergn=k
         if (radgmerg(k,iso4)<radgv) kmergv=k
        enddo

        k=kmergv
        kp1=k+1
        dmmerg=(fracmerg(k  )*(radgmerg(kp1,iso4)-radgv)  &
     &         -fracmerg(kp1)*(radgmerg(k  ,iso4)-radgv))  &
     &        /(radgmerg(kp1,iso4)-radgmerg(k,iso4))
        dmmerg=min(r1,dmmerg)*so4aer(imas)
        k=kmergn
        kp1=k+1
        dnmerg=(fracmerg(k  )*(radgmerg(kp1,iso4)-radg)  &
     &         -fracmerg(kp1)*(radgmerg(k  ,iso4)-radg))  &
     &        /(radgmerg(kp1,iso4)-radgmerg(k,iso4))
        dnmerg=min(r1,dnmerg)*so4aer(inum)
!
!---->add changes (simple Euler backward)
!
        so4aer(imas)  =so4aer(imas)  -dmmerg
        so4aer(immerg)=so4aer(immerg)+dmmerg

        if (so4aer(imas)==r0) dnmerg=so4aer(inum)
        so4aer(inum)  =so4aer(inum)  -dnmerg
        so4aer(inmerg)=so4aer(inmerg)+dnmerg
       endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine accumode(nv,so4aer,so4radv,so4non)
!
!---->calculate fraction of particles in accumulation mode
!
!
!---->load modules
!
      use precision
      use modpar_umaer
      use modbox_umaer
!
!---->define accuracy and input/output arrays
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      dimension so4aer(nvec,naer),so4radv(nv,nso4)
      dimension so4non(nvec,nnon)
!
!---->define local arrays and parameter
!
      parameter (r3h=1.5d0,r5=5.d0)
!
!---->add for non-sulfate aerosol that can act as CCN
#ifndef NONON
      parameter (r3q=0.75d0,r1td=1.d0/3.d0)
      parameter (fso4crit=10.d0)
      dimension xso4non(nv,nnon),fso4non(nv,nnon)
      dimension rso4min(nv)
#endif
!
!---->sulfate aerosol
!
      do iso4=1,nso4
       inum=iso4*nmomso4
       do iv=1,nv
        radg=so4radv(iv,iso4)*exp(-r3h*log(so4sig(iv,iso4))**2)
        if (r5*so4sig(iv,iso4)*radg>radaccu(1)) then
         if (radg<racculim) then
          do iaccu=1,naccu
           if (radgaccu(iaccu,iso4,1)>radg) exit
          enddo
          k=iaccu
          km1=k-1
          fracsrc(iv,iso4)=so4aer(iv,inum)  &
     &         *(fracaccu(km1,iso4)*(radgaccu(k  ,iso4,1)-radg)  &
     &          -fracaccu(k  ,iso4)*(radgaccu(km1,iso4,1)-radg))  &
     &         /(radgaccu(k  ,iso4,1)-radgaccu(km1,iso4,1))
         else
          do iaccu=1,naccu
           if (radgaccu(iaccu,iso4,2)<=radg) exit
          enddo
          k=iaccu
          km1=k-1
          fracsrc(iv,iso4)=so4aer(iv,inum)  &
     &         *(fracaccu(km1,iso4)*(radgaccu(k  ,iso4,2)-radg)  &
     &          -fracaccu(k  ,iso4)*(radgaccu(km1,iso4,2)-radg))  &
     &         /(radgaccu(k  ,iso4,2)-radgaccu(km1,iso4,2))
         endif
        else
         fracsrc(iv,iso4)=r0
        endif
       enddo
      enddo
#ifndef NONON
!
!---->non-sulfate aerosol
!
      do inon=1,nnon
       do iv=1,nv
        radg=radvnon(inon)*exp(-r3h*log(xnonsig(iv,inon))**2)
        if (xnonsig(iv,inon)<=r1) then
          if (radg>radaccu(1)) then
           fracsrcnon(iv,inon)=xnonum(iv,inon)
          else
           fracsrcnon(iv,inon)=r0
          endif
        else
          if (r5*xnonsig(iv,inon)*radg>radaccu(1)) then
            do iaccu=1,naccu
             if (radgaccnon(iaccu,inon)>radg) exit
            enddo
            k=iaccu
            km1=k-1
            fracsrcnon(iv,inon)=xnonum(iv,inon)  &
     &           *(fraccnon(km1)*(radgaccnon(k  ,inon)-radg)  &
     &            -fraccnon(k  )*(radgaccnon(km1,inon)-radg))  &
     &           /(radgaccnon(k,inon)-radgaccnon(km1,inon))
          else
           fracsrcnon(iv,inon)=r0
          endif
        endif
       enddo
      enddo
!
!---->so4 molecules on non-so4 aerosol surface (#molec/particle)
      do inon=1,nnon
        do iv=1,nv
          xso4non(iv,inon)=so4non(iv,inon)/  &
     &                     max(r1,xnonum(iv,inon))
        enddo
      enddo
!---->radius of each so4 molecule (m)
      do iv=1,nv
        rso4min(iv)=(so4min*r3q /  &
     &              (pi*so4mfrac(iv)*so4dens(iv)))**r1td
      enddo
!---->area fraction of so4 molecules on non-so4 aerosol surface
      do inon=1,nnon
        do iv=1,nv
          fso4non(iv,inon)=rso4min(iv)**2.d0*xso4non(iv,inon)  &
     &                    /(4.d0*(radvnon(inon)+rso4min(iv))**2.d0)
        enddo
      enddo
!---->non-sulfate aerosol that can act as CCN
      do inon=1,nnon
        do iv=1,nv
          if (inon<=nnon-4) then
            fracsrcnon(iv,inon)=fracsrcnon(iv,inon)  &
     &                         *min(fso4non(iv,inon)/fso4crit,r1)
          endif
        enddo
      enddo
#endif
      return
      end
!-----------------------------------------------------------------------
      function h2osat_f(temp)
!
!---->number of h2o molecules in saturation
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->Preining et al. (1981), report of University of Vienna
!     used in Kulmala et al. (1998) for 233K<T<298K
!
      parameter (a1=77.34491296d0,  &
     &           a2=7235.424651d0,  &
     &           a3=8.2d0,  &
     &           a4=5.7113d-3)
!
!---->Tabazadeh, Toon, Clegg, and Hamill (1997), GRL (24) pp 1931-1934
!     based on data from
!     Clegg and Brimblecombe (1995), J. Chem.Eng. Data (40) 43
!     185K<T<260
!
      parameter (b1=18.452406985d0+4.605170185988092d0,  &
     &           b2=3505.1578807d0,  &
     &           b3=330918.55082d0,  &
     &           b4=12725068.262d0)

      parameter (templim=229.d0)

      if (temp<templim) then
       temp2=temp*temp
       temp3=temp2*temp
       h2osat_f=exp(b1-b2/temp-b3/temp2+b4/temp3)/(boltz*temp)
      else
       h2osat_f=exp(a1-a2/temp-a3*log(temp)+a4*temp)/(boltz*temp)
      endif

      return
      end
!-----------------------------------------------------------------------
      function so4sat_f(temp)
!
!---->number of h2so4 molecules in saturation
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
#ifndef PARNEW
!
!---->Ayers et.al. (1980), GRL (7) pp 433-436
!
      parameter (a1=27.78492066d0,  &
     &           a2=10156.0d0)

      so4sat_f=exp(a1-a2/temp)/(boltz*temp)
#else
!
!---->Ayers et.al. (1980), GRL (7) pp 433-436
!     plus corrections for lower temperatures by Kulmala and Laaksonen (1990)
!     and Noppel et al. (1990)
!
!     ATTENTION! in Noppel et al. (1990) b2=11.94 which fits better
!                plot in Kulmala and Laaksonen (1990)

      parameter (b1=1.01325d5,  &
     &           b2=11.695d0,  &
     &           b3=1.0156d4,  &
     &           b4=0.38d0/545.d0,  &
     &           tref=360.15d0)
      so4sat_f=b1*exp(-b2+b3*(r1/tref-r1/temp  &
     &                        +b4*(r1+log(tref/temp)-tref/temp)))  &
     &        /(boltz*temp)
#endif
      return
      end
!-----------------------------------------------------------------------
      function h2osol_f(temp,h2osat,so4frac)
!
!---->calculate number concentration over the solution for h2o
!     Taleb et.al. (1996), JGR (101) pp 25967-25977
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (r10=10.0d0)
      parameter (a1=2989.d0,  &
     &           a2=2147000.d0,  &
     &           a3=233000000.d0,  &
     &           b11=0.527d0)
      parameter (tempmin=190.d0,tempmax=298.d0)

      tempnew=max(tempmin,min(tempmax,temp))
      h2ofrac=max(epsilon,r1-so4frac)
      a11=a1-a2/tempnew+a3/(tempnew*tempnew)
      h2ogam=r10**(a11*so4frac*so4frac  &
     &              /(tempnew*(so4frac+b11*h2ofrac)**2))
      h2osol_f=h2ogam*h2ofrac*h2osat

      return
      end
!-----------------------------------------------------------------------
      function so4sol_f(temp,so4sat,so4frac)
!
!---->calculate number concentration over the solution for h2so4
!     Taleb et al. (1996), JGR (101) pp 25967-25977
!     based on Zeleznik (1991) data
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (r10=10.0d0)
      parameter (a1=5672.d0,  &
     &           a2=4074000.d0,  &
     &           a3=4.421d8,  &
     &           b12=1.8975332d0)
      parameter (tempmin=190.d0,tempmax=298.d0)

      tempnew=max(tempmin,min(tempmax,temp))
      h2ofrac=r1-so4frac
      a12=a1-a2/tempnew+a3/(tempnew*tempnew)
      so4gam=r10**(a12*h2ofrac*h2ofrac  &
     &              /(tempnew*(h2ofrac+b12*so4frac)**2))
      so4sol_f=so4gam*so4frac*so4sat

      return
      end
!-----------------------------------------------------------------------
#ifndef KULNEW
      function so4fracnuc_f(temp,relhum,relacd,h2ogas,so4gas)
!
!---->calculate h2so4 molfraction of critical nucleus
!     Kulmala et.al. (1998), JGR (103) pp 8301-8307
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->Kulmala et.al. (1998), JGR (103) pp 8301-8307
!
      parameter (a1=1.6557d0,  &
     &           a2=0.0154d0,  &
     &           a3=0.0102d0,  &
     &           a4=0.0415d0,  &
     &           a5=0.0016d0)

      parameter (fracmin=0.000d0,fracmax=0.995d0)

      so4fracnuc_f=min(r1,max(r0,a1-a2*relacd/max(epsilon,relacd+relhum)  &
     &          +a3*log(max(epsilon,so4gas))-a4*log(max(epsilon,h2ogas))  &
     &          +a5*temp))
      so4fracnuc_f=max(fracmin,min(fracmax,so4fracnuc_f))
      return
      end
!-----------------------------------------------------------------------
#else
      function so4fracnuc_f(temp,relhum,so4gas)
!
!---->calculate h2so4 molfraction of critical nucleus
!     Vehkamaki (2002)
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
!
!---->Vehkamaeki et al. (2002)
!
      parameter (a1=7.40997d-1,  &
     &           a2=2.66379d-3,  &
     &           a3=3.49998d-3,  &
     &           a4=5.04022d-5,  &
     &           a5=2.01048d-3,  &
     &           a6=1.83289d-4,  &
     &           a7=1.57407d-3,  &
     &           a8=1.79059d-5,  &
     &           a9=1.84403d-4,  &
     &           a10=1.50345d-6)

      parameter (so4gmin=1.d10,so4gmax=1.d17)
      parameter (convfac=1.d-6)

      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas))*convfac)
      so4fracnuc_f=a1       -a2*temp  &
     &            -a3*so4ln +a4*temp*so4ln  &
     &            +a5*rhln  -a6*temp*rhln  &
     &            +a7*rhln2 -a8*temp*rhln2  &
     &            +a9*rhln3-a10*temp*rhln3

      return
      end
#endif
!-----------------------------------------------------------------------
      function so4crit_f(temp,relhum)
!
!
!---->calculate sulfuric acid concentration needed to produce 1 cm-3 s-1
!     so4crit in [#so4 molecules/m3]
!     Kulmala et.al. (1998), JGR (103) pp 8301-8307
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (a1=0.69699d0,  &
     &           a2=0.1335d0,  &
     &           a3=10.5462d0,  &
     &           a4=1958.4d0)

      parameter (tmin=233.d0,tmax=299.d0,  &
     &           rhmin=0.1d0,rhmax=1.d0)
!
!---->limit temperature and relative humidity
!     233K < temp   < 299K
!     0.1  < relhum < 1.0
!
!
      tempnuc=max(tmin,min(tmax,temp))
      relhnuc=max(rhmin,min(rhmax,relhum))
      so4crit_f=exp(-a1+a2*tempnuc-a3*relhnuc+a4*relhnuc/tempnuc)

      return
      end
!-----------------------------------------------------------------------
#ifndef KULNEW
      function rnucnum_f(temp,relhum,so4gas,so4frac,so4crit)
!
!---->calculate nucleation rate: rnucnum [#particles/m3]
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
!---->Kulmala et.al. (1998), JGR (103) pp 8301-8307
!
!---->WARNING! so4crit_f(temp,relhum) needs to be called first,
!     temperature and relative humidity should be limited to
!     233K < temp   < 299K
!     0.1  < relhum < 1.0
!
!
!---->local parameter
!
      parameter (b1=13.815511d0,  &
     &           b2=25.1289d0,  &
     &           b3=4890.80d0,  &
     &           b4=1743.3d0,  &
     &           b5=2.2479d0,  &
     &           b6=7643.3d0,  &
     &           b7=1.9712d0)
      parameter (t0=273.15d0)

      dnum=log(max(epsilon,so4gas)/so4crit)
      dtem=r1+(temp-t0)/t0

      rnucnum_f=exp(b1+b2*dnum-b3*dnum/temp-b4/temp  &
     &             -b5*dtem*dnum*relhum+b6*so4frac/temp  &
     &             -b7*so4frac*dtem/relhum)
      return
      end
!-----------------------------------------------------------------------
#else
      function rnucnum_f(temp,relhum,so4gas,so4frac)
!
!---->calculate nucleation rate: rnucnum [#particles/m3]
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
!---->Vehkamaeki et al. (2002)
!
      parameter (a1=0.14309d0,  &
     &           a2=2.21956d0,  &
     &           a3=2.73911d-2,  &
     &           a4=7.22811d-5,  &
     &           a5=5.91822d0)
      parameter (b1=0.117489d0,  &
     &           b2=0.462532d0,  &
     &           b3=1.18059d-2,  &
     &           b4=4.04196d-5,  &
     &           b5=15.7963d0)
      parameter (c1=0.215554d0,  &
     &           c2=8.10269d-2,  &
     &           c3=1.43581d-3,  &
     &           c4=4.7758d-6,  &
     &           c5=2.91297d0)
      parameter (d1=3.58856d0,  &
     &           d2=4.9508d-2,  &
     &           d3=2.1382d-4,  &
     &           d4=3.10801d-7,  &
     &           d5=2.93333d-2)
      parameter (e1=1.14598d0,  &
     &           e2=6.00796d-1,  &
     &           e3=8.64245d-3,  &
     &           e4=2.28947d-5,  &
     &           e5=8.44985)
      parameter (f1=2.15855d0,  &
     &           f2=8.08121d-2,  &
     &           f3=4.07382d-4,  &
     &           f4=4.01947d-7,  &
     &           f5=7.21326d-1)
      parameter (g1=1.6241d0,  &
     &           g2=1.60106d-2,  &
     &           g3=3.77124d-5,  &
     &           g4=3.21794d-8,  &
     &           g5=1.13255d-2)
      parameter (h1=9.71682d0,  &
     &           h2=1.15048d-1,  &
     &           h3=1.57098d-4,  &
     &           h4=4.00914d-7,  &
     &           h5=0.71186d0)
      parameter (p1=1.05611d0,  &
     &           p2=9.03378d-3,  &
     &           p3=1.98417d-5,  &
     &           p4=2.46048d-8,  &
     &           p5=5.79087d-2)
      parameter (q1=0.148712d0,  &
     &           q2=2.83508d-3,  &
     &           q3=9.24619d-6,  &
     &           q4=5.00427d-9,  &
     &           q5=1.27081d-2)
      parameter (convfac=1.d-6,convinv=1.d6)
      parameter (expmax=46.d0)

      sfracinv=r1/so4frac
      temp2=temp*temp
      temp3=temp*temp2

      a=+a1+a2*temp-a3*temp2+a4*temp3+a5*sfracinv
      b=+b1+b2*temp-b3*temp2+b4*temp3+b5*sfracinv
      c=-c1-c2*temp+c3*temp2-c4*temp3-c5*sfracinv
      d=-d1+d2*temp-d3*temp2+d4*temp3-d5*sfracinv
      e=+e1-e2*temp+e3*temp2-e4*temp3-e5*sfracinv
      f=+f1+f2*temp-f3*temp2-f4*temp3+f5*sfracinv
      g=+g1-g2*temp+g3*temp2+g4*temp3-g5*sfracinv
      h=+h1-h2*temp+h3*temp2+h4*temp3+h5*sfracinv
      p=-p1+p2*temp-p3*temp2+p4*temp3-p5*sfracinv
      q=-q1+q2*temp-q3*temp2+q4*temp3-q5*sfracinv


      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=max(r0,log(so4gas*convfac))
      so4ln2=so4ln*so4ln
      so4ln3=so4ln*so4ln2

      expon=min(expmax,a+b*rhln+c*rhln2+d*rhln3  &
     &                +e*so4ln+f*rhln*so4ln+g*rhln2*so4ln  &
     &                +h*so4ln2+p*rhln*so4ln2+q*so4ln3)
      rnucnum_f=convinv*exp(expon)

      return
      end
#endif
!-----------------------------------------------------------------------
      function rnucmas_f(temp,relhum,so4gas,so4frac)
!
!---->calculate total number of molecules in critical cluster [#molec/part]
!     Vehkamaeki et al. (2002)
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)

      parameter (a1=2.95413d-3,  &
     &           a2=9.76834d-2,  &
     &           a3=1.02485d-3,  &
     &           a4=2.18646d-6,  &
     &           a5=1.01717d-1)
      parameter (b1=2.05064d-3,  &
     &           b2=7.58504d-3,  &
     &           b3=1.92654d-4,  &
     &           b4=6.70430d-7,  &
     &           b5=2.55774d-1)
      parameter (c1=3.22308d-3,  &
     &           c2=8.52637d-4,  &
     &           c3=1.54757d-5,  &
     &           c4=5.66661d-8,  &
     &           c5=3.38444d-2)
      parameter (d1=4.74323d-2,  &
     &           d2=6.25104d-4,  &
     &           d3=2.65066d-6,  &
     &           d4=3.67471d-9,  &
     &           d5=2.67251d-4)
      parameter (e1=1.25211d-2,  &
     &           e2=5.80655d-3,  &
     &           e3=1.01674d-4,  &
     &           e4=2.88195d-7,  &
     &           e5=9.42243d-2)
      parameter (f1=3.85460d-2,  &
     &           f2=6.72316d-4,  &
     &           f3=2.60288d-6,  &
     &           f4=1.19416d-8,  &
     &           f5=8.51515d-3)
      parameter (g1=1.83749d-2,  &
     &           g2=1.72072d-4,  &
     &           g3=3.71766d-7,  &
     &           g4=5.14875d-10,  &
     &           g5=2.68660d-4)
      parameter (h1=6.19974d-2,  &
     &           h2=9.06958d-4,  &
     &           h3=9.11728d-7,  &
     &           h4=5.36796d-9,  &
     &           h5=7.74234d-3)
      parameter (p1=1.21827d-2,  &
     &           p2=1.06650d-4,  &
     &           p3=2.53460d-7,  &
     &           p4=3.63519d-10,  &
     &           p5=6.10065d-4)
      parameter (q1=3.20184d-4,  &
     &           q2=1.74762d-5,  &
     &           q3=6.06504d-8,  &
     &           q4=1.42177d-11,  &
     &           q5=1.35751d-4)


      parameter (so4gmin=1.d10,so4gmax=1.d17)
      parameter (convfac=1.d-6)

      sfracinv=r1/so4frac
      temp2=temp*temp
      temp3=temp*temp2

      a=-a1-a2*temp+a3*temp2-a4*temp3-a5*sfracinv
      b=-b1-b2*temp+b3*temp2-b4*temp3-b5*sfracinv
      c=+c1+c2*temp-c3*temp2+c4*temp3+c5*sfracinv
      d=+d1-d2*temp+d3*temp2-d4*temp3-d5*sfracinv
      e=-e1+e2*temp-e3*temp2+e4*temp3+e5*sfracinv
      f=-f1-f2*temp+f3*temp2+f4*temp3-f5*sfracinv
      g=-g1+g2*temp-g3*temp2-g4*temp3+g5*sfracinv
      h=-h1+h2*temp-h3*temp2-h4*temp3-h5*sfracinv
      p=+p1-p2*temp+p3*temp2-p4*temp3+p5*sfracinv
      q=+q1-q2*temp+q3*temp2-q4*temp3+q5*sfracinv


      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas))*convfac)
      so4ln2=so4ln*so4ln
      so4ln3=so4ln*so4ln2

      rnucmas_f=exp(a+b*rhln+c*rhln2+d*rhln3  &
     &             +e*so4ln+f*rhln*so4ln+g*rhln2*so4ln  &
     &             +h*so4ln2+p*rhln*so4ln2+q*so4ln3)

      return
      end
!-----------------------------------------------------------------------
      function so4dens_f(temp,so4mfrac)
!
!---->calculate density sulfate aerosol [kg/m3]
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
#ifndef PARNEW
!
!---->Jaecker-Voirol  (1988), J. Chem. Phys. (87) pp 4849-4857
!
      parameter (a10=998.94d0,  a20=982.99d0,  &
     &           a11=748.23d0,  a21=608.19d0,  &
     &           a12=4.07622d0, a22=233.26d0,  &
     &           a13=317.88d0,  a23=154.19d0)

      parameter (b10=473.52d0,  b20=250.52d0,  &
     &           b11=4903.99d0, b21=5733.14d0,  &
     &           b12=11916.5d0, b22=13138.14d0,  &
     &           b13=15057.6d0, b23=15565.78d0,  &
     &           b14=6668.37d0, b24=6618.70d0)

      parameter (t0=273.15d0)

      parameter (r60=60.d0)

      frac1=so4mfrac
      frac2=frac1*so4mfrac
      frac3=frac2*so4mfrac
      if (so4mfrac.lt.0.6d0) then
       rho1=a10+a11*frac1-a12*frac2+a13*frac3
       rho2=a20+a21*frac1+a22*frac2+a23*frac3
      else
       frac4=frac3*so4mfrac
       rho1=b10+b11*frac1-b12*frac2+b13*frac3-b14*frac4
       rho2=b20+b21*frac1-b22*frac2+b23*frac3-b24*frac4
      endif
      so4dens_f=rho1+(rho2-rho1)*(temp-t0)/r60
#else
!
!---->Vehkamaeki et al. (2002)
!
      parameter (a1= 0.7681724d0,  &
     &           a2= 2.184714d0,  &
     &           a3= 7.163002d0,  &
     &           a4=-44.31447d0,  &
     &           a5= 88.74606d0,  &
     &           a6=-75.73729d0,  &
     &           a7= 23.43228d0)
      parameter (b1= 1.808225d-3,  &
     &           b2=-9.294656d-3,  &
     &           b3=-3.742148d-2,  &
     &           b4= 2.565321d-1,  &
     &           b5=-5.362872d-1,  &
     &           b6= 4.857736d-1,  &
     &           b7=-1.629592d-1)
      parameter (c1=-3.478524d-6,  &
     &           c2= 1.335867d-5,  &
     &           c3= 5.195706d-5,  &
     &           c4=-3.717636d-4,  &
     &           c5= 7.990811d-4,  &
     &           c6=-7.458060d-4,  &
     &           c7= 2.581390d-4)
      parameter (convfac=1.d3)

      so4m2=so4mfrac*so4mfrac
      so4m3=so4mfrac*so4m2
      so4m4=so4mfrac*so4m3
      so4m5=so4mfrac*so4m4
      so4m6=so4mfrac*so4m5

      a=+a1+a2*so4mfrac+a3*so4m2+a4*so4m3  &
     &        +a5*so4m4+a6*so4m5+a7*so4m6
      b=+b1+b2*so4mfrac+b3*so4m2+b4*so4m3  &
     &        +b5*so4m4+b6*so4m5+b7*so4m6
      c=+c1+c2*so4mfrac+c3*so4m2+c4*so4m3  &
     &        +c5*so4m4+c6*so4m5+c7*so4m6
      so4dens_f=(a+b*temp+c*temp*temp)*convfac
#endif
      return
      end
!-----------------------------------------------------------------------
      function so4sig_f(temp,so4frac,so4mfrac)
!
!---->calculate surface tension of nucleus
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
#ifndef PARNEW
!
!---->Sabinina and Terbugow  (1935), Z. Phys. Chem (A173), pp 237-241
!
      parameter (a10=74.00296d0,    b10=77.04932d0,  &
     &           a11=7.68634d0,     b11=9.73321d0,  &
     &           a12=625.86132d0,   b12=59.42380d0,  &
     &           a13=5117.53488d0,  b13=5.65940d0,  &
     &           a14=10646.24244d0, b14=19.78486d0,  &
     &           a20=67.822d0,      b20=72.55489d0,  &
     &           a21=78.97377d0,    b21=32.99004d0,  &
     &           a22=207.81448d0,   b22=115.88314d0,  &
     &           a23=165.6474d0,    b23=62.03460d0,  &
     &           a24=654.16827d0)

      frac1=so4frac
      frac2=frac1*so4frac
      frac3=frac2*so4frac
      frac4=frac3*so4frac
      if (so4frac.lt.0.16d0) then
       sig1=a10+a11*frac1+a12*frac2-a13*frac3+a14*frac4
      else
       sig1=b10+b11*frac1-b12*frac2+b13*frac3+b14*frac4
      endif
      if (so4frac.lt.0.25d0) then
       sig2=a20+a21*frac1-a22*frac2-a23*frac3+a24*frac4
      else
       sig2=b20+b21*frac1-b22*frac2+b23*frac3
      endif
      so4sig_f=0.001d0*sig1+(sig2-sig1)*(temp-283.15d0)*0.000025d0
#else
!
!---->Vehkamaeki et al. (2002)
!
      parameter (a1= 0.11864d0,  &
     &           a2=-0.11651d0,  &
     &           a3= 0.76852d0,  &
     &           a4=-2.40909d0,  &
     &           a5= 2.95434d0,  &
     &           a6=-1.25852d0)
      parameter (b1=-1.5709d-4,  &
     &           b2= 4.0102d-4,  &
     &           b3=-2.3995d-3,  &
     &           b4= 7.611235d-3,  &
     &           b5=-9.37386d-3,  &
     &           b6= 3.89722d-3)

      so4m2=so4mfrac*so4mfrac
      so4m3=so4mfrac*so4m2
      so4m4=so4mfrac*so4m3
      so4m5=so4mfrac*so4m4

      a=+a1+a2*so4mfrac+a3*so4m2+a4*so4m3+a5*so4m4+a6*so4m5
      b=+b1+b2*so4mfrac+b3*so4m2+b4*so4m3+b5*so4m4+b6*so4m5
      so4sig_f=a+b*temp
#endif
      return
      end
!-----------------------------------------------------------------------
      function so4mfrac_f(temp,h2osat,h2ogas)
!
!---->calculate weight fraction h2so4 in aerosols as function
!     of temp and h2o partial pressure using parameterization of
!     Tabazadeh, Toon, Clegg, and Hamill (1997), GRL (24) pp 1931-1934
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (a11= 12.372089320d0,a12=-0.16125516114d0,  &
     &           a13=-30.490657554d0,a14=-2.1133114241d0,  &
     &           a21= 13.455394705d0,a22=-0.19213122550d0,  &
     &           a23=-34.285174607d0,a24=-1.7620073078d0)

      parameter (b11= 11.820654354d0,b12=-0.20786404244d0,  &
     &           b13=-4.8073063730d0,b14=-5.1727540348d0,  &
     &           b21= 12.891938068d0,b22=-0.23233847708d0,  &
     &           b23=-6.4261237757d0,b24=-4.9005471319d0)

      parameter (c11=-180.06541028d0 ,c12=-0.38601102592d0,  &
     &           c13=-93.317846778d0 ,c14=+273.88132245d0,  &
     &           c21=-176.95814097d0 ,c22=-0.36257048154d0,  &
     &           c23=-90.469744201d0 ,c24=+267.45509988d0)

      parameter (r70=70.d0,r190=190.d0)
      parameter (fracmin=0.002d0,fracmax=0.998d0)
      parameter (rhmin=0.01d0,rhmax=1.d0,  &
     &           rhlim1=0.05d0,rhlim2=0.85d0)
      parameter (tempmin=185.d0,tempmax=260.d0)
!
!---->water activity: aw
!
      aw=min(rhmax,max(rhmin,h2ogas/h2osat))
      if(aw.le.rhlim1) then
       y1=a11*aw**a12+a13*aw+a14
       y2=a21*aw**a22+a23*aw+a24
      else if(aw.le.rhlim2) then
       y1=b11*aw**b12+b13*aw+b14
       y2=b21*aw**b22+b23*aw+b24
      else
       y1=c11*aw**c12+c13*aw+c14
       y2=c21*aw**c22+c23*aw+c24
      end if
!
!---->sulfuric acid molality
!
      tempnew=min(tempmax,max(tempmin,temp))
      so4molal=y1+(tempnew-r190)*(y2-y1)/r70
!
!---->h2so4 weight fraction
!
      so4mfrac_f=98.d0*so4molal/(98.d0*so4molal+1000.d0)
      so4mfrac_f=max(fracmin,min(fracmax,so4mfrac_f))

      return
      end
!-----------------------------------------------------------------------
      function pathair_f(temp,press)
!
!---->calculate mean free path for air
!     following eq. 8.5 Seinfeld (1997)
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (dair=4.7149676d-10,sqr2=1.4142136d0)

      pathair_f=boltz*temp/(sqr2*pi*dair*dair*press)

      return
      end
!-----------------------------------------------------------------------
      function pathso4_f(temp,press)
!
!---->calculate mean free path for sulfate molecules in air
!     following Davis (1983) - eq. 8.11 Seinfeld (1997)
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (dso4=5.1209323d-10)

      pathso4_f=boltz*temp/(sqrt(r1+so4mol/drymol)*pi*dso4*dso4*press)

      return
      end
!-----------------------------------------------------------------------
      function diffso4_f(temp,press)
!
!---->calculate diffusion coefficient for h2so4 gas molec in air
!     following Chapman and Cowling (1970), Davis (1983)
!     eq. 8.13 Seinfeld (1997) / eq. 17.17 Jacobson (1998)
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (r3=3.d0,r8=8.d0)
      parameter (dso4=5.1209323d-10)

      diffso4_f=r3*boltz*temp/(r8*dso4*dso4*press)  &
     &         *sqrt(rgas*temp*(drymol+so4mol)/(r2*pi*drymol*so4mol))

      return
      end
!-----------------------------------------------------------------------
      function viscos_f(temp)
!
!---->calculate dynamic viscosity of air
!     Sutherland's equations (List, 1984) in Jacobson (1999) eq. 4.55
!
!
!---->load modules
!
      use precision
      use modpar_umaer
!
!---->define accuracy and local parameter
!
      implicit real(kdef) (a-h,o-z), integer (i-n)
      parameter (r3h=1.5d0)
      parameter (a1=1.8325d-5,  &
     &           a2=416.16d0,  &
     &           a3=120.d0,  &
     &           a4=296.16d0)

      viscos_f=a1*a2/(temp+a3)*(temp/a4)**r3h

      return
      end
