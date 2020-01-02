      subroutine psgprocess(ierr)
      use psgalloc
      implicit none
c
      integer*4 ierr
c
c     work space
c
      integer*4 i,l,ir,it,izs,nls,istp,isp
      integer*4 nr1,nr2,nlr,nprf,leninp,iunit
      real*8 am,rsmin,dratio,swap
c
      nfmin=64
      if(taumin.le.0.d0)then
        nfmax=nfmin
      else
        nfmax=nfmin
150     nfmax=2*nfmax
        if(dble(2+nfmax-1)*0.1d0*taumin.lt.twindow)goto 150
      endif
c
      allocate(tgrn(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: tgrn not allocated!'
      allocate(fgrn(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: fgrn not allocated!'
      allocate(dobs(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: obs not allocated!'
      allocate(dswap(4*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: dswap not allocated!'
c
      
      allocate(rs(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: rs not allocated!'
      allocate(geow(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: geow not allocated!'
      allocate(obs(nr,16,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: obs not allocated!'
      allocate(obs0(nr,16,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: obs0 not allocated!'
c
      allocate(idec(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: idec not allocated!'
      allocate(nout(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: nout not allocated!'
      allocate(du(-1:nfmax,nr,14,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: du not allocated!'
c
      lp=n0+2
      allocate(hp(lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: hp not allocated!'
      allocate(nno(lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: nno not allocated!'
      allocate(zp(lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: zp not allocated!'
c
      allocate(hkup(2,2,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: hkup not allocated!'
      allocate(hklw(2,2,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: hklw not allocated!'
      allocate(maup(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: maup not allocated!'
      allocate(maiup(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: maiup not allocated!'
      allocate(malw(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: malw not allocated!'
      allocate(mailw(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: mailw not allocated!'
c
      zs=0.d0
      call psglayer(ierr)
      nlr=nno(lzrec)
c
      leninp=index(inputfile,' ')-1
c
      stype(1)='explosion (M11=M22=M33=1*kappa)'
      stype(2)='strike-slip (M12=M21=1*mue)'
      stype(3)='dip-slip (M13=M31=1*mue)'
      stype(4)='clvd (M33=1*mue, M11=M22=-M33/2)'
c
      comptxt(1)='Uz (vertical displacement)'
      comptxt(2)='Ur (radial displacement)'
      comptxt(3)='Ut (tangential displacement)'
      comptxt(4)='Szz (linear stress)'
      comptxt(5)='Srr (linear stress)'
      comptxt(6)='Stt (linear stress)'
      comptxt(7)='Szr (shear stress)'
      comptxt(8)='Srt (shear stress)'
      comptxt(9)='Stz (shear stress)'
      comptxt(10)='Tr (tilt -dUr/dz)'
      comptxt(11)='Tt (tilt -dUt/dz)'
      comptxt(12)='Rot (rotation ar. z-axis)'
      comptxt(13)='Gd (geoid changes)'
      comptxt(14)='Gr (gravity changes)'
c
      iunit=10
      do istp=1,4
        do i=1,14
          if(select(i,istp))then
            iunit=iunit+1
            unit(i,istp)=iunit
            open(unit(i,istp),file=green(i,istp),status='unknown')
            write(unit(i,istp),'(a)')'################################'
            write(unit(i,istp),'(a)')'# The input file used: '
     &                        //inputfile(1:leninp)
            write(unit(i,istp),'(a)')'################################'
            write(unit(i,istp),'(a)')'# Greens function component: '
     &                        //comptxt(i)
            write(unit(i,istp),'(a)')'#(Okada solutions subtracted)'
            write(unit(i,istp),'(a)')'# Source type: '//stype(istp)
            write(unit(i,istp),'(a)')'# Observation distance sampling:'
            write(unit(i,istp),'(a)')'#    nr        r1[m]        r2[m]'
     &                             //'  samp_ratio'
            write(unit(i,istp),'(i7,2E14.6,f10.4)')nr,r1,r2,sampratio
            write(unit(i,istp),'(a)')'# Uniform obs. site parameters:'
            write(unit(i,istp),'(a)')'#    depth[m]       la[Pa]       '
     &       //'mu[Pa]  rho[kg/m^3]    etk[Pa*s]    etm[Pa*s]     alpha'
            write(unit(i,istp),'(7E13.6)')zrec,la(nlr),
     &           mu(nlr),rho(nlr),etk(nlr),etm(nlr),alf(nlr)
            write(unit(i,istp),'(a)')'# Source depth sampling:'
            write(unit(i,istp),'(a)')'#   nzs       zs1[m]       zs2[m]'
            write(unit(i,istp),'(i7,2d14.6)')nzs,zs1,zs2
            write(unit(i,istp),'(a)')'# Time sampling:'
            write(unit(i,istp),'(a)')'#    nt        t-window[s]'
            write(unit(i,istp),'(i7,E24.16)')nt,twindow
            write(unit(i,istp),'(a)')'# Data in each source depth block'
            write(unit(i,istp),'(a)')'# ==============================='
            write(unit(i,istp),'(a)')'# Line 1: source layer parameters'
            write(unit(i,istp),'(a)')'#  s_depth, la, mu, rho,'
     &                              //' etk, etm, alpha'
            write(unit(i,istp),'(a)')'# Line 2: coseismic responses '
     &                             //'(f(ir,it=1),ir=1,nr)'
            write(unit(i,istp),'(a)')'# Line 3: (idec(ir),ir=1,nr)'
     &                //'(decimal exponents for postseismic responses)'
            write(unit(i,istp),'(a)')'# Line 4: (f(ir,it=2),ir=1,nr)'
            write(unit(i,istp),'(a)')'# Line 5: (f(ir,it=3),ir=1,nr)'
            write(unit(i,istp),'(a)')'#  ...'
          endif
        enddo
      enddo
c
      call psgbsj(ierr)
c
      do izs=1,nzs
        zs=zs1+dble(izs-1)*dzs
        write(*,'(/,a,E13.4,a)')' Processing for the '
     &                 //'source at depth:',zs,' m.'
c
        call psglayer(ierr)
c
        do l=1,lp
          zp(l)=0.d0
          do i=1,l-1
            if(nno(i).eq.nno(l))zp(l)=zp(l)+hp(i)
          enddo
        enddo
        nls=nno(ls)
c
        do istp=1,4
          do i=1,14
            do ir=1,nr
              do it=1,nfmax
                du(it,ir,i,istp)=(0.d0,0.d0)
              enddo
            enddo
          enddo
        enddo
c
        zrs2=(zrec-zs)**2
        rsmin=0.5d0*dr
        do ir=1,nr
          rs(ir)=dmax1(rsmin,0.1d0*dsqrt(zrs2+r(ir)**2))
          geow(ir)=zrs2+(rs(ir)+r(ir))**2
        enddo
c
        swap=dsqrt(zrs2+(rs(nr)+r(nr))**2)/dsqrt(zrs2+(rs(1)+r(1))**2)
        nprf=1+idnint(dlog(swap)/dlog(2.5d0))
        if(nprf.gt.1)then
          dratio=swap**(1.d0/dble(nprf-1))
        else
          dratio=2.5d0
        endif
c
        isp=0
        nr2=0
200     isp=isp+1
        nr1=nr2+1
        nr2=nr1
        do ir=nr1+1,nr
          if(r(ir).le.dratio*dsqrt(zrs2+(rs(nr1)+r(nr1))**2))nr2=ir
        enddo
        call psgspec(isp,nr1,nr2)
        if(nr2.lt.nr)goto 200
c
        do istp=1,4
          do i=1,14
            if(.not.select(i,istp))goto 400
            write(unit(i,istp),'(a)')'#################################'
            write(unit(i,istp),'(a,i2,a)')'# the ',izs,'. source depth:'
            write(unit(i,istp),'(a)')'#################################'
            write(unit(i,istp),'(7E13.6)')zs,la(nls),
     &        mu(nls),rho(nls),etk(nls),etm(nls),alf(nls)
            do ir=1,nr-1
              write(unit(i,istp),'(E14.6,$)')dreal(du(1,ir,i,istp))
            enddo
            write(unit(i,istp),'(E14.6)')dreal(du(1,nr,i,istp))
            do ir=1,nr
              du(1,ir,i,istp)=dcmplx(0.d0,dimag(du(1,ir,i,istp)))
              am=0.d0
              do it=1,(nt+1)/2
                am=dmax1(am,dabs(dreal(du(it,ir,i,istp))),
     &                      dabs(dimag(du(it,ir,i,istp))))
              enddo
              if(am.le.0.d0)then
                idec(ir)=0
              else
                idec(ir)=idint(dlog10(am))-4
                do it=1,(nt+1)/2
                  du(it,ir,i,istp)=du(it,ir,i,istp)
     &                  *dcmplx(10.d0**(-idec(ir)),0.d0)
                enddo
              endif
            enddo
            call outint(unit(i,istp),idec,nr)
            do it=1,(nt+1)/2
              if(it.gt.1)then
                do ir=1,nr
                  nout(ir)=idnint(dreal(du(it,ir,i,istp)))
                enddo
                call outint(unit(i,istp),nout,nr)
              endif
              if(it*2.le.nt)then
                do ir=1,nr
                  nout(ir)=idnint(dimag(du(it,ir,i,istp)))
                enddo
                call outint(unit(i,istp),nout,nr)
              endif
            enddo
400         continue
          enddo
        enddo
      enddo
c
c     end of izs loop
c
      do istp=1,4
        do i=1,14
          if(select(i,istp))close(unit(i,istp))
        enddo
      enddo
c
      return
      end
