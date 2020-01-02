      subroutine psgspec(isp,nr1,nr2)
      use psgalloc
      implicit none
c
      integer*4 isp,nr1,nr2
c
      integer*4 i,ir,istp,l,lf,jf,istate,it,its,itga,itgb
      integer*4 ntg,nfg,nkc,ierr
      real*8 f,dfg,t,dtg,fc,fcut,rc,kc,dk0,dk
      real*8 a,b,alpha,beta,thick,dkmin
      real*8 fgrnabs,fgrnmax
      real*8 dspabs(4),stsabs(4),stnabs(4),potabs(4),grvabs(4)
      real*8 dspmax(4),stsmax(4),stnmax(4),potmax(4),grvmax(4)
      complex*16 cs
      complex*16 dum(14,4),dum0(14,4)
      logical*2 again,dfgfound,lowpass,tty
c
      integer*4 nkmax
      real*8 eps,pi2
      data nkmax/100000/
      data eps,pi2/0.05d0,6.28318530717959d0/

      print *,' ==============================================='
      write(*,'(i3,a,F12.4,a,f12.4,a)')isp,'. sub-profile: ',
     &        r(nr1)/km2m,' -> ',r(nr2)/km2m,' km'
c
c     determine wavenumber sampling rate
c
      rc=rs(nr1)+0.5d0*(r(nr1)+r(nr2))
      kc=pi2/rc
      dk0=eps*pi2/dsqrt((zrec-zs)**2+rc**2)
      thick=0.d0
      do l=1,lp-1
        thick=thick+hp(l)
      enddo
      dkmin=eps*accuracy*pi2/dsqrt(thick**2+rc**2)
c
c     determine frequency sampling rate
c
      dfg=1.d0/(2.d0*dble(nt)*dt)
      alpha=pi2*dfg
c
c     istate = -1: elastic case
c
      istate=-1
      cs=(0.d0,0.d0)
      dk=dk0
      do istp=1,4
        do i=1,14
          do ir=nr1,nr2
            du(0,ir,i,istp)=(0.d0,0.d0)
          enddo
        enddo
      enddo
100   nkc=idint(2.d0*kc/dk)
      call psgwvint(istate,1,nr1,nr2,cs,nkc,nkmax,dk,.false.)
      again=.false.
      do istp=1,4
        dspabs(istp)=0.d0
        stsabs(istp)=0.d0
        stnabs(istp)=0.d0
        potabs(istp)=0.d0
        grvabs(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspabs(istp)=dspabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=4,9
            stsabs(istp)=stsabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=10,12
            stnabs(istp)=stnabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          potabs(istp)=potabs(istp)
     &                +geow(ir)*cdabs(du(1,ir,13,istp)-du(0,ir,13,istp))
          grvabs(istp)=grvabs(istp)
     &                +geow(ir)*cdabs(du(1,ir,14,istp)-du(0,ir,14,istp))
        enddo
        dspmax(istp)=0.d0
        stsmax(istp)=0.d0
        stnmax(istp)=0.d0
        potmax(istp)=0.d0
        grvmax(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspmax(istp)=dspmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=4,9
            stsmax(istp)=stsmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=10,12
            stnmax(istp)=stnmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          potmax(istp)=potmax(istp)+geow(ir)*cdabs(du(1,ir,13,istp))
          grvmax(istp)=grvmax(istp)+geow(ir)*cdabs(du(1,ir,14,istp))
        enddo
        again=again
     &        .or.dspabs(istp).gt.accuracy*dspmax(istp)
     &        .or.stsabs(istp).gt.accuracy*stsmax(istp)
     &        .or.stnabs(istp).gt.accuracy*stnmax(istp)
     &        .or.potabs(istp).gt.accuracy*potmax(istp)
     &        .or.grvabs(istp).gt.accuracy*grvmax(istp)
      enddo
      if(again.and.dk.gt.dkmin)then
        do istp=1,4
          do i=1,14
            do ir=nr1,nr2
              du(0,ir,i,istp)=du(1,ir,i,istp)
            enddo
          enddo
        enddo
        dk=0.5d0*dk
        goto 100
      endif
      if(again)then
        print *,' Warning in psgspec: required acc. not satisfied!'
        nwarn=nwarn+1
      endif
c
      istate=1
      cs=dcmplx(alpha,dfg)
      dk=4.d0*dk
      do istp=1,4
        do i=1,14
          do ir=nr1,nr2
            du(0,ir,i,istp)=(0.d0,0.d0)
          enddo
        enddo
      enddo
200   nkc=idint(2.d0*kc/dk)
      call psgwvint(istate,1,nr1,nr2,cs,nkc,nkmax,dk,.false.)
      again=.false.
      do istp=1,4
        dspabs(istp)=0.d0
        stsabs(istp)=0.d0
        stnabs(istp)=0.d0
        potabs(istp)=0.d0
        grvabs(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspabs(istp)=dspabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=4,9
            stsabs(istp)=stsabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=10,12
            stnabs(istp)=stnabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          potabs(istp)=potabs(istp)
     &                +geow(ir)*cdabs(du(1,ir,13,istp)-du(0,ir,13,istp))
          grvabs(istp)=grvabs(istp)
     &                +geow(ir)*cdabs(du(1,ir,14,istp)-du(0,ir,14,istp))
        enddo
        dspmax(istp)=0.d0
        stsmax(istp)=0.d0
        stnmax(istp)=0.d0
        potmax(istp)=0.d0
        grvmax(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspmax(istp)=dspmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=4,9
            stsmax(istp)=stsmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=10,12
            stnmax(istp)=stnmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          potmax(istp)=potmax(istp)+geow(ir)*cdabs(du(1,ir,13,istp))
          grvmax(istp)=grvmax(istp)+geow(ir)*cdabs(du(1,ir,14,istp))
        enddo
        again=again
     &        .or.dspabs(istp).gt.accuracy*dspmax(istp)
     &        .or.stsabs(istp).gt.accuracy*stsmax(istp)
     &        .or.stnabs(istp).gt.accuracy*stnmax(istp)
     &        .or.potabs(istp).gt.accuracy*potmax(istp)
     &        .or.grvabs(istp).gt.accuracy*grvmax(istp)
      enddo
      if(again.and.dk.gt.dkmin)then
        do istp=1,4
          do i=1,14
            do ir=nr1,nr2
              du(0,ir,i,istp)=du(1,ir,i,istp)
            enddo
          enddo
        enddo
        dk=0.5d0*dk
        goto 200
      endif
      if(again)then
        print *,' Warning in psgspec: required acc. not satisfied!'
        nwarn=nwarn+1
      endif
c
      istate=-1
      cs=(0.d0,0.d0)
      call psgwvint(istate,-1,nr1,nr2,cs,nkc,nkmax,dk,.true.)
c
      if(nt.eq.1)then
        do istp=1,4
          do i=1,14
            do ir=nr1,nr2
              du(1,ir,i,istp)=dcmplx(dreal(du(-1,ir,i,istp)),0.d0)
            enddo
          enddo
          do ir=nr1,nr2
            du(1,ir,14,istp)=du(1,ir,14,istp)
     &                +dcmplx(2.d0*g0/rearth,0.d0)*du(1,ir,1,istp)
          enddo
        enddo
        return
      endif
c
      if(nt.eq.2)then
        istate=0
        cs=(0.d0,0.d0)
        call psgwvint(istate,0,nr1,nr2,cs,nkc,nkmax,dk,.true.)
        do istp=1,4
          do i=1,14
            do ir=nr1,nr2
              du(1,ir,i,istp)=dcmplx(dreal(du(-1,ir,i,istp)),
     &                dreal(du(0,ir,i,istp)-du(-1,ir,i,istp)))
            enddo
          enddo
          do ir=nr1,nr2
            du(1,ir,14,istp)=du(1,ir,14,istp)
     &                +dcmplx(2.d0*g0/rearth,0.d0)*du(1,ir,1,istp)
          enddo
        enddo
        return
      endif
c
      print *,' Finding upper cut-off frequency'
      do istp=1,4
        dspmax(istp)=0.d0
        stsmax(istp)=0.d0
        stnmax(istp)=0.d0
        potmax(istp)=0.d0
        grvmax(istp)=0.d0
      enddo
      istate=1
      lf=1
      nfg=1
300   fcut=dble(nfg-1)*dfg
      cs=dcmplx(alpha,pi2*fcut)
      call psgwvint(istate,lf,nr1,nr2,cs,nkc,nkmax,dk,.true.)
      do istp=1,4
        do i=1,14
          do ir=nr1,nr2
            du(lf,ir,i,istp)=du(lf,ir,i,istp)-du(-1,ir,i,istp)
          enddo
        enddo
      enddo
      again=.false.
      do istp=1,4
        dspabs(istp)=0.d0
        stsabs(istp)=0.d0
        stnabs(istp)=0.d0
        potabs(istp)=0.d0
        grvabs(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspabs(istp)=dspabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
          do i=4,9
            stsabs(istp)=stsabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
          do i=10,12
            stnabs(istp)=stnabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
          potabs(istp)=potabs(istp)+geow(ir)*cdabs(du(lf,ir,13,istp))
          grvabs(istp)=grvabs(istp)+geow(ir)*cdabs(du(lf,ir,14,istp))
        enddo
        dspmax(istp)=dmax1(dspmax(istp),dspabs(istp))
        stsmax(istp)=dmax1(stsmax(istp),stsabs(istp))
        stnmax(istp)=dmax1(stnmax(istp),stnabs(istp))
        potmax(istp)=dmax1(potmax(istp),potabs(istp))
        grvmax(istp)=dmax1(grvmax(istp),grvabs(istp))
        again=again
     &        .or.dspabs(istp).gt.eps*dspmax(istp)
     &        .or.stsabs(istp).gt.eps*stsmax(istp)
     &        .or.stnabs(istp).gt.eps*stnmax(istp)
     &        .or.potabs(istp).gt.eps*potmax(istp)
     &        .or.grvabs(istp).gt.eps*grvmax(istp)
      enddo
      if(again.and.nfg.le.nfmax.or.nfg.lt.max0(nfmin,nt/2))then
        lf=lf+1
        nfg=2*nfg
        goto 300
      endif
c
c     cut-off frequency found
c
      if(nfg.gt.nfmax)then
        nfg=nfmax
        lowpass=.true.
      else if(nfg.lt.nfmin)then
        nfg=nfmin
        lowpass=.false.
      else
        lowpass=.false.
      endif
      fcut=dble(nfg-1)*dfg
      write(*,'(i6,a,E12.5,a)')nfg,' frequencies will be used,'
     &                             //' cut-off: ',fcut,'Hz'
      print *,' ==============================================='
c
      istate=1
      jf=1+nfg/10
      do lf=1,nfg
        f=dble(lf-1)*dfg
        cs=dcmplx(alpha,pi2*f)
        tty=lf.eq.1.or.lf.eq.nfg.or.mod(lf,jf).eq.1
        call psgwvint(istate,lf,nr1,nr2,cs,nkc,nkmax,dk,tty)
      enddo
c
c     convert to time domain by FFT
c
      dtg=1.d0/(dble(2*nfg)*dfg)
      ntg=2*nfg
c
      do istp=1,4
        do i=1,14
          if(select(i,istp))then
            do ir=nr1,nr2
              do lf=1,nfg
                fgrn(lf)=du(lf,ir,i,istp)-du(-1,ir,i,istp)
              enddo
              if(lowpass)then
c
c               low-pass filter has to be used
c
c               find a proper corner frequency as high as possible
c
                fc=fcut/eps
600             fgrnmax=0.d0
                beta=1.d0+alpha/(pi2*fc)
                do lf=1,nfg
                  fgrnabs=cdabs(fgrn(lf))/(beta**2
     &                   +(dble(lf-1)*dfg/fc)**2)**1.5d0
                  fgrnmax=dmax1(fgrnmax,fgrnabs)
                enddo
                if(fgrnabs.gt.eps*fgrnmax)then
                  fc=0.75d0*fc
                  goto 600
                endif
c
c               use low-pass filter with the corner frequency found
c
                do lf=1,nfg
                  fgrn(lf)=fgrn(lf)/(dcmplx(beta,dble(lf-1)*dfg/fc))**3
                enddo
              endif
	        jf=1
	        do lf=2*nfg,nfg+2,-1
	          jf=jf+1
	          fgrn(lf)=dconjg(fgrn(jf))
	        enddo
	        fgrn(nfg+1)=(0.d0,0.d0)
c
c	        convention for Fourier transform:
c	        f(t)=\int F(f) exp(i2\pi f t) df
c
	        call four1w(fgrn,dswap,2*nfg,+1)
c
c             calculate response to the Heaviside source history
c
              tgrn(1)=0.d0
              do it=2,ntg
                tgrn(it)=tgrn(it-1)+dreal(fgrn(it-1))
     &                  *dtg*dfg*dexp(alpha*dble(it-1)*dtg)
              enddo
c
c             linear interpolation
c
c             t=0: co-seismic
c             t>0: postseismic (coseismic subtracted)
c
              dobs(1)=dreal(du(-1,ir,i,istp))
              do it=2,2*((1+nt)/2)
                t=dble(it-1)*dt
                b=dmod(t/dtg,1.d0)
                a=1.d0-b
                itga=min0(1+idint(t/dtg),ntg)
                itgb=min0(2+idint(t/dtg),ntg)
                dobs(it)=tgrn(itga)*a+tgrn(itgb)*b
              enddo
c
              do it=1,(1+nt)/2
                du(it,ir,i,istp)=dcmplx(dobs(2*it-1),dobs(2*it))
              enddo
            enddo
          endif
        enddo
c
c       gravity change measured by gravimeters fixed on the surface
c
        do ir=nr1,nr2
          do it=1,(1+nt)/2
            du(it,ir,14,istp)=du(it,ir,14,istp)
     &                +dcmplx(2.d0*g0/rearth,0.d0)*du(it,ir,1,istp)
          enddo
        enddo
      enddo
c
      return
      end