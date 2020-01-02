      subroutine pscdisc(npsum,tmax)
      use pscalloc
      implicit none
c
c     Last modified: Potsdam, Feb, 2019, by R. Wang
c
      integer*4 npsum
      real*8 tmax
c
c     returned outputs:
c     outputs through common blocks
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 i,j,is,izs,ix,iy,nx,ny,nxy,ips,iptch,npsmax,ierr
      real*8 dr,x,y,x0,y0,dx,dy,dzs,st,di,ra
      real*8 step_s,step_d,dmslp,dmopn_cl,dmopn_ep
      real*8 pz0,bga,bgc,sma,smb,smc,dmw
      real*8 sss,sss2,ss2s,ssd,ssd2,ss2d,ssr,ssr2,ss2r
      real*8 css,css2,cs2s,csd,csd2,cs2d,csr,csr2,cs2r
      real*8 sm(3,3),om(3,3)
c
      allocate(nps(nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: nps not allocated!'
      allocate(mtenexp(3,3,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: mtenexp not allocated!'
      allocate(mtenshr(3,3,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: mtenshr not allocated!'
      allocate(mscaexp(nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: mscaexp not allocated!'
      allocate(mscashr(nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: mscashr not allocated!'
c
      write(*,'(a)')' ... discretise rectangular fault segments:'
c
      do izs=1,nzs
        nps(izs)=0
      enddo
c
      dr=2.d0*(r2-r1)/dble(nr-1)/(1.d0+sampratio)
      if(nzs.gt.1)then
        dzs=(zs2-zs1)/dble(nzs-1)
      else
        dzs=dr
      endif
      dr=dmax1(dr,dzs)
c
c     determine npsmax
c
      do is=1,ns
        if(tstart(is).gt.tmax)goto 100
        step_s=length(is)/dble(nptch_s(is))
        step_d=width(is)/dble(nptch_d(is))
        do iptch=1,nptch_s(is)*nptch_d(is)
          dx=dr
          nx=max0(1,idnint(step_s/dx))
          dx=step_s/dble(nx)
c
          if(dabs(dzs*ssd).gt.0.d0)then
            dy=dmin1(dr,dzs/ssd)
          else
            dy=dr
          endif
          ny=max0(1,idnint(step_d/dy))
          dy=step_d/dble(ny)
c
          x0=ptch_s(is,iptch)-0.5d0*step_s
          y0=ptch_d(is,iptch)-0.5d0*step_d
c
          do iy=1,ny
            y=y0+dy*(dble(iy)-0.5d0)
            pz0=zref(is)+y*ssd
            if(dzs.gt.0.d0)then
              izs=idnint((pz0-zs1)/dzs)+1
              izs=max0(min0(izs,nzs),1)
            else
              izs=1
            endif
            do ix=1,nx
              x=x0+dx*(dble(ix)-0.5d0)
              nps(izs)=nps(izs)+1
            enddo
          enddo
        enddo
100     continue
      enddo
c
      npsmax=0
      do izs=1,nzs
        npsmax=max0(npsmax,nps(izs))
      enddo
c
      allocate(isno(npsmax,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: isno not allocated!'
      allocate(plat(npsmax,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: plat not allocated!'
      allocate(plon(npsmax,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: plon not allocated!'
      allocate(pz(npsmax,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: pz not allocated!'
      allocate(pmwei(6,npsmax,nzs),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscdisc: pmwei not allocated!'
c
      npsum=0
c
      do izs=1,nzs
        do i=1,3
          do j=1,3
            mtenexp(i,j,izs)=0.d0
            mtenshr(i,j,izs)=0.d0
          enddo
        enddo
        mscaexp(izs)=0.d0
        mscashr(izs)=0.d0
      enddo
c
      do izs=1,nzs
        nps(izs)=0
      enddo
c
      do is=1,ns
        if(tstart(is).gt.tmax)goto 200
        nxy=0
        st=strike(is)*DEG2RAD
        di=dip(is)*DEG2RAD
c
        sss=dsin(st)
        sss2=sss*sss
        ss2s=dsin(2.d0*st)
        ssd=dsin(di)
        ssd2=ssd*ssd
        ss2d=dsin(2.d0*di)
c
        css=dcos(st)
        css2=css*css
        cs2s=dcos(2.d0*st)
        csd=dcos(di)
        csd2=csd*csd
        cs2d=dcos(2.d0*di)
c
        step_s=length(is)/dble(nptch_s(is))
        step_d=width(is)/dble(nptch_d(is))
c
        do iptch=1,nptch_s(is)*nptch_d(is)
          ra=datan2(-slip_d(is,iptch),slip_s(is,iptch))
c
          ssr=dsin(ra)
          ssr2=ssr*ssr
          ss2r=dsin(2.d0*ra)
c
          csr=dcos(ra)
          csr2=csr*csr
          cs2r=dcos(2.d0*ra)
c
          sm(1,1)=-ssd*csr*ss2s-ss2d*ssr*sss2
          sm(2,2)= ssd*csr*ss2s-ss2d*ssr*css2
          sm(3,3)=-(sm(1,1)+sm(2,2))
          sm(1,2)= ssd*csr*cs2s+0.5d0*ss2d*ssr*ss2s
          sm(2,1)=sm(1,2)
          sm(2,3)=-csd*csr*sss+cs2d*ssr*css
          sm(3,2)=sm(2,3)
          sm(3,1)=-csd*csr*css-cs2d*ssr*sss
          sm(1,3)=sm(3,1)
c
c         openning => explosion (ep) + clvd (om)
c
          om(1,1)=-0.5d0+1.5d0*sss2*ssd2
          om(2,2)=-0.5d0+1.5d0*css2*ssd2
          om(3,3)=-(om(1,1)+om(2,2))
          om(1,2)=-1.5d0*sss*css*ssd2
          om(2,1)=om(1,2)
          om(2,3)=-1.5d0*css*ssd*csd
          om(3,2)=om(2,3)
          om(3,1)= 1.5d0*sss*ssd*csd
          om(1,3)=om(3,1)
c
          dx=dr
          nx=max0(1,idnint(step_s/dx))
          dx=step_s/dble(nx)
c
          if(dabs(dzs*ssd).gt.0.d0)then
            dy=dmin1(dr,dzs/ssd)
          else
            dy=dr
          endif
          ny=max0(1,idnint(step_d/dy))
          dy=step_d/dble(ny)
c
          dmslp=dsqrt(slip_s(is,iptch)**2+slip_d(is,iptch)**2)
     &         *step_s*step_d/dble(nx*ny)
          dmopn_ep=opening(is,iptch)
     &         *step_s*step_d/dble(nx*ny)
          dmopn_cl=dmopn_ep*4.d0/3.d0
c
          x0=ptch_s(is,iptch)-0.5d0*step_s
          y0=ptch_d(is,iptch)-0.5d0*step_d
c
          do iy=1,ny
            y=y0+dy*(dble(iy)-0.5d0)
            pz0=zref(is)+y*ssd
            if(pz0.lt.zs1-dzs.or.pz0.gt.zs2+dzs)then
              nwarn=nwarn+1
              print *,'Warning in pscdisc: the fault plane in the'
              print *,' depth range not covered by Green funcions!'
            endif
            if(dzs.gt.0.d0)then
              izs=idnint((pz0-zs1)/dzs)+1
              izs=max0(min0(izs,nzs),1)
            else
              izs=1
            endif
c
            do i=1,3
              mtenexp(i,i,izs)=mtenexp(i,i,izs)+dmopn_ep*dble(nx)
              do j=1,3
                mtenshr(i,j,izs)=mtenshr(i,j,izs)+dmslp*dble(nx)*sm(i,j)
     &                  +dmopn_cl*dble(nx)*om(i,j)
              enddo
            enddo
c
            mscaexp(izs)=mscaexp(izs)
     &                  +dabs(dmopn_ep)*dble(nx)*dsqrt(1.5d0)
            dmw=0.d0
            do i=1,3
              dmw=dmw+0.5d0*(dmslp*dble(nx)*sm(i,i)
     &                      +dmopn_cl*dble(nx)*om(i,i))**2
              do j=i+1,3
                dmw=dmw+(dmslp*dble(nx)*sm(i,j)
     &                  +dmopn_cl*dble(nx)*om(i,j))**2
              enddo
            enddo
            mscashr(izs)=mscashr(izs)+dsqrt(dmw)
c
            do ix=1,nx
              x=x0+dx*(dble(ix)-0.5d0)
              nps(izs)=nps(izs)+1
              ips=nps(izs)
c
c             spherical triangle:
c             A = pole, B = source position, C = reference position
c
              sma=dsqrt(x**2+(y*dcos(di))**2)/REARTH
              smb=0.5d0*PI-latref(is)*DEG2RAD
              bgc=st+datan2(y*dcos(di),x)
              smc=dacos(dcos(sma)*dcos(smb)
     &           +dsin(sma)*dsin(smb)*dcos(bgc))
              bga=dasin(dsin(sma)*dsin(bgc)/dsin(smc))
c
c             geographic coordinate of the point source
c
              plat(ips,izs)=90.d0-smc/DEG2RAD
              plon(ips,izs)=dmod(lonref(is)+bga/DEG2RAD,360.d0)
              pz(ips,izs)=pz0
              isno(ips,izs)=is
c
c             1 = weight for strike-slip: m12=m21=1;
c             2 = weight for dip-slip: m13=m31=1
c             3 = weight for clvd: m33=-(m11+m22)=1
c             4 = weight for 45 deg strike-slip: m11=-m22=1
c             5 = weight for 45 deg dip-slip: m23=m32=1
c             6 = weight for explosion: m11=m22=m33=1
c
              pmwei(1,ips,izs)=sm(1,2)*dmslp+om(1,2)*dmopn_cl
              pmwei(2,ips,izs)=sm(3,1)*dmslp+om(3,1)*dmopn_cl
              pmwei(3,ips,izs)=sm(3,3)*dmslp+om(3,3)*dmopn_cl
              pmwei(4,ips,izs)=0.5d0*(sm(1,1)-sm(2,2))*dmslp
     &                        +0.5d0*(om(1,1)-om(2,2))*dmopn_cl
              pmwei(5,ips,izs)=sm(2,3)*dmslp+om(2,3)*dmopn_cl
c
              pmwei(6,ips,izs)=dmopn_ep
            enddo
          enddo
          nxy=nxy+nx*ny
        enddo
        if(nxy.gt.1)then
          write(*,'(a,2(i4,a))')'  ',is,'. segment re-discretised to ',
     &                nxy,' patches.'
        endif
        npsum=npsum+nxy
200     continue
      enddo
      write(*,*)'------------------------------------------------'
     &        //'----------------'
      print *,'    total number of fault patches: ',npsum
      write(*,*)'------------------------------------------------'
     &        //'----------------'
      return
      end
