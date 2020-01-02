      subroutine psggetinp(ierr)
      use psgalloc
      implicit none
c
      integer*4 ierr
c
c     work space
c
      integer*4 i,j,l,ir,istp,lend,lenf
      real*8 dract
      real*8 kg,swap,vp,vs
c
c     read input file file
c
      open(10,file=inputfile,status='old')
c
c     parameters for source-observation array
c     =======================================
c
      call skipdoc(10)
      read(10,*)zrec,ioc
      zrec=zrec*km2m
      call skipdoc(10)
      read(10,*)nr,r1,r2,sampratio
c
      allocate(r(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: r not allocated!'
c
      if(sampratio.lt.1.d0)then
        stop 'Error: max. to min. sampling ratio < 1!'
      endif
      r1=r1*km2m
      r2=r2*km2m
      if(r1.gt.r2)then
        swap=r1
        r1=r2
        r2=swap
      endif
      if(r1.lt.0.d0.or.r2.lt.0.d0.or.nr.lt.1)then
        stop 'Error: wrong no of distance samples!'
      else if(nr.eq.1.or.r1.eq.r2)then
        r2=r1
        nr=1
        dr=0.d0
        r(1)=r1
      else if(nr.eq.2)then
        dr=r2-r1
        r(1)=r1
        r(2)=r2
      else
        dr=2.d0*(r2-r1)/dble(nr-1)/(1.d0+sampratio)
        r(1)=r1
        do i=2,nr
          dract=dr*(1.d0+(sampratio-1.d0)*dble(i-2)/dble(nr-2))
          r(i)=r(i-1)+dract
        enddo
      endif
c
      call skipdoc(10)
      read(10,*)nzs,zs1,zs2
      if(zs1.gt.zs2)then
        swap=zs1
        zs1=zs2
        zs2=swap
      endif
      zs1=zs1*km2m
      zs2=zs2*km2m
      if(zs1.lt.0.d0.or.zs2.lt.0.d0.or.nzs.lt.1)then
        stop 'Error: wrong no of source depths!'
      else if(nzs.eq.1.or.zs1.eq.zs2)then
        nzs=1
        dzs=0.d0
      else
        dzs=(zs2-zs1)/dble(nzs)
        zs1=zs1+0.5d0*dzs
        zs2=zs1+dble(nzs-1)*dzs
        if(zrec.ge.zs1.and.zrec.le.zs2)then
          swap=dmod(zrec,dzs)
          if(swap.lt.0.25d0*dzs)then
            swap=0.25d0*dzs-swap
            zs1=zs1-swap
            zs2=zs2-swap
          else if(swap.gt.0.75d0*dzs)then
            swap=swap-0.75d0*dzs
            zs1=zs1+swap
            zs2=zs2+swap
          endif
        endif
      endif
c
      call skipdoc(10)
      read(10,*)nt,twindow
      if(twindow.le.0.d0)then
        stop ' Error in input: wrong time window!'
      else if(nt.le.0)then
        stop ' Error in input: time sampling no <= 0!'
      endif
      twindow=twindow*day2sec
      if(nt.le.2)then
        dt=twindow
      else
        dt=twindow/dble(nt-1)
      endif
c
c     wavenumber integration parameters
c     =================================
c
      call skipdoc(10)
      read(10,*)accuracy
      if(accuracy.le.0.d0.or.accuracy.ge.1.d0)accuracy=0.1d0
c
      call skipdoc(10)
      read(10,*)grfac
      if(grfac.le.grfacmin)then
        grfac=0.d0
      endif
c
c     parameters for output files
c     ===========================
c
      call skipdoc(10)
      read(10,*)outdir
c
      do lend=80,1,-1
        if(outdir(lend:lend).ne.' ')goto 100
      enddo
100   continue
c
      if(lend.lt.1)then
        stop 'Error: wrong format for output directory!'
      endif
c
      call skipdoc(10)
      read(10,*)(fname(i),i=1,3)
      call skipdoc(10)
      read(10,*)(fname(i),i=4,9)
      call skipdoc(10)
      read(10,*)(fname(i),i=10,14)
      do i=1,14
        do lenf=80,1,-1
          if(fname(i)(lenf:lenf).ne.' ')goto 110
        enddo
110     continue
        green(i,1)=outdir(1:lend)//fname(i)(1:lenf)//'.ep'
        green(i,2)=outdir(1:lend)//fname(i)(1:lenf)//'.ss'
        green(i,3)=outdir(1:lend)//fname(i)(1:lenf)//'.ds'
        green(i,4)=outdir(1:lend)//fname(i)(1:lenf)//'.cl'
        do istp=1,4
          select(i,istp)=.true.
        enddo
      enddo
c
c     no tangential components for clvd sources
c
      select(3,1)=.false.
      select(8,1)=.false.
      select(9,1)=.false.
      select(11,1)=.false.
      select(12,1)=.false.
      select(3,4)=.false.
      select(8,4)=.false.
      select(9,4)=.false.
      select(11,4)=.false.
      select(12,4)=.false.
c
c     global model parameters
c     =======================
c
      call skipdoc(10)
      read(10,*)l
c
      allocate(h(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: h not allocated!'
      allocate(rho(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: rho not allocated!'
      allocate(la(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: la not allocated!'
      allocate(mu(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: mu not allocated!'
      allocate(etk(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: etk not allocated!'
      allocate(etm(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: etm not allocated!'
      allocate(alf(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: alf not allocated!'
c
c      multilayered model parameters
c      =============================
c
      kgmax=0.d0
      do i=1,l
        call skipdoc(10)
        read(10,*)j,h(i),vp,vs,rho(i),etk(i),etm(i),alf(i)
        if(alf(i).gt.1.d0.or.alf(i).le.0.d0)then
          stop 'Error in psggetinp: wrong value for parameter alpha!'
        endif
        h(i)=h(i)*km2m
        vp=vp*km2m
        vs=vs*km2m
        mu(i)=rho(i)*vs*vs
        la(i)=rho(i)*vp*vp-2.d0*mu(i)
        if(la(i).le.0.d0)then
          stop 'inconsistent Vp/Vs ratio!'
        endif
        if(etk(i).le.0.d0.or.alf(i).eq.1.d0)then
          etk(i)=0.d0
          alf(i)=1.d0
        endif
        if(etm(i).lt.0.d0)then
          etm(i)=0.d0
        endif
        kg=rho(i)*g0/(la(i)+2.d0*mu(i)/3.d0)
        kgmax=dmax1(kgmax,kg)
      enddo
      if(l.eq.1)h(l)=0.d0
c
c     end of inputs
c     =============
c
      close(10)
c
c     determine upper und lower parameter values of each layer
c
      allocate(z1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: z1 not allocated!'
      allocate(rho1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: rho1 not allocated!'
      allocate(la1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: la1 not allocated!'
      allocate(mu1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: mu1 not allocated!'
      allocate(etk1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: etk1 not allocated!'
      allocate(etm1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: etm1 not allocated!'
      allocate(alf1(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: alf1 not allocated!'
c
      allocate(z2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: z2 not allocated!'
      allocate(rho2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: rho2 not allocated!'
      allocate(la2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: la2 not allocated!'
      allocate(mu2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: mu2 not allocated!'
      allocate(etk2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: etk2 not allocated!'
      allocate(etm2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: etm2 not allocated!'
      allocate(alf2(l),stat=ierr)
      if(ierr.ne.0)stop ' Error in psggetinp: alf2 not allocated!'
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          la1(l0)=la(i-1)
          mu1(l0)=mu(i-1)
          rho1(l0)=rho(i-1)
          etk1(l0)=etk(i-1)
          etm1(l0)=etm(i-1)
          alf1(l0)=alf(i-1)
c
          z2(l0)=h(i)
          la2(l0)=la(i)
          mu2(l0)=mu(i)
          rho2(l0)=rho(i)
          etk2(l0)=etk(i)
          etm2(l0)=etm(i)
          alf2(l0)=alf(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          la1(l0)=la(i)
          mu1(l0)=mu(i)
          rho1(l0)=rho(i)
          etk1(l0)=etk(i)
          etm1(l0)=etm(i)
          alf1(l0)=alf(i)
        endif
      enddo
      z1(l0)=h(l)
      la1(l0)=la(l)
      mu1(l0)=mu(l)
      rho1(l0)=rho(l)
      etk1(l0)=etk(l)
      etm1(l0)=etm(l)
      alf1(l0)=alf(l)
c
      deallocate(h,rho,la,mu,etk,etm,alf)
c
      return
      end