      subroutine pscgetinp(ierr)
      use pscalloc
      implicit none
c
c     Last modified: Potsdam, Feb, 2019, by R. Wang
c
      integer*4 ierr
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 i,j,is,isnap,irec,ntout
      integer*4 iptch,ieq,itwas,nx,ny
      integer*4 ilatlon,ilatrec,ilonrec,nlatrec,nlonrec
      real*8 latrec1,latrec2,lonrec1,lonrec2,dlatrec,dlonrec
      real*8 sx,sy,sz,swap(10)
      complex*16 clatlon1,clatlon2
      logical*2 neweq
c
      integer ftell
c
      open(10,file=inputfile,status='old')
c
      call skipdoc(10)
      read(10,*)ilatlon
      if(ilatlon.eq.0)then
c
c       irregular observation positions
c
        call skipdoc(10)
        read(10,*)nrec
      else if(ilatlon.eq.1)then
c
c       1D observation profile
c
        call skipdoc(10)
        read(10,*)nrec
        call skipdoc(10)
        read(10,*)clatlon1,clatlon2
      else if(ilatlon.eq.2)then
c
c       2D rectanglar observation array
c
        call skipdoc(10)
        read(10,*)nlatrec,latrec1,latrec2
        call skipdoc(10)
        read(10,*)nlonrec,lonrec1,lonrec2
        nrec=nlatrec*nlonrec
      else
        stop' Error: wrong input for ilatlon!'
      endif
c
      if(nrec.lt.1)then
        stop ' Error: wrong input for nrec!'
      endif
c
      allocate(clatlon(nrec),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: clatlon not allocated!'
      allocate(latrec(nrec),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: latrec not allocated!'
      allocate(lonrec(nrec),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: lonrec not allocated!'
      allocate(obs1(nrec,14),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: obs1 not allocated!'
      allocate(obs2(nrec,14),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: obs2 not allocated!'
c
      if(ilatlon.eq.0)then
c
c       irregular observation positions
c
        read(10,*)(clatlon(irec),irec=1,nrec)
        do irec=1,nrec
          latrec(irec)=dreal(clatlon(irec))
          lonrec(irec)=dimag(clatlon(irec))
        enddo
      else if(ilatlon.eq.1)then
c
c        1D observation profile
c
        latrec(1)=dreal(clatlon1)
        lonrec(1)=dimag(clatlon1)
        if(nrec.gt.1)then
          dlatrec=dreal(clatlon2-clatlon1)/dble(nrec-1)
          dlonrec=dimag(clatlon2-clatlon1)/dble(nrec-1)
        else
          dlatrec=0.d0
          dlonrec=0.d0
        endif
        do irec=1,nrec
          latrec(irec)=dreal(clatlon1)+dlatrec*dble(irec-1)
          lonrec(irec)=dimag(clatlon1)+dlonrec*dble(irec-1)
        enddo
      else if(ilatlon.eq.2)then
c
c       2D rectanglar observation array
c
        if(nlatrec.gt.1)then
          dlatrec=(latrec2-latrec1)/dble(nlatrec-1)
        else
          dlatrec=0.d0
        endif
        if(nlonrec.gt.1)then
          dlonrec=(lonrec2-lonrec1)/dble(nlonrec-1)
        else
          dlonrec=0.d0
        endif
        irec=0
        do ilonrec=1,nlonrec
          do ilatrec=1,nlatrec
            irec=irec+1
            latrec(irec)=latrec1+dlatrec*dble(ilatrec-1)
            lonrec(irec)=lonrec1+dlonrec*dble(ilonrec-1)
          enddo
        enddo
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c      READ IN OUTPUT PARAMETERS
c      =========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(10)
      read(10,*)insar,xlos,ylos,zlos
      call skipdoc(10)
      read(10,*)icfs,friction,skempton,
     &               strike0,dip0,rake0,(swap(j),j=1,3)
      if(icfs.eq.1)then
        sigma0(1)=dmax1(swap(1),swap(2),swap(3))
        sigma0(3)=dmin1(swap(1),swap(2),swap(3))
        sigma0(2)=swap(1)+swap(2)+swap(3)-sigma0(1)-sigma0(3)
      endif
      call skipdoc(10)
      read(10,*)outdir
c
      call skipdoc(10)
      read(10,*)(itout(i),i=1,3)
      call skipdoc(10)
      read(10,*)(toutfile(i),i=1,3)
      call skipdoc(10)
      read(10,*)(itout(i),i=4,9)
      call skipdoc(10)
      read(10,*)(toutfile(i),i=4,9)
      call skipdoc(10)
      read(10,*)(itout(i),i=10,14)
      call skipdoc(10)
      read(10,*)(toutfile(i),i=10,14)
      call skipdoc(10)
      read(10,*)nsnap
c
      allocate(tsnap(nsnap),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: tsnap not allocated!'
      allocate(scoutfile(nsnap),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: tsnap not allocated!'
c
      do isnap=1,nsnap
        call skipdoc(10)
        read(10,*)tsnap(isnap),scoutfile(isnap)
        if(tsnap(isnap).lt.0.d0)then
          stop ' Error: wrong scenario time!'
        endif
        tsnap(isnap)=tsnap(isnap)*DAY2SEC
      enddo
      onlysc=.true.
      do i=1,14
        onlysc=onlysc.and.itout(i).ne.1
      enddo
      if(onlysc.and.nsnap.le.0)then
        stop ' No outputs have been selected!'
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR EARTH MODEL CHOICE
c     =========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(10)
      read(10,*)iesmodel
      if(iesmodel.eq.0)then
        call skipdoc(10)
        read(10,*)zrec,larec,murec
        zrec=zrec*KM2M
      else if(iesmodel.eq.1)then
        call skipdoc(10)
        read(10,*)grndir
        call skipdoc(10)
        read(10,*)(green(i),i=1,3)
        call skipdoc(10)
        read(10,*)(green(i),i=4,9)
        call skipdoc(10)
        read(10,*)(green(i),i=10,14)
      else
        stop ' Error: wrong selection of earth structure model!'
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR RECTANGULAR SOURCES
c     ==========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(10)
      read(10,*)ns
      if(ns.lt.1)then
        stop ' Error: wrong number of subfaults!'
      endif
c
      allocate(latref(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: latref not allocated!'
      allocate(lonref(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: lonref not allocated!'
      allocate(zref(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: zref not allocated!'
      allocate(length(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: length not allocated!'
      allocate(width(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: width not allocated!'
      allocate(strike(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: strike not allocated!'
      allocate(dip(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: dip not allocated!'
      allocate(nptch_s(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: nptch_s not allocated!'
      allocate(nptch_d(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: nptch_d not allocated!'
      allocate(tstart(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: tstart not allocated!'
      allocate(eqtime(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: eqtime not allocated!'
      allocate(ieqno(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: ieqno not allocated!'
c
c     current position of reading
c
      itwas=ftell(10)
c
      neq=0
      nptch=0
      do is=1,ns
        call skipdoc(10)
        read(10,*)i,(swap(j),j=1,7),nx,ny,tstart(is)
        nptch=max0(nptch,nx*ny)
        do iptch=1,nx*ny
          call skipdoc(10)
          read(10,*)(swap(j),j=1,5)
        enddo
        neweq=.true.
        do ieq=1,neq
          if(tstart(is).eq.eqtime(ieq))then
            neweq=.false.
            ieqno(is)=ieq
          endif
        enddo
        if(neweq)then
          neq=neq+1
          ieqno(is)=neq
          eqtime(neq)=tstart(is)
        endif
      enddo
c
      allocate(ptch_s(ns,nptch),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: ptch_s not allocated!'
      allocate(ptch_d(ns,nptch),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: ptch_d not allocated!'
      allocate(slip_s(ns,nptch),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: slip_s not allocated!'
      allocate(slip_d(ns,nptch),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: slip_d not allocated!'
      allocate(opening(ns,nptch),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscgetinp: opening not allocated!'
c
c     re-read finite fault data
c
      call fseek(10,itwas,0,i)
c
      do is=1,ns
        call skipdoc(10)
        read(10,*)i,latref(is),lonref(is),zref(is),
     &      length(is),width(is),strike(is),dip(is),
     &      nptch_s(is),nptch_d(is),tstart(is)
c
        zref(is)=zref(is)*KM2M
        length(is)=length(is)*KM2M
        width(is)=width(is)*KM2M
        tstart(is)=tstart(is)*DAY2SEC
        do iptch=1,nptch_s(is)*nptch_d(is)
          call skipdoc(10)
          read(10,*)ptch_s(is,iptch),ptch_d(is,iptch),sx,sy,sz
          ptch_s(is,iptch)=ptch_s(is,iptch)*KM2M
          ptch_d(is,iptch)=ptch_d(is,iptch)*KM2M
          slip_s(is,iptch)=sx
          slip_d(is,iptch)=sy
          opening(is,iptch)=sz
        enddo
        if(zref(is).lt.0.d0)then
          stop ' Error: source depth zs0 < 0!'
        endif
        if(dabs(strike(is)).gt.360.d0)then
          stop ' Error: wrong strike angle!'
        endif
        if(strike(is).lt.0.d0)then
          strike(is)=strike(is)+360.d0
        endif
        if(dip(is).gt.90.d0.or.dip(is).lt.0.d0)then
          stop ' Error: wrong dip angle!'
        endif
      enddo
c
      torigin=tstart(1)
      do is=2,ns
        torigin=dmin1(torigin,tstart(is))
      enddo
c
      do is=1,ns
        tstart(is)=tstart(is)-torigin
      enddo
      do isnap=1,nsnap
        tsnap(isnap)=tsnap(isnap)-torigin
      enddo
c
      close(10)
c
      if(iesmodel.eq.0)then
        ntout=0
        do i=1,14
          if(itout(i).eq.1)then
            ntout=ntout+1
            itout(i)=0
          endif
        enddo
        if(ntout.gt.0)then
          print *,'Warning: no time series will be calculated '
     &          //'for the selected homogeneous elastic model!'
        endif
        onlysc=.true.
      endif
c
      return
      end
