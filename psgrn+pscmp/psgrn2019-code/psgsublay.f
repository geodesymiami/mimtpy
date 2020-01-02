	subroutine psgsublay(ierr)
      use psgalloc
	implicit none
c
	integer*4 ierr
c
c	work space
c
	integer*4 i,l,n
	real*8 dh,dla,dmu,drho,detk,detm,dalf,z,dz,tau
c
      integer*4, allocatable:: i0(:)
c
      allocate(i0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: nout not allocated!'
c
      n0=0
      do l=1,l0-1
	  dla=2.d0*dabs(la2(l)-la1(l))/(la2(l)+la1(l))
	  dmu=2.d0*dabs(mu2(l)-mu1(l))/(mu2(l)+mu1(l))
        if(rho2(l)+rho1(l).gt.0.d0)then
	    drho=2.d0*dabs(rho2(l)-rho1(l))/(rho2(l)+rho1(l))
        else
          drho=0.d0
        endif
	  if(etk2(l)+etk1(l).gt.0.d0)then
          detk=2.d0*dabs(etk2(l)-etk1(l))/(etk2(l)+etk1(l))
        else
          detk=0.d0
        endif
	  if(etm2(l)+etm1(l).gt.0.d0)then
          detm=2.d0*dabs(etm2(l)-etm1(l))/(etm2(l)+etm1(l))
        else
          detm=0.d0
        endif
        if(alf2(l)+alf1(l).gt.0.d0)then
	    dalf=2.d0*dabs(alf2(l)-alf1(l))
     &            /(alf2(l)+alf1(l))
        else
          dalf=0.d0
        endif
	  i0(l)=idnint(dmax1(1.d0,dla/reslm,dmu/reslm,drho/resld,
     &                     detk/reslv,detm/reslv,dalf/reslv))
        n0=n0+i0(l)
      enddo
c
      n0=n0+1
c
      allocate(h(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: h not allocated!'
      allocate(rho(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: rho not allocated!'
      allocate(la(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: la not allocated!'
      allocate(mu(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: mu not allocated!'
      allocate(etk(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: etk not allocated!'
      allocate(etm(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: etm not allocated!'
      allocate(alf(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: alf not allocated!'
c
      allocate(cla(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: la not allocated!'
      allocate(cmu(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: mu not allocated!'
      allocate(elastic(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgsublay: elastic not allocated!'
c
	n=0
c
	do l=1,l0-1
	  dz=z2(l)-z1(l)
	  dla=(la2(l)-la1(l))/dz
	  dmu=(mu2(l)-mu1(l))/dz
	  drho=(rho2(l)-rho1(l))/dz
	  detk=(etk2(l)-etk1(l))/dz
	  detm=(etm2(l)-etm1(l))/dz
	  dalf=(alf2(l)-alf1(l))/dz
	  dh=dz/dble(i0(l))
	  do i=1,i0(l)
	    n=n+1
	    h(n)=dh
	    z=(dble(i)-0.5d0)*dh
	    la(n)=la1(l)+dla*z
	    mu(n)=mu1(l)+dmu*z
	    rho(n)=rho1(l)+drho*z
	    etk(n)=etk1(l)+detk*z
	    etm(n)=etm1(l)+detm*z
	    alf(n)=alf1(l)+dalf*z
	  enddo
	enddo
c
c	last layer is half-space
c
	n=n+1
	h(n)=0.d0
	la(n)=la1(l0)
	mu(n)=mu1(l0)
	rho(n)=rho1(l0)
	etk(n)=etk1(l0)
	etm(n)=etm1(l0)
	alf(n)=alf1(l0)
c
      write(*,*)'the multi-layered poroelastic model:'
	write(*,'(8a)')'  no',' thick(m)    ','  la(Pa)    ',
     &    '  mu(Pa)    ','rho(kg/m^3) ','  etk(Pa*s) ',
     &    '  etm(Pa*s) ','   alpha'
      taumin=0.d0
	do n=1,n0
	  write(*,1001)n,h(n),la(n),mu(n),
     &               rho(n),etk(n),etm(n),alf(n)
        elastic(n)=(etk(n).le.0.d0.or.alf(n).ge.1.d0).and.
     &             etm(n).le.0.d0
        if(.not.elastic(n))then
          if(etk(n).le.0.d0.or.alf(n).ge.1.d0)then
c           Maxwell body
            tau=etm(n)/mu(n)
          else
            tau=dmin1(etk(n)*(1.d0-alf(n))/alf(n),etm(n))/mu(n)
          endif
          if(taumin.le.0.d0.or.taumin.gt.tau)taumin=tau
        endif
	enddo
1001	format(i4,f11.4,6E12.4)
	ierr=0
c
      deallocate(i0)
c
	return
	end
