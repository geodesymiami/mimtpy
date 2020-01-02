      module psgalloc
c===================================================================
c     global constants
c===================================================================
      real*8 km2m,day2sec,relaxmin,grfacmin
      parameter(km2m=1.0d+03,day2sec=8.64d+04,relaxmin=1.0d-06,
     &          grfacmin=1.0d-03)
c
c     gravity, gravitational constant, earth radius
c     ==============================================================
c     gamma = 4*pi*G
c
      real*8 g0,gamma,rearth,denswater
      parameter(g0=9.82d+00,gamma=8.38579d-10,rearth=6.371d+06,
     &          denswater=1.0d+03)
c
c     parameters of Bessel functions
c     ==============================================================
c
      integer*4 dnx,nxmax,nbsjmax
      parameter(dnx=512,nxmax=1024,nbsjmax=dnx*nxmax)
c
c     resolution parameters for discetize layers with gradient
c     ==============================================================
c     reslm: for moduli
c     resld: for density
c     reslv: for viscosity
c
      real*8 reslm,resld,reslv
      parameter(reslm=0.5d-01,resld=0.5d-01,reslv=0.25d+00)
c
c     parameters of source functions
c     ==============================================================
c
      integer*8 ms(4),ics(4)
      real*8 cics(4),cms(4)
      complex*16 sfct0(8,4),sfct1(8,4),sfcths0(8,4),sfcths1(8,4)
      logical*2 select(14,4)
c===================================================================
c     global indices
c===================================================================
c     lzrec: layer index of receiver location
c     n0: number of homogeneous layers
c     lp: number of layers + pseudo layers (source and receiver)
c     nr: number of traces
c     nfmax: min. number of frequency samples
c     nfmax: max. number of frequency samples
c     nzs: number of source depths
c     nt: number of time samples
c
      integer*4 nfmin,nfmax,nzsmax
      integer*4 nr,nt,lp,lzrec,ioc,l0,n0,ls,nzs
c===================================================================
c     working space
c===================================================================
      logical*2 nongravity
      integer*4 nwarn
      integer*4 unit(14,4)
      real*8 twindow,taumin,grfac,zrec,kgmax,zs,dxbsj,dt,accuracy
      real*8 r1,r2,dr,zs1,zs2,dzs,zrs2,sampratio
      character*35 stype(4)
      character*35 comptxt(14)
      character*80 fname(14)
      character*163 green(14,4)
      character*80 inputfile,outdir
c
      integer*4, allocatable:: nno(:)
      real*8, allocatable:: zp(:),hp(:)
c===================================================================
c     original input model parameters
c===================================================================
      real*8, allocatable:: z1(:),z2(:),
     &        la1(:),la2(:),mu1(:),mu2(:),
     &        rho1(:),rho2(:),etk1(:),etk2(:),
     &        etm1(:),etm2(:),alf1(:),alf2(:)
c===================================================================
c     re-processed model parameters
c===================================================================
      logical*2, allocatable::  elastic(:)
      real*8, allocatable:: h(:),la(:),mu(:),rho(:),etk(:),etm(:),alf(:)
      real*8, allocatable:: bsj(:,:)
      complex*16, allocatable:: cla(:),cmu(:)
      complex*16, allocatable:: hkup(:,:,:),hklw(:,:,:)
      complex*16, allocatable:: maup(:,:,:),maiup(:,:,:),
     &                          malw(:,:,:),mailw(:,:,:)
c===================================================================
c     output data
c===================================================================
      integer*4, allocatable:: idec(:),nout(:)
      real*8, allocatable:: r(:),rs(:),geow(:)
      complex*16, allocatable:: obs(:,:,:),obs0(:,:,:),du(:,:,:,:)
      real*8, allocatable:: tgrn(:),dobs(:),dswap(:)
      complex*16, allocatable:: fgrn(:)
c
      end module