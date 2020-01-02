      module pscalloc
c
c     Last modified: Potsdam, Feb, 2019, by R. Wang
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real*8 DEG2RAD,KM2M,DAY2SEC,REARTH,G0,PI
      parameter(DEG2RAD=1.745329252d-02,KM2M=1.0d+03)
      parameter(DAY2SEC=8.64d+04,REARTH=6.371d+06,G0=9.82d+00)
      parameter(PI=3.14159265358979d0)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL INDICES
c     ==============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     nrec = number of observation positions
c     nzs = number of discrete source depths
c     nr = number of discrete radial diatances
c     ns = number of fault segments
c     neq = number of earthquakes
c     nptch = number of patches at a source rectangle
c     nps = number of discrete point sources per source depth
c     nt = number of time samples used for Green's functions
c     ntr = number of time samples of the outputs
c     nsnap = number of scenario outputs (<= ntr/2)
c     nwarn = total number of warnings
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 nzs,nr,neq,ns,nptch,nrec,nt,ntr,nsnap,ntrec
      integer*4 nwarn
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNCTION INFO
c     =====================
c     nzs,zs1,zs2 = number of depth samples, start and end depths used
c           in Green's functions
c     nr,r1,r2 = number of distance samples, start and end distances used
c           in Green's functions
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 iesmodel
      real*8 r1,r2,sampratio,zs1,zs2,twindow,torigin
      real*8 zrec,larec,murec,rhorec,etkrec,etmrec,alfrec
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     COSINES OF LOS TO INSAR ORBIT
c     =============================
c     insar = 1: output los displacements
c             0: not output los displacements
c     xlos, ylos, zlos = cosines of the los
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 insar,icfs
      real*8 xlos,ylos,zlos
      real*8 friction,skempton,strike0,dip0,rake0,mw,mwsca
      real*8 sigma0(3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     RECTANGULAR SOURCE PLANES
c     =========================
c     (latref,lonref) = geographic coordinates of the reference point
c     zref = depth of the reference point.
c     all angles in degree.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4, allocatable:: nptch_s(:),nptch_d(:),ieqno(:)
      real*8, allocatable:: latref(:),lonref(:),zref(:),
     &    length(:),width(:),strike(:),dip(:),tstart(:),eqtime(:),
     &    ptch_s(:,:),ptch_d(:,:),slip_s(:,:),slip_d(:,:),opening(:,:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     DISTRETE POINT SOURCES
c     ======================
c     (xs,ys,zs) = coordinates of the discrete point sources
c     with x = north, y = east, z = downward
c     angles in degree.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4, allocatable:: nps(:),isno(:,:)
      real*8, allocatable:: plat(:,:),plon(:,:),pz(:,:),pmwei(:,:,:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OBSERVATION POSITIONS AND OBSERVABLES
c     =====================================
c     (latrec(i),lonrec(i))=coordinates of the observation positions
c     the 3 displcement/velocity/acceleration components: ux,uy,uz
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*80 grndir,green(14)
      character*80 inputfile,outdir,outputfile,toutfile(14)
      logical*2 onlysc
      integer*4 itout(14)
c
      integer*4, allocatable:: idec(:),igrns(:,:),itsnap(:)
      real*8, allocatable:: latrec(:),lonrec(:),tsnap(:),obs(:,:,:),
     &    coobs(:,:,:),poobs(:,:,:),obs1(:,:),obs2(:,:),
     &    r(:),cogrns(:,:,:),grns(:,:,:,:),
     &    mtenexp(:,:,:),mtenshr(:,:,:),mscaexp(:),mscashr(:)
      complex*16, allocatable:: clatlon(:)
      character*7, allocatable:: rtxt(:)
      character*80, allocatable:: scoutfile(:)
c
      end module