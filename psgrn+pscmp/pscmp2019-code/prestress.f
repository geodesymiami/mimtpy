      subroutine prestress(s1,s2,s3,strike,dip,rake,p,f,
     &                     cmb,sxx,syy,szz,sxy,syz,szx)
      implicit none
c
c     updated by R. Wang on Feb 8, 2019
c
c     determine regional stress tensor using the known principal
c     stresses and master fault mechanism, assuming that the master
c     fault follows the Coulomb failure criterion.
c
c     input:
c     principal stresses, master fault strike, dip and rake,
c     pore pressure, friction coefficient
c
      real*8 s1,s2,s3,strike,dip,rake,p,f
c
c     output
c     max. preseismic Coulomb stress, prestress tensor
c
      real*8 cmb,sxx,syy,szz,sxy,syz,szx
c
c     local memories:
c
      integer*4 i,j,k
      real*8 pi,alpha,st,di,ra,deg2rad,cmb1,cmb2,cmb3
      real*8 ns(3),ts(3),rst(3),rdi(3)
      real*8 sig(3),s(3,3),rot(3,3)
      real*8 st1,di1,ra1,st2,di2,ra2,sigma
c
      cmb=0.d0
      sxx=0.d0
      syy=0.d0
      szz=0.d0
      sxy=0.d0
      syz=0.d0
      szx=0.d0
c
      if(s1.eq.0.d0.and.s2.eq.0.d0.and.s3.eq.0.d0)return
c
      pi=4.d0*datan(1.d0)
      deg2rad=pi/180.d0
c
      st=strike*deg2rad
      di=dip*deg2rad
      ra=rake*deg2rad
c
c     sort the principal stresses: sig1 = max, sig2 = min, sig3 = middle.
c
      sig(1)=dmax1(s1,s2,s3)
      sig(2)=dmin1(s1,s2,s3)
      sig(3)=s1+s2+s3-sig(1)-sig(2)
c
c     max. Coulomb stress (proof see below)
c     on a plane perpendicular to (axis-sig1, axis-sig2)
c
      cmb=0.5d0*(sig(1)-sig(2))*dsqrt(1.d0+f*f)
     &   +f*(0.5d0*(sig(1)+sig(2))+p)
c
c     determine principal stress orientations
c
c     alpha=0.5*arctan(1/f): angle between the axis-sig2 (minimum principal stress) and
c     the rake direction (maximum Coulomb stress orientation), positive from max(sig) to max(sig)
c
c     proof: cmb = 0.5*(sig1-sig2)*sin(2*alpha) + f*(sig1*cos^2(alpha)+sig2*sin^2(alpha)+p)
c            maximize cmb => alpha=0.5*arctan(1/f) or 0.5*arctan(1/f)+/-90° (here use of the first one)
c            and max(cmb) = 0.5*(sig1-sig2)*sqrt(1+f*f)+f*(0.5*(sig1+sig2)+p)
c
      alpha=0.5d0*datan2(1.d0,f)
c
c     ns: normal of the hanging plate
c     rst: unit vector in the strike direction
c     rdi: unit vector in the down-dip direction
c     ts: unit vector in the rake direction, i.e., orintation of sig1 = max(cmb)
c         in the local Cartesian coordinate system (x: north, y: east, z: downward)
c
      ns(1)=dsin(di)*dcos(st+0.5d0*pi)
      ns(2)=dsin(di)*dsin(st+0.5d0*pi)
      ns(3)=-dcos(di)
c
      rst(1)=dcos(st)
      rst(2)=dsin(st)
      rst(3)=0.d0
c
      rdi(1)=dcos(di)*dcos(st+0.5d0*pi)
      rdi(2)=dcos(di)*dsin(st+0.5d0*pi)
      rdi(3)=dsin(di)
c
      do i=1,3
        ts(i)=rst(i)*dcos(ra)-rdi(i)*dsin(ra)
      enddo
c
c     ns-ts_plane = sig1-sig2_plane
c     rot(i,1): unit vector of x' = axis-sig1,
c     rot(i,2): unit vector of y' = axis-sig2,
c     rot(i,3) = rot(i,1) x rot(i,2), equivalent to ns x ts = unit vector of z' = axis-sig3,
c     => rot = coordinate transformation matrix (x',y',z') => (x,y,z)
c     e.g., rot(1,2) = cos(x-axis,y'-axis)
c
      do i=1,3
        rot(i,1)=ts(i)*dsin(alpha)+ns(i)*dcos(alpha)
        rot(i,2)=ts(i)*dcos(alpha)-ns(i)*dsin(alpha)
      enddo
      rot(1,3)=rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2)
      rot(2,3)=rot(3,1)*rot(1,2)-rot(1,1)*rot(3,2)
      rot(3,3)=rot(1,1)*rot(2,2)-rot(2,1)*rot(1,2)
c
c     coordinate transformation from the principal stress system
c     (x' = sig1, y' = sig2, z' = sig3)to the local Cartesian coordinate
c     system (x = north, y = east and z = downward)
c
      do i=1,3
        do j=1,3
          s(i,j)=0.d0
          do k=1,3
            s(i,j)=s(i,j)+sig(k)*rot(i,k)*rot(j,k)
          enddo
        enddo
      enddo
c
      sxx=s(1,1)
      syy=s(2,2)
      szz=s(3,3)
      sxy=s(1,2)
      syz=s(2,3)
      szx=s(1,3)
c
      return
      end