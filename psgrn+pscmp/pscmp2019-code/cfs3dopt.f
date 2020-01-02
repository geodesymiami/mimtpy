      subroutine cfs3dopt(sxx,syy,szz,sxy,syz,szx,p,f,key,
     &                    st0,di0,ra0,
     &                    cfs,sig,st1,di1,ra1,st2,di2,ra2)
      implicit none
c
c     updated by Rongjiang Wang, GFZ, Feb 8, 2019
c
c     Coulomb failure stress with the optimal orientation
c
c     input:
c     sxx,syy,szz,sxy,syz,szx = stress tensor
c     p = pore pressure
c     f = friction coefficient
c     st0,di0,ra0 = strike, dip and rake [deg] of master fault
c     key = 0: determine optimal Coulomb stress only;
c           1: determine optimal Coulomb stress and orientations
c
      integer*4 key
      real*8 sxx,syy,szz,sxy,syz,szx,p,f,st0,di0,ra0
c
c     output
c     cfs = max. Coulomb failure stress, and
c     sig = normal stress on the two optimally oriented fault planes
c     st1,di1,ra1 = strike, dip and rake [deg] of the first optimally oriented fault
c                   which is the closest to the master fault
c     st2,di2,ra2 = strike, dip and rake [deg] of the 2. optimally oriented fault.
c                   Note:
c                   the two optimally oriented faults have an angle of arctan(1/f),
c                   or 180°-arctan(1/f). they become perpendicular only for f -> 0.
c
      real*8 cfs,sig,st1,di1,ra1,st2,di2,ra2
c
c     local memories:
c
      integer*4 i,j
      real*8 pi,rad2deg,b,c,d,s1,s2,s3,alpha,am,swap,det1,det2,det3,eps
      real*8 s(3),st(2),di(2),ra(2),rst(3),rdi(3)
      real*8 r(3,2),ns(3,2),ts(3,2)
      real*8 mscorr
c
      pi=4.d0*datan(1.d0)
      rad2deg=180.d0/pi
c
      if(sxy.eq.0.d0.and.syz.eq.0.d0.and.szx.eq.0.d0)then
        s(1)=sxx
        s(2)=syy
        s(3)=szz
      else
        b=-(sxx+syy+szz)
        c=sxx*syy+syy*szz+szz*sxx-sxy**2-syz**2-szx**2
        d=sxx*syz**2+syy*szx**2+szz*sxy**2-2.d0*sxy*syz*szx-sxx*syy*szz
        call roots3(b,c,d,s)
      endif
c
      s1=dmax1(s(1),s(2),s(3))
      s2=dmin1(s(1),s(2),s(3))
      s3=s(1)+s(2)+s(3)-s1-s2
c
      s(1)=s1
      s(2)=s2
      s(3)=s3
c
      eps=1.0d-10*dsqrt(s1**2+s2**2+s3**2)
c
      sig=0.5d0*((s1-s2)*f/dsqrt(1+f*f)+s1+s2)
      cfs=0.5d0*(s1-s2)*dsqrt(1+f*f)+f*(0.5d0*(s1+s2)+p)
c
      if(key.eq.0.or.dabs(s1-s2).le.eps)then
c
c       (nearly) isotropic stress state
c
        st1=0.d0
        di1=0.d0
        ra1=0.d0
c
        st2=0.d0
        di2=0.d0
        ra2=0.d0
        return
      endif
c
c     determine eigenvectors of the max. and min. eigenvalues (principal stresses)
c
      do j=1,2
        if(dabs(s(j)-s(3)).le.eps)then
c
c         two of the three principal stresses are identical
c         -> more than a pair of optimal oriented fault planes
c         here any arbitrary pair is chosen
c
          am=dmax1(dabs(sxx-s(j)),dabs(sxy),dabs(szx))
          if(dabs(sxx-s(j)).ge.am)then
            r(1,j)=-(sxy+szx)
            r(2,j)=sxx-s(j)
            r(3,j)=sxx-s(j)
          else if(dabs(sxy).ge.am)then
            r(1,j)=sxy
            r(2,j)=-(sxx-s(j)+szx)
            r(3,j)=sxy
          else
            r(1,j)=szx
            r(2,j)=-(sxx-s(j)+sxy)
            r(3,j)=sxy
          endif
        else
          det1=syz*syz-(syy-s(j))*(szz-s(j))
          det2=szx*szx-(sxx-s(j))*(szz-s(j))
          det3=sxy*sxy-(sxx-s(j))*(syy-s(j))
          am=dmax1(dabs(det1),dabs(det2),dabs(det3))
          if(dabs(det1).ge.am)then
            r(1,j)=det1
            r(2,j)=(szz-s(j))*sxy-syz*szx
            r(3,j)=(syy-s(j))*szx-syz*sxy
          else if(dabs(det2).ge.am)then
            r(1,j)=(szz-s(j))*sxy-szx*syz
            r(2,j)=det2
            r(3,j)=(sxx-s(j))*syz-szx*sxy
          else
            r(1,j)=(syy-s(j))*szx-sxy*syz
            r(2,j)=(sxx-s(j))*syz-sxy*szx
            r(3,j)=det3
          endif
        endif
        am=dsqrt(r(1,j)**2+r(2,j)**2+r(3,j)**2)
        do i=1,3
          r(i,j)=r(i,j)/am
        enddo
      enddo
c
      alpha=0.5d0*datan2(1.d0,f)
c
c     determine the two optimal fault-plane normals, ns,
c     their unit vectors in the rake direction, ts.
c
      do i=1,3
        ns(i,1)=r(i,1)*dcos(alpha)-r(i,2)*dsin(alpha)
        ts(i,1)=r(i,1)*dsin(alpha)+r(i,2)*dcos(alpha)
c
        ns(i,2)=r(i,1)*dcos(alpha)+r(i,2)*dsin(alpha)
        ts(i,2)=r(i,1)*dsin(alpha)-r(i,2)*dcos(alpha)
      enddo
c
      do j=1,2
        if(ns(3,j).gt.0.d0)then
          do i=1,3
            ns(i,j)=-ns(i,j)
            ts(i,j)=-ts(i,j)
          enddo
        endif
      enddo
c
c     determine strike, dip and rake
c
      do j=1,2
        st(j)=datan2(ns(2,j),ns(1,j))-0.5d0*pi
        di(j)=datan2(dsqrt(ns(1,j)**2+ns(2,j)**2),dabs(ns(3,j)))
c
c       rst = unit vector in the strike direction
c       rdi = unit vector in the down-dip direction
c
        rst(1)=dcos(st(j))
        rst(2)=dsin(st(j))
        rst(3)=0.d0
c
        rdi(1)=dcos(di(j))*dcos(st(j)+0.5d0*pi)
        rdi(2)=dcos(di(j))*dsin(st(j)+0.5d0*pi)
        rdi(3)=dsin(di(j))
c
        ra(j)=datan2(-(ts(1,j)*rdi(1)+ts(2,j)*rdi(2)+ts(3,j)*rdi(3)),
     &                (ts(1,j)*rst(1)+ts(2,j)*rst(2)+ts(3,j)*rst(3)))
      enddo
c
      st1=dmod(st(1)*rad2deg+360.d0,360.d0)
      di1=di(1)*rad2deg
      ra1=dmod(ra(1)*rad2deg+360.d0,360.d0)
c
      st2=dmod(st(2)*rad2deg+360.d0,360.d0)
      di2=di(2)*rad2deg
      ra2=dmod(ra(2)*rad2deg+360.d0,360.d0)
c
      if(mscorr(st0,di0,ra0,st1,di1,ra1).lt.
     &   mscorr(st0,di0,ra0,st2,di2,ra2))then
        swap=st1
        st1=st2
        st2=swap
c
        swap=di1
        di1=di2
        di2=swap
c
        swap=ra1
        ra1=ra2
        ra2=swap
      endif
      return
      end