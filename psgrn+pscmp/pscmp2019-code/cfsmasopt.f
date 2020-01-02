      subroutine cfsmasopt(sxx,syy,szz,sxy,syz,szx,
     &                     p,f,cfs,sig,strike,dip,rakeopt)
      implicit none
c
c     calculate Coulomb stress
c
c     input:
c     stress tensor, pore pressure, friction coefficient
c     fault parameters (strike, dip)
c
      real*8 sxx,syy,szz,sxy,syz,szx,p,f,cfs,sig,strike,dip,rakeopt
c
c     return:
c     max. Coulomb failure stress (cfs) in an optimal (rakeopt) orientation
c     on the given fault plane (strike, dip), normal stress (sig)
c
c     local memories:
c
      integer*4 i,j
      real*8 pi,deg2rad,st0,di0,tau
      real*8 s(3,3),ns(3),rst(3),rdi(3),sts(3)
c
      pi=4.d0*datan(1.d0)
      deg2rad=pi/180.d0
      st0=strike*deg2rad
      di0=dip*deg2rad
c
      s(1,1)=sxx
      s(1,2)=sxy
      s(1,3)=szx
      s(2,1)=sxy
      s(2,2)=syy
      s(2,3)=syz
      s(3,1)=szx
      s(3,2)=syz
      s(3,3)=szz
c
      ns(1)=dsin(di0)*dcos(st0+0.5d0*pi)
      ns(2)=dsin(di0)*dsin(st0+0.5d0*pi)
      ns(3)=-dcos(di0)
c
      rst(1)=dcos(st0)
      rst(2)=dsin(st0)
      rst(3)=0.d0
c
      rdi(1)=dcos(di0)*dcos(st0+0.5d0*pi)
      rdi(2)=dcos(di0)*dsin(st0+0.5d0*pi)
      rdi(3)=dsin(di0)
c
      do i=1,3
        sts(i)=0.d0
        do j=1,3
          sts(i)=sts(i)+s(i,j)*ns(j)
        enddo
      enddo
c
      sig=0.d0
      do i=1,3
        sig=sig+sts(i)*ns(i)
      enddo
c
      do i=1,3
        sts(i)=sts(i)-sig*ns(i)
      enddo
      tau=dsqrt(sts(1)**2+sts(2)**2+sts(3)**2)
c
      rakeopt=datan2(-sts(1)*rdi(1)-sts(2)*rdi(2)-sts(3)*rdi(3),
     &                sts(1)*rst(1)+sts(2)*rst(2)+sts(3)*rst(3))
     &       /deg2rad
      rakeopt=dmod(rakeopt+360.d0,360.d0)
c
      cfs=tau+f*(sig+p)
      return
      end