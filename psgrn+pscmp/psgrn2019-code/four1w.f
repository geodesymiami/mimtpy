      subroutine four1w(cdat,ddat,nn,isign)
      implicit none
      integer*4 nn,isign
      complex*16 cdat(nn)
      real*8 ddat(2*nn)
c
c     fast fourier transform (fft)
c     convention: f(t)=\int F(f)e^{-i2\pi ft} df 
c     replace ddat by its discrete fourier transform, if isign is input
c     as 1; or replace data by nn times its inverse discrete fourier
c     transform, if isign is input as -1. data is a double complex array of
c     length nn or, equivalently, a double precision array of length 2*nn.
c     nn must be an integer power of 2 (this is not checked!)
c
c     note for convention: f(t)=\int F(f)e^{i2\pi f t} df, t-domain ==>
c     f-domain, if isign=-1, and f-domain ==> t-domain, if isign=1.
c
      integer*4 i,j,n,m,mmax,istep
      real*8 tempr,tempi
      real*8 wr,wi,wpr,wpi,wtemp,theta
c
      do i=1,nn
        ddat(2*i-1)=dreal(cdat(i))
        ddat(2*i)=dimag(cdat(i))
      enddo
c
      n=2*nn
      j=1
      do i=1,n,2
        if(j.gt.i)then
          tempr=ddat(j)
          tempi=ddat(j+1)
          ddat(j)=ddat(i)
          ddat(j+1)=ddat(i+1)
          ddat(i)=tempr
          ddat(i+1)=tempi
        endif
        m=n/2
1       if((m.ge.2).and.(j.gt.m))then
          j=j-m
          m=m/2
          goto 1
        endif
        j=j+m
      enddo
      mmax=2
2     if(n.gt.mmax)then
        istep=2*mmax
        theta=6.28318530717959d0/dble(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do m=1,mmax,2
          do i=m,n,istep
            j=i+mmax
            tempr=wr*ddat(j)-wi*ddat(j+1)
            tempi=wr*ddat(j+1)+wi*ddat(j)
            ddat(j)=ddat(i)-tempr
            ddat(j+1)=ddat(i+1)-tempi
            ddat(i)=ddat(i)+tempr
            ddat(i+1)=ddat(i+1)+tempi
          enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
        enddo
        mmax=istep
        goto 2
      endif
c
      do i=1,nn
        cdat(i)=dcmplx(ddat(2*i-1),ddat(2*i))
      enddo
c
      return
      end
