      subroutine psghksh(hk,m,k,z,n)
      use psgalloc
      implicit none
c
      integer*4 m,n
      real*8 k,z
      complex*16 hk(m,m)
c
	complex*16 ck,c2x,cem,cch,csh
c
	ck=dcmplx(k,0.d0)
	c2x=dcmplx(2.d0*k*z,0.d0)
c
	if(m.eq.2)then
c
c	  haskell propagator matrix for SH waves
c
	if(z.gt.0.d0)then
	  cem=cdexp(-c2x)
	  cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
	  csh=(0.5d0,0.d0)*((1.d0,0.d0)-cem)
	else
	  cem=cdexp(c2x)
	  cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
	  csh=-(0.5d0,0.d0)*((1.d0,0.d0)-cem)
	endif
c
c	propagator matrix for SH waves
c
	hk(1,1)=cch
	hk(1,2)=csh/(cmu(n)*ck)
	hk(2,1)=csh*cmu(n)*ck
	hk(2,2)=cch
	else
	  print *,'error in peghask: m schould be 2!'
	  return
	endif
c
	return
	end
