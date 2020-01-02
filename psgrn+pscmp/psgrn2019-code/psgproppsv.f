      subroutine psgproppsv(ma,mai,nma,l1,l2,ypsv,ypsv0,k)
      use psgalloc
      implicit none
c
c     propagation of p-sv vectors
c
      integer*4 l1,l2,nma
      real*8 k
      complex*16 ypsv(6,3),ypsv0(6,3)
      complex*16 ma(6,6,nma),mai(6,6,nma)
c
c     work space
c
      integer*4 i,j,l
      real*8 h0
      complex*16 cdet,wave,wave2
      complex*16 c0(6,3),c1(6,3),y1(6,3),orth(3,3)
c
      if(l1.eq.l2)then
        return
      else if(l1.lt.l2)then
        do l=l1+1,l2
          h0=hp(l-1)
          wave=cdexp(dcmplx(-k*h0,0.d0))
          wave2=wave*wave
c
          call caxcb(mai(1,1,l-1),ypsv,6,6,3,c0)
c
c         orthonormalization of the p-sv modes
c
          cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &        +c0(3,1)*c0(5,2)*c0(1,3)
     &        +c0(5,1)*c0(1,2)*c0(3,3)
     &        -c0(5,1)*c0(3,2)*c0(1,3)
     &        -c0(3,1)*c0(1,2)*c0(5,3)
     &        -c0(1,1)*c0(5,2)*c0(3,3)
          orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
          orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
          orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
          orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
          orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
          orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
          orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
          orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
          orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
          call caxcb(c0,orth,6,3,3,c1)
          if(l.gt.lzrec)then
c
c           orthonormalization of the receiver vectors
c
            call caxcb(ypsv0,orth,6,3,3,y1)
            call cmemcpy(y1,ypsv0,18)
            do j=1,3
              do i=1,6
                ypsv0(i,j)=ypsv0(i,j)*wave
              enddo
            enddo
          endif
          c1(1,1)=(1.d0,0.d0)
          c1(2,1)=c1(2,1)*wave2
          c1(3,1)=(0.d0,0.d0)
          c1(4,1)=c1(4,1)*wave2
          c1(5,1)=(0.d0,0.d0)
          c1(6,1)=c1(6,1)*wave2
c
          c1(1,2)=(0.d0,0.d0)
          c1(2,2)=c1(2,2)*wave2
          c1(3,2)=(1.d0,0.d0)
          c1(4,2)=c1(4,2)*wave2
          c1(5,2)=(0.d0,0.d0)
          c1(6,2)=c1(6,2)*wave2
c
          c1(1,3)=(0.d0,0.d0)
          c1(2,3)=c1(2,3)*wave2
          c1(3,3)=(0.d0,0.d0)
          c1(4,3)=c1(4,3)*wave2
          c1(5,3)=(1.d0,0.d0)
          c1(6,3)=c1(6,3)*wave2
c
          call caxcb(ma(1,1,l-1),c1,6,6,3,ypsv)
          if(l.eq.lzrec)call cmemcpy(ypsv,ypsv0,18)
        enddo
      else
        do l=l1-1,l2,-1
          h0=hp(l)
          wave=cdexp(dcmplx(-k*h0,0.d0))
          wave2=wave*wave
c
          call caxcb(mai(1,1,l),ypsv,6,6,3,c0)
c
c         orthonormalization of the p-sv modes
c
          cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &        +c0(4,1)*c0(6,2)*c0(2,3)
     &        +c0(6,1)*c0(2,2)*c0(4,3)
     &        -c0(6,1)*c0(4,2)*c0(2,3)
     &        -c0(4,1)*c0(2,2)*c0(6,3)
     &        -c0(2,1)*c0(6,2)*c0(4,3)
          orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
          orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
          orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
          orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
          orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
          orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
          orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
          orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
          orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
c
          call caxcb(c0,orth,6,3,3,c1)
          if(l.lt.lzrec)then
c
c           orthonormalization of the receiver vectors
c
            call caxcb(ypsv0,orth,6,3,3,y1)
            call cmemcpy(y1,ypsv0,18)
            do j=1,3
              do i=1,6
                ypsv0(i,j)=ypsv0(i,j)*wave
              enddo
            enddo
          endif
c
	    c1(1,1)=c1(1,1)*wave2
 	    c1(2,1)=(1.d0,0.d0)
	    c1(3,1)=c1(3,1)*wave2
          c1(4,1)=(0.d0,0.d0)
	    c1(5,1)=c1(5,1)*wave2
          c1(6,1)=(0.d0,0.d0)
c
          c1(1,2)=c1(1,2)*wave2
          c1(2,2)=(0.d0,0.d0)
          c1(3,2)=c1(3,2)*wave2
          c1(4,2)=(1.d0,0.d0)
          c1(5,2)=c1(5,2)*wave2
          c1(6,2)=(0.d0,0.d0)
c
          c1(1,3)=c1(1,3)*wave2
          c1(2,3)=(0.d0,0.d0)
          c1(3,3)=c1(3,3)*wave2
          c1(4,3)=(0.d0,0.d0)
          c1(5,3)=c1(5,3)*wave2
          c1(6,3)=(1.d0,0.d0)
c
          call caxcb(ma(1,1,l),c1,6,6,3,ypsv)
          if(l.eq.lzrec)call cmemcpy(ypsv,ypsv0,18)
        enddo
      endif
      return
      end
