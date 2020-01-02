      subroutine psgkern(y,k,clahs,cmuhs)
      use psgalloc
      implicit none
c
c     calculation of response function in Laplace domain
c     y(8,4): solution vector (complex)
c     k: wave number (input)
c
      real*8 k
      complex*16 clahs,cmuhs,y(8,4)
c
      integer*4 i,istp
      complex*16 yhs(8,4)
c
      real*8 eps
      data eps/1.0d-06/
c
      do istp=1,4
        do i=1,8
          y(i,istp)=(0.d0,0.d0)
          yhs(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      nongravity=kgmax*grfac.lt.eps*k
c
      call psgpsv(y,k)
      call psgsh(y,k)
c
c     subtract the halfspace solution
c
      call psghskern(yhs,k,clahs,cmuhs)
c
      do istp=1,4
        do i=1,8
          y(i,istp)=y(i,istp)-yhs(i,istp)
        enddo
      enddo
c
      return
      end	  
