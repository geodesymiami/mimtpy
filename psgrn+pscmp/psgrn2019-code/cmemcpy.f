        subroutine cmemcpy(a,b,n)
        implicit none
c
c       copy complex array a to b
c
        integer*4 n
        complex*16 a(n),b(n)
c
        integer*4 i
c
        do i=1,n
          b(i)=a(i)
        enddo
        return
        end
