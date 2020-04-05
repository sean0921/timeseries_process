c --------------------- Common Library for TimeSeries Processing Tools ---------------------

c-----------------------------------------------------------------------

      subroutine bubble_sort(a,n)
      implicit none
      integer i,j,n
      real*8 a(n),tmp1

      do i=n-1,1,-1
        do j=1,i
          if(a(j).gt.a(j+1)) then
            tmp1=a(j)
            a(j)=a(j+1)
            a(j+1)=tmp1
          end if
        end do
      end do

      return
      end

c--------------------------------------------------------------------

      subroutine bubble_sort_3(a,b,c,n)
      implicit none
      integer n
      real*8 b(n),c(n),temp2,temp3
      character a(n)*4,temp1*4
      integer i,j

      do i=n-1,1,-1
        do j=1,i
          if(a(j).gt.a(j+1)) then
            temp1=a(j)
            temp2=b(j)
            temp3=c(j)
            a(j)=a(j+1)
            b(j)=b(j+1)
            c(j)=c(j+1)
            a(j+1)=temp1
            b(j+1)=temp2
            c(j+1)=temp3
          end if
        end do
      end do

      return
      end

c--------------------------------------------------------------------

      subroutine bubble_sort_10(a,b,c,d,e,f,g,h,k,l,n)
      implicit none
      integer n
      real*8 a(n),d(n),e(n),f(n),g(n),h(n),k(n),l(n)
      real*8 tmp1,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10
      integer i,j,b(n),c(n),tmp2,tmp3

      do i=n-1,1,-1
        do j=1,i
          if(a(j).gt.a(j+1)) then
            tmp1=a(j)
            tmp2=b(j)
            tmp3=c(j)
            tmp4=d(j)
            tmp5=e(j)
            tmp6=f(j)
            tmp7=g(j)
            tmp8=h(j)
            tmp9=k(j)
            tmp10=l(j)
            a(j)=a(j+1)
            b(j)=b(j+1)
            c(j)=c(j+1)
            d(j)=d(j+1)
            e(j)=e(j+1)
            f(j)=f(j+1)
            g(j)=g(j+1)
            h(j)=h(j+1)
            k(j)=k(j+1)
            l(j)=l(j+1)
            a(j+1)=tmp1
            b(j+1)=tmp2
            c(j+1)=tmp3
            d(j+1)=tmp4
            e(j+1)=tmp5
            f(j+1)=tmp6
            g(j+1)=tmp7
            h(j+1)=tmp8
            k(j+1)=tmp9
            l(j+1)=tmp10
          end if
        end do
      end do

      return
      end
