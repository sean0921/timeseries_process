c-----------------------------------------------------------------------

      subroutine lsqr(n,n1,g,d,m,w,s,sg,mis,rr)

      implicit none
      integer n,n1,i
      real*8 g(n,n1),d(n),gt(n1,n),gtg(n1,n1),invg(n1,n1),m(n1),rr
      real*8 w(n),nw(n,n),mg(n),avey,sg(n1),covm(n1,n1),sse,sst,ssr,mis
      real*8 s(n),covd(n,n)

      call diagm(n,w,nw)
      call diagm(n,s,covd)
      gt=transpose(g)
      gtg=matmul(matmul(gt,nw),g)
      call invert_matrix(gtg,invg,n1)
      m=matmul(matmul(matmul(invg,gt),nw),d)
      sse=0.  ! calculate "SSE"
c      mg=matmul(nw,d)-matmul(nw,matmul(g,m))
c      sig=sse/real(sum(nw)*(n-(n1-1))) ! Be careful!  Maybe mis=sqrt(sig)
      covm=matmul(matmul(matmul(matmul(invg,gt),nw),covd),
     +     transpose(matmul(matmul(invg,gt),nw)))
      do i=1,n1
        sg(i)=sqrt(covm(i,i))
      end do

      mg=d-matmul(g,m)
      do i=1,n
        sse=sse+mg(i)**2
      end do

      mis=sqrt(sse/(real(n)-n1))
      avey=sum(matmul(nw,d))/sum(nw)
      sst=0.  ! calculate "SST"
      do i=1,n
        sst=sst+(d(i)*w(i)-avey)**2
      end do
      mg=matmul(nw,d)-matmul(matmul(nw,g),m)
      ssr=0.  ! calculate "SSR"
      do i=1,n
        ssr=ssr+(d(i)*w(i)-avey)**2-mg(i)**2
      end do

      ! Calculate correlation coefficient
      rr=ssr/sst

      return
      end

c-----------------------------------------------------------------------

      subroutine invert_matrix(ma,inv,n) !求反矩陣

      integer i,j,n
      real*8 ma(n,n),temp(n,n),inv(n,n)

      do i=1,n
        do j=1,n
          temp(i,j)=ma(i,j)
          inv(i,j)=0.0
        end do
        inv(i,i)=1.
      end do

      call upper(temp,inv,n)  ! 做上三角矩陣處理
      call lower(temp,inv,n)  ! 做下三角矩陣處理

      do i=1,n
        do j=1,n
          inv(i,j)=inv(i,j)/temp(i,i)  ! 除上對角線元素
        end do
      end do

      return
      end

c---------------------------------

      subroutine upper(m,s,n)

      integer i,j,k,n
      real*8 e,m(n,n),s(n,n)

      do i=1,n-1
        do j=i+1,n
          e=m(j,i)/m(i,i)
          do k=1,n
            m(j,k)=m(j,k)-m(i,k)*e
            s(j,k)=s(j,k)-s(i,k)*e
          end do
        end do
      end do

      return
      end

c---------------------------------
      subroutine lower(m,s,n)

      integer i,j,k,n
      real*8 e,m(n,n),s(n,n)

      do i=n,2,-1
        do j=i-1,1,-1
          e=m(j,i)/m(i,i)
          do k=1,n
            m(j,k)=m(j,k)-m(i,k)*e
            s(j,k)=s(j,k)-s(i,k)*e
          end do
        end do
      end do

      return
      end

c-----------------------------------------------------------------------

      subroutine diagm(n,a,b)

      implicit none
      integer n,i,j
      real*8 a(n),b(n,n)

      do i=1,n
        do j=1,n
          b(i,j)=0.
        end do
        b(i,i)=a(i)
      end do

      return
      end
