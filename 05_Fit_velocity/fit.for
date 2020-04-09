      program fit

      implicit none
      integer i,j,s1,s2,m1,m2,m3,n,ns,stat,H,neq
      integer tyeq,tm,tm2,tm3,ty,statt,yy,day
      integer,allocatable::nn(:),tn(:)
      character infile*15,out*15,line*150,fmt*30,cm1*2,staa*4
      character,allocatable::sym(:)*1,sta(:)*4
      real*8 mis,rr,ve,vn,vh,se,sn,sh,lonn,latt,rt,c
      double precision,allocatable::t(:),g(:,:),d(:),w(:),m(:),sg(:),teq(:),ce(:),cn(:),ch(:),ces(:),cns(:),chs(:),pve(:),pvn(:),pvh(:),pves(:),pvns(:),pvhs(:),std(:),cal(:),t1(:),tf(:),dt(:),vee(:),vnn(:),vhh(:),see(:),snn(:),shh(:)

      open(1,file='fit.inp',status='old')
      ns=0
      stat=0
      do while (stat==0)
        read(1,*,iostat=stat)
        if(stat/=0)exit
        ns=ns+1
      end do
      close(1)
      allocate(sta(ns),tn(ns),t1(ns),tf(ns),dt(ns),nn(ns))
      allocate(vee(ns),vnn(ns),vhh(ns),see(ns),snn(ns),shh(ns))

      open(1,file='fit.para')
      open(2,file='fit.out')
      open(3,file='fit.inp',status='old')
      do s1=1,ns
        m1=2
        m2=0
        m3=0
        tyeq=0
        read(3,'(a)')line
        read(line,*)sta(s1),neq
        if(neq>0) then
          allocate(teq(neq),sym(neq))
          tyeq=1
          read(line,*)sta(s1),neq,(teq(i),sym(i),i=1,neq)
          do i=1,neq
            if(sym(i)=='C') then
              m1=m1+1
              m2=m2+1
            else if(sym(i)=='P') then
              m1=m1+1
            else if(sym(i)=='L') then
              m1=m1+2
              m2=m2+1
              m3=m3+1
            else if(sym(i)=='V') then
              m1=m1+1
              m3=m3+1
            else if(sym(i)=='T') then
              m1=m1+3
              m2=m2+1
              m3=m3+1
            end if
          end do
        end if
        allocate(ce(m2),cn(m2),ch(m2),ces(m2),cns(m2),chs(m2),pve(m3),pvn(m3),pvh(m3),pves(m3),pvns(m3),pvhs(m3))
        do s2=1,3
          if(s2==1) infile='ts_'//sta(s1)//'_e_o.dat'
          if(s2==2) infile='ts_'//sta(s1)//'_n_o.dat'
          if(s2==3) infile='ts_'//sta(s1)//'_u_o.dat'
          print*,'  Processing file: ',infile
          open(10,file=infile,status='old')
          n=0
          stat=0
          do while (stat==0)
            read(10,*,iostat=stat)
            if(stat/=0) exit
            n=n+1
          end do
          rewind(10)
          allocate(m(m1),sg(m1),t(n),g(n,m1),d(n),w(n),std(n),cal(n))
          do i=1,n
            read(10,*)t(i),d(i),std(i)
            tm=2
            g(i,1)=1
            g(i,2)=t(i)
            if(neq>0) then
              do j=1,neq
                if(t(i)>teq(j)) H=1
                if(t(i)<=teq(j))  H=0
                if(sym(j)=='C') then
                  tm=tm+1
                  g(i,tm)=H
                else if(sym(j)=='P') then
                  tm=tm+1
                  if(H==0) then
                    g(i,tm)=0.
                  else if(H==1.and.abs(t(i)-teq(j))>0.001) then
                    g(i,tm)=H*(log10(1+(t(i)-teq(j))))
                  end if
                else if(sym(j)=='L') then
                  tm=tm+1
                  g(i,tm)=H
                  tm=tm+1
                  g(i,tm)=H*(t(i)-teq(j))
                else if(sym(j)=='V') then
                  tm=tm+1
                  g(i,tm)=H*(t(i)-teq(j))
                else if(sym(j)=='T') then
                  tm=tm+1
                  g(i,tm)=H
                  tm=tm+1
                  g(i,tm)=H*(t(i)-teq(j))
                  tm=tm+1
                  if(H==0) then
                    g(i,tm)=0.
                  else if(H==1.and.abs(t(i)-teq(j))>0.001) then
                    g(i,tm)=H*(log10(1+(t(i)-teq(j))))
                  end if
                end if
              end do
            end if
            w(i)=1./real(n)
          end do
          close(10)
          call lsqr(n,m1,g,d,m,w,std,sg,mis,rr)
          if(s2==1) out='ts_'//sta(s1)//'_e_c.dat'
          if(s2==2) out='ts_'//sta(s1)//'_n_c.dat'
          if(s2==3) out='ts_'//sta(s1)//'_u_c.dat'
          do i=1,n
            cal(i)=0
            do j=1,m1
              cal(i)=cal(i)+g(i,j)*m(j)
            end do
          end do
          open(10,file=out)
          do i=1,n
            write(10,'(f10.5,f10.2)')t(i),cal(i)
          end do
          close(10)

          fmt='(a4,1x,a1,1x,10f19.9)'
          call int2char(m1,cm1)
          fmt(14:15)=cm1(1:2)
          if(s2==1) write(1,fmt)sta(s1),'e',(m(i),i=1,m1)
          if(s2==2) write(1,fmt)sta(s1),'n',(m(i),i=1,m1)
          if(s2==3) write(1,fmt)sta(s1),'u',(m(i),i=1,m1)
          fmt='(7x,10f19.9)'
          fmt(5:6)=cm1(1:2)
          write(1,fmt)(sg(i),i=1,m1)

          if(s2==1) then
            tn(s1)=n
            t1(s1)=t(1)
            tf(s1)=t(n)
            dt(s1)=t(n)-t(1)
            nn(s1)=1
            do i=2,n
              if((t(i)-t(i-1))>0.5) nn(s1)=nn(s1)+1
c              if((t(i)-t(i-1))<0.5) t(i)=t(i-1)
            end do
          end if

          ! strat to create trend of coordinate time-series
          if(s2==1) out='ts_'//sta(s1)//'_e_t.dat'
          if(s2==2) out='ts_'//sta(s1)//'_n_t.dat'
          if(s2==3) out='ts_'//sta(s1)//'_u_t.dat'
          open(11,file=out)
          yy=int(t(1))
          if(mod(yy,4)==0) day=nint((t(1)-real(yy))*366.+1.)-1
          if(mod(yy,4)/=0) day=nint((t(1)-real(yy))*365.+1.)-1
          do while (.true.)
            day=day+1                                 ! generate time table
            if(mod(yy,4)==0.and.day==367) then
              yy=yy+1
              day=day-366
            else if(mod(yy,4)/=0.and.day==366) then
              yy=yy+1
              day=day-365
            end if
            if(mod(yy,4)==0) rt=real(yy)+((real(day)-1.0)/366.0)
            if(mod(yy,4)/=0) rt=real(yy)+((real(day)-1.0)/365.0)
c            if(abs(rt-t(n))<.001) exit
            if(rt>t(n)) exit
            if(neq==0) c=m(1)+rt*m(2)
            if(neq>0) then
              c=m(1)+rt*m(2)
              tm=2
              do j=1,neq
                if(rt>teq(j).and.sym(j)=='C') then
                  tm=tm+1
                  c=c+m(tm)
                else if(rt>teq(j).and.sym(j)=='P') then
                  tm=tm+1
                  c=c+m(tm)*(log10(1+(rt-teq(j))))
                else if(rt>teq(j).and.sym(j)=='L') then
                  tm=tm+2
                  c=c+m(tm-1)+m(tm)*(rt-teq(j))
                else if(rt>teq(j).and.sym(j)=='V') then
                  tm=tm+1
                  c=c+m(tm)*(rt-teq(j))
                else if(rt>teq(j).and.sym(j)=='T') then
                  tm=tm+3
                  c=c+m(tm-2)+m(tm-1)*(rt-teq(j))+m(tm)*(log10(1+(rt-teq(j))))
                end if
              end do
            end if
            write(11,'(f9.4,f13.5)')rt,c
          end do
          close(11)                                   ! end of creating coordinate time-series

          if(s2==1) then
            ve=m(2)
            vee(s1)=m(2)
            se=sqrt(sg(2)*sg(2)*(mis/2.d0)*(mis/2.d0))
            see(s1)=se
            tm=2
            tm2=0
            tm3=0
            if(neq>0) then
              do i=1,neq
                if(sym(i)=='C') then
                  tm=tm+1
                  tm2=tm2+1
                  ce(tm2)=m(tm)
                  ces(tm2)=sg(tm)
                else if(sym(i)=='L') then
                  tm=tm+2
                  tm2=tm2+1
                  tm3=tm3+1
                  ce(tm2)=m(tm-1)
                  ces(tm2)=sg(tm-1)
                  if(tm3==1) then
                    pve(tm3)=m(tm)+m(2)
                    pves(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                  else
                    pve(tm3)=m(tm)+pve(tm3-1)
                    pves(tm3)=sqrt(sg(tm)**2+pves(tm3-1)**2)
                  end if
                else if(sym(i)=='V') then
                  tm=tm+1
                  tm3=tm3+1
                  if(tm3==1) then
                    pve(tm3)=m(tm)+m(2)
                    pves(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                  else
                    pve(tm3)=m(tm)+pve(tm3-1)
                    pves(tm3)=sqrt(sg(tm)**2+pves(tm3-1)**2)
                  end if
                else if(sym(i)=='T') then
                  tm=tm+3
                  tm2=tm2+1
                  tm3=tm3+1
                  ce(tm2)=m(tm-2)+m(tm)*log10(0.001)
                  ces(tm2)=sg(tm-2)
                  if(tm3==1) then
                    pve(tm3)=m(tm-1)+ve
                    pves(tm3)=sqrt(sg(tm-1)**2+se**2)
                  else
                    pve(tm3)=m(tm-1)+pve(tm3-1)
                    pves(tm3)=sqrt(sg(tm-1)**2+pves(tm3-1)**2)
                  end if
                end if
              end do
            end if
          else if(s2==2) then
            vn=m(2)
            vnn(s1)=m(2)
            sn=sqrt(sg(2)*sg(2)*(mis/2.d0)*(mis/2.d0))
            snn(s1)=sn
            tm=2
            tm2=0
            tm3=0
            if(neq>0) then
              do i=1,neq
                if(sym(i)=='C') then
                  tm=tm+1
                  tm2=tm2+1
                  cn(tm2)=m(tm)
                  cns(tm2)=sg(tm)
                else if(sym(i)=='L') then
                  tm=tm+2
                  tm2=tm2+1
                  tm3=tm3+1
                  cn(tm2)=m(tm-1)
                  cns(tm2)=sg(tm-1)
                  if(tm3==1) then
                    pvn(tm3)=m(tm)+m(2)
                    pvns(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                  else
                    pvn(tm3)=m(tm)+pvn(tm3-1)
                    pvns(tm3)=sqrt(sg(tm)**2+pvns(tm3-1)**2)
                  end if
                else if(sym(i)=='V') then
                  tm=tm+1
                  tm3=tm3+1
                  if(tm3==1) then
                    pvn(tm3)=m(tm)+m(2)
                    pvns(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                  else
                    pvn(tm3)=m(tm)+pvn(tm3-1)
                    pvns(tm3)=sqrt(sg(tm)**2+pvns(tm3-1)**2)
                  end if
                else if(sym(i)=='T') then
                  tm=tm+3
                  tm2=tm2+1
                  tm3=tm3+1
                  cn(tm2)=m(tm-2)+m(tm)*log10(0.001)
                  cns(tm2)=sg(tm-2)
                  if(tm3==1) then
                    pvn(tm3)=m(tm-1)+vn
                    pvns(tm3)=sqrt(sg(tm-1)**2+sn**2)
                  else
                    pvn(tm3)=m(tm-1)+pvn(tm3-1)
                    pvns(tm3)=sqrt(sg(tm-1)**2+pvns(tm3-1)**2)
                  end if
                end if
              end do
            end if
          else if(s2==3) then
            vh=m(2)
            vhh(s1)=m(2)
            sh=sqrt(sg(2)*sg(2)*(mis/2.d0)*(mis/2.d0))
            shh(s1)=sh
            tm=2
            tm2=0
            tm3=0
            if(neq>0) then
              do i=1,neq
                if(sym(i)=='C') then
                  tm=tm+1
                  tm2=tm2+1
                  ch(tm2)=m(tm)
                  chs(tm2)=sg(tm)
                else if(sym(i)=='L') then
                  tm=tm+2
                  tm2=tm2+1
                  tm3=tm3+1
                  ch(tm2)=m(tm-1)
                  chs(tm2)=sg(tm-1)
                  if(tm3==1) then
                    pvh(tm3)=m(tm)+m(2)
                    pvhs(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                  else
                    pvh(tm3)=m(tm)+pvh(tm3-1)
                    pvhs(tm3)=sqrt(sg(tm)**2+pvhs(tm3-1)**2)
                  end if
                else if(sym(i)=='V') then
                  tm=tm+1
                  tm3=tm3+1
                  if(tm3==1) then
                    pvh(tm3)=m(tm)+m(2)
                    pvhs(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                  else
                    pvh(tm3)=m(tm)+pvh(tm3-1)
                    pvhs(tm3)=sqrt(sg(tm)**2+pvhs(tm3-1)**2)
                  end if
                else if(sym(i)=='T') then
                  tm=tm+3
                  tm2=tm2+1
                  tm3=tm3+1
                  ch(tm2)=m(tm-2)+m(tm)*log10(0.001)
                  chs(tm2)=sg(tm-2)
                  if(tm3==1) then
                    pvh(tm3)=m(tm-1)+vh
                    pvhs(tm3)=sqrt(sg(tm-1)**2+sh**2)
                  else
                    pvh(tm3)=m(tm-1)+pvh(tm3-1)
                    pvhs(tm3)=sqrt(sg(tm-1)**2+pvhs(tm3-1)**2)
                  end if
                end if
              end do
            end if
          end if
          deallocate(m,sg,t,g,d,w,std,cal)
        end do
        ty=neq
        if(neq==0) ty=neq
        if(neq>0) then
          do i=1,neq
            if(sym(i)=='P') ty=ty-1
          end do
        end if
        write(2,'("# ",a4,"  Number of Events: ",i3)')sta(s1),ty
        write(2,'("    secular motion (mm/yr):")')
        write(2,'("      E component: ",f8.3," +-",f8.3)')ve,se
        write(2,'("      N component: ",f8.3," +-",f8.3)')vn,sn
        write(2,'("      U component: ",f8.3," +-",f8.3)')vh,sh
        if(neq>0) then
          tm2=0
          tm3=0
          do i=1,neq
            if(sym(i)=='C') then
              tm2=tm2+1
              write(2,'("  ! earthquake occurred time (yr):",f11.4)')teq(i)
              write(2,'("    coseismic displacements (mm):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')ce(tm2),ces(tm2)
              write(2,'("      N component: ",f8.3," +-",f8.3)')cn(tm2),cns(tm2)
              write(2,'("      U component: ",f8.3," +-",f8.3)')ch(tm2),chs(tm2)
            else if(sym(i)=='L') then
              tm2=tm2+1
              tm3=tm3+1
              write(2,'("  ! earthquake occurred time (yr):",f11.4)')teq(i)
              write(2,'("    coseismic displacements (mm):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')ce(tm2),ces(tm2)
              write(2,'("      N component: ",f8.3," +-",f8.3)')cn(tm2),cns(tm2)
              write(2,'("      U component: ",f8.3," +-",f8.3)')ch(tm2),chs(tm2)
              write(2,'("    linear velocity (mm/yr):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')pve(tm3),pves(tm3)
              write(2,'("      N component: ",f8.3," +-",f8.3)')pvn(tm3),pvns(tm3)
              write(2,'("      U component: ",f8.3," +-",f8.3)')pvh(tm3),pvhs(tm3)
            else if(sym(i)=='V') then
              tm3=tm3+1
              write(2,'("  ! velocity change time (yr):",f11.4)')teq(i)
              write(2,'("    linear velocity (mm/yr):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')pve(tm3),pves(tm3)
              write(2,'("      N component: ",f8.3," +-",f8.3)')pvn(tm3),pvns(tm3)
              write(2,'("      U component: ",f8.3," +-",f8.3)')pvh(tm3),pvhs(tm3)
            else if(sym(i)=='T') then
              tm2=tm2+1
              tm3=tm3+1
              write(2,'("  ! earthquake occurred time (yr):",f11.4)')teq(i)
              write(2,'("    coseismic displacements (mm):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')ce(tm2),ces(tm2)
              write(2,'("      N component: ",f8.3," +-",f8.3)')cn(tm2),cns(tm2)
              write(2,'("      U component: ",f8.3," +-",f8.3)')ch(tm2),chs(tm2)
              write(2,'("    linear velocity (mm/yr):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')pve(tm3),pves(tm3)
              write(2,'("      N component: ",f8.3," +-",f8.3)')pvn(tm3),pvns(tm3)
              write(2,'("      U component: ",f8.3," +-",f8.3)')pvh(tm3),pvhs(tm3)
           end if
         end do
        end if

        if(tyeq==1) deallocate(teq,sym)
        deallocate(ce,cn,ch,ces,cns,chs,pve,pvn,pvh,pves,pvns,pvhs)
      end do
      close(3)
      close(2)
      close(1)

      open(99,file='fit.fout')
      do i=1,ns
c        if(nn(i)<2) cycle
c        if(dt(i)<2) cycle
        open(98,file='getpl.gout',status='old')
        statt=0
        do while (statt==0)
          if(statt/=0) exit
          read(98,*,iostat=statt)staa,lonn,latt
          if(sta(i).eq.staa) then
             write(99,'(a5,2x,2f8.4,6f10.2,2i5,3f10.4)')sta(i),lonn,latt,vee(i),vnn(i),vhh(i),see(i),snn(i),shh(i),nn(i),tn(i),dt(i),t1(i),tf(i)
             exit
          end if
        end do
       rewind(98)
      end do
      close(99)
      close(98)

      stop
      end

c.......................................................................
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

c.......................................................................
      subroutine invert_matrix(ma,inv,n) ! 求反矩陣

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

c.......................................................................
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

c***********************************************************************
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
c***********************************************************************
      subroutine int2char(a,b)
c-------------------------------------------------------------
c     subroutine int2char : change integer to character (1-40)
c-------------------------------------------------------------
      integer a
      character b*2

      if (a==0)  b='00'
      if (a==1)  b='01'
      if (a==2)  b='02'
      if (a==3)  b='03'
      if (a==4)  b='04'
      if (a==5)  b='05'
      if (a==6)  b='06'
      if (a==7)  b='07'
      if (a==8)  b='08'
      if (a==9)  b='09'
      if (a==10) b='10'
      if (a==11) b='11'
      if (a==12) b='12'
      if (a==13) b='13'
      if (a==14) b='14'
      if (a==15) b='15'
      if (a==16) b='16'
      if (a==17) b='17'
      if (a==18) b='18'
      if (a==19) b='19'
      if (a==20) b='20'
      if (a==21) b='21'
      if (a==22) b='22'
      if (a==23) b='23'
      if (a==24) b='24'
      if (a==25) b='25'
      if (a==26) b='26'
      if (a==27) b='27'
      if (a==28) b='28'
      if (a==29) b='29'
      if (a==30) b='30'
      if (a==31) b='31'
      if (a==32) b='32'
      if (a==33) b='33'
      if (a==34) b='34'
      if (a==35) b='35'
      if (a==36) b='36'
      if (a==37) b='37'
      if (a==38) b='38'
      if (a==39) b='39'
      if (a==40) b='40'
      if (a==41) b='41'
      if (a==42) b='42'
      if (a==43) b='43'
      if (a==44) b='44'
      if (a==45) b='45'
      if (a==46) b='46'
      if (a==47) b='47'
      if (a==48) b='48'
      if (a==49) b='49'
      if (a==50) b='50'
      if (a==51) b='51'
      if (a==52) b='52'
      if (a==53) b='53'
      if (a==54) b='54'
      if (a==55) b='55'
      if (a==56) b='56'
      if (a==57) b='57'
      if (a==58) b='58'
      if (a==59) b='59'
      if (a==60) b='60'
      if (a==61) b='61'
      if (a==62) b='62'
      if (a==63) b='63'
      if (a==64) b='64'
      if (a==65) b='65'
      if (a==66) b='66'
      if (a==67) b='67'
      if (a==68) b='68'
      if (a==69) b='69'
      if (a==70) b='70'
      if (a==71) b='71'
      if (a==72) b='72'
      if (a==73) b='73'
      if (a==74) b='74'
      if (a==75) b='75'
      if (a==76) b='76'
      if (a==77) b='77'
      if (a==78) b='78'
      if (a==79) b='79'
      if (a==80) b='80'
      if (a==81) b='81'
      if (a==82) b='82'
      if (a==83) b='83'
      if (a==84) b='84'
      if (a==85) b='85'
      if (a==86) b='86'
      if (a==87) b='87'
      if (a==88) b='88'
      if (a==89) b='89'
      if (a==90) b='90'
      if (a==91) b='91'
      if (a==92) b='92'
      if (a==93) b='93'
      if (a==94) b='94'
      if (a==95) b='95'
      if (a==96) b='96'
      if (a==97) b='97'
      if (a==98) b='98'
      if (a==99) b='99'
      return
      end
