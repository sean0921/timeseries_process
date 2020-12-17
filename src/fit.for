      program fit

      implicit none
      integer i,j,s1,s2,m1,m2,m3,n,ns,stat,H,neq
      integer tyeq,tm,tm2,tm3,ty,statt,yy,day
      integer,allocatable::nn(:),tn(:)
      character infile*15,out*15,line*150,fmt*30,cm1*2,staa*4
      character,allocatable::sym(:)*1,sta(:)*4
      real*8 mis,rr,ve,vn,vh,se,sn,sh,lonn,latt,rt,c
      double precision,allocatable::t(:),g(:,:),d(:),w(:),m(:),sg(:),teq(:),ce(:),cn(:),ch(:),ces(:),cns(:),chs(:),pve(:),pvn(:),pvh(:),pves(:),pvns(:),pvhs(:),std(:),cal(:),t1(:),tf(:),dt(:),vee(:),vnn(:),vhh(:),see(:),snn(:),shh(:)

      print*,'                           PROGRAM FIT                                 '
      print*,'                                            by Kuo-En Ching 2005.01.11 '
      print*,'Seismology Lab, Department of Earth Sciences, NCKU   Update 2020.04.08+'
      print*,''
      print*,'source code: https://github.com/sean0921/timeseries_process'
      print*,''

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
          if(s2==1) infile='ts_'//sta(s1)//'_e_b.dat'
          if(s2==2) infile='ts_'//sta(s1)//'_n_b.dat'
          if(s2==3) infile='ts_'//sta(s1)//'_u_b.dat'
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
ccccccc  TODO: Figure out why this line are added, and its original usage
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
ccccccc  TODO: Figure out why those 2 lines are added, and those original usage
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
