c--     by Kuo-En Ching 2006.05.02
c--              Update 2006.05.07
c--     Considering that the post-seismic decay is LOG-type decay
c
c----comfilt.for---
c
      program comfilt

      implicit none
      character infile*15,outfile*15,outfil1*15,line*150,sta*4
      character,allocatable::sym(:)*1
      real*8 mis,rr,pi
      double precision,allocatable::t(:),std(:),g(:,:),d(:),w(:),cal(:)
      double precision,allocatable::m(:),sg(:),teq(:)
      integer i,j,s1,s2,m1,n,ns,stat,neq,H,tyeq,tm

      pi=atan(1.)*4

      open(1,file='comfilt.inp',status='old')
      ns=0
      stat=0
      do while (stat==0)
        read(1,*,iostat=stat)
        if(stat/=0)exit
        ns=ns+1
      end do
      rewind(1)
      do s1=1,ns
        m1=6
        tyeq=0
        read(1,'(a)')line
        read(line,*)sta,neq
        if(neq>0) then
          allocate(teq(neq),sym(neq))
          tyeq=1
          read(line,*)sta,neq,(teq(i),sym(i),i=1,neq)
          do i=1,neq
            if(sym(i)=='C') m1=m1+1
            if(sym(i)=='P') m1=m1+1
            if(sym(i)=='L') m1=m1+2
            if(sym(i)=='T') m1=m1+3
          end do
        end if
        do s2=1,3
          if(s2==1) infile='ts_'//sta//'_e_b.dat'
          if(s2==2) infile='ts_'//sta//'_n_b.dat'
          if(s2==3) infile='ts_'//sta//'_u_b.dat'
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
          allocate(m(m1),sg(m1),t(n),std(n),g(n,m1),d(n),w(n),cal(n))
          do i=1,n
            read(10,*)t(i),d(i),std(i)
            tm=6
            g(i,1)=1
            g(i,2)=t(i)
            g(i,3)=sin(2*pi*t(i))
            g(i,4)=cos(2*pi*t(i))
            g(i,5)=sin(4*pi*t(i))
            g(i,6)=cos(4*pi*t(i))
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
                    g(i,tm)=H*(log10(t(i)-teq(j))-log10(0.00273))
                  end if
                else if(sym(j)=='L') then
                  tm=tm+1
                  g(i,tm)=H
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
                    g(i,tm)=H*(log10(t(i)-teq(j))-log10(0.00273))
                  end if
                end if
              end do
            end if
            w(i)=1./std(i)**2
          end do
          close(10)

          call lsqr(n,m1,g,d,m,w,std,sg,mis,rr)
          if(s2==1) outfile='ts_'//sta//'_e_r.dat'
          if(s2==1) outfil1='ts_'//sta//'_e_c.dat'
          if(s2==2) outfile='ts_'//sta//'_n_r.dat'
          if(s2==2) outfil1='ts_'//sta//'_n_c.dat'
          if(s2==3) outfile='ts_'//sta//'_u_r.dat'
          if(s2==3) outfil1='ts_'//sta//'_u_c.dat'
          open(10,file=outfile)
          open(11,file=outfil1)
          cal=d-matmul(g,m)
          do i=1,n
            if(sta/='ILAN'.and.dabs(std(i))>3*mis) cycle
            if(dabs(cal(i))>3*mis) cycle
            write(10,'(f10.5,f10.2)')t(i),cal(i)
            write(11,'(f10.5,2f10.2)')t(i),d(i),std(i)
          end do
          close(10)
          close(11)
          deallocate(m,sg,t,std,g,d,w,cal)
        end do
        if(tyeq==1) deallocate(teq,sym)
      end do
      close(1)

      call filter(ns)
      call fit(ns)

      stop
      end

c.......................................................................
      subroutine filter(ns)

      implicit none
      integer i,j,c,k,tn,ne,ns,nt,stat,t1,n(ns)
      character infile*15,outfile*8,outfile1*15,sta(ns)*4
      real*8 err,rt1,std
      double precision,allocatable::t(:,:),res(:,:),et(:),e(:),rt(:)
      double precision,allocatable::tt(:)
      logical alive

      print*,'Start filter....'

#ifdef MINGW
      inquire(file='warning-e.sum',exist=alive)
      if(alive) call system('del warning-e.sum')
      inquire(file='warning-n.sum',exist=alive)
      if(alive) call system('del warning-n.sum')
      inquire(file='warning-u.sum',exist=alive)
      if(alive) call system('del warning-u.sum')
      inquire(file='fil.dat',exist=alive)
      if(alive) call system('del fil.dat')
      inquire(file='tmp.dat',exist=alive)
      if(alive) call system('del fil.dat')
#else
      inquire(file='warning-e.sum',exist=alive)
      if(alive) call system('rm -f warning-e.sum')
      inquire(file='warning-n.sum',exist=alive)
      if(alive) call system('rm -f warning-n.sum')
      inquire(file='warning-u.sum',exist=alive)
      if(alive) call system('rm -f warning-u.sum')
      inquire(file='fil.dat',exist=alive)
      if(alive) call system('rm -f fil.dat')
      inquire(file='tmp.dat',exist=alive)
      if(alive) call system('rm -f fil.dat')
#endif

      do c=1,3
        open(10,file='comfilt.inp',status='old')
        open(11,file='fil.dat')
        do i=1,ns
          read(10,'(a4)')sta(i)
          if(c==1) infile='ts_'//sta(i)//'_e_r.dat'
          if(c==2) infile='ts_'//sta(i)//'_n_r.dat'
          if(c==3) infile='ts_'//sta(i)//'_u_r.dat'
          write(11,'(a15)')infile
        end do
        close(10)
        close(11)
        open(10,file='fil.dat',status='old')
        do i=1,ns
          read(10,'(a15)')infile
          read(infile,'(3x,a4)')sta(i)
          open(11,file=infile,status='old')
          n(i)=0
          stat=0
          do while (stat==0)
            read(11,*,iostat=stat)
            if(stat/=0) exit
            n(i)=n(i)+1
          end do
          close(11)
        end do
        tn=0
        do i=1,ns
          tn=tn+n(i)
        end do
        rewind(10)
        allocate(tt(tn))
        tn=0
        do i=1,ns
          read(10,'(a15)')infile
          open(11,file=infile,status='old')
          do j=1,n(i)
            tn=tn+1
            read(11,*)tt(tn)
          end do
          close(11)
        end do
        close(10,status='delete')

        call bubble_sort(tt,tn)
        t1=1
        do i=2,tn
          if(tt(i)-tt(i-1)>0.0001) t1=t1+1
        end do
        allocate(rt(t1))
        t1=1
        rt(t1)=tt(1)
        do i=2,tn
          if(tt(i)-tt(i-1)>.0001) then
            t1=t1+1
            rt(t1)=tt(i)
          end if
        end do

        allocate(t(ns,tn),res(ns,tn))
        do i=1,ns
          if(c==1) infile='ts_'//sta(i)//'_e_r.dat'
          if(c==2) infile='ts_'//sta(i)//'_n_r.dat'
          if(c==3) infile='ts_'//sta(i)//'_u_r.dat'
          open(10,file=infile,status='old')
          do j=1,n(i)
            read(10,*)t(i,j),res(i,j)
          end do
          close(10)
        end do

        if(c==1) outfile='come.gmt'  ! step 2: stacking
        if(c==2) outfile='comn.gmt'
        if(c==3) outfile='comu.gmt'
        open(11,file=outfile)
        ne=0
        do k=1,t1
          err=0.
          nt=0
          do i=1,ns
            do j=1,n(i)
              if(abs(rt(k)-t(i,j))<.0001) then
                err=err+res(i,j)
                nt=nt+1
                exit
              end if
            end do
          end do
          if(nt<4) then
            if(c==1) open(10,file='warning-e.sum',position='append')
            if(c==2) open(10,file='warning-n.sum',position='append')
            if(c==3) open(10,file='warning-h.sum',position='append')
            write(10,'(f10.5," Number of stations: ",i1," < 4 !!")')
     +      rt(k),nt
            close(10)
          end if
          if(nt<2) cycle
          err=err/real(nt)
          write(11,'(f10.5,f9.2)')rt(k),err
          ne=ne+1
        end do
        close(11)

        open(10,file=outfile,status='old')  ! step 3: filtering
        allocate(et(ne),e(ne))
        do i=1,ne
          read(10,*)et(i),e(i)
        end do
        close(10)
        do i=1,ns
          if(c==1) infile='ts_'//sta(i)//'_e_c.dat'
          if(c==1) outfile1='ts_'//sta(i)//'_e_f.dat'
          if(c==2) infile='ts_'//sta(i)//'_n_c.dat'
          if(c==2) outfile1='ts_'//sta(i)//'_n_f.dat'
          if(c==3) infile='ts_'//sta(i)//'_u_c.dat'
          if(c==3) outfile1='ts_'//sta(i)//'_u_f.dat'
          open(10,file=infile,status='old')
          open(11,file=outfile1)
          do j=1,n(i)
            read(10,*)rt1,err,std
            do k=1,ne
              if(abs(rt1-et(k))<.0001)then
                err=err-e(k)
                write(11,'(f10.5,2f10.2)')rt1,err,std
              end if
            end do
          end do
          close(10)
          close(11)
        end do
        deallocate(t,res,et,e,rt,tt)
      end do

      print*,'End filter....'

      return
      end

c.......................................................................
      subroutine fit(ns)

      implicit none
      integer i,j,s1,s2,m1,m2,m3,n,ns,stat,H,neq
      integer tyeq,tm,tm2,tm3,ty
      character infile*15,out*15,sta(ns)*4,line*150,fmt*30,cm1*2
      character,allocatable::sym(:)*1
      real*8 mis,rr,ve,vn,vh,se,sn,sh,pi
      double precision,allocatable::t(:),g(:,:),d(:),w(:),m(:),sg(:),
     +teq(:),ce(:),cn(:),ch(:),ces(:),cns(:),chs(:),pve(:),pvn(:),pvh(:)
     +,pves(:),pvns(:),pvhs(:),std(:),cal(:)

      pi=atan(1.)*4

      print*,'Start fit....'

      open(1,file='comfilt.para')
      open(2,file='comfilt.out')
      open(3,file='comfilt.inp',status='old')
      write(2,'("Number of stations: ", i5)')ns
      do s1=1,ns
        m1=6
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
            else if(sym(i)=='T') then
              m1=m1+3
              m2=m2+1
              m3=m3+1
            end if
          end do
        end if
        allocate(ce(m2),cn(m2),ch(m2),ces(m2),cns(m2),chs(m2),pve(m3),
     +  pvn(m3),pvh(m3),pves(m3),pvns(m3),pvhs(m3))
        do s2=1,3
          if(s2==1) infile='ts_'//sta(s1)//'_e_f.dat'
          if(s2==2) infile='ts_'//sta(s1)//'_n_f.dat'
          if(s2==3) infile='ts_'//sta(s1)//'_u_f.dat'
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
            tm=6
            g(i,1)=1
            g(i,2)=t(i)
            g(i,3)=sin(2*pi*t(i))
            g(i,4)=cos(2*pi*t(i))
            g(i,5)=sin(4*pi*t(i))
            g(i,6)=cos(4*pi*t(i))
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
                    g(i,tm)=H*(log10(t(i)-teq(j))-log10(0.00273))
                  end if
                else if(sym(j)=='L') then
                  tm=tm+1
                  g(i,tm)=H
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
                    g(i,tm)=H*(log10(t(i)-teq(j))-log10(0.00273))
                  end if
                end if
              end do
            end if
            w(i)=1./real(n)
          end do
          close(10)
          call lsqr(n,m1,g,d,m,w,std,sg,mis,rr)
          if(s2==1) out='ts_'//sta(s1)//'_e_t.dat'
          if(s2==2) out='ts_'//sta(s1)//'_n_t.dat'
          if(s2==3) out='ts_'//sta(s1)//'_u_t.dat'
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
            ve=m(2)
            se=sg(2)
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
                  pve(tm3)=m(tm)+m(2)
                  pves(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                else if(sym(i)=='T') then
                  tm=tm+3
                  tm2=tm2+1
                  tm3=tm3+1
                  ce(tm2)=m(tm-2)+m(tm)*log10(0.001)
                  ces(tm2)=sg(tm-2)
                  if(i==1) then
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
            sn=sg(2)
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
                  pvn(tm3)=m(tm)+m(2)
                  pvns(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                else if(sym(i)=='T') then
                  tm=tm+3
                  tm2=tm2+1
                  tm3=tm3+1
                  cn(tm2)=m(tm-2)+m(tm)*log10(0.001)
                  cns(tm2)=sg(tm-2)
                  if(i==1) then
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
            sh=sg(2)
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
                  pvh(tm3)=m(tm)+m(2)
                  pvhs(tm3)=sqrt(sg(tm)**2+sg(2)**2)
                else if(sym(i)=='T') then
                  tm=tm+3
                  tm2=tm2+1
                  tm3=tm3+1
                  ch(tm2)=m(tm-2)+m(tm)*log10(0.001)
                  chs(tm2)=sg(tm-2)
                  if(i==1) then
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
        write(2,'("# ",a4,"  Recorded Number of EQ: ",i3)')sta(s1),ty
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
              write(2,'("  ! earthquake occurred time (yr):",f11.4)')
     +        teq(i)
              write(2,'("    coseismic displacements (mm):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')ce(tm2),
     +        ces(tm2)
              write(2,'("      N component: ",f8.3," +-",f8.3)')cn(tm2),
     +        cns(tm2)
              write(2,'("      U component: ",f8.3," +-",f8.3)')ch(tm2),
     +        chs(tm2)
            else if(sym(i)=='L') then
              tm2=tm2+1
              tm3=tm3+1
              write(2,'("  ! earthquake occurred time (yr):",f11.4)')
     +        teq(i)
              write(2,'("    coseismic displacements (mm):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')ce(tm2),
     +        ces(tm2)
              write(2,'("      N component: ",f8.3," +-",f8.3)')cn(tm2),
     +        cns(tm2)
              write(2,'("      U component: ",f8.3," +-",f8.3)')ch(tm2),
     +        chs(tm2)
               write(2,'("    linear velocity (mm/yr):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')pve(tm3)
     +        ,pves(tm3)
              write(2,'("      N component: ",f8.3," +-",f8.3)')pvn(tm3)
     +        ,pvns(tm3)
              write(2,'("      U component: ",f8.3," +-",f8.3)')pvh(tm3)
     +        ,pvhs(tm3)
            else if(sym(i)=='T') then
              tm2=tm2+1
              tm3=tm3+1
              write(2,'("  ! earthquake occurred time (yr):",f11.4)')
     +        teq(i)
              write(2,'("    coseismic displacements (mm):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')ce(tm2),
     +        ces(tm2)
              write(2,'("      N component: ",f8.3," +-",f8.3)')cn(tm2),
     +        cns(tm2)
              write(2,'("      U component: ",f8.3," +-",f8.3)')ch(tm2),
     +        chs(tm2)
               write(2,'("    linear velocity (mm/yr):")')
              write(2,'("      E component: ",f8.3," +-",f8.3)')pve(tm3)
     +        ,pves(tm3)
              write(2,'("      N component: ",f8.3," +-",f8.3)')pvn(tm3)
     +        ,pvns(tm3)
              write(2,'("      U component: ",f8.3," +-",f8.3)')pvh(tm3)
     +        ,pvhs(tm3)
           end if
         end do
        end if
        if(tyeq==1) deallocate(teq,sym)
        deallocate(ce,cn,ch,ces,cns,chs,pve,pvn,pvh,pves,pvns,pvhs)
      end do
      close(3)
      close(2)
      close(1)

      print*,'End fit....'

      return
      end
