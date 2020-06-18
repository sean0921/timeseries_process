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
        m1=2
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
          if(s2==1) infile='ts_'//sta//'_e.gmt'
          if(s2==2) infile='ts_'//sta//'_n.gmt'
          if(s2==3) infile='ts_'//sta//'_u.gmt'
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
            std(i)=abs(std(i))
            if(std(i)<0.0001) std(i)=0.001d0
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
          if(s2==1) outfile='ts_'//sta//'_e_r.gmt'
          if(s2==1) outfil1='ts_'//sta//'_e_c.gmt'
          if(s2==2) outfile='ts_'//sta//'_n_r.gmt'
          if(s2==2) outfil1='ts_'//sta//'_n_c.gmt'
          if(s2==3) outfile='ts_'//sta//'_u_r.gmt'
          if(s2==3) outfil1='ts_'//sta//'_u_c.gmt'
          open(10,file=outfile)
          open(11,file=outfil1)
          cal=d-matmul(g,m)
          do i=1,n
            if(dabs(std(i))>3*mis) cycle
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

      do c=1,3
        open(10,file='comfilt.inp',status='old')
        open(11,file='fil.dat')
        do i=1,ns
          read(10,'(a4)')sta(i)
          if(c==1) infile='ts_'//sta(i)//'_e_r.gmt'
          if(c==2) infile='ts_'//sta(i)//'_n_r.gmt'
          if(c==3) infile='ts_'//sta(i)//'_u_r.gmt'
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
          if(c==1) infile='ts_'//sta(i)//'_e_r.gmt'
          if(c==2) infile='ts_'//sta(i)//'_n_r.gmt'
          if(c==3) infile='ts_'//sta(i)//'_u_r.gmt'
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
          if(c==1) infile='ts_'//sta(i)//'_e_c.gmt'
          if(c==1) outfile1='ts_'//sta(i)//'_e_f.gmt'
          if(c==2) infile='ts_'//sta(i)//'_n_c.gmt'
          if(c==2) outfile1='ts_'//sta(i)//'_n_f.gmt'
          if(c==3) infile='ts_'//sta(i)//'_u_c.gmt'
          if(c==3) outfile1='ts_'//sta(i)//'_u_f.gmt'
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
          if(s2==1) infile='ts_'//sta(s1)//'_e_f.gmt'
          if(s2==2) infile='ts_'//sta(s1)//'_n_f.gmt'
          if(s2==3) infile='ts_'//sta(s1)//'_u_f.gmt'
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
          if(s2==1) out='ts_'//sta(s1)//'_e_t.gmt'
          if(s2==2) out='ts_'//sta(s1)//'_n_t.gmt'
          if(s2==3) out='ts_'//sta(s1)//'_u_t.gmt'
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
      if (a==51) b='61'
      if (a==52) b='62'
      if (a==53) b='63'
      if (a==54) b='64'
      if (a==55) b='65'
      if (a==56) b='66'
      if (a==57) b='67'
      if (a==58) b='68'
      if (a==59) b='69'
      if (a==60) b='70'
      if (a==51) b='71'
      if (a==52) b='72'
      if (a==53) b='73'
      if (a==54) b='74'
      if (a==55) b='75'
      if (a==56) b='76'
      if (a==57) b='77'
      if (a==58) b='78'
      if (a==59) b='79'
      if (a==60) b='80'
      if (a==51) b='81'
      if (a==52) b='82'
      if (a==53) b='83'
      if (a==54) b='84'
      if (a==55) b='85'
      if (a==56) b='86'
      if (a==57) b='87'
      if (a==58) b='88'
      if (a==59) b='89'
      if (a==60) b='90'
      if (a==51) b='91'
      if (a==52) b='92'
      if (a==53) b='93'
      if (a==54) b='94'
      if (a==55) b='95'
      if (a==56) b='96'
      if (a==57) b='97'
      if (a==58) b='98'
      if (a==59) b='99'
      return
      end
