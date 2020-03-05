      program remove

      implicit none
      character infile*15,outfile*15,line*500,sta*4,stap*4
      character rc*1,c*1
      character,allocatable::sym(:)*1
      real*8 cal,t,d,s,pi
      double precision,allocatable::m(:),teq(:),g(:)
      integer i,j,s1,s2,m1,H,ns,stat,neq,tyeq,tm,ty,nT

      pi=atan(1.)*4

      open(1,file='remove.inp',status='old')
      open(2,file='remove.para',status='old')
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
            if(sym(i)=='V') m1=m1+1
            if(sym(i)=='T') m1=m1+3
          end do
        end if
        do s2=1,3
          if(s2==1) infile='ts_'//sta//'_e_o.dat'
          if(s2==2) infile='ts_'//sta//'_n_o.dat'
          if(s2==3) infile='ts_'//sta//'_u_o.dat'
          if(s2==1) outfile='ts_'//sta//'_e_c.dat'
          if(s2==2) outfile='ts_'//sta//'_n_c.dat'
          if(s2==3) outfile='ts_'//sta//'_u_c.dat'
          if(s2==1) rc='e'
          if(s2==2) rc='n'
          if(s2==3) rc='u'
          allocate(m(m1),g(m1))
          do i=1,ns*6
            read(2,'(a)')line
            read(line,*)stap,c
            if(stap(1:4)==sta(1:4).and.c(1:1)==rc(1:1)) then
              read(line,*)stap,c,(m(j),j=1,m1)
              exit
            end if
          end do
          print*,'  Processing file: ',infile
          open(10,file=infile,status='old')
          open(11,file=outfile)
          stat=0
          do while (stat==0)
            read(10,*,iostat=stat)t,d,s
            if(stat/=0) exit
            tm=6
            g(1)=1
            g(2)=t
            g(3)=sin(2*pi*t)
            g(4)=cos(2*pi*t)
            g(5)=sin(4*pi*t)
            g(6)=cos(4*pi*t)
            if(neq>0) then
              do j=1,neq
                if(t>teq(j)) H=1
                if(t<=teq(j))  H=0
                if(sym(j)=='C') then
                  tm=tm+1
                  g(tm)=H
                else if(sym(j)=='P') then
                  tm=tm+1
                  if(H==0) then
                    g(tm)=0.
                  else if(H==1.and.abs(t-teq(j))>0.001) then
                    g(tm)=H*(log10(t-teq(j))-log10(0.00273))
                  end if
                else if(sym(j)=='L') then
                  tm=tm+1
                  g(tm)=H
                  tm=tm+1
                  g(tm)=H*(t-teq(j))
                else if(sym(j)=='V') then
                  tm=tm+1
                  g(tm)=H*(t-teq(j))
                else if(sym(j)=='T') then
                  tm=tm+1
                  g(tm)=H
                  tm=tm+1
                  g(tm)=H*(t-teq(j))
                  tm=tm+1
                  if(H==0) then
                    g(tm)=0.
                  else if(H==1.and.abs(t-teq(j))>0.001) then
                    g(tm)=H*(log10(t-teq(j))-log10(0.00273))
                  end if
                end if
              end do
            end if
            if(neq==0) then
              cal=0
              do j=3,6
                cal=cal+g(j)*m(j)
              end do
              write(11,'(f10.5,2f10.2)')t,d-cal,s
            elseif(neq>=1) then
              nT=0
              ty=0
              cal=0
              do j=3,m1
                if(ty==1) then
                  if(nT==0) then
                    ty=0
                    cycle
                  else if(nT==1) then
                    nT=0
                    ty=1
                    cycle
                  end if
                end if
                if(j<7)then
                  cal=cal+g(j)*m(j)
                else if(j>=7.and.abs(teq(j-6)-2010.16619)<0.003.and.sym(j-6)=='C')then
                  cycle
                else if(j>=7.and.abs(teq(j-6)-2010.16619)<0.003.and.sym(j-6)=='L')then
                  ty=1
                  cycle
                else if(j>=7.and.abs(teq(j-6)-2010.16619)<0.003.and.sym(j-6)=='T')then
                  nT=1
                  ty=1
                  cycle
                else
                  cal=cal+g(j)*m(j)
                end if
              end do
              write(11,'(f10.5,2f10.2)')t,d-cal,s
            end if
          end do
          close(10)
          close(11)
          deallocate(m,g)
        end do
        if(tyeq==1) deallocate(teq,sym)
      end do
      close(1)
      close(2)

      stop
      end
