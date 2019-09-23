      program plot_ts

      implicit none
      integer,parameter::m1=10000
      character site(m1)*4,file(16)*15,fmt*200,maping*11,ps*10,para*1
      double precision,allocatable::teq(:)
      real*8 xmin,ymin,xmax,ymax,xminb,yminb,xmaxb,ymaxb
      real*8 dx,dy,dxb,dyb,pxminb,pxmaxb
      integer i,j,k,n,stat,neq
      data dxb,dyb /9.,8./
      logical alive

      open(11,file='plot_ts.inp',status='old')
      read(11,*) pxminb,pxmaxb
      read(11,*) neq
      if(neq>0) allocate(teq(neq))
      if(neq>0) read(11,*)(teq(i),i=1,neq)
      if(neq==0) read(11,*)
      read(11,*)
      n=0
      stat=0
      do while(.true.)
        n=n+1
        read(11,'(a4)',iostat=stat) site(n)
        if(stat/=0) exit
      end do
      close(11)
      n=n-1

      do i=1,n
        print*,'Processing in station: ',site(i)
#ifdef MINGW
        maping='map'//site(i)(1:4)//'.bat'
#else
        maping='map'//site(i)(1:4)//'.sh'
#endif
        ps='map'//site(i)(1:4)//'.ps'

        file(1)='ts_'//site(i)//'_n_o.dat'
        file(2)='ts_'//site(i)//'_e_o.dat'
        file(3)='ts_'//site(i)//'_u_o.dat'
        file(4)='ts_'//site(i)//'_n_c.dat'
        file(5)='ts_'//site(i)//'_e_c.dat'
        file(6)='ts_'//site(i)//'_u_c.dat'
        file(7)='ts_'//site(i)//'_n_f.dat'
        file(8)='ts_'//site(i)//'_e_f.dat'
        file(9)='ts_'//site(i)//'_u_f.dat'
        file(10)='ts_'//site(i)//'_n_t.dat'
        file(11)='ts_'//site(i)//'_e_t.dat'
        file(12)='ts_'//site(i)//'_u_t.dat'

        open(11,file=maping)
#ifdef MINGW
        write(11,'("@echo off")')
#else
        write(11,'("#!/bin/sh")')
c        write(11,'("")')
c        write(11,'("set -eux")')
#endif
        write(11,'("")')
        write(11,'("gmt gmtset FONT_ANNOT 8p FONT_LABEL 9p MAP_TICK_",
     + "LENGTH -0.030i MAP_ANNOT_OFFSET 0.035i MAP_LABEL_OFFSET 0.",
     + "0525i FONT_TITLE 12 MAP_TITLE_OFFSET 0.01i")')

        do j=1,3

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find the minmum and maxmum values for boundaries                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          xmin=99999.
          ymin=99999.
          xmax=-99999.
          ymax=-99999.

          inquire(file=file(j),exist=alive)             ! find the minmum and maxmum values from "o.dat"
          if(alive) call minmax(file(j),xmin,xmax,ymin,ymax,0)
          inquire(file=file(j+3),exist=alive)           ! find the minmum and maxmum values from "c.dat"
          if(alive) call minmax(file(j+3),xmin,xmax,ymin,ymax,0)
          inquire(file=file(j+6),exist=alive)           ! find the minmum and maxmum values from "f.dat"
          if(alive) call minmax(file(j+6),xmin,xmax,ymin,ymax,0)
          inquire(file=file(j+9),exist=alive)           ! find the minmum and maxmum values from "t.dat"
          if(alive) call minmax(file(j+9),xmin,xmax,ymin,ymax,0)


          ! find x-axis and y-axis boundaries
          call bound(xmax,xmin,0.1,xmaxb,xminb)
          call bound(ymax,ymin,0.1,ymaxb,yminb)

          if(pxminb>-9999) xminb=pxminb
          if(pxmaxb< 9999) xmaxb=pxmaxb

          call tick(xmax,xmin,xmaxb,xminb,dxb,dx,1.5,0.05,1.,10.)
          call tick(ymax,ymin,ymaxb,yminb,dyb,dy,5.5,0.10,5.,10.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     north component                                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if(j==1)then
            fmt='("gmt psbasemap -Bxa",f7.3,"+l""Year"" -Bya",f7.3,"+l""
     +North (mm)"" -BWesn -JX6.5i/1.8i -R",f7.3,"/",f7.3,"/",f7.3,"/",f7
     +.3," -K -P -V2 -Y6.2i -X.8i > ",a10)'
            fmt(24:24)=para(dx)
            fmt(47:47)=para(dy)
            fmt(94:94)=para(xminb)
            fmt(103:103)=para(xmaxb)
            fmt(112:112)=para(yminb)
            fmt(121:121)=para(ymaxb)
            write(11,fmt)dx,dy,xminb,xmaxb,yminb,ymaxb,ps
            do k=1,neq
              write(11,'("echo "f9.4,f10.3," >  sym")')teq(k),yminb
              write(11,'("echo "f9.4,f10.3," >> sym")')teq(k),ymaxb
              write(11,'("gmt psxy sym -JX -R -W0.2/50 -O -P -V2 -K >> "
     +        ,a10)')ps
            end do
            write(11,'("echo ",2f8.3," 10 0 0 ML Site: ",a4," | gmt pste
     +xt -JX -R -N -O -K -D-.7i/.2i -P -V2 >> ",a10)')xmaxb,ymaxb,site(i
     +),ps
            inquire(file=file(1),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +255/0/0 -G250/250/250 -Ey.05i/0.2,178/178/178 -P -K -O -V2 >> ",a1
     +0)')file(1),ps
            inquire(file=file(1),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +255/0/0 -G250/250/250 -P -K -O -V2 >> ",a10)')file(1),ps
            inquire(file=file(4),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +0/0/255 -G250/250/250 -P -K -O -V2 >> ",a10)')file(4),ps
            inquire(file=file(7),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +0/0/0 -G250/250/250 -P -K -O -V2 >> ",a10)')file(7),ps
            inquire(file=file(10),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -W0.2,0/180/0
     + -P -K -O -V2 >> ",a10)')file(10),ps

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     east component                                                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          else if(j==2)then
            fmt='("gmt psbasemap -Bxa",f7.3,"+l""Year"" -Bya",f7.3,"+l""
     +East(mm)"" -BWesn -JX6.5i/1.8i -R",f7.3,"/",f7.3,"/",f7.3,"/",f7.3
     +," -K -O -P -V2 -Y-1.95i -X0i >> ",a10)'
            fmt(24:24)=para(dx)
            fmt(47:47)=para(dy)
            fmt(92:92)=para(xminb)
            fmt(101:101)=para(xmaxb)
            fmt(110:110)=para(yminb)
            fmt(119:119)=para(ymaxb)
            write(11,fmt)dx,dy,xminb,xmaxb,yminb,ymaxb,ps
            do k=1,neq
              write(11,'("echo "f9.4,f10.3," >  sym")')teq(k),yminb
              write(11,'("echo "f9.4,f10.3," >> sym")')teq(k),ymaxb
              write(11,'("gmt psxy sym -JX -R -W0.2/50 -O -P -V2 -K >> "
     +        ,a10)')ps
            end do
            inquire(file=file(2),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +255/0/0 -G250/250/250 -Ey.05i/0.2,178/178/178  -P -K -O -V2 >> ",a
     +10)')file(2),ps
            inquire(file=file(2),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +255/0/0 -G250/250/250 -P -K -O -V2 >> ",a10)')file(2),ps
            inquire(file=file(5),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +0/0/255 -G250/250/250 -P -K -O -V2 >> ",a10)')file(5),ps
            inquire(file=file(8),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +0/0/0 -G250/250/250 -P -K -O -V2 >> ",a10)')file(8),ps
            inquire(file=file(11),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -W0.2,0/180/0 
     + -P -K -O -V2 >> ",a10)')file(11),ps

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     up component                                                                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          else if(j==3)then
            fmt='("gmt psbasemap -Bxa",f7.3,"+l""Year"" -Bya",f7.3,"+l""
     +Up (mm)"" -BWeSn -JX6.5i/1.8i -R",f7.3,"/",f7.3,"/",f7.3,"/",f7.3,
     +" -K -O -P -V2 -Y-1.95i -X0i >> ",a10)'
            fmt(24:24)=para(dx)
            fmt(47:47)=para(dy)
            fmt(91:91)=para(xminb)
            fmt(100:100)=para(xmaxb)
            fmt(109:109)=para(yminb)
            fmt(118:118)=para(ymaxb)
            write(11,fmt)dx,dy,xminb,xmaxb,yminb,ymaxb,ps
            do k=1,neq
              write(11,'("echo "f9.4,f10.3," >  sym")')teq(k),yminb
              write(11,'("echo "f9.4,f10.3," >> sym")')teq(k),ymaxb
              write(11,'("gmt psxy sym -JX -R -W0.2/50 -O -P -V2 -K >> "
     +        ,a10)')ps
            end do
            inquire(file=file(3),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +255/0/0 -G250/250/250 -Ey.05i/0.2,178/178/178 -P -K -O -V2 >> ",a1
     +0)')file(3),ps
            inquire(file=file(3),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +255/0/0 -G250/250/250 -P -K -O -V2 >> ",a10)')file(3),ps
            inquire(file=file(6),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +0/0/255 -G250/250/250 -P -K -O -V2 >> ",a10)')file(6),ps
            inquire(file=file(9),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -Sc.02i -W0.2,
     +0/0/0 -G250/250/250 -P -K -O -V2 >> ",a10)')file(9),ps
            inquire(file=file(12),exist=alive)
            if(alive) write(11,'("gmt psxy ",a15," -JX -R -W0.2,0/180/0 
     + -P -K -O -V2 >> ",a10)')file(12),ps
          end if
        end do
        write(11,'("gmt psxy -R -J -O -T >> ",a10)')ps
        write(11,'("gmt psconvert -Tg -A ",a10)')ps
#ifdef MINGW
        write(11,'("del sym gmt.history gmt.conf")')
#else
        write(11,'("rm -f sym gmt.history gmt.conf")')
#endif
        close(11)
#ifdef MINGW
        call system(maping)
#else
        call system('sh '//maping)
#endif
        print*,'======== Completed! ========'
      end do
      if(neq>0) deallocate(teq)

      stop
      end

c***********************************************************************
      subroutine minmax(inp,x1,x2,y1,y2,ty)

      implicit none
      integer stat,ty
      character inp*15,line*100
      real*8 x,y,x1,x2,y1,y2

      stat=0
      open(12,file=inp,status='old')
      do while (stat==0)
        if(ty==1) read(12,'(a)',iostat=stat)line
        if(ty==0) read(12,*,iostat=stat)x,y
        if(stat/=0)exit
        if(ty==1.and.line(1:1)=='X')cycle
        if(ty==1) read(line,*)x,y
        if(x<x1)x1=x
        if(x>x2)x2=x
        if(y<y1)y1=y
        if(y>y2)y2=y
      end do
      close(12)

      return
      end

c***********************************************************************
      subroutine bound(max,min,time,maxb,minb)

      implicit none
      real time
      real*8 max,min,maxb,minb

      maxb=max+(max-min)*time
      minb=min-(max-min)*time

      return
      end

c***********************************************************************
      subroutine tick(max,min,maxb,minb,db,d,v1,v2,v3,v4)

      implicit none
      integer i,num
      real v1,v2,v3,v4
      real*8 max,min,maxb,minb,db,d,mtl,mtld,tmp

      if((max-min)<=v1) mtl=v2
      if((max-min)>=v1) mtl=v3
      if((maxb-minb)<=(mtl*db)) mtld=mtl/v4
      if((maxb-minb)>(mtl*db)) mtld=mtl
      do i=1,1000
        tmp=(maxb-minb)/(mtld*i)
        if(tmp<=db)then
          num=i
          if(tmp<(db/2)) num=num-1
          d=mtld*num
          exit
        end if
      end do

      return
      end

c***********************************************************************
      function para(val)

      implicit none
      character para*1
      real*8 val

      if((val<=-1000).and.(val>-10000)) para='9'
      if((val<=-100).and.(val>-1000)) para='8'
      if((val<=-10).and.(val>-100)) para='7'
      if((val<0).and.(val>-10)) para='6'
      if((val<10).and.(val>=0)) para='5'
      if((val<100).and.(val>=10)) para='6'
      if((val<1000).and.(val>=100)) para='7'
      if((val<10000).and.(val>=1000)) para='8'

      return
      end
