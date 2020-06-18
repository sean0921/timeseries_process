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
        maping='map'//site(i)(1:4)//'.bat'
        ps='map'//site(i)(1:4)//'.ps'
        
        file(1)='ts_'//site(i)//'_n.gmt'
        file(2)='ts_'//site(i)//'_e.gmt'
        file(3)='ts_'//site(i)//'_u.gmt'
        file(4)='ts_'//site(i)//'_n_c.gmt'
        file(5)='ts_'//site(i)//'_e_c.gmt'
        file(6)='ts_'//site(i)//'_u_c.gmt'
        file(7)='ts_'//site(i)//'_n_f.gmt'
        file(8)='ts_'//site(i)//'_e_f.gmt'
        file(9)='ts_'//site(i)//'_u_f.gmt'
        file(10)='ts_'//site(i)//'_n_t.gmt'
        file(11)='ts_'//site(i)//'_e_t.gmt'
        file(12)='ts_'//site(i)//'_u_t.gmt'
        
        open(11,file=maping)
        write(11,'("echo off")')
        write(11,'("gmtset ANNOT_FONT_SIZE 8p LABEL_FONT_SIZE 9p TICK_",
     +  "LENGTH -0.030i ANNOT_OFFSET 0.035i LABEL_OFFSET 0.0525i HEADE",
     +  "R_FONT_SIZE 12 HEADER_OFFSET 0.01i")')

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
            fmt='("psbasemap -B",f7.3,":""Year"":/",f7.3,":""North (mm)"
     +":Wesn -JX6.5i/1.8i -R",f7.3,"/",f7.3,"/",f7.3,"/",f7.3," -K -P -V
     + -Y6.2i -X.8i > ",a10)'
            fmt(18:18)=para(dx)
            fmt(37:37)=para(dy)
            fmt(81:81)=para(xminb)
            fmt(90:90)=para(xmaxb)
            fmt(99:99)=para(yminb)
            fmt(108:108)=para(ymaxb)
            write(11,fmt)dx,dy,xminb,xmaxb,yminb,ymaxb,ps
            do k=1,neq
              write(11,'("echo "f9.4,f10.3," >  sym")')teq(k),yminb
              write(11,'("echo "f9.4,f10.3," >> sym")')teq(k),ymaxb
              write(11,'("psxy sym -JX -R -W3/0ta -O -P -V -K >> ",a10)'
     +        )ps
            end do
            write(11,'("echo ",2f8.3," 10 0 0 ML Site: ",a4," | pstext -
     +JX -R -N -O -K -D-.7i/.2i -P -V >> ",a10)')xmaxb,ymaxb,site(i),ps
            inquire(file=file(1),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/255/0/
     +0 -G255 -Ey.05i/2/125 -P -K -O -V >> ",a10)')file(1),ps
            inquire(file=file(4),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/0/0/25
     +5 -G255 -P -K -O -V >> ",a10)')file(4),ps
            inquire(file=file(7),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/0/0/0 
     +-G255 -P -K -O -V >> ",a10)')file(7),ps
            inquire(file=file(10),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -W3/0/180/0 -P -K 
     +-O -V >> ",a10)')file(10),ps

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     east component                                                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          else if(j==2)then
            fmt='("psbasemap -B",f7.3,":""Year"":/",f7.3,":""East (mm)""
     +:Wesn -JX6.5i/1.8i -R",f7.3,"/",f7.3,"/",f7.3,"/",f7.3," -K -O -P 
     +-V -Y-1.95i -X0i >> ",a10)'
            fmt(18:18)=para(dx)
            fmt(37:37)=para(dy)
            fmt(80:80)=para(xminb)
            fmt(89:89)=para(xmaxb)
            fmt(98:98)=para(yminb)
            fmt(107:107)=para(ymaxb)
            write(11,fmt)dx,dy,xminb,xmaxb,yminb,ymaxb,ps
            do k=1,neq
              write(11,'("echo "f9.4,f10.3," >  sym")')teq(k),yminb
              write(11,'("echo "f9.4,f10.3," >> sym")')teq(k),ymaxb
              write(11,'("psxy sym -JX -R -W3/0ta -O -P -V -K >> ",a10)'
     +        )ps
            end do
            inquire(file=file(2),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/255/0/
     +0 -G255 -Ey.05i/2/125 -P -K -O -V >> ",a10)')file(2),ps
            inquire(file=file(5),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/0/0/25
     +5 -G255 -P -K -O -V >> ",a10)')file(5),ps
            inquire(file=file(8),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/0/0/0 
     +-G255 -P -K -O -V >> ",a10)')file(8),ps
            inquire(file=file(11),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -W3/0/180/0 -P -K 
     +-O -V >> ",a10)')file(11),ps

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     up component                                                                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          else if(j==3)then
            fmt='("psbasemap -B",f7.3,":""Year"":/",f7.3,":""Up (mm)"":W
     +eSn -JX6.5i/1.8i -R",f7.3,"/",f7.3,"/",f7.3,"/",f7.3," -K -O -P -V
     + -Y-1.95i -X0i >> ",a10)'
            fmt(18:18)=para(dx)
            fmt(37:37)=para(dy)
            fmt(78:78)=para(xminb)
            fmt(87:87)=para(xmaxb)
            fmt(96:96)=para(yminb)
            fmt(105:105)=para(ymaxb)
            write(11,fmt)dx,dy,xminb,xmaxb,yminb,ymaxb,ps
            do k=1,neq
              write(11,'("echo "f9.4,f10.3," >  sym")')teq(k),yminb
              write(11,'("echo "f9.4,f10.3," >> sym")')teq(k),ymaxb
              write(11,'("psxy sym -JX -R -W3/0ta -O -P -V -K >> ",a10)'
     +        )ps
            end do
            inquire(file=file(3),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/255/0/
     +0 -G255 -Ey.05i/2/125 -P -K -O -V >> ",a10)')file(3),ps
            inquire(file=file(6),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/0/0/25
     +5 -G255 -P -K -O -V >> ",a10)')file(6),ps
            inquire(file=file(9),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -Sc.03i -W2/0/0/0 
     +-G255 -P -O -V >> ",a10)')file(9),ps
            inquire(file=file(12),exist=alive)
            if(alive) write(11,'("psxy ",a15," -JX -R -W3/0/180/0 -P -O 
     +-V >> ",a10)')file(12),ps
          end if
        end do
        write(11,'("ps2raster ",a10," -A -Tj")')ps
        write(11,'("del sym .gmt*")')
        close(11)
        call system(maping)
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
