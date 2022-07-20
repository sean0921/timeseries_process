c--     by Kuo-En Ching 2005.01.13
c--              Update 2008.05.28
c
c----bern2time.for---
c
      program bern2time

      use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
      implicit none
      character line*200
      character sta1*4,inpfile*12
      character,allocatable::outfile(:)*10,sta(:)*4
      real*8 lats,lons,height,herr,laterr,lonerr
      real*8 sigh,sign,sige,fve,fvn,fvh,t,t1,t2
      integer i,n,lat,lon,latm,lonm,stat,ch,sn,hh,mm,ss
      logical alive
      common/fix/fve,fvn,fvh

      print*,'                         PROGRAM BERN2TIME                             '
      print*,'                                            by Kuo-En Ching 2005.01.11 '
      print*,'                                                     Update 2008.05.28 '
      print*,'Seismology Lab, Department of Earth Sciences, NCKU   Update 2020.12.17+'
      print*,''
      print*,'source code: https://github.com/sean0921/timeseries_process'
      print*,''
      print*,'Please choose: 1. Calculate velocities relative to one station (Fixed Station)'
      print*,'               2. Calculate station coordinate variations'
      print*,'(1 or 2) '
      read(*,'(i1)') ch
      if(ch==1) then
        print*,'Enter vel. of Fixed Stat in three components:'
        print*,'ve (mm/yr): '
        read(*,*)fve
        print*,'vn (mm/yr): '
        read(*,*)fvn
        print*,'vu (mm/yr): '
        read(*,*)fvh
      else if(ch==2) then
        fvn=0.
        fve=0.
        fvh=0.
      else
        print*,'Input Error!!'
        read(stdin,*)
        stop
      end if

      inquire(file='fil.dat',exist=alive)
      if(alive) open(8,file='fil.dat',status='old')
      if(alive) close(8,status='delete')
      inquire(file='error.msg',exist=alive)
      if(alive) open(9,file='error.msg',status='old')
      if(alive) close(9,status='delete')
#ifdef PLATFORM_IS_WINDOWS
      call system('del ts_????_b.dat ts_????_n_b.dat ts_????_e_b.dat ts_????_u_b.dat')
      call system('for %f in (FN??????.OUT) do echo %f >> fil.dat')
#else
      call system('rm -f ts_????_b.dat ts_????_n_b.dat ts_????_e_b.dat ts_????_u_b.dat')
      call system('for f in $(ls FN??????.OUT);do echo $f;done >> fil.dat')
#endif

      open(15,file='fil.dat',status='old')
      n=0
      stat=0
      do while(stat==0)
        read(15,*,iostat=stat)
        if(stat/=0) exit
        n=n+1
      end do
      rewind(15)
      allocate(outfile(n))
      do i=1,n
        read(15,'(a12)',iostat=stat)inpfile
        if(stat/=0) exit
        print*,'Processing '//inpfile//' ...'
        outfile(i)=inpfile(3:8)//'.plh'
        open(11,file=inpfile)
        open(12,file=outfile(i))
        stat=0
        do while (stat==0)
          read(11,'(a)',iostat=stat) line
          if(stat/=0) exit
          if(line(1:15)==' FILE TYP FREQ.') then
            t=0.
            t1=0.
            t2=24.
            read(11,'(/)')
            do while (.true.)
              read(11,'(a)')line
              if(line(5:5)==' ') exit
              read(line,'(65x,i3,1x,i2,1x,i2)')hh,mm,ss
              t1=real(hh)+real(mm)/60.+real(ss)/3600.
              if(t1<t2.and.t1>=1.0) t2=t1
            end do
            if(t2<24.) t=t2
            if(t2>=23.99999) t=t1
            write(12,'(f7.4)')t
          end if
          if(line(24:24)=='X') then
            read(line,'(6x,a4)') sta1
          else if(line(24:29)=='HEIGHT') then
            read(line,'(59x,f9.4,18x,f6.4)')height,herr
            sigh=herr*1000*5
          else if(line(24:31)=='LATITUDE') then
            read(line,'(54x,i3,1x,i2,1x,f9.6,16x,f6.4)')lat,latm,lats,laterr
            sign=laterr*1000*10
          else if(line(24:32)=='LONGITUDE') then
            read(line,'(54x,i3,1x,i2,1x,f9.6,16x,f6.4)')lon,lonm,lons,lonerr
            sige=lonerr*1000*10
            write(12,'(a4,2(i4,i3,f10.6),f10.4,3f10.2)')sta1,lon,lonm,lons,lat,latm,lats,height,sige,sign,sigh
          end if
        end do

        close(11)
        close(12)
      end do
      close(15,status='delete')

      open (11,file='sta-file')
      sn=0
      stat=0
      do while (stat==0)
        read(11,*,iostat=stat)
        if(stat/=0) exit
        sn=sn+1
      end do
      rewind(11)
      allocate(sta(sn))
      do i=1,sn
        read(11,'(a4)')sta(i)
      end do
      close(11)

      call time_anal(n,outfile,sn,sta)
      deallocate(outfile,sta)
      stop
      end

c***********************************************************************
      subroutine time_anal(n,plhfil,sn,sta)

      implicit none
      integer i,j,n,sn,stat,lon,lat,lonm,latm
      character plhfil(n)*10,sta(sn)*4,yr*4,dy*3,t*8,sta1*4,out*8
      real*8 lons,lats,hgh,en,ee,eh,tt

      print*,'Start time_anal'

      do i=1,n
        open(11,file=plhfil(i),status='old')
        read(11,'(f7.4)',iostat=stat)tt
        if(stat/=0)exit
        dy(1:3)=plhfil(i)(3:5)
        yr(3:4)=plhfil(i)(1:2)
        if(yr(3:3)=='9') yr(1:2)='19'
        if(yr(3:3)/='9') yr(1:2)='20'
        t=yr//' '//dy
        stat=0
        do while (stat==0)
          read(11,'(a4,2(i4,i3,f10.6),f10.4,3f10.2)',iostat=stat)sta1,lon,lonm,lons,lat,latm,lats,hgh,ee,en,eh
          if(stat/=0) exit
          do j=1,sn
            if(sta1(1:4)==sta(j)(1:4)) then
              out(1:8)=sta(j)//'.out'
              open (12,file=out,position='append')
              write(12,'(a8,f8.4,2(i4,i3,f10.6),f10.4,3f10.2)')t,tt,lon,lonm,lons,lat,latm,lats,hgh,ee,en,eh
              close(12)
            end if
          end do
        end do
        close(11,status='delete')
      end do

      print*,'End time_anal'
      call coordvar(sn,sta)

      return
      end

c***********************************************************************
      subroutine coordvar(stn,sta)

      use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
      implicit none
      integer i,k,n,stn,stat,jud,lon,lat,lonm,latm
      integer,allocatable::yr(:),dy(:)
      character sta(stn)*4,inp*8
      character out1*15,out2*15,out3*15,out4*13
      real*8 lons,lats,avgn,avge,avgh,londt,latdt,pi,adn,ade,adh
      real*8 fve,fvn,fvh
      double precision,allocatable::rlon(:),rlat(:),hgh(:),sn(:),se(:),sh(:),t(:),time(:),dn(:),de(:),dh(:)
      common/fix/fve,fvn,fvh
      data latdt /111.325/
      pi=atan(1.)*4

      print*,'Start coordvar'

      do i=1,stn
        inp=sta(i)//'.out'
        out1='ts_'//sta(i)//'_n_b.dat'
        out2='ts_'//sta(i)//'_e_b.dat'
        out3='ts_'//sta(i)//'_u_b.dat'
        out4='ts_'//sta(i)//'_b.dat'

        open(12,file=inp)
        open(13,file=out1)
        open(14,file=out2)
        open(15,file=out3)
        open(16,file=out4)

        write(16,'(a)')'year doy      hr       dn      sn      de      se      dh      sh (unit: mm)'

        n=0
        avgn=0.
        avge=0.
        avgh=0.
        stat=0
        do while (stat==0)
          read(12,*,iostat=stat)
          if(stat/=0) exit
          n=n+1
        end do
        if(n==1) then
          jud=1
          close(12,status='delete')
          close(13,status='delete')
          close(14,status='delete')
          close(15,status='delete')
          close(16,status='delete')
          open(11,file='error.msg',position='append')
          write(11,'("Station ",a4," is not found!")')sta(i)
          write(*,'(" Warning!! Station ",a4," is not found!")')sta(i)
          close(11)
          cycle
        end if
        rewind(12)
        allocate(yr(n),dy(n),t(n),time(n),hgh(n),se(n),sn(n),sh(n),rlon(n),rlat(n),dn(n),de(n),dh(n))
        do k=1,n
          read(12,'(2i4,f8.4,2(i4,i3,f10.6),f10.4,3f10.2)',iostat=stat)yr(k),dy(k),t(k),lon,lonm,lons,lat,latm,lats,hgh(k),se(k),sn(k),sh(k)
          if(stat/=0) exit
          rlon(k)=lon+lonm/60.+lons/3600.
          rlat(k)=lat+latm/60.+lats/3600.
          if(mod(yr(k),4)==0) then
            time(k)=real(yr(k))+(real(dy(k))-1.)/366.+t(k)/8784.
          else
            time(k)=real(yr(k))+(real(dy(k))-1.)/365.+t(k)/8760.
          end if
          avgn=avgn+rlat(k)
          avge=avge+rlon(k)
          avgh=avgh+hgh(k)
        end do

        avgn=avgn/real(n)
        avge=avge/real(n)
        avgh=avgh/real(n)

        call bubble_sort_10(time,yr,dy,t,rlon,rlat,hgh,se,sn,sh,n)

        adn=0
        ade=0
        adh=0
        do k=1,n
          londt=cos(avgn*pi/180.)*latdt
          dn(k)=(rlat(k)-avgn)*latdt*1000000.-(fvn*(time(k)-time(1)))
          de(k)=(rlon(k)-avge)*londt*1000000.-(fve*(time(k)-time(1)))
          dh(k)=(hgh(k)-avgh)*1000.-(fvh*(time(k)-time(1)))
          adn=adn+dn(k)
          ade=ade+de(k)
          adh=adh+dh(k)
        end do
        adn=adn/real(n)
        ade=ade/real(n)
        adh=adh/real(n)
        do k=1,n
          write(13,'(f10.5,1x,2f8.2)') time(k),dn(k)-adn,sn(k)
          write(14,'(f10.5,1x,2f8.2)') time(k),de(k)-ade,se(k)
          write(15,'(f10.5,1x,2f8.2)') time(k),dh(k)-adh,sh(k)
          write(16,'(2i4,f8.4,1x,6f8.2)') yr(k),dy(k),t(k),dn(k)-adn,sn(k),de(k)-ade,se(k),dh(k)-adh,sh(k)
        end do
        close(12,status='delete')
        close(13)
        close(14)
        close(15)
        close(16)

        deallocate(yr,dy,t,time,hgh,se,sn,sh,rlon,rlat,dn,de,dh)
      end do

      print*,'End coordvar'
      if(jud==1) read(stdin,*)

      return
      end
