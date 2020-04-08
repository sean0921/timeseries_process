c--     by Kuo-En Ching 2005.01.13
c--              Update 2005.09.19
c
c
      program getpl

      implicit none
      character line*200,sta1*4,inpfile*13
      character,allocatable::sta(:)*4
      real*8 lats,lons,hgh,lat,lon
      double precision,allocatable::phy(:),lda(:)
      integer i,j,n,latd,lond,latm,lonm,stat,ty
      logical alive

      print*,'                            PROGRAM GETPL                              '
      print*,'                                            by Kuo-En Ching 2005.01.13 '
      print*,'                                                     Update 2005.09.19 '
      print*,'Seismology Lab, Department of Earth Sciences, NCKU   Update 2020.04.08+'
      print*,''
      print*,'source code: https://github.com/sean0921/timeseries_process'
      print*,''

      inquire(file='fil.dat',exist=alive)
#ifdef PLATFORM_IS_WINDOWS
      if(alive) call system('del fil.dat')
      call system('for %f in (FN??????.OUT) do echo %f >> fil.dat')
#else
      if(alive) call system('rm -f fil.dat')
      call system('for f in $(ls FN??????.OUT);do echo $f;done>fil.dat')
#endif

      open(12,file='tmp')
      open(15,file='fil.dat',status='old')
      stat=0
      do while(stat==0)
        read(15,*,iostat=stat)
        if(stat/=0) exit
        n=n+1
      end do
      rewind(15)
      do i=1,n
        read(15,'(a13)',iostat=stat) inpfile
        if(stat/=0) exit
        print*,'Processing '//inpfile//' ...'
        open(11,file=inpfile)
        stat=0
        do while (stat==0)
          read(11,'(a)',iostat=stat) line
          if(line(24:24).eq.'X') then
            read(line,'(6x,a4)') sta1
          else if(line(24:29).eq.'HEIGHT') then
            read(line,'(59x,f9.4)')hgh
          else if(line(24:31).eq.'LATITUDE') then
            read(line,'(54x,i3,1x,i2,1x,f9.6)')latd,latm,lats
            lat=real(latd)+real(latm)/60.+lats/3600.
          else if(line(24:32).eq.'LONGITUDE') then
            read(line,'(54x,i3,1x,i2,1x,f9.6)')lond,lonm,lons
            lon=real(lond)+real(lonm)/60.+lons/3600.
            write(12,'(a4,1x,2f8.4,f10.4)')sta1,lon,lat,hgh
          end if
        end do
        close(11)
      end do
      close(15,status='delete')
      close(12)

      open(12,file='tmp',status='old')
      n=0
      stat=0
      do while (stat==0)
        read(12,*,iostat=stat)
        if(stat/=0) exit
        n=n+1
      end do
      rewind(12)
      allocate(sta(n),phy(n),lda(n))
      do i=1,n
        read(12,*)sta(i),lda(i),phy(i)
      end do
      close(12,status='delete')

      open(12,file='getpl_o.gout')
      write(12,'(a4,1x,2f8.4)')sta(1),lda(1),phy(1)
      do i=2,n
        ty=0
        do j=1,i-1
          if(sta(i)==sta(j)) then
            ty=1
            exit
          end if
        end do
        if(ty==1) cycle
        write(12,'(a4,1x,2f8.4)')sta(i),lda(i),phy(i)
      end do
      deallocate(sta,phy,lda)
      rewind(12)
      n=0
      stat=0
      do while (stat==0)
      read(12,*,iostat=stat)
        if(stat/=0) exit
        n=n+1
      end do
      rewind(12)
      allocate(sta(n),phy(n),lda(n))
      do i=1,n
        read(12,*)sta(i),lda(i),phy(i)
      end do
      close(12,status='delete')

      call bubble_sort_3(sta,phy,lda,n)

      open(12,file='getpl.gout')
      open(13,file='sta-file')
      do i=1,n
        write(12,'(a4,1x,2f8.4)')sta(i),lda(i),phy(i)
        write(13,'(a4)')sta(i)
      end do
      close(12)
      close(13)

      deallocate(sta,phy,lda)

      stop
      end
