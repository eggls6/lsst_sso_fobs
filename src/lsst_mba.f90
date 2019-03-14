program lsst_mba
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! LSST Solar System Object observation simulator
!
! needs: 1)LSST observtion campaign database such as kraken_2026.db (exported to ASCII for the moment)
!        2)SSO orbit file (OpenOrb cometary elements)
!        3)(optional) JPL planetary ephemeris and spice kernels (e.g. DE431 etc.) for Earth position and orientation
!
! written by S. Eggl 20181101
!
! last modified 20190311
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

!/////////////////////
!MODULES
  use orbdyn2_m, only:kema2rv,semia2n,epoch2ma,com2ke,pi,pix2,deg2rad,rad2deg,dist3d,lighttraveltime
  use trafo3_m, only:icrf2he,lat2cart,cart2lat,he2icrf,icrf2radec,rvhe2icrf,rvicrf2dradec,norm,angv1v2,angv1v2n
  use photometry_m
  use camera_m
  use fileio_m
  use earth_m
  use spice_m, only:getrvearth,load_spk_earth,unload_spk_earth ,rotearth
!/////////////////////
!VARIABLES
  implicit none
      ! aberration correction, use Ephemeris for the Earth, use LSST camera
      logical::abcorr,simplearth,cam
      
      ! object recorded/visible
      logical,dimension(:),allocatable::recorded,visible

      ! end of file, allocation status, file unit number
      integer::eof,astat, un
      
      ! dummy integers
      integer::ia,idum(6)
      
      ! number of header lines, length of orbit name
      integer::forbhl,fobshl,lenname
      
      ! number of SSOs in sample
      integer*8::nsmpl
      
      ! linked list of visible asteroids to original orbitid
      integer*8,dimension(:),allocatable::visid
      
      ! number of visible asteroids
      integer*8::nvis
      
      ! orbit identifaction numbers
      integer*8,dimension(:),allocatable::orbid,orbid0
      
      ! loop indices
      integer*8::i,j,k
       
      ! number of observations, number of orbits 
      integer*8::nobs,norbs
      
      ! observation time [JD], pericenter distance [au]
      real*8::tobsjd,qorb
      
      ! opening angle of FOV cone, same x 2
      real*8::fovoa,fovoax2 
      ! distances: Earth asteroid, Sun asteroid, Sun earth, phase angle, solar elongation
      real*8,dimension(:),allocatable::dea,dsa,dse,phase,solarelon
      
      ! visual magnitudes
      real*8,dimension(:),allocatable::Vmag
      real*8::Vism
      
      ! relativ position and velocity Earth asteroid, ecliptic distance, latitude and longitude
      real*8::rvrea(6),obsdeclip(3)
      
      ! center of FOV elevation and azimuth
      real*8,dimension(:),allocatable::elf,azf
      
      ! center of FOV right ascension and declination
      real*8,dimension(:),allocatable::raf,decf

      real*8,dimension(:),allocatable::lst !local sideral time
      real*8,dimension(:),allocatable::rot !FOV rotation angle in Azimuth Elevation frame (counter clockwise)
      !observation time, exposure time, airmass, five sigmal depth, effective seeing FWHM
      real*8,dimension(:),allocatable::tobs, &
                texp_s,airmass,fivesigdepth,seeingFwhmEff
    
      ! right ascension and declination, radec velocities, ecliptic coordinates
      real*8,dimension(:,:),allocatable::rradec,drradec,declip
      real*8::az,el !azimuth, elevation
      real*8,dimension(:),allocatable::camx,camy !camera coordinates (az el system)
      real*8::rdum(10)     !dummy variables
      real*8,dimension(1)::fsdmax,fsdmin !max Vmag for five sigma depth throughout the program

      ! time of ephemerides, new time of ephemerides
      real*8::teph,tenew
      
      ! FOV pointing vector and ICRF direction
      real*8::dx(3),dicrf(3)
   
      ! orbit variables, distances and velocities of asteroid, relative vector Earth asteroid
      real*8,dimension(:,:),allocatable::orbs,oast,rva,rea
      ! initial mean anomaly
      real*8,dimension(:),allocatable::ma0
      
      ! state vector of the Earth 
      real*8::rve(6)
      ! angle between FOV pointing direction and Earth asteroid vector, five sigma depth + 1 [mag]
      real*8::obsang,fsdp1
      
      !LSST latitude, xyz in various frames
      real*8::lsstlat,lsstxyz0(1:3),lsstxyz(1:3),lsstxyz1(1:3)
      
      ! orbital elements for Earth, epoch, mean anomaly, mean motion
      real*8::oearth(1:6),t0e,me0,nearth
      
      !filter
      character,dimension(:),allocatable::filter
      
      !file names: orbit, observations, output
      character(len=80)::forb,fobs,fout
      !heaer,orbit format, command line arguments, dummy string 
      character(len=80)::hdr,forbfmt,arg(4),dumc
      
      !name of orbits
      character(len=:),dimension(:),allocatable::orbname0,orbname
      
      !format strings
      character(len=200)::fmt
      character(len=3)::afmt
      

      
!/////////////////////
!PARAMETERS
      
      real*8,parameter::photG=0.15d0 !photometric G factor for Vmag calc (JPL)
      real*8,dimension(3),parameter::LSSTlatlonalt=(/-30.244639d0,-70.749417d0,2.663d0/) !latitude [deg], longitude(E) [deg] and altidue [km] of LSST observatory 
      real*8,parameter::earthapo=1.02d0   !Earth's maximum apocenter distance in au
      
      abcorr=.true.         !Aberration correction (requires an additional ephemeris calculation step)
      simplearth=.false.    !.true.:use 2body propagation rather than ephemeris to get the position of the Earth. All observations are geocentric
                            !.false.: use JPL Ephemeris + SPICE for Earth position, rotation and orientation of observer

      tenew=1.d0            !time between calculation of new ephemeris for the whole set of SSOs (once every day)
                            !visible SSOs are updated more often
      cam=.true.           !.true.:use LSST camera / .false. FOV circle with radius 1.75 deg
   

!///////////
! COMMAND LINE ARGUMENTS & I/O
! GET PROGRAM INPUT FILES AS COMMAND LINE ARGUMENTS 
! START PROGRAM WITH ./PROG.EXE INPUT1 INPUT2 OUTPUT
      
      ia=1
      DO
        CALL get_command_argument(ia, arg(ia))
        IF (LEN_TRIM(arg(ia)) == 0) EXIT 
        !WRITE (*,*)ia, TRIM(arg(ia))
        ia = ia+1
      END DO
      
      if(ia.lt.3) then
       write(*,*)'Input error, 3 input file names separated by blanks are required: '
       write(*,*)'1) LSST observation campaing file'
       write(*,*)'2) Orbit file containing orbit information for SSOs'
       write(*,*)'3) output file name'
       STOP
      end if

!///////////
! LSST OBSERVATION DATABASE (ASCII TEXT INPUT BECAUSE I AM LAZY)
! COLUMNS FROM SELECT * IN SQL DATABASE 
      !fobs='bs2018b2.txt'
      fobs=trim(arg(1))
      !number of header lines
      fobshl=1
      write(*,*)'LSST observation campaing file: ',fobs


!///////////
! ORBIT FILE AND ORBIT NAME FORMAT (ASCII TEXT)
! FORMAT FROM OPENORB
      forb=trim(arg(2))
      !forb='mba_H14_5to25_a0_3.orb'
      !forb='mba_H25to26_a0_3.orb'
      forbfmt='COM' 
     
      !number of header lines
      forbhl=1
      
      !FILE UNIT    
      un=21
      ! get length of OPENORB format name filed from first string character count
      open(un,file=forb,status='old')
      do i=1,forbhl+1
        read(un,*)dumc
      end do
      lenname=LEN_TRIM(ADJUSTL(dumc))
      close(un)

      write(*,*)'Orbit file: ',forb

!///////////
! OUTPUT FILE AND FORMAT      

      fout=trim(arg(3))
      write(*,*)'Output file: ', fout
      !output file format
      write(afmt,'(I3)')lenname
      fmt='(A'//TRIM(afmt)//',1X,F15.7,2(1X,F15.10),2(1X,F15.7),5(1X,F15.10),4(1X,F7.3),1X,I12,12(1X,F15.10),1X,A2)'

!///////////
! PROGRAM START
      fovoa=fovang !find the opening angle of observation cone (from camera_m module) [rad]
      fovoax2=fovoa*2.d0
      lsstlat=LSSTlatlonalt(1)*deg2rad

        
      if(simplearth) then
      ! Orbit of Earth at starting time for LSST (JD 2459853.016793982)
        oearth=(/1.000493878204624E+00,1.698422446885815E-02, &
        1.012458402718546E-03,3.116100091825618E+02,1.498214891834334E+02,2.675305461889174E+02/)
        oearth(4:6)=oearth(4:6)*deg2rad
        t0e=59853.016793982
        me0=oearth(6)
        nearth=9.848794429170397E-01*deg2rad !mean motion of the Earth
        else

        call latlonalt2r(LSSTlatlonalt(1),LSSTlatlonalt(2),LSSTlatlonalt(3),lsstxyz0(1:3))
        call load_spk_earth
       end if

      write(*,*)"obs input file opened: ",fobs
      write(*,*)"counting observations..."
      
!/////////                        
!READ LSST OBS SCHEDULE FILE  

      call get_file_stats(un,fobs,1,(/25/),fobshl,nobs,fsdmax,fsdmin)
      write(*,*)'number of LSST observations to be processed: ',nobs
      write(*,*)'maximum 5-sig depth during LSST campaign [mag]: ',fsdmax
      write(*,*)'minimum 5-sig depth during LSST campaign [mag]: ',fsdmin
      
      allocate(tobs(nobs),lst(nobs),filter(nobs),raf(nobs),decf(nobs),rot(nobs), &
                elf(nobs),azf(nobs),texp_s(nobs),airmass(nobs),seeingFwhmEff(nobs),fivesigdepth(nobs),stat=astat)
      if(astat.ne.0) then
        write(*,*)'ERROR: Could not allocate enough memory LSST observations. Reduce number of observations.'
        STOP
      end if
      
      open(un,file=fobs,status='old')
      if(fobshl.gt.0) call read_hdr(un,fobshl,hdr)
      eof=0
      teph=0.d0
        
      write(*,*)'reading LSST observations...' 
      do i=1,nobs 
       !read(19,*,iostat=eof)j,tobs,raf,decf,alt_deg,texp_s,airmass,fivesigdepth
       read(un,*,iostat=eof)j,idum(1:2),rdum(1),tobs(i),lst(i),idum(3:5),filter(i),raf(i),decf(i),rot(i), &
                            elf(i),azf(i),idum(6),rdum(2),texp_s(i),airmass(i),rdum(3:6),seeingFwhmEff(i),fivesigdepth(i)       
       if(eof.ne.0) exit
      end do
      
      !convert angular measurements from [deg] to [rad]
      lst(:)=lst(:)*deg2rad
      raf(:)=raf(:)*deg2rad
      decf(:)=decf(:)*deg2rad
      elf(:)=elf(:)*deg2rad
      azf(:)=azf(:)*deg2rad
      rot(:)=rot(:)*deg2rad
      
       !to make sure that the az el values are compatible with the ones calculated for SSO's calculate the az el values from RADEC 
      !with the same transformation - they conincide up to 10^-4 [deg] with those in the file (perhaps due to other latitude values for LSST)
      do i=1,nobs
       call radec2azel(raf(i),decf(i),lst(i),lsstlat,azf(i),elf(i))
      end do
      
      close(un)
      write(*,*)'number of LSST observations read:',j

!/////////////////////
! COUNT NUMBER OF SSOs IN INPUT FILE
      write(*,*)'counting number of orbits in file ',forb
      call file_cnt_ntrs(un,forb,forbhl,norbs)
      write(*,*)'number of orbits:',norbs
      write(*,*)'reading orbits...'
      
      allocate(character(lenname) :: orbname0(norbs))
      allocate(orbs(norbs,10),visid(norbs),orbid0(norbs),stat=astat)
      if(astat.ne.0) then
        write(*,*)'ERROR: Could not allocate enough memory for asteroid orbits. Reduce number of asteroids in sample.'
        STOP
      end if
      
!/////////
! READ ORBITS
      open(un,file=forb,status='old')
      !db
      !write(*,*)'a,e,i,w,node,m,n,epoch,H'
       !write(66,*)'q,a,e,i,w,node,m,n,epoch,H,qorb-earthapo,Vism,fsdmax(1)'
      !edb
      if (trim(forbfmt)=='COM') then        
       call readopenorb(un,forbhl,lenname,norbs,orbs,orbid0,orbname0)
       !oast(1:10)={a,e,i,w,node,mean_anomaly,n,epoch,Hmag,Vmag}
      end if
       nsmpl=0
       do i=1,norbs
       !FILTER those objects that can never be seen immediately
        qorb=orbs(i,1)*(1.d0-orbs(i,2))
        Vism=VismagHG( orbs(i,9),0.d0,qorb,qorb-earthapo,photG)
        if(Vism.le.fsdmax(1)) then
            nsmpl=nsmpl+1
            visid(nsmpl)=i  
        end if      
       end do
      close(un)
      
      write(*,*)"orbit input file read: ",forb

      allocate(character(lenname) :: orbname(nsmpl))
      allocate(oast(nsmpl,10),rva(nsmpl,6),ma0(nsmpl),rea(nsmpl,3),&
              orbid(nsmpl),visible(nsmpl),stat=astat)
              
      do i=1,nsmpl
        oast(i,1:2)=orbs(visid(i),1:2)
        oast(i,3:6)=orbs(visid(i),3:6)*deg2rad
        oast(i,7:10)=orbs(visid(i),7:10)
        ma0(i)=oast(i,6)
        orbid(i)=orbid0(visid(i))
        orbname(i)=orbname0(visid(i))
      end do
      
      deallocate(orbs,orbid0,orbname0,visid)
      !number of orbits after filtering
      write(*,*)'number of orbits after filtering',nsmpl
 
    
!/////////      
!OPEN OUTPUT FILE
      open(un,file=fout,status='replace')
      write(un,*)'objId ' ,'time ' ,'ra ' ,'dec ','dradt ','ddecdt ' ,'phase ','solarelon ' ,'helio_dist ', &
                 'geo_dist ','tar_geo_dist ','visitExposureTime ','fiveSigmaDepth ','Hmag ','Vmag ', &
                 'orbId ','ecLong ','ecLat ','x ','y ','z ', 'FOVra ','FOVdec ', &
                 'FOVaz ','FOVel ','camx ','camy ','seeing_FWHM_EFF ', 'filter '

     allocate(visid(1:nsmpl),rradec(1:nsmpl,1:3),drradec(1:nsmpl,1:3),&
             phase(1:nsmpl),solarelon(1:nsmpl),dsa(1:nsmpl),dea(1:nsmpl),dse(1:nsmpl),&
              Vmag(1:nsmpl),declip(1:nsmpl,1:3),camx(1:nsmpl),camy(1:nsmpl),recorded(1:nsmpl),stat=astat)
     if(astat.ne.0) then
        write(*,*)'ERROR: Could not allocate enough memory for asteroid orbits. Reduce number of asteroids in sample.'
        STOP
     end if
              
           
    write(*,*)'processing...' 
    recorded(:)=.false.

!////////////////
! MAIN PROGRAM LOOP
      do k=1,nobs
       !5 SIG DEPTH FILTER
       fsdp1=fivesigdepth(k)+1.d0
       !PROGRESS OUTPUT TO SCREEN 
       if(mod(k,500)==0) then
          write(*,*)'survey time',tobs(k),'# asteroids',nsmpl,'# visible',nvis
       end if

       !determine observation direction
       dicrf(:)=(/1.d0,decf(k),raf(k)/)
       !dx is the observation direction unit vector
       call lat2cart(dicrf,dx)
       call icrf2he(dx)
       call cart2lat(dx,obsdeclip(1:3))
      ! write(*,*)'ecliptic lat and long',declip(2:3)*rad2deg
 
       !if ecliptic latitude too large or too small to find MBA's cycle
!       if(declip(1,2)*rad2deg.lt.-80.d0) cycle
!       if(declip(1,2)*rad2deg.gt.80.d0) cycle
       
       if(simplearth) then
       !Heliocentric ecliptic position of the Earth through 2body propagation
        call epoch2ma(nearth,t0e,me0,tobs(k),oearth(6),.false.)  
        call kema2rv(oearth(1:6),rve(1:6))
        
       else
        !Heliocentric ecliptic position of the Earth through JPL ephemeris
        tobsjd=tobs(k)+2400000._wp
        call getrvearth(tobsjd,'ECLIPJ2000','solar system barycenter',rve(1:6))

        call rotearth(tobsjd,lsstxyz0,'ECLIPJ2000',lsstxyz)
        call rotearth(tobsjd-0.01_wp,lsstxyz0,'ECLIPJ2000',lsstxyz1)
    
        rve(1:3)=(rve(1:3)+lsstxyz(1:3))/aukm
        rve(4:6)=(rve(4:6)+(lsstxyz(1:3)-lsstxyz1(1:3))/(0.01_wp*86400._wp))/vau2vkm
       end if

       !Heliocentric ecliptic position of asteroids through 2body propagation 
       !Check if ephemeris has to be recalculated
       if(abs(tobs(k)-teph).gt.tenew) then
!       
!/////////////
!determine new positions for asteroids
        visible(:)=.false.
!$OMP PARALLEL DO PRIVATE(i)
       do i=1,nsmpl
         call epoch2ma(oast(i,7),oast(i,8),ma0(i),tobs(k),oast(i,6),.false.)  
         call kema2rv(oast(i,1:6),rva(i,1:6))
         call dist3d(rva(i,1:3),rve(1:3),rea(i,1:3),dea(i),dsa(i),dse(i),phase(i),solarelon(i))
         !Muinonen
         !oast(i,10)=Vismag(oast(i,9),phase,dsa,dea,g12)
         !JPL (photG=G)
         oast(i,10)=VismagHG(oast(i,9),phase(i),dsa(i),dea(i),photG)
         
         !normalize asteroid Earth vector for faster FOV check
         rea(i,1:3)=rea(i,1:3)/dea(i)
         
         !is the object visible?
         !if(oast(i,10).lt.fsdmax(1).and.abs(angv1v2n(3,rea(i,1:3),rve(1:3)))*rad2deg.lt.140.d0) then
         if(oast(i,10).lt.fsdmax(1)) then
           visible(i)=.true.
         end if !visible
       end do !nsmpl
!$OMP END PARALLEL DO  

       nvis=0
       do i=1,nsmpl
        if(visible(i)) then
            nvis=nvis+1
            visid(nvis)=i

        end if
       end do

       teph=tobs(k)
       end if !ephemeris update
    
!$OMP PARALLEL DO PRIVATE(i,j,obsang,rean,rvrea,az,el)
       !check whether or not asteroids are in FOV only for objects that are bright enough. 
do j=1,nvis
   i=visid(j)
   recorded(i)=.false.

!rough Vmag cut
   if(oast(i,10).lt.fsdp1) then       
     !see whether object is in the extended FOV (x2) via angle between pointing direction dx and Earth asteroid vector
     obsang=angv1v2n(3,rea(i,1:3),dx)
     
     !if in FOV x 2
       if(obsang.le.fovoax2) then
         !update ephemeris
         call epoch2ma(oast(i,7),oast(i,8),ma0(i),tobs(k),oast(i,6),.false.)  
         call kema2rv(oast(i,1:6),rva(i,1:6))
         !update distances and directions
         call dist3d(rva(i,1:3),rve(1:3),rea(i,1:3),dea(i),dsa(i),dse(i),phase(i),solarelon(i))
        
         !aberration correction
         if(abcorr) then
          call epoch2ma(oast(i,7),oast(i,8),ma0(i),tobs(k)-lighttraveltime(dea(i)),oast(i,6),.false.)  
          call kema2rv(oast(i,1:6),rva(i,1:6))
          call dist3d(rva(i,1:3),rve(1:3),rea(i,1:3),dea(i),dsa(i),dse(i),phase(i),solarelon(i))
         end if
         
         !Vmag: Vismag(H,phase,rs,re,g12)
         !Muinonen Vmag=Vismag(oast(i,9),phase,dsa,dea,g12)
         !JPL:
         Vmag(i)=VismagHG(oast(i,9),phase(i),dsa(i),dea(i),photG)
         
         !fine cut on Vmag < five sigma depth of LSST at that moment
         if(Vmag(i).le.fivesigdepth(k)) then
          ! renormalize Earth asteroid vector
          rea(i,1:3)=rea(i,1:3)/dea(i)
          ! get ecliptic longitude and latitude
          call cart2lat(rea(i,1:3),declip(i,1:3))
          rvrea(1:6)=rva(i,1:6)-rve(1:6)
          call rvhe2icrf(rvrea)
          call rvicrf2dradec(rvrea,rradec(i,:),drradec(i,:))
    
          call radec2azel(rradec(i,2),rradec(i,3),lst(k),lsstlat,az,el)
          call azel2camxy(az,el,azf(k),elf(k),rot(k),camx(i),camy(i))

          if(cam) then
          ! check whether the object actually is in the field of view (use camera)    
            call camcheck(camx(i),camy(i),recorded(i))
          else
            !see if it is still there
            obsang=angv1v2n(3,rea(i,1:3),dx)
            if(obsang.le.fovoa) then  
                recorded(i)=.true.
            end if
          end if

        end if !fine cut on Vmag
       end if !rough cut on FOV
      end if !rough cut on Vmag 
    end do !nvis
!$OMP END PARALLEL DO
    
!/////////      
!OUTPUT
 
      do j=1,nvis
           i=visid(j)
           if(recorded(i)) then
             write(un,fmt) &
              !objID, time of obs,     RightAscension [deg], Declination [deg], dRA/dt ["/hr], dDEC/dt ["/hr]
              trim(orbname(i)),tobs(k),rradec(i,2:3)*rad2deg,drradec(i,2:3)*rad2deg*3600.d0/24.d0,&
              !phase angle [deg], solar elongation [deg], distance Sun asteroid [au], distance Earth asteroid [au], distance Sun Earth [au]
              phase(i)*rad2deg,solarelon(i)*rad2deg,dsa(i),dea(i),dse(i),&
              !exposure time [s], five sigma depth [mag], H [mag], V [mag], orbit #, Ecliptic longitude and latitude [deg] 
              texp_s(k),fivesigdepth(k),oast(i,9),Vmag(i),orbid(i),declip(i,3)*rad2deg,declip(i,2)*rad2deg, &
              !heliocentric x,y,z of asteroid, FOV RA DEC [deg], FOV AZ EL [deg], 
              rva(i,1:3),raf(k)*rad2deg,decf(k)*rad2deg,azf(k)*rad2deg,elf(k)*rad2deg, &
              !camera x,y [m], seeing FWHM, filter
              camx(i),camy(i),seeingFwhmEff(k),filter(k)
            end if !recorded
     end do !nvis  
    end do !nobs
    close(un)

!/////////      
!END PROGRAM 
end program


