module spice_m
 use kind_m
 implicit none
!*************************************************************************************************
!   INTERFACES TO SPICE LIBARY ROUTINES
!*************************************************************************************************
 interface 
 
 subroutine FURNSH(spkname)
  character(len=*)::spkname
 end subroutine
 
 subroutine str2et(tstr,et)
  double precision ET
  character(len=*),intent(in)::tstr
 end subroutine
 
 subroutine spkezr(nstr,et,reffstr,abcorstr,centerstr,state,lt)
    character(len=*),intent(in)::nstr
    character(len=*),intent(in)::reffstr
    character(len=*),intent(in)::abcorstr
    character(len=*),intent(in)::centerstr
    double precision::et,state(1:6),lt
end subroutine

 end interface
!*************************************************************************************************
!   MODULE CONSTANTS (INCLUDES SPICE KERNEL LOCATION, EPHEMERIS AND PERTUBER LISTS)
!*************************************************************************************************
 PUBLIC
    !path to spice kernels on local machine
    character(len=*),parameter::pathtospk='/Users/eggl/Documents/src/spice/toolkit/eph/'
    
    !name of binary spk files to be loaded (includes tls files)
    character(len=*),parameter::spkde430='de430.bsp'            !JPL DIGITAL EPHEMERIS
    
    character(len=*),parameter::spkast='ast343de430.bsp'  ! 343 ASTEROID PERTURBERS USED IN EPHEMERIS
    
    character(len=*),parameter::spkt='naif0012.tls'          ! NAIF SPK COLLECTION INCLUDING LEAP SECONDS KERNEL
    character(len=*),parameter::spkgm='gm_de431.tpc'       !GM VALUES OF BODIES USED IN THE DE 431 EPHEMERIS
    character(len=*),parameter::spkrot='pck00010.tpc'      !PCK orientaion file
    
    character(len=*),parameter::spke1='earth_070425_370426_predict.bpc'      !PCK orientaion file
    character(len=*),parameter::spke2='earth_latest_high_prec.bpc'      !PCK orientaion file
    
    character(len=*),parameter::spkdidymos='a65803_1600-2500_jpl134.bsp'      !DIDYMOS SOLUMTION JPL 134
    
    !PLANETS
    integer,parameter::neph=11
    
    character(len=3),parameter::earthid='399'
    character(len=3),dimension(1:neph),parameter::planets=(/'009','001','301', &
                                                                         '008','007','004','006','002','399','005','010'/)     
    !PLUTO BARYC, MERCURY BARYC, MOON, NEPTUNE BC, URANUS BC, MARS BC, SATURN BC, VENUS BC, EARTH, JUPITER BC, SUN
    
    real(kind=wp),dimension(1:neph)::gmplanets
    
    !PERTURBERS
    integer,parameter::npert=4
    
    character(len=7),dimension(1:4),parameter::cpvh=(/'2000001', '2000002','2000004','2000010'/) 
    !CERES PALLAS VESTA HYGIEA
     
    
    character(len=7),dimension(1:16),parameter::big16=(/'2000001', &    !CERES
    '2000002','2000004','2000010', &                                       !PALLAS VESTA HYGIEA
    '2000003', '2000006','2000007','2000015','2000016', &     !JUNO, HEBE, IRIS, EUNOMIA, PSYCHE
    '2000029','2000052','2000065','2000087','2000088','2000511','2000704'/)  ! AMPHITRITE, EUROPA, CYBELE, SYLVIA, THISBE, DAVIDA, INTERAMNIA
    
    real(kind=wp),dimension(1:neph)::gmpert

    
    
!*************************************************************************************************
!   MODULE SUBROUTINES
!*************************************************************************************************  

    public::load_spk              ! LOAD ALL SPK FILES GIVEN IN SPK1-5 ABOVE
    public::unload_spk          ! UNLOAD ALL SPK FILES

    public::getgmplanets            ! GET GM OF PLANETS FROM SPICE FILE 
    public::getgmpert           ! GET GM OF PERTURBERS FROM SPICE FILE
    
    public::getreph                ! GET EME J2000 COORDINATES OF ALL PLANETS AND THE EARTH MOON
    public::getrveph              ! GET EME J2000 COORDINATES AND VELOCITIES OF ALL PLANETS AND THE EART MOON
    
    public::getrearth                ! GET EME J2000 COORDINATES OF THE EARTH
    public::getrvearth              ! GET EME J2000 COORDINATES AND VELOCITIES OF THE EARTH
    
    public::getrpert                ! GET EME J2000 COORDIANTES OF NPERT PERTURBERS
    public::getrvpert               ! GET EME J2000 COORDINATES AND VELOCITIES OF NPERT PERTURBERS
    
    public::rotearth                ! Rotate a topocentric vector with the Earth
    
    public::test_spk                ! TEST THE ROUNTINES IN THIS MODULE (UNDER DEVELOPMENT)
    
    

    contains

!****************************************************************************************
    subroutine load_spk
    implicit none
    !load spice kernels with ascending priority (priority spk4>spk3>...>spk1)
            CALL FURNSH ( trim(pathtospk)//trim(spkde430))
            CALL FURNSH ( trim(pathtospk)//trim(spkast))
            CALL FURNSH ( trim(pathtospk)//trim(spkt))
            CALL FURNSH ( trim(pathtospk)//trim(spkgm))
            CALL FURNSH ( trim(pathtospk)//trim(spkrot))
            CALL FURNSH ( trim(pathtospk)//trim(spke1))
            CALL FURNSH ( trim(pathtospk)//trim(spke2))
            CALL FURNSH ( trim(pathtospk)//trim(spkdidymos))
            
            
    return
    end subroutine
    
        !****************************************************************************************
    subroutine unload_spk
    implicit none
    !load spice kernels with ascending priority (priority spk4>spk3>...>spk1)
            CALL UNLOAD ( trim(pathtospk)//trim(spkde430))
            CALL UNLOAD ( trim(pathtospk)//trim(spkast))
            CALL UNLOAD ( trim(pathtospk)//trim(spkt))
            CALL UNLOAD ( trim(pathtospk)//trim(spkgm))
            CALL UNLOAD ( trim(pathtospk)//trim(spkrot))
            CALL UNLOAD ( trim(pathtospk)//trim(spke1))
            CALL UNLOAD ( trim(pathtospk)//trim(spke2))
            CALL UNLOAD ( trim(pathtospk)//trim(spkdidymos))
            
    return
    end subroutine
       !****************************************************************************************
    
    subroutine load_spk_earth
    implicit none
    !load spice kernels with ascending priority used for Earth orientation and position (priority spk4>spk3>...>spk1)
            CALL FURNSH ( trim(pathtospk)//trim(spkde430))
            CALL FURNSH ( trim(pathtospk)//trim(spkt))
            CALL FURNSH ( trim(pathtospk)//trim(spkrot))
            CALL FURNSH ( trim(pathtospk)//trim(spke1))
            CALL FURNSH ( trim(pathtospk)//trim(spke2))     
    return
    end subroutine
    !****************************************************************************************
    subroutine unload_spk_earth
    implicit none
    !load spice kernels with ascending priority (priority spk4>spk3>...>spk1)
            CALL UNLOAD ( trim(pathtospk)//trim(spkde430))
            CALL UNLOAD ( trim(pathtospk)//trim(spkt))
            CALL UNLOAD ( trim(pathtospk)//trim(spkrot))
            CALL UNLOAD ( trim(pathtospk)//trim(spke1))
            CALL UNLOAD ( trim(pathtospk)//trim(spke2))
            return
    end subroutine
    !****************************************************************************************
    subroutine rotearth(t,r,frame,rrot)
    !use spice routines to get the corect geocentric position for a station at coordinates r on the Earth
    
    !Important available frames are:
    !'J2000' (RADEC)
    !'ECLIPJ2000' (ecliptic)
    
    implicit none
    
    !INPUT
    real(kind=wp),intent(in)::t !time [JD]
    real(kind=wp),dimension(1:3),intent(in)::r !position on the Earth (WGS84)
    character(len=*),intent(in)::frame !SPK frame and center of coordinate origin
    
    !OUTPUT
    real(kind=wp),dimension(1:3),intent(out)::rrot !rotated position on the Earth

    real*8::LT,ET,ROTM(1:3,1:3)
   
    !integer::i
   
    call jd2et(t,ET)
    CALL PXFORM ('ITRF93', trim(frame), ET, ROTM )
    
    rrot=matmul(ROTM,r)
    
!    write(*,*)r
!    
!    do i=1,3
!     write(*,*)ROTM(i,:)
!    end do
!    
!    write(*,*)rrot
    
    return
    end subroutine
    
    !****************************************************************************************
    subroutine getrearth(t,frame,center,reph)
    !use spice routines to get J2000 position for the Earth
    
    !Important available frames are:
    !'J2000' (RADEC)
    !'ECLIPJ2000' (ecliptic)
    
    implicit none
    
    real(kind=wp),intent(in)::t !time [JD]
    character(len=*),intent(in)::frame,center !SPK frame and center of coordinate origin
    real(kind=wp),dimension(1:3),intent(out)::reph !position of the Earth

    real*8::LT,ET
   
    call jd2et(t,ET)
    CALL SPKPOS (earthid, et , trim(frame), 'NONE', trim(center) ,reph(1:3), LT  )
    return
    end subroutine
    !****************************************************************************************
    subroutine getrvearth(t,frame,center,rveph)
    !use spice routines to get J2000 position for the Earth
    
     !Important available frames are:
    !'J2000' (RADEC)
    !'ECLIPJ2000' (ecliptic)
    
    implicit none
    
    real(kind=wp),intent(in)::t !time [JD]
    character(len=*),intent(in)::frame,center !SPK frame and center of coordinate origin
    real(kind=wp),dimension(1:6),intent(out)::rveph !positions and velocities of the Earth

    real*8::LT,ET
   
    call jd2et(t,ET)
    CALL SPKEZR (earthid, et , trim(frame), 'NONE', trim(center) ,rveph(1:6), LT  )
    
    return
    end subroutine
    !****************************************************************************************
    subroutine getreph(t,reph)
    !use spice routines to get J2000 position for the planets
    implicit none
    
    real(kind=wp),intent(in)::t
    real(kind=wp),dimension(1:neph,1:3),intent(out)::reph
    
    integer*4::i
    real*8::LT, ET
    
    call jd2et(t,ET)

    write(*,*)'t',t,'ET',ET
    
    do i=1,neph
     CALL SPKPOS (planets(i), et , 'J2000', 'NONE', 'solar system barycenter',reph(i,1:3), LT  )
    end do
    return
    end subroutine
    
     !****************************************************************************************
    subroutine getrveph(t,rveph)
    !use spice routines to get J2000 position for the planets
    implicit none
    
    real(kind=wp),intent(in)::t
    real(kind=wp),dimension(1:neph,1:6),intent(out)::rveph
    
    integer*4::i
    real*8::LT,ET
    
    call jd2et(t,ET)
    
    do i=1,neph
     CALL SPKEZR (planets(i), et , 'J2000', 'NONE', 'solar system barycenter',rveph(i,1:6), LT  )
    end do
    return
    end subroutine
    
    !****************************************************************************************
    subroutine getrpert(t,rpert)
    !use spice routines to get J2000 position for the massive perturbers
    implicit none
    
    real(kind=wp),intent(in)::t
    real(kind=wp),dimension(1:npert,1:3),intent(out)::rpert
    
    integer*4::i
    real*8::LT,ET
    
    call jd2et(t,ET)
    
    do i=1,npert
     CALL SPKPOS (big16(i), ET , 'J2000', 'NONE', 'solar system barycenter',rpert(i,1:3), LT )
    end do
    return
    end subroutine
    
    
    !****************************************************************************************
    subroutine getrvpert(t,rvpert)
    !use spice routines to get J2000 position for the massive perturbers
    implicit none
    
    real(kind=wp),intent(in)::t
    real(kind=wp),dimension(1:npert,1:6),intent(out)::rvpert
    
    integer*4::i
    real*8::LT, ET
    
    call jd2et(t,ET)
      
    do i=1,npert
     CALL SPKPOS (big16(i), et , 'J2000', 'NONE', 'solar system barycenter',rvpert(i,1:6), LT )
    end do
    return
    end subroutine
    
    !****************************************************************************************
    subroutine getgmplanets
    !use spice routines to get GM for planets from corresponding spice kernel (GM_DE431.tpc)
    !read them into global constants
    implicit none
    
    !real(kind=wp),dimension(1:neph)::gmplanets
    integer*4::i,d
    
    do i=1,neph
     CALL BODVRD (planets(i), 'GM', 1, d, gmplanets(i))  
    end do 
    return
    end subroutine
    
    !****************************************************************************************
    subroutine getgmpert
        !use spice routines to get GM for planets from corresponding spice kernel (GM_DE431.tpc)
    !read them into global constants
    implicit none
    
    !real(kind=wp),dimension(1:neph)::gmplanets
    integer*4::i,d
    
    do i=1,npert
     CALL BODVRD (big16(i), 'GM', 1, d, gmpert(i))  
    end do 
    return
    end subroutine
    
    !***************************************************************************************
    subroutine ET2JD(et,jd)
    !routine to convert SPICE ephemeris time (TDB) given in seconds post epoch J2000 (2000 JAN 1, 12:00 TT) 
    ! to Julian Days
    implicit none
     real(kind=wp),intent(in)::et ![s]
     real(kind=wp),intent(out)::jd  ![julian days]
    
     jd=et/86400._wp+2451545._wp
     return
    end subroutine
    
    !***************************************************************************************
    subroutine JD2ET(jd,et)
    !routine to convert Julian Days to SPICE ephemeris time (TDB) given in seconds post epoch J2000 (2000 JAN 1, 12:00 TT) 
       implicit none
     real(kind=wp),intent(in)::jd   ![julian days]
     real(kind=wp),intent(out)::et ![s]
    
     et=(jd-2451545._wp)*86400._wp
     return
    end subroutine
    !***************************************************************************************
        

subroutine test_spk
            DOUBLE PRECISION   ET
            DOUBLE PRECISION   STATE (6),CERES(6),MOON(6),EARTH(6),Hygiea(6),LSST(3),LSSTrot(3)
            DOUBLE PRECISION   LT
            
            real(kind=wp),dimension(1:neph,1:3)::reph
            real(kind=wp),dimension(1:3)::dr    
            real(kind=wp)::gm,jd
                
            double precision, parameter::au=149597870.7d0
            integer*8::i,dd
            character(len=50)::iaun
            character(len=200)::hygieas='2000010'
            logical::fnd
            
            !Geocentric position of LSST
            LSST=(/1.81894d6, -5.20847d6, -3.19517d6/)
            
      !
      !  Convert UTC to ET 
      !
            CALL STR2ET ( '2022 OCT 5 15:03:06.272', ET )

            write(*,*)'ET:',ET
            call et2jd(et,jd)
            write(*,*)jd
            call jd2et(jd,et)
              write(*,*)'ET:',ET
            
               CALL STR2ET ( '2000 JAN 1 12:00:00.000', ET )

            write(*,*)'ET:',ET
            call et2jd(et,jd)
            write(*,*)jd
            call jd2et(jd,et)
            write(*,*)'ET:',ET
            
          !  STOP
      !
      !  Compute geometric state of MGS relative to Mars 
      !
      !      CALL SPKEZR ( '65803', ET, 'J2000', 'NONE', 'SUN',STATE, LT )
 
     open(21,file='tst.out',status='replace')
     do i=0,20,10
     
        jd=(et+real(i))/86400+2451545.0
        write(*,*)' '
        write(*,*)'ET',ET+real(i),'JD',jd
!       CALL SPKEZR ( '2065803', ET+real(i), 'J2000', 'NONE', 'solar system barycenter',STATE, LT )
!        write(*,*)'DIDYMOS ', jd, state(1:3)
!       CALL SPKEZR ( 'Ceres', ET+real(i), 'J2000', 'NONE', 'solar system barycenter',ceres, LT )
!        write(*,*)'CERES ',  jd, ceres(1:3)
!       CALL SPKEZR ( trim(hygieas), ET+real(i), 'J2000', 'NONE', 'solar system barycenter',hygiea, LT )
!          write(*,*)'HYGIEA ',  jd, ceres(1:3)
!       CALL SPKEZR ( 'MOON', ET+real(i), 'J2000', 'NONE', 'solar system barycenter',moon, LT )
!        write(*,*)'MOON ',  jd, moon(1:3)
       CALL SPKEZR ( 'EARTH', ET+real(i), 'J2000', 'NONE', 'solar system barycenter',earth(1:6), LT )
        write(*,*)'EARTH SPKEZR ', jd, earth(1:3),sqrt(dot_product(earth(1:3),earth(1:3)))/au
       CALL SPKPOS ('Earth', ET+real(i), 'J2000', 'NONE', 'solar system barycenter',earth(1:3), LT  )
          write(*,*)'EARTH SPKPOS ',  jd, earth(1:3),sqrt(dot_product(earth(1:3),earth(1:3)))/au
       Call getrearth (jd,'J2000','solar system barycenter',earth(1:3))               ! GET EME J2000 COORDINATES OF THE EARTH
        write(*,*)'GETREARTH  ',  jd, earth(1:3),sqrt(dot_product(earth(1:3),earth(1:3)))/au
      CALL getrvearth (jd,'J2000','solar system barycenter',earth(1:6))               ! GET EME J2000 COORDINATES AND VELOCITIES OF THE EARTH
          write(*,*)'GETRVEARTH  ',  jd, earth(1:6),sqrt(dot_product(earth(1:3),earth(1:3)))/au
          
          CALL SPKPOS ('Earth', ET+real(i), 'ECLIPJ2000', 'NONE', 'solar system barycenter',earth(1:3), LT  )
          write(*,*)'EclipJ2000  ',  jd, earth(1:6),sqrt(dot_product(earth(1:3),earth(1:3)))/au
          
        CALL BODVRD (planets(3), 'GM', 1, dd, GM)  
        write(*,*)'IAU NR GM',GM,dd
        

     end do
        call rotearth(jd+0.5d0,LSST,'J2000',LSSTrot)
        write(*,*)'geocentric position of LSST',LSST,' at ', jd
        write(*,*)'geocentric position of rotated LSST rot',LSSTrot,' at ', jd+0.5
     
!     call getreph(2451545._wp,reph)
!     
!     do i=1,11
!     dd=0
!     read(planets(i),*)dd
!     call BODC2N(dd,iaun,fnd)
!      write(*,*)trim(iaun),fnd,reph(i,:)/au
!      end do
!      
!      do i=1,11
!          dd=0
!         read(planets(i),*)dd
!        call BODC2N(dd,iaun,fnd)
!       CALL BODVRD (planets(i), 'GM', 1, dd, GM)  
!       dr=(reph(i,:)-reph(6,:))
!       write(*,*)iaun,GM/dot_product(dr,dr)
!     end do
!     
!     call getgmplanets
!     do i=1,11
!      read(planets(i),*)dd
!        call BODC2N(dd,iaun,fnd)
!      write(*,*)iaun,gmplanets(i)
!      end do
     
     close(21)
     return
    !****************************************************************************************
end subroutine

end module
