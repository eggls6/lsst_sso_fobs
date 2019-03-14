module camera_m
!//////////////////////////////////////////////////////////////////////////////////
! LSST camera module used to check whether RA DEC observations fall on the camera
!
! written by S.Eggl 20181201
!//////////////////////////////////////////////////////////////////////////////////

        use kind_m
        use const_m
        use trafo3_m
        
!       public::gapcheck
        public::camcheck
        public::azel2camxy

        public::radec2altaz
        public::altaz2radec

        Public
        
        real(kind=wp),parameter::focalLen=10.31_wp ![m],
        real(kind=wp),parameter::fovang=0.0305433 ! [rad], half of the 3.5 [deg] = 0.061 [rad] FOV diameter
        real(kind=wp),parameter::sensord=0.629998_wp ![m], sensor diameter
        real(kind=wp),parameter::sensorr=0.314999_wp ![m], sensor diameter
        real(kind=wp),parameter::raftd=0.126_wp ![m], raft diameter
        real(kind=wp),parameter::fillfactor=0.9089d0 !Veres et al. 2017

contains
!************************************************************************************************************************  
subroutine camcheck(x,y,recorded)
!check whether x,y coordinates of observed objects actually fall on the detector (includes statistical fill factor)
implicit none

real(kind=wp),intent(in)::x,y

logical,intent(out)::recorded

real(kind=wp)::fovrand

recorded=.false.

!CHECK WHETHER X,Y coordinates are on LSST SENSOR CONFIGURATION
    if(abs(x).lt.sensorr.and.abs(y).lt.sensorr) then
     if(abs(x).gt.(sensorr-raftd).and.abs(y).gt.(sensorr-raftd)) then
     else
!     !check for fill factor
!        call random_number(fovrand)
!          if(fovrand.lt.fillfactor) then 
            recorded=.true.
!          end if
     end if
    end if
    
return
end subroutine
    
!************************************************************************************************************************    
subroutine azel2camxy(az,el,azcen,elcen,rot,x,y)
! rotate a point with (az,el) direction around the center of the FOV (azcen,elcen) for an angle rot (clockwise!)
! the new position is given in LSST sensor coordinates (x,y) 
! all angles in [rad]!
! 
implicit none

real*8::az,el,azcen,elcen,rot
real*8::x,y

real*8,dimension(1:3)::p,px,pt,pxdx,fov,fovx,dx,x0
real*8::crot,srot
real*8::m(1:2,1:2)

x0=(/1.d0,0.d0,0.d0/)
!align system relative to FOV axis (azcen,elcen) pointing towards (1,0,0)
!daz=az-azcen
!if(daz.gt.pi) then
!    if(az.gt.azcen) then
!      daz=-(pix2-daz)
!    else
!      daz=pix2-daz
!    end if
!end if
!del=el-elcen

! not sure what I wanted here...
!cdaz=cos(daz)
!sdaz=sin(daz)
!cdel=cos(del)
!sdel=sin(del)
fov=(/1.d0,elcen,azcen/)
call lat2cart(fov,fovx)

p=(/1.d0,el,az/)
call lat2cart(p,px)

!write(*,*)'obsang ',angv1v2(3,fovx,px)*rad2deg

!write(*,*)'fovx',fovx
call roty(-elcen,fovx)
call rotz(-azcen,fovx)

!write(*,*)'fovx rot',fovx
!correct for small displacement to desired FOV normal at (1,0,0)
dx=x0-fovx
!write(*,*)'fovx rot +dx ',fovx+dx

call roty(-elcen,px)
call rotz(-azcen,px)
!write(*,*)'px',px
pxdx=px+dx
!write(*,*)'px+dx',pxdx,norm(pxdx-x0),norm(pxdx(2:3))

!write(*,*)'obsang rot',angv1v2(3,x0,pxdx)*rad2deg

!calculate relative vector of measurement to FOV axis
!p(1:3)=(/cdaz*cdel,sdaz*cdel,sdel/)

!d=px-fovx
!project relative vector onto plane orthogonal to FOV axis
!dp=px-dot_product(px,fovx)*fovx

!since FOV normal is (1,0,0) only y and z coordinates count
!project angles to distances
!call cart2lat(px,pt)

!write(*,*)'cart2lat px',pt(2:3)*rad2deg

p(2:3)=focalLen*tan(pxdx(2:3))


!correct for FOV rotation wrt to az el 
crot=cos(-rot)
srot=sin(-rot)

m(1,1:2)=(/crot,-srot/)
m(2,1:2)=(/srot,crot/)

pt(2:3)=matmul(m,p(2:3))

x=pt(2)
y=pt(3)


! write(*,*)'az,azcen,daz',az,azcen,daz
! write(*,*)'el,elcen,del',el,elcen,del
! write(*,*)'p,rot',p(:),rot
! STOP
end subroutine


!subroutine camxy2azel(x,y,azcen,elcen,rot,az,el)
!!inverse function of azel2camxy. not working properly at the moment
!! 
!implicit none
!
!real*8::az,el,azcen,elcen,rot
!real*8::x,y,r2
!
!real*8,dimension(1:3)::p,pt
!!real*8::caz,cel,sel,saz,sazcen,cazcen,selcen,crot,srot,
!real*8::daz,del,cdaz,sdaz,cdel,sdel,crot,srot,ryz2,rxy
!real*8::m(1:2,1:2)
!
!pt(2)=x/focallen
!pt(3)=y/focallen
!
!crot=cos(rot)
!srot=sin(rot)
!
!m(1,1:2)=(/crot,-srot/)
!m(2,1:2)=(/srot,crot/)
!
!p(2:3)=matmul(m,pt(2:3))
!
!ryz2=dot_product(p(2:3),p(2:3))
!
!!point must lie on the unit sphere -> x=sqrt(1-y^2-z^2)
!p(1)=sqrt(1.d0-ryz2)
!
!daz=atan2(p(2),p(1))
!
!!del=asin(p(3))
!rxy=sqrt(dot_product(p(1:2),p(1:2)))
!if(rxy==0.d0) then
! del=0.d0
!else
! del=atan2(p(3),sqrt(dot_product(p(1:2),p(1:2))))
!end if
!
!!rotate system to be aligned with FOV axis pointing towards (1,0,0)
!az=daz+azcen
!if(az.lt.0.d0) then
!    az=az+pix2
!end if
!el=del+elcen
!
!return
!end subroutine

!******************************************

SUBROUTINE azel2radec(az,el,lst,obslat,ra,dec,radecvec)
!verified

!transform horizontal coordinates (azimut, elevation) = (az, el) into
! equatorial coordinates (right ascesion, declination) = (ra,dec) 
! all angles in [rad]! 
! local solar time also in [rad]!
! observer latitude also in [rad]
! see IAU_SOFA routine iau_AE2HD for details: http://www.iausofa.org/
!

      IMPLICIT NONE

      real*8::az,el,lst,obslat
      real*8::ra,dec
      real*8,dimension(1:3),optional::radecvec
      
      real*8:: SA, CA, SE, CE, SP, CP
      real*8:: X, Y, Z, R,HA

      SA = SIN(AZ)
      CA = COS(AZ)
      SE = SIN(EL)
      CE = COS(EL)
      SP = SIN(obslat)
      CP = COS(obslat)

!azel unit vector
      X = - CA*CE*SP + SE*CP
      Y = - SA*CE
      Z = CA*CE*CP + SE*SP


      R = SQRT(X*X + Y*Y)
      IF ( R.EQ.0.D0 ) THEN
         HA = 0.D0
      ELSE
         HA = ATAN2(Y,X)
      END IF
      
      RA=LST-HA
      
      if(RA.lt.0.d0) then
       RA=RA+pix2
      end if
      
      RA=modulo(RA,pix2)
      
      DEC = ATAN2(Z,R)

     if (present(radecvec)) then
                radecvec(:)=(/x,y,z/)
     end if
    return
end subroutine
!***********************************************************************
subroutine radec2azel(ra,dec,lst,obslat,az,el,azelvec)
!verified

!transform equatorial coordinates (right ascesion, declination) = (ra,dec) into 
!horizontal coordinates (azimut, elevation) = (az, el)
! all angles in [rad]! 
! local solar time also in [rad]!
! observer latitude also in [rad]
! see IAU_SOFA routine iau_HD2AE for details: http://www.iausofa.org/
!

implicit none
      real*8::ra,dec,lst,obslat
      real*8::az,el
        
      real*8,dimension(1:3),optional::azelvec

      real*8::ha,sh,ch,sd,cd,sp,cp
      real*8::x,y,z,a,r
      
    !local hour angle
      ha=lst-ra
        
      SH = SIN(HA)
      CH = COS(HA)
      SD = SIN(DEC)
      CD = COS(DEC)
      SP = SIN(OBSLAT)
      CP = COS(OBSLAT)

!*  Az,El unit vector.
      X = - CH*CD*SP + SD*CP
      Y = - SH*CD
      Z = CH*CD*CP + SD*SP

      R = SQRT(X*X + Y*Y)
      IF ( R.EQ.0.D0 ) THEN
         A = 0.D0
      ELSE
         A = ATAN2(Y,X)
      END IF
      IF ( A.LT.0.D0 ) A = A+pix2
      AZ = modulo(A,pix2)
      EL = ATAN2(Z,R)

      if (present(azelvec)) then
         azelvec(:)=(/x,y,z/)
      end if
        
        return
end subroutine 


end module

