module earth_m
!contains constants and subroutines pretaining to the EARTH
!written by S.Eggl 20190127

    use kind_m, only:wp
    use const_m, only:pix2,deg2rad

implicit none

public::latlonalt2r !calculate geocentric x,y,z positions of WGS84 position in latitude, longitude and altitude on the Earth

public
    real(kind=wp),parameter::WGS84a=6378.1370_wp       !World Geodetic System 84 reference ellipsoid semimajor axis a [km] 
    real(kind=wp),parameter::WGS84b=6356.752314245_wp  !World Geodetic System 84 reference ellipsoid semiminor axis b  [km] 
    

contains
!####################################################################################
subroutine latlonalt2r(lat,lon,alt,xyz) 
!verified
    implicit none
    real(kind=wp),intent(in)::lat,lon !Latitude (N), Longitude (E) WGS84 [deg]
    real(kind=wp),intent(in)::alt !Altitude above sea level WGS84 [m]
    
    real(kind=wp),dimension(1:3),intent(out)::xyz !position of X,Y,Z in WGS84 rectangular coordinates
    
    real(kind=wp)::nphi!nphi is the radius of curvature in the prime vertical
    real(kind=wp)::a2,b2,clat,clon,slat,slon,latrad,lonrad
    
    a2=WGS84a*WGS84a
    b2=WGS84b*WGS84b
    
    latrad=lat*deg2rad    
    lonrad=lon*deg2rad
    
    clat=cos(latrad)
    slat=sin(latrad)
    clon=cos(lonrad)
    slon=sin(lonrad)
    
    nphi=a2/sqrt(a2*clat*clat+b2*slat*slat)
    
    xyz(1)=(nphi+alt)*clat*clon
    xyz(2)=(nphi+alt)*clat*slon
    xyz(3)=(b2/a2*nphi+alt)*slat
    
    return
end subroutine

end module