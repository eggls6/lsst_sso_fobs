module photometry_m
 use kind_m
 use const_m
 use trafo3_m, only:angv1v2,angv1v2n
!contains subroutines and functions relevant to asteroid photometry such
!as Vmag calulations

!implicit none
!
!PUBLIC
!!    real(kind=wp)::g1,g2,phi1,phi2,phi3
!
!public::g1
!public::g2
!public::phi1
!public::phi2
!public::phi3
!public::vred
!public::vismag
!public::vismaghg

!interface
! function angv1v2(d,v1,v2)
!!calculate angle between two vectors [rad]
!    integer,intent(in)::d !dimension
!    real(kind=wp),intent(in),dimension(d)::v1,v2 !two non unit vectors
!    real(kind=wp)::angv1v2 !angle between v1 and v2 in rad
! end function
!
! function angv1v2n(d,v1,v2)
!!calculate angle between two normalized vectors [rad]
!    integer,intent(in)::d !dimension
!    real(kind=wp),intent(in),dimension(d)::v1,v2 !two non unit vectors
!    real(kind=wp)::angv1v2n !angle between v1 and v2 in rad
! end function

! function phi1(y)
!  real(kind=wp),intent(in)::y !phase
!  real(kind=wp)::phi1
! end function
!
! function phi2(y)
!  real(kind=wp),intent(in)::y !phase
!  real(kind=wp)::phi2
! end function
!
! function phi3(phase)
!  real(kind=wp),intent(in)::phase
!  real(kind=wp)::phi3
! end function
!
! function g1(g12)
!  real(kind=wp),intent(in)::g12
!  real(kind=wp)::g1
! end function
!
! function g2(g12)
!  real(kind=wp),intent(in)::g12
!  real(kind=wp)::g2
! end function
! 
! function Vred(H,phase,g12)
!!reduced magnitude (phase angle correction)
!!following Muinonen et al. (2010)
!  real(kind=wp),intent(in)::H,phase,g12
!  real(kind=wp)::Vred
! end function
!
! function VismagHG(H,phase,rs,re,g)
!!Visual magnitude using HG model
!  real(kind=wp),intent(in)::H,phase,rs,re,g
!  real(kind=wp)::VismagHG
! end function
!
! function Vismag(H,phase,rs,re,g12)
!!Visual magnitude using G12 model (Muinonen et al., 2010) 
!  real(kind=wp),intent(in)::H,phase,rs,re,g12
!  real(kind=wp)::Vismag
! end function
!
!end interface

!PUBLIC

contains
!############################################################
!***********************************************
function phi1(y)
!verified
 implicit none
 real(kind=wp),intent(in)::y !phase
 real(kind=wp)::phi1
! real(kind=wp),parameter::pi=3.1415926535897932384626433d0
! phi1=1.d0-6.d0*phase/pi ! valid for phase angles < 30 deg
!Mathematica derived fit from Muinonen et al. (2010) paper data Table 3
 if(y.lt.5.d0/6.d0*pi) then
  phi1=0.993847-2.14146d0*y+2.30312d0*y**2-1.3907d0*y**3 &
      +0.43383d0*y**4-0.0536246*y**5
 else
  phi1=0.d0
 end if
 
return
end function
!***********************************************
function phi2(y)
!verified
 implicit none
 real(kind=wp),intent(in)::y !phase
 real(kind=wp)::phi2
! real(kind=wp),parameter::pi=3.1415926535897932384626433d0
! phi2=1.d0-9.d0/5.d0*phase/pi
 if(y.lt.5.d0/6.d0*pi) then
   phi2=1.02385d0-0.7355d0*y-0.19773d0*y**2+0.385224d0*y**3&
      -0.157956d0*y**4+0.0224944d0*y**5
  else
   phi2=0.d0
  end if
return
end function
!***********************************************
function phi3(phase)
!verified
 implicit none
 real(kind=wp),intent(in)::phase
 real(kind=wp)::phi3
 !real(kind=wp),parameter::pi=3.1415926535897932384626433d0
 if(phase<pi/6.d0) then
  phi3=exp(-4.d0*pi*tan(0.5d0*phase)**(2.d0/3.d0))
 else
  phi3=0.d0
 end if
return
end function
!***********************************************
function g1(g12)
!verified
 implicit none
  real(kind=wp)::g12
  real(kind=wp)::g1
  
  if(g12<0.2) then
   g1=0.7527d0*g12+0.06164d0 
  else
   g1=0.9529d0*g12+0.02162d0
  end if
return
end function
!***********************************************
function g2(g12)
!verified
 implicit none
  real(kind=wp)::g12
  real(kind=wp)::g2
  
  if(g12<0.2) then
   g2=-0.9612d0*g12+0.6270d0 
  else
   g2=-0.6125d0*g12+0.5572d0
  end if
return
end function
!**********************************************
function Vred(H,phase,g12)
!reduced magnitude (phase angle correction)
!following Muinonen et al. (2010)
implicit none
 real(kind=wp),intent(in)::H,phase,g12
! real(kind=wp)::phi1,phi2,phi3
! real(kind=wp)::g1,g2
 real(kind=wp),dimension(1:3)::g,phi

 real(kind=wp)::Vred

 g(1)=g1(g12) 
 g(2)=g2(g12)
 g(3)=1.d0-g(1)-g(2)

 phi(1)=phi1(phase)
 phi(2)=phi2(phase)
 phi(3)=phi3(phase)

 Vred=H-2.5d0*log10(dot_product(g,phi))
 return
end function
!***********************************************
function VismagHG(H,phase,rs,re,g)
!verified
implicit none
 real(kind=wp),intent(in)::H,phase,rs,re,g

 real(kind=wp)::VismagHG
 real(kind=wp)::q,phi(1:2)
 real(kind=wp),dimension(1:2),parameter::a=(/-3.332d0,-1.862d0/),b=(/0.631d0,1.218d0/)
 real(kind=wp),parameter::phaselim=2.094395102393195d0
 
 integer::i
 
 if(phase.gt.phaselim) then
    VismagHG=H+5.d0*log10(rs*re)
 else
    do i=1,2
     phi(i)=exp(a(i)*tan(phase*0.5d0)**b(i))
    end do
    q=(1.d0-g)*phi(1)+g*phi(2)
    VismagHG=H+5.d0*log10(rs*re)-2.5d0*log10(q)
end if
return
end function
!***********************************************
function Vismag(H,phase,rs,re,g12)
!verified
 implicit none
 real(kind=wp),intent(in)::H,phase,rs,re,g12
 !real(kind=wp)::Vflux,Vred
 real(kind=wp)::Vismag
!        Hflux=1.d0
!        Vflux=1.d0/(rs*rs)/(re*re)
!        Vmag=H-2.5d0*log10(Vflux/Hflux)
!        Vismag=Vred(H,phase,g12)-2.5d0*log10(Vflux)
!        Vismag=Vred(H,phase,g12)+5.d0*(log10(rs)+log10(re))
         Vismag=Vred(H,phase,g12)+5.d0*log10(rs*re)
return
end function
!***********************************************

!############################################################
end module