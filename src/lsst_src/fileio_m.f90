module fileio_m
    use kind_m
    use const_m
    use orbdyn2_m, only:com2ke !used by readopenorb only
    
!  interface
!    subroutine get_file_stats(un,fname,scol,col,hl,ne,maxe,mine,sume,ave)
!    integer,intent(in)::un                      !unit name for opening file
!    integer,intent(in)::hl                      !number of header lines
!    integer,intent(in)::scol                    !number of columns selected for stats
!    integer,intent(in),dimension(scol)::col     !which columns are selected for stats e.g. scol=3, col=(/1,3,8/)
!    character(len=80)::fname                    !input file name
!
!    integer*8,intent(out)::ne                   !number of entries
!    real(kind=wp),dimension(scol),intent(out)::mine,maxe
!    real(kind=wp),dimension(scol),intent(out),optional::sume,ave
!    end subroutine
!  end interface


public:: read_hdr       ! read headerlines from a textfile
public:: get_file_stats ! get basic stats from a file with multiple columns
public:: file_cnt_ntrs  ! count number of readable entries in file
public:: readopenorb    ! read orbital elements from openorb (requires orbdyn2_m)

contains
!############################################################

subroutine readopenorb(un,nhdr,lenname,nsmpl,oast,orbid,orbname)
!READ OPEN ORB OE INPUT FORMAT
use orbdyn2_m, only:com2ke
implicit none

integer,intent(in)::un,nhdr,lenname
integer*8,intent(in)::nsmpl

integer*8,intent(out)::orbid(1:nsmpl)
real(kind=wp),dimension(1:nsmpl,1:10),intent(out)::oast
character(len=lenname),dimension(1:nsmpl),intent(out)::orbname
           
integer*8::i,eof
real(kind=wp)::oe(1:6),ke(1:6),hmag,epoch,n
character(len=lenname)::dum1
character(len=3)::dum2

        
eof=0
if(nhdr.gt.0) then
 do i=1,nhdr
  read(un,*,iostat=eof)dum2
 end do
end if

do i=1,nsmpl
 read(un,*,iostat=eof)dum1,dum2,oe(1:6),hmag,epoch
 if(eof.ne.0) exit
 orbid(i)=i
 orbname(i)=dum1
 if(trim(dum2)=='COM') then
    call com2ke(oe(1:6),epoch,1._wp,0._wp,ke,n)
    oast(i,3:6)=ke(3:6)
    oast(i,1:2)=ke(1:2)
    oast(i,7)=n
    oast(i,8)=epoch
    oast(i,9)=hmag      
        
!    if(epoch.gt.2400000) then
!        oast(i,8)=epoch-2400000.d0
!    else
!        oast(i,8)=epoch
!    end if
 else
  write(*,*)'INPUT ERROR:orbit element format unknown'
 end if
end do
 
return
end subroutine

!***********************************************************************
subroutine get_file_stats(un,fname,scol,col,hl,ne,maxe,mine,sume,ave)
implicit none
integer,intent(in)::un                      !unit name for opening file
integer,intent(in)::hl                      !number of header lines
integer,intent(in)::scol                    !number of columns selected for stats
integer,intent(in),dimension(scol)::col     !which columns are selected for stats e.g. scol=3, col=(/1,3,8/)
character(len=80)::fname                    !input file name

integer*8,intent(out)::ne                   !number of entries
real*8,dimension(scol),intent(out)::mine,maxe
real*8,dimension(scol),intent(out),optional::sume,ave
real*8::kc,ky,kt                               !Kahan summation variables

integer*8::i,j
integer::eof
character::dum
character,dimension(:),allocatable::dum2
character(len=80)::hdr
real*8::x

maxe(:)=-1.d99
mine(:)=1.d99
sume(:)=0.d0
ave(:)=0.d0

open(un,file=trim(fname),status='old')

call read_hdr(un,hl,hdr)
ne=0
eof=0

do while(eof.eq.0) 
 read(un,*,iostat=eof)dum
 if(eof.ne.0) exit
 ne=ne+1
end do

do j=1,scol
 rewind(un)
 kc=0.d0
 kt=0.d0
 ky=0.d0
 call read_hdr(un,hl,hdr)
 
 if (col(j).eq.1) then
  do i=1,ne
   read(un,*)x
   !Min/Max
   if(x.gt.maxe(j)) then 
     maxe(j)=x
   end if
   if(x.lt.mine(j)) then
     mine(j)=x
   end if
if(present(sume)) then   
   !Sum
   ky=x-kc
   kt=sume(j)+ky
   kc=(kt-sume(j))-ky
   sume(j)=kt
end if
  end do
if(present(sume).and.present(ave)) then  
   !Average
   ave(j)=sume(j)/real(ne)
end if
 else
  allocate(dum2(1:col(j)-1))
  do i=1,ne
   read(un,*)dum2,x
   !Min/Max
   if(x.gt.maxe(j)) then 
     maxe(j)=x
   end if
   if(x.lt.mine(j)) then
     mine(j)=x
   end if
if(present(sume)) then    
   !Sum
   ky=x-kc
   kt=sume(j)+ky
   kc=(kt-sume(j))-ky
   sume(j)=kt
end if
  end do
if(present(sume).and.present(ave)) then  
   !Average
   ave(j)=sume(j)/real(ne)
end if
   deallocate(dum2)
 end if
 
end do !j number of columns to analyse
close(un)
return
end subroutine

!***********************************************************************
 
subroutine read_hdr(un,hl,hdr)
implicit none
integer,intent(in)::un !file unit number
integer,intent(in)::hl !number of header lines to read
character(len=80),intent(out)::hdr

integer::i

if (hl>0) then
 do i=1,hl
  read(un,'(A)') hdr
 end do
end if
return
end subroutine
!***********************************************************************
subroutine file_cnt_ntrs(un,fname,hl,ne)
!count number of entries in file
implicit none
integer,intent(in)::un                      !unit name for opening file
integer,intent(in)::hl                      !number of header lines
character(len=80)::fname                    !input file name

integer*8,intent(out)::ne               !number of entries
character::dum
character(len=80)::hdr
integer::eof

open(un,file=fname,status='old')

call read_hdr(un,hl,hdr)
ne=0
eof=0

do while(eof.eq.0) 
 read(un,*,iostat=eof)dum
 if(eof.ne.0) exit
 ne=ne+1
end do

rewind(un)

return
end subroutine
!############################################################
end module