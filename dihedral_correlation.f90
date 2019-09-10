module dihedral_correlation
implicit none
real(kind=8),allocatable,public::coord(:,:)
real(kind=8),allocatable,public::atom_coord(:,:)
integer(kind=8),allocatable::atomic_number(:),center_number(:)
integer(kind=8),parameter::natom=22,coor_dim=3
character(len=30)::rst_filename,log_filename

public sub_get_log_coord
public sub_get_rst_coord
public sub_dihedral

contains
subroutine sub_get_log_coord
implicit none
integer(kind=8)::j,atomic_type
character(len=110)::temp110

open(12,file='../log_data/'//log_filename)
do while(.true.)
  read(12,100) temp110
  100 format(a110)
  if(temp110(8:29).eq.'Stationary point found') then
    do while(.true.) 
      read(12,100) temp110
      if(temp110(26:45).eq.'Standard orientation') exit
    end do  
    read(12,*)
    read(12,*)
    read(12,*)
    read(12,*)
    do j=1,natom
      read(12,111) center_number(j),atomic_number(j),atomic_type,coord(1,j),coord(2,j),coord(3,j)
!      write(*,111) center_number(j),atomic_number(j),atomic_type,coord(1,j),coord(2,j),coord(3,j)
      111 format(5x,i2,10x,i1,11x,i1,7x,f9.6,3x,f9.6,3x,f9.6) 
    end do
  end if
  if(temp110(2:19).eq.'Normal termination') exit
end do
close(12)
end subroutine sub_get_log_coord

subroutine sub_get_rst_coord
implicit none
integer(kind=8)::i,istat2

!Get the molecular coordinates (x,y,z) from the rst_file bring by MM optimization.
open(14,file='../all_md_rst_file/'//rst_filename)
atom_coord=0.d0
do while(.true.)
  read(14,*,iostat=istat2)
  if(istat2.ne.0) exit
  read(14,*)
  read(14,444) (atom_coord(1,i),atom_coord(2,i),atom_coord(3,i),i=1,natom)
!  write(*,444) (atom_coord(1,i),atom_coord(1,i),atom_coord(1,i),i=1,natom)
  444 format(6f12.7)
end do 
close(14)
!if(allocated(atom_x)) deallocate(atom_x)
!if(allocated(atom_y)) deallocate(atom_y)
!if(allocated(atom_z)) deallocate(atom_z)
end subroutine sub_get_rst_coord

!The formula of calculating the dihedral.
subroutine sub_dihedral(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dihedral)
implicit none
real(kind=8),parameter::pai=atan(1.d0)*4.d0
real(kind=8),intent(in)::x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
real(kind=8),intent(out)::dihedral
real(kind=8)::ax,ay,az,bx,by,bz,cx,cy,cz
real(kind=8)::p1x,p1y,p1z,p2x,p2y,p2z
real(kind=8)::p1p2,p1_scalar,p2_scalar,cosphi
real(kind=8)::p2p1_cross_x,p2p1_cross_y,p2p1_cross_z,r2p2p1_cross,r2_scalar,sinphi
!http://www.bio.net/bionet/mm/molmodel/1999-March/001360.html
!              A             D
!           r1  \           /  r3
!                 B- - - -C 
!                    r2
! in which the points A(1), B(2), C(3) and D(4) are connected by the vectors
! r1(a), r2(b) and r3(c), we can obtain the torsion angle phi about r2 by applying:
!              p1= r1 x r2
!              p2= r2 x r3
!              p1 . p2 = |p1| |p2| cos (phi)
!              r2 . (p2 x p1) = |p1| |p2| |r2| sin(phi)
! and finally phi=atan[sin(phi)/cos(phi)] through a function like ATAN2 in Fortran.
ax=x2-x1
ay=y2-y1
az=z2-z1
bx=x3-x2
by=y3-y2
bz=z3-z2
cx=x4-x3
cy=y4-y3
cz=z4-z3
!    a x b                                           b x c   
!  |i  j  k |                                      |i  j  k |
!  |ax ay az|                                      |bx by bz|
!  |bx by bz|                                      |cx cy cz|
!  =(ay*bz-az*by,az*bx-ax*bz,ax*by-ay*bx)          =(by*cz-bz*cy,bz*cx-bx*cz,bx*cy-by*cx)
p1x=ay*bz-az*by
p1y=az*bx-ax*bz
p1z=ax*by-ay*bx
p2x=by*cz-bz*cy
p2y=bz*cx-bx*cz
p2z=bx*cy-by*cx
!p1 . p2
p1p2=p1x*p2x+p1y*p2y+p1z*p2z
p1_scalar=dsqrt(p1x*p1x+p1y*p1y+p1z*p1z)
p2_scalar=dsqrt(p2x*p2x+p2y*p2y+p2z*p2z)
cosphi=p1p2/(p1_scalar*p2_scalar)
!r2 . (p2 x p1)
!    p2 x p1
!   | i   j   k |
!   |p2x p2y p2z|
!   |p1x p1y p1z|
!   =(p2y*p1z-p2z*p1y,p2z*p1x-p2x*p1z,p2x*p1y-p2y*p1x)
p2p1_cross_x=p2y*p1z-p2z*p1y
p2p1_cross_y=p2z*p1x-p2x*p1z
p2p1_cross_z=p2x*p1y-p2y*p1x
r2p2p1_cross=bx*p2p1_cross_x+by*p2p1_cross_y+bz*p2p1_cross_z
r2_scalar=dsqrt(bx*bx+by*by+bz*bz)
sinphi=r2p2p1_cross/(p1_scalar*p2_scalar*r2_scalar)
!phi=atan[sin(phi)/cos(phi)]
dihedral=-atan2(sinphi,cosphi)*180.d0/pai
end subroutine sub_dihedral

end module dihedral_correlation
