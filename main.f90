program main
use dihedral_correlation
implicit none
integer(kind=8)::j,k 
integer(kind=8)::phi,psi
integer(kind=8),parameter::phi_num=25,psi_num=25 
integer(kind=8),parameter::phi_initial=-180,psi_initial=-180
integer(kind=8),parameter::increm=15      
character(len=4)::char_phi,char_psi 
real(kind=8),allocatable::qm_C_N_CA_C(:,:),qm_H_N_CA_C(:,:),qm_C_N_CA_CB(:,:),qm_H_N_CA_CB(:,:),qm_C_N_CA_H(:,:),qm_H_N_CA_H(:,:)
real(kind=8),allocatable::qm_N_CA_C_N(:,:),qm_HA_CA_C_N(:,:),qm_CB_CA_C_N(:,:),qm_N_CA_C_O(:,:),qm_HA_CA_C_O(:,:),qm_CB_CA_C_O(:,:)
real(kind=8),allocatable::mm_C_N_CA_C(:,:),mm_H_N_CA_C(:,:),mm_C_N_CA_CB(:,:),mm_H_N_CA_CB(:,:),mm_C_N_CA_H(:,:),mm_H_N_CA_H(:,:)
real(kind=8),allocatable::mm_N_CA_C_N(:,:),mm_HA_CA_C_N(:,:),mm_CB_CA_C_N(:,:),mm_N_CA_C_O(:,:),mm_HA_CA_C_O(:,:),mm_CB_CA_C_O(:,:)

allocate(qm_C_N_CA_C(phi_num,psi_num))
allocate(qm_H_N_CA_C(phi_num,psi_num))
allocate(qm_C_N_CA_CB(phi_num,psi_num))
allocate(qm_H_N_CA_CB(phi_num,psi_num))
allocate(qm_C_N_CA_H(phi_num,psi_num))
allocate(qm_H_N_CA_H(phi_num,psi_num))
allocate(qm_N_CA_C_N(phi_num,psi_num))
allocate(qm_HA_CA_C_N(phi_num,psi_num))
allocate(qm_CB_CA_C_N(phi_num,psi_num))
allocate(qm_N_CA_C_O(phi_num,psi_num))
allocate(qm_HA_CA_C_O(phi_num,psi_num))
allocate(qm_CB_CA_C_O(phi_num,psi_num))
allocate(mm_C_N_CA_C(phi_num,psi_num))
allocate(mm_H_N_CA_C(phi_num,psi_num))
allocate(mm_C_N_CA_CB(phi_num,psi_num))
allocate(mm_H_N_CA_CB(phi_num,psi_num))
allocate(mm_C_N_CA_H(phi_num,psi_num))
allocate(mm_H_N_CA_H(phi_num,psi_num))
allocate(mm_N_CA_C_N(phi_num,psi_num))
allocate(mm_HA_CA_C_N(phi_num,psi_num))
allocate(mm_CB_CA_C_N(phi_num,psi_num))
allocate(mm_N_CA_C_O(phi_num,psi_num))
allocate(mm_HA_CA_C_O(phi_num,psi_num))
allocate(mm_CB_CA_C_O(phi_num,psi_num))
allocate(atomic_number(natom)) 
allocate(center_number(natom))
allocate(coord(coor_dim,natom))
allocate(atom_coord(coor_dim,natom))

call system('mkdir qm_dih_cor_dat')
open(201,file='./qm_dih_cor_dat/qm_dih_cor01.dat')
open(202,file='./qm_dih_cor_dat/qm_dih_cor02.dat')
open(203,file='./qm_dih_cor_dat/qm_dih_cor03.dat')
open(204,file='./qm_dih_cor_dat/qm_dih_cor04.dat')
open(205,file='./qm_dih_cor_dat/qm_dih_cor05.dat')
open(206,file='./qm_dih_cor_dat/qm_dih_cor06.dat')
open(207,file='./qm_dih_cor_dat/qm_dih_cor07.dat')
open(208,file='./qm_dih_cor_dat/qm_dih_cor08.dat')
open(209,file='./qm_dih_cor_dat/qm_dih_cor09.dat')
open(210,file='./qm_dih_cor_dat/qm_dih_cor10.dat')
open(211,file='./qm_dih_cor_dat/qm_dih_cor11.dat')
open(212,file='./qm_dih_cor_dat/qm_dih_cor12.dat')
do j=1,phi_num
  phi=phi_initial+increm*(j-1) 
   write(char_phi,'(i4)') phi
  do k=1,psi_num  
    psi=psi_initial+increm*(k-1) 
     write(char_psi,'(i4)') psi
    log_filename='ala_dipeptide_'//trim(adjustl(char_phi))//'_'//trim(adjustl(char_psi))//'.log'
    rst_filename='ala_dipeptide_'//trim(adjustl(char_phi))//'_'//trim(adjustl(char_psi))//'.rst'
    call sub_get_log_coord
    call sub_get_rst_coord
    call sub_dihedral(coord(1,5),coord(2,5),coord(3,5),& 
                      coord(1,7),coord(2,7),coord(3,7),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      qm_C_N_CA_C(j,k))
     write(201,*) phi,psi,qm_C_N_CA_C(j,k)
    call sub_dihedral(coord(1,8),coord(2,8),coord(3,8),&
                      coord(1,7),coord(2,7),coord(3,7),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),& 
                      qm_H_N_CA_C(j,k))                   
     write(202,*) phi,psi,qm_H_N_CA_C(j,k)
    call sub_dihedral(coord(1,5),coord(2,5),coord(3,5),&                            
                      coord(1,7),coord(2,7),coord(3,7),&                            
                      coord(1,9),coord(2,9),coord(3,9),&                            
                      coord(1,11),coord(2,11),coord(3,11),&                         
                      qm_C_N_CA_CB(j,k))      
     write(203,*) phi,psi,qm_C_N_CA_CB(j,k)
    call sub_dihedral(coord(1,8),coord(2,8),coord(3,8),&                            
                      coord(1,7),coord(2,7),coord(3,7),&                            
                      coord(1,9),coord(2,9),coord(3,9),&                            
                      coord(1,11),coord(2,11),coord(3,11),&                         
                      qm_H_N_CA_CB(j,k))      
     write(204,*) phi,psi,qm_H_N_CA_CB(j,k)
    call sub_dihedral(coord(1,5),coord(2,5),coord(3,5),&                            
                      coord(1,7),coord(2,7),coord(3,7),&                            
                      coord(1,9),coord(2,9),coord(3,9),&                            
                      coord(1,10),coord(2,10),coord(3,10),&                         
                      qm_C_N_CA_H(j,k))      
     write(205,*) phi,psi,qm_C_N_CA_H(j,k)
    call sub_dihedral(coord(1,8),coord(2,8),coord(3,8),&                            
                      coord(1,7),coord(2,7),coord(3,7),&                            
                      coord(1,9),coord(2,9),coord(3,9),&                            
                      coord(1,10),coord(2,10),coord(3,10),&                         
                      qm_H_N_CA_H(j,k))         
     write(206,*) phi,psi,qm_H_N_CA_H(j,k)


    call sub_dihedral(coord(1,7),coord(2,7),coord(3,7),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      coord(1,17),coord(2,17),coord(3,17),&
                      qm_N_CA_C_N(j,k))
     write(207,*) phi,psi,qm_N_CA_C_N(j,k)
    call sub_dihedral(coord(1,10),coord(2,10),coord(3,10),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      coord(1,17),coord(2,17),coord(3,17),&
                      qm_HA_CA_C_N(j,k))
     write(208,*) phi,psi,qm_HA_CA_C_N(j,k)
    call sub_dihedral(coord(1,11),coord(2,11),coord(3,11),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      coord(1,17),coord(2,17),coord(3,17),&
                      qm_CB_CA_C_N(j,k))
     write(209,*) phi,psi,qm_CB_CA_C_N(j,k)
    call sub_dihedral(coord(1,7),coord(2,7),coord(3,7),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      coord(1,16),coord(2,16),coord(3,16),&
                      qm_N_CA_C_O(j,k))
     write(210,*) phi,psi,qm_N_CA_C_O(j,k)
    call sub_dihedral(coord(1,10),coord(2,10),coord(3,10),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      coord(1,16),coord(2,16),coord(3,16),&
                      qm_HA_CA_C_O(j,k))
     write(211,*) phi,psi,qm_HA_CA_C_O(j,k)
    call sub_dihedral(coord(1,11),coord(2,11),coord(3,11),&
                      coord(1,9),coord(2,9),coord(3,9),&
                      coord(1,15),coord(2,15),coord(3,15),&
                      coord(1,16),coord(2,16),coord(3,16),&
                      qm_CB_CA_C_O(j,k))
     write(212,*) phi,psi,qm_HA_CA_C_O(j,k)
  end do
end do

call system('mkdir qm_dih_d_value')
open(301,file='./qm_dih_d_value/qm_phi_02-01.dat')
open(302,file='./qm_dih_d_value/qm_phi_03-01.dat')
open(303,file='./qm_dih_d_value/qm_phi_04-01.dat')
open(304,file='./qm_dih_d_value/qm_phi_05-01.dat')
open(305,file='./qm_dih_d_value/qm_phi_06-01.dat')
open(306,file='./qm_dih_d_value/qm_psi_08-07.dat')
open(307,file='./qm_dih_d_value/qm_psi_09-07.dat')
open(308,file='./qm_dih_d_value/qm_psi_10-07.dat')
open(309,file='./qm_dih_d_value/qm_psi_11-07.dat')
open(310,file='./qm_dih_d_value/qm_psi_12-07.dat')
do j=1,phi_num
  phi=phi_initial+increm*(j-1)
  do k=1,psi_num
    psi=psi_initial+increm*(k-1)
    if(abs(qm_H_N_CA_C(j,k)-qm_C_N_CA_C(j,k)).ge.180.d0) then
      write(301,*) phi,psi,abs(qm_H_N_CA_C(j,k)-qm_C_N_CA_C(j,k))-360.d0
    else     
      write(301,*) phi,psi,abs(qm_H_N_CA_C(j,k)-qm_C_N_CA_C(j,k))
    end if

    if(abs(qm_C_N_CA_CB(j,k)-qm_C_N_CA_C(j,k)).ge.180.d0) then
      write(302,*) phi,psi,abs(qm_C_N_CA_CB(j,k)-qm_C_N_CA_C(j,k))-360.d0
    else
      write(302,*) phi,psi,abs(qm_C_N_CA_CB(j,k)-qm_C_N_CA_C(j,k))
    end if 
    
    if(abs(qm_H_N_CA_CB(j,k)-qm_C_N_CA_C(j,k)).ge.180.d0) then
      write(303,*) phi,psi,abs(qm_H_N_CA_CB(j,k)-qm_C_N_CA_C(j,k))-360.d0
    else
      write(303,*) phi,psi,abs(qm_H_N_CA_CB(j,k)-qm_C_N_CA_C(j,k))
    end if  
    
    if(abs(qm_C_N_CA_H(j,k)-qm_C_N_CA_C(j,k)).ge.180.d0) then
      write(304,*) phi,psi,abs(qm_C_N_CA_H(j,k)-qm_C_N_CA_C(j,k))-360.d0
    else
      write(304,*) phi,psi,abs(qm_C_N_CA_H(j,k)-qm_C_N_CA_C(j,k))
    end if
    
    if(abs(qm_H_N_CA_H(j,k)-qm_C_N_CA_C(j,k)).ge.180.d0) then
      write(305,*) phi,psi,abs(qm_H_N_CA_H(j,k)-qm_C_N_CA_C(j,k))-360.d0
    else
      write(305,*) phi,psi,abs(qm_H_N_CA_H(j,k)-qm_C_N_CA_C(j,k))
    end if  





    if(abs(qm_HA_CA_C_N(j,k)-qm_N_CA_C_N(j,k)).ge.180.d0) then
      write(306,*) phi,psi,abs(qm_HA_CA_C_N(j,k)-qm_N_CA_C_N(j,k))-360.d0
    else
      write(306,*) phi,psi,abs(qm_HA_CA_C_N(j,k)-qm_N_CA_C_N(j,k))
    end if 
    
    if(abs(qm_CB_CA_C_N(j,k)-qm_N_CA_C_N(j,k)).ge.180.d0) then
      write(307,*) phi,psi,abs(qm_CB_CA_C_N(j,k)-qm_N_CA_C_N(j,k))-360.d0
    else
      write(307,*) phi,psi,abs(qm_CB_CA_C_N(j,k)-qm_N_CA_C_N(j,k))
    end if
 
    if(abs(qm_N_CA_C_O(j,k)-qm_N_CA_C_N(j,k)).ge.180.d0) then 
      write(308,*) phi,psi,abs(qm_N_CA_C_O(j,k)-qm_N_CA_C_N(j,k))-360.d0
    else
      write(308,*) phi,psi,abs(qm_N_CA_C_O(j,k)-qm_N_CA_C_N(j,k))
    end if  
    
    if(abs(qm_HA_CA_C_O(j,k)-qm_N_CA_C_N(j,k)).ge.180.d0) then 
      write(309,*) phi,psi,abs(qm_HA_CA_C_O(j,k)-qm_N_CA_C_N(j,k))-360.d0
    else
      write(309,*) phi,psi,abs(qm_HA_CA_C_O(j,k)-qm_N_CA_C_N(j,k))
    end if  
  
    if(abs(qm_CB_CA_C_O(j,k)-qm_N_CA_C_N(j,k)).ge.180.d0) then
      write(310,*) phi,psi,abs(qm_CB_CA_C_O(j,k)-qm_N_CA_C_N(j,k))-360.d0
    else
      write(310,*) phi,psi,abs(qm_CB_CA_C_O(j,k)-qm_N_CA_C_N(j,k))
    end if
 
  end do 
end do

call system('mkdir mm_dih_cor_dat')
open(201,file='./mm_dih_cor_dat/mm_dih_cor01.dat')
open(202,file='./mm_dih_cor_dat/mm_dih_cor02.dat')
open(203,file='./mm_dih_cor_dat/mm_dih_cor03.dat')
open(204,file='./mm_dih_cor_dat/mm_dih_cor04.dat')
open(205,file='./mm_dih_cor_dat/mm_dih_cor05.dat')
open(206,file='./mm_dih_cor_dat/mm_dih_cor06.dat')
open(207,file='./mm_dih_cor_dat/mm_dih_cor07.dat')
open(208,file='./mm_dih_cor_dat/mm_dih_cor08.dat')
open(209,file='./mm_dih_cor_dat/mm_dih_cor09.dat')
open(210,file='./mm_dih_cor_dat/mm_dih_cor10.dat')
open(211,file='./mm_dih_cor_dat/mm_dih_cor11.dat')
open(212,file='./mm_dih_cor_dat/mm_dih_cor12.dat')
do j=1,phi_num
  phi=phi_initial+increm*(j-1) 
   write(char_phi,'(i4)') phi
  do k=1,psi_num  
    psi=psi_initial+increm*(k-1) 
     write(char_psi,'(i4)') psi
    rst_filename='ala_dipeptide_'//trim(adjustl(char_phi))//'_'//trim(adjustl(char_psi))//'.rst'
    call sub_get_rst_coord
    call sub_dihedral(atom_coord(1,5),atom_coord(2,5),atom_coord(3,5),& 
                      atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      mm_C_N_CA_C(j,k))
    write(201,*) phi,psi,mm_C_N_CA_C(j,k)
    call sub_dihedral(atom_coord(1,8),atom_coord(2,8),atom_coord(3,8),&
                      atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),& 
                      mm_H_N_CA_C(j,k))                   
    write(202,*) phi,psi,mm_H_N_CA_C(j,k)
    call sub_dihedral(atom_coord(1,5),atom_coord(2,5),atom_coord(3,5),&                            
                      atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&                            
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&                            
                      atom_coord(1,11),atom_coord(2,11),atom_coord(3,11),&                         
                      mm_C_N_CA_CB(j,k))      
    write(203,*) phi,psi,mm_C_N_CA_CB(j,k)
    call sub_dihedral(atom_coord(1,8),atom_coord(2,8),atom_coord(3,8),&                            
                      atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&                            
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&                            
                      atom_coord(1,11),atom_coord(2,11),atom_coord(3,11),&                         
                      mm_H_N_CA_CB(j,k))      
    write(204,*) phi,psi,mm_H_N_CA_CB(j,k)
    call sub_dihedral(atom_coord(1,5),atom_coord(2,5),atom_coord(3,5),&                            
                      atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&                            
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&                            
                      atom_coord(1,10),atom_coord(2,10),atom_coord(3,10),&                         
                      mm_C_N_CA_H(j,k))      
    write(205,*) phi,psi,mm_C_N_CA_H(j,k)
    call sub_dihedral(atom_coord(1,8),atom_coord(2,8),atom_coord(3,8),&                            
                      atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&                            
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&                            
                      atom_coord(1,10),atom_coord(2,10),atom_coord(3,10),&                         
                      mm_H_N_CA_H(j,k))         
    write(206,*) phi,psi,mm_H_N_CA_H(j,k)


    call sub_dihedral(atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      atom_coord(1,17),atom_coord(2,17),atom_coord(3,17),&
                      mm_N_CA_C_N(j,k))
    write(207,*) phi,psi,mm_N_CA_C_N(j,k)
    call sub_dihedral(atom_coord(1,10),atom_coord(2,10),atom_coord(3,10),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      atom_coord(1,17),atom_coord(2,17),atom_coord(3,17),&
                      mm_HA_CA_C_N(j,k))
    write(208,*) phi,psi,mm_HA_CA_C_N(j,k)
    call sub_dihedral(atom_coord(1,11),atom_coord(2,11),atom_coord(3,11),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      atom_coord(1,17),atom_coord(2,17),atom_coord(3,17),&
                      mm_CB_CA_C_N(j,k))
    write(209,*) phi,psi,mm_CB_CA_C_N(j,k)
    call sub_dihedral(atom_coord(1,7),atom_coord(2,7),atom_coord(3,7),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      atom_coord(1,16),atom_coord(2,16),atom_coord(3,16),&
                      mm_N_CA_C_O(j,k))
    write(210,*) phi,psi,mm_N_CA_C_O(j,k)
    call sub_dihedral(atom_coord(1,10),atom_coord(2,10),atom_coord(3,10),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      atom_coord(1,16),atom_coord(2,16),atom_coord(3,16),&
                      mm_HA_CA_C_O(j,k))
    write(211,*) phi,psi,mm_HA_CA_C_O(j,k)
    call sub_dihedral(atom_coord(1,11),atom_coord(2,11),atom_coord(3,11),&
                      atom_coord(1,9),atom_coord(2,9),atom_coord(3,9),&
                      atom_coord(1,15),atom_coord(2,15),atom_coord(3,15),&
                      atom_coord(1,16),atom_coord(2,16),atom_coord(3,16),&
                      mm_CB_CA_C_O(j,k))
    write(212,*) phi,psi,mm_HA_CA_C_O(j,k)
  end do
end do

call system('mkdir mm_dih_d_value')
open(301,file='./mm_dih_d_value/mm_phi_02-01.dat')
open(302,file='./mm_dih_d_value/mm_phi_03-01.dat')
open(303,file='./mm_dih_d_value/mm_phi_04-01.dat')
open(304,file='./mm_dih_d_value/mm_phi_05-01.dat')
open(305,file='./mm_dih_d_value/mm_phi_06-01.dat')
open(306,file='./mm_dih_d_value/mm_psi_08-07.dat')
open(307,file='./mm_dih_d_value/mm_psi_09-07.dat')
open(308,file='./mm_dih_d_value/mm_psi_10-07.dat')
open(309,file='./mm_dih_d_value/mm_psi_11-07.dat')
open(310,file='./mm_dih_d_value/mm_psi_12-07.dat')
do j=1,phi_num
  phi=phi_initial+increm*(j-1)
  do k=1,psi_num
    psi=psi_initial+increm*(k-1)
    if(abs(mm_H_N_CA_C(j,k)-mm_C_N_CA_C(j,k)).ge.180.d0) then
      write(301,*) phi,psi,abs(mm_H_N_CA_C(j,k)-mm_C_N_CA_C(j,k))-360.d0
    else     
      write(301,*) phi,psi,abs(mm_H_N_CA_C(j,k)-mm_C_N_CA_C(j,k))
    end if

    if(abs(mm_C_N_CA_CB(j,k)-mm_C_N_CA_C(j,k)).ge.180.d0) then
      write(302,*) phi,psi,abs(mm_C_N_CA_CB(j,k)-mm_C_N_CA_C(j,k))-360.d0
    else
      write(302,*) phi,psi,abs(mm_C_N_CA_CB(j,k)-mm_C_N_CA_C(j,k))
    end if 
    
    if(abs(mm_H_N_CA_CB(j,k)-mm_C_N_CA_C(j,k)).ge.180.d0) then
      write(303,*) phi,psi,abs(mm_H_N_CA_CB(j,k)-mm_C_N_CA_C(j,k))-360.d0
    else
      write(303,*) phi,psi,abs(mm_H_N_CA_CB(j,k)-mm_C_N_CA_C(j,k))
    end if  
    
    if(abs(mm_C_N_CA_H(j,k)-mm_C_N_CA_C(j,k)).ge.180.d0) then
      write(304,*) phi,psi,abs(mm_C_N_CA_H(j,k)-mm_C_N_CA_C(j,k))-360.d0
    else
       write(304,*) phi,psi,abs(mm_C_N_CA_H(j,k)-mm_C_N_CA_C(j,k))
    end if
    
    if(abs(mm_H_N_CA_H(j,k)-mm_C_N_CA_C(j,k)).ge.180.d0) then
      write(305,*) phi,psi,abs(mm_H_N_CA_H(j,k)-mm_C_N_CA_C(j,k))-360.d0
    else
      write(305,*) phi,psi,abs(mm_H_N_CA_H(j,k)-mm_C_N_CA_C(j,k))
    end if  





    if(abs(mm_HA_CA_C_N(j,k)-mm_N_CA_C_N(j,k)).ge.180.d0) then
      write(306,*) phi,psi,abs(mm_HA_CA_C_N(j,k)-mm_N_CA_C_N(j,k))-360.d0
    else
      write(306,*) phi,psi,abs(mm_HA_CA_C_N(j,k)-mm_N_CA_C_N(j,k))
    end if 
    
    if(abs(mm_CB_CA_C_N(j,k)-mm_N_CA_C_N(j,k)).ge.180.d0) then
      write(307,*) phi,psi,abs(mm_CB_CA_C_N(j,k)-mm_N_CA_C_N(j,k))-360.d0
    else
      write(307,*) phi,psi,abs(mm_CB_CA_C_N(j,k)-mm_N_CA_C_N(j,k))
    end if
 
    if(abs(mm_N_CA_C_O(j,k)-mm_N_CA_C_N(j,k)).ge.180.d0) then 
      write(308,*) phi,psi,abs(mm_N_CA_C_O(j,k)-mm_N_CA_C_N(j,k))-360.d0
    else
      write(308,*) phi,psi,abs(mm_N_CA_C_O(j,k)-mm_N_CA_C_N(j,k))
    end if  
    
    if(abs(mm_HA_CA_C_O(j,k)-mm_N_CA_C_N(j,k)).ge.180.d0) then 
      write(309,*) phi,psi,abs(mm_HA_CA_C_O(j,k)-mm_N_CA_C_N(j,k))-360.d0
    else
      write(309,*) phi,psi,abs(mm_HA_CA_C_O(j,k)-mm_N_CA_C_N(j,k))
    end if  
  
    if(abs(mm_CB_CA_C_O(j,k)-mm_N_CA_C_N(j,k)).ge.180.d0) then
      write(310,*) phi,psi,abs(mm_CB_CA_C_O(j,k)-mm_N_CA_C_N(j,k))-360.d0
    else
      write(310,*) phi,psi,abs(mm_CB_CA_C_O(j,k)-mm_N_CA_C_N(j,k))
    end if
 
  end do 
end do
end program
