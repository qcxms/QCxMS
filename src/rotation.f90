module qcxms_cid_rotation
   use diag3x3
   use qcxms_mdinit, only: ekinet
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_constants, only: kB, pi
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 

contains 

  subroutine euler_rotation(nuc, iat, xyz, velo, io_rotate)
    
    integer  :: nuc,iat(nuc)
    integer  :: io_rotate
    integer  :: i

    real(wp),intent(inout) :: xyz(3,nuc),velo(3,nuc)
    real(wp) :: aalpha,abeta,agamma
    real(wp) :: a,b,c
    real(wp) :: rotalpha(4,4),rotbeta(4,4),rotgamma(4,4),R(4,4)
    real(wp) :: dummat(4,4)
    real(wp) :: rpoint(4),vpoint(4)
    real(wp) :: nxyz(1:4,nuc),vnxyz(1:4,nuc)
    

    !Define Euler angles
    call random_number(a)
    call random_number(b)
    call random_number(c)
    aalpha = a*2*pi
    abeta  = b*2*pi
    agamma = c*pi

    !Define rotational matrices - y-axis-Convention
    rotalpha=reshape( (/1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,          &
        cos(aalpha),sin(aalpha),0.0d0,0.0d0,-1.0d0*sin(aalpha), &
        cos(aalpha),0.0d0,0.0d0,0.0d0,0.0d0,1.0d0/),(/4,4/) )
  
    rotbeta=reshape( (/COS(abeta),0.0d0,-1.0d0*SIN(abeta),0.0d0, &
         0.0d0,1.0d0,0.0d0,0.0d0,SIN(abeta),0.0d0,COS(abeta),    &
         0.0d0,0.0d0,0.0d0,0.0d0,1.0d0 /),(/4,4/) )
  
    rotgamma= reshape( (/ COS(agamma), SIN(agamma), 0.0d0, 0.0d0, &
         -1.0d0*SIN(agamma), COS(agamma), 0.0d0,0.0d0,            &
         0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0/),(/4,4/) )

    dummat = matmul(rotalpha,rotbeta)
    R      = matmul(dummat,rotgamma)
  
    !Perform rotation 
    do i=1,nuc
       rpoint(1) = xyz(1,i)
       rpoint(2) = xyz(2,i)
       rpoint(3) = xyz(3,i)
       rpoint(4) = 1.0d0
       vpoint(1) = velo(1,i)
       vpoint(2) = velo(2,i)
       vpoint(3) = velo(3,i)
       vpoint(4) = 1.0d0
    !nxyz = new xyz (rotated coord -nxyz DIM[4,atoms])
       nxyz(1:4,i) = matmul(R,rpoint)
       vnxyz(1:4,i) = matmul(R,vpoint)
    !re-initilize the coordinate point
       rpoint = 0.0d0
       vpoint = 0.0d0
    end do
    !print out new coords to rotate.xyz
    write(io_rotate,*) nuc
    write(io_rotate,*) ' '
    do i=1,nuc
       write(io_rotate,*) toSymbol(iat(i)),' ',nxyz(1,i)/aatoau  &
           ,' ',nxyz(2,i)/aatoau   &
           ,' ',nxyz(3,i)/aatoau
    end do
    close(io_rotate)
     
  !copy the new-coords to xyz
    xyz = 0.0d0
    do i =1,nuc
       xyz(1,i)  = nxyz(1,i)
       xyz(2,i)  = nxyz(2,i)
       xyz(3,i)  = nxyz(3,i)
       velo(1,i) = vnxyz(1,i) 
       velo(2,i) = vnxyz(2,i)
       velo(3,i) = vnxyz(3,i)
    end do

  end subroutine euler_rotation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rotation_velo(xyz, nuc, mass, velo, velo_rot, E_rot)

    integer  :: nuc, i

    real(wp) :: xyz(3,nuc),velo(3,nuc)
    real(wp) :: mass(nuc)
    real(wp) :: mat(3,3)
    real(wp) :: velo_rot(3,nuc),xyz_omega(3,3)
    real(wp) :: eigenvalue(3),eigenvector(3,3)
    real(wp) :: w_new(3)
    real(wp) :: E_Rot,E_Rot_Ang,ang_mom(3)
    real(wp) :: Ekin,Tinit


     ! calculate the main axes of inertia for rotational energy set-up    
     mat(:,:) = 0.0d0
     velo_rot(3,3) = 0.0d0
     xyz_omega = 0.0d0


     do i=1,nuc
        mat(1,1) = mat(1,1) + ( xyz(2,i)**2 + xyz(3,i)**2)*mass(i)
        mat(2,1) = mat(2,1) + (-xyz(1,i)    * xyz(2,i)   )*mass(i)
        mat(3,1) = mat(3,1) + (-xyz(1,i)    * xyz(3,i)   )*mass(i)
 
        mat(1,2) = mat(1,2) + (-xyz(2,i)    * xyz(1,i)   )*mass(i)
        mat(2,2) = mat(2,2) + ( xyz(1,i)**2 + xyz(3,i)**2)*mass(i)
        mat(3,2) = mat(3,2) + (-xyz(2,i)    * xyz(3,i)   )*mass(i)
                  
        mat(1,3) = mat(1,3) + (-xyz(3,i)    * xyz(1,i)   )*mass(i)
        mat(2,3) = mat(2,3) + (-xyz(3,i)    * xyz(2,i)   )*mass(i)
        mat(3,3) = mat(3,3) + ( xyz(1,i)**2 + xyz(2,i)**2)*mass(i)
     enddo 
 
     call eigvec3x3(mat,eigenvalue,eigenvector)
 
 !    do j=1,3
 !       write(*,*)  'Trägheitsmomente',j, eigenvalue(j)
 !       write(*,*)  'Haupt-Trägheitsachsen',j, eigenvector(j)
 !    enddo
 !    max_val = maxval(eigenvalue)
 !    min_val = minval(eigenvalue)
 
 !    verh = eigenvalue(1) / max_val
 !    write(*,*) 'Verhaltniss1', verh
 !    verh = eigenvalue(2) / max_val
 !    write(*,*) 'Verhaltniss2', verh
 !    verh = eigenvalue(3) / max_val
 !    write(*,*) 'Verhaltniss3', verh
 
 !    write(*,*)  'Hauptachsen     '
 !    write(*,*)  eigenvector
 !    write(*,*)  ''     

     ! Get Temperature and Ekin
     call ekinet(nuc,velo,mass,Ekin,Tinit)
 
  !   check = dot_product(eigenvector(1:3,3),eigenvector(1:3,2))
  !   write(*,*) 'Ortho', check
     w_new(1) = (sqrt((kB*Tinit) / eigenvalue(1)))  !eigenvealue * 2 ?!? wegen kB T/2 ! wo ist die 2?
     w_new(2) = (sqrt((kB*Tinit) / eigenvalue(2)))  
     w_new(3) = (sqrt((kB*Tinit) / eigenvalue(3)))
 
 !    write(*,*) 'wnew',w_new
 
     do i=1,3
       xyz_omega(1,i) = eigenvector(1,i) * w_new(i)  !* ((2*pi)/(800*fstoau)) ! * (0.5 * kB * 300) 
       xyz_omega(2,i) = eigenvector(2,i) * w_new(i)  !* ((2*pi)/(800*fstoau)) ! * (0.5 * kB * 300) 
       xyz_omega(3,i) = eigenvector(3,i) * w_new(i)  !* ((2*pi)/(800*fstoau)) ! * (0.5 * kB * 300) 
 
     enddo
 
     velo_rot = 0.0d0
     do i=1,nuc
        do j=1,3 !Along all Eigenvectors
          velo_rot(1,i) = velo_rot(1,i) + (xyz_omega(2,j)*xyz(3,i) - xyz_omega(3,j)*xyz(2,i)) 
          velo_rot(2,i) = velo_rot(2,i) + (xyz_omega(3,j)*xyz(1,i) - xyz_omega(1,j)*xyz(3,i)) 
          velo_rot(3,i) = velo_rot(3,i) + (xyz_omega(1,j)*xyz(2,i) - xyz_omega(2,j)*xyz(1,i)) 
        enddo
     enddo

    E_Rot= 0.0d0
    E_Rot_Ang= 0.0d0
    do i=1,3
       !ang_mom(i) = eigenvalue(i) * w_new(i)
       !E_Rot_Ang = E_Rot_ang + (ang_mom(i)**2 / (2*eigenvalue(i)))
       E_Rot = E_Rot + (0.5 * eigenvalue(i) *  (w_new(i)**2))
    enddo


  end subroutine rotation_velo

end module qcxms_cid_rotation
