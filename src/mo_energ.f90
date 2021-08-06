module qcxms_mo_energy
  use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
  implicit none

  contains 
    subroutine write_qmo(nat,nao,aoat,coeff,overlap,emo,focc,ihomo)

      integer  :: nat, nao, aoat(nao)
      integer  :: ia, n, j , k
      integer  :: ihomo
      integer  :: io_tmp, io_mspec

      real(wp) :: coeff(nao,nao), overlap(nao, nao)
      real(wp) :: qmo_sum, summa
      real(wp) :: qmo(nao, nat)
      real(wp) :: emo(nao),focc(nao)

      logical :: ex,nex

      ex = .false.
      nex = .false.

      inquire(file='tmp.mspec',exist=ex)
      if(ex) call system('rm tmp.mspec')
      open(file='tmp.mspec',newunit=io_tmp,status='new')

      inquire(file='qcxms.Mspec.tbxtb',exist=nex)
      if(nex) call system('rm qcxms.Mspec.tbxtb')
      open(file='qcxms.Mspec.tbxtb',newunit=io_mspec,status='new')

      do ia=1,nat
         do n=1, nao
            qmo_sum=0.0d0
            do j=1, nao
               if(aoat(j).eq.ia) then
                  do k=1, nao
                     qmo_sum=qmo_sum+coeff(j,n)*coeff(k,n)*overlap(j,k)
                  enddo
               endif
            enddo
            qmo_sum=qmo_sum+1.d-10
            qmo (n,ia)=qmo_sum
         enddo
      enddo


      do k=1,nao
         summa = 0.0d0
         do j=1,nat
            summa = summa+qmo(k,j)
         end do
         do j=1,nat
            qmo(k,j)=qmo(k,j)/summa
         end do
      end do
      
      write(io_mspec,*) nao, ihomo
      do  k = 1, nao
         write (io_mspec, *) 
         write (io_mspec,'(1(i3 ,1x, f10.3))') k, emo(k)*autoev 
         write (io_mspec,'(1(1x, f6.2))')     focc(k)
         write (io_mspec,'(10(1x, f6.2))')   (qmo(k,j)*100.0d0, j=1,nat)
      enddo
     
      write(io_tmp,*) nao, ihomo

      do  k = 1, nao
        write(io_tmp,*) emo(k)*autoev
        write(io_tmp,*) focc(k)
        do j=1,nat
           write (io_tmp, *) qmo(k,j)
        end do
      enddo


       close(io_tmp)
       close(io_mspec)

    end subroutine write_qmo

    subroutine readpopmo(nat,nao,ihomo,emo,focc,qmo)

       integer :: nat,ihomo,j,z,nao,io
       integer :: io_tmp
       !real*8 :: emo(2,5*nat),qmo(5*nat,nat)
       real(wp) :: emo(5*nat)
       real(wp) :: focc(5*nat)
       real(wp) :: qmo(5*nat,nat)

       OPEN(file='tmp.mspec',newunit=io_tmp,status='old')

       read(io_tmp,*) nao, ihomo

       do  j = 1, 5*nat
          read(io_tmp,*,IOSTAT=io) emo(j)
          read(io_tmp,*,IOSTAT=io) focc(j)
          if (io < 0) exit 
          do z = 1, nat
             read(io_tmp,*,IOSTAT=io) qmo(j,z)
             if (io < 0) exit  
          end do
          if (io < 0) exit  
       end do
       close(io_tmp, status='delete')
!       call system('rm tmp.mspec')
    end subroutine readpopmo

end module qcxms_mo_energy
