!     subroutine prints matrix r,which is supposed
!     to have dimension n,m  when m is nonzero and
!     ((n+1)*n)/2 when m is zero

!     only used for printing the fragment distance matrix (analyse.f90)

subroutine print_matrix(r,n,m)
   use xtb_mctc_accuracy, only: wp

   real(wp) :: r(*)

   nkpb=6
   if(m)10,10,80
!
10 continue
   ibl=n/nkpb
   ir=n-ibl*nkpb
   j1=1
   k1s=1
   kd=0
   if(ibl.eq.0) go to 50
   j2=nkpb
   do 40 i=1,ibl
      write(*,1002)(j,j=j1,j2)
      k1=k1s
      k2=k1
      kk=0
      do 20 j=j1,j2
         write(*,1003)j,(r(k),k=k1,k2)
         kk=kk+1
         k1=k1+kd+kk
20    k2=k1+kk
      j1=j1+nkpb
      if(j1.gt.n) return
      j2=j2+nkpb
      k2=k1-1
      k1=k2+1
      k2=k1+(nkpb-1)
      k1s=k2+1
      kk=kd+nkpb
      do 30 j=j1,n
         write(*,1003)j,(r(k),k=k1,k2)
         kk=kk+1
         k1=k1+kk
30    k2=k2+kk
40 kd=kd+nkpb
50 if(ir.eq.0) go to 70
   k1=k1s
   j2=j1+ir-1
   kk=0
   k2=k1
   write(*,1002)(j,j=j1,j2)
   write(*,1003)
   do 60 j=j1,j2
      write(*,1003)j,(r(k),k=k1,k2)
      kk=kk+1
      k1=k1+kd+kk
60 k2=k1+kk
70 return
80 ibl=m/nkpb
   ir=m-ibl*nkpb
   i2=0
   k2=0
   if(ibl.eq.0) go to 100
   do 90 i=1,ibl
      i1=(i-1)*n*nkpb+1
      i2=i1+(nkpb-1)*n
      k1=k2+1
      k2=k1+(nkpb-1)
      write(*,1002)(k,k=k1,k2)
      do 90 j=1,n
         write(*,1003)j,(r(ij),ij=i1,i2,n)
         i1=i1+1
90 i2=i1+(nkpb-1)*n
100 if(ir.eq.0) go to 120
   i1=ibl*n*nkpb+1
   i2=i1+(ir-1)*n
   k1=k2+1
   k2=m
   write(*,1002)(k,k=k1,k2)
   write(*,1003)
   do 110 j=1,n
      write(*,1003)j,(r(ij),ij=i1,i2,n)
      i1=i1+1
      i2=i1+(ir-1)*n
110 continue
120 write(*,1003)
   return
1002 format(/,' ',4x,6(3x,i4,3x),/)
1003 format(' ',i4,6f10.5)
end


