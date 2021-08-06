!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! copy from machine generated data statements inside pars.f
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine copyc6
   use d3common
   implicit none

   integer nlines
   integer iat,jat,iadr,jadr,nn,kk
   include 'pars.fh'

! define r2r4 an rcov at this point

   call setr2r4(r2r4)
   call setrcov(rcov)

   c6ab=-1
   mxc=0
! process all entries in pars.f
   kk=1
   do nn=1,nlines
      iat=int(pars(kk+1))
      jat=int(pars(kk+2))
      call limit(iat,jat,iadr,jadr)
      mxc(iat)=max(mxc(iat),iadr)
      mxc(jat)=max(mxc(jat),jadr)

      c6ab(iat,jat,iadr,jadr,1)=pars(kk)
      c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
      c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

      c6ab(jat,iat,jadr,iadr,1)=pars(kk)
      c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
      c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
      kk=(nn*5)+1
   enddo

end subroutine copyc6


subroutine limit(iat,jat,iadr,jadr)
   implicit none
   integer iat,jat,iadr,jadr,i
   iadr=1
   jadr=1
   i=100
10 if(iat.gt.100) then
      iat=iat-100
      iadr=iadr+1
      goto 10
   endif

   i=100
20 if(jat.gt.100) then
      jat=jat-100
      jadr=jadr+1
      goto 20
   endif

end subroutine limit

subroutine setrcov(rcov)
   use xtb_mctc_accuracy, only: wp
   real(wp) :: rcov(94)

   rcov( 1 )= 0.80628308_wp
   rcov( 2 )= 1.15903197_wp
   rcov( 3 )= 3.02356173_wp
   rcov( 4 )= 2.36845659_wp
   rcov( 5 )= 1.94011865_wp
   rcov( 6 )= 1.88972601_wp
   rcov( 7 )= 1.78894056_wp
   rcov( 8 )= 1.58736983_wp
   rcov( 9 )= 1.61256616_wp
   rcov( 10 )= 1.68815527_wp
   rcov( 11 )= 3.52748848_wp
   rcov( 12 )= 3.14954334_wp
   rcov( 13 )= 2.84718717_wp
   rcov( 14 )= 2.62041997_wp
   rcov( 15 )= 2.77159820_wp
   rcov( 16 )= 2.57002732_wp
   rcov( 17 )= 2.49443835_wp
   rcov( 18 )= 2.41884923_wp
   rcov( 19 )= 4.43455700_wp
   rcov( 20 )= 3.88023730_wp
   rcov( 21 )= 3.35111422_wp
   rcov( 22 )= 3.07395437_wp
   rcov( 23 )= 3.04875805_wp
   rcov( 24 )= 2.77159820_wp
   rcov( 25 )= 2.69600923_wp
   rcov( 26 )= 2.62041997_wp
   rcov( 27 )= 2.51963467_wp
   rcov( 28 )= 2.49443835_wp
   rcov( 29 )= 2.54483100_wp
   rcov( 30 )= 2.74640188_wp
   rcov( 31 )= 2.82199085_wp
   rcov( 32 )= 2.74640188_wp
   rcov( 33 )= 2.89757982_wp
   rcov( 34 )= 2.77159820_wp
   rcov( 35 )= 2.87238349_wp
   rcov( 36 )= 2.94797246_wp
   rcov( 37 )= 4.76210950_wp
   rcov( 38 )= 4.20778980_wp
   rcov( 39 )= 3.70386304_wp
   rcov( 40 )= 3.50229216_wp
   rcov( 41 )= 3.32591790_wp
   rcov( 42 )= 3.12434702_wp
   rcov( 43 )= 2.89757982_wp
   rcov( 44 )= 2.84718717_wp
   rcov( 45 )= 2.84718717_wp
   rcov( 46 )= 2.72120556_wp
   rcov( 47 )= 2.89757982_wp
   rcov( 48 )= 3.09915070_wp
   rcov( 49 )= 3.22513231_wp
   rcov( 50 )= 3.17473967_wp
   rcov( 51 )= 3.17473967_wp
   rcov( 52 )= 3.09915070_wp
   rcov( 53 )= 3.32591790_wp
   rcov( 54 )= 3.30072128_wp
   rcov( 55 )= 5.26603625_wp
   rcov( 56 )= 4.43455700_wp
   rcov( 57 )= 4.08180818_wp
   rcov( 58 )= 3.70386304_wp
   rcov( 59 )= 3.98102289_wp
   rcov( 60 )= 3.95582657_wp
   rcov( 61 )= 3.93062995_wp
   rcov( 62 )= 3.90543362_wp
   rcov( 63 )= 3.80464833_wp
   rcov( 64 )= 3.82984466_wp
   rcov( 65 )= 3.80464833_wp
   rcov( 66 )= 3.77945201_wp
   rcov( 67 )= 3.75425569_wp
   rcov( 68 )= 3.75425569_wp
   rcov( 69 )= 3.72905937_wp
   rcov( 70 )= 3.85504098_wp
   rcov( 71 )= 3.67866672_wp
   rcov( 72 )= 3.45189952_wp
   rcov( 73 )= 3.30072128_wp
   rcov( 74 )= 3.09915070_wp
   rcov( 75 )= 2.97316878_wp
   rcov( 76 )= 2.92277614_wp
   rcov( 77 )= 2.79679452_wp
   rcov( 78 )= 2.82199085_wp
   rcov( 79 )= 2.84718717_wp
   rcov( 80 )= 3.32591790_wp
   rcov( 81 )= 3.27552496_wp
   rcov( 82 )= 3.27552496_wp
   rcov( 83 )= 3.42670319_wp
   rcov( 84 )= 3.30072128_wp
   rcov( 85 )= 3.47709584_wp
   rcov( 86 )= 3.57788113_wp
   rcov( 87 )= 5.06446567_wp
   rcov( 88 )= 4.56053862_wp
   rcov( 89 )= 4.20778980_wp
   rcov( 90 )= 3.98102289_wp
   rcov( 91 )= 3.82984466_wp
   rcov( 92 )= 3.85504098_wp
   rcov( 93 )= 3.88023730_wp
   rcov( 94 )= 3.90543362_wp

   return

end subroutine setrcov

subroutine setr2r4(r2r4)
   use xtb_mctc_accuracy, only: wp
   real(wp) :: r2r4(94)

   r2r4( 1 )= 2.00734898_wp
   r2r4( 2 )= 1.56637132_wp
   r2r4( 3 )= 5.01986934_wp
   r2r4( 4 )= 3.85379032_wp
   r2r4( 5 )= 3.64446594_wp
   r2r4( 6 )= 3.10492822_wp
   r2r4( 7 )= 2.71175247_wp
   r2r4( 8 )= 2.59361680_wp
   r2r4( 9 )= 2.38825250_wp
   r2r4( 10 )= 2.21522516_wp
   r2r4( 11 )= 6.58585536_wp
   r2r4( 12 )= 5.46295967_wp
   r2r4( 13 )= 5.65216669_wp
   r2r4( 14 )= 4.88284902_wp
   r2r4( 15 )= 4.29727576_wp
   r2r4( 16 )= 4.04108902_wp
   r2r4( 17 )= 3.72932356_wp
   r2r4( 18 )= 3.44677275_wp
   r2r4( 19 )= 7.97762753_wp
   r2r4( 20 )= 7.07623947_wp
   r2r4( 21 )= 6.60844053_wp
   r2r4( 22 )= 6.28791364_wp
   r2r4( 23 )= 6.07728703_wp
   r2r4( 24 )= 5.54643096_wp
   r2r4( 25 )= 5.80491167_wp
   r2r4( 26 )= 5.58415602_wp
   r2r4( 27 )= 5.41374528_wp
   r2r4( 28 )= 5.28497229_wp
   r2r4( 29 )= 5.22592821_wp
   r2r4( 30 )= 5.09817141_wp
   r2r4( 31 )= 6.12149689_wp
   r2r4( 32 )= 5.54083734_wp
   r2r4( 33 )= 5.06696878_wp
   r2r4( 34 )= 4.87005108_wp
   r2r4( 35 )= 4.59089647_wp
   r2r4( 36 )= 4.31176304_wp
   r2r4( 37 )= 9.55461698_wp
   r2r4( 38 )= 8.67396077_wp
   r2r4( 39 )= 7.97210197_wp
   r2r4( 40 )= 7.43439917_wp
   r2r4( 41 )= 6.58711862_wp
   r2r4( 42 )= 6.19536215_wp
   r2r4( 43 )= 6.01517290_wp
   r2r4( 44 )= 5.81623410_wp
   r2r4( 45 )= 5.65710424_wp
   r2r4( 46 )= 5.52640661_wp
   r2r4( 47 )= 5.44263305_wp
   r2r4( 48 )= 5.58285373_wp
   r2r4( 49 )= 7.02081898_wp
   r2r4( 50 )= 6.46815523_wp
   r2r4( 51 )= 5.98089120_wp
   r2r4( 52 )= 5.81686657_wp
   r2r4( 53 )= 5.53321815_wp
   r2r4( 54 )= 5.25477007_wp
   r2r4( 55 )= 11.02204549_wp
   r2r4( 56 )= 10.15679528_wp
   r2r4( 57 )= 9.35167836_wp
   r2r4( 58 )= 9.06926079_wp
   r2r4( 59 )= 8.97241155_wp
   r2r4( 60 )= 8.90092807_wp
   r2r4( 61 )= 8.85984840_wp
   r2r4( 62 )= 8.81736827_wp
   r2r4( 63 )= 8.79317710_wp
   r2r4( 64 )= 7.89969626_wp
   r2r4( 65 )= 8.80588454_wp
   r2r4( 66 )= 8.42439218_wp
   r2r4( 67 )= 8.54289262_wp
   r2r4( 68 )= 8.47583370_wp
   r2r4( 69 )= 8.45090888_wp
   r2r4( 70 )= 8.47339339_wp
   r2r4( 71 )= 7.83525634_wp
   r2r4( 72 )= 8.20702843_wp
   r2r4( 73 )= 7.70559063_wp
   r2r4( 74 )= 7.32755997_wp
   r2r4( 75 )= 7.03887381_wp
   r2r4( 76 )= 6.68978720_wp
   r2r4( 77 )= 6.05450052_wp
   r2r4( 78 )= 5.88752022_wp
   r2r4( 79 )= 5.70661499_wp
   r2r4( 80 )= 5.78450695_wp
   r2r4( 81 )= 7.79780729_wp
   r2r4( 82 )= 7.26443867_wp
   r2r4( 83 )= 6.78151984_wp
   r2r4( 84 )= 6.67883169_wp
   r2r4( 85 )= 6.39024318_wp
   r2r4( 86 )= 6.09527958_wp
   r2r4( 87 )= 11.79156076_wp
   r2r4( 88 )= 11.10997644_wp
   r2r4( 89 )= 9.51377795_wp
   r2r4( 90 )= 8.67197068_wp
   r2r4( 91 )= 8.77140725_wp
   r2r4( 92 )= 8.65402716_wp
   r2r4( 93 )= 8.53923501_wp
   r2r4( 94 )= 8.85024712_wp

   return

end subroutine setr2r4

