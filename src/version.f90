module get_version
contains 

   subroutine version(i)
     implicit none

     integer,intent(in) :: i
     character(len=:), allocatable  :: line

     line ='Aug 27 10:30:00 CEST 2021 '

     if (i.eq.0)  write(*,' (22x,''*'',18x,''V5.1.2'',18x,'' *'')'  )
     if (i.eq.1)  write(*,' (22x,''*        '',(a)''         *'')'     ) line
     if (i.eq.2)  write(*,' (22x,''--- QCxMS V5.1.2 '',(a)'' ---'')') line
     if (i.eq.33) write(33,'(22x,''--- QCxMS V5.1.2 '',(a)'' ---'')') line

   end subroutine version

end module get_version
