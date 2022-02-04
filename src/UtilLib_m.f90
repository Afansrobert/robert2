MODULE UtilLib_m
USE NumParameters_m
IMPLICIT NONE


CONTAINS

  !!@description: Change the case of a string (default lowercase)
  !!@param: name_string the string
  !!@param: lower If the variable is present and its value is F,
  !!              the string will be converted into a uppercase string, otherwise,
  !!              it will be convert into a lowercase string.
  SUBROUTINE string_uppercase_TO_lowercase(name_string,lower)
  IMPLICIT NONE

   character (len=*), intent(inout)          :: name_string
   logical,           intent(in),  optional  :: lower

   logical  :: lower_loc
   integer  :: i,ascii_char

   IF (present(lower)) THEN
     lower_loc = lower
   ELSE
     lower_loc = .TRUE.
   END IF

   !write(out_unitp,*) 'name_string: ',name_string
   IF (lower_loc) THEN ! uppercase => lowercase
     DO i=1,len_trim(name_string)
       ascii_char = iachar(name_string(i:i))
       IF (ascii_char >= 65 .AND. ascii_char <= 90)                 &
                           name_string(i:i) = achar(ascii_char+32)
     END DO
   ELSE ! lowercase => uppercase
     DO i=1,len_trim(name_string)
       ascii_char = iachar(name_string(i:i))
       IF (ascii_char >= 97 .AND. ascii_char <= 122)                 &
                            name_string(i:i) = achar(ascii_char-32)

     END DO
   END IF
   !write(out_unitp,*) 'name_string: ',name_string


  END SUBROUTINE string_uppercase_TO_lowercase

  PURE FUNCTION strdup(string)
  IMPLICIT NONE
   character (len=:), allocatable  :: strdup

   character (len=*), intent(in)   :: string

   allocate(character(len=len_trim(string)) :: strdup)
   strdup = trim(string)

  END FUNCTION strdup
  PURE FUNCTION int_TO_char(i)

    character (len=:), allocatable  :: int_TO_char

    integer, intent(in)             :: i

    character (len=:), allocatable  :: name_int
    integer :: clen

    ! first approximated size of name_int
    IF (i == 0) THEN
      clen = 1
    ELSE IF (i < 0) THEN
      clen = int(log10(abs(real(i,kind=Rk))))+2
    ELSE
      clen = int(log10(real(i,kind=Rk)))+1
    END IF

    ! allocate name_int
    allocate(character(len=clen) :: name_int)

    ! write i in name_int
    write(name_int,'(i0)') i

    ! transfert name_int in int_TO_char
    int_TO_char = strdup(name_int)

    ! deallocate name_int
    deallocate(name_int)

  END FUNCTION int_TO_char
  FUNCTION string_IS_empty(String)
    logical                          :: string_IS_empty
    character(len=*), intent(in)     :: String

    string_IS_empty = (len_trim(String) == 0)

  END FUNCTION string_IS_empty

  SUBROUTINE sub_Format_OF_Line(wformat,nb_line,max_col,cplx,       &
                                Rformat,name_info)
  IMPLICIT NONE

   character (len=*), optional  :: Rformat
   character (len=*), optional  :: name_info

   character (len=Name_longlen) :: wformat,name_info_loc
   integer                      :: nb_line,max_col
   logical                      :: cplx


   character (len=Name_longlen) :: iformat,cformat
   character (len=Name_longlen) :: NMatformat


   IF (present(name_info)) THEN
     name_info_loc = '(2x,"' // trim(adjustl(name_info)) // ': ",'
   ELSE
     name_info_loc = '('
   END IF

   IF (present(Rformat)) THEN
     IF (len_trim(Rformat) > 10) THEN
       write(out_unitp,*) ' ERROR in sub_Format_OF_Line'
       write(out_unitp,*) ' The format (len_trim) in "Rformat" is too long',len_trim(Rformat)
       write(out_unitp,*) ' Rformat: ',Rformat
       STOP
     END IF
       IF (cplx) THEN
         NMatformat = "'('," // trim(adjustl(Rformat)) //           &
                      ",' +i'," // trim(adjustl(Rformat)) // ",')'"
       ELSE
         NMatformat = trim(adjustl(Rformat))
       END IF
   ELSE
       IF (cplx) THEN
         NMatformat = trim(adjustl(CMatIO_format))
       ELSE
         NMatformat = trim(adjustl(RMatIO_format))
       END IF
   END IF

   IF (nb_line > 0) THEN

       write(iformat,*) int(log10(real(nb_line,kind=Rk)))+1
       write(cformat,*) max_col

       wformat = trim(adjustl(name_info_loc)) // '1x,' //           &
                 'i' // trim(adjustl(iformat)) // ',2x,' //         &
                 trim(adjustl(cformat)) // '(' //                   &
                 trim(adjustl(NMatformat)) // ',1x))'
   ELSE
       write(cformat,*) max_col

         wformat = trim(adjustl(name_info_loc)) //                  &
                   trim(adjustl(cformat)) // '(' //                 &
                   trim(adjustl(NMatformat)) // ',1x))'


   END IF

  !write(out_unitp,*) 'format?: ',trim(wformat)
  END SUBROUTINE sub_Format_OF_Line

  SUBROUTINE Write_RMat(f,nio,nbcol1,Rformat,name_info)
  IMPLICIT NONE

     real(kind=Rk),     intent(in)           :: f(:,:)
     integer,           intent(in)           :: nio,nbcol1

     character (len=*), intent(in), optional :: Rformat
     character (len=*), intent(in), optional :: name_info



     integer         :: nl,nc
     integer i,j,nb,nbblocs,nfin,nbcol
     character (len=:), allocatable :: wformat

     nl = size(f,dim=1)
     nc = size(f,dim=2)
     !write(out_unitp,*) 'nl,nc,nbcol',nl,nc,nbcol
     nbcol = nbcol1
     IF (nbcol > 10) nbcol=10
     nbblocs=int(nc/nbcol)
     IF (nbblocs*nbcol == nc) nbblocs=nbblocs-1


     IF (present(name_info)) THEN
       wformat = strdup( '(" ' // trim(name_info) // ': ",' )
     ELSE
       wformat = strdup( '(' )
     END IF
     IF (present(Rformat)) THEN
       wformat = strdup( wformat // 'i0,1x,10' // trim(Rformat) // ')' )
     ELSE
       wformat = strdup( wformat // 'i0,1x,10f18.9)' )
     END IF

     DO nb=0,nbblocs-1
       DO j=1,nl
         write(nio,wformat) j,(f(j,i+nb*nbcol),i=1,nbcol)
       END DO
       IF (nl > 1 ) write(nio,*)
     END DO
     DO j=1,nl
       nfin=nc-nbcol*nbblocs
       write(nio,wformat) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
     END DO
     flush(nio)

  END SUBROUTINE Write_RMat
  SUBROUTINE Write_RVec(l,nio,nbcol1,Rformat,name_info)
  IMPLICIT NONE

     real(kind=Rk),     intent(in)           :: l(:)
     integer,           intent(in)           :: nio,nbcol1

     character (len=*), intent(in), optional :: Rformat
     character (len=*), intent(in), optional :: name_info


    integer           :: n,i,nb,nbblocs,nfin,nbcol
    character (len=Name_longlen) :: wformat

    n = size(l)
    !write(out_unitp,*) 'n,nbcol',n,nbcol
    nbcol = nbcol1
    IF (nbcol > 10) nbcol=10
    nbblocs=int(n/nbcol)
    IF (nbblocs*nbcol == n) nbblocs=nbblocs-1


    IF (present(Rformat)) THEN
      IF (present(name_info)) THEN
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.,Rformat,name_info)
       ELSE
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.,Rformat=Rformat)
       END IF
     ELSE
       IF (present(name_info)) THEN
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.,name_info=name_info)
       ELSE
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.)
       END IF
     END IF


     DO nb=0,nbblocs-1
       write(nio,wformat) (l(i+nb*nbcol),i=1,nbcol)
     END DO
     nfin=n-nbcol*nbblocs
     write(nio,wformat) (l(i+nbcol*nbblocs),i=1,nfin)
     flush(nio)

  END SUBROUTINE Write_RVec

  SUBROUTINE Read_RMat(f,nio,nbcol,err)

  real(kind=Rk),    intent(inout) :: f(:,:)
  integer,          intent(in)    :: nio,nbcol
  integer,          intent(inout) :: err


  integer i,j,jj,nb,nbblocs,nfin,nl,nc

!$OMP    CRITICAL (Read_RMat_CRIT)

  nl = size(f,dim=1)
  nc = size(f,dim=2)
! write(out_unitp,*) 'nl,nc,nbcol',nl,nc,nbcol


  nbblocs=int(nc/nbcol)

  IF (nbblocs*nbcol == nc) nbblocs=nbblocs-1
  err = 0

  DO nb=0,nbblocs-1

      DO j=1,nl
        read(nio,*,IOSTAT=err) jj,(f(j,i+nb*nbcol),i=1,nbcol)
        IF (err /= 0) EXIT
      END DO

      IF (err /= 0) EXIT

      IF (nl > 1) read(nio,*,IOSTAT=err)
      IF (err /= 0) EXIT

  END DO

  IF (err == 0) THEN
    DO j=1,nl
      nfin=nc-nbcol*nbblocs
      read(nio,*,IOSTAT=err) jj,(f(j,i+nbcol*nbblocs),i=1,nfin)
      IF (err /= 0) EXIT
    END DO
  END IF

  IF (err /= 0) THEN
    CALL Write_RMat(f,out_unitp,nbcol)
    write(out_unitp,*) ' ERROR in Read_RMat'
    write(out_unitp,*) '  while reading a matrix'
    write(out_unitp,*) '  end of file or end of record'
    write(out_unitp,*) '  The matrix paramters: nl,nc,nbcol',nl,nc,nbcol
    write(out_unitp,*) ' Check your data !!'
  END IF

!$OMP    END CRITICAL (Read_RMat_CRIT)

  END SUBROUTINE Read_RMat

  SUBROUTINE Read_RVec(l,nio,nbcol,err)

  real(kind=Rk),    intent(inout) :: l(:)
  integer,          intent(in)    :: nio,nbcol
  integer,          intent(inout) :: err

  integer :: n,i,nb,nbblocs,nfin

!$OMP    CRITICAL (Read_RVec_CRIT)

  n = size(l,dim=1)
  nbblocs=int(n/nbcol)
  err = 0


  IF (nbblocs*nbcol == n) nbblocs=nbblocs-1

  DO nb=0,nbblocs-1
    read(nio,*,IOSTAT=err) (l(i+nb*nbcol),i=1,nbcol)
    IF (err /= 0) EXIT
  END DO

  nfin=n-nbcol*nbblocs
  read(nio,*,IOSTAT=err) (l(i+nbcol*nbblocs),i=1,nfin)

  IF (err /= 0) THEN
    write(out_unitp,*) ' ERROR in Read_RVec'
    write(out_unitp,*) '  while reading a vector'
    write(out_unitp,*) '  end of file or end of record'
    write(out_unitp,*) '  The vector paramters: n,nbcol',n,nbcol
    write(out_unitp,*) ' Check your data !!'
  END IF

!$OMP    END CRITICAL (Read_RVec_CRIT)

  END SUBROUTINE Read_RVec

END MODULE UtilLib_m
