MODULE Basis_m
  USE NumParameters_m

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Basis_t,Read_Basis,Write_Basis,Basis_IS_allocated

  TYPE :: Basis_t
    integer                         :: nb         = 0
    integer                         :: nq         = 0

    character (len=:),  allocatable :: Basis_name

    real(kind=Rk),   allocatable :: x(:)
    real(kind=Rk),   allocatable :: w(:)
    real(kind=Rk),   allocatable :: d0gb(:,:)      ! basis functions d0gb(nq,nb)
    real(kind=Rk),   allocatable :: d1gb(:,:,:)    ! basis functions d2gb(nq,nb,1)
    real(kind=Rk),   allocatable :: d2gb(:,:,:,:)  ! basis functions d2gb(nq,nb,1,1)

  END TYPE Basis_t

CONTAINS
  FUNCTION Basis_IS_allocated(Basis) RESULT(alloc)

    TYPE(Basis_t),   intent(in)  :: Basis
    logical                      :: alloc

    alloc =             allocated(Basis%x)
    alloc = alloc .AND. allocated(Basis%w)
    alloc = alloc .AND. allocated(Basis%d0gb)
    alloc = alloc .AND. allocated(Basis%d1gb)
    alloc = alloc .AND. allocated(Basis%d2gb)

  END FUNCTION Basis_IS_allocated
  RECURSIVE SUBROUTINE Write_Basis(Basis)
  USE UtilLib_m

    TYPE(Basis_t),       intent(in)  :: Basis

    write(out_unitp,*) '-------------------------------------------------'
    write(out_unitp,*) 'Write_Basis'
    write(out_unitp,*) 'nb,nq',Basis%nb,Basis%nq
    IF (Basis_IS_allocated(Basis)) THEN
      CALL Write_RVec(Basis%x,out_unitp,5,name_info='x')
      write(out_unitp,*)
      CALL Write_RVec(Basis%w,out_unitp,5,name_info='w')
      write(out_unitp,*)
      CALL Write_RMat(Basis%d0gb,out_unitp,5,name_info='d0gb')
      write(out_unitp,*)
      CALL Write_RMat(Basis%d1gb(:,:,1),out_unitp,5,name_info='d1gb')
      write(out_unitp,*)
      CALL Write_RMat(Basis%d2gb(:,:,1,1),out_unitp,5,name_info='d2gb')
    ELSE
      write(out_unitp,*) ' Basis tables (x, w, dngb) are not allocated.'
    END IF

    write(out_unitp,*) '-------------------------------------------------'

  END SUBROUTINE Write_Basis
  RECURSIVE SUBROUTINE Read_Basis(Basis,nio)
    USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis
    integer,             intent(in)     :: nio



    integer                         :: err_io

    integer                         :: nb,nq,i,j
    character (len=Name_len)        :: name
    real(kind=Rk)                   :: A,B,scaleQ,Q0


    NAMELIST /basis_nD/ name,nb,nq,A,B,scaleQ,Q0


    nb        = 0
    nq        = 0
    A         = ZERO
    B         = ZERO
    Q0        = ZERO
    scaleQ    = ONE
    name      = '0'

    read(nio,nml=basis_nD,IOSTAT=err_io)
    write(out_unitp,nml=basis_nD)
    IF (err_io < 0) THEN
      write(out_unitp,basis_nD)
      write(out_unitp,*) ' ERROR in Read_Basis'
      write(out_unitp,*) '  while reading the namelist "basis_nD"'
      write(out_unitp,*) ' end of file or end of record'
      write(out_unitp,*) ' Probably, you forget a basis set ...'
      write(out_unitp,*) ' Check your data !!'
      STOP ' ERROR in Read_Basis: problems with the namelist.'
    END IF
    IF (err_io > 0) THEN
      write(out_unitp,basis_nD)
      write(out_unitp,*) ' ERROR in Read_Basis'
      write(out_unitp,*) '  while reading the namelist "basis_nD"'
      write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
      write(out_unitp,*) ' Check your data !!'
      STOP ' ERROR in Read_Basis: problems with the namelist.'
    END IF

    Basis%nb        = nb
    Basis%nq        = nq
    Basis%Basis_name      = trim(adjustl(name))
    CALL string_uppercase_TO_lowercase(Basis%Basis_name)

    SELECT CASE (Basis%Basis_name)
    CASE ('boxab')
     CALL Construct_Basis_Sin(Basis)
      Q0      = A
      scaleQ  = pi/(B-A)
      CASE ('herm','Ho')
     CALL Construct_Basis_Ho(Basis)
    CASE default
      STOP 'ERROR in Read_Basis: no default basis.'
    END SELECT

    CALL Scale_Basis(Basis,Q0,scaleQ)
    CALL CheckOrtho_Basis(Basis,nderiv=2)

    CALL Write_Basis(Basis)

  END SUBROUTINE Read_Basis
  SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
  USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis


    real(kind=Rk)          :: dx
    integer                :: ib,iq,nb,nq

    nb = Basis%nb
    nq = Basis%nq
    dx = pi/nq

    ! grid and weight
    Basis%x = [(dx*(iq-HALF),iq=1,nq)]
    Basis%w = [(dx,iq=1,nq)]

    allocate(Basis%d0gb(nq,nb))
    allocate(Basis%d1gb(nq,nb,1))
    allocate(Basis%d2gb(nq,nb,1,1))

    DO ib=1,nb
      Basis%d0gb(:,ib)     =          sin(Basis%x(:)*ib) / sqrt(pi*HALF)
      Basis%d1gb(:,ib,1)   =  ib    * cos(Basis%x(:)*ib) / sqrt(pi*HALF)
      Basis%d2gb(:,ib,1,1) = -ib**2 * Basis%d0gb(:,ib)
    END DO

    IF (nb == nq) THEN
      Basis%d0gb(:,nb)      = Basis%d0gb(:,nb)      / sqrt(TWO)
      Basis%d1gb(:,nb,:)    = Basis%d1gb(:,nb,:)    / sqrt(TWO)
      Basis%d2gb(:,nb,:,:)  = Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
    END IF

  END SUBROUTINE Construct_Basis_Sin


!!!!! la modification de Robert!!!!!!!!!!
 SUBROUTINE Construct_Basis_Ho(Basis) ! HO :
  USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis

    integer                :: i,j,nb,nq

    nb = Basis%nb
    nq = Basis%nq
    allocate(Basis%x(nq))
    allocate(Basis%w(nq))

    call hercom(nq, Basis%x(:), Basis%w(:))

    allocate(Basis%d0gb(nq,nb))
    allocate(Basis%d1gb(nq,nb,1))
    allocate(Basis%d2gb(nq,nb,1,1))

    DO i = 1, nq
          DO j = 1, nb
           CALL Construct_Basis_poly_Hermite_exp(Basis%x(i),Basis%d0gb(i,j), Basis%d1gb(i,j,1),Basis%d2gb(i,j,1,1), j-1,.TRUE.)
          END DO
      END DO


  END SUBROUTINE Construct_Basis_Ho

FUNCTION poly_Hermite(x,l)
    Implicit none
    real(kind = Rk):: poly_Hermite
    real(kind = Rk):: pl0,pl1,pl2,norme,x
    integer        :: i,l

      IF ( l .LT. 0 ) THEN
         Write(out_unitp,*) 'Bad arguments in poly_hermite :'
         Write(out_unitp,*) ' l < 0 : ',l
         STOP
        end if

       norme  =  sqrt(PI)

       IF (l .EQ. 0) THEN
         poly_Hermite = ONE/sqrt(norme)
       ELSE IF (l .EQ. 1) THEN
           norme = norme * TWO
         poly_Hermite = TWO * x/sqrt(norme)
       ELSE

         pl2 = ONE
         pl1 = TWO * x
         norme = norme * TWO

         DO i=2,l
           norme = norme * TWO * i
           pl0 = TWO*( x*pl1 - (i-1)*pl2 )
           pl2 = pl1
           pl1 = pl0
         END DO
         poly_Hermite = pl0/sqrt(norme)
       END IF

     END FUNCTION poly_Hermite

 FUNCTION gamma_perso(n)
    Implicit none

    real(kind = Rk)  :: gamma_perso
    real(kind = Rk)  :: a
    integer          :: i,n
       IF (n .LE. 0) THEN
           write(out_unitp,*) 'ERROR: gamma( n<=0)',n
           STOP
         END IF
         a = ONE
         DO i = 1,n-1
           a = a * dble (i)
      END DO
      gamma_perso = a
  END FUNCTION gamma_perso

   SUBROUTINE herrec ( p2, dp2, p1, x, nq )
      Implicit none
      integer       ::i
      integer       :: nq
      real(kind = Rk):: dp0
      real(kind = Rk):: dp1
      real(kind = RK):: dp2

      real(kind = Rk):: p0
      real(kind = Rk):: p1
      real(kind = Rk):: p2
      real(kind = Rk):: x

      p1  = ONE
      dp1 = ZERO

      p2  = x
      dp2 = ONE

      DO i = 2, nq

        p0  = p1
        dp0 = dp1

        p1  = p2
        dp1 = dp2

        p2  = x * p1 - HALF * ( dble ( i ) - ONE ) * p0
        dp2 = x * dp1 + p1 - HALF * ( dble ( i ) - ONE ) * dp0

    END DO
   END SUBROUTINE herrec

  SUBROUTINE herroot ( x, nq, dp2, p1 )
    Implicit none
      integer          :: i
      integer          :: nq
      real(kind = Rk),parameter  :: eps = TEN**(-TWELVE) ! 1.0d-12
      real(kind = Rk)  :: d
      real(kind = Rk)  :: dp2
      real(kind = Rk)  :: p1
      real(kind = Rk)  :: p2
      real(kind = Rk)  :: x

      DO i = 1, 10

        CALL herrec ( p2, dp2, p1, x, nq )

        d = p2 / dp2
        x = x - d

        IF ( ABS ( d ) .LE. eps * ( ABS ( x ) + ONE ) ) THEN
          RETURN
        END IF

      END DO
   END SUBROUTINE herroot

SUBROUTINE hercom (nq,xp,w)
   Implicit none

      integer        :: i,nq
      real(kind = Rk):: cc,dp2,p1,s,temp,x
      real(kind = Rk):: w(nq),xp(nq)



      cc = 1.7724538509_Rk * gamma_perso(nq ) / ( TWO**( nq-1) )

      s = ( TWO * dble (real(nq,Kind=Rk) ) + ONE )**( SIXTH )

      DO i = 1, ( nq + 1 ) / 2

       IF ( i .EQ. 1 ) THEN

          x = s**3 - 1.85575_Rk / s

        ELSE IF ( i .EQ. 2 ) THEN

          x = x - 1.14_Rk * ( ( dble ( nq ) )**0.426_Rk ) / x

        ELSE IF ( i .EQ. 3 ) THEN

          x = 1.86_Rk * x - 0.86_Rk * xp(1)

        ELSE IF ( i .EQ. 4 ) THEN

          x = 1.91_Rk * x - 0.91_Rk * xp(2)

        ELSE

          x = TWO * x - xp(i-2)

        END IF

        CALL herroot ( x,  nq, dp2, p1 )

        xp(i) = x
        W(i) = cc / dp2 / p1

        xp( nq-i+1) = - x
        w( nq-i+1) = w(i)

      END DO

      DO i = 1,  nq/2
        temp = xp(i)
        xp(i) = xp( nq+1-i)
        xp( nq+1-i) = temp
      END DO

      DO i = 1, nq
     w(i) = w(i)*exp(xp(i)*xp(i))
      END DO
 END SUBROUTINE hercom

 SUBROUTINE Construct_Basis_poly_Hermite_exp(x,d0gb,d1gb,d2gb,l,deriv)


      logical deriv

      integer        :: l

      real(kind = RK):: pexp,x,d0gb,d1gb,d2gb

       IF (deriv) THEN

          d0gb = poly_Hermite( x,l)
         IF (l .EQ. 0) THEN
          d1gb     = ZERO
          d2gb     = ZERO
         ELSE IF (l .EQ. 1) THEN
          d1gb = sqrt(TWO)*poly_Hermite( x,0)
          d2gb = ZERO

         ELSE IF (l .EQ. 2) THEN
          d1gb = sqrt(TWO*l) * poly_Hermite( x,l-1)
          d2gb = TWO*( x*d1gb-d0gb *l)
         ELSE
          d1gb = sqrt(TWO*l) * poly_Hermite( x,l-1)
          d2gb = TWO*( x* d1gb-d0gb*l)
         END IF

          pexp = exp(- HALF* x* x)
          d2gb = (d2gb-TWO*x*d1gb+( x* x-ONE)*d0gb)*pexp
          d1gb = (d1gb- x*d0gb)*pexp
          d0gb = d0gb*pexp

       ELSE
         d0gb = poly_Hermite(x ,l)*exp(-HALF* x* x)
         d1gb = ZERO
         d2gb = ZERO
       END IF

  END SUBROUTINE Construct_Basis_poly_Hermite_exp

!!!!!Fin de la modification de Robert!!!!!!!!!!

  SUBROUTINE CheckOrtho_Basis(Basis,nderiv)
  USE UtilLib_m

    TYPE(Basis_t),           intent(in)     :: Basis
    integer,                 intent(in)     :: nderiv

    integer                      :: ib
    real(kind=Rk), ALLOCATABLE   :: S(:,:)
    real(kind=Rk), ALLOCATABLE   :: d0bgw(:,:)
    real(kind=Rk)                :: Sii,Sij


    IF (Basis_IS_allocated(Basis)) THEN
      d0bgw = transpose(Basis%d0gb)
      DO ib=1,Basis%nb
        d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
      END DO

      S = matmul(d0bgw,Basis%d0gb)
      IF (nderiv > -1) CALL Write_RMat(S,out_unitp,5,name_info='S')
      Sii = ZERO
      Sij = ZERO
      DO ib=1,Basis%nb
        IF (abs(S(ib,ib)-ONE) > Sii) Sii = abs(S(ib,ib)-ONE)
        S(ib,ib) = ZERO
      END DO
      Sij = maxval(S)
      write(out_unitp,*) 'Sii,Sij',Sii,Sij

      IF (nderiv > 0) THEN
        write(out_unitp,*)
        S = matmul(d0bgw,Basis%d1gb(:,:,1))
        !CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>',Rformat='e13.4')
        CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
      END IF

      IF (nderiv > 1) THEN
        write(out_unitp,*)
        S = matmul(d0bgw,Basis%d2gb(:,:,1,1))
        !CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d2b>',Rformat='e13.4')
        CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
      END IF

    ELSE
      write(out_unitp,*) ' WARNNING in CheckOrtho_Basis'
      write(out_unitp,*) ' the basis is not allocated.'
    END IF

  END SUBROUTINE CheckOrtho_Basis
  SUBROUTINE Scale_Basis(Basis,x0,sx)
  USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis
    real(kind=Rk),       intent(in)     :: x0,sx

    IF (abs(sx) > ONETENTH**6 .AND. Basis_IS_allocated(Basis)) THEN

      Basis%x(:) = x0 + Basis%x(:) / sx
      Basis%w(:) =      Basis%w(:) / sx

      Basis%d0gb(:,:)     = Basis%d0gb(:,:)     * sqrt(sx)
      Basis%d1gb(:,:,:)   = Basis%d1gb(:,:,:)   * sqrt(sx)*sx
      Basis%d2gb(:,:,:,:) = Basis%d2gb(:,:,:,:) * sqrt(sx)*sx*sx
    ELSE
      write(out_unitp,*) ' ERROR in Scale_Basis'
      write(out_unitp,*) ' sx is too small  or ...'
      write(out_unitp,*) ' the basis is not allocated.'
      STOP 'ERROR in Scale_Basis'
    END IF

  END SUBROUTINE Scale_Basis

END MODULE Basis_m
