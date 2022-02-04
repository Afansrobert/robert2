  MODULE NumParameters_m
      USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64
      IMPLICIT NONE

      PUBLIC
      PRIVATE :: INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64

      integer, parameter :: Rk            = real64 ! 8
      integer, parameter :: Ik            = int32  ! 4
      integer, parameter :: ILk           = int64  ! 8
      integer, parameter :: Name_len      = 20
      integer, parameter :: Name_longlen  = 50
      integer, parameter :: Line_len      = 255


      real (kind=Rk), parameter :: ZERO    = 0._Rk
      real (kind=Rk), parameter :: ONE     = 1._Rk
      real (kind=Rk), parameter :: TWO     = 2._Rk
      real (kind=Rk), parameter :: THREE   = 3._Rk
      real (kind=Rk), parameter :: FOUR    = 4._Rk
      real (kind=Rk), parameter :: FIVE    = 5._Rk
      real (kind=Rk), parameter :: SIX     = 6._Rk
      real (kind=Rk), parameter :: SEVEN   = 7._Rk
      real (kind=Rk), parameter :: EIGHT   = 8._Rk
      real (kind=Rk), parameter :: NINE    = 9._Rk
      real (kind=Rk), parameter :: TEN     = 10._Rk
      real (kind=Rk), parameter :: ELEVEN  = 11._Rk
      real (kind=Rk), parameter :: TWELVE  = 12._Rk
      real (kind=Rk), parameter :: HUNDRED = 100._Rk

      real (kind=Rk), parameter :: HALF      = ONE/TWO
      real (kind=Rk), parameter :: THIRD     = ONE/THREE
      real (kind=Rk), parameter :: FOURTH    = ONE/FOUR
      real (kind=Rk), parameter :: QUARTER   = ONE/FOUR
      real (kind=Rk), parameter :: FIFTH     = ONE/FIVE
      real (kind=Rk), parameter :: SIXTH     = ONE/SIX
      real (kind=Rk), parameter :: ONETENTH  = ONE/TEN
      real (kind=Rk), parameter :: TWOTENTHS = TWO/TEN

      real (kind=Rk), parameter ::                                              &
                    PI = 3.14159265358979323846264338327950288419716939937511_Rk

      complex (kind=Rk), parameter :: EYE      = (0._Rk,1._Rk)
      complex (kind=Rk), parameter :: CZERO    = (0._Rk,0._Rk)
      complex (kind=Rk), parameter :: CONE     = (1._Rk,0._Rk)
      complex (kind=Rk), parameter :: CHALF    = (0.5_Rk,0._Rk)

      integer :: print_level = 0        ! 0 minimal, 1 default, 2 large, -1 nothing

      character (len=Name_longlen) :: EneIO_format = "f20.5"
      character (len=Name_longlen) :: RMatIO_format = "f18.10"
      character (len=Name_longlen) :: CMatIO_format = "'(',f15.7,' +i',f15.7,')'"

      integer :: in_unitp  = INPUT_UNIT  ! Unit for the ouptput files, with the ISO_FORTRAN_ENV
      integer :: out_unitp = OUTPUT_UNIT ! Unit for the input, with the ISO_FORTRAN_ENV

  END MODULE NumParameters_m
