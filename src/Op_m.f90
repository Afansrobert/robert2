module Op_m
  USE NumParameters_m
  USE diago_m
  USE Basis_m, only : Basis_t
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis

    real (kind=Rk),    allocatable :: RMat(:,:)
  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi

contains
  SUBROUTINE alloc_Op(Op,nb)
    TYPE(Op_t),  intent(inout) :: Op
    integer,     intent(in)    :: nb

    CALL dealloc_Op(Op)

    IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

    allocate(Op%RMat(nb,nb))

  END SUBROUTINE alloc_Op

  SUBROUTINE dealloc_Op(Op)
    TYPE(Op_t), intent(inout) :: Op


    nullify(Op%Basis)

    IF (allocated(Op%RMat)) THEN
      deallocate(Op%RMat)
    END IF

  END SUBROUTINE dealloc_Op

  SUBROUTINE write_Op(Op)
    USE UtilLib_m, only : Write_RMat

    TYPE(Op_t), intent(in) :: Op

    integer :: i

    IF (associated(Op%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to Op.'
    END IF

    IF (allocated(Op%RMat)) THEN
      write(out_unitp,*) 'Writing Op (real):'
      write(out_unitp,*)
      CALL Write_RMat(Op%RMat,out_unitp,5,name_info='Op%Rmat')
      write(out_unitp,*) 'END Writing Op'
    END IF

  END SUBROUTINE write_Op

  SUBROUTINE Set_Op(Op,Basis)
  USE Basis_m
  USE Molec_m

    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis
    
   ! character (len=*) :: VECPRO
    integer :: ib,iq

    real (kind=Rk), allocatable :: V(:),OpPsi_g(:),EigenVal(:),EigenVec(:,:)

    IF (.NOT. Basis_IS_allocated(Basis)) THEN
      STOP 'ERROR in Set_Op: the Basis is not initialized'
    END IF

    CALL alloc_Op(Op,Basis%nb)

    Op%Basis => Basis

    ! calculation of a potential on the grid
    allocate(V(Basis%nq))
    DO iq=1,Basis%nq
      V(iq) = Calc_pot(Basis%x(iq))
    END DO

    ! calculation of Op|b_i>
    DO ib=1,Basis%nb
      OpPsi_g = V * Basis%d0gb(:,ib) ! potential part
      OpPsi_g = OpPsi_g -HALF/mass * Basis%d2gb(:,ib,1,1) ! -1/2mass d2./dx2 part
      ! OpPsi_g is a vector on the grid. It must be projected on the basis (integration)
      OpPsi_g = OpPsi_g * Basis%w
      Op%RMat(:,ib) = matmul(transpose(Basis%d0gb),OpPsi_g)
    END DO


    CALL write_Op(Op)
    allocate(EigenVal(Basis%nb))
    allocate(EigenVec(Basis%nb,Basis%nb))
    CALL  diagonalization(Op%RMat,EigenVal,EigenVec,Basis%nb)
    
    Write(out_unitp,*)
    Write(out_unitp,*)
     Write(out_unitp,*) 'eigenvalues = '
     do ib=1,Basis%nb
        write(out_unitp,*) EigenVal(ib)
     end do
    Write(out_unitp,*)
    Write(out_unitp,*)
    
    do ib=1,Basis%nb
        write(*,*) (EigenVec(ib,iq),iq=1,Basis%nb)
     end do
    
    !CALL Write_RMat(EigenVec,out_unitp,5,name_info='Op%RMat')

  END SUBROUTINE Set_Op


  SUBROUTINE calc_OpPsi(Op,Psi,OpPsi)
    USE psi_m, ONLY : psi_t

    TYPE(Op_t),  intent(in)     :: Op
    TYPE(psi_t), intent(in)     :: Psi
    TYPE(psi_t), intent(inout)  :: OpPsi


    IF (.NOT. allocated(Op%RMat)) THEN
      STOP 'ERROR in calc_OpPsi: Op is not initialized!'
    END IF

    IF (allocated(Psi%RVec)) THEN
      OpPsi%RVec = matmul(Op%RMat,Psi%RVec)
    ELSE IF (allocated(Psi%CVec)) THEN
      OpPsi%CVec = matmul(Op%RMat,Psi%CVec)
    ELSE
      STOP 'ERROR in calc_OpPsi: Psi is not initialized!'
    END IF


  END SUBROUTINE calc_OpPsi

end module Op_m
