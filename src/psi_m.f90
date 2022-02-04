module psi_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none

  private

  TYPE :: psi_t
    TYPE (Basis_t),    pointer     :: Basis
    real (kind=Rk),    allocatable :: RVec(:)
    complex (kind=Rk), allocatable :: CVec(:)
  END TYPE psi_t

  public :: psi_t,write_psi,init_psi,dealloc_psi
  ! operation on psi has to be defined: psi=psi1, psi1+psi2, psi=psi1*cte ...
contains
  SUBROUTINE init_psi(psi,Basis,cplx)
  USE Basis_m

    TYPE(psi_t),    intent(inout)      :: psi
    TYPE (Basis_t), intent(in), target :: Basis
    logical,        intent(in)         :: cplx

    CALL dealloc_psi(psi)

    IF (.NOT. Basis_IS_allocated(Basis)) STOP 'ERROR in Set_Op: the Basis is not initialized'


    IF (Basis%nb < 1) STOP 'ERROR in init_psi: Basis%nb < 1!'

    psi%Basis => Basis

    IF (cplx) THEN
      allocate(psi%CVec(Basis%nb))
    ELSE
      allocate(psi%RVec(Basis%nb))
    END IF

  END SUBROUTINE init_psi

  SUBROUTINE dealloc_psi(psi)
    TYPE(psi_t), intent(inout) :: psi

    nullify(psi%Basis)

    IF (allocated(psi%RVec)) THEN
      deallocate(psi%RVec)
    END IF

    IF (allocated(psi%CVec)) THEN
      deallocate(psi%CVec)
    END IF

  END SUBROUTINE dealloc_psi

  SUBROUTINE write_psi(psi)
    TYPE(psi_t), intent(in) :: psi

    IF (associated(psi%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to psi.'
    END IF

    IF (allocated(psi%RVec)) THEN
      write(out_unitp,*) 'Writing psi (real):'
      write(out_unitp,*) psi%RVec
      write(out_unitp,*) 'END Writing psi'
    END IF
    IF (allocated(psi%CVec)) THEN
      write(out_unitp,*) 'Writing psi (complex):'
      write(out_unitp,*) psi%CVec
      write(out_unitp,*) 'END Writting psi'
    END IF
  END SUBROUTINE write_psi
end module psi_m
