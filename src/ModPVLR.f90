! ModPVLR.f90
!
! This modules provides C declarations for binding with C/C++.

module pvlr
    use, intrinsic :: iso_c_binding

    use ModAppr, only : appr_aca => aca, appr_maxvol => maxvol, maxvolproj
    use ModVec
    use ModMtrx

    implicit none

    type, bind(C) :: Matrix
        type(c_ptr) :: ptr
    end type

    type, bind(C) :: IntVector
        type(c_ptr) :: ptr
    end type

contains

    subroutine matrix_create(this, n, m) bind(C,  name="pvlr_matrix_create")
        type(Matrix), intent(inout) :: this
        integer(c_int), value, intent(in) :: n, m
        type(Mtrx), pointer :: ptr

        allocate(ptr)
        call ptr%init(n, m)
        this%ptr = c_loc(ptr)
    end subroutine

    subroutine matrix_destroy(this) bind(C, name="pvlr_matrix_destroy")
        type(Matrix), intent(inout) :: this
        type(Mtrx), pointer :: ptr

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            call ptr%deinit()
            deallocate(ptr)
            this%ptr = c_null_ptr
        end if
    end subroutine

    pure function matrix_at(i, j, mat) result(elem)
        integer(4), intent(in) :: i, j
        type(Mtrx), intent(in) :: mat
        double precision :: elem
        elem = mat%d(i,j)
    end function


    subroutine matrix_shape(this, rows, cols) bind(C, name="pvlr_matrix_shape")
        type(Matrix), intent(in) :: this
        integer(c_int), intent(out) :: rows, cols

        type(Mtrx), pointer ::ptr

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            rows = ptr%n
            cols = ptr%m
        else
            rows = 0
            cols = 0
        end if
    end subroutine

    subroutine matrix_from_buffer(           &
        this, buffer                         &
    ) bind(C,  name="pvlr_matrix_from_buffer")
        type(Matrix), intent(inout) :: this
        real(c_double), dimension(*), intent(in) :: buffer

        type(Mtrx), pointer :: ptr
        integer(4) :: i, j

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            do i = 1, ptr%n
                do j = 1, ptr%m
                    ! Coverts C-contingues to F-contingues arrays.
                    ptr%d(i, j) = buffer((i - 1) * ptr%m + j)
                end do
            end do
        end if
    end subroutine

    subroutine matrix_to_buffer(           &
        this, buffer                       &
    ) bind(C,  name="pvlr_matrix_to_buffer")
        type(Matrix), intent(inout) :: this
        real(c_double), dimension(*), intent(out) :: buffer

        type(Mtrx), pointer :: ptr
        integer(4) :: i, j

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            do i = 1, ptr%n
                do j = 1, ptr%m
                    ! Coverts F-contingues to C-contingues arrays.
                    buffer((i - 1) * ptr%m + j) = ptr%d(i, j)
                end do
            end do
        end if
    end subroutine

    subroutine int_vector_create(this, n                                 &
                                 ) bind(C,  name="pvlr_int_vector_create")
        type(IntVector), intent(inout) :: this
        integer(c_int), value, intent(in) :: n
        type(IntVec), pointer :: ptr

        allocate(ptr)
        call ptr%init(n)
        this%ptr = c_loc(ptr)
    end subroutine

    subroutine int_vector_destroy(this) bind(C, name="pvlr_int_vector_destroy")
        type(IntVector), intent(inout) :: this
        type(IntVec), pointer :: ptr

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            call ptr%deinit()
            deallocate(ptr)
            this%ptr = c_null_ptr
        end if
    end subroutine

    subroutine int_vector_shape(           &
        this, size                         &
    ) bind(C,  name="pvlr_int_vector_shape")
        type(IntVector), intent(in) :: this
        integer(c_int), intent(out) :: size

        type(IntVec), pointer :: ptr

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            size = ptr%n
        else
            size = 0
        end if
    end subroutine

    subroutine int_vector_from_buffer(           &
        this, buffer                             &
    ) bind(C,  name="pvlr_int_vector_from_buffer")
        type(IntVector), intent(inout) :: this
        integer(c_int), dimension(*), intent(in) :: buffer

        type(IntVec), pointer :: ptr
        integer(4) :: i

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            do i = 1, ptr%n
                ptr%d(i) = buffer(i)
            end do
        end if
    end subroutine

    subroutine int_vector_to_buffer(           &
        this, buffer                           &
    ) bind(C,  name="pvlr_int_vector_to_buffer")
        type(IntVector), intent(in) :: this
        integer(c_int), dimension(*), intent(out) :: buffer

        type(IntVec), pointer :: ptr
        integer(4) :: i

        if (c_associated(this%ptr)) then
            call c_f_pointer(this%ptr, ptr)
            do i = 1, ptr%n
                buffer(i) = ptr%d(i)
            end do
        end if
    end subroutine

    subroutine aca(A, Ni, Nj, MaxRank, jpmax_, U, V, per, rho_, rel_err, &
                   rel_err_fro, abs_err, iNs_, jNs_, maxsteps            &
                   ) bind(C, name="pvlr_aca")
        type(Matrix), intent(in) :: A
        integer(c_int), value, intent(in) :: Ni, Nj
        integer(c_int), value, intent(in) :: MaxRank
        integer(c_int), value, intent(in) :: jpmax_
        type(Matrix), intent(out) :: U, V
        type(IntVector), intent(in) :: per
        real(c_double), value, intent(in) :: rho_
        real(c_double), value, intent(in) :: rel_err
        real(c_double), value, intent(in) :: rel_err_fro
        real(c_double), value, intent(in) :: abs_err
        type(IntVector), intent(out) :: iNs_, jNs_
        integer(c_int), value, intent(in) :: maxsteps

        type(Mtrx), pointer :: ptr_A, ptr_U, ptr_V
        type(IntVec), pointer :: ptr_per, ptr_iNs_, ptr_jNs_

        if (.not. c_associated(A%ptr) .or. &
            .not. c_associated(per%ptr)) then
            return
        end if

        call c_f_pointer(A%ptr, ptr_A)
        call c_f_pointer(per%ptr, ptr_per)

        allocate(ptr_U)
        allocate(ptr_V)
        allocate(ptr_iNs_)
        allocate(ptr_jNs_)

        call appr_aca(matrix_at, ptr_A, Ni, Nj, MaxRank, jpmax_,         &
                      ptr_U, ptr_V, ptr_per, rho_, rel_err, rel_err_fro, &
                      abs_err, ptr_iNs_, ptr_jNs_, maxsteps)

        U%ptr = c_loc(ptr_U)
        V%ptr = c_loc(ptr_V)
        iNs_%ptr = c_loc(ptr_iNs_)
        jNs_%ptr = c_loc(ptr_jNs_)
    end subroutine

    subroutine maxvol(A, M, N, rank, per1, per2, C, UR, maxsteps, maxswaps, &
                      CA, premaxvol) bind(C, name="pvlr_maxvol")
        type(Matrix), intent(in) :: A
        integer(c_int), value, intent(in) :: M, N, rank
        type(IntVector) :: per1, per2
        type(Matrix), intent(out) :: C
        type(Matrix), intent(out) :: UR
        integer(c_int), value, intent(in) :: maxsteps, maxswaps
        type(Matrix), intent(out) :: CA
        logical(c_bool), value, intent(in) :: premaxvol

        type(Mtrx), pointer :: ptr_A, ptr_C, ptr_UR, ptr_CA
        type(IntVec), pointer :: ptr_per1, ptr_per2
        logical :: premaxvol_

        if (.not. c_associated(A%ptr) .or. &
            .not. c_associated(per1%ptr) .or. &
            .not. c_associated(per2%ptr)) then
            return
        end if

        call c_f_pointer(A%ptr, ptr_A)
        call c_f_pointer(per1%ptr, ptr_per1)
        call c_f_pointer(per2%ptr, ptr_per2)

        allocate(ptr_C)
        allocate(ptr_UR)
        allocate(ptr_CA)

        premaxvol_ = premaxvol

        call appr_maxvol(matrix_at, M, N, rank, ptr_per1, ptr_per2, ptr_A,    &
                         ptr_C, ptr_UR, maxsteps, maxswaps, ptr_CA, premaxvol_)

        C%ptr = c_loc(ptr_C)
        UR%ptr = c_loc(ptr_UR)
        CA%ptr = c_loc(ptr_CA)
    end subroutine

    subroutine maxvol_proj(A, M, N, rank, k, l, per1, per2, C, UR, &
                           maxsteps, maxswaps, acc_type            &
                           ) bind(C, name="pvlr_maxvol_proj")
        type(Matrix), intent(in) :: A
        integer(c_int), intent(in) :: M, N
        integer(c_int), intent(in) :: rank
        integer(c_int), intent(in) :: k, l
        type(IntVector) :: per1, per2
        type(Matrix), intent(out) :: C
        type(Matrix), intent(out) :: UR
        integer(c_int), intent(in) :: maxsteps, maxswaps
        integer(c_int), intent(in) :: acc_type

        type(Mtrx), pointer :: ptr_A, ptr_C, ptr_UR
        type(IntVec), pointer :: ptr_per1, ptr_per2


        if (.not. c_associated(A%ptr) .or. &
            .not. c_associated(per1%ptr) .or. &
            .not. c_associated(per2%ptr)) then
            return
        end if

        call c_f_pointer(A%ptr, ptr_A)
        call c_f_pointer(per1%ptr, ptr_per1)
        call c_f_pointer(per2%ptr, ptr_per2)

        call c_f_pointer(C%ptr, ptr_C)
        call c_f_pointer(UR%ptr, ptr_UR)

        call maxvolproj(matrix_at, M, N, rank, k, l, ptr_per1, ptr_per2,  &
                        ptr_A, ptr_C, ptr_UR, maxsteps, maxswaps, acc_type)

        C%ptr = c_loc(ptr_C)
        UR%ptr = c_loc(ptr_UR)
  end subroutine

end module pvlr
