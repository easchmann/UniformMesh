! Fortran bindings for the uniform mesh (thin layer over mesh.c / mesh_f.c, see BINDINGS.md).
!
! Conventions (Flash-X style):
!   - block IDs, variables and axes are 1-based
!   - cell indices are global and 1-based: interior of the domain is 1..N_global per axis,
!     unused axes are 1..1, halo cells extend the block limits by MESH_NHALO
!   - mesh_f_data_ptr returns a pointer that aliases the C storage (no copy):
!       AoS layout: U(1:nvar, iloGC:ihiGC, jloGC:jhiGC, kloGC:khiGC)   (like Flash-X UG)
!       SoA layout: U(iloGC:ihiGC, jloGC:jhiGC, kloGC:khiGC, 1:nvar)   (like Flash-X UGReordered)
#include "mesh_config.h"

module mesh_f
    use, intrinsic :: iso_c_binding
    implicit none
    private

    ! kind of the mesh data, matches MESH_TYPE
    integer, parameter, public :: mesh_rk = merge(c_double, c_float, MESH_TYPE_BYTES == 8)

    integer, parameter, public :: MESH_LOW = 1, MESH_HIGH = 2
    integer, parameter, public :: MESH_LEFT_EDGE = 1, MESH_CENTER = 2, MESH_RIGHT_EDGE = 3
    integer, parameter, public :: MESH_LAYOUT_AOS = 1, MESH_LAYOUT_SOA = 2

    ! error codes of mesh_f_create
    integer, parameter, public :: MESH_OK = 0, MESH_ERR_CREATE = -1, MESH_ERR_TYPE = -2, MESH_ERR_LAYOUT = -3

    ! opaque handle; sizes are cached because they never change after creation
    type, public :: mesh_t
        private
        type(c_ptr) :: p = c_null_ptr
        integer :: layout = 0
        integer(c_int) :: n(3) = 0, nh(3) = 0, nb(3) = 0
        integer(c_int) :: nblock = 0, nvar = 0
    end type mesh_t

    public :: mesh_f_create, mesh_f_remove
    public :: mesh_f_ndim, mesh_f_nblocks, mesh_f_nvars, mesh_f_nhalo, mesh_f_layout
    public :: mesh_f_limits, mesh_f_data_ptr
    public :: mesh_f_deltas, mesh_f_domain, mesh_f_cell_coords
    public :: mesh_f_fill_halo, mesh_f_fill_halo_vars
    public :: mesh_f_c_ptr

    interface
        ! mesh.h
        function c_mesh_create(nx, ny, nz, bx, by, bz, n_var) bind(C, name="mesh_create") result(p)
            import :: c_int, c_ptr
            integer(c_int), value :: nx, ny, nz, bx, by, bz, n_var
            type(c_ptr) :: p
        end function c_mesh_create

        subroutine c_mesh_remove(p) bind(C, name="mesh_remove")
            import :: c_ptr
            type(c_ptr), value :: p
        end subroutine c_mesh_remove

        subroutine c_mesh_fill_halo(p) bind(C, name="mesh_fill_halo")
            import :: c_ptr
            type(c_ptr), value :: p
        end subroutine c_mesh_fill_halo

        function c_mesh_fill_halo_vars(p, v0, nv) bind(C, name="mesh_fill_halo_vars") result(ierr)
            import :: c_int, c_ptr
            type(c_ptr), value :: p
            integer(c_int), value :: v0, nv
            integer(c_int) :: ierr
        end function c_mesh_fill_halo_vars

        ! mesh_f.c
        function c_type_bytes() bind(C, name="mesh_c_type_bytes") result(n)
            import :: c_int
            integer(c_int) :: n
        end function c_type_bytes

        function c_ndim() bind(C, name="mesh_c_ndim") result(n)
            import :: c_int
            integer(c_int) :: n
        end function c_ndim

        subroutine c_sizes(p, n, nh, nb, n_block, n_var) bind(C, name="mesh_c_sizes")
            import :: c_int, c_ptr
            type(c_ptr), value :: p
            integer(c_int) :: n(3), nh(3), nb(3), n_block, n_var
        end subroutine c_sizes

        subroutine c_domain(p, low, dx) bind(C, name="mesh_c_domain")
            import :: c_double, c_ptr
            type(c_ptr), value :: p
            real(c_double) :: low(3), dx(3)
        end subroutine c_domain

        subroutine c_block_limits(p, b, lo, hi) bind(C, name="mesh_c_block_limits")
            import :: c_int, c_ptr
            type(c_ptr), value :: p
            integer(c_int), value :: b
            integer(c_int) :: lo(3), hi(3)
        end subroutine c_block_limits

        subroutine c_block_layout(p, b, base, stride) bind(C, name="mesh_c_block_layout")
            import :: c_int, c_int64_t, c_ptr
            type(c_ptr), value :: p
            integer(c_int), value :: b
            type(c_ptr) :: base
            integer(c_int64_t) :: stride(4)
        end subroutine c_block_layout

        subroutine c_coords(p, b, a, edge, x) bind(C, name="mesh_c_coords")
            import :: c_int, c_double, c_ptr
            type(c_ptr), value :: p
            integer(c_int), value :: b, a, edge
            real(c_double) :: x(*)
        end subroutine c_coords
    end interface

contains

    ! nx, ny, nz: global interior cells, bx, by, bz: blocks per axis, nvar: variables.
    ! ierr = MESH_OK on success, otherwise one of MESH_ERR_* and m stays unassociated.
    subroutine mesh_f_create(m, nx, ny, nz, bx, by, bz, nvar, ierr)
        type(mesh_t), intent(out) :: m
        integer, intent(in) :: nx, ny, nz, bx, by, bz, nvar
        integer, intent(out) :: ierr
        type(c_ptr) :: base
        integer(c_int64_t) :: s(4)
        integer(c_int64_t) :: nt(3)

        if (c_type_bytes() /= storage_size(1.0_mesh_rk) / 8) then
            ierr = MESH_ERR_TYPE
            return
        end if

        m%p = c_mesh_create(int(nx, c_int), int(ny, c_int), int(nz, c_int), &
                            int(bx, c_int), int(by, c_int), int(bz, c_int), int(nvar, c_int))
        if (.not. c_associated(m%p)) then
            ierr = MESH_ERR_CREATE
            return
        end if
        call c_sizes(m%p, m%n, m%nh, m%nb, m%nblock, m%nvar)

        ! the layout is decided by MESH_IDX in C; accept it if it is one of the two
        ! orders a contiguous Fortran array can represent
        call c_block_layout(m%p, 0_c_int, base, s)
        nt = m%n + 2 * m%nh
        if (s(1) == 1 .and. s(2) == m%nvar .and. s(3) == m%nvar * nt(1) &
                .and. s(4) == m%nvar * nt(1) * nt(2)) then
            m%layout = MESH_LAYOUT_AOS
        else if (s(2) == 1 .and. s(3) == nt(1) .and. s(4) == nt(1) * nt(2) &
                .and. s(1) == nt(1) * nt(2) * nt(3)) then
            m%layout = MESH_LAYOUT_SOA
        else
            call mesh_f_remove(m)
            ierr = MESH_ERR_LAYOUT
            return
        end if
        ierr = MESH_OK
    end subroutine mesh_f_create

    subroutine mesh_f_remove(m)
        type(mesh_t), intent(inout) :: m
        if (c_associated(m%p)) call c_mesh_remove(m%p)
        m%p = c_null_ptr
        m%layout = 0
    end subroutine mesh_f_remove

    ! compile-time dimension of the C side (must match Flash-X NDIM)
    integer function mesh_f_ndim()
        mesh_f_ndim = c_ndim()
    end function mesh_f_ndim

    integer function mesh_f_nblocks(m)
        type(mesh_t), intent(in) :: m
        mesh_f_nblocks = m%nblock
    end function mesh_f_nblocks

    integer function mesh_f_nvars(m)
        type(mesh_t), intent(in) :: m
        mesh_f_nvars = m%nvar
    end function mesh_f_nvars

    ! halo width on axis 1 (must match Flash-X NGUARD)
    integer function mesh_f_nhalo(m)
        type(mesh_t), intent(in) :: m
        mesh_f_nhalo = m%nh(1)
    end function mesh_f_nhalo

    integer function mesh_f_layout(m)
        type(mesh_t), intent(in) :: m
        mesh_f_layout = m%layout
    end function mesh_f_layout

    ! interior and guard-cell index limits of a block, like Flash-X blkLimits/blkLimitsGC
    subroutine mesh_f_limits(m, blockID, limits, limitsGC)
        type(mesh_t), intent(in) :: m
        integer, intent(in) :: blockID
        integer, intent(out) :: limits(MESH_LOW:MESH_HIGH, 3)
        integer, intent(out), optional :: limitsGC(MESH_LOW:MESH_HIGH, 3)
        integer(c_int) :: lo(3), hi(3)

        call c_block_limits(m%p, int(blockID - 1, c_int), lo, hi)
        limits(MESH_LOW, :) = lo + 1
        limits(MESH_HIGH, :) = hi + 1
        if (present(limitsGC)) then
            limitsGC(MESH_LOW, :) = limits(MESH_LOW, :) - m%nh
            limitsGC(MESH_HIGH, :) = limits(MESH_HIGH, :) + m%nh
        end if
    end subroutine mesh_f_limits

    ! pointer to all data of a block incl. halo, indexed with global indices (see top)
    subroutine mesh_f_data_ptr(m, blockID, U)
        type(mesh_t), intent(in) :: m
        integer, intent(in) :: blockID
        real(mesh_rk), pointer, intent(out) :: U(:, :, :, :)
        real(mesh_rk), pointer :: flat(:, :, :, :)
        type(c_ptr) :: base
        integer(c_int64_t) :: s(4)
        integer :: lim(MESH_LOW:MESH_HIGH, 3), limGC(MESH_LOW:MESH_HIGH, 3), nt(3)

        call c_block_layout(m%p, int(blockID - 1, c_int), base, s)
        call mesh_f_limits(m, blockID, lim, limGC)
        nt = limGC(MESH_HIGH, :) - limGC(MESH_LOW, :) + 1

        if (m%layout == MESH_LAYOUT_AOS) then
            call c_f_pointer(base, flat, [int(m%nvar), nt(1), nt(2), nt(3)])
            U(1:, limGC(MESH_LOW, 1):, limGC(MESH_LOW, 2):, limGC(MESH_LOW, 3):) => flat
        else
            call c_f_pointer(base, flat, [nt(1), nt(2), nt(3), int(m%nvar)])
            U(limGC(MESH_LOW, 1):, limGC(MESH_LOW, 2):, limGC(MESH_LOW, 3):, 1:) => flat
        end if
    end subroutine mesh_f_data_ptr

    subroutine mesh_f_deltas(m, deltas)
        type(mesh_t), intent(in) :: m
        real(c_double), intent(out) :: deltas(3)
        real(c_double) :: low(3)
        call c_domain(m%p, low, deltas)
    end subroutine mesh_f_deltas

    ! low and high corner of the global domain
    subroutine mesh_f_domain(m, low, high)
        type(mesh_t), intent(in) :: m
        real(c_double), intent(out) :: low(3), high(3)
        real(c_double) :: dx(3)
        call c_domain(m%p, low, dx)
        high = low + real(m%n * m%nb, c_double) * dx
    end subroutine mesh_f_domain

    ! positions on axis of all cells of a block incl. halo, like Grid_getCellCoords.
    ! coords must have limitsGC(HIGH,axis) - limitsGC(LOW,axis) + 1 elements;
    ! edge is MESH_LEFT_EDGE, MESH_CENTER or MESH_RIGHT_EDGE
    subroutine mesh_f_cell_coords(m, blockID, axis, edge, coords)
        type(mesh_t), intent(in) :: m
        integer, intent(in) :: blockID, axis, edge
        real(c_double), intent(out) :: coords(:)

        if (size(coords) /= m%n(axis) + 2 * m%nh(axis)) then
            error stop "mesh_f_cell_coords: coords has the wrong size"
        end if
        call c_coords(m%p, int(blockID - 1, c_int), int(axis - 1, c_int), int(edge - MESH_CENTER, c_int), coords)
    end subroutine mesh_f_cell_coords

    subroutine mesh_f_fill_halo(m)
        type(mesh_t), intent(in) :: m
        call c_mesh_fill_halo(m%p)
    end subroutine mesh_f_fill_halo

    ! fill the halo of variables var_first .. var_first+nvars-1 (1-based);
    ! ierr = 0, or -1 for an invalid range (nothing filled)
    subroutine mesh_f_fill_halo_vars(m, var_first, nvars, ierr)
        type(mesh_t), intent(in) :: m
        integer, intent(in) :: var_first, nvars
        integer, intent(out) :: ierr
        ierr = c_mesh_fill_halo_vars(m%p, int(var_first - 1, c_int), int(nvars, c_int))
    end subroutine mesh_f_fill_halo_vars

    ! raw C handle, for calling other C functions of the mesh directly
    type(c_ptr) function mesh_f_c_ptr(m)
        type(mesh_t), intent(in) :: m
        mesh_f_c_ptr = m%p
    end function mesh_f_c_ptr

end module mesh_f
