! Fortran bindings for the uniform mesh (thin layer over mesh.c / mesh_bind.c, see BINDINGS.md).
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
    public :: mesh_f_limits, mesh_f_data_ptr, mesh_f_field_ptr
    public :: mesh_f_neighbor
    public :: mesh_f_set_domain, mesh_f_deltas, mesh_f_domain, mesh_f_cell_coords
    public :: mesh_f_fill_halo, mesh_f_fill_halo_vars
    public :: mesh_f_write_dump
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

        function c_mesh_set_domain(p, low, high) bind(C, name="mesh_set_domain") result(ierr)
            import :: c_int, c_double, c_ptr
            type(c_ptr), value :: p
            real(c_double), intent(in) :: low(3), high(3)
            integer(c_int) :: ierr
        end function c_mesh_set_domain

        ! mesh_bind.c
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

        function c_mesh_block_neighbor(p, b, a, dir) bind(C, name="mesh_c_block_neighbor") result(n)
            import :: c_int, c_ptr
            type(c_ptr), value :: p
            integer(c_int), value :: b, a, dir
            integer(c_int) :: n
        end function c_mesh_block_neighbor
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

    ! AoS only: Fortran pointer to all blocks of the whole mesh at once, so a solver
    ! can view it as a single array U(var, 1-nh1:, 1-nh2:, 1-nh3:, block) with block-local
    ! halo lower bounds. This is the layout NewImpl's block-loop kernels index directly.
    subroutine mesh_f_field_ptr(m, U)
        type(mesh_t), intent(in) :: m
        real(mesh_rk), pointer, intent(out) :: U(:, :, :, :, :)
        real(mesh_rk), pointer :: flat(:, :, :, :, :)
        type(c_ptr) :: base
        integer(c_int64_t) :: s(4)
        integer :: lim(MESH_LOW:MESH_HIGH, 3), limGC(MESH_LOW:MESH_HIGH, 3), nt(3)

        if (m%layout /= MESH_LAYOUT_AOS) error stop "mesh_f_field_ptr: AoS layout required"
        call c_block_layout(m%p, 0_c_int, base, s)
        call mesh_f_limits(m, 1, lim, limGC)
        nt = limGC(MESH_HIGH, :) - limGC(MESH_LOW, :) + 1
        call c_f_pointer(base, flat, [int(m%nvar), nt(1), nt(2), nt(3), m%nblock])
        U(1:, limGC(MESH_LOW, 1):, limGC(MESH_LOW, 2):, limGC(MESH_LOW, 3):, 1:) => flat
    end subroutine mesh_f_field_ptr

    ! neighbor block of blockID along axis (1..3) in direction dir (-1 low, +1 high);
    ! the periodic block-neighbor logic lives in mesh.c (mesh_block_neighbor)
    integer function mesh_f_neighbor(m, blockID, axis, dir)
        type(mesh_t), intent(in) :: m
        integer, intent(in) :: blockID, axis, dir
        mesh_f_neighbor = int(c_mesh_block_neighbor(m%p, int(blockID, c_int), &
                                                    int(axis - 1, c_int), int(dir, c_int)), kind(blockID))
    end function mesh_f_neighbor

    ! set the global domain [low, high] per axis (default [0,1]); dx follows from it.
    ! ierr = 0, or -1 if high <= low on a used axis (nothing changed)
    subroutine mesh_f_set_domain(m, low, high, ierr)
        type(mesh_t), intent(in) :: m
        real(c_double), intent(in) :: low(3), high(3)
        integer, intent(out) :: ierr
        ierr = c_mesh_set_domain(m%p, low, high)
    end subroutine mesh_f_set_domain

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

    ! Binary dump of the whole mesh (all blocks incl. halo), read by flashx/plot_mesh.py.
    ! Stream access, native endianness:
    !   char(8) "UMDUMP01"
    !   int32   ndim, nvar, nblock, step, layout (1 AoS, 2 SoA), bytes per value
    !   int32   global cells(3), blocks per axis(3), halo(3)
    !   float64 time, domain low(3), domain high(3)
    !   char(8) name of each variable (blank padded)
    !   per block: int32 interior limits (low,high) x 3 axes, global 1-based,
    !              mesh data of the block incl. halo, first index fastest
    !              (AoS: (nvar, i, j, k), SoA: (i, j, k, nvar); MESH_TYPE_BYTES per value)
    ! ierr = 0, or the iostat of the failed open/write
    subroutine mesh_f_write_dump(m, filename, names, step, time, ierr)
        type(mesh_t), intent(in) :: m
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: names(:)
        integer, intent(in) :: step
        real(c_double), intent(in) :: time
        integer, intent(out) :: ierr
        real(mesh_rk), pointer :: U(:, :, :, :)
        real(c_double) :: low(3), high(3)
        integer :: lim(MESH_LOW:MESH_HIGH, 3), iu, b, v
        character(len=8) :: name8

        call mesh_f_domain(m, low, high)
        open(newunit=iu, file=filename, access="stream", form="unformatted", &
             status="replace", action="write", iostat=ierr)
        if (ierr /= 0) return
        write(iu, iostat=ierr) "UMDUMP01", &
            int([mesh_f_ndim(), int(m%nvar), int(m%nblock), step, m%layout, MESH_TYPE_BYTES], c_int32_t), &
            int(m%n * m%nb, c_int32_t), int(m%nb, c_int32_t), int(m%nh, c_int32_t), &
            time, low, high
        do v = 1, m%nvar
            name8 = ""
            if (v <= size(names)) name8 = names(v)
            if (ierr == 0) write(iu, iostat=ierr) name8
        end do
        do b = 1, m%nblock
            call mesh_f_limits(m, b, lim)
            call mesh_f_data_ptr(m, b, U)
            if (ierr == 0) write(iu, iostat=ierr) int(lim, c_int32_t), U
        end do
        close(iu)
    end subroutine mesh_f_write_dump

    ! raw C handle, for calling other C functions of the mesh directly
    type(c_ptr) function mesh_f_c_ptr(m)
        type(mesh_t), intent(in) :: m
        mesh_f_c_ptr = m%p
    end function mesh_f_c_ptr

end module mesh_f
