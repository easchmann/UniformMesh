! Tests of the Fortran bindings (mesh_f). Build with -DMESH_SOA for the SoA layout, so
! the array order here matches (Flash-X does the same with its !!Reorder tool).
#include "mesh_config.h"

#ifdef MESH_SOA
#define AT(U,v,i,j,k) U(i,j,k,v)
#define DIMS4(lo) dimension(lo(1):, lo(2):, lo(3):, 1:)
#define EXPECTED_LAYOUT MESH_LAYOUT_SOA
#else
#define AT(U,v,i,j,k) U(v,i,j,k)
#define DIMS4(lo) dimension(1:, lo(1):, lo(2):, lo(3):)
#define EXPECTED_LAYOUT MESH_LAYOUT_AOS
#endif

module test_util
    use, intrinsic :: iso_c_binding
    use mesh_f
    implicit none

    integer :: n_fail = 0

    ! C accessor, to check that Fortran and C see the same storage
    interface
        function c_mesh_get(p, b, v, i, j, k) bind(C, name="mesh_get") result(x)
            import :: c_ptr, c_int, mesh_rk
            type(c_ptr), value :: p
            integer(c_int), value :: b, v, i, j, k
            real(mesh_rk) :: x
        end function c_mesh_get
    end interface

contains

    subroutine check(cond, msg)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        if (.not. cond) then
            write(*, '(a,a)') "  FAIL: ", msg
            n_fail = n_fail + 1
        end if
    end subroutine check

    ! encoded value from global indices, exact in float for |i|,|j|,|k| < 50, v < 8
    real(mesh_rk) function code(v, i, j, k)
        integer, intent(in) :: v, i, j, k
        code = real(v * 1000000 + (i + 50) * 10000 + (j + 50) * 100 + (k + 50), mesh_rk)
    end function code

    ! periodic partner of 1-based global index g on an axis with n cells
    integer function wrap(g, n)
        integer, intent(in) :: g, n
        wrap = modulo(g - 1, n) + 1
    end function wrap

    ! Spark-shaped kernel: assumed-shape dummy with lower bounds from loGC, loop over the
    ! interior grown by g (only used axes), read the +-1 neighbours along every used axis.
    ! Built with -fcheck=bounds, so any access outside the block's storage aborts.
    subroutine neighbour_sum(Uin, limits, loGC, g, ndim, out)
        integer, intent(in) :: limits(MESH_LOW:MESH_HIGH, 3), loGC(3), g, ndim
        real(mesh_rk), DIMS4(loGC), intent(in) :: Uin
        real(mesh_rk), dimension(loGC(1):, loGC(2):, loGC(3):), intent(out) :: out
        integer :: i, j, k, gi(3)

        gi = 0
        gi(1:ndim) = g
        do k = limits(MESH_LOW, 3) - gi(3), limits(MESH_HIGH, 3) + gi(3)
        do j = limits(MESH_LOW, 2) - gi(2), limits(MESH_HIGH, 2) + gi(2)
        do i = limits(MESH_LOW, 1) - gi(1), limits(MESH_HIGH, 1) + gi(1)
            out(i, j, k) = AT(Uin, 1, i - 1, j, k) + AT(Uin, 1, i + 1, j, k)
            if (ndim > 1) out(i, j, k) = out(i, j, k) + AT(Uin, 1, i, j - 1, k) + AT(Uin, 1, i, j + 1, k)
            if (ndim > 2) out(i, j, k) = out(i, j, k) + AT(Uin, 1, i, j, k - 1) + AT(Uin, 1, i, j, k + 1)
        end do
        end do
        end do
    end subroutine neighbour_sum

end module test_util


program test_mesh_f
    use, intrinsic :: iso_c_binding
    use mesh_f
    use test_util
    implicit none

    integer, parameter :: H = MESH_NHALO
    ! global sizes, non-square, every block (3 x 3 x 2) has >= H cells
    integer, parameter :: NG(3) = [3 * (H + 2), 3 * (H + 1), 2 * (H + 1)]
    integer, parameter :: NB(3) = [3, 3, 2]
    integer, parameter :: NVAR = 3

    type(mesh_t) :: m
    real(mesh_rk), pointer :: U(:, :, :, :)
    real(mesh_rk), allocatable :: out(:, :, :)
    real(c_double), allocatable :: xc(:), xl(:), xr(:)
    real(c_double) :: deltas(3), low(3), high(3), x
    integer :: lim(MESH_LOW:MESH_HIGH, 3), limGC(MESH_LOW:MESH_HIGH, 3), loGC(3)
    integer :: ierr, ndim, b, v, i, j, k, a, g, nused(3), ngl(3), ix
    integer, allocatable :: cover(:, :, :)
    logical :: inside
    real(mesh_rk), parameter :: sentinel = -1
    real(c_double), parameter :: dom_low(3) = [-2.0_c_double, 1.0_c_double, 0.5_c_double]
    real(c_double), parameter :: dom_high(3) = [8.0_c_double, 4.0_c_double, 2.5_c_double]

    ndim = mesh_f_ndim()
    write(*, '(a,i0,a,i0,a,i0)') "variant: ndim=", ndim, " halo=", H, " layout=", EXPECTED_LAYOUT

    ! global cell counts actually used (unused axes collapse to 1)
    nused = 1
    nused(1:ndim) = NB(1:ndim)
    ngl = 1
    ngl(1:ndim) = NG(1:ndim)

    ! --- create
    write(*, '(a)') "test_create"
    call mesh_f_create(m, NG(1), NG(2), NG(3), NB(1), NB(2), NB(3), NVAR, ierr)
    call check(ierr == MESH_OK, "mesh_f_create")
    if (ierr /= MESH_OK) stop 1
    call check(mesh_f_nblocks(m) == product(nused), "nblocks")
    call check(mesh_f_nvars(m) == NVAR, "nvars")
    call check(mesh_f_nhalo(m) == H, "nhalo")
    call check(mesh_f_layout(m) == EXPECTED_LAYOUT, "layout")

    block
        type(mesh_t) :: bad
        call mesh_f_create(bad, NG(1) + 1, NG(2), NG(3), 3, 1, 1, NVAR, ierr)
        call check(ierr == MESH_ERR_CREATE, "invalid sizes rejected")
    end block

    ! --- limits tile the global index space exactly
    write(*, '(a)') "test_limits"
    allocate(cover(ngl(1), ngl(2), ngl(3)))
    cover = 0
    do b = 1, mesh_f_nblocks(m)
        call mesh_f_limits(m, b, lim, limGC)
        do a = 1, 3
            call check(lim(MESH_HIGH, a) - lim(MESH_LOW, a) + 1 == ngl(a) / nused(a), "block size")
            call check(limGC(MESH_LOW, a) == lim(MESH_LOW, a) - merge(H, 0, a <= ndim), "limitsGC low")
            call check(limGC(MESH_HIGH, a) == lim(MESH_HIGH, a) + merge(H, 0, a <= ndim), "limitsGC high")
        end do
        cover(lim(1, 1):lim(2, 1), lim(1, 2):lim(2, 2), lim(1, 3):lim(2, 3)) = &
            cover(lim(1, 1):lim(2, 1), lim(1, 2):lim(2, 2), lim(1, 3):lim(2, 3)) + 1
    end do
    call check(all(cover == 1), "blocks cover the domain exactly once")

    ! --- data pointer: bounds, and Fortran writes are what C reads
    write(*, '(a)') "test_data_ptr"
    do b = 1, mesh_f_nblocks(m)
        call mesh_f_limits(m, b, lim, limGC)
        call mesh_f_data_ptr(m, b, U)
#ifdef MESH_SOA
        call check(all(lbound(U) == [limGC(MESH_LOW, :), 1]), "lbound")
        call check(all(ubound(U) == [limGC(MESH_HIGH, :), NVAR]), "ubound")
#else
        call check(all(lbound(U) == [1, limGC(MESH_LOW, :)]), "lbound")
        call check(all(ubound(U) == [NVAR, limGC(MESH_HIGH, :)]), "ubound")
#endif
        do k = limGC(1, 3), limGC(2, 3)
        do j = limGC(1, 2), limGC(2, 2)
        do i = limGC(1, 1), limGC(2, 1)
        do v = 1, NVAR
            AT(U, v, i, j, k) = code(v, i, j, k)
        end do
        end do
        end do
        end do
        ! C local index = global - first interior global index
        do k = limGC(1, 3), limGC(2, 3)
        do j = limGC(1, 2), limGC(2, 2)
        do i = limGC(1, 1), limGC(2, 1)
        do v = 1, NVAR
            x = c_mesh_get(mesh_f_c_ptr(m), int(b - 1, c_int), int(v - 1, c_int), int(i - lim(1, 1), c_int), &
                           int(j - lim(1, 2), c_int), int(k - lim(1, 3), c_int))
            call check(x == code(v, i, j, k), "Fortran write seen by C mesh_get")
        end do
        end do
        end do
        end do
    end do

    ! --- halo fill through the binding: halos hold the periodic partner (global wrap)
    write(*, '(a)') "test_fill_halo"
    call set_halo(sentinel)
    call mesh_f_fill_halo(m)
    call check_halo(1, NVAR)

    write(*, '(a)') "test_fill_halo_vars"
    call mesh_f_fill_halo_vars(m, 0, 1, ierr)
    call check(ierr == -1, "var_first = 0 rejected")
    call mesh_f_fill_halo_vars(m, NVAR, 2, ierr)
    call check(ierr == -1, "range past nvar rejected")
    call set_halo(sentinel)
    call mesh_f_fill_halo_vars(m, 2, 1, ierr)
    call check(ierr == 0, "fill var 2")
    call check_halo(2, 1)

    ! --- geometry
    write(*, '(a)') "test_coords"
    ! non-unit domain, so dx and positions are not just 1/N
    call mesh_f_set_domain(m, dom_low, dom_high, ierr)
    call check(ierr == 0, "mesh_f_set_domain")
    call mesh_f_set_domain(m, [0.0_c_double, 0.0_c_double, 0.0_c_double], &
                           [-1.0_c_double, 1.0_c_double, 1.0_c_double], ierr)
    call check(ierr == -1, "invalid domain rejected")
    call mesh_f_deltas(m, deltas)
    call mesh_f_domain(m, low, high)
    do a = 1, 3
        call check(abs(deltas(a) - (high(a) - low(a)) / ngl(a)) < 1e-12_c_double, "deltas")
        call check(abs(low(a) - dom_low(a)) < 1e-12_c_double, "domain low")
        call check(abs(high(a) - dom_high(a)) < 1e-12_c_double, "domain high")
    end do
    do b = 1, mesh_f_nblocks(m)
        call mesh_f_limits(m, b, lim, limGC)
        do a = 1, 3
            allocate(xc(limGC(1, a):limGC(2, a)), xl(limGC(1, a):limGC(2, a)), xr(limGC(1, a):limGC(2, a)))
            call mesh_f_cell_coords(m, b, a, MESH_CENTER, xc)
            call mesh_f_cell_coords(m, b, a, MESH_LEFT_EDGE, xl)
            call mesh_f_cell_coords(m, b, a, MESH_RIGHT_EDGE, xr)
            do ix = limGC(1, a), limGC(2, a)
                ! Flash-X formulas (Grid_getCenterCoords etc.) with 1-based global index
                call check(abs(xc(ix) - (low(a) + (ix - 0.5_c_double) * deltas(a))) < 1e-12_c_double, "center")
                call check(abs(xl(ix) - (low(a) + (ix - 1) * deltas(a))) < 1e-12_c_double, "left edge")
                call check(abs(xr(ix) - (low(a) + ix * deltas(a))) < 1e-12_c_double, "right edge")
            end do
            deallocate(xc, xl, xr)
        end do
    end do

    ! --- Spark-shaped kernel over grown boxes, bounds-checked
    write(*, '(a)') "test_kernel"
    call mesh_f_fill_halo(m)
    do b = 1, mesh_f_nblocks(m)
        call mesh_f_limits(m, b, lim, limGC)
        call mesh_f_data_ptr(m, b, U)
        loGC = limGC(MESH_LOW, :) ! as in Spark's driver
        allocate(out(limGC(1, 1):limGC(2, 1), limGC(1, 2):limGC(2, 2), limGC(1, 3):limGC(2, 3)))
        do g = 0, H - 1
            out = 0
            call neighbour_sum(U, lim, loGC, g, ndim, out)
            do k = lim(1, 3) - merge(g, 0, ndim > 2), lim(2, 3) + merge(g, 0, ndim > 2)
            do j = lim(1, 2) - g, lim(2, 2) + g
            do i = lim(1, 1) - g, lim(2, 1) + g
                x = expected_sum(i, j, k)
                call check(out(i, j, k) == x, "neighbour sum")
            end do
            end do
            end do
        end do
        deallocate(out)
    end do

    call mesh_f_remove(m)

    if (n_fail == 0) then
        write(*, '(a)') "ALL PASSED"
    else
        write(*, '(i0,a)') n_fail, " CHECK(S) FAILED"
        stop 1
    end if

contains

    ! write the code into every interior cell and val into every halo cell
    subroutine set_halo(val)
        real(mesh_rk), intent(in) :: val
        do b = 1, mesh_f_nblocks(m)
            call mesh_f_limits(m, b, lim, limGC)
            call mesh_f_data_ptr(m, b, U)
            do k = limGC(1, 3), limGC(2, 3)
            do j = limGC(1, 2), limGC(2, 2)
            do i = limGC(1, 1), limGC(2, 1)
                inside = all([i, j, k] >= lim(1, :)) .and. all([i, j, k] <= lim(2, :))
                do v = 1, NVAR
                    AT(U, v, i, j, k) = merge(code(v, i, j, k), val, inside)
                end do
            end do
            end do
            end do
        end do
    end subroutine set_halo

    ! halo of vars v0..v0+nv-1 holds the partner's code, all other halos still hold the sentinel
    subroutine check_halo(v0, nv)
        integer, intent(in) :: v0, nv
        do b = 1, mesh_f_nblocks(m)
            call mesh_f_limits(m, b, lim, limGC)
            call mesh_f_data_ptr(m, b, U)
            do k = limGC(1, 3), limGC(2, 3)
            do j = limGC(1, 2), limGC(2, 2)
            do i = limGC(1, 1), limGC(2, 1)
                inside = all([i, j, k] >= lim(1, :)) .and. all([i, j, k] <= lim(2, :))
                do v = 1, NVAR
                    if (inside .or. (v0 <= v .and. v < v0 + nv)) then
                        call check(AT(U, v, i, j, k) == code(v, wrap(i, ngl(1)), wrap(j, ngl(2)), wrap(k, ngl(3))), &
                                   "halo holds periodic partner")
                    else
                        call check(AT(U, v, i, j, k) == sentinel, "halo of unfilled variable untouched")
                    end if
                end do
            end do
            end do
            end do
        end do
    end subroutine check_halo

    real(mesh_rk) function expected_sum(i, j, k)
        integer, intent(in) :: i, j, k
        expected_sum = code(1, wrap(i - 1, ngl(1)), wrap(j, ngl(2)), wrap(k, ngl(3))) &
                     + code(1, wrap(i + 1, ngl(1)), wrap(j, ngl(2)), wrap(k, ngl(3)))
        if (ndim > 1) expected_sum = expected_sum &
                     + code(1, wrap(i, ngl(1)), wrap(j - 1, ngl(2)), wrap(k, ngl(3))) &
                     + code(1, wrap(i, ngl(1)), wrap(j + 1, ngl(2)), wrap(k, ngl(3)))
        if (ndim > 2) expected_sum = expected_sum &
                     + code(1, wrap(i, ngl(1)), wrap(j, ngl(2)), wrap(k - 1, ngl(3))) &
                     + code(1, wrap(i, ngl(1)), wrap(j, ngl(2)), wrap(k + 1, ngl(3)))
    end function expected_sum

end program test_mesh_f
