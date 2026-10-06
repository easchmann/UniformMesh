!> Glue between Flash-X and the UniformMesh for the NewImpl sub-unit
!! (see flashx/README.md in the UniformMesh repo and NewImpl/README.md on MOL).
!!
!! Flash-X keeps owning the solution (unk of the UG Grid), which initialization, IO and
!! the timestep use. Hydro copies the interior into mesh 1 (state), the NewImpl whole-step
!! driver works on meshes 1..3 (state/reference/stage) and Hydro copies the result back.
!!
!! Requirements, checked in um_init / um_copy:
!!   UG Grid on one MPI rank, periodic boundaries on all used axes, NDIM equal to
!!   MESH_NDIM of mesh_config.h, AoS layout (Flash-X without -index-reorder).
!!   Unlike the wrapper-path Spark/UniformMesh, MESH_NHALO only needs to be >= 3
!!   (NewImpl's NGUARD>=3 stencil), not equal to NGUARD.
!!
!! With um_dumpInterval = N > 0, um_dump writes the mesh (all blocks incl. halo, all variables)
!! to um_dump_<step>.bin at the start and every N steps, for flashx/plot_mesh.py.

#include "Simulation.h"
#include "constants.h"

module um_flashx
    use mesh_f
    implicit none
    private

    type(mesh_t), save, public :: um_mesh(3)   ! 1=state, 2=reference, 3=stage
    logical, save :: um_ready = .false.
    integer, save :: um_gsize(MDIM)   ! global interior cells per axis
    integer, save :: um_dumpInterval = 0
    integer, save :: um_nSteps = 0    ! Hydro calls so far

    public :: um_init, um_copyIn, um_copyOut, um_output

contains

    ! create the 3 meshes on the first call; aborts if the Flash-X setup does not fit
    subroutine um_init()
        use Grid_interface, ONLY: Grid_getGlobalIndexLimits, Grid_getDomainBoundBox
        use RuntimeParameters_interface, ONLY: RuntimeParameters_get
        use Driver_interface, ONLY: Driver_abort

        character(len=16), parameter :: bcNames(2*MDIM) = [character(len=16) :: &
            "xl_boundary_type", "xr_boundary_type", "yl_boundary_type", &
            "yr_boundary_type", "zl_boundary_type", "zr_boundary_type"]
        character(len=MAX_STRING_LENGTH) :: bc
        real :: bbox(LOW:HIGH, MDIM)
        integer :: nb(MDIM), ierr, i, m

        if (um_ready) return

        if (mesh_f_ndim() /= NDIM) then
            call Driver_abort("[UniformMeshNewImpl] mesh_config.h was generated for another NDIM")
        end if
        do i = 1, 2*NDIM
            call RuntimeParameters_get(trim(bcNames(i)), bc)
            if (trim(bc) /= "periodic") then
                call Driver_abort("[UniformMeshNewImpl] only periodic boundaries are supported: " // trim(bcNames(i)))
            end if
        end do

        call Grid_getGlobalIndexLimits(um_gsize)
        call RuntimeParameters_get("um_nblockx", nb(IAXIS))
        call RuntimeParameters_get("um_nblocky", nb(JAXIS))
        call RuntimeParameters_get("um_nblockz", nb(KAXIS))
        call RuntimeParameters_get("um_dumpInterval", um_dumpInterval)

        ! the three meshes are identical geometry; state(1) holds the Flash-X solution,
        ! reference(2)/stage(3) are NewImpl RK scratch filled by hy_copyBlockState
        do m = 1, 3
            call mesh_f_create(um_mesh(m), um_gsize(IAXIS), um_gsize(JAXIS), um_gsize(KAXIS), &
                               nb(IAXIS), nb(JAXIS), nb(KAXIS), NUNK_VARS, ierr)
            if (ierr /= MESH_OK) then
                call Driver_abort("[UniformMeshNewImpl] mesh_f_create failed: um_nblock[xyz] must divide the " // &
                                  "global cell counts and every block needs >= 3 cells per axis")
            end if
            ! same physical domain as Flash-X, so the mesh's dx (used by the solver) is Flash-X's dx
            call Grid_getDomainBoundBox(bbox)
            call mesh_f_set_domain(um_mesh(m), bbox(LOW, :), bbox(HIGH, :), ierr)
            if (ierr /= 0) then
                call Driver_abort("[UniformMeshNewImpl] invalid domain bounding box")
            end if
        end do

        ! NewImpl needs at least 3 guard cells per axis (the driver hard-asserts
        ! guardLo <= lo-3 and guardHi >= hi+3); mesh_config.h is built with halo >= 3.
        ! The wrapper-path UniformMesh requires MESH_NHALO == NGUARD, but this sub-unit
        ! accepts any halo >= 3 (the mesh still fills all of its own halo cells).
        if (mesh_f_nhalo(um_mesh(1)) < 3) then
            call Driver_abort("[UniformMeshNewImpl] MESH_NHALO < 3, regenerate mesh_config.h (install_newimpl.sh)")
        end if
        if (mesh_f_layout(um_mesh(1)) /= MESH_LAYOUT_AOS) then
            call Driver_abort("[UniformMeshNewImpl] needs the AoS layout, like Flash-X's unk(var,i,j,k)")
        end if

        write(*, '(a,3i6,a,3i4,a,i3)') " [UniformMeshNewImpl] global cells", um_gsize, &
            ", blocks", nb, ", halo", mesh_f_nhalo(um_mesh(1))
        um_ready = .true.
    end subroutine um_init

    ! called by Hydro before (afterStep = .false.) and after (.true.) each step;
    ! writes a dump at the start and every um_dumpInterval steps
    subroutine um_output(time, afterStep)
        real, intent(in) :: time
        logical, intent(in) :: afterStep

        if (um_dumpInterval <= 0) return
        if (.not. afterStep) then
            if (um_nSteps == 0) call um_dump(0, time)     ! initial state, halo already filled
        else
            um_nSteps = um_nSteps + 1
            if (mod(um_nSteps, um_dumpInterval) == 0) then
                call mesh_f_fill_halo(um_mesh(1))         ! so the dumped halo matches the interior
                call um_dump(um_nSteps, time)
            end if
        end if
    end subroutine um_output

    ! dump of the mesh (format: mesh_f_write_dump in mesh_f.F90) with Flash-X's variable names
    subroutine um_dump(step, time)
        use iso_c_binding, ONLY: c_double
        use Simulation_interface, ONLY: Simulation_mapIntToStr
        use Driver_interface, ONLY: Driver_abort

        integer, intent(in) :: step
        real, intent(in) :: time
        character(len=8) :: names(NUNK_VARS)
        character(len=MAX_STRING_LENGTH) :: name
        character(len=32) :: fname
        integer :: v, ierr

        do v = 1, NUNK_VARS
            name = ""
            call Simulation_mapIntToStr(v, name, MAPBLOCK_UNK)
            names(v) = name
        end do
        write(fname, '(a,i6.6,a)') "um_dump_", step, ".bin"
        call mesh_f_write_dump(um_mesh(1), trim(fname), names, step, real(time, c_double), ierr)
        if (ierr /= 0) call Driver_abort("[UniformMeshNewImpl] could not write " // trim(fname))
        write(*, '(a,a,a,es12.5)') " [UniformMeshNewImpl] wrote ", trim(fname), ", t =", time
    end subroutine um_dump

    ! Flash-X Grid -> mesh 1 (state), interior cells of all variables
    subroutine um_copyIn()
        call um_copy(.true.)
    end subroutine um_copyIn

    ! mesh 1 (state) -> Flash-X Grid, interior cells of all variables
    subroutine um_copyOut()
        call um_copy(.false.)
    end subroutine um_copyOut

    subroutine um_copy(toMesh)
        use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
        use Grid_iterator, ONLY: Grid_iterator_t
        use Grid_tile, ONLY: Grid_tile_t
        use Driver_interface, ONLY: Driver_abort

        logical, intent(in) :: toMesh
        type(Grid_iterator_t) :: itor
        type(Grid_tile_t) :: tileDesc
        real, dimension(:,:,:,:), pointer :: solnData
        real(mesh_rk), dimension(:,:,:,:), pointer :: U
        real :: gridDel(MDIM), meshDel(MDIM)
        integer :: lim(LOW:HIGH, MDIM), off(MDIM), lo(MDIM), hi(MDIM)
        integer :: b, ntiles

        nullify(solnData)
        ntiles = 0
        call Grid_getTileIterator(itor, LEAF, tiling=.false.)
        do while (itor%isValid())
            call itor%currentTile(tileDesc)
            ntiles = ntiles + 1

            ! one rank: the single Grid block is the whole domain, global index = local + off
            off = 1 - tileDesc%limits(LOW, :)
            if (any(tileDesc%limits(HIGH, :) + off /= um_gsize)) then
                call Driver_abort("[UniformMeshNewImpl] the Grid block is not the whole domain (run on 1 MPI rank)")
            end if

            ! Spark takes dx from the mesh, so it must be Flash-X's dx
            call tileDesc%deltas(gridDel)
            call mesh_f_deltas(um_mesh(1), meshDel)
            if (any(abs(meshDel(1:NDIM) - gridDel(1:NDIM)) > 1e-12 * gridDel(1:NDIM))) then
                call Driver_abort("[UniformMeshNewImpl] mesh dx differs from the Grid's dx")
            end if

            call tileDesc%getDataPtr(solnData, CENTER)
            do b = 1, mesh_f_nblocks(um_mesh(1))
                call mesh_f_limits(um_mesh(1), b, lim)
                call mesh_f_data_ptr(um_mesh(1), b, U)
                lo = lim(LOW, :) - off
                hi = lim(HIGH, :) - off
                if (toMesh) then
                    U(:, lim(LOW,1):lim(HIGH,1), lim(LOW,2):lim(HIGH,2), lim(LOW,3):lim(HIGH,3)) = &
                        solnData(:, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
                else
                    solnData(:, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = &
                        U(:, lim(LOW,1):lim(HIGH,1), lim(LOW,2):lim(HIGH,2), lim(LOW,3):lim(HIGH,3))
                end if
            end do
            call tileDesc%releaseDataPtr(solnData, CENTER)
            call itor%next()
        end do
        call Grid_releaseTileIterator(itor)

        if (ntiles /= 1) then
            call Driver_abort("[UniformMeshNewImpl] expected exactly one Grid block (UG on 1 MPI rank)")
        end if
    end subroutine um_copy

end module um_flashx
