!> Glue between Flash-X and the UniformMesh (see flashx/README.md in the UniformMesh repo).
!!
!! Flash-X keeps owning the solution (unk of the UG Grid), which initialization, IO and
!! the timestep use. Hydro copies the interior into the mesh, lets Spark work on the mesh
!! blocks and copies the result back.
!!
!! Requirements, checked in um_init / um_copy:
!!   UG Grid on one MPI rank, periodic boundaries on all used axes, NDIM and NGUARD equal
!!   to MESH_NDIM and MESH_NHALO of mesh_config.h, AoS layout (Flash-X without -index-reorder).

#include "Simulation.h"
#include "constants.h"

module um_flashx
    use mesh_f
    implicit none
    private

    type(mesh_t), save, public :: um_mesh
    logical, save :: um_ready = .false.
    integer, save :: um_gsize(MDIM)   ! global interior cells per axis

    public :: um_init, um_copyIn, um_copyOut

contains

    ! create the mesh on the first call; aborts if the Flash-X setup does not fit the mesh
    subroutine um_init()
        use Grid_interface, ONLY: Grid_getGlobalIndexLimits
        use RuntimeParameters_interface, ONLY: RuntimeParameters_get
        use Driver_interface, ONLY: Driver_abort

        character(len=16), parameter :: bcNames(2*MDIM) = [character(len=16) :: &
            "xl_boundary_type", "xr_boundary_type", "yl_boundary_type", &
            "yr_boundary_type", "zl_boundary_type", "zr_boundary_type"]
        character(len=MAX_STRING_LENGTH) :: bc
        integer :: nb(MDIM), ierr, i

        if (um_ready) return

        if (mesh_f_ndim() /= NDIM) then
            call Driver_abort("[UniformMesh] mesh_config.h was generated for another NDIM")
        end if
        do i = 1, 2*NDIM
            call RuntimeParameters_get(trim(bcNames(i)), bc)
            if (trim(bc) /= "periodic") then
                call Driver_abort("[UniformMesh] only periodic boundaries are supported: " // trim(bcNames(i)))
            end if
        end do

        call Grid_getGlobalIndexLimits(um_gsize)
        call RuntimeParameters_get("um_nblockx", nb(IAXIS))
        call RuntimeParameters_get("um_nblocky", nb(JAXIS))
        call RuntimeParameters_get("um_nblockz", nb(KAXIS))

        call mesh_f_create(um_mesh, um_gsize(IAXIS), um_gsize(JAXIS), um_gsize(KAXIS), &
                           nb(IAXIS), nb(JAXIS), nb(KAXIS), NUNK_VARS, ierr)
        if (ierr /= MESH_OK) then
            call Driver_abort("[UniformMesh] mesh_f_create failed: um_nblock[xyz] must divide the " // &
                              "global cell counts and every block needs >= NGUARD cells per axis")
        end if
        if (mesh_f_nhalo(um_mesh) /= NGUARD) then
            call Driver_abort("[UniformMesh] MESH_NHALO /= NGUARD, regenerate mesh_config.h (install.sh)")
        end if
        if (mesh_f_layout(um_mesh) /= MESH_LAYOUT_AOS) then
            call Driver_abort("[UniformMesh] needs the AoS layout, like Flash-X's unk(var,i,j,k)")
        end if

        write(*, '(a,3i6,a,3i4,a,i3)') " [UniformMesh] global cells", um_gsize, &
            ", blocks", nb, ", halo", mesh_f_nhalo(um_mesh)
        um_ready = .true.
    end subroutine um_init

    ! Flash-X Grid -> mesh, interior cells of all variables
    subroutine um_copyIn()
        call um_copy(.true.)
    end subroutine um_copyIn

    ! mesh -> Flash-X Grid, interior cells of all variables
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
                call Driver_abort("[UniformMesh] the Grid block is not the whole domain (run on 1 MPI rank)")
            end if

            call tileDesc%getDataPtr(solnData, CENTER)
            do b = 1, mesh_f_nblocks(um_mesh)
                call mesh_f_limits(um_mesh, b, lim)
                call mesh_f_data_ptr(um_mesh, b, U)
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
            call Driver_abort("[UniformMesh] expected exactly one Grid block (UG on 1 MPI rank)")
        end if
    end subroutine um_copy

end module um_flashx
