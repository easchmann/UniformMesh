// C side of the Fortran bindings: plain C ABI functions around the mesh macros,
// so Fortran never sees the Mesh struct or the macros (see BINDINGS.md).
// Indices here are the C ones (0-based block, variable and local cell index);
// mesh_f.F90 converts to the 1-based Flash-X conventions.

#include <stdint.h>
#include "mesh.h"

_Static_assert(sizeof(type_t) == MESH_TYPE_BYTES, "MESH_TYPE_BYTES does not match MESH_TYPE");

int mesh_c_type_bytes(void)
{
    return (int)sizeof(type_t);
}

int mesh_c_ndim(void)
{
    return MESH_NDIM;
}

void mesh_c_sizes(const Mesh *m, int n[3], int nh[3], int nb[3], int *n_block, int *n_var)
{
    for (int a = 0; a < 3; ++a){
        n[a] = MESH_N(m, a);
        nh[a] = MESH_NH(m, a);
        nb[a] = MESH_NB(m, a);
    }
    *n_block = MESH_NBLOCK(m);
    *n_var = MESH_NVAR(m);
}

void mesh_c_domain(const Mesh *m, double low[3], double dx[3])
{
    for (int a = 0; a < 3; ++a){
        low[a] = MESH_LOW(m, a);
        dx[a] = MESH_DX(m, a);
    }
}

// global (0-based) index of the first and last interior cell of block b per axis
void mesh_c_block_limits(const Mesh *m, int b, int lo[3], int hi[3])
{
    for (int a = 0; a < 3; ++a){
        lo[a] = MESH_GIDX(m, b, a, 0);
        hi[a] = MESH_GIDX(m, b, a, MESH_N(m, a) - 1);
    }
}

// address of the first (lowest guard) cell of block b and the distance in elements
// between neighbours in v, i, j, k. Fortran builds its array pointer from this.
void mesh_c_block_layout(Mesh *m, int b, type_t **base, int64_t stride[4])
{
    int i0 = -MESH_NH(m, 0), j0 = -MESH_NH(m, 1), k0 = -MESH_NH(m, 2);
    int64_t o = (int64_t)MESH_IDX(m, b, 0, i0, j0, k0);

    *base = &MESH_AT(m, b, 0, i0, j0, k0);
    stride[0] = (int64_t)MESH_IDX(m, b, 1, i0, j0, k0) - o;
    stride[1] = (int64_t)MESH_IDX(m, b, 0, i0 + 1, j0, k0) - o;
    stride[2] = (int64_t)MESH_IDX(m, b, 0, i0, j0 + 1, k0) - o;
    stride[3] = (int64_t)MESH_IDX(m, b, 0, i0, j0, k0 + 1) - o;
}

// cell positions on axis a for all cells of block b incl. halo (MESH_NT(m,a) values)
// edge < 0: low edge, 0: center, > 0: high edge
void mesh_c_coords(const Mesh *m, int b, int a, int edge, double *x)
{
    int nh = MESH_NH(m, a);
    for (int idx = -nh; idx < MESH_N(m, a) + nh; ++idx){
        x[idx + nh] = (edge < 0) ? MESH_XL(m, b, a, idx)
                    : (edge > 0) ? MESH_XR(m, b, a, idx)
                    : MESH_X(m, b, a, idx);
    }
}
