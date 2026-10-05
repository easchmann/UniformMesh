#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

Mesh *mesh_create(int nx, int ny, int nz, int bx, int by, int bz, int n_var)
{
    if (n_var < 1) return NULL;

    Mesh *m = calloc(1, sizeof *m);
    if (!m) return NULL;

    int size[3] = {nx, ny, nz};
    int blocks[3] = {bx, by, bz};

    for (int a =0; a < 3; ++a){
        int used = (a < MESH_NDIM);
        m->nb[a] = used ? blocks[a] : 1;
        if (!used) size[a] = 1;

        if (size[a] < 1 || m->nb[a] < 1 || size[a] % m->nb[a] != 0){
            fprintf(stderr, "mesh_create: axis %d: %d cells not divisible into %d blocks\n", a, size[a], m->nb[a]);
            free(m);
            return NULL;
        }

        m->n[a] = size[a] / m->nb[a];
        m->n_h[a] = used ? MESH_NHALO : 0;
        m->n_t[a] = m->n[a] + 2*m->n_h[a];

        m->low[a] = 0.0;
        m->high[a] = 1.0;
        m->dx[a] = (m->high[a] - m->low[a]) / size[a];
    }

    m->n_block = m->nb[0] * m->nb[1] * m->nb[2];
    m->n_var = n_var;

    for (int a=0; a < MESH_NDIM; ++a){
        if (m->n[a] < m->n_h[a]){
            fprintf(stderr, "mesh_create: axis %d: block smaller than halo\n", a);
            free(m);
            return NULL;
        }
    }

    m->n_cell = (size_t)m->n_t[0] * (size_t)m->n_t[1] * (size_t)m->n_t[2];
    size_t n_bytes = (size_t)m->n_block * (size_t)n_var * m->n_cell * sizeof(type_t);

    void *buf = NULL;

    if (MESH_ALLOC(buf, n_bytes) != 0){
        fprintf(stderr, "mesh_create: allocation failed\n");
        free(m);
        return NULL;
    }
    m->data = buf;

    // set all to 0
    MESH_LOOP_BLOCKS(m, b)
        for (int v=0; v < MESH_NVAR(m); ++v){
            MESH_LOOP_ALL(m, i, j, k)
                MESH_AT(m, b, v, i, j, k) = 0;
            MESH_LOOP_3D_END
        }
    MESH_LOOP_END

    return m;
}

int mesh_set_domain(Mesh *m, const double low[3], const double high[3])
{
    for (int a = 0; a < MESH_NDIM; ++a){
        if (!(high[a] > low[a])){
            return -1;
        }
    }
    for (int a = 0; a < 3; ++a){
        m->low[a] = low[a];
        m->high[a] = high[a];
        m->dx[a] = (high[a] - low[a]) / (MESH_N(m, a) * MESH_NB(m, a));
    }
    return 0;
}

void mesh_remove(Mesh *m)
{
    if (!m){
        return;
    }
    MESH_FREE(m->data);
    free(m);
}

// id of the block next to b along axis a (dir = -1 low, +1 high), periodic
static int block_neighbor(const Mesh *m, int b, int a, int dir)
{
    int c[3];
    for (int x = 0; x < 3; ++x){
        c[x] = MESH_BCOORD(m, b, x);
    }
    c[a] = (c[a] + dir + MESH_NB(m, a)) % MESH_NB(m, a);
    return MESH_BLOCK_ID(m, c[0], c[1], c[2]);
}

void mesh_fill_halo(Mesh *m)
{
    mesh_fill_halo_vars(m, 0, MESH_NVAR(m));
}

int mesh_fill_halo_vars(Mesh *m, int v0, int nv)
{
    if (v0 < 0 || nv < 0 || v0 + nv > MESH_NVAR(m)){
        return -1;
    }

    //for each axis fill the low and high halo strip/box (depending on halo size)
    // axis a fixed. All blocks finish axis a before axis a+1, so corners are
    // copied from halos the neighbour already filled
    for (int a = 0; a < MESH_NDIM; ++a){
        int N = MESH_N(m,a);
        int NH = MESH_NH(m,a);
        int low[3];
        int high[3];

        // extend the strip/box along the other axes
        for (int b = 0; b < 3; ++b){
            // interior low and high bounds
            low[b] = 0;
            high[b] = MESH_N(m,b) - 1;
            // extend bounds if necessary
            if (b < a){ // axis b filled in earlier step, include its halo
                low[b] -= MESH_NH(m, b);
                high[b] += MESH_NH(m,b);
            }
        }

        MESH_LOOP_BLOCKS(m, blk)
            for (int side = 0; side < 2; ++side){
                // extend along a
                low[a] = (side == 0) ? -NH : N;
                high[a] = (side == 0) ? -1 : N+NH-1;
                int shift = (side==0)? N : -N; // distance from halo cell to the interior cell of the neighbour it copies
                int src = block_neighbor(m, blk, a, (side == 0) ? -1 : 1);

                for (int v= v0; v < v0 + nv; ++v){
                    MESH_LOOP_3D(low, high, i, j, k)
                        int position[3] = {i, j, k};
                        position[a] += shift;
                        MESH_AT(m, blk, v, i, j, k) = MESH_AT(m, src, v, position[0], position[1], position[2]);
                    MESH_LOOP_3D_END
                }
            }
        MESH_LOOP_END
    }
    return 0;
}

type_t mesh_get(const Mesh *m, int b, int v, int i, int j, int k){
    return MESH_AT(m, b, v, i, j, k);
}

void   mesh_set(Mesh *m, int b, int v, int i, int j, int k, type_t x){
    MESH_AT(m, b, v, i, j, k) = x;
}

void mesh_print(const Mesh *m, int b, int v, int k)
{
    int n_x = MESH_N(m, 0);
    int h_x = MESH_NH(m, 0);
    int n_y = MESH_N(m, 1);
    int h_y = MESH_NH(m, 1);

    printf("block %d (%d,%d,%d), var %d, k = %d   (halo values in brackets)\n",
           b, MESH_BCOORD(m, b, 0), MESH_BCOORD(m, b, 1), MESH_BCOORD(m, b, 2), v, k);
    for (int j = n_y + h_y - 1; j >= -h_y; --j) { /* top row first, so y points up */
        printf("y=%6.3f |", MESH_X(m, b, 1, j));
        for (int i = -h_x; i < n_x + h_x; ++i) {
            int halo = (i < 0 || i >= n_x || j < 0 || j >= n_y);
            printf(halo ? " [%6.3g]" : "  %6.3g ", (double)MESH_AT(m, b, v, i, j, k));
        }
        printf("\n");
    }
    printf("x=      |");
    for (int i = -h_x; i < n_x + h_x; ++i) {
        printf("  %6.3f ", MESH_X(m, b, 0, i));
    }
    printf("\n");
}

int mesh_write(const Mesh *m, const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    // fixed order b, v, k, j, i incl. halo, independent of the layout
    for (int b = 0; b < MESH_NBLOCK(m); ++b)
        for (int v = 0; v < MESH_NVAR(m); ++v)
            for (int k = -MESH_NH(m,2); k < MESH_N(m,2) + MESH_NH(m,2); ++k)
                for (int j = -MESH_NH(m,1); j < MESH_N(m,1) + MESH_NH(m,1); ++j)
                    for (int i = -MESH_NH(m,0); i < MESH_N(m,0) + MESH_NH(m,0); ++i) {
                        type_t x = MESH_AT(m, b, v, i, j, k);
                        fwrite(&x, sizeof x, 1, f);
                    }

    return fclose(f) == 0 ? 0 : -1;
}
