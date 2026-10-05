#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

Mesh *mesh_create(int nx, int ny, int nz, int n_var)
{
    if (n_var < 1) return NULL;

    Mesh *m = calloc(1, sizeof *m);
    if (!m) return NULL;

    int size[3] = {nx, ny, nz};

    for (int a =0; a < 3; ++a){
        int used = (a < MESH_NDIM);
        m->n[a] = used ? size[a] : 1;
        m->n_h[a] = used ? MESH_NHALO : 0;
        m->n_t[a] = m->n[a] + 2*m->n_h[a];

        m->low[a] = 0.0;
        m->high[a] = 1.0;
        m->dx[a] =1.0/m->n[a];
    }

    m->n_var = n_var;

    for (int a=0; a < MESH_NDIM; ++a){
        if (m->n[a] < 1 || m->n[a] < m->n_h[a]){
            fprintf(stderr, "mesh_create: axis %d too small\n", a);
            free(m);
            return NULL;
        }
    }

    m->n_cell = (size_t)m->n_t[0] * (size_t)m->n_t[1] * (size_t)m->n_t[2];
    size_t n_bytes = (size_t)n_var * m->n_cell * sizeof(type_t);

    void *buf = NULL;

    if (MESH_ALLOC(buf, n_bytes) != 0){
        fprintf(stderr, "mesh_create: allocation failed\n");
        free(m);
        return NULL;
    }
    m->data = buf;

    // set all to 0
    for (int v=0; v < MESH_NVAR(m); ++v){
        MESH_LOOP_ALL(m, i, j, k)
            MESH_AT(m, v, i, j, k) = 0;
        MESH_LOOP_3D_END
    }

    return m;
}

void mesh_remove(Mesh *m)
{
    if (!m){
        return;
    }
    MESH_FREE(m->data);
    free(m);
}

void mesh_fill_halo(Mesh *m)
{   //for each axis fill the low and high halo strip/box (depending on halo size)
    // axis a fixed
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
            else {
                
            }
        }

        for (int side = 0; side < 2; ++side){
            // extend along a
            low[a] = (side == 0) ? -NH : N;
            high[a] = (side == 0) ? -1 : N+NH-1;
            int shift = (side==0)? N : -N; // distance from halo cell to the interior cell it needs to copy

            for (int v= 0; v < MESH_NVAR(m); ++v){
                MESH_LOOP_3D(low, high, i, j, k)
                    int position[3] = {i, j, k};
                    position[a] += shift;
                    MESH_AT(m, v, i, j, k) = MESH_AT(m, v, position[0], position[1], position[2]);
                MESH_LOOP_3D_END

            }
        }
        
    }
}

type_t mesh_get(const Mesh *m, int v, int i, int j, int k){
    return MESH_AT(m, v, i, j, k);
}

void   mesh_set(Mesh *m, int v, int i, int j, int k, type_t x){
    MESH_AT(m, v, i, j, k) = x;
}

void mesh_print(const Mesh *m, int v, int k)
{
    int n_x = MESH_N(m, 0); 
    int h_x = MESH_NH(m, 0);
    int n_y = MESH_N(m, 1);
    int h_y = MESH_NH(m, 1);

    printf("var %d, k = %d   (halo values in brackets)\n", v, k);
    for (int j = n_y + h_y - 1; j >= -h_y; --j) { /* top row first, so y points up */
        printf("j=%3d |", j);
        for (int i = -h_x; i < n_x + h_x; ++i) {
            int halo = (i < 0 || i >= n_x || j < 0 || j >= n_y);
            printf(halo ? " [%6.3g]" : "  %6.3g ", (double)MESH_AT(m, v, i, j, k));
        }
        printf("\n");
    }


}

int mesh_write(const Mesh *m, const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    // fixed order v, k, j, i incl. halo, independent of the layout
    for (int v = 0; v < MESH_NVAR(m); ++v)
        for (int k = -MESH_NH(m,2); k < MESH_N(m,2) + MESH_NH(m,2); ++k)
            for (int j = -MESH_NH(m,1); j < MESH_N(m,1) + MESH_NH(m,1); ++j)
                for (int i = -MESH_NH(m,0); i < MESH_N(m,0) + MESH_NH(m,0); ++i) {
                    type_t x = MESH_AT(m, v, i, j, k);
                    fwrite(&x, sizeof x, 1, f);
                }

    return fclose(f) == 0 ? 0 : -1;
}
