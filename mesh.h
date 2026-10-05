#ifndef MESH_H
#define MESH_H

#include <stddef.h>
#include "mesh_config.h"

typedef MESH_TYPE type_t; // type name physics uses

typedef struct {
    int n[3], n_h[3], n_t[3]; // per block (all blocks have the same size)
    int nb[3]; // blocks per axis
    int n_block;
    int n_var;
    size_t n_cell; // cells per block incl. halo
    double low[3], high[3], dx[3]; // global domain
    type_t *data;
} Mesh;

// declarations
// nx, ny, nz: global interior cells, bx, by, bz: blocks per axis (must divide n)
Mesh *mesh_create(int nx, int ny, int nz, int bx, int by, int bz, int n_var);
void mesh_remove(Mesh *m);
// global domain [low, high] per axis (default [0,1]), recomputes dx; returns -1 if high <= low on a used axis
int mesh_set_domain(Mesh *m, const double low[3], const double high[3]);

void mesh_fill_halo(Mesh *m);
// fill the halo of variables v0 .. v0+nv-1 only; returns -1 (and fills nothing) for an invalid range
int mesh_fill_halo_vars(Mesh *m, int v0, int nv);

// analysis/check purpose (e.g. print)
type_t mesh_get(const Mesh *m, int b, int v, int i, int j, int k);
void mesh_set(Mesh *m, int b, int v, int i, int j, int k, type_t x);

// I/O
void mesh_print(const Mesh *m, int b, int v, int k); //text
int mesh_write(const Mesh *m, const char *path); //write to a binary file

#endif
