#ifndef MESH_H
#define MESH_H

#include <stddef.h>
#include "mesh_config.h"

typedef MESH_TYPE type_t; // type name physics uses

typedef struct {
    int n[3], n_h[3], n_t[3];
    int n_var;
    size_t n_cell;
    double low[3], high[3], dx[3];
    type_t *data;
} Mesh;

// declarations
Mesh *mesh_create(int nx, int ny, int nz, int n_var);
void mesh_remove(Mesh *m);

void mesh_fill_halo(Mesh *m);

// analysis/check purpose (e.g. print)
type_t mesh_get(const Mesh *m, int v, int i, int j, int k);
void mesh_set(Mesh *m, int v, int i, int j, int k, type_t x);

// I/O
void mesh_print(const Mesh *m, int v, int k); //text
int mesh_write(const Mesh *m, const char *path); //write to a binary file

#endif
