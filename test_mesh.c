#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"


static int n_fail = 0;

/* report a failed condition with file:line, keep running */
#define CHECK(cond) do { if (!(cond)) {fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++n_fail;} } while (0)

/* Test mesh: non-square on purpose, so swapped axes show up as failures.
 * NX, NY, NZ are global cell counts; the block counts are set per run. */
enum { NX = 12, NY = 6, NZ = 4, NVAR = 3 };
static int BX = 1, BY = 1, BZ = 1;

static Mesh *create(void)
{
    return mesh_create(NX, NY, NZ, BX, BY, BZ, NVAR);
}

/* Encoded value from global indices: tells where a value came from. Exactly
 * representable in float for |i|,|j|,|k| < 50 and v < 8. */
static type_t code(int v, int i, int j, int k)
{
    return (type_t)(v * 1000000 + (i + 50) * 10000 + (j + 50) * 100 + (k + 50));
}

/* Periodic partner of index g on an axis with n interior cells. */
static int wrap(int g, int n)
{
    return ((g % n) + n) % n;
}

/* fill the interior of every block and variable with code() of the global index */
static void fill_with_code(Mesh *m)
{
    MESH_LOOP_BLOCKS(m, b)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_LOOP_INTERIOR(m, i, j, k)
                MESH_AT(m, b, v, i, j, k) = code(v, MESH_GIDX(m, b, 0, i), MESH_GIDX(m, b, 1, j), MESH_GIDX(m, b, 2, k));
            MESH_LOOP_3D_END
        }
    MESH_LOOP_END
}


/* 1. sizes after mesh_create follow the logical model */
static void test_sizes(void)
{
    printf("test_sizes\n");
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m) {
        return;
    }

    int size[3] = { NX, NY, NZ };
    int blocks[3] = { BX, BY, BZ };
    int n_block = 1;
    for (int a = 0; a < 3; ++a) {
        int used = (a < MESH_NDIM);
        CHECK(MESH_NB(m, a) == (used ? blocks[a] : 1));
        CHECK(MESH_N(m, a)  == (used ? size[a] / blocks[a] : 1));
        n_block *= MESH_NB(m, a);
        CHECK(MESH_NH(m, a) == (used ? MESH_NHALO : 0));
        CHECK(MESH_NT(m, a) == MESH_N(m, a) + 2 * MESH_NH(m, a));
    }
    CHECK(MESH_NBLOCK(m) == n_block);
    CHECK(MESH_NVAR(m) == NVAR);
    CHECK(MESH_NCELL(m) == (size_t)MESH_NT(m, 0) * MESH_NT(m, 1) * MESH_NT(m, 2));

    mesh_remove(m);
}

/* 2. MESH_IDX hits every storage slot exactly once */
static void test_bijection(void)
{
    printf("test_bijection\n");
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    size_t n_slot = (size_t)MESH_NBLOCK(m) * (size_t)MESH_NVAR(m) * MESH_NCELL(m);
    int *counter = calloc(n_slot, sizeof(*counter));
    if (!counter){
        mesh_remove(m);
        return;
    }

    MESH_LOOP_BLOCKS(m, b)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_LOOP_ALL(m, i, j, k)
                size_t idx = MESH_IDX(m, b, v, i, j, k);
                CHECK(idx < n_slot);
                if (idx < n_slot) ++counter[idx];
            MESH_LOOP_3D_END
        }
    MESH_LOOP_END

    for (size_t count = 0; count < n_slot; count++){
            CHECK(counter[count] == 1);
    }

    free(counter);
    mesh_remove(m);

}

/* 3. each loop macro visits each cell of its region exactly once */
static void test_loops(void)
{
    printf("test_loops\n");

    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    MESH_LOOP_BLOCKS(m, b)
    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, b, v, i, j, k) = 0;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_INTERIOR(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, b, v, i, j, k) += 1;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        int val = 0;
        if (0 <= i && i < MESH_N(m, 0) && 0 <= j && j < MESH_N(m, 1) && 0 <= k && k< MESH_N(m, 2)){
            val = 1;
        }
        for (int v = 0; v < MESH_NVAR(m); ++v){
            CHECK(MESH_AT(m, b, v, i, j, k) == val); 
        }
    MESH_LOOP_3D_END
    MESH_LOOP_END

    
    // part 2
    MESH_LOOP_BLOCKS(m, b)
    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, b, v, i, j, k) = 0;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, b, v, i, j, k) = 1;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            CHECK(MESH_AT(m, b, v, i, j, k) == 1); 
        }
    MESH_LOOP_3D_END
    MESH_LOOP_END

    mesh_remove(m);




}

/* 4. values written are read by mesh_get */
static void test_get_set(void)
{
    printf("test_get_set\n");

    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    fill_with_code(m);
    
    MESH_LOOP_BLOCKS(m, b)
    MESH_LOOP_INTERIOR(m, i, j, k)
        for (int v =0; v < MESH_NVAR(m); ++v){
            CHECK(mesh_get(m, b, v, i, j, k) == code(v, MESH_GIDX(m, b, 0, i), MESH_GIDX(m, b, 1, j), MESH_GIDX(m, b, 2, k)));
        }
    MESH_LOOP_3D_END
    MESH_LOOP_END

    mesh_remove(m);

}

/* 5. after mesh_fill_halo every cell equals the code of its periodic partner
 * in global index space, i.e. halos come from the neighbouring block */
static void test_halo(void)
{
    printf("test_halo\n");

    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    fill_with_code(m);
    mesh_fill_halo(m);

    MESH_LOOP_BLOCKS(m, b)
    for (int v =0; v < MESH_NVAR(m); ++v){
        MESH_LOOP_ALL(m, i, j, k)
            int p1 = wrap(MESH_GIDX(m, b, 0, i), MESH_N(m, 0) * MESH_NB(m, 0));
            int p2 = wrap(MESH_GIDX(m, b, 1, j), MESH_N(m, 1) * MESH_NB(m, 1));
            int p3 = wrap(MESH_GIDX(m, b, 2, k), MESH_N(m, 2) * MESH_NB(m, 2));
            CHECK(MESH_AT(m, b, v, i, j, k) == code(v, p1, p2, p3));
        MESH_LOOP_3D_END
    }
    MESH_LOOP_END

    mesh_print(m, MESH_NBLOCK(m) - 1, 0, 0);
    mesh_remove(m);


}

/* 6. block id <-> block coordinates is a bijection */
static void test_blocks(void)
{
    printf("test_blocks\n");
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    int *counter = calloc((size_t)MESH_NBLOCK(m), sizeof(*counter));
    if (!counter){
        mesh_remove(m);
        return;
    }

    for (int bk = 0; bk < MESH_NB(m, 2); ++bk)
        for (int bj = 0; bj < MESH_NB(m, 1); ++bj)
            for (int bi = 0; bi < MESH_NB(m, 0); ++bi){
                int b = MESH_BLOCK_ID(m, bi, bj, bk);
                CHECK(0 <= b && b < MESH_NBLOCK(m));
                if (b < 0 || b >= MESH_NBLOCK(m)) continue;
                ++counter[b];
                CHECK(MESH_BCOORD(m, b, 0) == bi);
                CHECK(MESH_BCOORD(m, b, 1) == bj);
                CHECK(MESH_BCOORD(m, b, 2) == bk);
            }

    for (int b = 0; b < MESH_NBLOCK(m); ++b){
        CHECK(counter[b] == 1);
    }

    free(counter);
    mesh_remove(m);
}

/* 7. cell centers: uniform spacing, continuous across blocks, inside the domain */
static void test_centers(void)
{
    printf("test_centers\n");
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    const double eps = 1e-12;
    MESH_LOOP_BLOCKS(m, b)
        for (int a = 0; a < 3; ++a){
            for (int idx = -MESH_NH(m, a); idx < MESH_N(m, a) + MESH_NH(m, a); ++idx){
                int g = MESH_BCOORD(m, b, a) * MESH_N(m, a) + idx;
                double x = MESH_X(m, b, a, idx);
                CHECK(fabs(x - (MESH_LOW(m, a) + (g + 0.5) * MESH_DX(m, a))) < eps);
                // neighbouring cells are dx apart, also across halo/interior
                if (idx > -MESH_NH(m, a)){
                    CHECK(fabs(x - MESH_X(m, b, a, idx - 1) - MESH_DX(m, a)) < eps);
                }
                if (0 <= idx && idx < MESH_N(m, a)){
                    CHECK(MESH_LOW(m, a) < x && x < m->high[a]);
                }
            }
            // first interior center of the next block is dx after this block's last one
            if (MESH_BCOORD(m, b, a) + 1 < MESH_NB(m, a)){
                int c[3] = { MESH_BCOORD(m, b, 0), MESH_BCOORD(m, b, 1), MESH_BCOORD(m, b, 2) };
                ++c[a];
                int nb = MESH_BLOCK_ID(m, c[0], c[1], c[2]);
                double gap = MESH_X(m, nb, a, 0) - MESH_X(m, b, a, MESH_N(m, a) - 1);
                CHECK(fabs(gap - MESH_DX(m, a)) < eps);
                // the high halo of b sits on the interior of its neighbour
                CHECK(fabs(MESH_X(m, b, a, MESH_N(m, a)) - MESH_X(m, nb, a, 0)) < eps);
            }
        }
    MESH_LOOP_END

    mesh_remove(m);
}

/* 8. invalid block counts are rejected */
static void test_invalid(void)
{
    printf("test_invalid\n");
    Mesh *m = mesh_create(NX, NY, NZ, 5, 1, 1, NVAR); // 12 % 5 != 0
    CHECK(m == NULL);
    mesh_remove(m);
    m = mesh_create(NX, NY, NZ, 0, 1, 1, NVAR);
    CHECK(m == NULL);
    mesh_remove(m);
}

/* 9. write a deterministic mesh to a file, for cmp across variants */
static void test_write(const char *path)
{
    printf("test_write -> %s\n", path);
    
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    fill_with_code(m);
    mesh_fill_halo(m);

    CHECK(mesh_write(m, path) == 0);

    mesh_remove(m);

}




int main(int argc, char **argv)
{
    const char *out = (argc > 1) ? argv[1] : "out.bin";

    printf("variant: ndim=%d halo=%d\n", MESH_NDIM, MESH_NHALO);

    /* single block and a tiled domain; out is written for the tiled one */
    int configs[2][3] = { {1, 1, 1}, {3, 3, 2} };
    for (int c = 0; c < 2; ++c){
        BX = configs[c][0]; BY = configs[c][1]; BZ = configs[c][2];
        printf("--- blocks %d x %d x %d\n", BX, BY, BZ);
        test_sizes();
        test_bijection();
        test_loops();
        test_get_set();
        test_halo();
        test_blocks();
        test_centers();
    }
    test_invalid();
    test_write(out);

    if (n_fail == 0) printf("ALL PASSED\n");
    else printf("%d CHECK(S) FAILED\n", n_fail);
    return n_fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
