#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"


static int n_fail = 0;

/* report a failed condition with file:line, keep running */
#define CHECK(cond) do { if (!(cond)) {fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++n_fail;} } while (0)

/* Test mesh: non-square on purpose, so swapped axes show up as failures.
 * NX, NY, NZ are global cell counts; the block counts are set per run.
 * Sizes grow with the halo so every block (up to 3 x 3 x 2) has >= MESH_NHALO cells. */
enum { NX = 3 * (MESH_NHALO + 2), NY = 3 * (MESH_NHALO + 1), NZ = 2 * (MESH_NHALO + 1), NVAR = 3 };
static int BX = 1, BY = 1, BZ = 1;

static Mesh *create(void)
{
    return mesh_create(NX, NY, NZ, BX, BY, BZ, NVAR);
}

/* Encoded value from global indices: tells where a value came from. Exactly
 * representable in float for |i|,|j|,|k| < 50 and v < 8 (holds up to MESH_NHALO = 10). */
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

/* 5b. mesh_fill_halo_vars fills exactly the requested variables, others keep their halo */
static void test_halo_vars(void)
{
    printf("test_halo_vars\n");

    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    // invalid ranges are rejected
    CHECK(mesh_fill_halo_vars(m, -1, 1) == -1);
    CHECK(mesh_fill_halo_vars(m, 0, -1) == -1);
    CHECK(mesh_fill_halo_vars(m, 0, NVAR + 1) == -1);
    CHECK(mesh_fill_halo_vars(m, NVAR - 1, 2) == -1);
    CHECK(mesh_fill_halo_vars(m, NVAR, 0) == 0); // empty range is fine

    const type_t sentinel = -1;
    const int v0 = 1, nv = 1; // middle variable only
    MESH_LOOP_BLOCKS(m, b)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_LOOP_ALL(m, i, j, k)
                MESH_AT(m, b, v, i, j, k) = sentinel;
            MESH_LOOP_3D_END
        }
    MESH_LOOP_END
    fill_with_code(m);
    CHECK(mesh_fill_halo_vars(m, v0, nv) == 0);

    MESH_LOOP_BLOCKS(m, b)
    for (int v = 0; v < MESH_NVAR(m); ++v){
        MESH_LOOP_ALL(m, i, j, k)
            int interior = 0 <= i && i < MESH_N(m, 0) && 0 <= j && j < MESH_N(m, 1) && 0 <= k && k < MESH_N(m, 2);
            int p1 = wrap(MESH_GIDX(m, b, 0, i), MESH_N(m, 0) * MESH_NB(m, 0));
            int p2 = wrap(MESH_GIDX(m, b, 1, j), MESH_N(m, 1) * MESH_NB(m, 1));
            int p3 = wrap(MESH_GIDX(m, b, 2, k), MESH_N(m, 2) * MESH_NB(m, 2));
            if (interior || (v0 <= v && v < v0 + nv)){
                CHECK(MESH_AT(m, b, v, i, j, k) == code(v, p1, p2, p3));
            }
            else {
                CHECK(MESH_AT(m, b, v, i, j, k) == sentinel);
            }
        MESH_LOOP_3D_END
    }
    MESH_LOOP_END

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
                // edges: dx apart, center in the middle, shared exactly with the next cell
                double xl = MESH_XL(m, b, a, idx);
                double xr = MESH_XR(m, b, a, idx);
                CHECK(fabs(xr - xl - MESH_DX(m, a)) < eps);
                CHECK(fabs(0.5 * (xl + xr) - x) < eps);
                CHECK(xl < x && x < xr);
                CHECK(xr == MESH_XL(m, b, a, idx + 1));
            }
            // domain boundaries coincide with the outer edges of the first/last block
            if (MESH_BCOORD(m, b, a) == 0){
                CHECK(MESH_XL(m, b, a, 0) == MESH_LOW(m, a));
            }
            if (MESH_BCOORD(m, b, a) == MESH_NB(m, a) - 1){
                CHECK(fabs(MESH_XR(m, b, a, MESH_N(m, a) - 1) - m->high[a]) < eps);
            }
            // first interior center of the next block is dx after this block's last one
            if (MESH_BCOORD(m, b, a) + 1 < MESH_NB(m, a)){
                int c[3] = { MESH_BCOORD(m, b, 0), MESH_BCOORD(m, b, 1), MESH_BCOORD(m, b, 2) };
                ++c[a];
                int nb = MESH_BLOCK_ID(m, c[0], c[1], c[2]);
                // blocks share their boundary edge exactly
                CHECK(MESH_XR(m, b, a, MESH_N(m, a) - 1) == MESH_XL(m, nb, a, 0));
                double gap = MESH_X(m, nb, a, 0) - MESH_X(m, b, a, MESH_N(m, a) - 1);
                CHECK(fabs(gap - MESH_DX(m, a)) < eps);
                // the high halo of b sits on the interior of its neighbour
                CHECK(fabs(MESH_X(m, b, a, MESH_N(m, a)) - MESH_X(m, nb, a, 0)) < eps);
            }
        }
    MESH_LOOP_END

    mesh_remove(m);
}

/* zero var 0 everywhere, add 1 on the box, check that exactly the box was hit */
static void check_box(Mesh *m, int b, const int low[3], const int high[3])
{
    MESH_LOOP_ALL(m, i, j, k)
        MESH_AT(m, b, 0, i, j, k) = 0;
    MESH_LOOP_3D_END

    MESH_LOOP_3D(low, high, i, j, k)
        MESH_AT(m, b, 0, i, j, k) += 1;
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        int inside = low[0] <= i && i <= high[0] && low[1] <= j && j <= high[1] && low[2] <= k && k <= high[2];
        CHECK(MESH_AT(m, b, 0, i, j, k) == (inside ? 1 : 0));
    MESH_LOOP_3D_END
}

/* 8. box macros: grown interior and face boxes, as used by telescoping RK
 * (stage s works on the interior grown by (MAXSTAGE-s)*NSTENCIL) */
static void test_boxes(void)
{
    printf("test_boxes\n");
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    int low[3], high[3];
    int b = MESH_NBLOCK(m) - 1;

    for (int g = 0; g <= MESH_NHALO; ++g){
        MESH_BOX_GROWN(m, low, high, g);
        for (int a = 0; a < 3; ++a){
            int ga = (a < MESH_NDIM) ? g : 0;
            CHECK(low[a] == -ga);
            CHECK(high[a] == MESH_N(m, a) - 1 + ga);
            // stays inside the storage for g <= halo
            CHECK(low[a] >= -MESH_NH(m, a));
            CHECK(high[a] <= MESH_N(m, a) + MESH_NH(m, a) - 1);
        }
        check_box(m, b, low, high);

        for (int d = 0; d < 3; ++d){
            // reconstruction box: one extra cell on both sides along d
            if (g < MESH_NHALO){
                MESH_BOX_GROWN(m, low, high, g);
                MESH_BOX_EXTEND(low, high, d, 1, 1);
                for (int a = 0; a < 3; ++a){
                    int ext = (a == d && a < MESH_NDIM) ? 1 : 0;
                    int ga = (a < MESH_NDIM) ? g : 0;
                    CHECK(low[a] == -ga - ext);
                    CHECK(high[a] == MESH_N(m, a) - 1 + ga + ext);
                }
                check_box(m, b, low, high);

                // faces normal to d: N + 2g + 1 faces along d (unused axis: unchanged)
                MESH_BOX_FACES(m, low, high, g, d);
                for (int a = 0; a < 3; ++a){
                    int ext = (a == d && a < MESH_NDIM) ? 1 : 0;
                    int ga = (a < MESH_NDIM) ? g : 0;
                    CHECK(low[a] == -ga);
                    CHECK(high[a] == MESH_N(m, a) - 1 + ga + ext);
                }
                check_box(m, b, low, high);
            }
        }
    }

    mesh_remove(m);
}

/* 8b. mesh_set_domain: dx and cell positions follow a non-unit domain, bad domains are rejected */
static void test_domain(void)
{
    printf("test_domain\n");
    Mesh *m = create();
    CHECK(m != NULL);
    if (!m){
        return;
    }

    const double low[3] = { -2.0, 1.0, 0.5 }, high[3] = { 8.0, 4.0, 2.5 };
    const double eps = 1e-12;
    CHECK(mesh_set_domain(m, low, high) == 0);
    MESH_LOOP_BLOCKS(m, b)
        for (int a = 0; a < 3; ++a){
            int n_global = MESH_N(m, a) * MESH_NB(m, a);
            CHECK(fabs(MESH_DX(m, a) - (high[a] - low[a]) / n_global) < eps);
            CHECK(MESH_LOW(m, a) == low[a]);
            if (MESH_BCOORD(m, b, a) == 0){
                CHECK(MESH_XL(m, b, a, 0) == low[a]);
            }
            if (MESH_BCOORD(m, b, a) == MESH_NB(m, a) - 1){
                CHECK(fabs(MESH_XR(m, b, a, MESH_N(m, a) - 1) - high[a]) < eps);
            }
        }
    MESH_LOOP_END

    // high <= low on a used axis: rejected, domain unchanged
    const double bad_high[3] = { -2.0, 4.0, 2.5 };
    CHECK(mesh_set_domain(m, low, bad_high) == -1);
    CHECK(fabs(MESH_DX(m, 0) - 10.0 / (MESH_N(m, 0) * MESH_NB(m, 0))) < eps);

    mesh_remove(m);
}

/* 9. invalid block counts and blocks smaller than the halo are rejected */
static void test_invalid(void)
{
    printf("test_invalid\n");
    Mesh *m = mesh_create(NX + 1, NY, NZ, 3, 1, 1, NVAR); // NX + 1 is not divisible by 3
    CHECK(m == NULL);
    mesh_remove(m);
    m = mesh_create(NX, NY, NZ, 0, 1, 1, NVAR);
    CHECK(m == NULL);
    mesh_remove(m);
    m = mesh_create(2 * (MESH_NHALO - 1), NY, NZ, 2, 1, 1, NVAR); // blocks of NHALO-1 cells
    CHECK(m == NULL);
    mesh_remove(m);
}

/* 10. write a deterministic mesh to a file, for cmp across variants */
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
        test_halo_vars();
        test_blocks();
        test_centers();
        test_boxes();
    }
    test_domain();
    test_invalid();
    test_write(out);

    if (n_fail == 0) printf("ALL PASSED\n");
    else printf("%d CHECK(S) FAILED\n", n_fail);
    return n_fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
