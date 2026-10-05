#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"


static int n_fail = 0;

/* report a failed condition with file:line, keep running */
#define CHECK(cond) do { if (!(cond)) {fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++n_fail;} } while (0)

/* Test mesh: non-square on purpose, so swapped axes show up as failures. */
enum { NX = 7, NY = 5, NZ = 4, NVAR = 3 };

/* Encoded value: tells where a value came from. Exactly representable in float
 * for |i|,|j|,|k| < 50 and v < 8. */
static type_t code(int v, int i, int j, int k)
{
    return (type_t)(v * 1000000 + (i + 50) * 10000 + (j + 50) * 100 + (k + 50));
}

/* Periodic partner of index g on an axis with n interior cells. */
static int wrap(int g, int n)
{
    return ((g % n) + n) % n;
}

/* fill the interior of every variable with code() */
static void fill_with_code(Mesh *m)
{
    for (int v = 0; v < MESH_NVAR(m); ++v){
        MESH_LOOP_INTERIOR(m, i, j, k)
            MESH_AT(m, v, i, j, k) = code(v, i, j, k);
        MESH_LOOP_3D_END
    }
}


/* 1. sizes after mesh_create follow the logical model */
static void test_sizes(void)
{
    printf("test_sizes\n");
    Mesh *m = mesh_create(NX, NY, NZ, NVAR);
    CHECK(m != NULL);
    if (!m) {
        return;
    }

    int size[3] = { NX, NY, NZ };
    for (int a = 0; a < 3; ++a) {
        int used = (a < MESH_NDIM);
        CHECK(MESH_N(m, a)  == (used ? size[a] : 1));
        CHECK(MESH_NH(m, a) == (used ? MESH_NHALO : 0));
        CHECK(MESH_NT(m, a) == MESH_N(m, a) + 2 * MESH_NH(m, a));
    }
    CHECK(MESH_NVAR(m) == NVAR);
    CHECK(MESH_NCELL(m) == (size_t)MESH_NT(m, 0) * MESH_NT(m, 1) * MESH_NT(m, 2));

    mesh_remove(m);
}

/* 2. MESH_IDX hits every storage slot exactly once */
static void test_bijection(void)
{
    printf("test_bijection\n");
    Mesh *m = mesh_create(NX,NY, NZ, NVAR);
    CHECK(m != NULL);
    if (!m){
        return;
    }

    int *counter = calloc((size_t)(MESH_NVAR(m) * MESH_NCELL(m)), sizeof(*counter));
    if (!counter){
        mesh_remove(m);
        return;
    }

    for (int v = 0; v < MESH_NVAR(m); ++v){
        MESH_LOOP_ALL(m, i, j, k)
            int idx = MESH_IDX(m, v, i, j, k);
            CHECK((size_t)idx <  (size_t)MESH_NVAR(m)*(size_t)MESH_NCELL(m));
            ++counter[idx];
        MESH_LOOP_3D_END

    }

    for (size_t count = 0; count < (MESH_NVAR(m) * MESH_NCELL(m)); count++){
            CHECK(counter[count] == 1);
    }

    free(counter);
    mesh_remove(m);

}

/* 3. each loop macro visits each cell of its region exactly once */
static void test_loops(void)
{
    printf("test_loops\n");

    Mesh *m = mesh_create(NX,NY, NZ, NVAR);
    CHECK(m != NULL);
    if (!m){
        return;
    }

    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, v, i, j, k) = 0;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_INTERIOR(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, v, i, j, k) += 1;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        int val = 0;
        if (0 <= i && i < MESH_N(m, 0) && 0 <= j && j < MESH_N(m, 1) && 0 <= k && k< MESH_N(m, 2)){
            val = 1;
        }
        for (int v = 0; v < MESH_NVAR(m); ++v){
            CHECK(MESH_AT(m, v, i, j, k) == val); 
        }
    MESH_LOOP_3D_END

    
    // part 2
    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, v, i, j, k) = 0;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            MESH_AT(m, v, i, j, k) = 1;
        }
    MESH_LOOP_3D_END

    MESH_LOOP_ALL(m, i, j, k)
        for (int v = 0; v < MESH_NVAR(m); ++v){
            CHECK(MESH_AT(m, v, i, j, k) == 1); 
        }
    MESH_LOOP_3D_END

    mesh_remove(m);




}

/* 4. values written are read by mesh_get */
static void test_get_set(void)
{
    printf("test_get_set\n");

    Mesh *m = mesh_create(NX,NY, NZ, NVAR);
    CHECK(m != NULL);
    if (!m){
        return;
    }

    fill_with_code(m);
    
    MESH_LOOP_INTERIOR(m, i, j, k)
        for (int v =0; v < MESH_NVAR(m); ++v){
            CHECK(mesh_get(m, v, i, j, k) == code(v, i, j, k));
        }
    MESH_LOOP_3D_END

    mesh_remove(m);

}

/* 5. after mesh_fill_halo every cell equals the code of its periodic partner */
static void test_halo(void)
{
    printf("test_halo\n");

    Mesh *m = mesh_create(NX,NY, NZ, NVAR);
    CHECK(m != NULL);
    if (!m){
        return;
    }

    fill_with_code(m);
    mesh_fill_halo(m);

    for (int v =0; v < MESH_NVAR(m); ++v){
        MESH_LOOP_ALL(m, i, j, k)
            int p1 = wrap(i, MESH_N(m, 0));
            int p2 = wrap(j, MESH_N(m, 1));
            int p3 = wrap(k, MESH_N(m, 2));
            CHECK(MESH_AT(m, v, i, j, k) == code(v, p1, p2, p3));
        MESH_LOOP_3D_END
    }

    mesh_print(m, 0, 0);
    mesh_remove(m);


}

/* 6. write a deterministic mesh to a file, for cmp across variants */
static void test_write(const char *path)
{
    printf("test_write -> %s\n", path);
    
    Mesh *m = mesh_create(NX, NY, NZ, NVAR);
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

    test_sizes();
    test_bijection();
    test_loops();
    test_get_set();
    test_halo();
    test_write(out);

    if (n_fail == 0) printf("ALL PASSED\n");
    else printf("%d CHECK(S) FAILED\n", n_fail);
    return n_fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
