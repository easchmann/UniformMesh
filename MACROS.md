# Macros

This file lists the macros the glue code must provide and that the physics and mesh may use. 
A glue will be built form a set of `.ini` fragments (like in  `macro_example.ini`) that `ini_to_h.py` (WIP) turns into `mesh_config.h`.

A glue variant has to define every macro below and pass the property listed under **Check** for each one.

---

## Conventions

### Argument names
I use this in every macro with the same meaning:

| Arg | Meaning |
|---|---|
| `m` | mesh handle (`Mesh *`). Not hardcoded, so several blocks can exist at once |
| `v` | logical variable id, `0 <= v < n_var` |
| `i, j, k` | logical cell index per axis |
| `low, high` | index bounds |
| `p` | a pointer value |
| `n_bytes` | size in bytes (`size_t`) |

### Logical index space
`N(a) = MESH_N(m,a)`, `H(a) = MESH_NH(m,a)`. a is the axis, `0 <= a < 3`.
- Interior cells: `0 .. N(a)-1` on axis a.
- halo cells: `-H(a) .. -1` (low) and `N(a) .. N(a)+H(a)-1` (high).
- For an unused axis a >= MESH_NDIM: `N(a) = 1`, `H(a) = 0` whcih means the index is always 0.
- Physics should only use logical indices and not arithmetic/offsets

### Mesh representation
We don't want to fix the mesh representation, i.e. how the mesh stores its sizes or data. 
everything reads them thorugh the accessor macros below.
this allows the glue to swtch representations without touching the rest of the code

Representation:

| What | Accessor|
|---|---|
| interior cells per axis | MESH_N(m, a) |
| halo cells on each side per axis | MESH_NH(m, a) |
| total cells per axis = interior + 2*halo | MESH_NT(m, a) |
| number of variables | MESH_NVAR(m) |
| number of cells | MESH_NCELL(m) | 
| number of blocks | MESH_NBLOCK(m) |
| pointer to the storage | MESH_DATA(m) |

### Rules for macro body
1. Parenthesise every argument and the whole expression.
2. Function-like macros may evaluate an argument more than once -> don't pass expressions with side effects like i++.
3. Use _Pragma("") for pragmas instead of #pragma to avoid it being interpreted as a comment

---

## Mesh accessors

- These macros are the only way to get sizes and storage from a mesh handle.
- no side effects and return always the same result for the same m.


### MESH_N(m, a), MESH_NH(m, a), MESH_NT(m, a)
- **Kind:** expression (`int`)
- **Meaning:** number of interior cells (`MESH_N`), halo cells per side (`MESH_NH`)
  and total cells (`MESH_NT`) on axis `a`.
- **Implementation**:
  ```ini
  [MESH_N]
  args = m,a
  definition = ((m)->n[a])

  [MESH_NH]
  args = m,a
  definition = ((m)->n_h[a])

  [MESH_NT]
  args = m,a
  definition = ((m)->n_t[a])
  ```



### MESH_NVAR(m)
- **Kind:** expression (`int`)
- **Meaning:** number of variables stored in the mesh.
- **Implementation**
  ```ini
  [MESH_NVAR]
  args = m
  definition = ((m)->n_var)
  ```

### MESH_NCELL(m) 
- **Kind:** expression (`size_t`)
- **Meaning:** total number of cells, including halo cells.
- **Properties:** `== MESH_NT(m,0) * MESH_NT(m,1) * MESH_NT(m,2)`.
- **Open:** stored field or computed from `MESH_NT`? Storing would save two multiplications for each `MESH_IDX` in a structure of arrays
- **Implementation**:
  ```ini
  [MESH_NCELL]
  args = m
  definition = ((m)->n_cell)

  ; orif  computed
  [MESH_NCELL]
  args = m
  definition = ((size_t)MESH_NT(m,0) * (size_t)MESH_NT(m,1) * (size_t)MESH_NT(m,2))
  ```

### MESH_DATA(m)
- **Kind:** expression (`MESH_TYPE *`)
- **Meaning:** base pointer of the storage that `MESH_IDX` offsets into.
- **Implementation:** 
  ```ini
  [MESH_DATA]
  args = m
  definition = ((m)->data)
  ```

---

## Types

### MESH_TYPE 
- **Kind:** type
- **Meaning:** the floating-point storage type of all mesh variables. Used by physics in mesh.h (?)
- **Check:** `sizeof(MESH_TYPE)` should match the size implied by `MESH_TYPE_NAME`
- **Implementation:**
  ```ini
  [MESH_TYPE]
  definition = double
  ```
  ```ini
  [MESH_TYPE]
  definition = float
  ```

### MESH_TYPE_NAME
- **Kind:** constant (string literal)
- **Meaning:** language independent name for `MESH_TYPE`, written to headers so other variants and tools can read the data.
- **Check:** check if consistent with `sizeof(MESH_TYPE)`.
- **Implemention:**
  ```ini
  [MESH_TYPE_NAME]
  definition = "float64"
  ```


---

## Dimensionality and halo

### MESH_NDIM
- **Kind:** constant (integer)
- **Meaning:** number of used dimensions. All axes `0..MESH_NDIM-1` are used.
- **Properties:** Must be a plain integer literal so that `#if MESH_NDIM == 2` works.
- **Implementation (for now):**
  ```ini
  [MESH_NDIM]
  definition = 2
  ```

### MESH_IS2D, MESH_IS3D
- **Kind:** constant (0 or 1)
- **Meaning:** `1` if axis 1 (resp. axis 2) is used, else `0`. Physics multiplies stencil offsets by these, so one can work inmultiple dimensions:
  `MESH_AT(m, v, i, j + MESH_IS2D, k)`.
- **Implementation:**
  ```ini
  [MESH_IS2D]
  definition = (MESH_NDIM >= 2)
  ```
  ```ini
  [MESH_IS3D]
  definition = (MESH_NDIM >= 3)
  ```

### MESH_NHALO
- **Kind:** constant
- **Meaning:** halo width on each used axis. used during creation of mesh to set
  `MESH_NH(m,a)` (`MESH_NHALO` for active axes, 0 otherwise).
- **Open:** physics knows its stencil size. Should physics set a minimum (`MESH_NHALO_MIN`) and the glue pick `MESH_NHALO >= MESH_NHALO_MIN`? (My reasoning is that maybe setting a halo x*NHALO_MIN would need to fill the halo only every x steps?)
- **Implementation**
  ```ini
  [MESH_NHALO]
  definition = 1
  ```

---

## Layout

### MESH_CELL(m, i, j, k)
- **Kind:** expression (`size_t`)
- **Meaning:** linear position of cell `(i,j,k)` among all the `MESH_NCELL(m)` cells, (including halo cells). This is needed by `MESH_IDX`.
- **Properties:** projects the multi dimensional cells on the range `[0, MESH_NCELL(m))`.
  can use the previously defined accessor macros (e.g. `MESH_NH`, `MESH_NT`)
- **Implementation:** 
  ```ini
  [MESH_CELL]
  args = m, i, j, k
  definition = ( (size_t)((i) + MESH_NH(m,0))
    + (size_t)MESH_NT(m,0) * ( (size_t)((j) + MESH_NH(m,1))
    + (size_t)MESH_NT(m,1) *   (size_t)((k) + MESH_NH(m,2)) ) )
  ```

### MESH_IDX(m, v, i, j, k)
- **Kind:** expression (`size_t`)
- **Meaning:** offset into `MESH_DATA(m)` of variable `v` in cell `(i,j,k)`.
  **This macro is the layout.** Physics never calls it directly; it goes through
  `MESH_AT`.

- **Implementation:**
  ```ini
  ; structure of arrays
  [MESH_IDX]
  args = m,v,i,j,k
  definition = ((size_t)(v) * MESH_NCELL(m) + MESH_CELL(m,i,j,k))
  ```
  ```ini
  ; array of structures: all vars of cell 0, then cell 1 etc
  [MESH_IDX]
  args = m,v,i,j,k
  definition = (MESH_CELL(m,i,j,k) * (size_t)MESH_NVAR(m) + (size_t)(v))
  ```

---

## Access

### MESH_AT(m, v, i, j, k)  
- **Meaning:** the value of variable `v` in cell `(i,j,k)`.
- **IMplementation*:**
  ```ini
  [MESH_AT]
  args = m,v,i,j,k
  definition = (MESH_DATA(m)[MESH_IDX(m,v,i,j,k)])
  ```

---

## Loops

### MESH_LOOP_INTERIOR(m, i, j, k)
- **Kind:** statement-opening
- **Meaning:** iterate over all interior cells of `m`.
- **Implementation:**
  ```ini
  [MESH_LOOP_INTERIOR]
  args = m,i,j,k
  definition =
      for (int k = 0; k < MESH_N(m,2); ++k)
      for (int j = 0; j < MESH_N(m,1); ++j)
      for (int i = 0; i < MESH_N(m,0); ++i) {
  ```

### MESH_LOOP_END
- **Kind:** statement-closing
- **Meaning:** closes the most recent `MESH_LOOP_*` opener.
- **Implementation:**
  ```ini
  [MESH_LOOP_END]
  definition = }
  ```

### MESH_LOOP_ALL(m, i, j, k)
- **Meaning:** iterate over interior and halo cells.
- **Implementation**:
  [MESH_LOOP_ALL]
  args = m,i,j,k
  definition =
      for (int k = -MESH_NH(m,2); k < MESH_N(m,2)+MESH_NH(m,2); ++k)
      for (int j = -MESH_NH(m,1); j < MESH_N(m,2)+MESH_NH(m,1); ++j)
      for (int i = -MESH_NH(m,0); i < MESH_N(m,2)+MESH_NH(m,0); ++i) {
  ```

### MESH_LOOP_3D(low, high, i, j, k)
- **Meaning:** iterate over a box `low[a] <= idx <= high[a]`.
- **Open:** mesh-independent (no `m`), so the glue can't use mesh
  information here. Is that acceptable?
- **Implementation**:
  ```ini
  for (int k = (low)[2]; k <= (high)[2]; ++k)
  for (int j = (low)[1]; j <= (high)[1]; ++j)
  for (int i = (low)[0]; i <= (high)[0]; ++i) {
  ```

### MESH_LOOP_3D_END
- **Kind:** statement-closing
- **Meaning:** closes MESH_LOOP_3D opener.
- **Implementation:**
  ```ini
  [MESH_LOOP_END]
  definition = }}}
  ```

---

## Memory

### MESH_ALIGNMENT 
- **Kind:** constant (bytes, power of two)
- **Implementation:**
  ```ini
  [MESH_ALIGNMENT]
  definition = 64

### MESH_ALLOC(p, n_bytes) 
- **Kind:** expression (`int`, 0 on success)
- **Meaning:** allocate `n_bytes` aligned to `MESH_ALIGNMENT` and store the
  pointer in lvalue `p`.
- **Open:** return a status code or the pointer?
- **Implementation**: 
  ```ini
  args = p,n_bytes
  definition = posix_memalign((void **)&(p), MESH_ALIGNMENT, (n_bytes))
  ```

### MESH_FREE(p)
- **Kind:** statement
- **Implementation**:
  ```ini
  [MESH_FREE]
  args = p
  definition = free(p)

---

## Functions (not macros)

Functions in `mesh.c`, written in terms of the macros above.

| Function | Purpose |
|---|---|
| `mesh_create`, `mesh_remove` | create and remove a mesh |
| `mesh_fill_halo(m)`, `mesh_fill_halo_vars(m, v0, nv)` | periodic wraparound halo fill (exposed to physics) |
| `mesh_write`, `mesh_read` | I/O |

