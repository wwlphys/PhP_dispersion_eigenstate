# Phonon Polariton Solver

A MATLAB implementation of the **generalized Huang's equation** — a parameter-free,
microscopic theory of phonon polaritons (PhPs) in arbitrary polyatomic polar
crystals and finite-thickness van der Waals (vdW) slabs (bulk, monolayer, and
finite-thickness geometries).

This code accompanies the manuscript

> *Unified Parameter-Free Theory of Phonon Polaritons Beyond Diatomic Crystals*

All material parameters — second-order force constants Φ, Born effective charge
tensor **Z**, electronic dielectric tensor ε∞ — are obtained from first-principles
calculations (VASP + Phonopy); **no fitting** to polariton measurements is used.
The software solves the coupled lattice-vibration / electromagnetic-field
eigenproblem of the generalized Huang's equation.

---

## What it computes

For any wave vector `q = (qx, qy, qz)` and any crystal, the program yields:

| Output | Description |
|---|---|
| PhP dispersion ω(**q**) and isofrequency surfaces (IFSs) | 1D dispersion lines, 2D isofrequency contours, 3D isofrequency surfaces |
| Full eigenstates | Complete eigenvectors of ionic motion **and** electromagnetic field |
| Longitudinal/transverse character of each branch | Electric-field polarization of every branch |
| Phonon composition (projection weights) | Which phonon modes make up each PhP eigenstate |
| Out-of-plane circular polarization (phonon angular momentum) s_z | Chirality analysis in slab-confined geometries |

---

## Quick start (3-minute worked example)

On Linux with MATLAB:

```bash
# 1. get the code
cd public_code

# 2. in MATLAB, add the code folder (with all sub-folders) to the search path
addpath(genpath('/path/to/your/public_code'))
```

```matlab
% 3. run the bundled hBN bulk example (inputs shipped in examples/bulk_hBN)
cd('examples/bulk_hBN')
calculate_dispersion          % bulk dispersion -> dispersion.txt (~minutes)
```

```bash
# 4. (only for the eigenvector step) two mandatory preparations:
truncate -s -1 dispersion.txt          # strip the trailing newline (the reader requires it)
mkdir -p nime displacement             # output dirs; missing them aborts with "Invalid file identifier"
```

```matlab
calculate_eigenvector         % -> ovlp.txt, Edispersion.txt, dispeigvec.txt …
```

> 💡 The bundled `examples/bulk_hBN` inputs are self-contained — **no VASP/Phonopy
> needed** to verify the whole pipeline. Expected results: see [§6 Examples](#6-examples-no-vasp-required).

---

## Repository layout

```
public_code/
├── calculate_dispersion.m, main.m …       # bulk dispersion drivers + helpers
├── bulk_eigenvector/                       # bulk eigenvector / phonon-composition module
├── chiral/                                 # chirality / phonon angular momentum module
├── finite_thickness_based_on_unitcell/     # finite-thickness slab dispersion module
│   └── eigenvector/                        # slab eigenvector / field-plot sub-module
├── examples/
│   ├── bulk_hBN/                           # bulk-hBN input set (§6.1)
│   ├── vasp_hBN_441/                       # raw VASP/Phonopy outputs for hBN (§3.1.7 tutorial)
│   └── slab_hBN/                           # finite-thickness slab-hBN input set (§6.3)
└── read_* / *_preprocess / plot_* .m       # input readers, preprocessors, plotting helpers
```

---

## 1. System requirements

| Dependency | Role | Notes |
|---|---|---|
| **MATLAB** | runtime | R2020b or later recommended (`parpool`/`parfor`, `digits(64)`, standard linear algebra) |
| **Parallel Computing Toolbox** | multi-core root finding | optional for small runs, recommended for production |
| **Symbolic Math Toolbox** | high-precision arithmetic | `digits` used in `findeigenvector3.m` |
| **Phonopy ≥ 2.x** | upstream phonon inputs (`phonopy.yaml`, `qpoints.yaml`, `band.yaml`, `FORCE_CONSTANTS`) | only needed to *generate* inputs |
| **VASP ≥ 5.x** | DFPT inputs: Born effective charges, dielectric tensors | only needed to *generate* inputs |

- **Platform**: Linux x86-64 is primary (developed on Red-Hat-family cluster nodes).
  Plain MATLAB scripts should also run on macOS/Windows with two caveats: use
  forward-slash paths (already normalized), and reduce `ncpu` in `parameters.input`
  (e.g. `4`) on laptops.
- **Hardware**: none required. Multi-core paths scale with CPU count (typical runs
  use 8–28 cores; a single core works). Optional `gpuArray` paths exist behind
  comment toggles — **not** the default, GPU hardware is **not** needed.
- **No installation**: interpreted MATLAB scripts, no build step. The bundled
  `examples/` demo inputs are self-contained.

---

## 2. Installation & path setup

There is no build step — only three steps:

1. **Get the code.** Copy or `git clone` the `public_code` directory. From here
   on, `<repo>` denotes its absolute path on *your* machine.

2. **Add the code folder (with all sub-folders) to the MATLAB path** (in `startup.m`
   or at the prompt):
   ```matlab
   addpath(genpath('<repo>'));
   ```
   After this, all helper functions are reachable from any working directory.


**Typical installation time:** a few minutes (copy + path fix + verifying §6.1);
no compilation.

---

## 3. Input files

**Single most important convention** 🔑: all solvers read **all** inputs and write
**all** outputs as relative paths in the **current MATLAB working directory**.
`public_code` only supplies **functions** via the search path — never `cd` into
`bulk_eigenvector/`, `chiral/`, `eigenvector/`, … before running, or the run fails
with "file not found". **Every run must start from the directory that contains
your input files.**

### 3.1 Raw inputs you must generate with Phonopy / VASP

#### `OUTCAR.txt` (all calculation types)

Run a VASP **static** (DFPT) calculation of the unit cell to obtain the Born
effective charges and the dielectric tensor, then save the following two blocks
of the resulting `OUTCAR` as `OUTCAR.txt`:

```
MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)
 -------------------------------------
           5.447    -0.000     0.000
          -0.000     6.208     0.000
           0.000     0.000     4.248
 -------------------------------------

 BORN EFFECTIVE CHARGES (including local field effects)
 ------------------------------------------------------
 ion    1
     1     6.51576    -0.00000    -0.28531
 ...
```

⚠️ Delete any `PIEZOELECTRIC TENSOR` block before saving — the preprocessor
would mis-read it.

#### Bulk dispersion / isofrequency: `qpoints.yaml` + `phonopy.yaml`

Phonopy `band.conf` (with `WRITEDM` and `QPOINTS`):

```
ATOM_NAME = ***
DIM = ***
FORCE_SETS = READ
EIGENVECTORS = .TRUE.
WRITEDM = .TRUE.
QPOINTS = .TRUE.
```

`QPOINTS` file: **the first line is the number of q points**, followed by one
`qx qy qz` per line — fractional coordinates in Phonopy convention (the solver
multiplies them by the reciprocal-lattice vectors to obtain absolute q in m⁻¹).
Example: 201 q points along the second reciprocal-lattice direction (b*, the y
axis), from Γ (0,0,0) to (0,0.001,0):

```
201
0 0 0
0 5E-6 0
...
0 0.001 0
```

Then run:

```bash
phonopy -c POSCAR -p -s band.conf
```

- **1D dispersion**: use a 1D `QPOINTS` line → read `dispersion.txt` (iq, q [m⁻¹], freq [Hz]).
- **2D isofrequency**: use a 2D `QPOINTS` grid → plot with `plot2Disofreq.m`.
- **3D isofrequency surface**: see the parallel workflow in §4.1.

#### Bulk eigenvectors / phonon composition: `band.yaml`

Run Phonopy again on a **Γ-origin** `BAND` (only the first q point is used, so
only one wave-vector direction can be compared per run):

```
ATOM_NAME = ***
DIM = ***
FORCE_SETS = READ
EIGENVECTORS = .TRUE.
BAND = 0 0 0  0.5 0 0.0
```

```bash
phonopy -c POSCAR -p -s band.conf --nac
```

⚠️ **This run rewrites `phonopy.yaml`** (the `configuration:` block records the
`BAND` settings). If you need `phonopy.yaml` for the dispersion run afterwards,
regenerate it with the §3.1.2 `QPOINTS` command once more.

#### Slab dispersion: `qpoints.yaml` + `phonopy.yaml` + `FORCE_CONSTANTS`

Same `band.conf` as §3.1.2 but with `FORCE_CONSTANTS = WRITE` added, then run
with `--full-fc`:

```bash
phonopy -c POSCAR -p -s band.conf --full-fc
```

Also create `iq.txt` yourself — a single integer selecting the q-point slice
(e.g. `1`). ⚠️ `iq.txt` is **not shipped** with the examples; you must create it
before any slab run.

#### Chirality: `band.yaml`

Use a `FORCE_CONSTANTS = READ` + `BAND_POINTS` configuration (example in §6.2) and
run `phonopy -c POSCAR-unitcell -p -s band.conf` to obtain `band.yaml`.

#### From raw VASP/Phonopy outputs to solver inputs (worked tutorial)

`examples/vasp_hBN_441/` is a complete raw input set for the same hBN system
(`POSCAR-unitcell`, `SPOSCAR`, `phonopy_disp.yaml`, `FORCE_SETS`, `BORN`,
`OUTCAR.txt`) — the upstream data from which `examples/bulk_hBN/` was derived
(identical force constants and Z*/ε∞). To reproduce the bulk pipeline from scratch:

1. Create a `QPOINTS` file (§3.1.2), then with a `FORCE_SETS = READ` + `WRITEDM` +
   `QPOINTS` `band.conf` run `phonopy -c POSCAR-unitcell -p -s band.conf` →
   `qpoints.yaml` + `phonopy.yaml`;
2. Switch to a `BAND = 0 0 0  0.5 0 0.0` + `BAND_POINTS = 2` config and run
   `phonopy ... --nac` → `band.yaml` (keep it; re-run step 1 if you want the
   QPOINTS version of `phonopy.yaml` back);
3. `OUTCAR.txt` is already in the required format (§3.1.1); `BORN` is the
   Phonopy/DFPT equivalent and is *not* read by the solvers;
4. Create `parameters.input` (§3.3) and run `calculate_dispersion` from this
   directory.

### 3.2 Preprocessed files (generated automatically)

The drivers run the `*_preprocess.m` helpers automatically at the top of every
run; they strip alphabetical characters from the raw YAML text so `fscanf` can
read the numeric fields directly. The `clean_*.txt` intermediates need no manual
preparation:

| File | Written by | Read by |
|---|---|---|
| `clean_qpoints.txt` | `dm_preprocess` | `read_dm` |
| `clean_phonopy_yaml.txt` | `phonopy_yaml_preprocess` | `read_VnM` |
| `clean_OUTCAR.txt` | `OUTCAR_preprocess` | `read_Z` |
| `clean_band_yaml.txt` | `band_yaml_preprocess` | `read_displace` / `read_displaceq` |

The raw `qpoints.yaml` / `phonopy.yaml` / `OUTCAR.txt` / `band.yaml` must be
present in the current working directory. (Slab drivers have a `firsttime` flag
at the top; keep the default — it controls whether the preprocessors re-run.)

### 3.3 `parameters.input` formats

The solver reads its run parameters from `parameters.input`. **The exact field
order matters**; the bulk and slab formats differ.

**Bulk dispersion / IFS** (4 lines):

```
ncpu            % number of CPU cores (parpool size)
nw              % number of frequency steps (⚠️ plain decimal digits only — no scientific notation; make it as large as feasible, double until results stop changing)
minow           % lower frequency bound (Hz)
maxow           % upper frequency bound (Hz)
```

Example (first wide scan):

```
28
5000000
0
5e13
```

**Finite-thickness slab dispersion** (15 lines):

```
thick           % slab thickness d (m)
ncpu            % number of CPU cores
qxi             % imaginary in-plane q component (real-space decay)
minlamda        % lower out-of-plane decay length (qzi) bound
lamdastep       % decay-length step
maxlamda        % upper decay-length bound
minqzr          % lower real out-of-plane q component
qzrstep         % step
maxqzr          % upper real out-of-plane q component
nw              % number of frequency steps (plain digits only)
minw            % lower frequency bound (Hz)
maxw            % upper frequency bound (Hz)
outputunit      % 1:THz, 2:cm^-1, 3:meV
is_slab         % 0 or 1 (whether to use the finite-thickness dielectric ratio)
niter           % iterative refinement count
polarization    % 1:p, 2:s (q along a principal axis, _ps solver); 4:arbitrary q p+s (arb_q solver)
```

Example (q along a principal axis, p polarization):

```
4.9e-9      % thick (m)
28          % ncpu
0           % qxi
-1e6        % minlamda
1e5         % lamdastep
1e7         % maxlamda
-1e6        % minqzr
1e5         % qzrstep
1e7         % maxqzr
2000        % nw
3.5e13      % minw (Hz)
5.5e13      % maxw (Hz)
3           % outputunit: meV
1           % is_slab
4           % niter
1           % polarization: p
```

Frequency precision ≈ `1/nw/(ncpu+1)^niter`. Start with a wide λ/qzr window to
locate where solutions concentrate, then refine; too-large windows produce many
higher-order roots and slow the solver by an order of magnitude (see §5).

### 3.4 Auxiliary control files

| File | Purpose |
|---|---|
| `iq.txt` | a single integer selecting the q-point slice (slab runs; create it yourself) |
| `target_k.txt` | `kx ky kz` (3 floats) to focus the 3D-IFS root search near a given k point |
| `path.txt` | optional; if present, `read_dm` reads the q-point file from the path in its first line |
| `iq_done.txt` / `../all_iq_done.txt` | marker files (list of finished q-indices) for the random-q 3D-IFS workflow (§4.1) |

---

## 4. How to use (four calculation types)

> All commands below assume §2 is done (paths repointed, `addpath` run) and the
> input files sit in your current MATLAB working directory.

### 4.1 Bulk dispersion + isofrequency surfaces

```matlab
calculate_dispersion
```

`calculate_dispersion.m` runs the preprocessors, loads the data, then calls a
root-finding solver. Choose the solver by editing the solver line(s) at the
bottom of the driver:

- **1D dispersion / 2D IFS**: `findzeropointrealallq_insert` (sequential q order).
- **3D isofrequency surface**: `findzeropointRANDq_insert` (random q order, so many
  nodes can work simultaneously) — spread the inputs into several job
  sub-directories, `touch ../all_iq_done.txt` in the parent, run in each
  sub-directory, periodically refresh the marker file with
  `cat */iq_done.txt > ../all_iq_done.txt`, and finally concatenate
  `cat */dispersion2d.txt > dispersion2d.txt`.

Plot: `plot3Disofreq.m` (general anisotropic) or `plot3D_rotate_isofreq.m`
(in-plane isotropic).

### 4.2 Bulk eigenvectors / phonon composition

Two mandatory preparations:

```bash
truncate -s -1 dispersion.txt     # strip trailing newline (the reader loops `while ~feof`; a final \n makes the last fscanf return empty and abort)
mkdir -p nime displacement        # output dirs; missing them aborts with "Invalid file identifier"
```

Then run **from the input directory** (not from `bulk_eigenvector/`):

```matlab
calculate_eigenvector
```

Outputs: `ovlp.txt` (q, ω, phonon branch, polarization–phonon overlap projection
weight — only PhP branches kept), `Edispersion.txt`, `dispeigvec.txt` (per PhP
branch: 12 phonon displacement components + 3 electric-field components).
`findeigenvector3.m` / `plot_displacement_*.m` produce atomic-vibration
animations/plots.

### 4.3 Finite-thickness slab dispersion (unit-cell-based)

⚠️ Two solver families exist for slab runs — **choose by wave-vector direction,
not by filename recency**:

| Family | Entry point | Handles | When to use |
|---|---|---|---|
| `_ps` | `calculate_surface_dispersion_ps` | **p (1) or s (2) polarization, q along a principal axis** (q is projected onto x: `qx=sqrt(qx²+qy²)`); does not accept `polarization=4` | q along a crystal axis (the common case) |
| `arb_q` | `calculate_surface_dispersion_arb_q` | **full p+s, arbitrary q directions** (`polarization=4`, 6N+18 unknowns, matching SI Eq. (S66)) | q **not** along a principal axis |

```matlab
calculate_surface_dispersion_ps       % q along a principal axis: polarization 1 (p) or 2 (s)
% calculate_surface_dispersion_arb_q  % q NOT along a principal axis: polarization 4 (p+s)
% material-specific variants (same solver, tuned for a specific crystal):
%   calculate_surface_dispersion_hBN, _MoO3, _MoO3_w1
```

The run writes `surface_dispersionq*.txt` in the current directory, then:

```matlab
pick3D_minw_arb_q
cat surface_dispersion_picked_q* > dispersion.txt
```

extracts the lowest-`mfrfdiff` points per q slice. A slight staircase in the
curves is an unavoidable Phonopy `qpoints.yaml` output precision limit.

> ℹ️ Run `pick3D_minw_arb_q` **only after** the scan of a q slice has finished:
> on a partially written file it may emit empty outputs, because it discards
> boundary points and requires a full neighbourhood around each candidate.

Slab eigenvectors (run directly; the `finite_thickness_based_on_unitcell`
folder is on the MATLAB path via the §2 `genpath` setup):

```matlab
calculate_surface_eigenvector_ps     % -> LargeDet_minw_eigenvector_ps.m
```

Run from the input directory (`eigenvector/` only supplies helpers via the search
path). The atomic-vibration plot is triggered in `LargeDet_minw_eigenvector_ps.m`
around line ~357; the preceding `if` controls which q point is drawn.

### 4.4 Chirality / angular momentum

```matlab
calculate_chirality       % -> SpinAlong_*.txt (q/x/y/z projection axes)
```

The working directory needs `qpoints.yaml` and `band.yaml`; `chiral/` is only on
the search path — **do not `cd` into it**.

---

## 5. Recommended frequency-range workflow (time savers)

- Bulk: start with a wide `maxow ≈ 5×10¹³` Hz, inspect where solutions
  concentrate, then zoom in.
- `nw` should be large (≥ 5×10⁶); double it until results stop changing. The
  quick check in §6.1 uses `nw=5×10⁵` just to locate branches fast.
- Slabs: compute a small λ/qzr window first and enlarge ×3 per step — large
  windows introduce many higher-order roots and slow the solver by an order of
  magnitude.
- A slight staircase in slab curves is a Phonopy `qpoints.yaml` precision limit.

---

## 6. Examples (no VASP required)

### 6.1 Bulk hBN

`examples/bulk_hBN/` is a complete, verified bulk-dispersion input set for
hexagonal BN (`DIM = 4 4 1`): `qpoints.yaml`, `phonopy.yaml`, `OUTCAR.txt`,
`band.yaml` — all raw (unpreprocessed) inputs.

**a) Dispersion.** `parameters.input` must contain four lines (§3.3), then:

```matlab
cd('examples/bulk_hBN')
calculate_dispersion
```

Expected result (28 cores, `nw=5e5`, `0–5e13` Hz, ~20–40 min on a shared cluster
node; `nw=5e4` locates all branches in ~10 min as a quick check):
`dispersion.txt` (201 q points, hundreds to ~1300 branch points, 1.5–49.8 THz),
`neig.txt` (3–6 branches per q; the Γ point contributes the minimum of 3),
`fw.txt`. The Γ-point frequencies ~22.45 / 24.53 / 48.29 THz are the hBN PhP
branches and match the known Reststrahlen bands.

> ⚠️ **Do not launch MATLAB with `-nojvm`**: the solver opens a process-based
> `parpool`, which requires the JVM. Start with `matlab -nodisplay -nosplash`
> (no `-nojvm`).

**b) Eigenvectors / phonon composition.**

```bash
truncate -s -1 dispersion.txt
mkdir -p nime displacement
```

```matlab
calculate_eigenvector
```

Expected result (~3 min, 28 cores): `ovlp.txt` rows with `overlap ≈ 1` and
`minerror ≲1e-7` are genuine PhP eigenstates; at Γ each branch folds into
specific phonons (e.g. 22.45 THz over branches 7/8, 48.29 THz over 11/12). A few
rows show larger `minerror` (~0.18) and NaN E-field components — those are the
non-polariton branches. `displacement/*.jpg` and `nime/*.ascii` are written into
the pre-created directories.

### 6.2 Chirality demo (seconds to a few minutes)

From the `bulk_hBN` directory (which already contains `qpoints.yaml` and
`band.yaml`):

```matlab
calculate_chirality
% -> SpinAlong_q/x/y/z.txt   each row: distq  freq  s (phonon angular momentum)
```

The four output files are written in under a minute. Exact numeric values are
machine/data dependent.

### 6.3 Slab hBN (inputs shipped, runs on SLURM)

`examples/slab_hBN/` is a complete finite-thickness slab input set: 10 nm thick,
28 cores, p polarization (`polarization=1`), meV output (`outputunit=3`), a large
λ/qzr window, plus `ncompare.in` for `pick3D_minw_arb_q`. q is along a principal
axis, so it runs with the `_ps` solver.

⚠️ Create `iq.txt` yourself (`echo 1 > iq.txt`) before running, then submit:

```bash
sbatch -N 1 --ntasks=1 --cpus-per-task=28 -p batch <<'EOF'
#!/bin/bash
cd /path/to/examples/slab_hBN
export PATH=/path/to/matlab/bin:$PATH
matlab -nodisplay -nosplash -r "calculate_surface_dispersion_ps" > slab.log 2>&1
EOF
```

Expected first row (iqzr=1, ilam=1): `1  1.4566e+06  1  1  65  171.12  170.9178  -1500000000  -148000000  ~3e-5`
— a ~171 meV branch with small `mfrfdiff` (~3×10⁻⁵).

---

## 7. Output file reference

| File | Produced by | Columns |
|---|---|---|
| `dispersion.txt` | bulk dispersion §4.1 | `iq  dist[m⁻¹]  freq[Hz]` — one line per PhP branch point |
| `dispersion2d.txt` | bulk dispersion §4.1 | `qx qy qz freq[Hz]` — isofrequency lines; ends with the sentinel line `-1 -1 -1 -1` |
| `neig.txt` | bulk dispersion §4.1 | number of branches per q-point (whitespace-separated) |
| `fw.txt` | bulk dispersion §4.1 | `freq[Hz] Re f Im f` diagnostic of the last q-point |
| `ovlp.txt` | eigenvectors §4.2 | `\|q\| freq[Hz] phonon_branch overlap minerror` |
| `Edispersion.txt` | eigenvectors §4.2 | `\|q\| freq[Hz] Re/Im Ex,Ey,Ez E_par E_perp minerror` |
| `dispeigvec.txt` | eigenvectors §4.2 | header `nband`, `nq`; per branch: `iq qx qy qz freq[Hz]`, 12 phonon displacements (Re/Im), 3 E-field comps |
| `phase.txt` | eigenvectors §4.2 | atomic-displacement & E-field phases |
| `PPP_CrystalEnergy.txt` | eigenvectors §4.2 | polarization & time-dependent energy densities (frequency subset) |
| `Spinz_a/Spinx_a/Spiny_a.txt` | eigenvectors §4.2 | phonon angular momentum projections (filled only for converged eigenstates) |
| `SpinAlong_{q,x,y,z}.txt` | chirality §4.4 | `distq freq s [projection]` |
| `surface_dispersionq*.txt` | slab dispersion §4.3 | `iq diste iqzr ilam iw ω/2π qzr qzi mfrfdiff`; column 6 unit follows `outputunit` |
| `surface_dispersion2dq*.txt` | slab dispersion §4.3 | `qx qy qz freq`; the freq column is **always** in cm⁻¹ (×0.03e12), regardless of `outputunit` |
| `surface_dispersion_picked_q*.txt` | `pick3D_minw_arb_q` §4.3 | lowest-`mfrfdiff` points extracted per q slice |

---

## 8. FAQ

**Q: MATLAB reports "file not found"?**
A: You `cd`'d into `bulk_eigenvector/`, `chiral/`, etc. All solvers read/write
from the **current working directory**; sub-directories only supply functions
(the §3 convention).

**Q: "Invalid file identifier" error?**
A: You skipped the eigenvector preparations: create `nime/` and `displacement/`
directories first (§4.2).

**Q: "Java must be initialized … parpool" error?**
A: MATLAB was launched with `-nojvm`. Use `matlab -nodisplay -nosplash` instead.

**Q: `nw = 5e6` fails to parse?**
A: In `parameters.input`, `nw` accepts plain decimal digits only — no scientific
notation (§3.3).

**Q: Staircase in slab dispersion curves?**
A: An unavoidable Phonopy `qpoints.yaml` output precision limit — normal.


---

## 9. Known limitations

- Primarily developed/tested on Linux/MATLAB; cross-platform caveats in §1.
- Some `.m` files have non-UTF-8 comments; this does not affect execution.

---

## 10. License & citation

  - **License**: This code is released under the MIT License.                                                                                                                                  
    See the `LICENSE` file in the repository root for the full text.                                                                                                                           
  - **Citation**: if you use this code in published work, please cite the                                                                                                                      
    accompanying manuscript.                                                                                                                                                                   
  - **Contact**: open an issue in the repository, or contact the first                                                                                                                 
    author of the manuscript (W. Wang, Sun Yat-sen University).                                                                 