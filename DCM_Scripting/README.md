# DCM Scripting

Scripts for applying spatially-varying predefined fields to spinal cord FE models, modelling compression sites in degenerative cervical myelopathy (DCM) simulations.

## Scripts

- **sets_scaled_inpmod_gm_script.py** ⭐ **PRIMARY — USE THIS** ⭐ — Extends `sets_scaled_inpmod_overlap.py` with `CORD_SET_NAME` filtering for models where the cord shares a combined part instance with other anatomy. Requires a `Cord` assembly node set (GM + WM combined). Includes all features: preload scaling, overlap detection and resolution, summary file output (Abaqus Python 2.7).
- **sets_scaled_inpmod_script.py** - Extends `sets_inpmod_script.py` with automatic per-site preload scaling based on sagittal cord diameter measurements. Use only if the cord IS its own separate assembly instance (Abaqus Python 2.7).
- **sets_scaled_inpmod_overlap.py** - Extends `sets_scaled_inpmod_script.py` with automatic detection and resolution of overlapping node assignments across multiple compression sites, plus a structured summary output file. Does **not** have `CORD_SET_NAME` filtering — use `sets_scaled_inpmod_gm_script.py` instead (Abaqus Python 2.7).
- **sets_inpmod_script.py** - Writes node sets and predefined fields directly into the `.inp` file without preload scaling. Superceded by `sets_scaled_inpmod_script.py` (Abaqus Python 2.7).
- **sets_script.py** - Original version using the Abaqus CAE API. **Not recommended** — creates assembly-level sets via the API which can corrupt the assembly tree (see below). Retained for reference only (Abaqus Python 2.7).
- **field_band_plot.py** - Visualises the raised cosine field distribution (Python 3, matplotlib)
- **coordinates.csv** - Compression site coordinates (center, upper, lower points per site) — basic format without cord diameters
- **coordinates_scaled.csv** - Compression site coordinates with sagittal cord diameter measurements for preload scaling and overlap resolution

---

## Usage: sets_scaled_inpmod_gm_script.py ⭐ PRIMARY

Use this script when the spinal cord **shares a combined part instance** with other anatomy (bone, disc, ligaments) — i.e. the cord is NOT its own separate assembly instance.

**Prerequisites:**
1. Write a `.inp` from your model (Job Manager or `mdb.jobs['Job-1'].writeInput()`)
2. In CAE, create an assembly node set called `Cord` combining your GM and WM sets (e.g. `PART-1-1.P44;GM` + `PART-1-1.P45;WM`)
3. Check your `.inp` for the amplitude name (search `*Amplitude`) — the default assumed is `Amp-1-preload`. If your model uses a different name (e.g. `AMP-1`), the job will fail with `THERE IS NO AMPLITUDE BY THE NAME Amp-1-preload`. To fix: open the output `.inp` and search-replace `Amp-1-preload` with your actual amplitude name.

Run in Abaqus CAE kernel:

```python
import os
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')
```

**Option 1 — Single site (manual peak value, no scaling):**
```python
MODEL_NAME       = 'Model-1'
INSTANCE_NAME    = 'PART-1-1'
INP_FILE         = 'Job-212.inp'
PEAK_FIELD_VALUE = 0.3
CORD_SET_NAME    = 'Cord'
CENTER_POINT     = (x, y, z)
UPPER_POINT      = (x, y, z)
LOWER_POINT      = (x, y, z)
execfile('sets_scaled_inpmod_gm_script.py')
```

**Option 2 — Multiple sites from CSV (auto-scaled, overlap detection applied automatically):**
```python
MODEL_NAME    = 'Model-1'
INSTANCE_NAME = 'PART-1-1'
INP_FILE      = 'Job-212.inp'
CORD_SET_NAME = 'Cord'
COORDS_FILE   = 'coordinates_scaled.csv'
execfile('sets_scaled_inpmod_gm_script.py')
```

If `CORD_SET_NAME` is not set, all nodes in `INSTANCE_NAME` are classified.

**Point placement guidance:** upper, center, and lower points define the axis vector (`normalize(upper - lower)`). Place all three on the **same face and same mesh layer** of the cord surface. Mixing mesh layers on the same face introduces an artificial tilt into the axis vector, creating non-uniform temperature across the cord cross-section. Ensure `upper_cord_sag_dist` > `indent_cord_sag_dist` — swapping these produces negative field values.
fff
For overlap detection behaviour, output file naming, and summary file format, see the `sets_scaled_inpmod_overlap.py` section below — `sets_scaled_inpmod_gm_script.py` uses the same three-pass approach with identical output.

---

---

## Other scripts (secondary / legacy)

---

## Usage: sets_script.py (CAE API — not recommended)

Run in Abaqus CAE kernel:

```python
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')

MODEL_NAME = 'Model-1'
INSTANCE_NAME = 'PART-1_1-1'

# Option 1 - Single site
CENTER_POINT = (x, y, z)
UPPER_POINT = (x, y, z)
LOWER_POINT = (x, y, z)
execfile('sets_script.py')

# Option 2 - Multiple sites from CSV (takes priority if set)
COORDS_FILE = 'coordinates.csv'
execfile('sets_script.py')
```

## Usage: sets_inpmod_script.py (direct .inp modification — superceded)

This script exists because creating assembly-level sets via the Abaqus Python API (`assembly.Set`) can corrupt the assembly tree into sub-assemblies, causing the `.inp` writer to silently drop the sets and predefined fields. This script bypasses the CAE entirely by inserting `*Nset` and `*Temperature` keywords directly into an existing `.inp` file.

**Prerequisites:** You must have already written a `.inp` file from your model (e.g. via Job Manager or `mdb.jobs['Job-1'].writeInput()`). The script reads node coordinates from the CAE model but writes all output to a new `.inp` file (e.g. `Job-212.inp` -> `Job-212_modified.inp`). The original file is not modified.

Run in Abaqus CAE kernel:

```python
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')

MODEL_NAME = 'Model-1'
INSTANCE_NAME = 'PART-1_1-1'
INP_FILE = 'Job-212.inp'
PEAK_FIELD_VALUE = 0.15  # optional, defaults to 0.15

# Option 1 - Single site
CENTER_POINT = (x, y, z)
UPPER_POINT = (x, y, z)
LOWER_POINT = (x, y, z)
execfile('sets_inpmod_script.py')

# Option 2 - Multiple sites from CSV (takes priority if set)
COORDS_FILE = 'coordinates.csv'
execfile('sets_inpmod_script.py')
```

The script inserts:
- `*Nset` blocks before `*End Assembly`
- `*Temperature` blocks after `** PREDEFINED FIELDS` in the step section

It will skip writing if the sets already exist in the `.inp` file to avoid duplicates.

CSV format:
```
site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z
Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0
```

## Usage: sets_scaled_inpmod_script.py (direct .inp modification with preload scaling — use gm_script instead if cord shares an instance)

This script extends `sets_inpmod_script.py` with automatic per-site preload scaling. It is intended for use across multiple indent sites or patients where cord geometry varies. All `.inp` file handling is identical to `sets_inpmod_script.py`; the only difference is how the peak field value (preload) is determined.

### Preload scaling

The predefined temperature field represents swelling of the cord back towards its pre-compression size. Because the field drives thermal **strain** (fractional deformation) rather than displacement directly, the displacement produced is:

```
displacement = strain × current_cord_diameter
```

where `current_cord_diameter` is the cord's **compressed** diameter (`indent_cord_sag_dist`) — the size the cord elements are at when the preload is applied. To close the full compression gap, the required thermal strain is therefore:

```
required strain = (upper_cord_sag_dist - indent_cord_sag_dist) / indent_cord_sag_dist
```

where:
- `upper_cord_sag_dist` — sagittal (AP) cord diameter at the unaffected level (mm)
- `indent_cord_sag_dist` — sagittal (AP) cord diameter at the indent/compression site (mm)

This is the **compression ratio** (gap relative to compressed size), not the compression fraction (gap relative to uncompressed size). Using the uncompressed diameter as denominator systematically underestimates the required preload, with the error growing with compression severity.

The peak field value (preload) is scaled proportionally to this required strain, normalised to a known reference calibration point:

```
compression_ratio_site = (upper_cord_sag_dist - indent_cord_sag_dist) / indent_cord_sag_dist

compression_ratio_ref  = (REF_UPPER_CORD_SAG_DIST - REF_INDENT_CORD_SAG_DIST) / REF_INDENT_CORD_SAG_DIST

peak_field_value = PEAK_FIELD_VALUE × (compression_ratio_site / compression_ratio_ref)
```

`PEAK_FIELD_VALUE` is the desired preload at the reference site geometry — set this before `execfile()` to run different preload levels (e.g. 0.3, 0.4, 0.5). All other sites scale proportionally from it. It defaults to `REFERENCE_PRELOAD` (0.3) if not set.

Reference calibration constants (fixed — define the geometry at which preload was manually validated):

| Constant | Value | Description |
|---|---|---|
| `REFERENCE_PRELOAD` | 0.3 | Default value of `PEAK_FIELD_VALUE` if not set; preload validated at the reference site |
| `REF_UPPER_CORD_SAG_DIST` | 6.82259 mm | Healthy cord AP diameter at reference site |
| `REF_INDENT_CORD_SAG_DIST` | 4.24591 mm | Compressed cord AP diameter at reference site |
| `compression_ratio_ref` | 0.6068 | (6.82259 − 4.24591) / 4.24591 |

The calibration constants only need overriding if the reference site itself changes.

### Usage

Run in Abaqus CAE kernel:

```python
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')

MODEL_NAME    = 'Model-1'
INSTANCE_NAME = 'PART-1_1-1'
INP_FILE      = 'Job-212.inp'

# Option 1 - Single site (manual PEAK_FIELD_VALUE, no scaling)
PEAK_FIELD_VALUE = 0.5
CENTER_POINT = (x, y, z)
UPPER_POINT  = (x, y, z)
LOWER_POINT  = (x, y, z)
execfile('sets_scaled_inpmod_script.py')

# Option 2 - Multiple sites from CSV (scaling applied automatically per site)
PEAK_FIELD_VALUE = 0.3    # desired preload at reference site; all others scale from this (default)
COORDS_FILE = 'coordinates.csv'
execfile('sets_scaled_inpmod_script.py')
```

CSV format:
```
site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z,upper_cord_sag_dist,indent_cord_sag_dist
Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0,6.82259,4.24591
```

The columns `upper_cord_sag_dist` and `indent_cord_sag_dist` are optional. If absent, `PEAK_FIELD_VALUE` is applied directly with no scaling. `PEAK_FIELD_VALUE` defaults to `REFERENCE_PRELOAD` (0.3) if not set.

---

## Usage: sets_scaled_inpmod_overlap.py (multi-site with overlap detection)

This script extends `sets_scaled_inpmod_script.py` for multi-level DCM cases where compression sites at adjacent vertebral levels may have overlapping band regions. It uses a three-pass approach: classify all sites, resolve overlaps, then write the `.inp` file.

### When to use this script vs `sets_scaled_inpmod_gm_script.py`

| Scenario | Script to use |
|---|---|
| Single compression site | Either script (single-site path is identical in both) |
| Multiple non-overlapping sites | Either script |
| Multiple sites at adjacent levels (e.g. C4/5 and C5/6) | `sets_scaled_inpmod_gm_script.py` (primary) or this script if cord is its own instance |
| Anterior + posterior compression at the same level | Define as **one site** in the CSV (see below), either script |

### Overlap detection and resolution rule

When classifying nodes into bands, a node near the boundary between two adjacent compression levels may fall inside the band regions of **both** sites. Without resolution, both sites would write a predefined temperature field to that node in the `.inp` file — leading to conflicting or double-loaded boundary conditions.

`sets_scaled_inpmod_overlap.py` resolves this automatically using the **minimum-field-value rule**:

> For each contested node, compute the field value each competing site would assign it based on which band it falls in and that site's raised-cosine profile. Assign the node exclusively to the site giving it the **lower** field value. Remove it from all other sites' bands.

**Physical rationale:** The transition zone between two adjacent compression levels (e.g. the disc space between C4/5 and C5/6) should receive the *lower* of the two competing preloads — not the higher, and not a double contribution. The minimum-field-value rule is conservative: it avoids over-loading the boundary region while still ensuring the bands of the winning site extend continuously through the full cord cross-section at that level.

**Significant overlap warning:** If a contested node falls in an inner band (bands 1 or 2 — the highest field value region) of either competing site, the script prints a `*** WARNING ***` and flags the pair in the summary file. This indicates that the upper/lower extents defined in the CSV for those two sites are probably too broad and are encroaching on each other's high-field regions. The recommended action is to reduce the extent coordinates in the CSV.

### Same-level anterior + posterior compressions

If a patient has both anterior and posterior compression at the same vertebral level, **do not define these as two separate sites**. Nearest-centre or minimum-value splitting would split the cord cross-section artificially — anterior nodes going to one site, posterior nodes to another — which does not reflect the continuous loading through the cord.

Instead, define the combined anterior + posterior effect as a **single site** in the CSV, with the centre, upper, and lower points and cord diameter measurements representing the combined compression. This is both physically correct and simpler.

### Three-pass execution

1. **Pass 1 — Classify:** All sites are classified into bands without writing anything. Each node is tentatively assigned to one or more sites based on its axial projection relative to each site's centre.
2. **Pass 2 — Resolve:** Contested nodes (assigned to more than one site) are resolved using the minimum-field-value rule. `all_site_band_nodes` is updated in-place. The overlap table is printed.
3. **Pass 3 — Write:** The resolved node sets and predefined fields are written sequentially to the output `.inp` file.

### Output files

Two files are written to the same directory as the input `.inp` file:

- **Modified `.inp` file** — name is derived from the input filename with the following transformations:
  - `TEMPLATE` (case-insensitive) is stripped along with its adjacent separator characters
  - `PEAK_FIELD_VALUE` is appended as e.g. `0pt50` (decimal point replaced with `pt`)
  - Site names from the CSV (in order) are appended, separated by `_`
  - Example: `Job-N01-015-TEMPLATE_2STEP.inp` → `Job-N01-015_2STEP_0pt50_Site1_Site2.inp`

- **Summary `.txt` file** — same base name as the output `.inp`, with `_overlap_summary` suffix:
  - Run configuration (model, instance, file paths)
  - Overlap resolution table with plain-English explanation of what was adjusted and why
  - Per-site summary tables showing node set names, node counts, field values, and predefined field names for each band

### Usage

Run in Abaqus CAE kernel:

```python
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')

MODEL_NAME    = 'N01-015_2026-02-09-scripting'
INSTANCE_NAME = 'PART-1_1-1'
INP_FILE      = 'D:\\path\\to\\Job-N01-015-TEMPLATE_2STEP.inp'

# Option 1 - Single site (manual PEAK_FIELD_VALUE, no scaling or overlap detection)
PEAK_FIELD_VALUE = 0.5
CENTER_POINT = (x, y, z)
UPPER_POINT  = (x, y, z)
LOWER_POINT  = (x, y, z)
execfile('sets_scaled_inpmod_overlap.py')

# Option 2 - Multiple sites from CSV (scaling + overlap detection applied automatically)
PEAK_FIELD_VALUE = 0.3    # desired preload at reference site; all others scale from this (default)
COORDS_FILE = 'coordinates_scaled.csv'
execfile('sets_scaled_inpmod_overlap.py')
```

CSV format (`coordinates_scaled.csv`):
```
site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z,upper_cord_sag_dist,indent_cord_sag_dist
Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0,6.82259,4.24591
Site2,0.0,10.0,0.0,0.0,10.0,10.0,0.0,10.0,-10.0,7.10,5.50
```

The columns `upper_cord_sag_dist` and `indent_cord_sag_dist` are optional. If absent, `PEAK_FIELD_VALUE` is applied directly with no scaling for that site. `PEAK_FIELD_VALUE` defaults to `REFERENCE_PRELOAD` (0.3) if not set.

---

## Predefined Field Magnitude Profile Calculation

Field values follow a **raised cosine** (smooth step) profile along the compression axis:

```
value = (peak - min) * (1 + cos(pi * x)) / 2 + min
```

where `x = band_index / num_bands`. For 5 bands this gives 6 evenly spaced points from `x = 0` to `x = 1.0`, mapping to 0% to 120% of the physical compression region, set by input coordinates.  Only the first 5 points get node sets and predefined fields; the 6th is a virtual zero-crossing beyond the region boundary.

The raised cosine has zero gradient at centre point of indentation and at the edge ($x$ intercept) at 120% of the distance. This gives a smooth step akin to the Abaqus smooth step amplitude definition.

With `PEAK_FIELD_VALUE = 0.3` (default), `min = 0.0`:

| Band | Distance from centre | Field value |
|------|---------------------|-------------|
| 1    | 0%                  | 0.300       |
| 2    | 24%                 | 0.271       |
| 3    | 48%                 | 0.196       |
| 4    | 72%                 | 0.104       |
| 5    | 96%                 | 0.029       |
| (virtual) | 120%           | 0.000       |

## Output

For each band, the script creates:
- An assembly-level node set (e.g. `SITE1_BAND_1`)
- A predefined temperature field (e.g. `predefinedfield-1-fieldband1`)

Naming convention: `predefinedfield-<site_index>-fieldband<band_number>`

`sets_scaled_inpmod_gm_script.py` and `sets_scaled_inpmod_overlap.py` both write a `_overlap_summary.txt` file alongside the `.inp` output, documenting the run configuration, overlap resolution results, and per-site band summaries. See the [Output files](#output-files) section above for naming details.

---

## Cord Oscillation Modeling: oscillation.py

Models cord motion driven by CSF (cerebrospinal fluid) pulsation, as an alternative to the flexion/extension loading used elsewhere in this README. Rather than applying a moment to bend the spine, it moves everything *except* the cord back and forth along the spine's own axis defined by reference points at top and bottom of spine, leaving the cord itself stationary. 

It does this by taking a completed job's `.inp` (which normally has `Step-1` for compression-site preload, then `Step-2` for the flexion/extension moment) and producing a new `.inp` where `Step-2` is dropped entirely and replaced with a `Step-3` that oscillates the vertebral column instead. 

### How it works

Every patient's `.inp` has two reference points (RPs) coupled to the top and bottom of the modeled vertebral segment — one to the C2 (or whichever level sits at the top) surface, one to the C7 (or whichever level sits at the bottom) surface. The lower RP is the one held fixed for the preload step and is kinematically coupled to the whole vertebral/ligament surface, so driving the lower RP moves the bony/ligamentous anatomy while the cord, which isn't tied to that coupling, stays where it is.

The script finds these two RPs automatically by scanning the `.inp` for `*Coupling` lines and matching node set names that look like `c<number>-top` / `c<number>-base` (e.g. `c2-top`, `c3-base {typo}`, `c7-base`). The exact `*Nset` name for each RP (`m_Set-3`, `m_Set-4`, `m_Set-5`...) also varies patient to patient, so it's never hardcoded.

The oscillation direction is the 3D vector between these two RPs' coordinates (upper minus lower) computed patient-by-patient from `.inp`. This is necessary because not every patient's model is oriented identically.

**RP scripting note**: an RP's coordinates can live in two different places in the `.inp`, and getting this wrong produces a plausible-looking but wrong answer rather than an error (oscillation mainly in $x$ direction rather than majority $z$). Some patients have their RPs inside a small, dedicated reference-point instance (their own tiny `*Instance` block). However, most of cohort have them as bare nodes sitting directly in the `*Assembly` block, with no `instance=` on their `*Nset`. For that second case, the search has to specifically exclude every `*Node` block that belongs to a `*Part` definition, because a `*Part`'s own internal node numbering restarts at 1 too, and a naive whole-file search would return the main anatomy mesh's node 1 instead of the real upper RP.

### Oscillation profile and magnitude

The time-varying displacement (a roughly cardiac-cycle-shaped curve, `Amp-3-osc`) and its 0.76mm peak magnitude are copied from Sam Schaefer's work, which itself uses the oscillatory cervical cord motion profile described by Mikulis et al. [1,2]. The magnitude is decomposed into three components along the computed direction vector (`dx, dy, dz`) and written as three `*Boundary` lines under one shared `*Boundary, ..., amplitude=Amp-3-osc` block, so all three move together in time and the net motion stays purely along that one direction, rather than tracing out some other path.

### Usage

Pure text processing — no Abaqus API is used, so it runs equally well via `abaqus python`, `execfile()` in the CAE kernel, or a plain Python interpreter.

```python
# Single file:
INP_PATH = r'D:\path\to\Job-xxx.inp'
execfile('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting\\oscillation.py')

# Batch (INP_PATH = None): every id_map.csv row where loading_condition == 'Flexion'
INP_PATH = None
execfile('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting\\oscillation.py')
```

### Output naming

The source `.inp` is only ever opened read-only — it's never modified — and the result is written into a new, renumbered sibling job folder, never into an existing one. The job number's leading digit is replaced with `3` (so it's visually distinct from the original run), and `removed`/`osc` markers are added: `removed` sits right before the compression-site suffix (`_0pt30_Site1_Site2...`) if the filename has one, or at the end if it doesn't (e.g. most PostOp jobs, which have no compression-site fields); `osc` is always appended at the very end.

```
Job-020-N01-011-PreOp-BC0pt35wEVOL/
  Job-020-N01-011-PreOp-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp
→
Job-320-N01-011-PreOp-BC0pt35wEVOL_removed_osc/
  Job-320-N01-011-PreOp-BC0pt35wEVOL_removed_0pt30_Site1_Site2_Site3_Site4_osc.inp
```

Re-running the script overwrites its *own* previously-generated output (so a bug fix can regenerate everyone's files in one go), but it will never overwrite the source `.inp` — that check is unconditional, not something a flag can bypass.

### References

[1] Mikulis DJ, Wood ML, Zerdoner OAM, Poncelet BP. Oscillatory motion of the normal cervical spinal cord. Radiology. 1994 Jul;192(1):117–121. doi: 10.1148/radiology.192.1.8208922. PMID: 8208922.

[2] Schaefer SD, Davies BM, Newcombe VFJ, Sutcliffe MPF. Could spinal cord oscillation contribute to spinal cord injury in degenerative cervical myelopathy? Brain and Spine. 2023;3:101743. doi: 10.1016/j.bas.2023.101743. PMID: 37383476; PMCID: PMC10293319.

---

## Results Extraction and Analysis

Post-processing pipeline for extracting maximum principal strain (MPS) and volume (EVOL) from the `Cord` set across completed jobs, and comparing strain-exceedance patterns across patients, loading conditions (flexion/extension), and surgical state (pre-op/post-op).

### Scripts

- **Alex_results_extraction.py** (Abaqus kernel, Python 2.7) — extracts max-MPS-over-integration-points and EVOL for every element in the `Cord` set, for every frame of the step. `Cord` is a **node** set, not an element set, so the script derives the element set itself: every element whose full connectivity (all nodes) lies within the `Cord` node set. Writes two files next to the `.odb`:
  - `<basename>_mps.csv` — one row per element per frame (`frame_index, frame_value, element_label, mps, volume`)
  - `<basename>_topology.csv` — element connectivity (`element_label, nodes`, semicolon-separated node labels), used later for spatial clustering ("blob") analysis

  Also auto-updates `id_map.csv`: after writing, it finds the row whose `odb_path` matches the `ODB_PATH` you set and fills in `csv_path` — no manual copy-paste needed. If no matching row exists yet, it prints a warning instead of failing.

- **mps_common.py** (Python 3) — shared helpers used by every script below:
  - `volume_weighted_percentile(df, p)` — the MPS value below which fraction `p` of total *volume* (not element count) lies. Weighting by volume rather than counting elements equally matters when comparing regions/patients with different mesh densities.
  - `pct_volume_above(df, threshold)` — % of total volume with `mps >= threshold`.
  - `PLOT_STYLE` — shared matplotlib rcParams.

- **Alex_results_plotting.py** (Python 3) — single-patient time-history diagnostic. Point `CSV_PATH` at one patient's `_mps.csv`. Computes T95/T99 from the **last frame**, then tracks % volume above those thresholds across **every** frame, to check whether exceedance peaks at the final frame or earlier in the simulation. Prints and logs the result to `id_map.csv` (matched by `csv_path`) in three new columns: `peak_frame_t95`, `peak_frame_t99`, `last_frame_idx`. Run this once per patient/condition to build a peak-frame log across the cohort — it's what decides `FRAME_MODE` below.

- **Alex_results_multipatient_plot - compare_threshold.py** (Python 3) — the main cross-patient comparison script. Reads `id_map.csv`, loads every patient's `_mps.csv`, and produces five plots in one run (see [Output plots](#output-plots) below).

### `id_map.csv`

Gitignored — maps an anonymized participant number to the real patient ID and file paths, so no real ID ever appears in anything committed to the repo.

| Column | Meaning |
|---|---|
| `participant` | Anonymized number (e.g. `6` → labeled `P6` in all plots/summaries) |
| `id` | Real patient/hospital ID — never leaves this file |
| `mJOA` | Pre-operative mJOA score |
| `loading_condition` | `Flexion` or `Extension` |
| `State` | `PreOp` or `PostOp` |
| `odb_path` | Full path to the job's `.odb` |
| `csv_path` | Full path to the `_mps.csv` (auto-filled by `Alex_results_extraction.py`) |
| `peak_frame_t95`, `peak_frame_t99`, `last_frame_idx` | Logged by `Alex_results_plotting.py`'s diagnostic |

**Anonymization boundary:** `id_map.csv` and the raw per-frame `_mps.csv`/`_topology.csv` files (which live next to the `.odb`, typically on `D:\`, outside the repo) are never committed. Only the aggregated summary CSVs and plots produced by the multipatient script — keyed solely by `P{n}` — are safe to commit.

### Usage

1. Add a row to `id_map.csv` for the new patient/condition/state (`participant, id, mJOA, loading_condition, State, odb_path` — leave `csv_path` blank).
2. In the Abaqus kernel:
   ```python
   ODB_PATH = r'D:\path\to\Job-xxx.odb'
   execfile('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting\\Alex_results_extraction.py')
   ```
   `csv_path` gets filled in automatically once this finishes.
3. Run `Alex_results_plotting.py` on at least one representative patient's CSV to check whether strain exceedance peaks at the final frame or earlier. If it peaks earlier and relaxes by the end, use `FRAME_MODE = 'peak'` (each element's max-ever MPS across all frames); if it peaks at the final frame, `'last'` is adequate.
4. Set `FRAME_MODE` at the top of `Alex_results_multipatient_plot - compare_threshold.py` accordingly, then run it (`python "Alex_results_multipatient_plot - compare_threshold.py"`).

### Output plots

| File | Shows |
|---|---|
| `multipatient_mps_plot_compare_thresholds.pdf` | % cord volume above T90/T95/T99 (cohort-pooled volume-weighted percentiles) per patient, flexion vs extension (PreOp only) |
| `multipatient_mps_plot_compare_thresholds_manual_thresholds.pdf` | Same style, fixed MPS thresholds 0.05/0.10/0.15/0.20 instead of percentiles (PreOp only) |
| `multipatient_mps_plot_blob_distribution.pdf` | Cumulative % of total cord volume above each threshold, by spatial cluster ("blob") size — pooled across all patients per (threshold, condition). Distinguishes a few large contiguous high-strain regions from many small scattered ones |
| `multipatient_mps_plot_blob_distribution_perpatient_*.pdf` (×4) + `..._perpatient_grid.pdf` | Same blob analysis, faceted per patient (one curve per patient, normalized to their own total cord volume) instead of pooled — one plot per threshold, plus a combined 2×2 grid |
| `multipatient_mps_plot_compare_thresholds_manual_thresholds_prepost.pdf` | Manual-threshold comparison (0.10/0.15 shown) **including PostOp** — solid markers = PreOp, hollow = PostOp |

Every plot has a matching `_summary.csv` written alongside it.

### Methodology notes

- **Volume-weighted, not element-count-weighted**: percentiles and thresholds are always computed by volume, so a coarse-meshed region can't be outvoted by a fine-meshed one just because it has more elements.
- **Blob clustering**: elements exceeding a threshold are grouped into connected components using face-sharing adjacency between elements (approximated as **≥4 shared nodes** — exact for this mesh since it's all `C3D8` hex elements, but not a formal face check against each element's specific face-node groups). Each blob's volume is converted to an effective radius via `r = (3V / 4π)^(1/3)` (equivalent-sphere radius). Cumulative distributions are drawn as step functions (`drawstyle='steps-post'`), not diagonally-interpolated lines, since nothing actually accumulates between one blob's size and the next.
- **FRAME_MODE**: `'last'` uses each element's MPS at the final frame; `'peak'` uses its max-ever MPS across all frames. Strain can spike mid-simulation and relax by the end, in which case `'last'` would underestimate true exposure — this is why the single-patient diagnostic exists.
