# DCM Scripting

Scripts for applying spatially-varying predefined fields to spinal cord FE models, modelling compression sites in degenerative cervical myelopathy (DCM) simulations.

## Scripts

- **sets_script.py** - Creates node sets and predefined temperature fields along a compression axis via the Abaqus CAE API (Abaqus Python 2.7)
- **sets_inpmod_script.py** - Alternative that writes node sets and predefined fields directly into the `.inp` file (Abaqus Python 2.7). Use this when the CAE API corrupts the assembly tree (see below).
- **sets_scaled_inpmod_script.py** - Extends `sets_inpmod_script.py` with automatic per-site preload scaling based on sagittal cord diameter measurements. Use this when applying preloads across multiple indent sites or patients with varying cord sizes (Abaqus Python 2.7).
- **sets_scaled_inpmod_overlap.py** - Extends `sets_scaled_inpmod_script.py` with automatic detection and resolution of overlapping node assignments across multiple compression sites, plus a structured summary output file. Use this for multi-level DCM cases (Abaqus Python 2.7).
- **field_band_plot.py** - Visualises the raised cosine field distribution (Python 3, matplotlib)
- **coordinates.csv** - Compression site coordinates (center, upper, lower points per site) — basic format without cord diameters
- **coordinates_scaled.csv** - Compression site coordinates with sagittal cord diameter measurements for preload scaling and overlap resolution

## Usage: sets_script.py (CAE API)

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

## Usage: sets_inpmod_script.py (direct .inp modification)

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

## Usage: sets_scaled_inpmod_script.py (direct .inp modification with preload scaling)

This script extends `sets_inpmod_script.py` with automatic per-site preload scaling. It is intended for use across multiple indent sites or patients where cord geometry varies. All `.inp` file handling is identical to `sets_inpmod_script.py`; the only difference is how the peak field value (preload) is determined.

### Preload scaling

The predefined temperature field represents swelling of the cord back towards its pre-compression size. Because the field drives thermal **strain** (fractional deformation) rather than displacement directly, the displacement produced is:

```
displacement = strain × cord diameter
```

To produce a target displacement equal to the compression at a given site, the required strain is therefore:

```
required strain = (upper_cord_sag_dist - indent_cord_sag_dist) / upper_cord_sag_dist
```

where:
- `upper_cord_sag_dist` — sagittal (AP) cord diameter measured at the unaffected level at the upper inflection point (mm)
- `indent_cord_sag_dist` — sagittal (AP) cord diameter measured at the indent/compression site (mm)

The peak field value (preload) is scaled proportionally to this required strain, normalised to a known reference calibration point:

```
compression_fraction_site = (upper_cord_sag_dist - indent_cord_sag_dist) / upper_cord_sag_dist

compression_fraction_ref  = (REF_UPPER_CORD_SAG_DIST - REF_INDENT_CORD_SAG_DIST) / REF_UPPER_CORD_SAG_DIST

peak_field_value = PEAK_FIELD_VALUE × (compression_fraction_site / compression_fraction_ref)
```

`PEAK_FIELD_VALUE` is the desired preload at the reference site geometry — set this before `execfile()` to run different preload levels (e.g. 0.3, 0.4, 0.5). All other sites scale proportionally from it. It defaults to `REFERENCE_PRELOAD` (0.5) if not set.

Reference calibration constants (fixed — define the geometry at which preload was manually validated):

| Constant | Value | Description |
|---|---|---|
| `REFERENCE_PRELOAD` | 0.3 | Default value of `PEAK_FIELD_VALUE` if not set; preload validated at the reference site |
| `REF_UPPER_CORD_SAG_DIST` | 6.82259 mm | Healthy cord AP diameter at reference site |
| `REF_INDENT_CORD_SAG_DIST` | 4.24591 mm | Compressed cord AP diameter at reference site |
| `compression_fraction_ref` | 0.3778 | (6.82259 − 4.24591) / 6.82259 |

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

The columns `upper_cord_sag_dist` and `indent_cord_sag_dist` are optional. If absent, `PEAK_FIELD_VALUE` is applied directly with no scaling. `PEAK_FIELD_VALUE` defaults to `REFERENCE_PRELOAD` (0.5) if not set.

---

## Usage: sets_scaled_inpmod_overlap.py (multi-site with overlap detection)

This script extends `sets_scaled_inpmod_script.py` for multi-level DCM cases where compression sites at adjacent vertebral levels may have overlapping band regions. It uses a three-pass approach: classify all sites, resolve overlaps, then write the `.inp` file.

### When to use this script vs `sets_scaled_inpmod_script.py`

| Scenario | Script to use |
|---|---|
| Single compression site | Either script (single-site path is identical in both) |
| Multiple non-overlapping sites | Either script |
| Multiple sites at adjacent levels (e.g. C4/5 and C5/6) | `sets_scaled_inpmod_overlap.py` |
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

The columns `upper_cord_sag_dist` and `indent_cord_sag_dist` are optional. If absent, `PEAK_FIELD_VALUE` is applied directly with no scaling for that site. `PEAK_FIELD_VALUE` defaults to `REFERENCE_PRELOAD` (0.5) if not set.

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

`sets_scaled_inpmod_overlap.py` additionally writes a `_overlap_summary.txt` file alongside the `.inp` output, documenting the run configuration, overlap resolution results, and per-site band summaries. See the [Output files](#output-files) section above for naming details.
