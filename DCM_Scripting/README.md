# DCM Scripting

Scripts for applying spatially-varying predefined fields to spinal cord FE models, modelling compression sites in degenerative cervical myelopathy (DCM) simulations.

## Scripts

- **sets_script.py** - Creates node sets and predefined temperature fields along a compression axis via the Abaqus CAE API (Abaqus Python 2.7)
- **sets_inpmod_script.py** - Alternative that writes node sets and predefined fields directly into the `.inp` file (Abaqus Python 2.7). Use this when the CAE API corrupts the assembly tree (see below).
- **sets_scaled_inpmod_script.py** - Extends `sets_inpmod_script.py` with automatic per-site preload scaling based on sagittal cord diameter measurements. Use this when applying preloads across multiple indent sites or patients with varying cord sizes (Abaqus Python 2.7).
- **field_band_plot.py** - Visualises the raised cosine field distribution (Python 3, matplotlib)
- **coordinates.csv** - Compression site coordinates (center, upper, lower points per site)

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

peak_field_value = REFERENCE_PRELOAD × (compression_fraction_site / compression_fraction_ref)
```

Default reference calibration constants (the site at which preload was manually validated):

| Constant | Value | Description |
|---|---|---|
| `REFERENCE_PRELOAD` | 0.5 | Preload validated at the reference site |
| `REF_UPPER_CORD_SAG_DIST` | 6.82259 mm | Healthy cord AP diameter at reference site |
| `REF_INDENT_CORD_SAG_DIST` | 4.24591 mm | Compressed cord AP diameter at reference site |
| `compression_fraction_ref` | 0.3778 | (6.82259 − 4.24591) / 6.82259 |

These can be overridden before `execfile()` if the reference site changes.

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
COORDS_FILE = 'coordinates.csv'
execfile('sets_scaled_inpmod_script.py')

# Optionally override reference calibration before execfile()
REFERENCE_PRELOAD        = 0.5
REF_UPPER_CORD_SAG_DIST  = 6.82259
REF_INDENT_CORD_SAG_DIST = 4.24591
```

CSV format:
```
site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z,upper_cord_sag_dist,indent_cord_sag_dist
Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0,6.82259,4.24591
```

The columns `upper_cord_sag_dist` and `indent_cord_sag_dist` are optional. If absent, the script falls back to `PEAK_FIELD_VALUE` (default 0.15).

---

## Predefined Field Magnitude Profile Calculation

Field values follow a **raised cosine** (smooth step) profile along the compression axis:

```
value = (peak - min) * (1 + cos(pi * x)) / 2 + min
```

where `x = band_index / num_bands`. For 5 bands this gives 6 evenly spaced points from `x = 0` to `x = 1.0`, mapping to 0% to 120% of the physical compression region, set by input coordinates.  Only the first 5 points get node sets and predefined fields; the 6th is a virtual zero-crossing beyond the region boundary.

The raised cosine has zero gradient at centre point of indentation and at the edge ($x$ intercept) at 120% of the distance. This gives a smooth step akin to the Abaqus smooth step amplitude definition.

With `peak = 0.15`, `min = 0.0`:

| Band | Distance from centre | Field value |
|------|---------------------|-------------|
| 1    | 0%                  | 0.150       |
| 2    | 24%                 | 0.136       |
| 3    | 48%                 | 0.098       |
| 4    | 72%                 | 0.052       |
| 5    | 96%                 | 0.014       |
| (virtual) | 120%           | 0.000       |

## Output

For each band, the script creates:
- An assembly-level node set (e.g. `SITE1_BAND_1`)
- A predefined temperature field (e.g. `predefinedfield-1-fieldband1`)

Naming convention: `predefinedfield-<site_index>-fieldband<band_number>`
