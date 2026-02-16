# DCM Scripting

Scripts for applying spatially-varying predefined fields to spinal cord FE models, modelling compression sites in degenerative cervical myelopathy (DCM) simulations.

## Scripts

- **sets_script.py** - Creates node sets and predefined temperature fields along a compression axis via the Abaqus CAE API (Abaqus Python 2.7)
- **sets_inpmod_script.py** - Alternative that writes node sets and predefined fields directly into the `.inp` file (Abaqus Python 2.7). Use this when the CAE API corrupts the assembly tree (see below).
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

**Prerequisites:** You must have already written a `.inp` file from your model (e.g. via Job Manager or `mdb.jobs['Job-1'].writeInput()`). The script reads node coordinates from the CAE model but writes all output to the `.inp` file.

Run in Abaqus CAE kernel:

```python
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')

MODEL_NAME = 'Model-1'
INSTANCE_NAME = 'PART-1_1-1'
INP_FILE = 'Job-212.inp'

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
