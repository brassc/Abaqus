# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Python Version

**All Abaqus scripts must use Python 2.7** - this is the Python version bundled with Abaqus. When writing or modifying scripts that use Abaqus APIs (`odbAccess`, `abaqusConstants`, `abaqus`):
- Use `raw_input()` not `input()` for user prompts
- Use `.format()` string formatting, not f-strings
- Use `print()` function syntax (works in Python 2.7)
- Use `pickle.dump(..., protocol=2)` for pickle compatibility

## Project Overview

This is an Abaqus FEA post-processing codebase for biomechanical simulations (brain/cerebrum modeling). It contains Python scripts for:
1. **Node set generation** - Creating spatial filters (spherical regions) within Abaqus CAE
2. **Results analysis** - Extracting and plotting strain, displacement, and volume evolution data from ODB files

## Key Commands

### Running scripts in Abaqus CAE kernel
```python
# Set working directory first
new_dir='C:\\Users\\cmb247\\repos\\Abaqus\\Scripting'
os.chdir(new_dir)

# Execute script
execfile('main.py')
execfile('evol_calcs_main.py')
```

### Running Abaqus Python scripts without GUI
```bash
abaqus cae nogui=roi.py
```

### ODB data extraction
```bash
abaqus python odb_extract.py <input.odb> <output.pkl> <variables>
# Example:
abaqus python odb_extract.py C:\Users\cmb247\ABAQUS\K_DC_FALX\k-nofalx-001\helmet-brain-0.5.odb Output\nofalxnohemi0.5.pkl LE
```

## Architecture

### DCM_Scripting/
Spinal cord compression site modelling - see [DCM_Scripting/README.md](DCM_Scripting/README.md)

### Scripting/
- **main.py** - Orchestrator for sphere-based node set creation
- **node_set_export.py** - Export assembly node sets to CSV (`create_file_from_node_set`)
- **node_set_import.py** - Import CSV to create node sets (`create_node_set_from_file`)
- **nodes_within_sphere.py** - Filter nodes within spherical radius (`create_node_set_within_sphere`)
- **list_check.py** - Set difference operations between node lists
- **evol_calcs_main.py** - Extract EVOL field data from ODB files
- **evol_plots.py** - Plot cerebrum volume vs time and mesh size
- **roi.py** - Get max strain at ROI within sphere (r=3mm)

### Scripting/Thesis_Extractions/
- **odb_extract.py** - Frame-by-frame ODB data extraction to pickle files
- **odb_extract_plot.py** - Visualization and analysis of extracted data
- **odb_plot.py** - Direct ODB plotting

## Key Conventions

- Reference point for sphere center is named 'RP-SC' (or set manually via center_x, center_y, center_z)
- Scripts use `execfile()` pattern for Abaqus kernel execution
- Node lists are stored as comma-separated values in .txt files
- EVOL data files follow naming: `evol_<N>mm.csv` where N is mesh size
- Matplotlib styling uses `set_publication_style()` for consistent output (Times New Roman, 300dpi)

## Dependencies

**Abaqus Python 2.7 scripts** (run via `abaqus python` or `execfile()`):
- Abaqus Python API: `odbAccess`, `abaqusConstants`, `abaqus`
- `numpy` (included with Abaqus)

**Standalone Python 3 scripts** (e.g., `evol_plots.py`, `odb_extract_plot.py`):
- pandas, numpy, matplotlib
- pathlib for file handling
