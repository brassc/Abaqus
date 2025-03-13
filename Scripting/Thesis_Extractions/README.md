# README.md for `odb_extract.py`
`odb_extract.py` extracts the field data from a particular odb. 

## Implementation 
the python script is implemented by calling `abaqus python odb_extract.py odb_to_extract.odb output_data.pkl LE,U,S`
Specifying variables to extract is optional - if not specified, the script will extract all of them (this is slow). 

### Example Command
`abaqus python odb_extract.py C:\Users\cmb247\ABAQUS\K_DC_FALX\k-nofalx-001\helmet-brain-0.5.odb Output\nofalxnohemi0.5.pkl LE`

## Outputs
The output of the script is a directory, stored in `Outputs` directory in repo. This can then be worked with directly using python3 visualisations. 

```
output_directory/
├── model_structure.pkl         (contains node coordinates and element info)
├── manifest.pkl                (contains summary information)
├── Step-1/
│   ├── frame_0000.pkl
│   ├── frame_0001.pkl
│   └── ...
└── Step-2/
    ├── frame_0000.pkl
    └── ...
```
