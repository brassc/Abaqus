For DCM Scripting, see README.md in DCM scripting folder. 

# Scripting for node set generation (selection of nodes inside a spherical area)
`nodes*` or `list*` `file*` or `main.py` scripts.
This set of programs exports a node set from Abaqus as comma separated list in a .txt file. 
From this exported list, node sets can be compared and a new node set created. 

The `main.py` is run in the Abaqus kernel via the `execfile('<filename.py>')` command. Set work directory in Abaqus to repo location. 
The `node_set_export.py` `create_file_from_node_set` function extracts nodes from assembly set and exports to .txt as comma separated list


To create a new node set from a .txt comma separated list in Abaqus, use `create_node_set_from_file` function in `file node_set_import.py`
This requires a text file called with nodes as comma separated values, filepath = 'filename.txt', and a new_node_set_name as a string

To find the nodes that aren't in a sub set of a larger set and write them to file, use `list_check.py`

To find nodes that are within a sphere and create node set, use `create_node_set_within_sphere` function in `nodes_within_sphere.py`
To find nodes external to this sphere and create complementary node set, use `list_check` then `create_node_set_from_file` functions. 

Centre point of sphere is set by reference point called 'RP-SC' for reference point sphere center. Centre point may also be set manually
using center_x, center_y, center_z variables in `main.py`. 

# Results Processing 
## EVOL evaluation and plotting
How to get EVOL data out of Abaqus:
1. Opn relevant CAE file (denoted by stated mesh size)
2. Open matching output database (.odb)
3. Navigate to Abaqus repo in kernel:
    - `new_dir='C:\\Users\\cmb247\\repos\\Abaqus\\Scripting'`
    - `os.chdir(new_dir)`
4. In kernel, run `execfile('evol_calcs_main.py')`
    - Within the script, check that .odb file name and output `'w'` filename is correct
    - This outputs a `evol_<N>mm.csv` file
5. Push new file to repo

The program `evol_plots.py` uses `evol_<N>mm.csv` to plot EVOL vs. time for multiple mesh sizes and to plot peak cerebrum volume vs. mesh size. 
    - Plots are saved as `peak_cerebrum_vol_vs_mesh_size.png` and `cerebrum_vol_vs_time.png`. 
Other plots are created directly from the `phd-thesis-template-2.4` repo using `plot_scripts/mesh_sensitivity/mesh_sensitivity_script.py` and `plot_scripts/mesh_sensitivity/mesh_sensitivity.csv`. 


## Max strain evaluation
The `roi.py` script gets maximum strain output from different mesh size output databases (.odb) at a particular $(x, y, z)$ location with radius $r=3$ and records in .csv file, `max_strain.csv`. The `roi.py` script is ran using the following command inside terminal from the directory containing the `roi.py` script: `abaqus cae nogui=roi.py`. This allows use of proprietary abaqus libraries, but avoids having to open cae GUI. 
Other plots are created directly from the `phd-thesis-template-2.4` repo using `plot_scripts/mesh_sensitivity/mesh_sensitivity_script.py` and `plot_scripts/mesh_sensitivity/mesh_sensitivity.csv`.





