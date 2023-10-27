This set of programs exports a node set from Abaqus as comma separated list in a .txt file. 
From this exported list, node sets can be compared and a new node set created. 

The main.py is run in the Abaqus kernel via the execfile('<filename.py>') command.
The node_set_export.py create_file_from_node_set function extracts nodes from assembly set and exports to .txt as comma separated list


To create a new node set from a .txt comma separated list in Abaqus, use create_node_set_from_file function in file node_set_import.py 
This requires a text file called with nodes as comma separated values, filepath = 'filename.txt', and a new_node_set_name as a string

To find the nodes that aren't in a sub set of a larger set and write them to file, use list_check.py

To find nodes that are within a sphere and create node set, use create_node_set_within_sphere function in nodes_within_sphere.py
To find nodes external to this sphere and create complementary node set, use list_check then create_node_set_from_file functions. 

Centre point of sphere is set by reference point called 'RP-SC' for reference point sphere center. Centre point may also be set manually
using center_x, center_y, center_z variables in 'main.py'. 

