from abaqus import *
from abaqusConstants import *
from nodes_within_sphere import create_node_set_within_sphere
from node_set_export import create_file_from_node_set
from node_set_import import create_node_set_from_file
from list_check import list_check

# Get the model database
modelDB = mdb.models['Model-1']

#export set node labels
create_file_from_node_set(set_name='GM_Node_Set')
create_file_from_node_set(set_name='Swelling_Region_DC_Side')

#compare the two lists
list_check(list1='GM_Node_Set_NodeList.txt', list2='Swelling_Region_DC_Side_NodeList.txt', filename_list='node_numbers.txt')

#create new set from new list of the nodes that are not in both lists
create_node_set_from_file(file_path='node_numbers.txt', new_node_set_name='Non_Swelling')


#import function parameters
# Specify the path to the .txt file containing node numbers
#file_path = 'node_numbers.txt'
#new_node_set_name = 'Non_Swelling'
#create_node_set_from_file(file_path, new_node_set_name)

# Specify the existing node set and the radius of the sphere
existing_node_set_name = 'GM_Node_Set'
radius = 30 # in mm

#specify filenames and set names
nodes_outside_sphere_filename = 'nodes_outside_sphere.txt'
nodes_outside_sphere_set_name='NodesOutsideSphere'

# Create node set within sphere of specified radius
create_node_set_within_sphere(existing_node_set_name,
                              radius,
                              filename='nodes_within_sphere.txt',
                              center_x=50,  # X-coordinate of the sphere's center
                              center_y=50,  # Y-coordinate of the sphere's center
                              center_z=50   # Z-coordinate of the sphere's center
                              ) 


#list_check(list1='GM_Node_Set_NodeList.txt', list2='Swelling_Region_DC_Side_NodeList.txt', filename_list='node_numbers.txt')


#Create node set outside sphere
list_check(list1='GM_Node_Set_NodeList.txt', list2='nodes_within_sphere.txt', filename_list='nodes_outside_sphere.txt')
create_node_set_from_file(nodes_outside_sphere_filename, nodes_outside_sphere_set_name)


