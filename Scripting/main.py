from abaqus import *
from abaqusConstants import *
from nodes_within_sphere import create_node_set_within_sphere
from node_set_export import create_file_from_node_set
from node_set_import import create_node_set_from_file
from list_check import list_check
import sys

# Get the model database
modelDB = mdb.models['Model-1']


##Create half and half swelling
#export set node labels
create_file_from_node_set(set_name='GM_Node_Set')
#create_file_from_node_set(set_name='Swelling_Region_DC_Side')
print('file(s) created from node set(s)')

#compare the two lists
#print('comparing node lists...')
#list_check(list1='GM_Node_Set_NodeList.txt', list2='Swelling_Region_DC_Side_NodeList.txt', filename_list='node_numbers.txt')

#create new set from new list of the nodes that are not in both lists
#create_node_set_from_file(file_path='node_numbers.txt', new_node_set_name='Non_Swelling')



##Create swelling sphere and non-swelling rest
# Specify the existing node set and the radius of the sphere
existing_node_set_name = 'GM_Node_Set'
radius = 90 # in mm

#specify filenames and set names
nodes_outside_sphere_filename = 'nodes_outside_sphereBIG.txt'
nodes_outside_sphere_set_name='NodesOutsideSphereBIG'
nodes_inside_sphere_filename = 'nodes_inside_sphereBIG.txt'
nodes_inside_sphere_set_name='NodesInsideSphereBIG'

# Get the sphere centre reference point by name
reference_point_name = "RP-SCALT"

# Get the coordinates of sphere center reference point
center_x = int(mdb.models['Model-1'].rootAssembly.features[reference_point_name].xValue)
center_y = int(mdb.models['Model-1'].rootAssembly.features[reference_point_name].yValue)
center_z = int(mdb.models['Model-1'].rootAssembly.features[reference_point_name].zValue)
print(center_x, center_y, center_z)

# Create node set within sphere of specified radius
create_node_set_within_sphere(existing_node_set_name, 
                              nodes_inside_sphere_set_name,
                              radius,
                              nodes_inside_sphere_filename,
                              center_x,#=center_x,  # X-coordinate of the sphere's center
                              center_y,#=center_y,  # Y-coordinate of the sphere's center
                              center_z,#=center_z   # Z-coordinate of the sphere's center
                              ) 


#Create node set outside sphere
list_check(list1='GM_Node_Set_NodeList.txt', list2=nodes_inside_sphere_filename, filename_list=nodes_outside_sphere_filename)
create_node_set_from_file(nodes_outside_sphere_filename, nodes_outside_sphere_set_name)


