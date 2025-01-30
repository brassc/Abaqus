from odbAccess import *
from abaqusConstants import *
import math
import numpy as np
from abaqus import *
import sys
import os


# Force terminal output
sys.stdout = sys.__stdout__

def euclidean_distance_3d(x1, y1, z1, x2, y2, z2):
    """
    Calculate Euclidean distance between two 3D points using individual coordinates.
    
    Parameters:
        x1, y1, z1: Coordinates of first point
        x2, y2, z2: Coordinates of second point
        
    Returns:
        float: Euclidean distance between the points
    """
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

def get_max_strain(odb, step_name, frame_number, center_point, radius, mesh_size):
    """
    Get the maximum strain value within a given radius of a centerpoint.
    
    Parameters:
        centerpoint: Tuple of (x, y, z) coordinates of the centerpoint
        odb: Path to the .odb file
        radius: Radius around the centerpoint to search for maximum strain
        frame_number: Frame number (int) to extract strain data
        step_name: Step number to extract strain data
        
    Returns:
        float: Maximum strain value within the radius of the centerpoint
    """
    # Open the database
    print('opening database...')
    odb = openOdb(path=odb, readOnly=True)
    print('opened odb')

    
    try: 
        # Get the step and frame
        print('get step')
        step = odb.steps[step_name]
        print('get frame')
        frame = step.frames[frame_number]
        print('got frame')
        
        # Get the field output
        print('get field output')
        field_output = frame.fieldOutputs['LE']
        print('got field output')
        
        # Get the nodes within the radius
        print('get nodes')
        nodes = odb.rootAssembly.instances['E_GMVC-1'].nodes
        print('got nodes')
        print(type(nodes))
        print(nodes[0])
        
        
        nodes_within_radius = []
        for node in nodes:
            distance = euclidean_distance_3d(center_point[0], center_point[1], center_point[2], node.coordinates[0], node.coordinates[1], node.coordinates[2])
            if distance <= radius:
                nodes_within_radius.append(node.label)    
        print('got nodes within radius')
        print('Nodes within radius:', nodes_within_radius)
        print(type(nodes_within_radius))

        mesh=odb.rootAssembly.instances['E_GMVC-1'].elements
        print('got mesh')
        node_to_elements = {}
        for element in mesh:
            for node in element.connectivity:
                if node in nodes_within_radius:
                    if node not in node_to_elements:
                        node_to_elements[node]=[element.label]
                    else:
                        node_to_elements[node].append(element.label)

        print('got node to elements')
        #print(node_to_elements)

        # First, get unique element labels from the dictionary
        element_list = []
        for node_elements in node_to_elements.values():
            for elem in node_elements:
               if elem not in element_list:
                    element_list.append(elem)
        print('got element list')
        print(element_list)

        #print("Total unique elements: ", len(element_list))
        # Get values at integration points for these specific elements
        elementSet = odb.rootAssembly.instances['E_GMVC-1'].ElementSetFromElementLabels(
            name='MyElements', 
            elementLabels=element_list
        )
        LE_at_elements = field_output.getSubset(region=elementSet, position=INTEGRATION_POINT)
        print('got LE at elements')
        print(LE_at_elements)
        
        # counters and containers
        max_strain_array=[]
        i=0
        # loop through the LE values at integration points
        for v in LE_at_elements.values:
            i+=1
            element_label = v.elementLabel
            integration_point = v.integrationPoint
            max_principal_strain=v.maxPrincipal
            max_strain_array.append(max_principal_strain)
            if i < 3:
                #print('v.data:', v.data) [LE11, LE22, LE33, LE12, LE23, LE13]
                print('element_label:', element_label)
                print('integration_point:', integration_point)
                print('principal_strain:', max_principal_strain)
       
        # Search for the maximum strain value
        max_strain = max(max_strain_array)
        print('Max strain:', max_strain)
        # append to file
        is_new_file = not os.path.exists('max_strain.csv') or os.path.getsize('max_strain.csv') == 0
        with open('max_strain.csv', 'a') as f:
            # Write header
            if is_new_file:
                f.write('mesh_size,max_strain,radius\n')
            
            # Write data based on odb path
            print('odb:', odb)
            f.write('{0}, {1}, {2}\n'.format(mesh_size, max_strain, radius))
            
        
        return nodes_within_radius, max_strain
    
    except Exception as e:
        print('Error:', e)
        return None
    
    finally:
        odb.close()
        print('closed odb')


# Test the function
odb_filepath='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-099\Job-99.odb' 
odb_filepath2='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-100\Job-100.odb' 
odb_filepath3='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-102\Job-102.odb' 
odb_filepath4='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-109\Job-109.odb' #7mm mesh
odb_filepath5='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-108\Job-108.odb' #4mm mesh
step_name='Step-2'
frame_number=20
center_point = (34, -52, 18)
radius = 3
nodes_within_radius, max_strain = get_max_strain(odb_filepath, step_name, frame_number, center_point, radius, mesh_size=10)
nodes_within_radius, max_strain = get_max_strain(odb_filepath4, step_name, frame_number, center_point, radius, mesh_size=7)
nodes_within_radius, max_strain = get_max_strain(odb_filepath2, step_name, frame_number, center_point, radius, mesh_size=5)
nodes_within_radius, max_strain = get_max_strain(odb_filepath5, step_name, frame_number, center_point, radius, mesh_size=4)
nodes_within_radius, max_strain = get_max_strain(odb_filepath3, step_name, frame_number, center_point, radius, mesh_size=3)



