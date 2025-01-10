from odbAccess import *
from abaqusConstants import *
import math
import numpy as np
from abaqus import *
import sys

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



def get_nodes_from_element(odb, element_id, instance):
    try:
        element = instance.elements[element_id]
        node_indices = element.connectivity
        node_ids = []
        for node_idx in node_indices:
            if node_idx >= len(instance.nodes):
                print("Warning: Node index %d out of range" % node_idx)
                print("element_id: ", element_id)
                continue
            try:
                node = instance.nodes[node_idx]
                node_ids.append(node.label)
            except:
                print("Failed to get node at index %d" % node_idx)
                
        return node_ids
        
    except Exception as e:
        print("Error processing element %d: %s" % (element_id, str(e)))
        return []

def get_strain_at_point(odb_path, x, y, z):
    print('opening odb...')
    #print('reload test')
    odb = openOdb(path=odb_path)
    print('odb opened')
    interestFrame = odb.steps['Step-2'].frames[24]  # For frame 25
    print('frame of interest selected')
    # Get strain field
    print('getting strain field...')
    strain = interestFrame.fieldOutputs['LE']  # 'LE' for logarithmic strain
    print("Available field outputs in frame:", interestFrame.fieldOutputs.keys())
    print('strain field obtained')
    
    # Find the closest node to your coordinates
    print('finding closest node to coordinates...')
    assembly = odb.rootAssembly
    all_nodes = assembly.instances['E_GMVC-1'].nodes  
    
    min_distance = float('inf')
    #closest_node = None
    node_distances = []
    
    for node in all_nodes:
        coords = node.coordinates
        distance = ((coords[0]-x)**2 + (coords[1]-y)**2 + (coords[2]-z)**2)**0.5
        
        # Calculate direction (vector from node to target point)
        direction = [x - coords[0], y - coords[1], z - coords[2]]
        
        # Normalize the direction vector (unit vector)
        magnitude = math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2)
        if magnitude != 0:
            direction = [d / magnitude for d in direction]  # Normalize
            
            
        node_distances.append((node, distance, direction))
        #if distance < min_distance:
        #    min_distance = distance
        #    closest_node = node
    
    # Sort nodes by distance and get the top 5 closest
    node_distances.sort(key=lambda pair: pair[1])
    closest_nodes = node_distances[:5]  # Get the top 5 closest nodes
    
    
    print('closest nodes found.')
    #print(closest_node)
    """
    for i, (node, distance, direction) in enumerate(closest_nodes):
        print('Node {}: Label {}, Distance {:.6f}, Direction: ({:.6f}, {:.6f}, {:.6f})'.format(
            i+1, node.label, distance, direction[0], direction[1], direction[2]
        ))
    """
    
    """
    # Get strain value at closest node
    print('getting strain values at closest node...')
    strain_values = strain.getSubset(position=NODAL)
    print('strain_values')
    print(strain_values)
    print(strain_values.validInvariants)
    #max_principal_values = strain_values.getSubset(invariant=MAX_PRINCIPAL, position=NODAL)
    #print(max_principal_values)
    print('hello')
    #for value in strain_values.validInvariants:
    print('Type of strain_values.values:', type(strain_values.values))

    # Check if strain_values is None before proceeding
    if strain_values is None:
        print('strain_values is None. Please check the getSubset() method call.')
    else:
        print('Type of strain_values.values:', type(strain_values.values))

        # Check if strain_values.values is an instance of FieldValueArray by type string
        if str(type(strain_values.values)) == "<type 'FieldValueArray'>":
            print('inside if')
            print('Number of values in strain_values:', len(strain_values.values))

            # Loop through the FieldValueArray using indexing
            for i in range(len(strain_values.values)):
                print('inside for')
                value = strain_values.values[i]  # Access individual FieldValue
                print('Node Label:', value.nodeLabel)
                print('Strain Data:', value.data)
        else:
            print('strain_values.values is not a FieldValueArray or iterable')
    
    for value in strain_values.values:
        print('inside for')
        print(value)
        break
        if value.nodeLabel == closest_node:
            print('inside if')
            print('strain data:')
            print(value.data)
            break
            #return value.data  # Returns strain components
        print('no strain data found for closest node')
    """
    odb.close()
    

# def get_max_strain_in_region(odb_path, step_name, frame_number, center_point, radius):
#     """
#     Extract maximum strain value within a spherical region of interest.
    
#     Parameters:
#     odb_path (str): Path to the .odb file
#     step_name (str): Name of the step to analyze
#     frame_number (int): Frame number to analyze
#     center_point (tuple): (x, y, z) coordinates of region center
#     radius (float): Radius of spherical region of interest
    
#     Returns:
#     tuple: (max_strain_value, node_id)
#     """
#     # Open the database
#     print('opening database...')
#     odb = openOdb(path=odb_path, readOnly=True)
#     print('opened odb')
#     try:
#         # Get the step and frame
#         print('get step')
#         step = odb.steps[step_name]
#         print('get frame')
#         frame = step.frames[frame_number]
#         print('step:', step)
#         print('frame:', frame)
        
#         # Get the logarithmic strain field output and instance
#         print('get LE strain field')
#         strain_field = frame.fieldOutputs['LE']
#         print(strain_field)
#         instance = odb.rootAssembly.instances['E_GMVC-1']
#         instance_strain = strain_field.getSubset(region=instance)
        
#         print('total strain values: ', len(instance_strain.values))
        
        
#         # Convert center_point to numpy array
        
#         #print("center point received:", center_point)
        
#         center = np.array(center_point, dtype=float)
#         print("center point np array: ", center)

#         # Initialize tracking variables
#         max_strain = float('-inf')
#         max_strain_element = None
#         max_int_point = None
        
#         # Keep track of processed elements to get coordinates only once
#         element_centers = {}
        
#         # Debug counters
#         processed_elements = 0
#         elements_in_range = 0
#         processed_element_ids = set()
        
#         # Iterate through strain values
#         for i, value in enumerate(instance_strain.values):
#             element_id = value.elementLabel
#             processed_elements += 1
            
#             # Skip if we've already processed this element
#             if element_id in processed_element_ids:
#                 continue
        
#             processed_element_ids.add(element_id)
            
#             # Calculate element center only once per element
#             if element_id not in element_centers:
#                 element = instance.elements[element_id - 1]
#                 coords = np.zeros(3)
#                 valid_nodes = 0
                
#                 node_ids = get_nodes_from_element(odb_path, element_id, instance)
#                 #print(node_ids)
#                 element_center = np.zeros(3)
#                 valid_count = 0

#                 for node_id in node_ids:
#                     try:
#                         node = instance.getNodeFromLabel(node_id)
#                         coords = np.array(node.coordinates, dtype=float)
#                         #print("Node", node_id)
#                         #print("coords",coords)
#                         element_center += coords
#                         valid_count += 1
#                     except:
#                         print("Failed to get node %d", node_id)
        
#                 if valid_count > 0:
#                     element_center = element_center / valid_count
#                     #print("Calculated center: ",actual_center)
#                     #print(type(actual_center))
#                     element_centers[element_id] = element_center #  element name , actual_center
                        
                
             
            

            
#             # Get cached element center and check distance
#             if element_id in element_centers:
#                 element_center = element_centers[element_id]
#                 #if element_id == 21595:
#                 #    print("center before distance calc: ", center)
#                 #distance=euclidean_distance_3d(element_center, center)
#                 distance = euclidean_distance_3d(
#                         element_center[0], element_center[1], element_center[2],
#                         center[0], center[1], center[2]
#                 )
#                 #if element_id == 21595:
#                 #    print("center after distance calc : ", center)
                  
#                 # Debug print for this specific element
#                 if element_id == 21595:
#                     print("\nDistance calculation debug:")
#                     print("Original element center coordinates were around (13, 8, 5)")
#                     print("Stored element center: ", element_centers[element_id])
#                     print("Target center: ", center)
#                     print("center0", center[0])
#                     print("center1", center[1])
#                     print("center2", center[2])
#                     print("Raw distance:", distance)# np.linalg.norm(element_center.astype(float) - center))           
#                     print("distance type: ", type(distance))
                
                
                
#                 if np.abs(distance) <= radius:
#                     elements_in_range += 1
#                     # Get maximum absolute strain component
#                     strain_components = value.data
#                     max_component = max(abs(comp) for comp in strain_components)
#                     """
#                     print("Element within radius:")
#                     print("Element:", element_id)
#                     print("Center:", element_center)
#                     print("Distance: ", distance)
#                     print("Max component: ", max_component)
#                     print("\n")
#                     """
#                     if max_component > max_strain:
#                         max_strain = max_component
#                         max_strain_element = element_id
#                         max_int_point = i % 4  # C3D10M typically has 4 integration points
#                         print("Max strain updated to: ", max_strain)
#                         print("Element: ", max_strain_element)


#         print("\nProcessing Summary:")
#         print("Total elements processed: %d" % processed_elements)
#         print("Elements within radius: %d" % elements_in_range)\
        
#         print("max strain: ", max_strain)
#         print("max strain element:", max_strain_element)
        
        
#         return max_strain, max_strain_element, max_int_point
    
#     finally:
#         # Always close the database
#         odb.close()


def get_max_strain_in_region(odb_path, step_name, frame_number, center_point, radius):
    """
    Extract maximum strain value within a spherical region of interest.
    Added enhanced debugging and flexible radius.
    """
    print('Opening database...')
    odb = openOdb(path=odb_path, readOnly=True)
    print('Opened odb')
    
    try:
        step = odb.steps[step_name]
        frame = step.frames[frame_number]
        strain_field = frame.fieldOutputs['LE']
        instance = odb.rootAssembly.instances['E_GMVC-1']
        instance_strain = strain_field.getSubset(region=instance)
        
        # Sample node coordinates near target point
        center = np.array(center_point, dtype=float)
        print("\nSearch Parameters:")
        print("Center point: %s" % str(center))
        print("Search radius: %s" % str(radius))
        
        # Initialize tracking variables
        max_strain = float('-inf')
        max_strain_element = None
        max_int_point = None
        element_centers = {}
        
        # Debug counters and storage
        processed_elements = 0
        elements_in_range = 0
        processed_element_ids = set()
        closest_element_distance = float('inf')
        closest_element_id = None
        closest_element_coords = None
        
        # Distance tracking
        distances = []
        
        # Modified element processing loop
        for i, value in enumerate(instance_strain.values):
            element_id = value.elementLabel
            
            if element_id in processed_element_ids:
                continue
                
            processed_element_ids.add(element_id)
            processed_elements += 1
            
            if element_id not in element_centers:
                element = instance.elements[element_id - 1]
                node_ids = get_nodes_from_element(odb_path, element_id, instance)
                element_center = np.zeros(3)
                valid_count = 0
                
                for node_id in node_ids:
                    try:
                        node = instance.getNodeFromLabel(node_id)
                        coords = np.array(node.coordinates, dtype=float)
                        element_center += coords
                        valid_count += 1
                    except:
                        continue
                
                if valid_count > 0:
                    element_center = element_center / valid_count
                    element_centers[element_id] = element_center
            
            # Distance calculation and tracking
            if element_id in element_centers:
                element_center = element_centers[element_id]
                distance = euclidean_distance_3d(
                    element_center[0], element_center[1], element_center[2],
                    center[0], center[1], center[2]
                )
                
                # Track closest element
                if distance < closest_element_distance:
                    closest_element_distance = distance
                    closest_element_id = element_id
                    closest_element_coords = element_center
                
                distances.append(distance)
                
                if np.abs(distance) <= radius:
                    elements_in_range += 1
                    strain_components = value.data
                    max_component = max(abs(comp) for comp in strain_components)
                    
                    if max_component > max_strain:
                        max_strain = max_component
                        max_strain_element = element_id
                        max_int_point = i % 4
        
        # Print distance statistics
        distances = np.array(distances)
        print("\nDistance Statistics:")
        print("Minimum distance to any element: %.2f" % np.min(distances))
        print("Mean distance to elements: %.2f" % np.mean(distances))
        print("Closest element:")
        print("- ID: %s" % str(closest_element_id))
        print("- Coordinates: %s" % str(closest_element_coords))
        print("- Distance: %.2f" % closest_element_distance)
        
        # Suggest adjusted radius
        if elements_in_range == 0:
            suggested_radius = np.percentile(distances, 1)  # radius that would capture 1% of elements
            print("\nNo elements found in current radius (%.2f)" % radius)
            print("Suggested minimum radius to capture elements: %.2f" % suggested_radius)
        
        print("\nProcessing Summary:")
        print("Total elements processed: %d" % processed_elements)
        print("Elements within radius: %d" % elements_in_range)
        print("Max strain: %s" % str(max_strain))
        print("Max strain element: %s" % str(max_strain_element))
        
        return max_strain, max_strain_element, max_int_point
    
    finally:
        odb.close()


        
odb_filepath='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-099\Job-99.odb'  
#odb_filepath='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-100\Job-100.odb'
#odb_filepath='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-102\Job-102.odb'

step_name='Step-2'
frame_number=25
#center_point=(34, -52, 18)
center_point=(34, -52, 18)
radius=5
#strain_values=get_strain_at_point(odb_filepath, 35, -52, 20)

get_max_strain_in_region(odb_filepath, step_name, frame_number, center_point, radius)


