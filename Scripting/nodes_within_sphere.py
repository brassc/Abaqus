from abaqus import *
from abaqusConstants import *
from node_set_import import create_node_set_from_file


def create_node_set_within_sphere(existing_node_set_name, radius, filename, center_x, center_y, center_z):
    # Get the model database
    modelDB = mdb.models['Model-1']
    
    # Create a new node set for nodes within the sphere
    new_node_set_name = 'NodesWithinSphere'
    nodes_within_sphere = []

    # Get the coordinates of the nodes in the existing node set
    existing_node_set = modelDB.rootAssembly.sets[existing_node_set_name]
    existing_nodes = existing_node_set.nodes

    # Open file into which to read nodes
    with open(filename, 'w') as txt_file:
        # Check each node in the existing set if it falls within the sphere
        for node in existing_node_set.nodes:
            x, y, z = node.coordinates
            distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5  # Euclidean distance
            if distance <= radius:
                nodes_within_sphere.append(node)
        txt_file.write(','.join(str(node.label) for node in nodes_within_sphere))
    # Print a confirmation message
    

    

    # Write the node labels to a .txt file separated by commas
    #node_labels = [str(node.label) for node in nodes_within_sphere]
    #with open(filename, 'w') as txt_file:
    #    txt_file.write(','.join(node_labels))

    print("Node labels within the sphere have been saved to '{}'.".format(filename))
    print("Creating new set...")
    create_node_set_from_file(filename, new_node_set_name)
    message = "Created a new node set '{}' with nodes within a sphere of radius {} mm at center ({}, {}, {})."
    print(message.format(new_node_set_name, radius, center_x, center_y, center_z))

