from abaqus import *
from abaqusConstants import *


def create_node_set_from_file(file_path, new_node_set_name):

    try:
        # Open the file and read the node numbers
        with open(file_path, 'r') as file:
            node_numbers = file.read().split(',')

        ## Convert the node numbers to integers
        node_numbers = [int(node_number) for node_number in node_numbers]
        #print(node_numbers)

        # Get the model database
        modelDB = mdb.models['Model-1']
    
        # Create node set from node_numbers
        # Create the assembly node set directly from the list of node labels
        
        # Create the assembly node set from node labels
        print('Creating new assembly node set from node labels...)
        modelDB.rootAssembly.SetFromNodeLabels(name=new_node_set_name, nodeLabels=(('E_GMVC-1', node_numbers), 
        ))
    

    # Print a confirmation message
        print("Created a new node set '" + new_node_set_name + "' with " + str(len(node_numbers)) + " nodes.")
        print('Saving model database...')
        mdb.save()
        print('Save complete.')
        
    except IOError:
        print("File '" + file_path + "' not found or an error occurred while reading it. Please ensure the file exists.")
    except ValueError:
        print("Error converting node numbers to integers. Please ensure the file contains valid integers separated by commas.")
    except KeyError:
        print("Error creating the new node set. Please check your model and script.")
     