#This function extracts nodes from assembly set and exports to .txt as comma separated list


from abaqus import mdb

# Get the model database
modelDB = mdb.models['Model-1']

def create_file_from_node_set(set_name):

    nodeSet = modelDB.rootAssembly.sets[set_name]

    # Get the list of nodes in the node set
    nodes = nodeSet.nodes

    # Specify the output file path (in the working directory)
    output_file_path = str(set_name) + '_NodeList.txt'

    #write comma separated node labels to file
    with open(output_file_path, 'w') as output_file:
        node_labels = [str(node.label) for node in nodes]
        output_file.write(','.join(node_labels))

    # Print a confirmation message
    print("Node labels have been exported to " + str(output_file_path))
