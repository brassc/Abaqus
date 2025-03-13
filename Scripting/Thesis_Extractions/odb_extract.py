#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comprehensive extraction of data from an Abaqus .odb file.
Saves all results to a Python pickle file to preserve full numerical precision.

Run this script using Abaqus Python (Python 2.7):
abaqus python extract_all_odb_data.py model.odb [output.pkl]
"""

import os
import sys
import pickle
import numpy as np
from abaqusConstants import *
from odbAccess import *

def extract_all_field_data(odb_file, output_file=None, field_filter=None):
    """
    Extract all available field data from an Abaqus .odb file and save to a pickle file.
    
    Parameters:
    -----------
    odb_file : str
        Path to the Abaqus .odb file
    output_file : str
        Path to save the extracted data as pickle file
    """
    print("Opening ODB file: {}".format(odb_file))
    odb = openOdb(path=odb_file, readOnly=True)
    
    # If no output file is specified, create one based on the input file
    if output_file is None:
        output_file = os.path.splitext(odb_file)[0] + "_data.pkl"
    
    # Container for all the data
    data = {
        'model_name': odb.name,
        'steps': {},
        'node_coordinates': {},
        'element_connectivity': {},
        'node_sets': {},
        'element_sets': {}
    }
    
    # Extract node coordinates for all instances
    print("Extracting node coordinates...")
    for instance_name, instance in odb.rootAssembly.instances.items():
        # Store node coordinates
        node_coords = {}
        for node in instance.nodes:
            node_coords[node.label] = {
                'coordinates': node.coordinates
            }
        data['node_coordinates'][instance_name] = node_coords
        
        # Store element connectivity
        element_conn = {}
        for element in instance.elements:
            element_conn[element.label] = {
                'connectivity': element.connectivity,
                'type': element.type
            }
        data['element_connectivity'][instance_name] = element_conn
    
    # Extract node and element sets
    print("Extracting node and element sets...")
    for set_name, node_set in odb.rootAssembly.nodeSets.items():
        nodes = []
        for node in node_set.nodes:
            for n in node:
                nodes.append((n.instanceName, n.label))
        data['node_sets'][set_name] = nodes
    
    for set_name, element_set in odb.rootAssembly.elementSets.items():
        elements = []
        for element in element_set.elements:
            for e in element:
                elements.append((e.instanceName, e.label))
        data['element_sets'][set_name] = elements
    
    # Process each step in the model
    for step_name, step in odb.steps.items():
        print("Processing step: {}".format(step_name))
        
        step_data = {
            'frames': {},
            'description': step.description,
            'time_period': step.timePeriod
        }
        
        # Process each frame in the step
        for frame_idx, frame in enumerate(step.frames):
            print("  Processing frame {}/{}: time={}".format(
                frame_idx + 1, len(step.frames), frame.frameValue))
            
            frame_data = {
                'frame_value': frame.frameValue,
                'description': frame.description,
                'field_outputs': {}
            }
            
            # Process each field output in the frame
            for field_name, field in frame.fieldOutputs.items():
                if field_filter and field_name not in field_filter:
                    continue
                
                print("    Extracting field: {}".format(field_name))
                
                field_data = {
                    'description': field.description,
                    'type': str(field.type),
                    'values': []
                }
                
                # Extract field values
                for value in field.values:
                    # Handle instance name safely
                    instance_name = None
                    if hasattr(value, 'instance') and value.instance is not None:
                        instance_name = value.instance.name
                    
                    val_data = {
                        'node_label': getattr(value, 'nodeLabel', None),
                        'element_label': getattr(value, 'elementLabel', None),
                        'instance_name': instance_name,
                        'integration_point': getattr(value, 'integrationPoint', None),
                        'data': value.data
                    }
                    # If the field has component labels, include them
                    if hasattr(field, 'componentLabels'):
                        field_data['component_labels'] = field.componentLabels
                    
                    field_data['values'].append(val_data)
                    
                # Add field data to frame
                frame_data['field_outputs'][field_name] = field_data
            
            # Add frame data to step
            step_data['frames'][frame_idx] = frame_data
        
        # Add step data to overall data
        data['steps'][step_name] = step_data
    
    # Close the ODB file
    odb.close()
    
    # Save data to pickle file
    print("Saving data to: {}".format(output_file))
    with open(output_file, 'wb') as f:
        pickle.dump(data, f, protocol=2)  # Use protocol 2 for Python 2/3 compatibility
    
    # Print summary
    print("\nExtraction Summary:")
    print("  Model: {}".format(data['model_name']))
    print("  Steps: {}".format(len(data['steps'])))
    
    for step_name, step in data['steps'].items():
        print("  Step '{}': {} frames".format(step_name, len(step['frames'])))
        
        # Count field outputs in the last frame
        last_frame_idx = max(step['frames'].keys())
        last_frame = step['frames'][last_frame_idx]
        print("    Fields in last frame: {}".format(len(last_frame['field_outputs'])))
        
        # List all fields
        print("    Available fields: {}".format(
            ", ".join(last_frame['field_outputs'].keys())))
    
    return output_file

def main():
    """Main function to extract all data from an Abaqus ODB file."""
    # Get parameters from command line or prompt for input
    if len(sys.argv) > 1:
        odb_file = sys.argv[1]
    else:
        odb_file = raw_input("Enter the path to the .odb file: ")
    
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        default_output = os.path.splitext(odb_file)[0] + "_data.pkl"
        output_file = raw_input("Enter the output file path (default: {}): ".format(default_output)) or default_output
    
    # Get field filter if provided
    field_filter = None
    if len(sys.argv) > 3:
        field_filter = sys.argv[3].split(',')
    else:
        field_filter_input = raw_input("Enter field variables to extract (comma-separated, e.g., 'LE,U' or blank for all): ")
        if field_filter_input.strip():
            field_filter = [f.strip() for f in field_filter_input.split(',')]
    
    if field_filter:
        print("Extracting only these fields: {}".format(", ".join(field_filter)))
    else:
        print("Extracting all available fields")
    
    
    # Extract data
    extract_all_field_data(odb_file, output_file, field_filter)

if __name__ == "__main__":
    main()
