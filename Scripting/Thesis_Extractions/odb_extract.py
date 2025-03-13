#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comprehensive extraction of data from an Abaqus .odb file.
Saves data frame-by-frame to compressed pickle files (.pkl.gz) to preserve full
numerical precision while minimizing file size. Creates a directory structure
with model information and individual frame files for efficient storage and access.

Run this script using Abaqus Python (Python 2.7):
abaqus python extract_odb_data.py model.odb output_directory LE
"""

import os
import sys
import pickle
import gzip
import time
from abaqusConstants import *
from odbAccess import *

def extract_frame_data(odb_file, output_dir=None, field_filter=None, max_frames=None):
    """
    Extract data from an Abaqus .odb file one frame at a time.
    
    Parameters:
    -----------
    odb_file : str
        Path to the Abaqus .odb file
    output_dir : str
        Directory to save the extracted frame data
    field_filter : list
        List of field names to extract (e.g., ['LE']). If None, extract all fields.
    max_frames : int
        Maximum number of frames to extract (useful for testing)
    """
    print("Opening ODB file: {}".format(odb_file))
    start_time = time.time()
    odb = openOdb(path=odb_file, readOnly=True)
    
    # Create output directory if not specified
    if output_dir is None:
        output_dir = os.path.splitext(odb_file)[0] + "_frames"
    
    # Create directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract basic model data (nodes, elements, sets)
    print("Extracting model structure...")
    model_data = {
        'model_name': odb.name,
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
        model_data['node_coordinates'][instance_name] = node_coords
        
        # Store element connectivity
        element_conn = {}
        for element in instance.elements:
            element_conn[element.label] = {
                'connectivity': element.connectivity,
                'type': element.type
            }
        model_data['element_connectivity'][instance_name] = element_conn
    
    # Extract node and element sets
    print("Extracting node and element sets...")
    for set_name, node_set in odb.rootAssembly.nodeSets.items():
        nodes = []
        for node in node_set.nodes:
            for n in node:
                nodes.append((n.instanceName, n.label))
        model_data['node_sets'][set_name] = nodes
    
    for set_name, element_set in odb.rootAssembly.elementSets.items():
        elements = []
        for element in element_set.elements:
            for e in element:
                elements.append((e.instanceName, e.label))
        model_data['element_sets'][set_name] = elements
    
    # Save model structure
    model_file = os.path.join(output_dir, "model_structure.pkl.gz")
    print("Saving model structure to: {}".format(model_file))
    with gzip.open(model_file, 'wb', compresslevel=9) as f:
        pickle.dump(model_data, f, protocol=2)
    
    # Process each step in the model
    total_frames = 0
    extracted_frames = 0
    
    for step_name, step in odb.steps.items():
        print("\nProcessing step: {}".format(step_name))
        print("  Time period: {}".format(step.timePeriod))
        print("  Number of frames: {}".format(len(step.frames)))
        
        # Create step directory
        step_dir = os.path.join(output_dir, step_name)
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        
        # Process each frame in the step
        for frame_idx, frame in enumerate(step.frames):
            # Check if we've reached max frames (if specified)
            if max_frames is not None and extracted_frames >= max_frames:
                print("Reached maximum number of frames ({}). Stopping.".format(max_frames))
                break
            
            frame_start_time = time.time()
            print("  Processing frame {}/{}: time={}".format(
                frame_idx + 1, len(step.frames), frame.frameValue))
            
            # Frame metadata
            frame_data = {
                'frame_value': frame.frameValue,
                'description': frame.description,
                'field_outputs': {}
            }
            
            # Process each field output in the frame
            has_fields = False
            for field_name, field in frame.fieldOutputs.items():
                # Skip fields not in the filter if a filter is provided
                if field_filter and field_name not in field_filter:
                    continue
                
                print("    Extracting field: {}".format(field_name))
                has_fields = True
                
                field_data = {
                    'description': field.description,
                    'type': str(field.type),
                    'values': []
                }
                
                # Store component labels if available
                if hasattr(field, 'componentLabels'):
                    field_data['component_labels'] = field.componentLabels
                
                # Extract field values
                for value in field.values:
                    try:
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
                        
                        field_data['values'].append(val_data)
                    except Exception as e:
                        print("      Warning: Skipping a value due to error: {}".format(e))
                        continue
                
                # Add field data to frame
                frame_data['field_outputs'][field_name] = field_data
            
            # Skip saving if no fields were extracted
            if not has_fields:
                print("    No requested fields found in this frame. Skipping.")
                continue
            
            # Save frame data to pickle file
            frame_file = os.path.join(step_dir, "frame_{:04d}.pkl.gz".format(frame_idx))
            print("    Saving frame data to: {}".format(frame_file))
            frame_save_start = time.time()
            with gzip.open(frame_file, 'wb', compresslevel=9) as f:
                pickle.dump(frame_data, f, protocol=2)
            
            frame_save_time = time.time() - frame_save_start
            frame_time = time.time() - frame_start_time
            print("    Frame extracted in {:.2f} seconds (save: {:.2f} s)".format(frame_time, frame_save_time))
            
            extracted_frames += 1
            total_frames += 1
    
    # Close the ODB file
    odb.close()
    
    # Print summary
    total_time = time.time() - start_time
    print("\nExtraction Summary:")
    print("  Model: {}".format(model_data['model_name']))
    print("  Total frames extracted: {}".format(total_frames))
    print("  Total time: {:.2f} seconds ({:.2f} minutes)".format(total_time, total_time/60))
    
    if field_filter:
        print("  Fields extracted: {}".format(", ".join(field_filter)))
    
    print("  Data saved to: {}".format(output_dir))
    
    # Create a manifest file with extraction information
    manifest = {
        'model_name': model_data['model_name'],
        'extraction_time': total_time,
        'total_frames': total_frames,
        'field_filter': field_filter,
        'steps': {}
    }
    
    # Add step information to manifest
    for step_name in os.listdir(output_dir):
        step_dir = os.path.join(output_dir, step_name)
        if os.path.isdir(step_dir) and step_name != '__pycache__':
            frame_files = [f for f in os.listdir(step_dir) if f.endswith('.pkl.gz')]
            manifest['steps'][step_name] = len(frame_files)
    
    # Save manifest
    manifest_file = os.path.join(output_dir, "manifest.pkl.gz")
    with gzip.open(manifest_file, 'wb', compresslevel=9) as f:
        pickle.dump(manifest, f, protocol=2)
    
    return output_dir

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
    
    # Get max frames if needed (for testing)
    max_frames = None
    if len(sys.argv) > 4:
        try:
            max_frames = int(sys.argv[4])
        except ValueError:
            pass
    else:
        max_frames_input = raw_input("Enter maximum number of frames to extract (blank for all): ")
        if max_frames_input.strip():
            try:
                max_frames = int(max_frames_input)
            except ValueError:
                pass
    
    
    if max_frames:
        print("Extracting at most {} frames".format(max_frames))
        
        
    # Extract data
    extract_frame_data(odb_file, output_file, field_filter, max_frames)

if __name__ == "__main__":
    main()
