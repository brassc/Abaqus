from odbAccess import *
from abaqusConstants import *
import sys

# Open the .odb file
odb = openOdb(path='C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-109\Job-109.odb')

# Access the step and frame information
step = odb.steps['Step-2']  # Change this to your specific step name

# Check all instances and print element sets in each instance
for instance_name, instance in odb.rootAssembly.instances.items():
    print("Instance: {}, Element sets: {}".format(instance_name, list(instance.elementSets.keys())))


instance_name = 'E_GMVC-1'  
element_set_name = 'SET-1'  

# Access the instance and element set by name
instance = odb.rootAssembly.instances[instance_name]
elem_set = instance.elementSets[element_set_name]



# Access the element set by name
#elem_set = odb.rootAssembly.elementSets['E_GMVC-1.SET-1']  

# Initialize a list to hold the total EVOL values for each frame (time step)
total_evol_over_time = []

# Loop over each frame in the step (i.e., each time increment)
for frame in step.frames:
    evol_field = frame.fieldOutputs['EVOL']  # Extract EVOL field data
    total_evol = 0.0

    # Sum the EVOL values only for elements in the specified element set
    for value in evol_field.getSubset(region=elem_set).values:
        total_evol += value.data

    # Store the sum for this frame (time step)
    total_evol_over_time.append(total_evol)
    

# Print the total EVOL values for each frame (time step)
#for i, evol in enumerate(total_evol_over_time):
#    print(f"Frame {i+1}, EVOL sum: {evol}")
  
with open('C:\\Users\\cmb247\\repos\\Abaqus\\Scripting\\evol_7mm.csv', 'w') as f:
    f.write('Frame,EVOL Sum\n')
    for i, evol in enumerate(total_evol_over_time):
        f.write("{},{}\n".format(i+1, evol))


# Close the odb file
odb.close()

