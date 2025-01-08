from abaqus import *
from abaqusConstants import *
from node_set_export import create_file_from_node_set
import sys


print('*************')
class Face:
    def __init__(self, featureName, index, instanceName, isReferenceRep, pointOn):
        self.featureName = featureName
        self.index = index
        self.instanceName = instanceName
        self.isReferenceRep = isReferenceRep
        self.pointOn = pointOn

    def __repr__(self):
        return 'Face(featureName={}, index={}, instanceName={}, isReferenceRep={}, pointOn={})'.format(
            self.featureName, self.index, self.instanceName, self.isReferenceRep, self.pointOn
        )
    # Add any additional methods as needed

# Create a custom face object with specified attributes
face_attributes = {
    'featureName': None,
    'index': 1002,
    'instanceName': 'DC_Skull-1',
    'isReferenceRep': False,
    'pointOn': ((48.738136, 71.184865, 4.793374), (0.615942, 0.787508, 0.021136))
}

custom_face = Face(**face_attributes)

# Print the custom face object to verify its attributes
print(custom_face)
print('********')


# Get the model database
modelDB = mdb.models['Model-1']
assembly=modelDB.rootAssembly

# surface names
large_surf_name = 'test-big'#'Skull_Surf'
small_surf_name = 'test-small'#'sticky_skull'

# get surfaces
large_surf=assembly.surfaces[large_surf_name]
small_surf=assembly.surfaces[small_surf_name]

for face in large_surf.faces[:]:
    if face in small_surf.faces:
        # remove it 
        large_surf.faces.remove(face)



exit()
# find faces in large surf that are NOT in the small surf
faces_to_include=[]

for face in large_surf.faces:
    if face not in small_surf.faces:
        custom_face = Face(
        featureName=face.featureName,
        index=face.index,
        instanceName=face.instanceName,
        isReferenceRep=face.isReferenceRep,
        pointOn=face.pointOn
        )
        faces_to_include.append(custom_face)

print('************')
print(faces_to_include)


# Create new surf for the non-overlapping region at the assembly level
new_surf_name='not_sticky_skull'
# Print list of faces in the large surface

"""
print("Faces in large surface:")
for face in large_surf.faces:
    print(face)
"""
# Print list of faces in the small surface
print("\nFaces in small surface:")
for face in small_surf.faces:
    print(face)


# Print faces to include
print("\nFaces to include:")
for face in faces_to_include:
    #print(f"ID: {face.id}, PointOn: {face.pointOn}")
    print(face)
   

assembly.Surface(name=new_surf_name, side1Faces=faces_to_include)










