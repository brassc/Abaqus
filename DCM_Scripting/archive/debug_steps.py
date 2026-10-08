# Debug script - run each step, then test manual set creation
# After each step, try creating a manual set to see if model is corrupted

print("Step 1: Get model")
modelDB = mdb.models[MODEL_NAME]
print("Done - test manual set creation now")
# STOP HERE FIRST - comment out below, test, then uncomment next section

print("Step 2: Get assembly")
assembly = modelDB.rootAssembly
print("Done - test manual set creation now")
# STOP HERE - comment out below

print("Step 3: Get instance")
instance = assembly.instances[INSTANCE_NAME]
print("Done - test manual set creation now")
# # STOP HERE - comment out below

print("Step 4: Get all_nodes")
all_nodes = instance.nodes
print("Found {} nodes".format(len(all_nodes)))
print("Done - test manual set creation now")
# # STOP HERE - comment out below

# print("Step 5: Iterate through nodes")
count = 0
for node in all_nodes:
    coords = node.coordinates
    count += 1
print("Iterated through {} nodes".format(count))
print("Done - test manual set creation now")
