# Minimal test - does this create flat sets?
# Run: execfile('test_set_creation.py')

print("test_set_creation.py starting")

m = mdb.models[MODEL_NAME]
inst = m.rootAssembly.instances[INSTANCE_NAME]

# Create 2 test sets with 100 nodes each
labels1 = [inst.nodes[i].label for i in range(100)]
labels2 = [inst.nodes[i].label for i in range(100, 200)]

m.rootAssembly.SetFromNodeLabels(name='SCRIPT_TEST_1', nodeLabels=((inst.name, labels1),))
m.rootAssembly.SetFromNodeLabels(name='SCRIPT_TEST_2', nodeLabels=((inst.name, labels2),))

print("Done - check if SCRIPT_TEST_1 and SCRIPT_TEST_2 are flat")
