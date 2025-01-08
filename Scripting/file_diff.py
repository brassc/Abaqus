def compare_files(file1_path, file2_path, output_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        # Read lines from both files
        file1_lines = file1.readlines()
        file2_lines = file2.readlines()

    # Calculate differences
    differences = []
    for line in file1_lines:
        if line not in file2_lines:
            differences.append(f"- {line}")
    
    for line in file2_lines:
        if line not in file1_lines:
            differences.append(f"+ {line}")
    
    print(differences)
    
    # Write differences to output file
    with open(output_path, 'w') as output_file:
        output_file.writelines(differences)

    print(f"Differences written to {output_path}")
    
    
    

# Example usage
file1_path = r'C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-045\Job-45.inp'
file2_path = r'C:\Users\cmb247\ABAQUS\K_DC_FALX\K-DCBH-046\Job-46.inp'
output_path = r'C:\Users\cmb247\repos\Abaqus\Scripting\differences.txt'

compare_files(file1_path, file2_path, output_path)
