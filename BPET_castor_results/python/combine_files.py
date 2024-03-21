from tqdm import tqdm

# Path to input files
path = "../Files/"
# List of input files
input_files = [path+"intersect_1th_points_with_time.txt", path + "intersect_2th_points_with_time.txt",path+"intersect_3th_points_with_time.txt", path + "intersect_4th_points_with_time.txt"]
output_file = path + 'LOR_all_points_time.txt'
print("Input files: ", input_files)
print("Output file: ", output_file)

# Initialize a variable to track if the header has been written
header_written = False

# Count the total number of lines across all files
total_files = len(input_files)

print("Combining files...")
# Open the output file in write mode
with open(output_file, 'w') as out_file:
    for input_file in tqdm(input_files, unit='file', desc="Combining Files"):
        with open(input_file, 'r') as in_file:
            lines = in_file.readlines()
            if not header_written:
                # Write the header from the first file
                out_file.writelines(lines)  
                header_written = True
            else:
                # Skip the header for subsequent files
                out_file.writelines(lines[1:])  

print("Files combined successfully.")
