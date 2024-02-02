import re

# Function to parse the result.txt file
def parse_result_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Start data extraction after the specified line
    data_start_line = 'Perform docking|================done=================|'.rstrip()
    start_index = content.find(data_start_line) + len(data_start_line) + 1
    content = content[start_index:]

    # Extract ligand information
    ligand_matches = re.findall(r'Refining ligand \./test/set17_1_out/(\d+) results(.*?)(?=(?:Refining ligand \./test/set17_1_out/|\Z))', content, re.DOTALL)

    # Organize the data into a list of dictionaries
    ligand_data = [{'ligand_id': match[0], 'result_data': match[1].strip()} for match in ligand_matches]

    return ligand_data

# Function to get the top N ligands based on docking score
def get_top_ligands(ligand_data, top_n=10):
    processed_ligands = 0  # To count the number of ligands processed

    # Extracting docking score from result_data with support for decimal and negative scores
    for ligand in ligand_data:
        scores = re.findall(r'\d+\s+(-?\d+(?:\.\d+)?)\s+\d+\.\d+\s+\d+\.\d+', ligand['result_data'])
        if scores:
            ligand['docking_score'] = float(scores[0])
            processed_ligands += 1

    # Remove ligands without docking score
    ligand_data = [ligand for ligand in ligand_data if 'docking_score' in ligand]

    # Sort ligands based on docking score
    sorted_ligands = sorted(ligand_data, key=lambda x: x['docking_score'], reverse=False)

    # Print the number of ligands processed
    print(f"\nNumber of ligands processed: {processed_ligands}")

    return sorted_ligands[:top_n]

# Specify the path to the result.txt file
result_file_path = './result.txt'

# Parse the result file
ligand_data = parse_result_file(result_file_path)

# Get user input for the number of top ligands
top_n = int(input("Enter the number of top ligands to display (default 10): ") or 10)

# Get the top N ligands
top_ligands = get_top_ligands(ligand_data, top_n)

# Display the results
print("\nTop {} Ligands:".format(top_n))
for i, ligand in enumerate(top_ligands, start=1):
    print(f"{i}. Ligand {ligand['ligand_id']} - Docking Score: {ligand['docking_score']:.2f}")
