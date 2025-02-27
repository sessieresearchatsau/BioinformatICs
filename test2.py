import requests
import nglview as nv
from IPython.display import display
import py3Dmol

# Function to search for a protein structure in AlphaFold DB
def fetch_alphafold_structure(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        pdb_filename = f"{uniprot_id}.pdb"
        with open(pdb_filename, "w") as pdb_file:
            pdb_file.write(response.text)
        return pdb_filename
    else:
        print("Error: Protein not found in AlphaFold DB.")
        return None

# Function to visualize the protein structure
def visualize_protein(pdb_file):
    viewer = py3Dmol.view(width=600, height=400)
    with open(pdb_file, "r") as f:
        pdb_data = f.read()
    viewer.addModel(pdb_data, "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.zoomTo()
    return viewer.show()

# Example: Search for a protein (e.g., UniProt ID "P12345")
protein_id = input("Enter UniProt ID: ")
pdb_file = fetch_alphafold_structure(protein_id)

if pdb_file:
    visualize_protein(pdb_file)