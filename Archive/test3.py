import requests
import re
import py3Dmol
from IPython.display import display, HTML

def fetch_alphafold_pdb(uniprot_id):
    """
    Fetches the AlphaFold-predicted PDB for the given UniProt ID from the AlphaFold database.
    Returns the text content of the PDB or None if not found.
    """
    # AlphaFold DB PDB URL format
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"[Error] Could not find PDB for UniProt ID '{uniprot_id}' on AlphaFold.")
        return None

def parse_pdb_for_amino_acids(pdb_text):
    """
    Parses the PDB text to extract the sequence of amino acids in order.
    Returns a list of tuples [(resName, resSeq), ...] e.g. [("ALA", 1), ("ARG", 2), ...].
    
    We look at ATOM lines and extract:
      - residue name (resName)
      - residue sequence number (resSeq)
    We'll track the first occurrence of each (resSeq) to build the chain.
    """
    amino_acid_list = []
    seen_residues = set()
    
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            # Columns from PDB format:
            # Residue name is in columns 17-19 (in python, line[17:20])
            # Residue sequence number is in columns 22-26 (line[22:26]) 
            # but we need to strip spaces
            res_name = line[17:20].strip()
            res_seq = line[22:26].strip()
            
            # Use an integer for res_seq
            res_seq = int(res_seq)
            
            # Check if we already added this residue
            if (res_name, res_seq) not in seen_residues:
                seen_residues.add((res_name, res_seq))
                amino_acid_list.append((res_name, res_seq))
                
    # Sort by residue number (in case they're not strictly in order)
    amino_acid_list.sort(key=lambda x: x[1])
    
    return amino_acid_list

def build_amino_acid_graph(amino_acid_list):
    """
    Builds a list of edges representing the consecutive amino acids.
    Each edge is formatted as '{AA# -> AA#}' for printing.
    
    We assume consecutive residues in the list are bonded (like a backbone).
    """
    edges = []
    for i in range(len(amino_acid_list) - 1):
        current_aa, current_num = amino_acid_list[i]
        next_aa, next_num = amino_acid_list[i + 1]
        
        # Format: {ALA1 -> ARG2}
        edge_str = f"{{{current_aa}{current_num} -> {next_aa}{next_num}}}"
        edges.append(edge_str)
    return edges

def display_protein_structure(pdb_text):
    """
    Displays the protein 3D structure using py3Dmol.
    """
    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_text, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    return view.show()

def main():
    print("=== Protein Graph & 3D Structure Viewer ===")
    uniprot_id = input("Enter a UniProt ID (e.g. P69905 for human hemoglobin subunit alpha): ").strip()
    
    # 1) Fetch PDB from AlphaFold
    pdb_data = fetch_alphafold_pdb(uniprot_id)
    if not pdb_data:
        return  # Stop if not found
    
    # 2) Parse PDB to extract amino acid sequence
    amino_acids = parse_pdb_for_amino_acids(pdb_data)
    if not amino_acids:
        print("[Error] Could not parse amino acids from the PDB file.")
        return
    
    # 3) Build graph edges
    edges = build_amino_acid_graph(amino_acids)
    
    # 4) Print Output 1: connections in your requested format
    print("\nOutput 1: Amino acid connections")
    print("=================================")
    print(", ".join(edges))
    
    # 5) Output 2: show an image of the protein (3D structure)
    print("\nOutput 2: 3D Protein Structure")
    print("=================================")
    display_protein_structure(pdb_data)

# Run the script (in a Jupyter environment, call main())
if __name__ == "__main__":
    main()