import requests
import io
import networkx as nx
from Bio import PDB
from Bio.Data import IUPACData
import py3Dmol

###############################################################
# 1. FETCH PDB FROM RCSB
###############################################################

def fetch_rcsb_pdb(pdb_id):
    """
    Fetch the official PDB file from RCSB by PDB ID.
    Returns the PDB text or None if not found.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"[Error] Could not find PDB file for '{pdb_id}' on RCSB.")
        return None

###############################################################
# 2. PARSE PDB WITH BIOPYTHON
###############################################################

def parse_structure_from_pdb_text(pdb_text, structure_id="PDB_model"):
    """
    Parse the PDB text into a Biopython Structure object.
    """
    pdb_io = io.StringIO(pdb_text)
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, pdb_io)
    return structure

###############################################################
# 3. BOND / INTERACTION DETECTION
###############################################################

POS_CHARGED = {"ARG", "LYS", "HIS"}  
NEG_CHARGED = {"ASP", "GLU"}
HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO"}

def residue_label(res):
    """
    Returns a label like ALA12 for a Biopython Residue object,
    using resname + residue number.
    """
    return f"{res.get_resname().strip()}{res.get_id()[1]}"

def is_consecutive(r1, r2):
    """
    Returns True if r2 is exactly the next residue number of r1.
    """
    return (r2.get_id()[1] == r1.get_id()[1] + 1)

def parse_disulfide_bonds(pdb_text):
    """
    Returns a list of ((chain1, resnum1), (chain2, resnum2)) from SSBOND lines.
    If the PDB is missing these lines, you'll get an empty list.
    """
    ssbonds = []
    for line in pdb_text.splitlines():
        if line.startswith("SSBOND"):
            # Example line: "SSBOND   1 CYS A    6    CYS A  127"
            chain1 = line[15]
            resnum1 = line[17:21].strip()
            chain2 = line[29]
            resnum2 = line[31:35].strip()
            try:
                resnum1 = int(resnum1)
                resnum2 = int(resnum2)
                ssbonds.append(((chain1, resnum1), (chain2, resnum2)))
            except ValueError:
                pass
    return ssbonds

def is_hbond(atom1, atom2, max_dist=3.5):
    """
    Checks if there's a potential H-bond between atom1 and atom2 by distance
    and basic element check (N or O). Simplified approach.
    """
    e1 = atom1.element.upper()
    e2 = atom2.element.upper()
    if (e1 in ["N","O"] and e2 in ["N","O"]):
        dist = atom1 - atom2
        if dist <= max_dist:
            return True
    return False

def is_ionic(res1, res2, max_dist=5.0):
    """
    Check if res1 is positive and res2 is negative (or vice versa),
    and if any of their side-chain atoms are within max_dist.
    """
    r1name = res1.get_resname()
    r2name = res2.get_resname()
    is_r1_pos = (r1name in POS_CHARGED)
    is_r2_pos = (r2name in POS_CHARGED)
    is_r1_neg = (r1name in NEG_CHARGED)
    is_r2_neg = (r2name in NEG_CHARGED)
    
    # Must be pos <-> neg
    if not ((is_r1_pos and is_r2_neg) or (is_r1_neg and is_r2_pos)):
        return False
    
    backbone_atoms = {"N", "CA", "C", "O"}
    for a1 in res1.get_atoms():
        if a1.get_name() not in backbone_atoms:
            for a2 in res2.get_atoms():
                if a2.get_name() not in backbone_atoms:
                    dist = a1 - a2
                    if dist <= max_dist:
                        return True
    return False

def is_hydrophobic(res1, res2, max_dist=4.5):
    """
    Returns True if both residues are in the HYDROPHOBIC set
    and any side-chain atoms are within max_dist.
    """
    if (res1.get_resname() not in HYDROPHOBIC) or (res2.get_resname() not in HYDROPHOBIC):
        return False
    
    backbone_atoms = {"N", "CA", "C", "O"}
    for a1 in res1.get_atoms():
        if a1.get_name() not in backbone_atoms:
            for a2 in res2.get_atoms():
                if a2.get_name() not in backbone_atoms:
                    dist = a1 - a2
                    if dist <= max_dist:
                        return True
    return False

###############################################################
# 4. BUILD GRAPH & CAPTURE EDGES BY CATEGORY
###############################################################

def build_full_graph(structure, pdb_text):
    """
    Build a graph of residues (nodes) and determine edges for each interaction type.
    Returns a dict of edge-lists keyed by category, e.g.:
      edges_by_category["peptide"] = [(MET1, ALA2), ...]
      edges_by_category["disulfide"] = [(CYS6, CYS127), ...]
      ...
    """
    model = structure[0]
    
    # Collect all standard residues (ignore water, hetero)
    all_residues = []
    for chain in model:
        for residue in chain.get_residues():
            if residue.id[0] == " ":  # standard residue
                all_residues.append((chain.id, residue))
    
    # Sort by chain + residue number
    all_residues.sort(key=lambda x: (x[0], x[1].get_id()[1]))
    
    edges_by_category = {
        "peptide": [],
        "disulfide": [],
        "hydrogen_bond": [],
        "ionic_interaction": [],
        "hydrophobic_interaction": []
    }
    
    # PEPTIDE BONDS: consecutive residues in the same chain
    for i in range(len(all_residues) - 1):
        chain_id1, r1 = all_residues[i]
        chain_id2, r2 = all_residues[i+1]
        if chain_id1 == chain_id2 and is_consecutive(r1, r2):
            label1 = residue_label(r1)
            label2 = residue_label(r2)
            edges_by_category["peptide"].append((label1, label2))
    
    # DISULFIDE BONDS: from SSBOND lines
    ssbonds = parse_disulfide_bonds(pdb_text)
    for (c1, n1), (c2, n2) in ssbonds:
        rA = None
        rB = None
        for (chain_id, res) in all_residues:
            if chain_id == c1 and res.get_id()[1] == n1 and res.get_resname() == "CYS":
                rA = res
            if chain_id == c2 and res.get_id()[1] == n2 and res.get_resname() == "CYS":
                rB = res
        if rA and rB:
            edges_by_category["disulfide"].append((residue_label(rA), residue_label(rB)))
    
    # For other interactions, do an O(N^2) check of residue pairs
    for i in range(len(all_residues)):
        _, r1 = all_residues[i]
        label1 = residue_label(r1)
        for j in range(i+1, len(all_residues)):
            _, r2 = all_residues[j]
            label2 = residue_label(r2)
            
            # HYDROGEN BONDS
            hbond_found = False
            for a1 in r1.get_atoms():
                for a2 in r2.get_atoms():
                    if is_hbond(a1, a2):
                        edges_by_category["hydrogen_bond"].append((label1, label2))
                        hbond_found = True
                        break
                if hbond_found:
                    break
            
            # IONIC
            if is_ionic(r1, r2):
                edges_by_category["ionic_interaction"].append((label1, label2))
            
            # HYDROPHOBIC
            if is_hydrophobic(r1, r2):
                edges_by_category["hydrophobic_interaction"].append((label1, label2))
    
    return edges_by_category

###############################################################
# 5A. GATHER SEQUENCES (SINGLE-LETTER) FOR CHAINS
###############################################################

def get_chain_sequences(structure):
    """
    Returns a dict: { chain_id: single_letter_sequence } 
    for each chain in the first model,
    plus a concatenated sequence of all chains.
    """
    three_to_one = IUPACData.protein_letters_3to1  # e.g., {'Ala': 'A', 'Arg': 'R', ...}
    
    chain_seq_map = {}
    model = structure[0]
    
    for chain in model:
        seq_str = ""
        for residue in chain.get_residues():
            if residue.id[0] == " ":  # standard residue
                # Convert 3-letter to 1-letter
                rname = residue.get_resname().strip().title()  # e.g. "Ala"
                one_letter = three_to_one.get(rname, 'X')  # fallback to 'X'
                seq_str += one_letter
        chain_seq_map[chain.id] = seq_str
    
    # Concatenate sequences in chain order: A, B, C...
    all_chain_ids = sorted(chain_seq_map.keys())
    concatenated_seq = "".join(chain_seq_map[ch] for ch in all_chain_ids)
    
    return chain_seq_map, concatenated_seq

###############################################################
# 5B. WRITE OUTPUT FILE (TXT) WITH TWO SECTIONS
###############################################################

def write_txt_file(pdb_id, edges_by_category):
    """
    Writes a text file: <pdb_id>_interactions.txt with two main sections:
    1) AMINO ACID CONNECTIONS (resname + index)
    2) RESIDUE INDEX CONNECTIONS (index only)
    
    Each section lists PEPTIDE, DISULFIDE, HYDROGEN_BOND, etc. plus ALL INTERACTIONS.
    """
    txt_filename = f"{pdb_id}_interactions.txt"
    
    # We'll store ALL edges in a single list for each approach
    # so we can make an "ALL INTERACTIONS" at the end.
    all_edges_amino = []
    
    # Convert our edges into the two labeling schemes:
    # 1) Amino acid + index => already in edges_by_category
    # 2) Index only => parse out numeric part from each label
    
    # Gather edges in a new dict for numeric-only
    edges_by_category_numeric = {}
    for category, edges in edges_by_category.items():
        numeric_pairs = []
        for (labelA, labelB) in edges:
            # Extract just the integer part from each label, e.g. "MET15" -> "15"
            A_num = "".join([ch for ch in labelA if ch.isdigit()])
            B_num = "".join([ch for ch in labelB if ch.isdigit()])
            numeric_pairs.append((A_num, B_num))
        edges_by_category_numeric[category] = numeric_pairs
    
    with open(txt_filename, "w") as f:
        # ===== Section 1: AMINO ACID CONNECTIONS =====
        f.write("AMINO ACID CONNECTIONS:\n")
        
        for category, edges in edges_by_category.items():
            # Build a string like {MET1->ALA2,ALA2->VAL3}
            connections = []
            for (a, b) in edges:
                connections.append(f"{a}->{b}")
                all_edges_amino.append((a, b))  # track for ALL INTERACTIONS
            line = "{" + ",".join(connections) + "}"
            
            f.write(f"{category.upper().replace('_',' ')}:\n")
            f.write(line + "\n\n")
        
        # ALL INTERACTIONS in the amino acid labeling
        all_connections_amino = []
        for (a, b) in all_edges_amino:
            all_connections_amino.append(f"{a}->{b}")
        all_line_amino = "{" + ",".join(all_connections_amino) + "}"
        
        f.write("ALL INTERACTIONS:\n")
        f.write(all_line_amino + "\n\n")
        
        # ===== Section 2: RESIDUE INDEX CONNECTIONS =====
        f.write("\nRESIDUE INDEX CONNECTIONS:\n")
        
        # We'll keep track for a numeric ALL INTERACTIONS
        all_edges_numeric = []
        
        for category, edges in edges_by_category_numeric.items():
            # Build a string like {1->2,2->3}
            connections_numeric = []
            for (numA, numB) in edges:
                connections_numeric.append(f"{numA}->{numB}")
                all_edges_numeric.append((numA, numB))
            line_numeric = "{" + ",".join(connections_numeric) + "}"
            
            f.write(f"{category.upper().replace('_',' ')}:\n")
            f.write(line_numeric + "\n\n")
        
        # ALL INTERACTIONS in numeric labeling
        all_connections_numeric = []
        for (x, y) in all_edges_numeric:
            all_connections_numeric.append(f"{x}->{y}")
        all_line_numeric = "{" + ",".join(all_connections_numeric) + "}"
        
        f.write("ALL INTERACTIONS:\n")
        f.write(all_line_numeric + "\n\n")
    
    return txt_filename

###############################################################
# 5C. WRITE HTML FILE WITH 3D VIEW + INFO
###############################################################

def write_html_file(pdb_id, pdb_text, chain_seq_map, concatenated_seq):
    """
    Exports an HTML file named <pdb_id>_view.html with:
      - The PDB ID at the top
      - Link to RCSB PDB
      - Single-letter sequences (concatenated + chain-by-chain)
      - An interactive 3D model (py3Dmol)
    """
    html_filename = f"{pdb_id}_view.html"
    
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    
    # 3Dmol JS
    startjs = view.startjs
    endjs = view.endjs
    
    # Simple HTML with link to RCSB, sequences, 3D viewer
    rcsb_link = f"https://www.rcsb.org/structure/{pdb_id}"
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Protein {pdb_id} 3D View</title>
</head>
<body>
    <h1>PDB ID: {pdb_id}</h1>
    <p><b>RCSB Link:</b> 
       <a href="{rcsb_link}" target="_blank">{rcsb_link}</a>
    </p>
    <p>
      <b>Amino Acid Sequence (Concatenated):</b><br/>
      {concatenated_seq}
    </p>
"""
    # Chain-by-chain
    for chain_id in sorted(chain_seq_map.keys()):
        seq = chain_seq_map[chain_id]
        html_content += f"""    <p><b>Chain {chain_id}:</b><br/>
        {seq}
    </p>
"""
    # 3D viewer
    html_content += f"""
    <div id="3dmolviewer">
    {startjs}
    {endjs}
    </div>
</body>
</html>
"""
    
    with open(html_filename, "w") as f:
        f.write(html_content)
    
    return html_filename

###############################################################
# 6. MAIN
###############################################################

def main():
    print("=== Protein Interaction Graph (RCSB) ===")
    pdb_id = input("Enter a PDB ID (e.g. 1CRN or 4HHB): ").strip()
    
    # 1) Fetch PDB
    pdb_data = fetch_rcsb_pdb(pdb_id)
    if not pdb_data:
        return  # Could not fetch the file, so exit.
    
    # 2) Parse structure
    structure = parse_structure_from_pdb_text(pdb_data, structure_id=pdb_id)
    
    # 3) Build graph of interactions
    edges_by_category = build_full_graph(structure, pdb_data)
    
    # 4) Write .txt file (with both amino acid & numeric labeling)
    txt_file = write_txt_file(pdb_id, edges_by_category)
    print(f"[Info] Wrote interactions to '{txt_file}'")
    
    # 5) Gather sequences
    chain_seq_map, concatenated_seq = get_chain_sequences(structure)
    
    # 6) Create interactive HTML
    html_file = write_html_file(pdb_id, pdb_data, chain_seq_map, concatenated_seq)
    print(f"[Info] Wrote interactive 3D view to '{html_file}'")
    print(f"Open '{html_file}' in your browser to explore the protein.")

if __name__ == "__main__":
    main()