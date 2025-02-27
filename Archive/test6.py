import requests
import io
import networkx as nx
from Bio import PDB
import py3Dmol

###############################################################
# 1. FETCH PDB FROM ALPHAFOLD
###############################################################

def fetch_alphafold_pdb(uniprot_id):
    """
    Fetch the AlphaFold-predicted PDB for the given UniProt ID.
    Returns the PDB text or None if not found.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"[Error] Could not find PDB for UniProt ID '{uniprot_id}' on AlphaFold.")
        return None

###############################################################
# 2. PARSE PDB WITH BIOPYTHON
###############################################################

def parse_structure_from_pdb_text(pdb_text, structure_id="AF_model"):
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

POS_CHARGED = {"ARG", "LYS", "HIS"}  # Simplified
NEG_CHARGED = {"ASP", "GLU"}
HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO"}

def residue_label(res):
    """
    Returns a label like ALA12 for a Biopython Residue object.
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
    """
    ssbonds = []
    for line in pdb_text.splitlines():
        if line.startswith("SSBOND"):
            # e.g.: "SSBOND   1 CYS A    6    CYS A  127"
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
    and basic element check (N or O). This is very simplified!
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
    
    # Check side-chain distances
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
      edges_by_category["peptide"] = [(ALA12, ARG13), ...]
      edges_by_category["disulfide"] = [(CYS6, CYS127), ...]
      ...
    """
    model = structure[0]
    
    # Collect all standard residues (ignoring water, hetero)
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
    
    # PEPTIDE BONDS: consecutive residues in same chain
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
    
    # For other interactions, we do an O(N^2) check of residue pairs
    for i in range(len(all_residues)):
        _, r1 = all_residues[i]
        label1 = residue_label(r1)
        for j in range(i+1, len(all_residues)):
            _, r2 = all_residues[j]
            label2 = residue_label(r2)
            
            # HYDROGEN BONDS (check each atom pair)
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
# 5. WRITE OUTPUT FILES (TXT + HTML)
###############################################################

def write_txt_file(uniprot_id, edges_by_category):
    """
    Writes a text file named <uniprot_id>_interactions.txt grouping edges by interaction type.
    """
    txt_filename = f"{uniprot_id}_interactions.txt"
    with open(txt_filename, "w") as f:
        for category, edges in edges_by_category.items():
            # Convert edges into a single line inside { ... } with commas
            connections = []
            for (a, b) in edges:
                connections.append(f"{a}->{b}")
            
            line = "{" + ",".join(connections) + "}"
            
            f.write(f"{category.upper().replace('_',' ')}:\n")
            f.write(line + "\n\n")  # extra newline for clarity
    return txt_filename

def write_html_file(uniprot_id, pdb_text):
    """
    Exports an HTML file named <uniprot_id>_view.html with an interactive 3D model using py3Dmol.
    """
    html_filename = f"{uniprot_id}_view.html"
    
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    
    # Generate standalone HTML
    html_content = view._make_html()
    
    with open(html_filename, "w") as f:
        f.write(html_content)
    
    return html_filename

###############################################################
# 6. MAIN
###############################################################

def main():
    print("=== Protein Interaction Graph ===")
    uniprot_id = input("Enter a UniProt ID (e.g. P69905 for Hemoglobin subunit alpha): ").strip()
    
    # 1) Fetch PDB
    pdb_data = fetch_alphafold_pdb(uniprot_id)
    if not pdb_data:
        return
    
    # 2) Parse Structure
    structure = parse_structure_from_pdb_text(pdb_data, structure_id=uniprot_id)
    
    # 3) Build Graph of Interactions
    edges_by_category = build_full_graph(structure, pdb_data)
    
    # 4) Write .txt with grouped interactions
    txt_file = write_txt_file(uniprot_id, edges_by_category)
    print(f"[Info] Wrote interactions to '{txt_file}'")
    
    # 5) Create interactive HTML with py3Dmol
    html_file = write_html_file(uniprot_id, pdb_data)
    print(f"[Info] Wrote interactive 3D view to '{html_file}'")
    print(f"Open '{html_file}' in your browser to explore the protein.")

if __name__ == "__main__":
    main()