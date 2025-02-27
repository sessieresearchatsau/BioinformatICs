import requests
import io
import networkx as nx
from Bio import PDB
import py3Dmol
from IPython.display import display

###############################################################
# 1. FETCH PDB FROM ALPHAFOLD
###############################################################

def fetch_alphafold_pdb(uniprot_id):
    """
    Fetch the AlphaFold-predicted PDB for the given UniProt ID.
    Returns the PDB text as a string or None if not found.
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
    Uses Biopython's PDBParser to parse the structure from raw PDB text.
    Returns a Biopython Structure object.
    """
    # We can parse from an in-memory string by wrapping in StringIO
    pdb_io = io.StringIO(pdb_text)
    parser = PDB.PDBParser(QUIET=True)  # QUIET suppresses warnings
    structure = parser.get_structure(structure_id, pdb_io)
    return structure


###############################################################
# 3. IDENTIFY INTERACTIONS
###############################################################

# 3.1 HELPER: Identify residue categories
POS_CHARGED = {"ARG", "LYS", "HIS"}  # Histidine can be partially charged, but let's keep it simple
NEG_CHARGED = {"ASP", "GLU"}
HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO"}  
# (GLY and TYR can also be considered partially hydrophobic, but let's keep it simpler)

def residue_label(res):
    """
    Creates a consistent label for a Biopython Residue object (e.g. "ARG57").
    """
    resname = res.get_resname().strip()  # e.g. "ARG"
    # Biopython residue ID tuple looks like: (hetero flag, 'A', resseq, iCode)
    # res.get_id() -> (' ', 57, ' ')
    # Let's ignore chain letter for simplicity, or you can include it if you want.
    resseq = res.get_id()[1]
    return f"{resname}{resseq}"

def is_consecutive(r1, r2):
    """
    Returns True if r2 is exactly the next residue of r1 in terms of sequence number.
    """
    return (r2.get_id()[1] == r1.get_id()[1] + 1)

# 3.2 Parse SSBOND lines for disulfides (AlphaFold may not have them, but standard PDBs might)
def parse_disulfide_bonds(pdb_text):
    """
    Returns a list of ((chain1, resnum1), (chain2, resnum2)) for each SSBOND record in the PDB text.
    """
    ssbonds = []
    for line in pdb_text.splitlines():
        if line.startswith("SSBOND"):
            # Per PDB format: 
            # Columns: SSBOND  ****  CYS chain resnum  ****  CYS chain resnum
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

# 3.3 Detect hydrogen bonds (naïve distance-based; no angle check)
#    We check for donor/acceptor pairs (N/O) within 3.5 Å
def is_hbond(atom1, atom2, max_dist=3.5):
    """
    Naïve check if there's a potential H-bond between atom1 and atom2.
    We'll assume any N/O pair under 3.5Å is an H-bond (very rough).
    """
    e1 = atom1.element.upper()
    e2 = atom2.element.upper()
    # Donor ~ N, Acceptors ~ O (this is simplified, ignoring S, etc.)
    if (e1 in ["N","O"] and e2 in ["N","O"]):
        distance = atom1 - atom2  # Biopython distance
        if distance <= max_dist:
            return True
    return False

# 3.4 Detect ionic interactions (pos vs neg charged side chains within threshold)
def is_ionic(res1, res2, max_dist=5.0):
    """
    Check if res1 is positive and res2 is negative (or vice versa),
    and if any side-chain atoms are within `max_dist`.
    """
    resname1 = res1.get_resname()
    resname2 = res2.get_resname()
    
    is_res1_pos = resname1 in POS_CHARGED
    is_res2_pos = resname2 in POS_CHARGED
    is_res1_neg = resname1 in NEG_CHARGED
    is_res2_neg = resname2 in NEG_CHARGED
    
    # Must be pos <-> neg
    if not ((is_res1_pos and is_res2_neg) or (is_res1_neg and is_res2_pos)):
        return False
    
    # Check distance between side-chain atoms only
    # We'll exclude the backbone atoms: N, CA, C, O
    backbone_names = {"N", "CA", "C", "O"}
    
    for atom1 in res1.get_atoms():
        if atom1.get_name() not in backbone_names:
            for atom2 in res2.get_atoms():
                if atom2.get_name() not in backbone_names:
                    dist = atom1 - atom2
                    if dist <= max_dist:
                        return True
    return False

# 3.5 Detect hydrophobic interactions (any pair of hydrophobic side-chain atoms within threshold)
def is_hydrophobic(res1, res2, max_dist=4.5):
    """
    Returns True if both residues are hydrophobic and have side-chain atoms within 4.5 Å.
    """
    if (res1.get_resname() not in HYDROPHOBIC) or (res2.get_resname() not in HYDROPHOBIC):
        return False
    
    backbone_names = {"N", "CA", "C", "O"}
    for atom1 in res1.get_atoms():
        if atom1.get_name() not in backbone_names:
            for atom2 in res2.get_atoms():
                if atom2.get_name() not in backbone_names:
                    dist = atom1 - atom2
                    if dist <= max_dist:
                        return True
    return False

###############################################################
# 4. BUILD THE GRAPH USING NETWORKX
###############################################################
def build_full_graph(structure, pdb_text):
    """
    Builds a NetworkX graph where each residue is a node, 
    and edges exist for:
      - Peptide bonds (consecutive residues)
      - Disulfide bonds
      - Hydrogen bonds
      - Ionic interactions
      - Hydrophobic interactions
    
    Returns the graph (G) and an edge list for printing.
    """
    import networkx as nx
    G = nx.Graph()
    
    # 4.1 Parse out the model(s) - let’s just take model[0] for simplicity
    model = structure[0]
    
    # 4.2 Collect all residues across all chains
    all_residues = []
    for chain in model:
        for residue in chain.get_residues():
            # Exclude water or hetero-residues if you want
            if residue.id[0] == " ":  # ' ' = standard residue
                all_residues.append((chain.id, residue))
    
    # Sort by chain + residue number (optional, to keep them in order)
    all_residues.sort(key=lambda x: (x[0], x[1].get_id()[1]))
    
    # 4.3 Add nodes to the graph
    for chain_id, residue in all_residues:
        node_label = residue_label(residue)
        G.add_node(node_label)
    
    edges_list = []
    
    # 4.4 Peptide bonds: consecutive residues in the same chain
    for i in range(len(all_residues) - 1):
        chain_id1, r1 = all_residues[i]
        chain_id2, r2 = all_residues[i+1]
        if chain_id1 == chain_id2 and is_consecutive(r1, r2):
            n1 = residue_label(r1)
            n2 = residue_label(r2)
            G.add_edge(n1, n2, interaction="peptide_bond")
            edges_list.append(f"{{{n1} -> {n2}}} [peptide_bond]")
    
    # 4.5 Disulfide bonds (from SSBOND lines in PDB)
    ssbonds = parse_disulfide_bonds(pdb_text)
    if ssbonds:
        # We'll find the matching residue objects by chain+resnum
        for (c1, num1), (c2, num2) in ssbonds:
            # Locate the actual residue objects
            rA = None
            rB = None
            for chain_id, r in all_residues:
                if chain_id == c1 and r.get_id()[1] == num1 and r.get_resname() == "CYS":
                    rA = r
                if chain_id == c2 and r.get_id()[1] == num2 and r.get_resname() == "CYS":
                    rB = r
            if rA and rB:
                nA = residue_label(rA)
                nB = residue_label(rB)
                G.add_edge(nA, nB, interaction="disulfide_bond")
                edges_list.append(f"{{{nA} -> {nB}}} [disulfide_bond]")
    
    # 4.6 For each pair of residues, check hydrogen, ionic, hydrophobic
    #     (Naïve O(N^2) approach; can be slow for very large proteins)
    for i in range(len(all_residues)):
        chain_id1, r1 = all_residues[i]
        n1 = residue_label(r1)
        for j in range(i+1, len(all_residues)):
            chain_id2, r2 = all_residues[j]
            n2 = residue_label(r2)
            
            # Hydrogen bond?
            # We'll check every atom pair (N,O) distance
            # A more advanced approach uses angles and partial charges
            hbond_found = False
            for atom1 in r1.get_atoms():
                for atom2 in r2.get_atoms():
                    if is_hbond(atom1, atom2):
                        G.add_edge(n1, n2, interaction="hydrogen_bond")
                        edges_list.append(f"{{{n1} -> {n2}}} [hydrogen_bond]")
                        hbond_found = True
                        break
                if hbond_found:
                    break
            
            # Ionic?
            # Check if at least one pair of side-chain atoms is close and + vs -
            if is_ionic(r1, r2):
                G.add_edge(n1, n2, interaction="ionic_interaction")
                edges_list.append(f"{{{n1} -> {n2}}} [ionic_interaction]")
            
            # Hydrophobic?
            if is_hydrophobic(r1, r2):
                G.add_edge(n1, n2, interaction="hydrophobic_interaction")
                edges_list.append(f"{{{n1} -> {n2}}} [hydrophobic_interaction]")
    
    return G, edges_list


###############################################################
# 5. VISUALIZE PROTEIN IN 3D (py3Dmol)
###############################################################
def show_protein_3d(pdb_text):
    """
    Displays the protein 3D structure using py3Dmol (cartoon representation).
    """
    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_text, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    display(view.show())


###############################################################
# 6. MAIN RUN (IN A JUPYTER NOTEBOOK OR SIMILAR)
###############################################################
def main():
    print("=== Protein Graph & 3D Structure (All Interactions) ===")
    uniprot_id = input("Enter a UniProt ID (e.g. P69905 for Hemoglobin subunit alpha): ").strip()
    
    # 1) Fetch from AlphaFold
    pdb_data = fetch_alphafold_pdb(uniprot_id)
    if not pdb_data:
        return
    
    # 2) Parse with Biopython
    structure = parse_structure_from_pdb_text(pdb_data, structure_id=uniprot_id)
    
    # 3) Build full graph
    G, edges = build_full_graph(structure, pdb_data)
    
    # 4) Output 1: Print edges
    print("\nOutput 1: All connections (edges) in the protein graph")
    print("=======================================================")
    for e in edges:
        print(e)
    
    # 5) Output 2: Show 3D structure
    print("\nOutput 2: 3D Protein Structure (cartoon, colored by spectrum)")
    print("=============================================================")
    show_protein_3d(pdb_data)

# If running in a notebook, call main() in a cell. 
# If running in a standalone script, uncomment the following:
if __name__ == "__main__":
    main()