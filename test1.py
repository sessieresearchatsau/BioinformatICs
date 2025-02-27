import networkx as nx
import matplotlib.pyplot as plt

# Example: Simple protein backbone graph with 5 residues
G = nx.Graph()
residues = ["Ala1", "Val2", "Gly3", "Ser4", "Leu5"]

# Add nodes
for res in residues:
    G.add_node(res)

# Add edges (simple linear backbone connectivity)
for i in range(len(residues) - 1):
    G.add_edge(residues[i], residues[i + 1])

# Draw graph
plt.figure(figsize=(5, 5))
nx.draw(G, with_labels=True, node_color="lightblue", edge_color="gray", node_size=2000, font_size=10)
plt.show()