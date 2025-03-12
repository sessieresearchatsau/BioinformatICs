import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 1) Create an empty graph
G = nx.Graph()

# 2) Define 3D coordinates for each node
#    (Below are just approximate positions to mimic a tapering tower shape.)

pos_3d = {
    # Level 0 (base, z=0)
    "A0": ( 1.0,  1.0, 0),
    "B0": ( 1.0, -1.0, 0),
    "C0": (-1.0, -1.0, 0),
    "D0": (-1.0,  1.0, 0),
    
    # Level 1 (z=1)
    "A1": ( 0.8,  0.8, 1),
    "B1": ( 0.8, -0.8, 1),
    "C1": (-0.8, -0.8, 1),
    "D1": (-0.8,  0.8, 1),
    
    # Level 2 (z=2)
    "A2": ( 0.4,  0.4, 2),
    "B2": ( 0.4, -0.4, 2),
    "C2": (-0.4, -0.4, 2),
    "D2": (-0.4,  0.4, 2),
    
    # Top spire (z=3)
    "T":  ( 0.0,  0.0, 3)
}

# Add these nodes to our graph
G.add_nodes_from(pos_3d.keys())

# 3) Define the edges (same pattern as the 2D example)

# -- Horizontal (Perimeter) --
perimeter_level0 = [("A0","B0"), ("B0","C0"), ("C0","D0"), ("D0","A0")]
perimeter_level1 = [("A1","B1"), ("B1","C1"), ("C1","D1"), ("D1","A1")]
perimeter_level2 = [("A2","B2"), ("B2","C2"), ("C2","D2"), ("D2","A2")]

# -- Vertical --
vertical_0_1 = [("A0","A1"), ("B0","B1"), ("C0","C1"), ("D0","D1")]
vertical_1_2 = [("A1","A2"), ("B1","B2"), ("C1","C2"), ("D1","D2")]

# -- Diagonal --
diag_0_1 = [("A0","B1"), ("B0","C1"), ("C0","D1"), ("D0","A1")]
diag_1_2 = [("A1","B2"), ("B1","C2"), ("C1","D2"), ("D1","A2")]

# -- Top spire connections --
top_connections = [("A2","T"), ("B2","T"), ("C2","T"), ("D2","T")]

# Combine all edges
all_edges = (
    perimeter_level0
    + perimeter_level1
    + perimeter_level2
    + vertical_0_1
    + vertical_1_2
    + diag_0_1
    + diag_1_2
    + top_connections
)

G.add_edges_from(all_edges)

# 4) Print adjacency sets (optional, for inspection)
print("Adjacency Sets (Node -> Neighbors):")
for node in G.nodes():
    print(f"{node} -> {sorted(G.adj[node])}")

# 5) Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot each node as a point, and label it
for node in G.nodes():
    x, y, z = pos_3d[node]
    ax.scatter(x, y, z, s=40)  # s=40 sets a moderate node size
    ax.text(x, y, z, node, fontsize=9, zorder=10)

# Draw edges by connecting the node coordinates
for edge in G.edges():
    x_vals = [pos_3d[edge[0]][0], pos_3d[edge[1]][0]]
    y_vals = [pos_3d[edge[0]][1], pos_3d[edge[1]][1]]
    z_vals = [pos_3d[edge[0]][2], pos_3d[edge[1]][2]]
    ax.plot(x_vals, y_vals, z_vals)

ax.set_title("Simplified 3D Eiffel Tower Graph")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()