import math
import networkx as nx
import matplotlib.pyplot as plt

def interpolate(a, b, t):
    """Linear interpolation between scalars a and b, by fraction t."""
    return a + (b - a)*t

def make_arc(center_x, center_y, radius, start_deg, end_deg, n_points, label_prefix):
    """
    Approximate an arc (from start_deg to end_deg in degrees)
    with n_points along the curve. Returns a list of (node_name, (x, y)).
    """
    nodes = []
    for i in range(n_points+1):
        frac = i / n_points
        angle_deg = interpolate(start_deg, end_deg, frac)
        angle_rad = math.radians(angle_deg)
        x = center_x + radius*math.cos(angle_rad)
        y = center_y + radius*math.sin(angle_rad)
        node_name = f"{label_prefix}_arc{i}"
        nodes.append( (node_name, (x, y)) )
    return nodes

def make_square(center_x, center_y, half_size, label_prefix):
    """
    Makes a square with center (center_x, center_y) and
    'half_size' from center to each edge. Returns 4 corner nodes.
    """
    nodes = []
    corners = [
        (center_x - half_size, center_y + half_size),
        (center_x + half_size, center_y + half_size),
        (center_x + half_size, center_y - half_size),
        (center_x - half_size, center_y - half_size),
    ]
    for i, (x, y) in enumerate(corners):
        node_name = f"{label_prefix}_corner{i}"
        nodes.append( (node_name, (x, y)) )
    return nodes

def make_cross_brace(x1, y1, x2, y2, label_prefix, n_points=2):
    """
    Makes a simple 'line' or 'brace' from (x1,y1) to (x2,y2),
    optionally subdivided into n_points segments.
    """
    nodes = []
    for i in range(n_points+1):
        t = i / n_points
        x = interpolate(x1, x2, t)
        y = interpolate(y1, y2, t)
        node_name = f"{label_prefix}_brace{i}"
        nodes.append( (node_name, (x, y)) )
    return nodes

# ---------------------------------------------------------------------
# BUILD THE 'TOWER KIT' GRAPH
# ---------------------------------------------------------------------

G = nx.Graph()
pos = {}  # Will map node_name -> (x, y) for plotting

# 1) Example arcs for the tower legs (like the base arches).
#    Suppose we want two arcs forming one 'leg' shape, joined at endpoints.
leg_arc_1 = make_arc(center_x=0, center_y=0, radius=10,
                     start_deg=180, end_deg=90,
                     n_points=6, label_prefix="legA_upper")
leg_arc_2 = make_arc(center_x=0, center_y=0, radius=10,
                     start_deg=90, end_deg=0,
                     n_points=6, label_prefix="legA_lower")

# Add them to the graph & connect consecutive points
previous_node = None
for (n, (x,y)) in leg_arc_1:
    G.add_node(n)
    pos[n] = (x, y)
    if previous_node:
        G.add_edge(previous_node, n)
    previous_node = n

# Connect arc_1’s end to arc_2’s start
G.add_node(leg_arc_2[0][0])
pos[leg_arc_2[0][0]] = leg_arc_2[0][1]
G.add_edge(previous_node, leg_arc_2[0][0])

previous_node = leg_arc_2[0][0]
for (n, (x,y)) in leg_arc_2[1:]:
    G.add_node(n)
    pos[n] = (x, y)
    G.add_edge(previous_node, n)
    previous_node = n

# 2) Example “platform square” at some location
platform_nodes = make_square(center_x=25, center_y=10,
                             half_size=5, label_prefix="platform1")
for i, (n, (x,y)) in enumerate(platform_nodes):
    G.add_node(n)
    pos[n] = (x, y)
# Connect the square corners in a loop
for i in range(len(platform_nodes)):
    n1 = platform_nodes[i][0]
    n2 = platform_nodes[(i+1) % len(platform_nodes)][0]
    G.add_edge(n1, n2)

# 3) Example cross-brace
brace_nodes = make_cross_brace(x1=20, y1=0, x2=30, y2=20,
                               label_prefix="brace1", n_points=3)
prev_n = None
for (n, (x,y)) in brace_nodes:
    G.add_node(n)
    pos[n] = (x, y)
    if prev_n:
        G.add_edge(prev_n, n)
    prev_n = n

# 4) You could create multiple arcs, squares, braces, and so on,
#    just like your laser-cut layout.

# ---------------------------------------------------------------------
# VISUALIZE
# ---------------------------------------------------------------------
plt.figure(figsize=(8,8))

# Draw the graph in the “pos” coordinates we assigned
nx.draw_networkx(
    G, 
    pos, 
    with_labels=True, 
    node_size=200, 
    font_size=8
)

plt.gca().set_aspect("equal", adjustable="datalim")
plt.title("Eiffel-Tower-Like Parts as a Graph of Connection Points")
plt.axis("off")
plt.show()

# Optionally, print some adjacency info:
print("\nGraph summary:")
print(f"  Nodes: {len(G.nodes())}")
print(f"  Edges: {len(G.edges())}")