import igraph as ig
import networkx as nx
import math
import re
import itertools
import numpy as np
import random

comments_re = re.compile(r'#.*')
default_fixed_pos = [(0,0), (1,0), (0,1)]
U_INDEX = 0
V_INDEX = 1
W_INDEX = 2

def find_x_embedding_triconnected(g, fixed_vertices=None, shortcut=False):
    print("Finding X-Embedding...")
    
    if fixed_vertices is None:
        edge_list = g.get_edgelist()
        nxg = nx.Graph(edge_list)
        # Worst-case timing: O(V^3)
        base_cycle_list = list(nx.cycle_basis(nxg))
        tri_cycle_list = [fixed_vertices for fixed_vertices in base_cycle_list if len(fixed_vertices) == 3]
        # If there are triangles
        if tri_cycle_list:
            layout, best_cycle, min_dist = find_x_embedding_triangle(g, tri_cycle_list, shortcut)
        # Otherwise
        else:
            non_separating_cycle = get_non_separating_cycle(g, base_cycle_list)
            layout, best_cycle, min_dist = find_x_embedding_no_triangle(g, non_separating_cycle, shortcut)
    else:
        layout, best_cycle, min_dist = find_x_embedding_fixed(g, fixed_vertices)
        
    return layout, best_cycle, min_dist

def find_x_embedding_fixed(g, fixed_vertices):
    positions = get_rubber_band_positions(g, fixed_vertices)
    min_dist = get_min_dist_between_v(positions)
    print("Minimum distance between points: {0} (X={1})".format(min_dist, fixed_vertices))
    
    return positions, fixed_vertices, min_dist

def find_x_embedding_triangle(g, tri_cycle_list, shortcut=False):
    print("(Triangles found)")
    all_cycle_variations = []
    # If shortcut is selected, just pick a random triangle from the
    # cycle list.
    if shortcut:
        rand_cycle_index = int(random.random() * len(tri_cycle_list))
        all_cycle_variations = [tri_cycle_list[rand_cycle_index]]
    else:
        for fixed_vertices in tri_cycle_list:
            for comb in itertools.permutations(fixed_vertices,3):
                all_cycle_variations.append(comb)
    
    max_min_dist = -1
    best_cycle = None
    for fixed_vertices in all_cycle_variations:
        #print("Fixed vertices: ", fixed_vertices)
        positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos)
        
        min_dist = get_min_dist_between_v(positions)
        if min_dist > max_min_dist:
            max_min_dist = min_dist
            max_min_dist_pos = positions
            best_cycle = fixed_vertices
        
        #print("Minimum distance between points: {0}\n".format(min_dist))
    
    print("Maximum min distance between points: {0} (X={1})".format(max_min_dist, best_cycle))
        
    return max_min_dist_pos, best_cycle, max_min_dist

def find_x_embedding_no_triangle(g, non_separating_cycles, shortcut=False):
    print("(No triangles found)")
    if shortcut:
        rand_cycle_index = int(random.random() * len(non_separating_cycles))
        fixed_vertex_details = get_details_for_no_triangle(g, non_separating_cycles[rand_cycle_index], shortcut)
    else:
        fixed_vertex_details = get_details_for_no_triangle(g, non_separating_cycles[0])
        for i in range(1, len(non_separating_cycles)):
            fixed_vertex_details.extend(get_details_for_no_triangle(g, non_separating_cycles[i]))
    
    max_min_dist = -1
    best_cycle = None
    for i in range(len(fixed_vertex_details)):
        fixed_vertices = fixed_vertex_details[i][0]
        # Fixed for now, will improve upon once non-separating cycles work
        V1 = fixed_vertex_details[i][1]
        V2 = fixed_vertex_details[i][2]
        #print("Fixed vertices: ", fixed_vertices)
        #print("V1: ", V1)
        #print("V2: ", V2)
        
        gravity = 1
        attempts = 0
        # First phase - get to the general area
        positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos, V1, V2, gravity)
        while not points_are_in_area(positions, fixed_vertices, V1, V2):
            gravity *= 2
            positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos, V1, V2, gravity)
            attempts += 1
        
        # Second Phase - use a semi-binary search to locate the best location,
        # trying to get the lowest node in V1 close to 0.5
        min_dist_from_area = 0.01
        if get_lowest_in_V1(positions, V1) - 0.5 > min_dist_from_area:
            # Use as the base last gravity that put the vertices outside the desired area
            base_gravity = gravity / 2
            multiplier = 1.5
            mult_delta = 0.25
            while True:
                gravity = base_gravity * multiplier
                positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos, V1, V2, gravity)
                
                if points_are_in_area(positions, fixed_vertices, V1, V2):
                    if get_lowest_in_V1(positions, V1) - 0.5 <= min_dist_from_area:
                        break
                    else:
                        multiplier -= mult_delta
                else:
                    multiplier += mult_delta
                mult_delta /= 2
                    
        min_dist = get_min_dist_between_v(positions)
        if min_dist > max_min_dist:
            max_min_dist = min_dist
            max_min_dist_pos = positions
            best_cycle = fixed_vertices
        #print("Minimum distance between points: {0}\n".format(min_dist))
            
    
    print("Maximum min distance between points: {0} (X={1})".format(max_min_dist, best_cycle))
            
    return max_min_dist_pos, best_cycle, max_min_dist

def get_lowest_in_V1(positions, V1):
    lowest_loc = 1
    for vertex in V1:
        if positions[vertex][1] < lowest_loc:
            lowest_loc = positions[vertex][1]
    return lowest_loc

def plot_x_embedding(g, positions, filename=None):
    # Set visual style for graph
    visual_style={}
    visual_style['bbox']=(1200,1200)
    visual_style['margin']=50
    visual_style['vertex_label']=g.vs.indices
    visual_style['vertex_size']=20
    visual_style["vertex_color"]="blue" 
    visual_style["vertex_label_size"] = 15
    visual_style["vertex_label_color"] = "white"
    
    layout = [ [p[0],-p[1]] for p in positions]
    if filename:
        ig.plot(g, layout=layout, **visual_style).save(filename)
    else:
        ig.plot(g, layout=layout, **visual_style)

def get_graph_matrices(g, fixed_vertices, fixed_pos_map, V1=[], V2=[], gravity=1):
    vertices = list(g.vs.indices)
    unfixed_vertices = [vertex for vertex in vertices if vertex not in fixed_vertices]
    num_unfixed = len(unfixed_vertices)
    
    v_to_i = {k:v for v,k in enumerate(unfixed_vertices)}
    
    # Temporarily going to hard code this distinction
    if V1:
        V1_subgraph_edges = get_induced_subgraph_edges(g, V1)
        V2_subgraph_edges = get_induced_subgraph_edges(g, V2)
    else:
        V1_subgraph_edges = []
        V2_subgraph_edges = []
    
    c = {}
    for e in g.es:
        if e.tuple in V1_subgraph_edges or e.tuple in V2_subgraph_edges:
            c[e.tuple] = gravity
        else:
            c[e.tuple] = 1
    
    # We're going to solve for the following equation with matrices:
    # f_v = (1 / sum_for_each_neighbor_u(c_uv)) * 
    #            sum_for_each_neighbor_u(c_uv * f_u)
    #     for all unfixed vertices, where c_uv is the elasticity coefficient
    #     and f_v and f_u are vertex positions
    pos_matrix_x = np.zeros((num_unfixed,num_unfixed))
    pos_matrix_y = np.zeros((num_unfixed,num_unfixed))
    const_x = [0 for _ in range(num_unfixed)]
    const_y = [0 for _ in range(num_unfixed)]
    for vertex in unfixed_vertices:
        v_index = v_to_i[vertex]
        nb_list = g.neighbors(vertex)
        denominator = 0
        # In the equation, we want one side to be a constant, so we
        # subtract the f_v variable from both sides
        pos_matrix_x[v_index][v_index] = -1
        pos_matrix_y[v_index][v_index] = -1
        # Get the denominator for the equation - the sum of the elasticity coefficients
        # We divide everything by this except for f_v 
        for nb in nb_list:
            c_index = tuple(sorted([vertex,nb]))
            denominator += c[c_index]
        for nb in nb_list:
            c_index = tuple(sorted([vertex,nb]))
            # If the neighbor is a fixed vertex, then we already know it's position, so it's not
            # a variable we have to solve for. Therefore, it goes over to the constant side
            # (subtracting from both sides, so it's negative)
            if nb in fixed_vertices:
                const_x[v_index] += (-c[c_index]*fixed_pos_map[nb][0]) / denominator
                const_y[v_index] += (-c[c_index]*fixed_pos_map[nb][1]) / denominator
            # Otherwise, it's a variable we have to solve for. We're multiplying the variable
            # by its spring value, so that's what shows up in the matrix
            else:
                nb_index = v_to_i[nb]
                pos_matrix_x[v_index][nb_index] = c[c_index] / denominator
                pos_matrix_y[v_index][nb_index] = c[c_index] / denominator
                
    return pos_matrix_x, pos_matrix_y, np.array(const_x), np.array(const_y)

# Current fixed vertex order: u, v, w
def get_rubber_band_positions(g, fixed_vertices, fixed_pos=None, V1=[], V2=[], gravity=1):
    
    if fixed_pos and not len(fixed_pos) == len(fixed_vertices):
        print("ERROR: The length of the fixed positions must be equal to the length of the fixed vertices")
        print("Using automated positions")
        fixed_pos = None
    if not fixed_pos:
        layout = list(g.layout_auto())
        fixed_pos = []
        for i in range(len(fixed_vertices)):
            fixed_pos.append(tuple(layout[fixed_vertices[i]]))
        
    fixed_pos_map = dict(zip(fixed_vertices, fixed_pos))
    
    pos_matrix_x, pos_matrix_y, const_x, const_y = get_graph_matrices(g, fixed_vertices, fixed_pos_map, V1, V2, gravity)

    ans_x = list(np.linalg.solve(pos_matrix_x,const_x))
    ans_y = list(np.linalg.solve(pos_matrix_y,const_y))
    
    for vertex in sorted(fixed_vertices):
        ans_x.insert(vertex,fixed_pos_map[vertex][0])
        ans_y.insert(vertex,fixed_pos_map[vertex][1])
    
    return list(zip(ans_x,ans_y))

# Not using built-in induced subgraph since it doesn't preserve
# the vertex/edge indices, which we need
# Alternately, we could name the edges, then induce the subgraph
# and select the edges from the original graph
def get_induced_subgraph_edges(g, vertices):
    return [e.tuple for e in g.es if e.target in vertices and e.source in vertices]
            
def get_dist(pos1, pos2):
    x1, y1 = pos1
    x2, y2 = pos2    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def points_are_in_area(v_positions, fixed_vertices, V1, V2):
    in_area = True
    for i in range(len(v_positions)):  
        if i not in fixed_vertices:     
            x,y = v_positions[i]
            if i in V1 and y < 0.5:
                return False
            elif i in V2 and (y > x or x > 1-(2*y)):
                return False
    return in_area

def get_details_for_no_triangle(g, cycle, shortcut=False):
    cycle_len = len(cycle)
    num_vertices = len(g.vs)
    vertices = g.vs.indices;
    details = []
    # If the cycle is less than half the length of the graph, then we
    # can pick any consecutive u,v,w. This finds all possible
    # u, v, w simply by starting at each index for u and counting
    # forward from there. It also gets them in reverse.
    # V2 is just every vertex in the cycle except for W, and V1
    # is everything else
    # Note: The order looks weird because I always use the vertices in the order
    # u, v, w, but they need to be in consecutive order u, w, v
    if cycle_len <= (num_vertices / 2) + 1:
        for i in range(cycle_len):
            uvw = [cycle[i], cycle[(i+2) % cycle_len], cycle[(i+1) % cycle_len]]
            V2 = [v for v in cycle if not v == uvw[W_INDEX]]
            V1 = [v for v in vertices if not v in V2]
            details.append([uvw, V1, V2])
            if shortcut:
                break
    # Here, w is a vertex not in the cycle (C) that has two neighbors in C.
    # U and V are those neighbors. 
    else:
        # Find w
        for vertex in vertices:
            # We skip the vertices in the cycle since none of them can be w
            if vertex in cycle:
                continue
            uvw = []
            nbh = g.neighbors(vertex)
            # Get the neighbors of the selected vertex that are in the cycle
            nb_in_cycle = list(set(cycle) & set(nbh))
            nb_size = len(nb_in_cycle)
            # We are only interested in vertices with two or more neighbors in the cycle
            if nb_size < 2:
                continue
            w = vertex
            # Find all the possible v,u combinations and get V1 and V2
            # V2 is the path between u and v (must be < |V|/2 - 1). V1 is everything else
            vu = list(itertools.permutations(nb_in_cycle, 2))
            for v,u in vu:
                u_index = cycle.index(u)
                v_index = cycle.index(v)
                path1 = [u]
                path2 = [v]
                # Get the shorter path between u and v (will definitely be less than |V|/2-1)
                for i in range(1,cycle_len):
                    if (u_index + i) % cycle_len == v_index:
                        path1.append(v)
                        V2 = path1
                        break
                    elif (v_index + i) % cycle_len == u_index:
                        path2.append(u)
                        V2 = path2
                        break
                    path1.append(cycle[(u_index + i) % cycle_len])
                    path2.append(cycle[(v_index + i) % cycle_len])
                V1 = [vertex for vertex in vertices if not vertex in V2]
                details.append([[u, v, w], V1, V2])
                if shortcut:
                    break
            if shortcut:
                break

    return details

# Brute force method - tried to find another method but haven't gotten it working yet
# (the other way is divide and conquer and is faster)
def get_min_dist_between_v(v_positions):
    V = list(range(len(v_positions)))
    min_dist = float("inf")
    vertex_pairs = list(itertools.combinations(V,2))
    for pair in vertex_pairs:
        v1,v2 = pair
        dist = get_dist(v_positions[v1], v_positions[v2])
        if dist < min_dist:
            min_dist = dist
            
    return min_dist

def all_vertices_neighbor_graph(g, cycle):
    all_vertices_connect = True
    non_cycle_vertices = [i for i in g.vs.indices if not i in cycle]
    for v in cycle:
        nbh = g.neighbors(v)
        nb_in_graph = [nb for nb in nbh if nb in non_cycle_vertices]
        if not nb_in_graph:
            all_vertices_connect = False
            break
    
    return all_vertices_connect

def cycle_has_chords(g, cycle):
    edges_in_cycle = [e for e in g.es if e.target in cycle and e.source in cycle]
    if not len(edges_in_cycle) == len(cycle):
        return True
    else:
        return False

def cycle_is_non_separating(g, cycle):
    g_copy = g.copy()
    g_copy.delete_vertices(cycle)
    if g_copy.is_connected():
        return True
    else:
        return False

def cycle_is_acceptable(g, cycle, min_length=4):
    if (len(cycle) >= min_length and len(cycle) < len(g.es) 
        and all_vertices_neighbor_graph(g, cycle)
        and not cycle_has_chords(g, cycle)
        and cycle_is_non_separating(g, cycle)):
        
        return True
    else:
        return False

def get_non_separating_cycle(g, cycle_basis):
    count = 0
    # Start by checking the cycle basis, as
    # it is much faster to check than
    # going through every possible cycle
    non_separating_cycles = []
    for cycle in cycle_basis:
        if (cycle_is_acceptable(g, cycle)):
            non_separating_cycles.append(cycle)
        count += 1
    
    # If acceptable cycles are found in the cycle basis,
    # search through all the cycles and pick the first
    # cycle found. This can take a significant amount of time
    # depending on the size of the graph
    if not non_separating_cycles:
        print("No non-separating cycles in cycle basis")
        g_directed = g.copy()
        g_directed.to_directed(mutual=True)
        edge_list = g_directed.get_edgelist()
        nxg = nx.DiGraph(edge_list)
        count = 0
        for cycle in nx.simple_cycles(nxg):
            if (count % 100000 == 0):
                print("Examined {0} cycles".format(count))
            if (cycle_is_acceptable(g, cycle)):
                non_separating_cycles.append[cycle]
                break
            count += 1
    
    #if non_separating_cycles:
    #    print("Found cycles: ", non_separating_cycles)
    #else:
    #    print("Did not find cycle.")
    
    return non_separating_cycles

def load_graph(data_file_name, directed=False):
    if (data_file_name.endswith(".gml")):
        return ig.Graph.Read_GML(data_file_name)
    else:
        V = set()
        E = set()
    
        with open(data_file_name) as data_file:
            for line in data_file:
                line = comments_re.sub('', line.strip()).strip()
                if line:
                    source, target = line.split()
                    source = int(source)
                    target = int(target)
                    V.add(source)
                    V.add(target)
                    E.add((source, target))
        
        g = ig.Graph(directed=directed)
        g.add_vertices(sorted(list(V)))
        g.add_edges(list(E))
        return g
    
def main(filename, fixed_vertices=None, plot_best=False, shortcut=False):
    g = load_graph("data/" + filename, False)
    #g = ig.Graph.Full_Bipartite(6, 6)
    positions,_,_ = find_x_embedding_triconnected(g, shortcut=shortcut)
    
    if plot_best:
        plot_x_embedding(g, positions)

if __name__=="__main__":
    example_to_use = 0
    filenames = ["square_graph.txt", "power_3_connected_subset_0047_18.gml", "32-node-random-triconnected.gml",
                 "GMLFile2.gml", "GMLFile3.gml", "GMLFile4.gml", "ring_graph.txt"]
    
    main(filenames[example_to_use], plot_best=True, shortcut=False)