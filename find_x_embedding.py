import igraph as ig
import networkx as nx
import math
import re
import itertools
import numpy as np

comments_re = re.compile(r'#.*')
default_fixed_pos = [(0,0), (1,0), (0.5,1)]
U_INDEX = 0
V_INDEX = 1
W_INDEX = 2

def main(filename, fixed_vertices=None, plot_best=False):
    g = load_graph("data/" + filename, False)
    positions,_ = find_x_embedding(g)
    
    if plot_best:
        plot_x_embedding(g, positions)

def find_x_embedding(g, fixed_vertices=None):
    if fixed_vertices is None:
        edge_list = g.get_edgelist()
        nxg = nx.Graph(edge_list)
        full_cycle_list = list(nx.cycle_basis(nxg))
        tri_cycle_list = [fixed_vertices for fixed_vertices in full_cycle_list if len(fixed_vertices) == 3]
        # For testing ------------------
        #tri_cycle_list = []
        #full_cycle_list = [[0,1,2]]
        #-------------------------------
        if tri_cycle_list:
            fixed_vertex_details = []
            for fixed_vertices in tri_cycle_list:
                for comb in itertools.permutations(fixed_vertices,3):
                    fixed_vertex_details.append(comb)
            triangle_exists = True
        else:
            # TODO: Will update to find non-separating cycle
            # This might be just one cycle
            non_separating_cycles = full_cycle_list
            triangle_exists = False
            
            fixed_vertex_details = get_details_for_no_triangle(g, non_separating_cycles[0])
            for i in range(1, len(full_cycle_list)):
                fixed_vertex_details.extend(get_details_for_no_triangle(g, non_separating_cycles[i]))
        
        # Get positions and plot the graph with the points that are farthest
        # from each other
        max_min_dist = 0
        best_cycle = None
        for i in range(len(fixed_vertex_details)):            
            if triangle_exists:
                fixed_vertices = fixed_vertex_details[i]
                print("Fixed vertices: ", fixed_vertices)
                positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos)
            else:
                fixed_vertices = fixed_vertex_details[i][0]
                # Fixed for now, will improve upon once non-separating cycles work
                V1 = fixed_vertex_details[i][1]
                V2 = fixed_vertex_details[i][2]
                print("Fixed vertices: ", fixed_vertices)
                
                positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos, V1, V2, 1)
                
                # Note, there is a max gravity mentioned in the paper, but it's huge and tends to just pull everything
                # straight to one corner
                if not points_are_in_area(positions, fixed_vertices, V1, V2):
                    gravity = 1
                    attempts = 0
                    while not points_are_in_area(positions, fixed_vertices, V1, V2) and attempts < 30:
                        gravity += 0.5
                        positions = get_rubber_band_positions(g, fixed_vertices, default_fixed_pos, V1, V2, gravity)
                        attempts += 1
                    print("Attempts: ", attempts)
                    
            min_dist = get_min_dist_between_v(positions)
            if (min_dist > max_min_dist):
                max_min_dist = min_dist
                max_min_dist_pos = positions
                best_cycle = fixed_vertices
            print("Minimum distance between points: {0}\n".format(min_dist))
            
        print("Maximum min distance between points: {0} (fixed_vertices {1})".format(max_min_dist, best_cycle))
        # Flip so the y orientation is correct (with 0,0 at the lower left)
        print(list(max_min_dist_pos))
            
        layout = max_min_dist_pos
                   
    else:
        print("Fixed positions: {0}".format(fixed_vertices))
        layout = get_rubber_band_positions(g, fixed_vertices)
        best_cycle = fixed_vertices
        # Flip so the y orientation is correct (with 0,0 at the lower left)
        print("Minimum distance between points: {0}".format(get_min_dist_between_v(positions)))
        
    return layout, best_cycle
    
def plot_x_embedding(g, positions, filename=None):
    # Set visual style for graph
    visual_style={}
    visual_style['bbox']=(1200,1200)
    visual_style['margin']=50
    visual_style['vertex_label']=g.vs.indices
    visual_style['vertex_size']=15
    
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
    
    if (fixed_pos == None):
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
    subgraph_edges = []
    for e in g.es:
        v1, v2 = e.tuple
        if v1 in vertices and v2 in vertices:
            subgraph_edges.append(e.tuple)
    return subgraph_edges
        
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

def get_details_for_no_triangle(g, cycle):
    cycle_len = len(cycle)
    num_vertices = len(g.vs)
    vertices = g.vs.indices;
    
    details = []
    if cycle_len <= (num_vertices / 2):
        for i in range(cycle_len):
            uvw = [cycle[i], cycle[(i+1) % cycle_len], cycle[(i+2) % cycle_len]]
            V2 = [v for v in cycle if not v == uvw[W_INDEX]]
            V1 = [v for v in vertices if not v in V2]
            details.append([uvw, V1, V2])
            # Reverse the list
            uvw = [cycle[(i+2) % cycle_len], cycle[(i+1) % cycle_len], cycle[i]]
            V2 = [v for v in cycle if not v == uvw[W_INDEX]]
            V1 = [v for v in vertices if not v in V2]
            details.append([uvw, V1, V2])
    else:
        for v in vertices:
            # Had this in an if statement before: Which is neater, continue
            # or a large portion of this in an if statement?
            if v in cycle:
                continue
            uvw = []
            nbh = g.neighbors(v)
            nb_in_cycle = list(set(cycle) & set(nbh))
            nb_size = len(nb_in_cycle)
            # If the vertex has fewer than two neighbors, it can't be w
            # Same question as above about continue
            if nb_size < 2:
                continue
            w = v
            # We go through all the neighbors and pick possible v and u
            for i in range(nb_size):
                u = nb_in_cycle[i]
                v = nb_in_cycle[(i+1) % nb_size]
                u_index = cycle.index(u)
                v_index = cycle.index(v)
                path1 = [u]
                path2 = [v]
                # Get the shorter path between u and v (will definitely be less than |V|/2-1)
                for i in range(1,nb_size):
                    if (u_index + i) % nb_size == v_index:
                        path1.append(v)
                        V2 = path1
                        break
                    elif (v_index + i) % nb_size == u_index:
                        path2.append(u)
                        V2 = path2
                        break
                    path1.append(cycle[(u_index + i) % nb_size])
                    path2.append(cycle[(v_index + i) % nb_size])
                V1 = [vertex for vertex in vertices if not vertex in V2]
                details.append([[u, v, w], V1, V2])
                # Reverse u and v and add that as well, unless this
                # neighborhood is only two vertices, in which case 
                # it will reverse anyway
                if (nb_size > 2):
                    details.append([[v, u, w], V1, V2])
    
    print("Details:", details)
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

if __name__=="__main__":
    example_to_use = 1
    filenames = ["triangle_graph.txt", "power_3_connected_subset_0047_18.gml", "GMLFile.gml",
                 "GMLFile2.gml", "GMLFile3.gml", "GMLFile4.gml"]
    
    main(filenames[example_to_use], plot_best=True)