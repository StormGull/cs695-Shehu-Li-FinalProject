#!/usr/bin/python3
import igraph as ig
import networkx as nx
import math
import re
import itertools
import numpy as np
import time

comments_re = re.compile(r'#.*')
default_fixed_pos = [(0,0), (1,0), (0,1)]
V_INDEX = 0
U_INDEX = 1
W_INDEX = 2

def compute_positions(g, fixed_vertices, pull_to_vertices):
    if fixed_vertices is None:
        edge_list = g.get_edgelist()
        nxg = nx.Graph(edge_list)
        cycle_list = list(nx.cycle_basis(nxg))
        temp_tri_cycle_list = [cycle for cycle in cycle_list if len(cycle) == 3]
        tri_cycle_list = []
        for cycle in temp_tri_cycle_list:
            for comb in itertools.permutations(cycle,3):
                tri_cycle_list.append(comb)
        print(tri_cycle_list)
        
        # Get positions and plot the graph with the points that are farthest
        # from each other
        max_min_dist = 0
        best_cycle = None
        for cycle in tri_cycle_list:
            print("Cycle: ", cycle)
            # Fixed for now, will improve upon once non-separating cycles work        
            u = cycle[U_INDEX]
            v = cycle[V_INDEX]
            V2 = [u, v]
            V1 = [vertex for vertex in g.vs.indices if (not vertex == u and not vertex == v)]
            
            positions = get_rubber_band_positions(g, cycle, default_fixed_pos, V1, V2, 1)
            
            # Note, there is a max g mentioned in the paper, but it's huge and tends to just pull everything
            # straight to one corner
            if (pull_to_vertices and not points_are_in_area(positions, cycle, V1, V2)):
                gravity = 1
                attempts = 0
                while not points_are_in_area(positions, cycle, V1, V2) and attempts < 30:
                    gravity += 1
                    positions = get_rubber_band_positions(g, cycle, default_fixed_pos, V1, V2, gravity)
                    attempts += 1
                print("Attempts: ", attempts)
                    
            min_dist = get_min_dist_between_v(positions)
            if (min_dist > max_min_dist):
                max_min_dist = min_dist
                max_min_dist_pos = positions
                best_cycle = cycle
            print("Minimum distance between points: {0}\n".format(min_dist))
            
        print("Maximum min distance between points: {0} (cycle {1})".format(max_min_dist, best_cycle))
        # Flip so the y orientation is correct (with 0,0 at the lower left)
        layout = [ [p[0],-p[1]] for p in max_min_dist_pos ]
        print(list(max_min_dist_pos))
#        end = time.clock()
#        ig.plot(g, layout=layout, **visual_style)
                   
    else:
        positions = get_rubber_band_positions(g, fixed_vertices)
        layout = positions
        best_cycle = fixed_vertices
#        end = time.clock()
#        ig.plot(g, layout=positions, **visual_style)
        print("Fixed positions: {0}".format(fixed_vertices))
        print("Minimum distance between points: {0}".format(get_min_dist_between_v(positions)))

    return layout, best_cycle


def main(filename, fixed_vertices=None, pull_to_vertices=False):
    start = time.clock()
    g = load_graph("data/" + filename, False)
    
    # Set visual style for graph
    visual_style={}
    visual_style['bbox']=(1200,1200)
    visual_style['margin']=50
    visual_style['vertex_label']=g.vs.indices
    visual_style['vertex_size']=15
    
    layout, fixed_vertices = compute_positions(g, fixed_vertices, pull_to_vertices)
    
    end = time.clock()
    ig.plot(g, layout=layout, **visual_style)

    print("Time: ", end - start)

    return

def get_graph_matrices(g, fixed_vertices, fixed_pos_map, V1=[], V2=[], edge_gravity=1):
    vertices = list(g.vs.indices)
    unfixed_vertices = [vertex for vertex in vertices if vertex not in fixed_vertices]
    num_unfixed = len(unfixed_vertices)
    
    v_to_i = {k:v for v,k in enumerate(unfixed_vertices)}
    
    # Temporarily going to hard code this distinction
    V1_subgraph_edges = get_induced_subgraph_edges(g, V1)
    V2_subgraph_edges = get_induced_subgraph_edges(g, V2)
    
    c = {}
    for e in g.es:
        if e.tuple in V1_subgraph_edges or e.tuple in V2_subgraph_edges:
            c[e.tuple] = edge_gravity
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
        #print("Denominator: ", denominator)
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

def get_rubber_band_positions(g, fixed_vertices, fixed_pos=None, V1=[], V2=[], edge_gravity=1):
    
    if (fixed_pos == None):
        layout = list(g.layout_random())
        fixed_pos = []
        for i in range(len(fixed_vertices)):
            fixed_pos.append(tuple(layout[fixed_vertices[i]]))
        
    fixed_pos_map = dict(zip(fixed_vertices, fixed_pos))
    
    pos_matrix_x, pos_matrix_y, const_x, const_y = get_graph_matrices(g, fixed_vertices, fixed_pos_map, V1, V2, edge_gravity)

    ans_x = list(np.linalg.solve(pos_matrix_x,const_x))
    ans_y = list(np.linalg.solve(pos_matrix_y,const_y))
    
    for vertex in sorted(fixed_vertices):
        ans_x.insert(vertex,fixed_pos_map[vertex][0])
        ans_y.insert(vertex,fixed_pos_map[vertex][1])
    
    return list(zip(ans_x,ans_y))

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
    example_to_use = 4
    #filenames = ["square_graph.txt", "pentagon_graph.txt", "triangle_graph.txt", "GMLFile.gml", "GMLFile3.gml"]
    filenames = ["triangle_graph.txt", "GMLFile.gml", "GMLFile2.gml", "GMLFile3.gml", "GMLFile4.gml"]
    #fixed_vertices = [[0,1,2,3],[0,1,2,3,4],[3,4,5],[2,4,10]]
    
    main(filenames[example_to_use], pull_to_vertices=True)
