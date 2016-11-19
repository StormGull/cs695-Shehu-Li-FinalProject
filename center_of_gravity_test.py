import igraph as ig
import networkx as nx
import math
import re
import itertools

comments_re = re.compile(r'#.*')

def main(filename, fixed_vertices=None):
     
    g = load_graph("data/" + filename, False)
    
    # Find all cycles of length 3
    if fixed_vertices is None:
        edge_list = g.get_edgelist()
        nxg = nx.Graph(edge_list)
        cycle_list = list(nx.cycle_basis(nxg))
        tri_cycle_list = [cycle for cycle in cycle_list if len(cycle) == 3]

    # Set visual style for graph
    visual_style={}
    visual_style['bbox']=(1200,1200)
    visual_style['margin']=50
    #visual_style['vertex_label']=g.vs.indices
    visual_style['vertex_size']=5
    
    # Get positions and plot graphs
    if fixed_vertices is None:
        max_min_dist = 0
        best_cycle = None
        for cycle in tri_cycle_list:
            #print("Fixed positions: {0}".format(cycle))
            positions = get_rubber_band_positions(g, cycle, True)
            #ig.plot(g, layout=positions, **visual_style)
            min_dist = get_min_dist_between_v(positions)
            if (min_dist > max_min_dist):
                max_min_dist = min_dist
                max_min_dist_pos = positions
                best_cycle = cycle
            #print("Minimum distance between points: {0}".format(min_dist))
        print("Maximum min distance between points: {0} (cycle {1})".format(max_min_dist, best_cycle))
        ig.plot(g, layout=max_min_dist_pos, **visual_style)
                   
    else:
        positions = get_rubber_band_positions(g, fixed_vertices)
        #ig.plot(g, layout=positions, **visual_style)
        print("Fixed positions: {0}".format(fixed_vertices))
        print("Minimum distance between points: {0}".format(get_min_dist_between_v(positions)))

    return

def get_rubber_band_positions(g, fixed_vertices, right_triangle_layout=False):
    V = list(g.vs.indices)
    positions = list(g.layout_kamada_kawai())
    loose_vertices = [v for v in V if v not in fixed_vertices]
    
    if right_triangle_layout:
        if len(fixed_vertices) != 3:
            print("Must have three fixed vertices for the right triangle layout")
        else:
            positions[fixed_vertices[0]] = [0,0]
            positions[fixed_vertices[1]] = [0,1]
            positions[fixed_vertices[2]] = [1,0]
    
    d_stop = 0.0000005;
    max_distance_moved = d_stop + 1
    num_iterations = 0
    while max_distance_moved > d_stop:
        max_distance_moved = 0
        for v_index in loose_vertices:
            x_sum = 0
            y_sum = 0
            nh_list = g.neighbors(v_index)
            for node in nh_list:
                x_sum += positions[node][0]
                y_sum += positions[node][1]
            old_pos = [positions[v_index][0], positions[v_index][1]]
            new_pos = [x_sum / len(nh_list), y_sum / len(nh_list)]
            distance_diff = get_dist(old_pos, new_pos)
            if distance_diff > max_distance_moved:
                max_distance_moved = distance_diff
            positions[v_index] = new_pos
        num_iterations += 1
    
    #print("Number of iterations: {0}".format(num_iterations))
    return positions
    
def get_dist(pos1, pos2):
    x1, y1 = pos1
    x2, y2 = pos2    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

# Brute force method - tried to find another method but haven't gotten it working yet
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
    
    main(filenames[example_to_use])