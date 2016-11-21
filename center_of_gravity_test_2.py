import igraph as ig
import networkx as nx
import math
import re
import itertools

comments_re = re.compile(r'#.*')

def main(filename, fixed_vertices=None, pull_to_vertices=False):
     
    g = load_graph("data/" + filename, False)
    
    # Set visual style for graph
    visual_style={}
    visual_style['bbox']=(1200,1200)
    visual_style['margin']=50
    visual_style['vertex_label']=g.vs.indices
    visual_style['vertex_size']=15
        
    # Find all cycles of length 3
    if fixed_vertices is None:
        edge_list = g.get_edgelist()
        nxg = nx.Graph(edge_list)
        cycle_list = list(nx.cycle_basis(nxg))
        tri_cycle_list = [cycle for cycle in cycle_list if len(cycle) == 3]
    
        # Get positions and plot the graph with the points that are farthest
        # from each other
        max_min_dist = 0
        best_cycle = None
        for cycle in tri_cycle_list:
            print("Cycle: ", cycle)
            positions = get_rubber_band_positions(g, cycle, True, pull_to_vertices)
            #ig.plot(g, layout=positions, **visual_style)
            min_dist = get_min_dist_between_v(positions)
            if (min_dist > max_min_dist):
                max_min_dist = min_dist
                max_min_dist_pos = positions
                best_cycle = cycle
            print("Minimum distance between points: {0}\n".format(min_dist))
        print("Maximum min distance between points: {0} (cycle {1})".format(max_min_dist, best_cycle))
        # Flip so the y orientation is correct (with 0,0 at the lower left)
        layout = [ [p[0],-p[1]] for p in max_min_dist_pos ]
        ig.plot(g, layout=layout, **visual_style)
                   
    else:
        positions = get_rubber_band_positions(g, fixed_vertices, False)
        ig.plot(g, layout=positions, **visual_style)
        print("Fixed positions: {0}".format(fixed_vertices))
        print("Minimum distance between points: {0}".format(get_min_dist_between_v(positions)))

    return

def get_rubber_band_positions(g, fixed_vertices, right_triangle_layout=True, pull_to_vertices=False):
    vertices = list(g.vs.indices)
    positions = g.layout_random()
    loose_vertices = [vertex for vertex in vertices if vertex not in fixed_vertices]
    
    if right_triangle_layout:
        if len(fixed_vertices) != 3:
            print("Must have three fixed vertices for the right triangle layout")
        else:
            positions[fixed_vertices[0]] = [0,0]
            positions[fixed_vertices[1]] = [0,1]
            positions[fixed_vertices[2]] = [1,0]
            u = fixed_vertices[0]
            w = fixed_vertices[1]
            v = fixed_vertices[2]
    
    d_stop = 0.0000005;
    max_distance_moved = d_stop + 1
    num_iterations = 0
    settled = False
    num_attempts = 0
    vertex_gravity = 1
    first_attempt = True
    V1 = []
    V2 = []
    while not settled:
        max_distance_moved = 0
        for v_index in loose_vertices:
            x_sum = 0
            y_sum = 0
            nh_list = g.neighbors(v_index)
            denominator_multiplier = 0;
            for node in nh_list:
                # When two vertices are in the same group, they pull towards each other
                # This simulating the stronger springs inside the induced subgraphs (that is,
                # a graph that includes the subset of plus edges where both endpoints are
                # in the selected vertices) of V1 and V2
                if node in V2 and v_index in V2:
                    weight = vertex_gravity
                    denominator_multiplier += 1
                elif node in V1 and v_index in V1:
                    weight = vertex_gravity
                    denominator_multiplier += 1
                else:
                    weight = 1
                x_sum += positions[node][0] * weight
                y_sum += positions[node][1] * weight
            old_pos = [positions[v_index][0], positions[v_index][1]]
            # Get weighted denominator
            # The formula is essentially: cx = (v1x*m1 + v2x*m2 + ... vnx*mn) / (m1 + m2 .... mn)
            # (same thing for y)
            denominator = len(nh_list)+denominator_multiplier*(vertex_gravity-1)
            new_pos = [x_sum / denominator, y_sum / denominator]
            distance_diff = get_dist(old_pos, new_pos)
            if distance_diff > max_distance_moved:
                max_distance_moved = distance_diff
            positions[v_index] = new_pos
        num_iterations += 1
        settled = max_distance_moved < d_stop
        # If we're pulling into the vertices, check to see if the unfixed vertices
        # are within the correct area. If not, increase the gravity of the vertices
        # and let the loose vertices re-settle
        if pull_to_vertices and settled:
            positions_without_fixed = [p for v, p in enumerate(positions) if v not in fixed_vertices]
            if not points_are_in_area(positions_without_fixed):
                if first_attempt:
                    # This is a temporary hard-coding of V1 and V2. Later will make it so it pulls out all the
                    # nodes of the cycle that are not w (will need to give this or whatever function this goes
                    # into both the full cycle and the fixed vertices
                    V1 = [u,v]
                    V2 = [vertex for vertex in vertices if (not vertex == u and not vertex == v)]
                    first_attempt = False
                vertex_gravity += 0.5
                num_attempts += 1
                settled = False
    
    print("Number of iterations: {0}".format(num_iterations))
    print("Number of attempts: {0}".format(num_attempts)) # To settle in the correct locations
    return positions
    
def get_dist(pos1, pos2):
    x1, y1 = pos1
    x2, y2 = pos2    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def points_are_in_area(v_positions):
    in_area = True
    for v in v_positions:
        x,y = v
        if (y < 0.5 and (y > x or x > 1-(2*y))):
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