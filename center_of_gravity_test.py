import igraph as ig
import math
import re

comments_re = re.compile(r'#.*')

def main(filename, fixed_vertices):
     
    g = load_tsv_edges("data/" + filename, False)
    V = list(g.vs.indices)    
    positions = list(g.layout_kamada_kawai())
    loose_vertices = [v for v in V if v not in fixed_vertices]
    
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
    
    print("Number of iterations: {0}".format(num_iterations))

    # Display graphs
    visual_style={}
    visual_style['bbox']=(600,600)
    visual_style['margin']=50
    visual_style['vertex_label']=g.vs.indices

    ig.plot(g, layout=positions, **visual_style)

    return

def get_dist(pos1, pos2):
    x1, y1 = pos1
    x2, y2 = pos2    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def load_tsv_edges(data_file_name, directed=False):

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
    filenames = ["square_graph.txt", "pentagon_graph.txt"]
    fixed_vertices = [[0,1,2,3],[0,1,2,3,4]]
    
    main(filenames[example_to_use],fixed_vertices[example_to_use])