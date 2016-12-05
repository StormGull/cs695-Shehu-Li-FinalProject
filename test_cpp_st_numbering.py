#!/usr/bin/python3

import center_of_gravity_test_3 as cog
import json
import matplotlib.pyplot as plt
import os
import pdb
import random

def render_line(line, **kw):
    xs = [line[0][0], line[1][0]]
    ys = [line[0][1], line[1][1]]
    return plt.Line2D(xs, ys, **kw)


def assign_power(g):
    # Assign power supply 1, sink -1 randomly (but evenly)
    power = []
    for i, _ in enumerate(g.vs):
        power.append(i%2 and 1 or -1)
    random.shuffle(power)
    return power

def slope(line):
    return (line[1][1] - line[0][1])/(line[1][0] - line[0][0])

def split_network(solution):
    # order vertices by y position if slope > 1, otherwise
    # by x position
    st_numbering = solution['st_numbering']
    tangent_line = solution['tangent']
    tangent_points = [solution['pair_indices']['p1'], solution['pair_indices']['p2']]

    if abs(slope(tangent_line)) > 1:
        key_func = lambda x: x[1][1]
    else:
        key_func = lambda x: x[1][0]
    ordering = [(v, p) for v, p in enumerate(st_numbering)]

    try:
        ordering.sort(key=key_func)
    except Exception as e:
        pdb.set_trace()
        foo = 1
    
    # Add all vertices to one set or the other, except those used to
    # create the tangent line. Then add them later.

    V1, V2 = set(), set()
    found_tangent_vertices = False
    for v, _ in ordering:
        if found_tangent_vertices:
            if v not in tangent_points:
                V2.add(v)
        elif v in tangent_points:
            found_tangent_vertices = True
        else:
            V1.add(v)
    
    # The theory says we can always adjust the line that runs through
    # the two vertices that created the tangent in a way that puts
    # reverses which community they are in. So just do that without
    # fiddling with the geometry.
    V1a, V2a = set(V1), set(V2)
    V1b, V2b = set(V1), set(V2)

    V1a.add(tangent_points[0])
    V2a.add(tangent_points[1])
    V1b.add(tangent_points[1])
    V2b.add(tangent_points[0])

    return [{'V1': V1a, 'V2': V2a},
            {'V1': V1b, 'V2': V2b}]


def evaluate_split(V1, V2, power):
    p1 = abs(sum(power[v] for v in V1))
    p2 = abs(sum(power[v] for v in V2))
    size_ratio = max([len(V1)/len(V2), len(V2)/len(V1)])
    return {'power': max(p1, p2), 'size': size_ratio}

def draw_st_numbering(graph, solution, points, power=None):

    plt.gca().add_patch(plt.Circle(tuple(solution['circle']['center']), radius=solution['circle']['radius'], fill=False))

    for v, p in enumerate(points):
        color = 'r'
        if power and power[v] < 0:
            color = 'blue'
        plt.gca().add_patch(plt.Circle(tuple(p), radius=0.01, fc=color))

    for e in graph.es:
        plt.gca().add_line(render_line([points[e.source], points[e.target]], lw=0.5, color='red'))

    plt.gca().add_line(render_line(solution["tangent"], color='blue', ls='--'))

    for v, p in enumerate(solution["st_numbering"]):
        color = 'black'
        plt.gca().add_patch(plt.Circle(tuple(p), radius=0.01, fc=color))

    plt.axis('scaled')
    plt.show()

def ensure_dir(dir_name):
    if os.path.exists(dir_name):
        if not os.path.isdir(dir_name):
            raise Exception(dir_name + " exists but is not a directory.")
        return
    os.makedirs(dir_name)

def main(filename):
    ensure_dir("work")
    g = cog.load_graph("data/" + filename, False)

    power = assign_power(g)
    layout, best_cycle = cog.compute_positions(g, fixed_vertices=None, pull_to_vertices=False)

    vertex_file_name       = os.path.join("work", filename + '.' + 'positions')
    st_numbering_file_name = os.path.join("work", filename + '.' + 'st_numbering')

    with open(vertex_file_name, 'w') as vertex_file:
        for x, y in layout:
            vertex_file.write("{0}\n{1}\n".format(x, y))

    os.system('wykobi-master/make-st-numbering {0} {1}'.format(vertex_file_name, st_numbering_file_name))

    st_numbering = None
    with open(st_numbering_file_name) as st_numbering_file:
        # -nan stuff is hack to deal with division by zero or something errors
        # coming out of the C++ program.
        st_numbering = json.loads(''.join([l.replace('-nan', '"-nan"') for l in st_numbering_file.readlines()]))

    results = []
    for i, solution in enumerate(st_numbering["solutions"]):
        if( any(x[0] == '-nan' for x in solution["st_numbering"]) or
            any(x[1] == '-nan' for x in solution["st_numbering"])):
            print("Found a solution with -nan values in st_numbering")
            continue
        try:
            splits = split_network(solution)

            for j in range(2):
                splits[j]['quality'] = evaluate_split(splits[j]['V1'],
                                                      splits[j]['V2'],
                                                      power)
                splits[j]['solution'] = solution
                results.append(splits[j])
        except ZeroDivisionError:
            print("Got a division by zero error for solution {0}".format(i))
        
    best_power_value  = min([x['quality']['power'] for x in results])
    best_power_result = min([x for x in results if x['quality']['power'] == best_power_value], 
                            key=lambda y:y['quality']['size'])

    print("Best power result. Power = {power}. Size = {size}".format(**(best_power_result['quality']))) 
    draw_st_numbering(g, best_power_result['solution'], st_numbering["input_points"])

    best_size_value  = min([x['quality']['size'] for x in results])
    best_size_result = min([x for x in results if x['quality']['size'] == best_size_value], 
                            key=lambda y:y['quality']['power'])
    print("Best size result. Size = {size}. Power = {power}".format(**(best_size_result['quality']))) 
    draw_st_numbering(g, best_size_result['solution'], st_numbering["input_points"])

    return

    
if __name__=="__main__":
    example_to_use = 5
    filenames = ["triangle_graph.txt", 
                 "GMLFile.gml", 
                 "GMLFile2.gml", 
                 "GMLFile3.gml", 
                 "GMLFile4.gml", 
                 "power_subset.gml",]
    main(filenames[example_to_use])