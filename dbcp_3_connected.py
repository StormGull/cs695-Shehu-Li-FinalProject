#!/usr/bin/python3

import argparse
import find_x_embedding as fx
import json
import matplotlib.pyplot as plt
import os
import numpy as np
import pdb
import random
import time

class SubgraphNotConnectedException(Exception):
    pass

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
    try:
        return (line[1][1] - line[0][1])/(line[1][0] - line[0][0])
    except ZeroDivisionError:
        return float('inf')

def sort_st_numbering(tangent_line, st_numbering):
    if abs(slope(tangent_line)) > 1:
        key_func = lambda x: x[1][1]
    else:
        key_func = lambda x: x[1][0]
    ordering = [(v, p) for v, p in enumerate(st_numbering)]
    ordering.sort(key=key_func)
    return ordering

def split_network(solution):
    # order vertices by y position if slope > 1, otherwise
    # by x position
    st_numbering = solution['st_numbering']
    tangent_line = solution['tangent']
    tangent_points = [solution['pair_indices']['p1'], solution['pair_indices']['p2']]

    ordering = sort_st_numbering(tangent_line, st_numbering)
    
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

def draw_st_numbering(graph, solution, points, power=None, partitions=None, pair=None, file_name = None):

    st_numbering = solution["st_numbering"]

    plt.gca().add_patch(plt.Circle(tuple(solution['circle']['center']), radius=solution['circle']['radius'], fill=False))

    for v, p in enumerate(points):
        radius = 0.01
        if pair and p in pair:
            radius = .02
        color = 'r'
        if power and power[v] < 0:
            color = 'blue'
        plt.gca().add_patch(plt.Circle(tuple(p), radius=radius, fc=color))

    for e in graph.es:
        p1 = e.source
        p2 = e.target
        ls='-'

        if (partitions and
            ((p1 in partitions[0] and p2 in partitions[1]) or
             (p1 in partitions[1] and p2 in partitions[0]))):
            ls = ':'
        plt.gca().add_line(render_line([points[p1], points[p2]], lw=0.5, color='red', ls=ls))

    ordering = [x[1] for x in sort_st_numbering(solution['tangent'], st_numbering)]
    plt.gca().add_line(render_line([ordering[0],ordering[-1]], lw=0.5, color='blue', ls='--'))

    for v, p in enumerate(st_numbering):
        if v in partitions[0]:
            color = 'white'
        else:
            color = 'black'
        plt.gca().add_patch(plt.Circle(tuple(p), radius=0.02, fc=color))

    plt.axis('scaled')
    if file_name:
        plt.savefig(file_name)
    else:
        plt.show()
    plt.clf()

def ensure_dir(dir_name):
    if os.path.exists(dir_name):
        if not os.path.isdir(dir_name):
            raise Exception(dir_name + " exists but is not a directory.")
        return
    os.makedirs(dir_name)

def check_connectedness(g, V1, V2):
    for nodes in (V1, V2):
        subgraph = g.subgraph(nodes)
        if len(subgraph.components()) > 1:
            return False
    return True

#Scales a 2D point
def scale_point(point, scale):
    point_matrix = np.array([point[0], point[1], 1])
    scale_matrix = np.array([[scale, 0, 0],
                             [0, scale, 0],
                             [0,0,1]])
    result = np.matmul(scale_matrix, point_matrix)
    
    return (result[0], result[1])

def scale_positions(positions, scale):
    return [scale_point(point,scale) for point in positions]

def main(filename, write=False, save_images=False):
    ensure_dir("work")
    ensure_dir("reports")
    ensure_dir("images")

    g = fx.load_graph(filename, False)

    power = assign_power(g)
    layout, best_cycle, min_dist = fx.find_x_embedding_triconnected(g, fixed_vertices=None)
    
    # Scale graph so that the minimum distance between the points will not be too close
    min_acceptable_dist = 0.00001
    print("min_dist: {}".format(min_dist))
    if min_dist < min_acceptable_dist:
        scale_factor = min_acceptable_dist / min_dist
        layout = scale_positions(layout, scale_factor)

    vertex_file_name       = os.path.join(filename.replace('data/', 'work/') + '.' + 'positions')
    st_numbering_file_name = os.path.join(filename.replace('data/', 'work/') + '.' + 'st_numbering')

    with open(vertex_file_name, 'w') as vertex_file:
        for x, y in layout:
            vertex_file.write("{0}\n{1}\n".format(x, y))

    exit_code = os.system('make-st-numbering {0} {1}'.format(vertex_file_name, st_numbering_file_name))
    if exit_code:
        raise Exception("Running make-st-numbering failed:{}".format(exit_code))

    st_numbering = None
    with open(st_numbering_file_name) as st_numbering_file:

        # The -nan stuff is hack to deal with errors related to floating point precision
        # coming out of the C++ program. We had to hack the Boost C++ json-writing code
        # to write numbers without quotes, so it causes nan to not have quotes, which the
        # json parser doesn't like.
        #
        # However, we think we've solved the problems leading to the nan values in the first place,
        # but just in case...
        st_numbering = json.loads(''.join([l.replace('-nan', '"-nan"').replace('-1.#IND', '"-1.#IND"')
                                            for l in st_numbering_file.readlines()]))

    results = []
    for i, solution in enumerate(st_numbering["solutions"]):
        if( any(x[0] == '-nan' for x in solution["st_numbering"]) or
            any(x[1] == '-nan' for x in solution["st_numbering"]) or
            any(x[1] == '-1.#IND' for x in solution["st_numbering"]) or
            any(x[1] == '-1.#IND' for x in solution["st_numbering"])):
            print("Found a solution with -nan values in st_numbering")
            continue
        splits = split_network(solution)

        for j in range(2):
            # we are only interested in partitions where |V1|=|V2|
            # when |V| is even, or abs(|V1|-|V2|) = 1, when |V| is
            # odd.
            if len(g.vs)%2 == 0:
                if len(splits[0]['V1']) != len(splits[0]['V2']):
                    continue
            elif abs(len(splits[0]['V1']) - len(splits[0]['V2'])) != 1:
                continue

            # Check that this split produces two connected subgraphs.
            if not check_connectedness(g, splits[j]['V1'], splits[j]['V2']):
                print("ERROR: subgraphs not connected with split points {}".format(solution['pair_points']))
                # draw_st_numbering(g, solution,
                #                   st_numbering["input_points"], 
                #                   power=power, 
                #                   partitions=[splits[j]['V1'],
                #                               splits[j]['V2']],
                #                   pair=[solution['pair_points']['p1'],
                #                         solution['pair_points']['p2']],
                # )
            else:
                splits[j]['quality'] = evaluate_split(splits[j]['V1'],
                                                      splits[j]['V2'],
                                                      power)
                splits[j]['solution'] = solution
                results.append(splits[j])
        
    best_power_value  = min([x['quality']['power'] for x in results])
    best_power_result = min([x for x in results if x['quality']['power'] == best_power_value], 
                            key=lambda y:y['quality']['size'])

    best_size_value  = min([x['quality']['size'] for x in results])
    best_size_result = min([x for x in results if x['quality']['size'] == best_size_value], 
                            key=lambda y:y['quality']['power'])

    save_report(filename, g, power, st_numbering, best_power_result, best_size_result, write, save_images)

    if not save_images:
        print("Best power result. Power = {power}. Size = {size}".format(**(best_power_result['quality']))) 
        draw_st_numbering(g, best_power_result['solution'], 
                          st_numbering["input_points"], 
                          power=power, 
                          partitions=[best_power_result['V1'],
                                      best_power_result['V2']],
                          pair=[best_power_result['solution']['pair_points']['p1'],
                                best_power_result['solution']['pair_points']['p2'],]
        )

        if best_size_result != best_power_result:
            print("Best size result. Size = {size}. Power = {power}".format(**(best_size_result['quality']))) 
            draw_st_numbering(g, best_size_result['solution'],
                              st_numbering["input_points"], 
                              power=power, 
                              partitions=[best_size_result['V1'],
                                          best_size_result['V2']],
                              pair=[best_size_result['solution']['pair_points']['p1'],
                                    best_size_result['solution']['pair_points']['p2'],]
            )
        else:
            print("Best size result is the same as the best power result.")

def save_report(graph_file_name, g, power, st_numbering, best_power_result, best_size_result, write, save_images):
    image_file_names = []
    ensure_dir('images')
    ensure_dir('reports')
        
    file_name_base = time.strftime("algo_3-conn_{}_%Y-%m-%d_%H-%M-%S".format(graph_file_name.replace('data/','').replace('.gml', '').replace('.txt', '')))

    if save_images:
            image_file_name = "images/{}_best_power.png".format(file_name_base)
            draw_st_numbering(g, best_power_result['solution'], 
                              st_numbering["input_points"], 
                              power=power, 
                              partitions=[best_power_result['V1'],
                                          best_power_result['V2']],
                              pair=[best_power_result['solution']['pair_points']['p1'],
                                    best_power_result['solution']['pair_points']['p2'],],
                              file_name=image_file_name,
            )
            image_file_names.append(image_file_name)

            if best_power_result != best_size_result:
                image_file_name = "images/{}_best_size.png".format(file_name_base)
                draw_st_numbering(g, best_size_result['solution'], 
                                  st_numbering["input_points"], 
                                  power=power, 
                                  partitions=[best_size_result['V1'],
                                              best_size_result['V2']],
                                  pair=[best_size_result['solution']['pair_points']['p1'],
                                        best_size_result['solution']['pair_points']['p2'],],
                                  file_name=image_file_name,
                )
                image_file_names.append(image_file_name)
            else:
                image_file_names.append("Best size result is same as best power result.")
    if write:
        report_file_name = "reports/{}.json".format(file_name_base)
        print("Writing {0}".format(report_file_name))
        report_dict = {
            'graph_file': graph_file_name,
            'best_power_quality': best_power_result['quality'],
            'best_size_quality': best_size_result['quality'],
            'power': power,
        }

        if image_file_names:
            report_dict['image_file_names'] = image_file_names

        with open(report_file_name, 'w') as report_file:
            report_file.write(json.dumps(report_dict, indent=2, sort_keys=True))

    return


    
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                        help='Graph file to process.',
                        default=None)
    parser.add_argument('-w',
                        action='store_true',
                        help='Write results to disk.',
                        default=False)
    parser.add_argument('-i',
                        action='store_true',
                        help='Save images of solutions to disk.',
                        default=False,)
    args = parser.parse_args()

    main(args.f, write=args.w, save_images=args.i)
