#!/usr/bin/python3

import center_of_gravity_test_3 as cog
import json
import matplotlib.pyplot as plt
import os
import pdb

def render_line(line, **kw):
    xs = [line[0][0], line[1][0]]
    ys = [line[0][1], line[1][1]]
    return plt.Line2D(xs, ys, **kw)


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
    layout, best_cycle = cog.compute_positions(g, fixed_vertices=None, pull_to_vertices=True)

    vertex_file_name       = os.path.join("work", filename + '.' + 'positions')
    st_numbering_file_name = os.path.join("work", filename + '.' + 'st_numbering')

    with open(vertex_file_name, 'w') as vertex_file:
        for x, y in layout:
            vertex_file.write("{0}\n{1}\n".format(x, y))

    os.system('wykobi-master/make-st-numbering {0} {1}'.format(vertex_file_name, st_numbering_file_name))

    st_numbering = None
    with open(st_numbering_file_name) as st_numbering_file:
        st_numbering = json.loads(st_numbering_file.read())

    for solution in st_numbering["solutions"]:
        draw_st_numbering(g, solution, st_numbering["input_points"])

    return

    
if __name__=="__main__":
    example_to_use = 4
    filenames = ["triangle_graph.txt", "GMLFile.gml", "GMLFile2.gml", "GMLFile3.gml", "GMLFile4.gml"]
    main(filenames[example_to_use])
