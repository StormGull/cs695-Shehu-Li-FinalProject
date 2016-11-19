#!/usr/bin/python3
import igraph as ig
import math
import re
import random
import matplotlib.pyplot as plt
import pdb
from sympy import geometry as geo
import cs695util
from collections import namedtuple

def assign_initial_positions(g, fixed_vertices):
    return list(g.layout_kamada_kawai())

# I'm not sure why I redefined this function. Made
# sense at the time :-)
def assign_initial_positions(g, fixed_vertices):
    V = g.vs.indices
    positions = [None for v in V]

    # Just picked these since they're the same
    # as the triangle points in the
    # non-triangle version
    positions[fixed_vertices[0]] = [0.0, 0.0]
    positions[fixed_vertices[1]] = [1.0, 0.0]
    positions[fixed_vertices[2]] = [0.0, 1.0]

    seen = set()
    for v in V:
        if v not in fixed_vertices:
            while True:
                p = (random.random(), random.random())
                if p not in seen:
                    seen.add(p)
                    positions[v] = list(p)
                    break
    return positions

def make_convex_x_embedding(g, fixed_vertices):
    V = list(g.vs.indices)    
    positions = assign_initial_positions(g, fixed_vertices)
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
    return positions

def render_line(line, **kw):
    xs = [line.p1.x, line.p2.x]
    ys = [line.p1.y, line.p2.y]
    return plt.Line2D(xs, ys, **kw)

def draw_st_numbering(circle, points, edges, tangent, tangent_points, lines_to_tangent, st_numbering, positions):
    pc = lambda p: [p.x, p.y]

    plt.gca().add_patch(plt.Circle((circle.center.x, circle.center.y), radius=circle.radius, fill=False))

    for p in points:
        color = 'r'
        if p in tangent_points:
            color = 'blue'
        plt.gca().add_patch(plt.Circle(pc(p), radius=0.03, fc=color))

    for e in edges:
        plt.gca().add_line(render_line(e, lw=2.0, color='red'))

    if tangent:
        plt.gca().add_line(render_line(tangent, lw=2.0, color='blue'))

    for p in st_numbering:
        plt.gca().add_patch(plt.Circle(pc(p), radius=0.03, fc='green'))

    plt.axis('scaled')
    plt.show()

    return

def main(filename, fixed_vertices):
    g = cs695util.load_tsv_edges("data/" + filename, False)
    positions = make_convex_x_embedding(g, fixed_vertices)

    fixed_points = [geo.Point(*positions[x]) for x in fixed_vertices]
    base_circle  = geo.Circle(*fixed_points)
    base_center  = base_circle.center
    base_radius  = base_circle.radius
    points       = [geo.Point(*v) for v in positions]
    edges        = [[e.source, e.target] for e in g.es]
    sym_edges    = [geo.Segment(points[e.source], points[e.target]) for e in g.es]

    # Find the line formed by each pair of vertices
    lines = []
    seen = set()
    for i, v1 in enumerate(points):
        for j, v2 in enumerate(points):
            if j != i and (i,j) not in seen:
                line = geo.Line(v1, v2)
                seen.add((i, j))
                seen.add((j, i))
                if line.slope >= 0.0:
                    #                lines.append(line_with_vertices(line, i, j))
                    lines.append(line)
                
                    center = geo.Point((line.p1.x + line.p2.x)/2.0, (line.p1.y + line.p2.y)/2.0)
                    center_distance = base_center.distance(center)
                    circle = geo.Circle(center, base_radius + center_distance)

                    # Find intersection points for each line. There will be two for each, 
                    # one of which we'll filter later.
                    intersection_points = circle.intersection(line)

                    # Find tangent lines at the points of intersection
                    line_tangents = circle.tangent_lines(intersection_points[0])
                    if not line_tangents:
                        print("Found no line_tangents")
                        draw_st_numbering(circle, points, sym_edges, None, line, None, [], [])
                        pdb.set_trace()
                        continue

                    line_tangent_1 = circle.tangent_lines(intersection_points[1])[0]

                    # Find intersection point for each 
                    # Just use one...not sure how to pick best, or if there is a best.
                    # Pick first for now.
                    st_numbering = []
                    # Find intersection for every point
                    lines_to_tangent = []
                    for p in points:
                        line_to_tangent = line_tangent_1.perpendicular_line(p)
                        lines_to_tangent.append(line_to_tangent)
                        projection_point = line_tangent_1.intersection(line_to_tangent)
                        st_numbering.append(projection_point[0])
    
                    draw_st_numbering(circle, points, sym_edges, line_tangent_1, line, lines_to_tangent, st_numbering, positions)
    return

def get_dist(pos1, pos2):
    x1, y1 = pos1
    x2, y2 = pos2    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

if __name__=="__main__":
    example_to_use = 2
    filenames = ["square_graph.txt", "pentagon_graph.txt", 'triangle_graph.txt', 'triangle_graph-3.txt']
    fixed_vertices = [[0,1,2,3],[0,1,2,3,4], [0, 1, 2], [0, 1, 2]]
    
    main(filenames[example_to_use],fixed_vertices[example_to_use])
