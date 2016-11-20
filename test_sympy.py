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

def draw_st_numbering(circle, points, edges, tangent, tangent_points, lines_to_tangent, st_numbering, V1=None, V2=None, power=None):
    pc = lambda p: [p.x, p.y]

    plt.gca().add_patch(plt.Circle((circle.center.x, circle.center.y), radius=circle.radius, fill=False))

    for v, p in enumerate(points):
        color = 'r'
        if power and power[v] < 0:
            color = 'blue'
        plt.gca().add_patch(plt.Circle(pc(p), radius=0.04, fc=color))

    for e in edges:
        plt.gca().add_line(render_line(e, lw=2.0, color='red'))

    plt.gca().add_line(render_line(tangent_points, color='blue', ls='--'))

    if tangent:
        plt.gca().add_line(render_line(tangent, lw=2.0, color='blue'))

    for v, p in enumerate(st_numbering):
        color = 'black'
        if V2 and v in V2:
            color = 'white'
        plt.gca().add_patch(plt.Circle(pc(p), radius=0.03, fc=color))

    plt.axis('scaled')
    plt.show()

    return

def assign_power(g):
    # Assign power supply 1, sink -1 randomly (but evenly)
    power = []
    half  = len(g.vs)/2
    num_supply = 0
    num_sink = 0
    for v in g.vs.indices:
        if num_supply == half:
            power.append(-1)
            num_sink += 1
        elif num_sink == half:
            power.append(1)
            num_supply += 1
        elif random.random() < .5:
            power.append(-1)
            num_sink += 1
        else:
            power.append(1)
            num_supply += 1
    return power

def split_network(st_numbering, tangent_line, tangent_points):
    # order vertices by y position if slope > 1, otherwise
    # by x position
    if abs(tangent_line.slope) > 1:
        key_func = lambda x: x[1].y
    else:
        key_func = lambda x: x[1].x
    ordering = [(v, p) for v, p in enumerate(st_numbering)]
    ordering.sort(key=key_func)
    
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

    return V1a, V2a, V1b, V2b

def evaluate_split(V1, V2, power):
    p1 = sum(power[v] for v in V1)
    p2 = sum(power[v] for v in V2)
    size_ratio = max([len(V1)/len(V2), len(V2)/len(V1)])
    return [p1, p2, size_ratio]

def report_partition(V1, V2, p1, p2, size_ratio):
    print("""Split: size ratio={0}
    V1: {1} power={2}
    V2: {3} power={4}
""".format(size_ratio, 
           ",".join(map(str,sorted(V1))), p1,
           ",".join(map(str,sorted(V2))), p2))

def main(filename, fixed_vertices, draw):
    g = cs695util.load_tsv_edges("data/" + filename, False)
    positions = make_convex_x_embedding(g, fixed_vertices)

    fixed_points = [geo.Point(*positions[x]) for x in fixed_vertices]
    base_circle  = geo.Circle(*fixed_points)
    base_center  = base_circle.center
    base_radius  = base_circle.radius
    points       = [geo.Point(*v) for v in positions]
    edges        = [[e.source, e.target] for e in g.es]
    sym_edges    = [geo.Segment(points[e.source], points[e.target]) for e in g.es]
    power        = assign_power(g)

    partitions = []

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
                
                    # The line from the two vertices to the circle needs to be perpendicular
                    # to the resulting tangent line, so center the circle at the midpoint of
                    # the two vertices.
                    center = geo.Point((line.p1.x + line.p2.x)/2.0, (line.p1.y + line.p2.y)/2.0)
                    center_distance = base_center.distance(center)
                    circle = geo.Circle(center, base_radius + center_distance)

                    # Find intersection points for each line. There will be two for each, 
                    # one of which we'll filter later.
                    intersection_points = circle.intersection(line)

                    # Find tangent lines at the points of intersection
                    tangent_lines = circle.tangent_lines(intersection_points[0])
                    if not tangent_lines:
                        print("Found no tangent_lines")
                        draw_st_numbering(circle, points, sym_edges, None, line, None, [])
                        pdb.set_trace()
                        continue

                    tangent_line_1 = circle.tangent_lines(intersection_points[1])[0]

                    # Find intersection point for each 
                    # Just use one...not sure how to pick best, or if there is a best.
                    # Pick first for now.
                    st_numbering = []
                    # Find intersection for every point
                    lines_to_tangent = []
                    for p in points:
                        line_to_tangent = tangent_line_1.perpendicular_line(p)
                        lines_to_tangent.append(line_to_tangent)
                        projection_point = tangent_line_1.intersection(line_to_tangent)
                        st_numbering.append(projection_point[0])

                    V1a, V2a, V1b, V2b = split_network(st_numbering, tangent_line_1, (i, j))

                    partitions.append([V1a, V2a] + evaluate_split(V1a, V2a, power))
                    report_partition(*partitions[-1])
                    if draw:
                        draw_st_numbering(circle, points, sym_edges, tangent_line_1, line, lines_to_tangent, st_numbering, V1a, V2a, power)

                    partitions.append([V1b, V2b] + evaluate_split(V1b, V2b, power))
                    report_partition(*partitions[-1])
                    if draw:
                        draw_st_numbering(circle, points, sym_edges, tangent_line_1, line, lines_to_tangent, st_numbering, V1b, V2b, power)
    return

def get_dist(pos1, pos2):
    x1, y1 = pos1
    x2, y2 = pos2    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

if __name__=="__main__":
    draw = False
    example_to_use = 2
    filenames = ["square_graph.txt", "pentagon_graph.txt", 'triangle_graph.txt', 'triangle_graph-3.txt']
    fixed_vertices = [[0,1,2,3],[0,1,2,3,4], [0, 1, 2], [0, 1, 2]]
    
    main(filenames[example_to_use],fixed_vertices[example_to_use], draw)
