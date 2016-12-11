#!/usr/bin/python3
import argparse
from   collections import defaultdict, namedtuple, deque
import display_graph as dg
import functools
import igraph as ig
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import os
import pdb
import pseudopath
import random
import sys

def show_graph(g, color_nodes=None, **kw):
    default_color = 'white'
    node_color = defaultdict(lambda: default_color)

    if color_nodes:
        for color, ids in color_nodes.items():
            for id in ids:
                node_color[id] = color
    
    nx.draw(g,with_labels=True,
            node_color=[node_color[id] for id in g.nodes()],
            **kw,
    )
    plt.show()
    
def show_graph_with_pseudo(g, numbering):
    labels = dict([(n, "{} ({})".format(n, numbering[n])) for n in g.nodes()])
    show_graph(g, labels=labels, node_size=700)

def assign_power(g):
    # Assign power supply 1, sink -1 randomly (but evenly)
    power = []
    for i, _ in enumerate(g):
        power.append(i%2 and 1 or -1)
    random.shuffle(power)
    return power

def ig_to_nx(g):
    nxg = nx.Graph()
    nxg.add_nodes_from(g.vs.indices)
    nxg.add_edges_from([(x.source, x.target) for x in g.es])
    return nxg

def nx_to_ig(g):
    node_map = dict([(n, i) for i, n in enumerate(sorted(g.nodes()))])
    edges = [(node_map[x], node_map[y]) for x,y in g.edges()]
    return ig.Graph(edges=edges)

def remove_degree_x_vertices(g, x, verbose=False):
    """Return a graph with all vertices of degree <= x. Since removing
    nodes will likely reduce the degree of other nodes, we have to do
    this iteratively until no nodes of degree <= k remain.

    """

    work_graph = nx.Graph(g)
    deg_dict = nx.degree(work_graph)
    iteration = 0
    while deg_dict and min(deg_dict.values()) <= x:
        ok_nodes = [n for n,d in deg_dict.items() if d > x]
        work_graph  = nx.Graph(work_graph.subgraph(ok_nodes))
        deg_dict = nx.degree(work_graph)
        iteration += 1
        if verbose:
            print("Iteration {0}. length deg_dict = {1}".format(iteration, len(deg_dict)))
    return work_graph

class SubgraphTooBigException(Exception):
    pass

class GraphNot2ConnectedException(Exception):
    pass

def find_two_connected_subgraphs(g, verbose=False):
    """Find the largest 2-connected subgraph of the graph passed in."""
    work_graph = remove_degree_x_vertices(g, 1, verbose)

    nodes_to_remove = set()
    iteration = 0
    for nodes in nx.all_node_cuts(work_graph, 1):
        if verbose:
            sys.stdout.write("{0}\r".format(iteration))
        iteration += 1
        nodes_to_remove = nodes_to_remove.union(nodes)

    if verbose:
        sys.stdout.write("\n")

    if nodes_to_remove:
        nodes_to_keep = set(work_graph.nodes()) - nodes_to_remove
        work_graph  = nx.Graph(work_graph.subgraph(list(nodes_to_keep)))
        components = list(nx.connected_components(work_graph))
        for component in components:
            yield(nx.Graph(work_graph.subgraph(list(component))))
    return

def find_case_1_two_connected_subgraphs(g, q, verbose=False):
    """Find subgraphs of g that has a single separating pair (ie, that is
    suitable for the first 2-connected case in section 3.2 of the DBCP
    paper. Each subgraph Gi must have |Vi| < |V|(q-1)/q. Make a
    generating function that returns valid subgraphs and nodes.

    """
    g = remove_degree_x_vertices(g, 1, verbose)

    all_nodes = g.nodes()
    size = len(all_nodes)
    q = float(q)
    max_subgraph_size = (q-1.0)/q*size
    
    if any(nx.all_node_cuts(g, 1)):
        raise GraphNot2ConnectedException()

    min_size = float('inf')
    found_one = False
    for nodes in nx.all_node_cuts(g, 2):
        sub_graph_nodes = [n for n in all_nodes if n not in nodes]
        disconnected_graph = nx.Graph(g.subgraph(sub_graph_nodes))
        try:
            for component in nx.connected_components(disconnected_graph):
                if verbose:
                    print("Subcomponent has size {0}".format(len(component)))
                if len(component) >= max_subgraph_size:
                    print("{} is too big!".format(len(component)))
                    min_size = min(min_size, len(component))
                    raise SubgraphTooBigException()
        except SubgraphTooBigException:
#            if verbose:
            print("For nodes {0} found a subgraph > (q-1)/q*size = {1}".format(nodes, max_subgraph_size))
            continue
        found_one=True
        yield(nodes, disconnected_graph)
    if not found_one:
        print("Found no suitable components. Smallest one was {}".format(min_size))
    return

def partition_pseudopaths(g, q, pseudopaths):
    S1 = []
    S2 = []
    S1_len = 0
    if len(pseudopaths[-1]) < len(g)/q:
        list_to_extend = S1
        for path in pseudopaths:
            list_to_extend.append(path)
            if list_to_extend == S1:
                S1_len += len(path)
            if list_to_extend == S1 and S1_len >= len(g)/(q-1):
                list_to_extend = S2
    else:
        S1 = pseudopaths[:-1]
        S2 = [pseudopaths[-1]]
    return S1, S2

def set_power(power, st):
    return sum([power[x] for x in st])

def validate_solution(g, power, V1, V2):
    print("V1 = {}. V2 = {}".format(V1, V2))

    p1 = abs(set_power(power, V1))
    p2 = abs(set_power(power, V2))
    size_ratio = max([len(V1)/len(V2), len(V2)/len(V1)])

    print("Power balance: {}. Size ratio: {}".format(max([p1, p2]), size_ratio))

    if not V1.intersection(V2):
        print("V1 and V2 do not intersect")
    else:
        print("ERROR: V1 and V2 intersect")

    if set(g.nodes()) == V1.union(V2):
        print("V1 and V2 have all nodes")
    else:
        print("ERROR: V1 and V2 do not contain all nodes")
    
    v1_subgraph = g.subgraph(V1)
    if len(list(nx.connected_components(v1_subgraph)))==1:
        print("Subgraph induced by V1 is connected.")
    else:
        print("ERROR: Subgraph induced by V1 is not connected.")
        
    v2_subgraph = g.subgraph(V2)
    if len(list(nx.connected_components(v2_subgraph)))==1:
        print("Subgraph induced by V2 is connected.")
    else:
        print("ERROR: Subgraph induced by V2 is not connected.")
        
    print("")

def do_case_1_two_connected_subgraphs(g, q, power, verbose=False):
    q *= 1.0 # make float
    V = set(g.nodes())
    size_V = len(V)

    print("Power: {}".format(power))

    find_set_power = functools.partial(set_power, power)

    found_count = 0
    solutions = []
    for nodes, subgraph in find_case_1_two_connected_subgraphs(g, q, verbose):
        found_count += 1
        nodes_list = list(nodes)
        pseudopaths = make_pseudopaths(g, nodes_list, subgraph, verbose)
        S1, S2 = partition_pseudopaths(g, q, pseudopaths)
        Q1 = [x for l in S1 for x in l]
        Q2 = [x for l in S2 for x in l]
        
        if len(Q2) > len(Q1):
            Q1, Q2 = Q2, Q1

        V1 = set([nodes_list[0]])
        V2 = set()

        if verbose:
            print("S1: {}".format(S1))
            print("S2: {}".format(S2))
            print("Q1: {}".format(Q1))
            print("Q2: {}".format(Q2))

        done = False
        for i, v in enumerate(Q1):
            # First just add |V|/q nodes from Q1 into V1
            if i + 1 < size_V/q:
                V1.add(v)
            else:
                # If we've add |V|/q nodes to V1 (so V2 is still empty)
                # check if we're done
                if not V2 and find_set_power(V1) == 0:
                    solutions.append([V1, V - V1])
                    done = True
                    break
                # Otherwise, start filling V2
                if not V2:
                    V2 = V1.copy()
                V2.add(v)

                # See if we're done
                if find_set_power(V2) == 0:
                    solutions.append([V2, V - V2])
                    done = True
                    break

        if done:
            continue
        # We got here, so we're not done yet
        # These are ordered from u -> v, but
        # now we want to work back from v.
        
        for i, v in enumerate([nodes_list[1]] + list(reversed(Q2))):
            V2.add(v)
            if find_set_power(V2) == 0:
                solutions.append([V2.difference([nodes_list[0]]), V.difference(V2).union([nodes_list[0]])])
                done = True
                break
            if len(V2) >= size_V * (q-1)/q:
                break

        if done:
            continue

        # Still no solution with p = 0.
        V3 = V.difference(V2)
        if abs(find_set_power(V1)) >= abs(find_set_power(V3)):
            V11 = set([nodes_list[0]])
            q1_idx = 0
            for v in Q1:
                if abs(find_set_power(V11)) == abs(find_set_power(V3)):
                    break
                V11.add(v)
            V_1 = V11.union(V3)
            solutions.append([V_1, V.difference(V_1)])
            done = True
            continue

        if abs(find_set_power(V1)) < abs(find_set_power(V3)):
            for v in Q2:
                V1.add(v)
                if find_set_power(V1) == 0:
                    break
            solutions.append(V1, V.difference(V1))
            done = True
            continue
                
    for solution in solutions:
        if verbose:
            show_graph(g, {'red': solution[0], 'yellow': solution[1]})
        validate_solution(g, power, solution[0], solution[1])

    return solutions

class InvalidPseudoPathException(Exception):
    pass

def make_pseudopaths(graph, separating_nodes, subgraph, verbose=False):
    numbering = pseudopath.find_st_ordering(graph, *separating_nodes)
    if not pseudopath.verify_pseudo_path(graph, numbering, *separating_nodes):
        print("{}".format(numbering))
        raise InvalidPseudoPathException()

    # lists of lists of nodes. 
    pseudopaths = []

    # Do a BFS through each subcomponent and add nodes to list of pseudopaths.
    seen_nodes = set()
    for neighbor in graph[separating_nodes[0]]:
        if neighbor not in seen_nodes:
            seen_nodes.add(neighbor)
            next_path = [neighbor]
            for s, t in nx.bfs_edges(subgraph, source=neighbor):
                if t not in seen_nodes:
                    seen_nodes.add(t)
                    next_path.append(t)
            next_path.sort(key=lambda x: numbering[x])
            pseudopaths.append(next_path)

    # Need to sort in increasing order by size:
    pseudopaths.sort(key=lambda x: len(x))
    return pseudopaths

def extract_power_k_connected(verbose=False):
    """Utility function to break Strogatz Watts power grid data into n-connected graphs."""
    nxg = ig_to_nx(ig.Graph.Read_GML('data/power.gml'))
    k_connected = nx.k_components(nxg)
    for k, components in k_connected.items():
        for i, component in enumerate(components):
            subgraph = nx.Graph(nxg.subgraph(component))
            file_name = "data/power_{}_connected_subset_{:>04}_{:>02}.gml".format(k, len(component), i)
            if verbose:
                print("Writing {}".format(file_name))
                nx_to_ig(subgraph).write_gml(file_name)

def main(file_name, q, power=None, verbose=False):
    g = ig_to_nx(ig.Graph.Read_GML(file_name))
    if not power:
        power = assign_power(g)
    do_case_1_two_connected_subgraphs(g, q, power, verbose=verbose)
    return

def test_1():
    main("data/power_2_connected_subset_0016_75.gml", 
         4, 
         power=[1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1],
         verbose=True,
    )

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                        help='Graph GML file to process',
                        default=None)
    parser.add_argument('-v',
                        action='store_true',
                        help='Verbose mode',
                        default=False)
    parser.add_argument('-q',
                        type=int,
                        help='Q constant to determine how large of subcomponents we should accept. Largest component can be (Q-1)/Q * |V|',
                        default=4)
    args = parser.parse_args()

    main(args.f, args.q, power=None, verbose=args.v)

    
