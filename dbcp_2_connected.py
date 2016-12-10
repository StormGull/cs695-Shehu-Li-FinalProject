#!/usr/bin/python3
from   collections import defaultdict, namedtuple, deque
import display_graph as dg
import functools
import igraph as ig
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import os
import pdb
import random
import sys

def show_graph(g, color_nodes=None):
    default_color = 'white'
    node_color = defaultdict(lambda: default_color)

    if color_nodes:
        for color, ids in color_nodes.items():
            for id in ids:
                node_color[id] = color
    
    nx.draw(g,with_labels=True,
            node_color=[node_color[id] for id in g.nodes()],
    )
    plt.show()
    
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

    for nodes in nx.all_node_cuts(g, 2):
        sub_graph_nodes = [n for n in all_nodes if n not in nodes]
        disconnected_graph = nx.Graph(g.subgraph(sub_graph_nodes))
        try:
            for component in nx.connected_components(disconnected_graph):
                if verbose:
                    print("Subcomponent has size {0}".format(len(component)))
                if len(component) >= max_subgraph_size:
                    if verbose:
                        print("Too big!")
                    raise SubgraphTooBigException()
        except SubgraphTooBigException:
            if verbose:
                print("For nodes {0} found a subgraph > (q-1)/q*size = {1}".format(nodes, max_subgraph_size))
            continue
        yield(nodes, disconnected_graph)
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

    print("Power balance: {}. Size ration: {}".format(max([p1, p2]), size_ratio))

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
        pseudopaths = make_pseudopaths(g, nodes_list, subgraph)
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
            show_graph(g)
            show_graph(subgraph)

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
        
        show_graph(g, {'red': V2})

        for i, v in enumerate([nodes_list[1]] + list(reversed(Q2))):
            V2.add(v)
            if find_set_power(V2) == 0:
                solutions.append([V2.difference(nodes_list[0]), V.minus(V2).union(nodes_list[0])])
                done = True
                break
            if len(V2) >= size_V * (q-1)/q:
                break

        if done:
            continue

        show_graph(g, {'red': V2})

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
                
    pdb.set_trace()
    for solution in solutions:
        show_graph(g, {'red': solution[0], 'yellow': solution[1]})
        validate_solution(g, power, solution[0], solution[1])

    return solutions

def make_pseudopaths(graph, separating_nodes, subgraph):
    # lists of lists of nodes. Each sub
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
            pseudopaths.append(next_path)

    # Need to sort in increasing order by size:
    pseudopaths.sort(key=lambda x: len(x))
    return pseudopaths

class NoNumberingFoundException(Exception):
    pass

class DupliciateNodeAssignmentException(Exception):
    pass

def pseudo_bfs(g, source, target):
    """Finding a pseudopath appears to be a sort of constraint
    satisfaction problem. Our solution is to use a BFS, but whenever
    we find a node n whose numbering would be N, but neighbors have
    all been numbered already (with values < N), we add a constraint
    that the neighbor of n with the max depth be given a numbering of
    at least N + 1.

    """

    todo_entry = namedtuple("todo_entry", ['node', 'depth'])

    forced_numbering = {}
    restart_required = True
    restart_nodes = defaultdict(set)
    result = None

    while restart_required :
        restart_required = False
        counter = 0
        numbering = {}
        todo = deque([todo_entry(source,0)])
        node_depth = defaultdict(lambda: -1)
        while todo:
            entry = todo.popleft()
            node_depth[entry.node] = entry.depth

            if entry.node in numbering:
                continue
            if entry.node in forced_numbering:
                numbering[entry.node] = forced_numbering[entry.node]
            else:
                neighbor_numberings = [numbering[p] for p in g.neighbors(entry.node) if p in numbering]
                if (not neighbor_numberings or
                    any([n < counter for n in neighbor_numberings])):
                    numbering[entry.node] = counter
                    counter += 1
                else:
                    numbering[entry.node] = max(neighbor_numberings) + 1

            if entry.node == target:
                continue
            found_suitable_neighbor = False
            neighbors = sorted(g.neighbors(entry.node))
            for neighbor in neighbors:
                if neighbor not in numbering:
                    found_suitable_neighbor = True
                    todo.append(todo_entry(neighbor, entry.depth + 1))
                elif (neighbor in forced_numbering and 
                      forced_numbering[neighbor] > numbering[entry.node]):
                    found_suitable_neighbor = True
            
            if not found_suitable_neighbor:
                selected_suitable_neighbor = False
                neighbors = sorted(g.neighbors(entry.node), 
                                   reverse=True, 
                                   key=lambda x: node_depth[x])
                for neighbor in neighbors:
                    # if we've already tried this one we won't try again
                    # but we can remove contraints on it.
                    if neighbor in restart_nodes[entry.node]:
                        if neighbor in forced_numbering:
                            del(forced_numbering[neighbor])
                        continue
                    else:
                        selected_suitable_neighbor = True
                        forced_numbering[neighbor] = counter
                        counter += 1
                        restart_nodes[entry.node].add(neighbor)
                        break
                if not selected_suitable_neighbor:
                    raise NoNumberingFoundException
                else:
                    restart_required = True
                    break
        result = numbering
    return result

def verify_pseudo_path(g, numbering, source, target):
    correct = True
    for node in g.nodes():
        found_a_successor = False
        found_a_predecessor = False
        for neighbor in g.neighbors(node):
            if numbering[neighbor] > numbering[node]:
                found_a_successor = True
            elif numbering[neighbor] < numbering[node]:
                found_a_predecessor = True
            else:
                raise DupliciateNodeAssignmentException()
        if not node == source and not found_a_predecessor:
            correct = False
        if not node == target and not found_a_successor:
            correct = False

    return correct

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
    if verbose:
        show_graph(g)
    if not power:
        power = assign_power(g)
    do_case_1_two_connected_subgraphs(g, q, power, verbose=verbose)
    return



def test_1():
    pdb.set_trace()
    main("data/power_2_connected_subset_0016_75.gml", 
         4, 
         power=[1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1],
         verbose=True,
    )
    

def test_make_pseudo_path_1():
    g = nx.Graph()
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(1,3)
    g.add_edge(2,4)
    g.add_edge(3,5)
    g.add_edge(5,4)
    g.add_edge(4,6)
#    show_graph(g)
#    pdb.set_trace()
    results = pseudo_bfs(g, 0, 6)    
    if (verify_pseudo_path(g, results, 0, 6) and 
        tuple(sorted(results.items())) == ((0, 0), (1, 1), (2, 2), (3, 3), (4, 6), (5, 4), (6, 7))):
        print("Test 1 succeeded")
    else:
        print("Test 1 failed")

    return

def test_make_pseudo_path_2():
    g = nx.Graph()
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(2,3)
    g.add_edge(3,4)
#    show_graph(g)
#    pdb.set_trace()
    results = pseudo_bfs(g, 0, 4)    
    if (verify_pseudo_path(g, results, 0, 4) and 
        tuple(sorted(results.items())) == ((0, 0), (1, 1), (2, 2), (3, 3), (4, 4))):
        print("Test 2 succeeded")
    else:
        print("Test 2 failed")
        
    return

def test_make_pseudo_path_3():
    g = nx.Graph()
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(1,3)
    g.add_edge(2,4)
    g.add_edge(3,4)
    g.add_edge(4,5)
#    show_graph(g)
#    pdb.set_trace()
    results = pseudo_bfs(g, 0, 5)    
    if (verify_pseudo_path(g, results, 0, 5) and 
        tuple(sorted(results.items())) == ((0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5))):
        print("Test 3 succeeded")
    else:
        print("Test 3 failed")

    return

def test_make_pseudo_path_4():
    g = nx.Graph()
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(1,3)
    g.add_edge(2,4)
    g.add_edge(3,5)
    g.add_edge(5,6)
    g.add_edge(6,4)
    g.add_edge(4,7)
#    show_graph(g)
#    pdb.set_trace()
    results = pseudo_bfs(g, 0, 7)    
    if (verify_pseudo_path(g, results, 0, 7) and 
        tuple(sorted(results.items())) == ((0, 0), (1, 1), (2, 2), (3, 3), (4, 7), (5, 4), (6, 5), (7, 8))):
        print("Test 4 succeeded")
    else:
        print("Test 4 failed")

    return

def test_make_pseudo_path_5():
    g = nx.Graph()
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(1,3)
    g.add_edge(2,4)
    g.add_edge(3,5)
    g.add_edge(3,6)
    g.add_edge(5,6)
    g.add_edge(6,4)
    g.add_edge(4,7)
#    show_graph(g)
#    pdb.set_trace()
    results = pseudo_bfs(g, 0, 7)    
    if (verify_pseudo_path(g, results, 0, 7) and 
        tuple(sorted(results.items())) == ((0, 0), (1, 1), (2, 2), (3, 3), (4, 7), (5, 4), (6, 5), (7, 8))):
        print("Test 5 succeeded")
    else:
        print("Test 5 failed")

    return


if __name__=='__main__':
    if len(sys.argv) > 1:
        file_name = sys.argv[1]
        verbose = len(sys.argv) > 2
        main(file_name, 4, power=None, verbose=verbose)
    else:
        test_1()

    
