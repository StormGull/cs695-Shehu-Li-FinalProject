#!/usr/bin/python3

import igraph as ig
import networkx as nx
import pdb
import sys

def main():
    powerg = ig.Graph.Read_GML('data/power.gml')
    powerx = nx.Graph()
    powerx.add_nodes_from(powerg.vs.indices)
    powerx.add_edges_from([(x.source, x.target) for x in powerg.es])

    work_graph = nx.Graph(powerx)
    deg_dict = nx.degree(work_graph)
    iteration = 0
    k = 2
    while min(deg_dict.values()) < k:
        deg2_nodes = [n for n,d in deg_dict.items() if d > k - 1]
        work_graph  = nx.Graph(work_graph.subgraph(deg2_nodes))
        deg_dict = nx.degree(work_graph)
        iteration += 1
        print("Interation {0}. length deg_dict = {1}".format(iteration, len(deg_dict)))
        
    nodes_to_remove = set()
    iteration = 0
    for nodes in nx.all_node_cuts(work_graph, 1):
        sys.stdout.write("{0}\r".format(iteration))
        iteration += 1
        nodes_to_remove = nodes_to_remove.union(nodes)

    sys.stdout.write("\n")
    nodes_to_keep = set(work_graph.nodes()) - nodes_to_remove
    work_graph  = nx.Graph(work_graph.subgraph(list(nodes_to_keep)))
    components = list(nx.connected_components(work_graph))
    largest_component = max(components, key=lambda x: len(x))
    work_graph = nx.Graph(work_graph.subgraph(list(largest_component)))
    
    k = 3
    iteration = 0
    while min(deg_dict.values()) < k:
        deg3_nodes = [n for n,d in deg_dict.items() if d > k - 1]
        work_graph  = nx.Graph(work_graph.subgraph(deg3_nodes))
        deg_dict = nx.degree(work_graph)
        iteration += 1
        print("Interation {0}. length deg_dict = {1}".format(iteration, len(deg_dict)))

    nodes_to_keep = set(work_graph.nodes()) - nodes_to_remove
    work_graph  = nx.Graph(work_graph.subgraph(list(nodes_to_keep)))
    components = list(nx.connected_components(work_graph))
    largest_component = max(components, key=lambda x: len(x))
    work_graph = nx.Graph(work_graph.subgraph(list(largest_component)))


    iteration = 0
    nodes_to_remove = set()
    for nodes in nx.all_node_cuts(work_graph, 2):
        sys.stdout.write("{0}\r".format(iteration))
        iteration += 1
        nodes_to_remove = nodes_to_remove.union(nodes)

    if not nodes_to_remove:
        node_map = dict([(n, i) for i, n in enumerate(sorted(work_graph.nodes()))])
        edges = [(node_map[x], node_map[y]) for x,y in work_graph.edges()]
        g = ig.Graph(edges=edges)
        g.write_gml("data/power_subset.gml")

    return

if __name__=='__main__':
    main()
    
