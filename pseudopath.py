#!/usr/bin/python3
from llist import dllist
from collections import defaultdict
import matplotlib.pyplot as plt
import networkx as nx
import pdb
import unittest

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
    
def show_graph_with_pseudo(g, numbering, labels=None):
    if not labels:
        l = dict([(n, "{} ({})".format(n, numbering[n])) for n in g.nodes()])
    else:
        l = dict([(n, "{} ({})".format(labels[n], numbering[n])) for n in g.nodes()])

    show_graph(g, labels=l, node_size=700)

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
            print("Node {} has no valid predecessor".format(node))
            correct = False
        if not node == target and not found_a_successor:
            print("Node {} has no valid successor".format(node))
            correct = False

    return correct

def paper_ordering(ordering):
    l = list('stghfbcaed')
    ordered_list = [l[x[1]] for x in sorted([(o,n) for n, o in ordering.items()])]
    return ordered_list

def find_st_ordering(graph, source, target):
    """The st-ordering algorithm is from "Eager st-Ordering" by Ulrik
    Brandes.
    """

    # edge dicts are child to parent
    tree_edge              = {source: source, target: source}
    ear_dependency         = defaultdict(set)
    current_child          = {source:target}
    ordering               = dllist()
    ordering_nodes         = {}
    ordering_nodes[source] = ordering.appendleft(source)
    ordering_nodes[target] = ordering.append(target)

    # orientation - True is in direction of nodes in edge
    #             - False is in opposite direction of nodes in edge
    
    orientation            = {(source,target): True}

    def find_ear_path(x, v, w):
        """Find u on path from x to v that is in the ordering and
        is closest to v"""
        path = [v]
        while path[-1] != x and path[-1] not in ordering_nodes:
            path.append(tree_edge[path[-1]])
        path.reverse()
        path.append(w)
        return path

    def process_ears(tree_edge):
        w, x = tree_edge
        if tree_edge in ear_dependency:
            for v, w in ear_dependency[tree_edge]:
                ear_path = find_ear_path(x, v, w)
                for i, j  in zip(ear_path, ear_path[1:]):
                    orientation[(i,j)] = not orientation[(w,x)]
                u = ear_path[0]
                if orientation[(w,x)]:
                    insert_node = ordering_nodes[u]
                    for i in ear_path[1:-1]:
                        insert_node = ordering_nodes[i] = ordering.insert(i, insert_node)
                else:
                    insert_node = ordering_nodes[u].next
                    for i in ear_path[1:-1]:
                        insert_node = ordering_nodes[i] = ordering.insert(i, insert_node)
                        insert_node = insert_node.next
                for i, j  in zip(ear_path, ear_path[1:]):
                    process_ears((i, j))

            del(ear_dependency[tree_edge])

    def dfs(v):
        for w in sorted(graph.neighbors(v)):
            if tree_edge[v] == w:
                continue

            if w not in tree_edge:
                tree_edge[w] = v
                current_child[v] = w
                dfs(w)
            else:
                if w in current_child:
                    x = current_child[w]
                    ear_dependency[(w,x)].add((v,w))
                    if x in ordering_nodes:
                        process_ears((w,x))

    dfs(target)
    ordered_list = []
    node = ordering.first
    while node:
        ordered_list.append(node.value)
        node = node.next

    return dict([(x, i) for i, x in enumerate(ordered_list)])

class TestSTOrdering(unittest.TestCase):
    def test_bandes_example(self):
        g = nx.Graph()
        g.add_edge(0,1)
        g.add_edge(0,7)
        g.add_edge(0,5)
        g.add_edge(0,9)
        g.add_edge(1,2)
        g.add_edge(1,3)
        g.add_edge(2,3)
        g.add_edge(2,4)
        g.add_edge(2,8)
        g.add_edge(2,9)
        g.add_edge(4,5)
        g.add_edge(4,6)
        g.add_edge(5,7)
        g.add_edge(5,8)
        g.add_edge(6,9)

        ordering = find_st_ordering(g, 0, 1)
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 1))

    def test_make_pseudo_path_1(self):
        g = nx.Graph()
        g.add_edge(0,1)
        g.add_edge(1,2)
        g.add_edge(1,3)
        g.add_edge(2,4)
        g.add_edge(3,5)
        g.add_edge(5,4)
        g.add_edge(4,6)

        ordering = find_st_ordering(g, 0, 6)    
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 6))


    def test_make_pseudo_path_2(self):
        g = nx.Graph()
        g.add_edge(0,1)
        g.add_edge(1,2)
        g.add_edge(2,3)
        g.add_edge(3,4)
        ordering = find_st_ordering(g, 0, 4)    
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 4))

    def test_make_pseudo_path_3(self):
        g = nx.Graph()
        g.add_edge(0,1)
        g.add_edge(1,2)
        g.add_edge(1,3)
        g.add_edge(2,4)
        g.add_edge(3,4)
        g.add_edge(4,5)
        ordering = find_st_ordering(g, 0, 5)    
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 5))

    def test_make_pseudo_path_4(self):
        g = nx.Graph()
        g.add_edge(0,1)
        g.add_edge(1,2)
        g.add_edge(1,3)
        g.add_edge(2,4)
        g.add_edge(3,5)
        g.add_edge(5,6)
        g.add_edge(6,4)
        g.add_edge(4,7)
        ordering = find_st_ordering(g, 0, 7)    
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 7))

    def test_make_pseudo_path_5(self):
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
        ordering = find_st_ordering(g, 0, 7)    
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 7))

    def test_make_pseudo_path_6(self):
        g = nx.Graph()
        g.add_edge(0,1)
        g.add_edge(0,2)
        g.add_edge(1,3)
        g.add_edge(3,4)
        g.add_edge(4,5)
        g.add_edge(4,7)
        g.add_edge(2,6)
        g.add_edge(2,7)
        g.add_edge(6,5)
        ordering = find_st_ordering(g, 0, 7)    
        self.assertTrue(verify_pseudo_path(g, ordering, 0, 7))

        return

def test_bandes():
    g = nx.Graph()
    g.add_edge(0,1)
    g.add_edge(0,7)
    g.add_edge(0,5)
    g.add_edge(0,9)
    g.add_edge(1,2)
    g.add_edge(1,3)
    g.add_edge(2,3)
    g.add_edge(2,4)
    g.add_edge(2,8)
    g.add_edge(2,9)
    g.add_edge(4,5)
    g.add_edge(4,6)
    g.add_edge(5,7)
    g.add_edge(5,8)
    g.add_edge(6,9)

    pdb.set_trace()
    ordering = find_st_ordering(g, 0, 1)
    verify_pseudo_path(g, ordering, 0, 1)
    print(paper_ordering(ordering))
    show_graph_with_pseudo(g, ordering, labels=list('stghfbcaed'))
    return

if __name__ == '__main__':
    unittest.main()
        
