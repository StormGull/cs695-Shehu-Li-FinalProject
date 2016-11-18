import igraph as ig
import re

comments_re = re.compile(r'#.*')


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
