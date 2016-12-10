#!/usr/bin/python3
import igraph as ig
import sys

def show_graph(g, write_to_file=False, **kw):
    visual_style={}
    visual_style['bbox']=(1200,1200)
    visual_style['margin']=50
    visual_style['vertex_label']=g.vs.indices
    visual_style['label_dist']=1

    layout = g.layout(len(g.vs) > 100 and 'lgl' or 'kk')

    visual_style['vertex_color']='blue'
    visual_style.update(kw)
    if write_to_file:
        image_file_name = file_name.replace('.gml', '.png').replace('data/','images/')
        ig.plot(g, image_file_name, layout=layout, **visual_style)
    else:
        ig.plot(g, layout=layout, **visual_style)

def show_file(file_name,  write_to_file=False, **kw):
    g = ig.Graph.Read_GML(file_name)
    show_graph(g, write_to_file, **kw)

if __name__=='__main__':
    show_file(sys.argv[1], len(sys.argv) > 2)
