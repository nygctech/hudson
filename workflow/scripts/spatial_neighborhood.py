import numpy as np
import skimage
import pandas as pd
from collections import namedtuple
import networkx as nx
import math
import sys
import getopt
from voronoi import *
import pickle
from pathlib import Path
import mudata as md

mdata = md.read(Path(snakemake.input[0]))
x_coord = mdata['morphological'].obsm['spatial'][:,0]
y_coord = mdata['morphological'].obsm['spatial'][:,1]

Point = namedtuple('Point', ['x','y'])
point_list = []
for x,y in zip(x_coord,y_coord):
    point_list.append(Point(x,y))
v,l,e = computeVoronoiDiagram(point_list)

#Building Graph
gr = nx.Graph()
for pt in v:
    gr.add_node(pt)
    
edge_list= []
for ei in range(len(e)):
    edge_list.append(e[ei][1:])
    
for ed in edge_list:
    gr.add_edge(ed[0],ed[1])
    
df = pd.read_csv(Path(snakemake.input[1]))
cell_type = {}
for i,ty in zip(df.id,df.type):
    cell_type.update({i: ty})
    
 
def get_neighbors(graph, node_val, rad):
    subgraph = nx.ego_graph(G = graph, n= node_val, radius=rad, center=True, undirected=False, distance=None)
    node_list = list(subgraph.nodes())
    return(node_list)
    

def get_cell_type(node_list):
    cell_frame = pd.DataFrame()
    cell_list = []
    final_node_list = []
    try:
        node_list.remove(0)
    except: 
        pass
    try:
        node_list.remove(-1)
    except:
        pass
    for n in node_list:
        try:
            cell_list.append(cell_type[n])
            final_node_list.append(n)
            
        except:
            node_list.remove(n)
            pass
    
    cell_frame['id'] = final_node_list
    cell_frame['type'] = cell_list
   
    return cell_frame 
        

def cell_pixel_intensity_gfap(frame):
    intensity_list = []
    for idx in frame.id:
        intensity_list.append(mean_intensities[idx - 1]) #minus 1 as we removed intensity for label 0
        
    frame['gfap_intensity'] = intensity_list
    
    
    return frame

def get_nodes(dt):
    selected_nodes = []
    for k in dt.keys():
        if dt[k] == 2:
            selected_nodes.append(k)
    try:
         selected_nodes.remove(-1) or selected_nodes.remove(0)
        
    except: 
          pass
    return(selected_nodes)

def n_th_neighbors(graph,node,limit):
    spatial_dict = nx.single_source_shortest_path_length(graph, node, cutoff=limit)
    return spatial_dict

neighbor_dict = {}
cell_type_per_node = {}

for g in range(len(v)):
    neighbor_dict.update({g:get_nodes(n_th_neighbors(gr,g,2))})
    
for g in range(len(v)):
    neigh = get_cell_type(get_nodes(n_th_neighbors(gr,node = g,limit = 2)))
    sub_df = neigh['type'].value_counts(normalize=True)
    cell_type_per_node.update({g:sub_df})
   
with open(Path(snakemake.output[0]), 'wb') as f:
     pickle.dump(cell_type_per_node, f)   

with open(Path(snakemake.output[1]), 'wb') as f:
    pickle.dump(neighbor_dict, f)  
        
        
    
