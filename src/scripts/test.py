import pandas as pd
import random
from graph_creator import BioGraph
import os
import networkx as nx
import statistics

edges = []
densities = []
clusterings = []

df = pd.read_csv('data/integrated/integrated_data_all.csv')

genes = df['hgnc_symbol'].tolist()
for i in range(30):
    rand_genes = random.choices(genes, k=100)

    complete_file_list = os.listdir('data/FantomV/hs.FANTOM.annotated/')

    file_dict = {}
    for f in complete_file_list:
        curr=f.split('@')[1].split('.')[0]
        if curr in rand_genes:
            if curr not in file_dict.keys():
                file_dict[curr]=[f]
            else:
                file_dict[curr].append(f)

    my_graph = BioGraph()
    my_graph.make_graph(rand_genes, file_dict)

    edges.append(my_graph._G.number_of_edges())
    densities.append(nx.density(my_graph._G))
    clusterings.append(nx.average_clustering(my_graph._G))

print('Mean of graph densities: ' + str(statistics.mean(densities)))
print('Mean of average clustering coefficients: ' + str(statistics.mean(clusterings)))
print('Standard deviation of graph densities: ' + str(statistics.pstdev(densities)))
print('Standard deviation of average clustering coefficients: ' + str(statistics.pstdev(clusterings)))
