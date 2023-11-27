import pandas as pd
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt


class BioGraph():
    def __init__(self, gene_list = None, file_list = None) -> None:
        self._G = nx.Graph()
        self._gene_list = gene_list
        self._file_list = file_list
        
    def set_gene_and_connections(self, df_list_index):
        file = self._file_list[df_list_index]
        gene_name = str(file.split('@')[1].split('.')[0])
        gene_df = pd.read_csv(f'data/FantomV/test/{file}',  header= 1, index_col=0)
        gene_connections = gene_df[['Frel','gene_name']]
        gene_connections = gene_connections[gene_connections['Frel'] > 0.5]
        gene_connections = gene_connections.to_numpy()
        
        self.add_gene(gene_name)
        for gene, weight in gene_connections:
            print(gene, weight)
            self.add_gene(gene)
            self.add_gene_link(gene_name, gene, weight)
        print(f"Finished to add {gene_name} and its connections")
        
    def add_gene(self, gene):
        self._G.add_node(gene)
        
    def add_gene_link(self, gene1, gene2, weight):
        self._G.add_edge(gene1, gene2, weight=weight)
        
    def get_gene_graph(self):
        return self._G
    
    def plot_graph(self):
        subax1 = plt.subplot(121)
        nx.draw(self._G, with_labels=True, font_weight='bold')
        subax2 = plt.subplot(122)
        nx.draw_shell(self._G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')
        plt.show()
        
    def analyxe_graph(self):
        
        list(nx.connected_components(self._G))
        sorted(d for n, d in self._G.degree())
        cluster = nx.clustering(self._G)
        print(cluster)
        
        
    def __quick_test__(self, df):
        pass
        
        
        

if __name__ == '__main__':
    file_list = os.listdir('data/FantomV/test')
    test = BioGraph(file_list=file_list)
    test.set_gene_and_connections(0)
    test.analyxe_graph()
    test.plot_graph()
    