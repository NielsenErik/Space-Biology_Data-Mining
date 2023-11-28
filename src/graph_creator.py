import pandas as pd
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt


class BioGraph():
    def __init__(self, gene_list = None, file_list = None, gene_list_to_discover = None) -> None:
        self._G = nx.Graph()
        self._aprox_G = nx.Graph()
        self._gene_list = gene_list
        self._file_list = file_list
        self._gene_to_discover = gene_list_to_discover
    
    def set_gene_and_connections(self, df_list_index):
        file = self._file_list[df_list_index]
        gene_name = str(file.split('@')[1].split('.')[0])
        gene_df = pd.read_csv(f'data/FantomV/test/{file}',  header= 1, index_col=0)
        gene_connections = gene_df[['Frel','gene_name']]
        gene_connections = gene_connections[gene_connections['Frel'] > 0.995]
        gene_connections = gene_connections.to_numpy()
        self.add_gene(gene_name)
        for weight, gene in gene_connections:
            self.add_gene(gene)
            self.add_gene_link(gene_name, gene, weight)
        self._G.remove_edges_from(nx.selfloop_edges(self._G))
        print(f"Gene set {gene_name} added to graph")
        
        
    def set_all_genes_connections(self):
        if self._gene_to_discover is None:
            print("No genes list to discover provided")
            print(self._file_list)
            for i, file in enumerate(self._file_list):
                self.set_gene_and_connections(i)
                print(f"Processing {file}")
        else:
            for index, file in enumerate(self._file_list):
                if str(file.split('@')[1].split('.')[0]) in self._gene_to_discover:
                    print(f"Processing {file}")
                    self.set_gene_and_connections(index)
    
        
    def check_gene_position(self, gene):
        for i in self._gene_position:
            if i[0] == gene:
                return i[1]
        return None
        
    def add_gene(self, gene):
        self._G.add_node(gene)
        
    def add_gene_link(self, gene1, gene2, weight):
        self._G.add_edge(gene1, gene2, weight=weight)
        
    def get_gene_graph(self):
        return self._G
    
    def plot_graph(self):
        
        print(self._G.number_of_edges())
        print(self._G.number_of_nodes())
        
        pos=nx.spring_layout(self._G)
        subax1 = plt.subplot(121)
        nx.draw(self._G, pos, with_labels=True)
        plt.title('Dependencies Graph of Genes')
        plt.savefig("figures/graph/dependencies_graph.png")
        # subax2 = plt.subplot(122)
        # nx.draw_networkx_edges(self._G, pos=nx.circular_layout(self._G))
        plt.show()
        nx.draw_networkx_nodes(self._G, pos, node_size=700)
        plt.title('Nodes representation of Genes')
        plt.savefig("figures/graph/nodes_graph.png")
        plt.show()
        nx.draw_networkx_edges(self._G, pos=pos, edgelist=self._G.edges(), width=6)
        plt.title('Edges representation of dependencies between Genes')
        plt.savefig("figures/graph/edges_graph.png")
        plt.show()
        
    def analyxe_graph(self):
        
        list(nx.connected_components(self._G))
        sorted(d for n, d in self._G.degree())
        cluster = nx.clustering(self._G)
    
    def approximize_graph(self):
        from networkx.algorithms import approximation
        self._aprox_G = approximation.k_components(self._G)
        print(self._aprox_G)
        self._aprox_G = approximation.average_clustering(self._G)
        print(self._aprox_G)
        self._aprox_G = approximation.node_connectivity(self._G)
        print(self._aprox_G)
        
        
    def __quick_test__(self, df):
        pass
        
        
        

if __name__ == '__main__':
    file_list = os.listdir('data/FantomV/test/')
    test = BioGraph(file_list=file_list)
    test.set_all_genes_connections()
    # test.approximize_graph()
    test.plot_graph()
    