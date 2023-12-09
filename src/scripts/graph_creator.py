import pandas as pd
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt
import time


class BioGraph():
    def __init__(self, gene_list = None, file_list = None, gene_list_to_discover = None, cluster = None) -> None:
        self._G = nx.Graph()
        self._aprox_G = nx.Graph()
        self._gene_list = gene_list
        self._file_list = file_list
        self._gene_to_discover = gene_list_to_discover
        self._cluster = cluster
    
    def set_gene_and_connections(self, df_list_index):
        file = self._file_list[df_list_index]
        gene_name = str(file.split('@')[1].split('.')[0])
        gene_df = pd.read_csv(f'data/FantomV/hs.FANTOM.annotated/{file}',  header= 1, index_col=0)
        gene_connections = gene_df[['Frel','gene_name']]
        gene_connections = gene_connections[gene_connections['Frel'] > 0.9]
        gene_connections = gene_connections.to_numpy()
        self.add_gene(gene_name)
        for weight, gene in gene_connections:
            mylist = self._gene_list
        #    self.add_gene(gene)
        #    if gene in self._gene_list:
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
        with open(f"figures/graph/{self._cluster}/graph_info.txt", "w") as f:
            f.write(f"Number of edges: {self._G.number_of_edges()}\n")
            f.write(f"Number of nodes: {self._G.number_of_nodes()}\n")
            f.write(f"Graph density: {nx.density(self._G)}\n")
            f.write(f"Graph average clustering: {nx.average_clustering(self._G)}\n")
            f.write(f"Graph average degree connectivity: {nx.average_degree_connectivity(self._G)}\n")
            # f.write(f"Graph average node connectivity: {nx.average_node_connectivity(self._G)}\n")
            # f.write(f"Graph average nearest neighbors degree: {nx.average_neighbor_degree(self._G)}\n")
            f.close()
        pos=nx.spring_layout(self._G, k=2)
        subax1 = plt.subplot(121)
        plt.figure(figsize=(10,10))
        nx.draw(self._G, pos, with_labels=True)
        plt.title(f'Dependencies Graph of Genes on {self._cluster}')

        plt.savefig(f"figures/graph/{self._cluster}/dependencies_graph_{self._cluster}.png")
        # subax2 = plt.subplot(122)
        # nx.draw_networkx_edges(self._G, pos=nx.circular_layout(self._G))
        plt.show()
        plt.figure(figsize=(10,10))
        nx.draw_networkx_nodes(self._G, pos, node_size=700)

        plt.title(f'Nodes representation of Genes {self._cluster}')
        plt.savefig(f"figures/graph/{self._cluster}/nodes_graph_{self._cluster}.png")
        plt.show()
        plt.figure(figsize=(10,10))
        nx.draw_networkx_edges(self._G, pos=pos, edgelist=self._G.edges(), width=6)
        # plt.show()
        # nx.draw_networkx_nodes(self._G, pos, node_size=700)
        # plt.title(f'Nodes representation of Genes {self._cluster}')
        # plt.savefig(f"C:/Users/ingma/PycharmProjects/Space-Biology_Data-Mining/figures/graph/{self._cluster}/nodes_graph_{self._cluster}.png")
        # plt.show()
        nx.draw_networkx_edges(self._G, pos=pos, edgelist=self._G.edges(), width=2)
        plt.title(f'Edges representation of dependencies between Genes on {self._cluster}')

        plt.savefig(f"figures/graph/{self._cluster}/edges_graph_{self._cluster}.png")
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

    def make_graph(self, gene_list, file_dict):
        self._G.add_nodes_from(gene_list)
        for gene in gene_list:
            try:
                for file in file_dict[gene]:
                    gene_df = pd.read_csv(
                        f'C:/Users/ingma/PycharmProjects/Space-Biology_Data-Mining/data/FantomV/hs.FANTOM.annotated/{file}',
                        header=1, index_col=0)
                    gene_connections = gene_df[['Frel', 'gene_name']]
                    gene_connections = gene_connections[gene_connections['Frel'] > 0.9]
                    gene_connections = gene_connections.to_numpy()
                    for weight, gene2 in gene_connections:
                        if gene2 in gene_list and not gene == gene2:
                            self.add_gene_link(gene, gene2, weight)
            except KeyError:
                # print(f'{gene} not in dataset')
                pass




def get_list(file):
    gene_list = open(file, 'r')
    return gene_list.read().split('\n')
        

if __name__ == '__main__':
    
    file = 'data/FantomV/cluster_4.txt'
    cluster = file.split('/')[-1].split('.')[0]
    print(cluster)
    gene_list = get_list(file)
    print(gene_list)
    
    complete_file_list = os.listdir('data/FantomV/hs.FANTOM.annotated/')
    
    file_dict = {}
    for f in complete_file_list:
        curr=f.split('@')[1].split('.')[0]
        if curr in gene_list:
            if curr not in file_dict.keys():
                file_dict[curr]=[f]
            else:
                file_dict[curr].append(f)

    test = BioGraph(cluster = cluster)
    tic = time.time()
    test.make_graph(gene_list, file_dict)
    toc = time.time()
    print(toc-tic)
    # test.approximize_graph()
    test.plot_graph()
    