import pandas as pd
import numpy as np
import gget
import biomart
import json
from joblib import Parallel, delayed
import multiprocessing
import os

n_cores = multiprocessing.cpu_count()
print(n_cores)

class Translator():
    def __init__(self, file):
        self._file = file
        self._df = pd.read_csv(file)
        self._data_list = self._df['ensembl_id'].tolist()
        
    def get_ensembl_mappings(self, id = 0):
        # Not working
        
        # Set up connection to server
        
        server = biomart.BiomartServer('http://www.ensembl.org/biomart')
        server.verbose = True
        mart = server.datasets['mmusculus_gene_ensembl']
        # mart.show_attributes()
        # server.show_databases()
        # mart.show_filters()  # uses pprint
        # mart.show_attributes()
        # data_list = data_list[0:20]
        attributes = ['ensembl_gene_id', 'uniprot_gn_id', ]# 'uniprot_isoform'
        response = mart.search({
            'filters': {
                'ensembl_gene_id': self._data_list,
            },
            'attributes': attributes
            }, header = 1 )
        print(response.text)
        save = open(f'data/rna-seq/ensembl_uniprot_mapping_{id}.txt', 'w')
        save.write(response.text)
        save.close()
        
        #response = mart.search({'attributes': ['ensembl_gene_id']})
                                                                                    
        return 0

    def parallel_gget(gene_id_list, id = 0):
        # Best working solution
        
        # gene_id_list = gene_id_list[0:2000]
        # splitting list in 10 parts
        print(f'Processing {id}')
        gene_id_list = np.array_split(gene_id_list, n_cores)
        data = Parallel(n_jobs=n_cores)(delayed(gget.info)(gene_id_list[i]) for i in range(n_cores))
        df = pd.concat(data)
        results = df[['ensembl_id', 'uniprot_id']]
        results.to_csv(f'data/rna-seq/ensembl_uniprot_mapping_{id}.csv', index=False)
        print(f'Finished {id}')
        
    def biomaoy_test(gene_id_list):
        # Not working
        import biomapy as bp
        #  Uniprot, entrez, ensembl, and Gene_symbol.

        gene_symbol_list=gene_id_list[0:100]
        gene_entrez_list=bp.gene_mapping_many(gene_symbol_list,'ensembl','entrez')
        print(gene_entrez_list)
    
    def start_translator(self):
        split_list = np.array_split(self._data_list, n_cores*100)
        #get_ensembl_mappings(gene_id_list)
        for i in range(634, n_cores*100):
            self.parallel_gget(split_list[i], id=i)
        print('Finished')
        
class TranslatorConcat():
    def __init__(self, file):
        self._file = file
        self.file_list = os.listdir(file)
        self.df = self.start_concatenator()
        
    def start_concatenator(self):
        df = pd.concat([pd.read_csv(f'{self._file}/{file}') for file in self.file_list])
        df = df.drop_duplicates(subset='ensembl_id', keep='first')
        return df
        


if __name__ == '__main__':
    # print(os.getcwd())
    # dir = 'data/translator_mapping'
    # concat_df = TranslatorConcat(dir)
    # concat_df = concat_df.start_concatenator()
    # save_df = concat_df.to_csv('data/translator_mapping.csv', index=False)
    
    prot_df = pd.read_csv('data/proteins/150929_KChatacharty_NASA_GeneLab_GroupA_CASIS_1_9_Fr1_TargetProtein.csv')
    translator = pd.read_csv('data/translator_mapping.csv')
    prot_list = prot_df['Accession'].tolist()
    translator_list = translator['uniprot_id'].tolist()
    
    count = 0
    for prot in prot_list:
        if prot in translator_list:
            count += 1
            
    print("Number of proteins in translator: ", count)