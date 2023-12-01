import pandas as pd
import numpy as np
import gget
import biomart
import json
from joblib import Parallel, delayed
import multiprocessing
import os
from pybiomart import Server

n_cores = multiprocessing.cpu_count()
print(n_cores)

class Translator():
    def __init__(self, file):
        self._file = file
        self._df = pd.read_csv(file)
        self._data_list = self._df['human_ensembl_gene_id'].tolist()
        
    def get_ensembl_mappings(self, gene_id_list):
        server = Server(host='http://www.ensembl.org')
        df = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'].query(attributes=['ensembl_gene_id', 'uniprotswissprot'], filters={'link_ensembl_gene_id': gene_id_list})
        print('Finished query')
        return df
        
    def parallel_biomart(self, gene_id_list):
        splitter = 150
        gene_id_list_ = np.array_split(gene_id_list, splitter)
        counter = len(gene_id_list)//splitter
        print(counter)
        for i in range(0, splitter):
            sublist = gene_id_list[counter*i:counter*(i+1)]
            df = self.get_ensembl_mappings(sublist)
            save = df.to_csv(f'data/translator_mapping/translator_mapping_human_ensembl_uniprot_{i}.csv', index=False)
            
        # df = pd.concat(data)
        # save = df.to_csv('data/translator_mapping_human_ensembl_uniprot.csv', index=False)

    def parallel_gget(self, gene_id_list):
        # Best working solution
        
        # gene_id_list = gene_id_list[0:2000]
        # splitting list in 10 parts
        gene_id_list = np.array_split(gene_id_list, n_cores*10)
        dfs = []
        data = Parallel(n_jobs=n_cores)(delayed(gget.info)(gene_id_list[i]) for i in range(n_cores*10))
        df = pd.concat(data)
        print(df.head())
        results = df[['ensembl_id', 'uniprot_id']]
        results.to_csv(f'data/rna-seq/ensembl_uniprot_mapping_{id}.csv', index=False)
        print(f'Finished {id}')
    
    def start_translator(self):
        # split_list = np.array_split(self._data_list, n_cores*100)
        # #get_ensembl_mappings(gene_id_list)
        # for i in range(0, n_cores*100):
        #     print(i)
        #     self.parallel_gget(split_list[i], id=i)
        self.parallel_biomart(self._data_list)
        
        print('Finished')
        
class Translator_Gene_Protein():
    def __init__(self, file, df = None):
        self._file = file
        self._file_list = os.listdir(file)
        self._df = df
        
    def start_concatenator(self):
        df = pd.concat([pd.read_csv(f'{self._file}/{file}') for file in self._file_list])
        df = df.dropna()
        return df
    
    def translations_counter(self, compare_list):
        if self._df is None:
            return 'No dataframe'
        count = 0
        self._df['Accession'] = self._df['Accession'].astype(str)
        for item in compare_list:
            if item in self._df['Accession'].tolist():
                count += 1
        return count
    
    def translations_merger(self, list_df):
        
        df = list_df[1]
        df['ensembl_id'] = df['ensembl_id'].astype(str).str.split('.').str[0]
        print(df.head())
        df = pd.merge(df, list_df[0], on='ensembl_id', how='left')
        df = pd.merge(df, list_df[2], on='Accession', how='inner')
        df = pd.merge(df, list_df[3], on='Accession', how='inner')
        
        return df
        
class Label_Translator():
    def __init__(self, list_df, mapper_path) -> None:
        self._list_df = list_df
        self._mapper = self.get_mapper(mapper_path)
    
    def get_mapper(self, mapper_path):
        mapper = pd.read_csv(mapper_path)
        mapper = mapper[['Sample Name', 'Label']]
        return mapper
        
    def from_tmt_to_label(self, df):
        print(self._mapper)
        for col in df.columns:
            col_ = col.split(' ')[-1]
            for row in self._mapper.iterrows():
                label = row[1]['Label']
                sample = row[1]['Sample Name']
                # print(col_, label, sample)
                if col_ == label:
                    print(col_, label, sample)
                    df = df.rename(columns={col: sample})
        return df
    
    def label_df_list(self):
        df_list = []
        for df in self._list_df:
            tmp = self.from_tmt_to_label(df)
            df_list.append(tmp)
        return df_list
    
    def save_dfs(self):
        df_list = self.label_df_list()
        counter = 0
        for df in df_list:
            df.to_csv(f'data/proteins/ProtonDiscoverer/renamed_labels_Proto_{counter}.csv', index=False)
            counter += 1
        return 0
    

class DataIntegrator():
    def __init__(self, rna_df, prot_df, mapper) -> None:
        self._rna_df = rna_df
        self._prot_df = prot_df
        self._mapper = mapper
        self._integrated_df = None
            
    def df_integration(self):
        
        for col in self._prot_df.columns.tolist():
            name = col + '_prot'
            if col != 'human_accession':
                self._prot_df = self._prot_df.rename(columns = {col : name})
        
        for col in self._rna_df.columns.tolist():
            name = col + '_rna'
            if col != 'mouse_ensembl_gene_id':
                self._rna_df = self._rna_df.rename(columns = {col : name})
        
        df = pd.merge(self._mapper, self._prot_df, on='human_accession', how='inner')
        df = pd.merge(df, self._rna_df, on='mouse_ensembl_gene_id', how='inner')
        
        print(df.head())
        print(df.shape)
        
        save_df = df.to_csv('data/integrated/integrated_data_1.csv', index=False)
        self._integrated_df = df
        
    def individuals_integration(self):
        for rna_col in self._rna_df.columns.tolist():
            for prot_cols in self._prot_df.columns.tolist():
                rna_col_ = rna_col.replace('_rna', '')
                prot_cols_ = prot_cols.replace('_prot', '')
                if rna_col_ == prot_cols_:
                    print(rna_col, prot_cols)
                    df = pd.DataFrame(self._integrated_df[['human_ensembl_gene_id','human_accession','mouse_ensembl_gene_id',rna_col, prot_cols]])
                    print(df.head())
                    print(df.shape)
                    save_df = df.to_csv(f'data/integrated/integrated_data_{rna_col_}_0.csv', index=False)
    
if __name__ == '__main__':
    # print(os.getcwd())
    mapper = pd.read_csv('data/mappers/mapper_human_mouse.csv')
    rna_df = pd.read_csv('data/rna_seq/GLDS-48_rna_seq_Normalized_Counts.csv')
    prot_df = pd.read_csv('data/proteins/ProtonDiscoverer/renamed_labels_Proto_0.csv')
    
    integrator = DataIntegrator(rna_df, prot_df, mapper)
    integrator.df_integration()
    integrator.individuals_integration()
    

    