import pandas as pd
import numpy as np
import gget
import biomart
import json
from joblib import Parallel, delayed
import multiprocessing

n_cores = multiprocessing.cpu_count()
print(n_cores)

def get_ensembl_mappings(data_list):
    # Set up connection to server
    server = biomart.BiomartServer('http://www.ensembl.org/biomart')
    mart = server.datasets['mmusculus_gene_ensembl']
    # mart.show_filters()  # uses pprint
    # mart.show_attributes()
    data_list = ['ENSMUSG00000064372', 'ENSMUSG00000064371']
    response = mart.search({
        'filters': {
            'ensembl_gene_id': data_list[0],
        },
        #'attributes': ['ensembl_gene_id', 'accession']
        }, header = 1 )
    print(response.text)
    
    print(response)
    
    #response = mart.search({'attributes': ['ensembl_gene_id']})
                                                                                
    return 0

def parallel_gget(gene_id_list):
    gene_id_list = gene_id_list[0:2000]
    #splitting list in 10 parts
    gene_id_list = np.array_split(gene_id_list, n_cores)

    data = Parallel(n_jobs=n_cores)(delayed(gget.info)(gene_id_list[i]) for i in range(n_cores))

    df = pd.concat(data)
    results = df[['ensembl_id', 'uniprot_id']]
    results.to_csv('data/rna-seq/ensembl_uniprot_mapping.csv', index=False)


rna_norm_counts = pd.read_csv('data/rna-seq/GLDS-48_rna_seq_Normalized_Counts.csv')
gene_id_list = rna_norm_counts['Unnamed: 0'].tolist()

# get_ensembl_mappings(gene_id_list)

parallel_gget(gene_id_list)

    

# results = gget.info(gene_id_list)


# gens_names_to_convert = rna_norm_counts['Unnamed: 0'].to_dict()


