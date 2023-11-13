import pandas as pd
import numpy as np
import biomart
import json


def get_ensembl_mappings(data_dict):
    # Set up connection to server
    server = biomart.BiomartServer('http://www.ensembl.org/biomart')
    mart = server.datasets['mmusculus_gene_ensembl']    
    # mart.show_filters()  # uses pprint
    # mart.show_attributes()
    data_dict = {'ensembl_gene_id': 'ENSMUSG00000064372',}
    response = mart.search({
        'filters': {
            'ensembl_gene_id': data_dict['ensembl_gene_id'],
        },
        'attributes': ['ensembl_gene_id', 'accession']
        }, header = 1 )
    print(response.text)
    
    print(response)
    
    #response = mart.search({'attributes': ['ensembl_gene_id']})
                                                                                
    return 0

rna_norm_counts = pd.read_csv('data/rna-seq/GLDS-48_rna_seq_Normalized_Counts.csv')

gens_names_to_convert = rna_norm_counts['Unnamed: 0'].to_dict()
tmp = get_ensembl_mappings(gens_names_to_convert)


