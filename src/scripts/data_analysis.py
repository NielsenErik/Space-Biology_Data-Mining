import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import math
from sklearn.ensemble import IsolationForest
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
from scipy.stats import norm

import warnings
warnings.filterwarnings('ignore')

class DataAnalyzer():
    def __init__(self, df) -> None:
        self._df = df
        self._gc_df = self.df_splitter()[0]
        self._flt_df = self.df_splitter()[1]
        self._rna_df = None
        self._prot_df = None
        self._rna_cols = []
    
    # human_ensembl_gene_id,human_accession,hgnc_symbol,mouse_ensembl_gene_id,external_gene_name,mouse_accession
    def df_splitter(self):
        gc_df = pd.DataFrame()
        flt_df = pd.DataFrame()
        gc_df['mouse_ensembl_gene_id'] = self._df['mouse_ensembl_gene_id']
        flt_df['mouse_ensembl_gene_id'] = self._df['mouse_ensembl_gene_id']
        gc_df['hgnc_symbol'] = self._df['hgnc_symbol']
        flt_df['hgnc_symbol'] = self._df['hgnc_symbol']
        gc_df['mouse_accession'] = self._df['mouse_accession']
        flt_df['mouse_accession'] = self._df['mouse_accession']
        
        
        for cols in self._df.columns.to_list():
            if '_GC_' in cols:
                gc_df[cols] = self._df[cols]
            elif '_FLT_' in cols:
                flt_df[cols] = self._df[cols]
            else:
                continue
        return gc_df, flt_df
    
    def compute_means(self, df):
        rna_cols = []
        prot_cols = []
        for cols in df.columns.to_list():
            if '_rna' in cols:
                rna_cols.append(cols)
            elif '_prot' in cols:
                prot_cols.append(cols)
        df['MEAN_rna'] = df[rna_cols].mean(axis=1)
        df['MEAN_prot'] = df[prot_cols].mean(axis=1)
        df['STD_rna'] = df[rna_cols].std(axis=1)
        df['STD_prot'] = df[prot_cols].std(axis=1)
        return df
    
    def df_exploration(self, df, sample):
       
        for col in df.columns.tolist():
            if 'MEAN_rna' in col:
                df = df.rename(columns = {col : 'RNA-Seq'})
            if 'MEAN_prot' in col:
                df = df.rename(columns = {col : 'Proteomics'})
        
        df['Protein-Rna_Ratios'] = df['Proteomics'] / df['RNA-Seq']
                    
        sns.set_theme(style="whitegrid")
        sns.heatmap(df.isnull(), cbar=False, cmap='viridis', xticklabels=False)
        plt.title(f'Missing values on sample\n{sample}')
        plt.savefig(f'figures/plots/rna-prot_means/missing_values_{sample}.png')
        plt.close()

        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna()
        df = df.drop_duplicates(keep = 'first')
        
        # ois = IsolationForest(random_state=0).fit(df[['RNA-Seq', 'Proteomics']])
        # df['outlier'] = ois.predict(df[['RNA-Seq', 'Proteomics']])
        #
        # clr = df['outlier'].map({1: 'blue', -1: 'red'})
        # plt.scatter(df['RNA-Seq'], df['Proteomics'], c=clr, s=2)
        # plt.xscale('log')
        # plt.title(f'Outliers on sample\n{sample}')
        # plt.xlabel('RNA-Seq')
        # plt.ylabel('Proteomics')
        # plt.savefig(f'figures/plots/rna-prot_means/outliers_{sample}.png')
        # plt.close()
        #
        # df = df[df['outlier'] == 1]
        
        plt.scatter(df['RNA-Seq'], df['Proteomics'], s=2)
        plt.xscale('log')
        plt.title(f'scatterplot sample\n{sample}')
        plt.xlabel('RNA-Seq')
        plt.ylabel('Proteomics')
        plt.savefig(f'figures/plots/rna-prot_means/outliers_removed_{sample}.png')
        plt.close()
        
        # df = df.drop(['outlier'], axis=1)
        return df
    
    def student_t_test_on_rna_prot_means(self):
        rna_gc = self._gc_df[['hgnc_symbol','mouse_ensembl_gene_id', 'RNA-Seq', 'STD_rna']]
        rna_flt = self._flt_df[['hgnc_symbol','mouse_ensembl_gene_id', 'RNA-Seq', 'STD_rna']]
        prot_gc = self._gc_df[['hgnc_symbol','mouse_ensembl_gene_id', 'Proteomics', 'STD_prot']]
        prot_flt = self._flt_df[['hgnc_symbol','mouse_ensembl_gene_id', 'Proteomics', 'STD_prot']]
        
        rna_gc = rna_gc.rename(columns = {'RNA-Seq' : 'RNA-Seq_GC'})
        rna_gc = rna_gc.rename(columns = {'STD_rna' : 'STD_rna_GC'})
        rna_flt = rna_flt.rename(columns = {'RNA-Seq' : 'RNA-Seq_FLT'})
        rna_flt = rna_flt.rename(columns = {'STD_rna' : 'STD_rna_FLT'})
        prot_gc = prot_gc.rename(columns = {'Proteomics' : 'Proteomics_GC'})
        prot_gc = prot_gc.rename(columns = {'STD_prot' : 'STD_prot_GC'})
        prot_flt = prot_flt.rename(columns = {'Proteomics' : 'Proteomics_FLT'})
        prot_flt = prot_flt.rename(columns = {'STD_prot' : 'STD_prot_FLT'})
        
        rna_merge = pd.merge(rna_gc, rna_flt, on=['hgnc_symbol','mouse_ensembl_gene_id'], how='inner')
        prot_merge = pd.merge(prot_gc, prot_flt, on=['hgnc_symbol','mouse_ensembl_gene_id'], how='inner')
        
        print(rna_merge.head())
        
        rna_merge['rna-t-test-statistic'] = (rna_merge['RNA-Seq_GC']-rna_merge['RNA-Seq_FLT'])/np.sqrt((rna_merge['STD_rna_GC']**2+rna_merge['STD_rna_FLT']**2)/2)
        rna_merge['rna-t-test-p-value'] = 2*(1-norm.cdf(abs(rna_merge['rna-t-test-statistic'])))
        prot_merge['prot-t-test-statistic'] = (prot_merge['Proteomics_GC']-prot_merge['Proteomics_FLT'])/np.sqrt((prot_merge['STD_prot_GC']**2+prot_merge['STD_prot_FLT']**2)/2)
        prot_merge['prot-t-test-p-value'] = 2*(1-norm.cdf(abs(prot_merge['prot-t-test-statistic'])))
        
        self.analyze_t_test(rna_merge, 'RNA-Seq')
        self.analyze_t_test(prot_merge, 'Proteomics')
        
        self.cluster_analysis(rna_merge, 'RNA-Seq')
        self.cluster_analysis(prot_merge, 'Proteomics')
        
        save = rna_merge.to_csv('data/processed/rna_means_t_test.csv', index=False)
        save = prot_merge.to_csv('data/processed/prot_means_t_test.csv', index=False)
        return rna_merge, prot_merge
    
    def analyze_t_test(self, df, title):
        # sns.histplot(x='rna-t-test-p-value', data=df)
        cols = ['rna-t-test-statistic', 'rna-t-test-p-value', 'prot-t-test-statistic', 'prot-t-test-p-value']
        if title == 'RNA-Seq':
            cols = ['rna-t-test-statistic', 'rna-t-test-p-value']
        elif title == 'Proteomics':
            cols = ['prot-t-test-statistic', 'prot-t-test-p-value']
        sns.histplot(x = cols[0], data=df)
        plt.title(f'{title} t-test statistic distribution')
        plt.savefig(f'figures/plots/t-test/{title}_t-test_statistic_distribution.png')
        # plt.show()
        sns.histplot(x = cols[1], data=df)
        plt.title(f'{title} t-test p-value distribution')
        plt.savefig(f'figures/plots/t-test/{title}_t-test_p-value_distribution.png')
        # plt.show()
    
    def cluster_eval(self, df, title):
        if title == 'RNA-Seq':
            col = 'rna-t-test-statistic'
        elif title == 'Proteomics':
            col = 'prot-t-test-statistic'
        scaler = StandardScaler()
        df_ = scaler.fit_transform(df[[col]])
        inertia = []
        silhouette = []
        for k in range(2, 10):
            kmeans = KMeans(n_clusters=k, random_state=0).fit(df_)
            inertia.append(kmeans.inertia_)
            silhouette.append(silhouette_score(df_, kmeans.labels_))
        figure, ax = plt.subplots(1, 2, figsize=(15, 5))
        ax[0].plot(range(2, 10), inertia)
        ax[0].set_xlabel('Number of clusters')
        ax[0].set_ylabel('Inertia')
        ax[0].set_title(f'Inertia for {title}')
        ax[1].plot(range(2, 10), silhouette)
        ax[1].set_xlabel('Number of clusters')
        ax[1].set_ylabel('Silhouette score')
        ax[1].set_title(f'Silhouette score for {title}')
        plt.savefig(f'figures/plots/clusters_analysis/{title}_cluster_eval_w_outliers.png')
        # plt.show()
        plt.close()
    
    def cluster_analysis(self, df, title):
        df_ = df.drop(['mouse_ensembl_gene_id'], axis=1)
        self.cluster_eval(df_, title)
        scaler = StandardScaler()
        
        if title == 'RNA-Seq':
            x = 'RNA-Seq_GC'
            y = 'RNA-Seq_FLT'
            col = 'rna-t-test-statistic'
        elif title == 'Proteomics':
            x = 'Proteomics_GC'
            y = 'Proteomics_FLT'
            col = 'prot-t-test-statistic'
        
        df_ = scaler.fit_transform(df_[[col]])
        
        kmeans = KMeans(n_clusters=5, random_state=0).fit(df_)
        df['cluster'] = kmeans.labels_

        sns.scatterplot(x=x, y=y, hue='cluster', data=df)
        if title=='RNA-Seq':
            plt.xscale('log')
            plt.yscale('log')
        plt.title(f'evaluating number of clusters')
        plt.savefig(f'figures/plots/clusters_analysis/{title}_clusters_w_outliers.png')
        # plt.show()
        plt.close()

        print(f'Silhouette score: {silhouette_score(df_, kmeans.labels_)}')
        
        return df
    
    def student_t_test_on_rna_prot_means_on_all_samples(self):
        
        gc_rna_cols = []
        gc_prot_cols = []
        flt_rna_cols = []
        flt_prot_cols = []
        
        gc_rna_cols.append('mouse_ensembl_gene_id')
        gc_prot_cols.append('mouse_ensembl_gene_id')
        flt_rna_cols.append('mouse_ensembl_gene_id')
        flt_prot_cols.append('mouse_ensembl_gene_id')

        for cols in self._gc_df.columns.tolist():
            if '_rna' in cols and cols != 'STD_rna':
                gc_rna_cols.append(cols)
            if '_prot' in cols and cols != 'STD_prot':
                gc_prot_cols.append(cols)
        for cols in self._flt_df.columns.tolist():
            if '_rna' in cols and cols != 'STD_rna':
                flt_rna_cols.append(cols)
            if '_prot' in cols and cols != 'STD_prot':
                flt_prot_cols.append(cols)
               
        gc_rna_df = self._gc_df[gc_rna_cols]
        gc_prot_df = self._gc_df[gc_prot_cols]
        flt_rna_df = self._flt_df[flt_rna_cols]
        flt_prot_df = self._flt_df[flt_prot_cols]
        
        rna_merge = pd.merge(gc_rna_df, flt_rna_df, on='mouse_ensembl_gene_id', how='inner')
        prot_merge = pd.merge(gc_prot_df, flt_prot_df, on='mouse_ensembl_gene_id', how='inner')
        
        sns.heatmap(rna_merge.isnull(), cbar=False, cmap='viridis', xticklabels=False)
        plt.title(f'Missing values on sample\nRNA-Seq both GC and FLT')
        # plt.show()
        plt.close()
        
        sns.heatmap(prot_merge.isnull(), cbar=False, cmap='viridis', xticklabels=False)
        plt.title(f'Missing values on sample\nProteomics both GC and FLT')
        # plt.show()
        plt.close()
        
        rna_merge = rna_merge.dropna()
        prot_merge = prot_merge.dropna()
        
        rna_merge_ = rna_merge.drop(['mouse_ensembl_gene_id'], axis=1)
        
        rna_res = []
        for row in rna_merge_.iterrows():
            res = ttest_ind(row[1][:len(gc_rna_cols)], row[1][len(gc_rna_cols):], equal_var=False)
            res = list(res)
            rna_res.append(res)
        rna_res = np.array(rna_res)
        rna_merge['t-test-statistic'] = rna_res[:, 0]
        rna_merge['t-test-p-value'] = rna_res[:, 1]
        
        prot_merge_ = prot_merge.drop(['mouse_ensembl_gene_id'], axis=1)
        prot_res = []
        for row in prot_merge_.iterrows():
            res = ttest_ind(row[1][:len(gc_prot_cols)], row[1][len(gc_prot_cols):], equal_var=False)
            prot_res.append(res)
        prot_res = np.array(prot_res)
        prot_merge['t-test-statistic'] = prot_res[:, 0]
        prot_merge['t-test-p-value'] = prot_res[:, 1]
        
        save = rna_merge.to_csv('data/processed/rna_t_test.csv', index=False)
        save = prot_merge.to_csv('data/processed/prot_t_test.csv', index=False)
        
        return rna_merge, prot_merge
    
    def genes_picker(self, df):
        for cluster in df['cluster_rna'].unique():
            
            cluster_df = df[df['cluster_rna'] == cluster]
            tmp = df.sort_values(by=['rna-t-test-p-value'], ascending=False)
            
            with open(f'data/FantomV/cluster_{cluster}.txt', 'w') as f:
                for row in range(100):
                    if len(tmp)>row:
                        f.write(tmp.iloc[row][0] + '\n')
                    else:
                        break
            
    
    def run(self):
        self._gc_df = self.compute_means(self._gc_df)
        self._flt_df = self.compute_means(self._flt_df)
        
        self._gc_df = self.df_exploration(self._gc_df, 'GC')
        self._flt_df = self.df_exploration(self._flt_df, 'FLT')
        
        self.student_t_test_on_rna_prot_means_on_all_samples()
        self._rna_df, self._prot_df = self.student_t_test_on_rna_prot_means()
        
        self._rna_df = self._rna_df.rename(columns = {'cluster' : 'cluster_rna'})
        self._prot_df = self._prot_df.rename(columns = {'cluster' : 'cluster_prot'})
        
        merged_df = pd.merge(self._rna_df, self._prot_df, on=['hgnc_symbol','mouse_ensembl_gene_id'], how='inner')
        print(merged_df)
        genes_pool = merged_df[((merged_df['cluster_rna'] == 1) & (merged_df['cluster_prot'] == 1)) |
                               ((merged_df['cluster_rna'] == 3) & (merged_df['cluster_prot'] == 2))]
        print(genes_pool)

        # with open(f'data/FantomV/outliers.txt',
        #           'w') as f:
        #     for row in genes_pool.index:
        #         f.write(genes_pool['hgnc_symbol'][row] + '\n')

        self.genes_picker(genes_pool)
        
        return self._gc_df, self._flt_df


if __name__ == '__main__':

    df = pd.read_csv('data/integrated/integrated_data_all.csv')
    analyzer = DataAnalyzer(df)
    gc_df, flt_df = analyzer.run()
    gc_df.to_csv('data/processed/gc_df.csv', index=False)
    flt_df.to_csv('data/processed/flt_df.csv', index=False)