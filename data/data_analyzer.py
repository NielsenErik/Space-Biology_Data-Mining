import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.ensemble import IsolationForest
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


class DataAnalyzer():
    def __init__(self, df_gc_list, df_flt_list, gc_sample, flt_sample, gc_all = None, flt_all = None) -> None:
        self._df_gc_list = df_gc_list
        self._df_flt_list = df_flt_list
        self._gc_sample = gc_sample
        self._flt_sample = flt_sample
        self._gc_mean_df = gc_all
        self._flt_mean_df = flt_all
        
    def df_exploration(self, df, sample, all = False):
        if all == False:
            for col in df.columns.tolist():
                if '_rna' in col:
                    df = df.rename(columns = {col : 'RNA-Seq'})
                if '_prot' in col:
                    df = df.rename(columns = {col : 'Proteomics'})
        else:
            for col in df.columns.tolist():
                if 'MEAN_rna' in col:
                    df = df.rename(columns = {col : 'RNA-Seq'})
                if 'MEAN_prot' in col:
                    df = df.rename(columns = {col : 'Proteomics'})
                    
        sns.set_theme(style="whitegrid")
        sns.heatmap(df.isnull(), cbar=False, cmap='viridis', xticklabels=False)
        plt.title(f'Missing values on sample\n{sample}')
        plt.savefig(f'figures/plots/rna-prot_per_samples/missing_values_{sample}.png')
        plt.close()
        
        df = df.dropna()
        df = df.drop_duplicates(keep = 'first')
        
        ois = IsolationForest(random_state=0).fit(df[['RNA-Seq', 'Proteomics']])
        df['outlier'] = ois.predict(df[['RNA-Seq', 'Proteomics']])
        
        clr = df['outlier'].map({1: 'blue', -1: 'red'})
        plt.scatter(df['RNA-Seq'], df['Proteomics'], c=clr)
        plt.title(f'Outliers on sample\n{sample}')
        plt.xlabel('RNA-Seq')
        plt.ylabel('Proteomics')
        plt.savefig(f'figures/plots/rna-prot_per_samples/outliers_{sample}.png')
        plt.close()
        
        df = df[df['outlier'] == 1]
        
        plt.scatter(df['RNA-Seq'], df['Proteomics'])
        plt.title(f'Outliers removed on sample\n{sample}')
        plt.xlabel('RNA-Seq')
        plt.ylabel('Proteomics')
        plt.savefig(f'figures/plots/rna-prot_per_samples/outliers_removed_{sample}.png')
        plt.close()
        return df
    
    def df_list_exploration(self):
        df_list = []
        for i in range(len(self._df_gc_list)):
            df = self._df_gc_list[i]
            sample = self._gc_sample[i]
            tmp = self.df_exploration(df, sample, all=False)
            tmp = self.get_cluster(tmp, sample)
            self._df_gc_list[i] = tmp
        for i in range(len(self._df_flt_list)):
            df = self._df_flt_list[i]
            sample = self._flt_sample[i]
            tmp = self.df_exploration(df, sample, all=False)
            tmp = self.get_cluster(tmp, sample)
            self._df_flt_list[i] = tmp
        return df_list
    
    def eval_n_clusters(self, df, sample):
        inertia = []
        silhouette = []
        for i in range(2, 11):
            kmeans = KMeans(n_clusters=i, random_state=0)
            clusters = kmeans.fit_predict(df[['RNA-Seq', 'Proteomics']])
            silhouette.append(silhouette_score(df[['RNA-Seq', 'Proteomics']], clusters))
            inertia.append(kmeans.inertia_)
        figure, ax = plt.subplots(1, 2, figsize=(16,8))
        ax[0].plot(range(2, 11), inertia)
        ax[0].set_title('Elbow Method')
        ax[0].set_xlabel('Number of clusters')
        ax[0].set_ylabel('Inertia')
        ax[1].plot(range(2, 11), silhouette)
        ax[1].set_title('Silhouette Method')
        ax[1].set_xlabel('Number of clusters')
        ax[1].set_ylabel('Silhouette')
        plt.savefig(f'figures/plots/rna-prot_per_samples/elbow_silhouette{sample}.png')
        plt.close()
        
        
    def get_cluster(self, df, sample):
        print("CLUSTERING")
        print(df.columns.tolist())
        
        self.eval_n_clusters(df, sample=sample)
        
        cluster = KMeans(n_clusters=5, random_state=0)
        
        X = df[['RNA-Seq', 'Proteomics']]
        
        df['cluster'] = cluster.fit_predict(X)
        plt.figure(figsize=(10, 10))
        
        plt.scatter(X['RNA-Seq'], X['Proteomics'], c=df['cluster'], cmap='viridis')
        plt.scatter(cluster.cluster_centers_[:, 0], cluster.cluster_centers_[:, 1], c='red', marker='x')
        
        plt.title(f'Cluster on sample\n{sample}')
        plt.xlabel('RNA-Seq')
        plt.ylabel('Cluster')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'figures/plots/rna-prot_per_samples/cluster_{sample}.png')
        plt.close()
        
        return df
    
    def get_means_df(self):
        gc_df = self._df_gc_list[0]
        flt_df = self._df_flt_list[0]
        for i in range(1, len(self._df_gc_list)):
            gc_df = pd.merge(gc_df, self._df_gc_list[i], on=['human_ensembl_gene_id','human_accession','mouse_ensembl_gene_id'], how='inner')
        for i in range(1, len(self._df_flt_list)):
            flt_df = pd.merge(flt_df, self._df_flt_list[i], on=['human_ensembl_gene_id','human_accession','mouse_ensembl_gene_id'], how='inner')
        gc_rna_cols = []
        for col in gc_df.columns.tolist():
            if '_rna' in col:
                gc_rna_cols.append(col)
        gc_prot_cols = []
        for col in gc_df.columns.tolist():
            if '_prot' in col:
                gc_prot_cols.append(col)
        flt_rna_cols = []
        for col in flt_df.columns.tolist():
            if '_rna' in col:
                flt_rna_cols.append(col)
        flt_prot_cols = []
        for col in flt_df.columns.tolist():
            if '_prot' in col:
                flt_prot_cols.append(col)
        
        gc_df = gc_df.drop_duplicates(keep = 'first')
        flt_df = flt_df.drop_duplicates(keep = 'first')
        print("CALCULATING MEANS")
        
        gc_df['MEAN_rna'] = gc_df[gc_rna_cols].mean(axis=1)
        gc_df['STD_rna'] = gc_df[gc_rna_cols].std(axis=1)
        gc_df['MEAN_prot'] = gc_df[gc_prot_cols].mean(axis=1)
        gc_df['STD_prot'] = gc_df[gc_prot_cols].std(axis=1)
        
        flt_df['MEAN_rna'] = flt_df[flt_rna_cols].mean(axis=1)
        flt_df['STD_rna'] = flt_df[flt_rna_cols].std(axis=1)
        flt_df['MEAN_prot'] = flt_df[flt_prot_cols].mean(axis=1)
        flt_df['STD_prot'] = flt_df[flt_prot_cols].std(axis=1)
        
        for cols in gc_df.columns.tolist():
            if 'Mmus_' in cols:
                gc_df = gc_df.drop([cols], axis=1)
        for cols in flt_df.columns.tolist():
            if 'Mmus_' in cols:
                flt_df = flt_df.drop([cols], axis=1)
        
        save_gc = gc_df.to_csv('data/integrated/all_gc.csv', index=False)
        save_flt = flt_df.to_csv('data/integrated/all_flt.csv', index=False)
        
        self._gc_mean_df = gc_df
        self._flt_mean_df = flt_df
            
        return gc_df, flt_df
    
    def analyze_means_exploration(self, gc_df, flt_df):
        self._gc_mean_df = gc_df
        self._flt_mean_df = flt_df
        print("EXPLORATION")
        self._gc_mean_df = self.df_exploration(self._gc_mean_df, 'GC', all=True)
        self._flt_mean_df = self.df_exploration(self._flt_mean_df, 'FLT', all=True)
        _ = self.get_cluster(self._gc_mean_df, sample='GC')
        _ = self.get_cluster(self._flt_mean_df, sample='FLT')
        
if __name__ == '__main__':
    file_list = os.listdir('data/integrated/')
    gc_file = []
    flt_file = []
    gc_sample = []
    flt_sample = []
    for file in file_list:
        if '_GC_' in file:
            df = pd.read_csv(f'data/integrated/{file}')
            gc_file.append(df)
            file = file.replace('integrated_data_', '')
            file = file.replace('.csv', '')
            gc_sample.append(file)
        if '_FLT_' in file:
            df = pd.read_csv(f'data/integrated/{file}')
            flt_file.append(df)
            file = file.replace('integrated_data_', '')
            file = file.replace('.csv', '')
            flt_sample.append(file)
    gc_df = pd.read_csv('data/integrated/all_gc.csv')
    flt_df = pd.read_csv('data/integrated/all_flt.csv')
    analyzer = DataAnalyzer(gc_file, flt_file, gc_sample, flt_sample)
    analyzer.df_list_exploration()
    analyzer.analyze_means_exploration(gc_df=gc_df, flt_df=flt_df)
    
    pass
        
        
        
        