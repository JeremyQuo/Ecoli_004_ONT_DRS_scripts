import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing
# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
gene_bed = "gene.bed"
# gene_bed="/t1/zhguo/Data/Saureus_data/ref/gene.bed"
gene_bed = pd.read_csv(gene_bed,sep='\t',header=None)

gene_bed = gene_bed.sort_values(1)
input_path="xxxxx/"

path_dict={
    'new_IVT':"IVT_neg_004",
    'new_WT':"ss&rd_004",
}

pbar = tqdm(total=gene_bed.shape[0]*len(path_dict.keys()), position=0, leave=True, unit='gene')
def update(*a):
    pbar.update(1)

final_df=pd.DataFrame()
cpu = 1
result_df = pd.DataFrame()
for key,value in path_dict.items():
    print(key)
    df_plus = pd.read_csv(input_path  + key + '/genomic.plus_strand.per.site.csv')
    df_minus = pd.read_csv(input_path  + key + '/genomic.minus_strand.per.site.csv')
    df = pd.concat([df_plus, df_minus],axis=0,ignore_index=True)
    df = df[df['cov']>=5]
    df['acc'] = 1 - df['mis'] - df['ins'] - df['del']
    final_temp_df = df[['acc']]
    final_temp_df.loc[:,'type']=None
    final_temp_df.loc[:,'group'] = value
    for index, line in gene_bed.iterrows():
        filtered_index = (df['strand']==line[5])&(df['pos']>=line[1])&(df['pos']<=line[2])
        final_temp_df.loc[filtered_index,'type'] = line[6]
        pbar.update(1)
    final_temp_df.dropna(axis=0,how='any',inplace=True)
    final_df = pd.concat([final_df,final_temp_df],axis=0,ignore_index=True)
    final_df.to_csv(input_path + key + '/site_acc.csv', index=False)
    print(1)
pbar.close()



