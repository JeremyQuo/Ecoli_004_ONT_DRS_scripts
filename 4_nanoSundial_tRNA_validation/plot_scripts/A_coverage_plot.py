import pandas as pd
import numpy as np
# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
gene_bed = "gene.bed"
gene_bed = pd.read_csv(gene_bed,sep='\t',header=None)

gene_bed = gene_bed.sort_values(1)
input_path="xxxxx/"

path_dict={
    'new_IVT':"IVT_neg_004",
    'new_WT':"ss&rd_004",
}

df=pd.DataFrame()
for key,value in path_dict.items():
    coverage_file = pd.read_csv(input_path + key + '/coverage.txt', sep='\t', header=None)
    for idx,line in gene_bed.iterrows():
        start_pos = line[1]
        end_pos = line[2]
        type = line[6]

        coverage_file.loc[(coverage_file[1] >= start_pos) & (coverage_file[1] <= end_pos), 4] = type
    coverage_file.dropna(axis=0,inplace=True)
    coverage_file.drop([0, 1], axis=1,inplace=True)
    coverage_file[5]=value
    df=pd.concat([df,coverage_file],axis=0)
df.to_csv('coverage.txt',header=False,index=False)