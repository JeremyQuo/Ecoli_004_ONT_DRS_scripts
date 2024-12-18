import numpy as np
import pandas as pd

df_a = pd.read_csv('merged_positive_region_sub1_gene.bed',sep='\t',header=None)
df_b = pd.read_csv('merged_positive_region_sub2_gene.bed',sep='\t',header=None)
def calculate_intersection(df_a,df_b):
    idx_list=[]
    for idx,line in df_a.iterrows():
        temp_line = df_b[(df_b[1]<=line[2])&(df_b[2]>=line[1])]
        if not temp_line.empty:
            idx_list.append(idx)
    return idx_list
inter_a=calculate_intersection(df_a,df_b)
inter_b=calculate_intersection(df_b,df_a)
if len(inter_a) <= len(inter_b):
    result=df_a.loc[inter_a]
else:
    result=df_b.loc[inter_b]

rna_type='mRNA'
df_a_mrna = df_a[df_a[12]==rna_type]
df_b_mrna = df_b[df_b[12]==rna_type]
result_mrna = result[result[12]==rna_type]

