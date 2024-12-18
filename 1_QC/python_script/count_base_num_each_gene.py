import pysam
import pandas as pd
import numpy as np

# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------

def intersection(a, b):
    start = max(a[0], b[0])
    end = min(a[1], b[1])
    if start <= end:
        return end - start
    else:
        return 0
gene_bed = "gene.bed"
gene_bed = pd.read_csv(gene_bed,sep='\t',header=None)

gene_bed = gene_bed.sort_values(1)

# 合并区间
merged_data = []
start = gene_bed[1].iloc[0]
end = gene_bed[2].iloc[0]

for i in range(1, len(gene_bed)):
    if gene_bed[1].iloc[i] <= end:
        end = max(end, gene_bed[2].iloc[i])
    else:
        merged_data.append([gene_bed[0].iloc[i-1], start, end,gene_bed[3].iloc[i-1],gene_bed[4].iloc[i-1],gene_bed[5].iloc[i-1],gene_bed[6].iloc[i-1]])
        start = gene_bed[1].iloc[i]
        end = gene_bed[2].iloc[i]

merged_data.append([gene_bed[0].iloc[-1], start, end,gene_bed[3].iloc[i],gene_bed[4].iloc[i],gene_bed[5].iloc[i],gene_bed[6].iloc[i]])

gene_df = pd.DataFrame(merged_data)

input_path="xxxxx/"

path_dict={
    'new_IVT':"IVT_neg_004",
    'new_WT':"ss&rd_004",
    'new_WT_rep':"ss&rd_004_rep",
}

result_df=pd.DataFrame()
for key,value in path_dict.items():

    print(key)
    bam_file=input_path+key+"/genomic.bam"

    bam_file= pysam.AlignmentFile(bam_file,'rb')
    result_dict={}
    gene_dict={}
    for idx,line in gene_df.iterrows():

        chrom = line[0]
        start = line[1]
        end = line[2]
        strand = line[5]
        gene_name = line[3]
        gene_interval = [start, end]
        for read in bam_file.fetch(chrom,start,end):
            if read.is_secondary or read.is_supplementary:
                continue
            if (strand == '+' and read.is_reverse) or (strand == '-' and not read.is_reverse):
                continue
            if gene_name not in result_dict:
                result_dict[gene_name] = 0
                gene_dict[gene_name] = 0
            pairs = pd.DataFrame(read.aligned_pairs)

            pairs = pairs[(pairs[1]<=end)&(pairs[1]>=start)]
            pairs = pairs.dropna(axis=0)
            result_dict[gene_name] = result_dict[gene_name] + pairs.shape[0]

    df = pd.DataFrame.from_dict(result_dict, orient='index')
    df.reset_index(inplace=True)
    df.columns=['gene','coverage_num']
    temp_df = gene_bed[[3,6]]
    temp_df.columns=['gene','type']
    # temp_df['type'] = temp_df['type'].apply(lambda x: 'ncRNA' if x == 'tRNA' else x)
    df = pd.merge(df,temp_df,on='gene',how='left')
    df.drop('gene',axis=1,inplace=True)
    grouped = df.groupby('type').sum()
    grouped.reset_index(inplace=True)
    rna_coverage = grouped['coverage_num'].sum()
    print(grouped)




