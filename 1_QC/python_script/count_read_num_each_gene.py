import pysam
import pandas as pd
from tqdm import tqdm

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
gene_bed['length']=gene_bed[2]-gene_bed[1]

input_path="xxxxx/"

path_dict={
    'new_IVT':"IVT_neg_004",
    'new_WT':"ss&rd_004",
    'new_WT_rep':"ss&rd_004_rep",
}

result_df=pd.DataFrame()
pbar = tqdm(total=gene_bed.shape[0]*len(path_dict.keys()), position=0, leave=True, unit='gene')

for key, value in path_dict.items():
    bam_file=input_path + key +'/'+"genomic.bam"
    bam_file= pysam.AlignmentFile(bam_file,'rb')
    result_dict={}
    for idx,line in gene_bed.iterrows():
        chrom = line[0]
        start = line[1]
        end = line[2]
        strand = line[5]
        gene_name = line[3]
        gene_interval = [start, end]
        for read in bam_file.fetch(chrom,start,end):
            if (strand == '+' and read.is_reverse) or (strand == '-' and not read.is_reverse):
                continue
            read_interval = [read.reference_start, read.reference_end]
            inter_len = intersection(gene_interval, read_interval)
            if  inter_len / (end - start) >= 0.5 or inter_len>100:
                if gene_name not in result_dict:
                    result_dict[gene_name] = 0
                result_dict[gene_name] = result_dict[gene_name]+1
        # if idx>100:
        #     break
        pbar.update(1)
    # print(item,gene_num)
    df = pd.DataFrame.from_dict(result_dict, orient='index')
    df.reset_index(inplace=True,drop=False)
    df.columns=['gene_name',value]
    if result_df.shape[0]==0:
        result_df=df
    else:
        result_df = pd.merge(result_df, df, how='outer', on='gene_name')
gene_bed = gene_bed[[3,6,'length']]
gene_bed.columns=['gene_name','type','length']
df = pd.merge(result_df,gene_bed,how='left',on='gene_name')
df.to_csv("gene.count",index=False)
