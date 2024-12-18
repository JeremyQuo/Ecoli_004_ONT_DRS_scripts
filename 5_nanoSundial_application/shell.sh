# find the rna name and type for each region
bedtools intersect -a merged_positive_region.bed -b gene.bed -wb -s >merged_positive_region_gene.bed
bedtools intersect -a merged_positive_region.bed -b cds_expand.bed -wb -s >merged_positive_region_gene_expand.bed

# find TU information
python extract_TU.py --bam genomic.bam --bed gene.bed --is_drs --cutoff 10 --output TU_result