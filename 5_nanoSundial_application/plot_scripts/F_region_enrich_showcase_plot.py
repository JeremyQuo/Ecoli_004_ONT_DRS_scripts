import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

gene_name= 'fkpA'
bed_file = pd.read_csv('gene.bed',header=None,sep='\t')
gene_info = bed_file[bed_file[3]==gene_name]
gene_info[1] = gene_info[1] - 100
gene_info[2] = gene_info[2] + 100
gene_info_arr = gene_info.values[0]
gene_info_arr[-1]=2
start = gene_info_arr[1]
end = gene_info_arr[2]
rep1_result = pd.read_csv('merged_positive_region_rep1_gene.bed',sep='\t',header=None)
rep1_result['group'] = 1
rep2_result = pd.read_csv('merged_positive_region_rep2_gene.bed',sep='\t',header=None)
rep2_result['group'] = 0
df = pd.concat([rep1_result,rep2_result])
df = df[(df[1]>=start)&(df[2]<=end)].values.tolist()
gene_info_arr[1] = gene_info_arr[1] + 100
gene_info_arr[2] = gene_info_arr[2] - 100
df.append(gene_info_arr)
df = pd.DataFrame(df)
df.columns=['chr','start','end','.','..','strand','group']
df['ymax'] = df['group'] + 0.25
df['ymin'] = df['ymax'] - 0.50

plot = p9.ggplot(df) \
       + p9.geom_rect(color='black',
                      mapping=p9.aes(xmin='start', ymin='ymin', ymax='ymax', xmax='end', fill='group')) \
       + p9.theme_bw() \
       + p9.xlim(start - 110, end + 110) \
       + p9.theme(
    figure_size=(8, 2),
    axis_text=p9.element_text(size=12, family='Arial'),
    axis_title=p9.element_text(size=12, family='Arial'),
    panel_grid_minor=p9.element_blank(),
    title=p9.element_text(size=12, family='Arial'),
    strip_background=p9.element_rect(alpha=0),
    strip_text=p9.element_text(size=12, family='Arial'),
    legend_position='bottom',
    legend_text=p9.element_text(size=8, family='Arial'),
    legend_title=p9.element_blank(),
    legend_key_size=10
)
plot.save(gene_name+'_enrich.pdf', dpi=300)
print(1)