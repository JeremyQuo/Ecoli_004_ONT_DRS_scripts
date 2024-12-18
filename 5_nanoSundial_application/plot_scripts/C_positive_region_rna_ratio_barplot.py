import pandas as pd
import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

def read_result_gene_bed(file_name,label,count_num):

    df = pd.read_csv( file_name   ,sep='\t',header=None)
    df.sort_values(by=[1,12],inplace=True)
    df.drop_duplicates(subset=[1],inplace=True,keep='last')
    temp = df[12].value_counts().reset_index()
    new_row = {12: 'Others', 'count': count_num-temp['count'].sum()}
    temp.loc[len(temp)] = new_row
    temp['label'] = label
    return temp
df_1 = read_result_gene_bed('merged_positive_region_rep1_gene.bed','rep1',1157)
df_2 = read_result_gene_bed('merged_positive_region_rep2_gene.bed','rep2',1797)

df = pd.concat([df_1,df_2])
df.columns=['type','count','label']
pp=  p9.ggplot(df, p9.aes(x='label',y='count',fill='type'))\
    + p9.geom_bar(stat='identity',width=0.5,color='black')\
    + p9.geom_text(p9.aes(label='count'),position=p9.position_stack,size=12)\
    + p9.theme_bw() \
    + p9.labs(y="", x='')\
    + p9.scale_fill_manual(values={'tRNA':'#CDFADB',
                               'ncRNA':'#F6FDC3',
                               'mRNA':'#FFC0B7',
                               'rRNA':'#FFCF96',
                               'Others':'#E3E3E3'})\
    + p9.theme(
        figure_size=(3,8),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )
pp.save('rna_ratio_rep.pdf',dpi=300)
print(pp)