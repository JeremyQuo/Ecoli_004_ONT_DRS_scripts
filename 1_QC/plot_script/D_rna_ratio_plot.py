import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
df = pd.read_csv("../genomic_cov.count")
df['coverage_num'] = round(df['coverage_num'] / 1000000, 1)
df['group'] = df['group'].replace({'ss&rd_002_2':'ss&rd_002'})
df =df[df['group']!='ss&rd_002_1']
name_list=['ss&rd_004','IVT_neg_004','ss&rd_002','IVT_neg_002']
rna_list=['tRNA','ncRNA','rRNA','mRNA','Others']
rna_list.reverse()
color_list = [ '#f8ac8c', '#f8ac8c', '#9Ac9DB', '#FF8884', '#2878B5']
name_list.reverse()
color_list.reverse()
category = pd.api.types.CategoricalDtype(categories=name_list, ordered=True)
df['group'] = df['group'].astype(category)

category = pd.api.types.CategoricalDtype(categories=rna_list, ordered=True)
df['type'] = df['type'].astype(category)
pp= p9.ggplot(df, p9.aes(x='group',y='coverage_num',fill='type'))\
    + p9.geom_bar(stat='identity',width=0.6,color='black',position='fill')\
    + p9.geom_text(p9.aes(label='coverage_num'),position=p9.position_fill(0.5),size=9)\
    + p9.theme_bw() \
    + p9.labs(y="", x='')\
    + p9.scale_fill_manual(values={'tRNA':'#CDFADB',
                               'ncRNA':'#F6FDC3',
                               'mRNA':'#FFC0B7',
                               'rRNA':'#FFCF96',
                               'Others':'#E3E3E3'})\
    + p9.theme(
        figure_size=(6,3),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )\
    + p9.coord_flip()


print(pp)
pp.save("rna_number.pdf",dpi=300)
print(pp)
print(1)