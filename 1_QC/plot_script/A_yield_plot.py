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
df = pd.read_csv("aligned_num.csv")
df['aligned_read'] = round(df['aligned_read'] / 1000, 2)
df['not_aligned_read'] = round(df['not_aligned_read'] / 1000, 2)
df['aligned_base'] = round(df['aligned_base'] / 1000000, 2)
df['not_aligned_base'] = round(df['not_aligned_base'] / 1000000, 2)
df.drop(['total_read','total_base'],axis=1,inplace=True)
df = pd.melt(df, id_vars=['group'], value_vars=['aligned_read', 'not_aligned_read', 'aligned_base', 'not_aligned_base'])
df['Sample'] = df.apply(lambda x:x['group'] if 'not' not in x['variable'] else 'Unmap',axis=1)
df['Type'] = df.apply(lambda x:'Reads num (K)' if 'base' not in x['variable'] else 'Bases num (Mb)',axis=1)
# name_list=[ "Total RNA","rRNA dep-1","rRNA dep-2","rRNA dep-3","rRNA dep+SS",'IVT_neg',"IVT_pos"]

name_list=['ss&rd_004','IVT_neg_004','ss&rd_002','IVT_neg_002','Unmap']

# color_list = [ '#ffd988', '#FFF2EC', '#ffc6b1', '#FFA47A', '#FDFB7D', '#d8ecec', '#ffc2d1']
# name_list=[ "E Total RNA","E rd RNA","E rd+ss RNA","S rd+ss RNA"]
color_list = [ '#FF8884', '#2878B5','#f8ac8c', '#9Ac9DB', '#FFFFFF']
name_list.reverse()
color_list.reverse()
category = pd.api.types.CategoricalDtype(categories=name_list, ordered=True)
df['Sample'] = df['Sample'].astype(category)
category = pd.api.types.CategoricalDtype(categories=name_list, ordered=True)
df['group'] = df['group'].astype(category)
category = pd.api.types.CategoricalDtype(categories=['Reads num (K)','Bases num (Mb)'], ordered=True)
df['Type'] = df['Type'].astype(category)

pp= p9.ggplot(df, p9.aes(x='group',y='value',fill='Sample'))\
    + p9.geom_bar(stat='identity',width=0.6,color='black')\
    + p9.theme_bw() \
    + p9.labs(y="", x='')\
    + p9.scale_fill_manual(values=color_list)\
    + p9.theme(
        figure_size=(6,2),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='none'
            )\
    + p9.facet_wrap('~ Type', ncol=2,scales='free_x') \
    + p9.coord_flip()

pp.save("yield_info.pdf",dpi=300)
print(pp)

print(pp)
print(1)