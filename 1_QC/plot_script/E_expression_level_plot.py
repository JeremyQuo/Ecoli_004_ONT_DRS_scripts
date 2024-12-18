import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.neighbors import KernelDensity
import numpy as np

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
axis_x='ss_rd_004'
axis_y='IVT_neg_004'
df = pd.read_csv('gene.count')
df_new=pd.read_csv('gene_wt.count')

df.columns=['gene','ss_rd_004','type','length']
columns_list=['gene','IVT_neg_004','type','length']
df_new.columns=columns_list
df = pd.merge(df,df_new,on=['gene','type','length'])



df.fillna(value=1, inplace=True)
df=df[df['type']=='mRNA']
df.dropna(subset=[axis_x ,axis_y],inplace=True)
df['Gene length (Kb)'] = df['length'].apply(lambda x: 1 if x < 1000 else 2 if x<2000 else 3 if x<3000 else 4)
# df = df[(df != 0).all(axis=1)]

# CPM normalization
df[axis_x]=df[axis_x]/df[axis_x].sum()*1000000
df[axis_y]=df[axis_y]/df[axis_y].sum()*1000000

pearson_correlation = df[axis_x].corr(df[axis_y])
spearman_correlation_ivt = df[axis_x].corr(df[axis_y], method='spearman').round(3)

df[axis_x]=np.log10(df[axis_x])
df[axis_y]=np.log10(df[axis_y])


kde = KernelDensity(bandwidth=0.1).fit(df[[axis_x,axis_y]])
df['Density']=np.exp(kde.score_samples(df[[axis_x,axis_y]]))
pp = p9.ggplot(df, p9.aes(x=axis_x,y=axis_y,fill='Density'))\
 + p9.theme_bw() \
    + p9.ylim(-0.1, 5) \
    + p9.xlim(-0.1, 5) \
 + p9.geom_point(p9.aes(size='Gene length (Kb)'),stroke=0.3) \
 + p9.scale_fill_gradient(low='#ffffff', high='#1856B6')\
 + p9.theme(
            figure_size=(5,11),
            legend_key_width=8,
            legend_position='bottom',
            axis_text=p9.element_text(size=12,family='Arial'),
            axis_title=p9.element_text(size=12,family='Arial'),
            panel_grid_minor=p9.element_blank(),
            title=p9.element_text(size=12,family='Arial'),
            legend_text=p9.element_text(size=10,family='Arial')
                      )\
 +p9.facet_wrap('~ group ',ncol=1,scales='free')

print(pp)
pp.save("mRNA_expression.pdf",dpi=300)
