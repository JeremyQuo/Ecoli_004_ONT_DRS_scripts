import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.neighbors import KernelDensity
import numpy as np
#
from scipy.stats import spearmanr

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

axis_x='ss_rd_004'
axis_y='IVT_neg_004'
df = pd.read_csv('gene.count')
df_new=pd.read_csv('gene_wt.count')

df.columns=['gene','ss_rd_004','type','length']
columns_list=['gene','IVT_neg_004','type','length']
df_new.columns=columns_list
df = pd.merge(df,df_new,on=['gene','type','length'])



df.fillna(value=1, inplace=True)
df=df[df['type']=='tRNA']

df['gene']=df['gene_name'].apply(lambda x:x[0].upper()+x[1:3])
wt_df=df.groupby('gene')['ss&rd_004'].sum().reset_index()
ivt_df=df.groupby('gene')['IVT_neg_004'].sum().reset_index()
df = pd.merge(wt_df,ivt_df,on='gene')

corr, _ = spearmanr(df['ss&rd_004'], df['IVT_neg_004'])
df['ss&rd_004'] = np.log10(df['ss&rd_004'])
df['IVT_neg_004'] = np.log10(df['IVT_neg_004'])

pp = p9.ggplot(df, p9.aes(x='ss&rd_004',y='IVT_neg_004',fill='gene')) \
 + p9.ylim(2, 6) \
 + p9.xlim(2, 6) \
 + p9.theme_bw() \
 + p9.geom_point(size=3,stroke=0.3,color='none',alpha=0.8) \
 + p9.theme(
            figure_size=(4, 3),
            legend_key_width=8,
            legend_position='right',
            axis_text=p9.element_text(size=12,family='Arial'),
            axis_title=p9.element_text(size=12,family='Arial'),
            panel_grid_minor=p9.element_blank(),
            title=p9.element_text(size=12,family='Arial'),
            legend_text=p9.element_text(size=8,family='Arial'),
            strip_background=p9.element_rect(alpha=0),
            strip_text=p9.element_text(size=12,family='Arial'),
                      )
print(pp)
pp.save('tRNA_dotplot.pdf',dpi=300)
print(1)