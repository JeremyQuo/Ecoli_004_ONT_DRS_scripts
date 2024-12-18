import pandas as pd
import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
from sklearn import metrics
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
from statsmodels.stats.multitest import multipletests

# tRNA_sites.csv is the positive sites on tRNA from nanoSundial
prd_df = pd.read_csv('tRNA_sites.csv',header=None)
prd_df.columns=['index','Result']

# tRNA_mod.csv is the file from MODOMICS
tp_df = pd.read_csv('tRNA_mod.csv',header=None,sep='\t').reset_index()
tp_df = tp_df[['index',7]]
tp_df.columns=['index','Type']
result_df = pd.merge(prd_df,tp_df,on='index')

oder_seq=result_df.value_counts(subset='Type').reset_index()['Type'][:13]
result_df = pd.merge(result_df,pd.DataFrame(oder_seq,columns=['Type']),on='Type')

category = pd.api.types.CategoricalDtype(categories=oder_seq, ordered=True)
result_df['Type'] = result_df['Type'].astype(category)

category = pd.api.types.CategoricalDtype(categories=['low coverage','not detected','detected'], ordered=True)
result_df['Result'] = result_df['Result'].astype(category)

pp = p9.ggplot(result_df, p9.aes(x='Type',fill='Result'))\
    + p9.geom_bar(stat = "count",width=0.5)\
    + p9.scale_fill_manual(values={
    'detected':"#D1D1F8",
    'not detected':'#BEFCFE',
    'low coverage':"#D7D7D7"
})\
    + p9.theme_bw() \
    + p9.theme(
        figure_size=(8,4),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )\

print(pp)
pp.save('tRNA_overall.pdf',dpi=300)
print(1)
