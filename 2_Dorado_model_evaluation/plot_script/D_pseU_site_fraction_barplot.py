import pandas as pd
import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

input_path="xxxxx/"

path_dict={
    'new_IVT':"IVT_neg_004",
    'new_WT':"ss&rd_004",
}


df =pd.DataFrame()
for key,value in path_dict.items():
    tem_df = pd.read_csv(input_path+key+"/basecall_pseU.bed", sep='\t',header=None)
    tem_df = tem_df[tem_df[4]>=20]
    tem_df = tem_df[[0,1,2,5,10]]
    print(tem_df[10].mean())
    tem_df.columns=[0,1,2,5,value]

    if df.shape[0]==0:
        df = tem_df
    else:
        df = pd.merge(df,tem_df,how='inner',on=[0,1,2,5])

rRNA_pos = pd.read_csv('ground_truth.csv',header=None)
rRNA_pos['pos']=rRNA_pos[1]-rRNA_pos[7]
rRNA_pos =rRNA_pos[[1,2,'pos']]
rRNA_pos.columns=['1','2','pos']

desired_values = rRNA_pos['1'].values
desired_rows = df[df['1'].isin(desired_values)]
desired_rows = pd.merge(desired_rows,rRNA_pos,on='1')
df = pd.melt(desired_rows,id_vars='pos',value_vars=['ss&rd_004','IVT_neg_004'])
df['pos'] = df['pos'].astype(str)
df=df[['pos','value','variable']]

order_list=[ 515,745, 954,1911, 1917, 2457, 2504, 2580, 2604,2605, 1210,743,1915,2028,2075, 2554]
order_list = [str(x) for x in order_list]
category = pd.api.types.CategoricalDtype(categories=order_list, ordered=True)
df['pos'] = df['pos'].astype(category)
category = pd.api.types.CategoricalDtype(categories=['ss&rd_004','IVT_neg_004'], ordered=True)
df['variable'] = df['variable'].astype(category)
pp= p9.ggplot(df, p9.aes(x='pos',y='value',fill='variable'))\
    + p9.geom_col(width=0.6,position='dodge')\
    + p9.theme_bw() \
    + p9.theme(
        figure_size=(9,3),
        axis_text=p9.element_text(size=10,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )
print(pp)
pp.save("pseU_showcase.pdf",dpi=300)
print(1)