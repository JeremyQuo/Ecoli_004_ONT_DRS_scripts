import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

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

# df.to_csv('mod_U_merged.csv',index=False)
df['differ'] = df['ss&rd_004'] - df['IVT_neg_004']

pp= p9.ggplot(df, p9.aes(x='differ'))\
    + p9.geom_histogram(binwidth=10,fill='#D7D7D7')\
    + p9.xlim(-100,100)\
    + p9.scale_y_log10()\
    + p9.theme_bw() \
    + p9.theme(
        figure_size=(4,3),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial')
            )
print(pp)
pp.save("pseU_differ_distribution.pdf",dpi=300)
print(1)