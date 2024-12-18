
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
from sklearn import metrics
import numpy as np
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
def read_fasta_to_dic(filename):
    """
    function used to parser small fasta
    still effective for genome level file
    """
    fa_dic = OrderedDict()

    with open(filename, "r") as f:
        for n, line in enumerate(f.readlines()):
            if line.startswith(">"):
                if n > 0:
                    fa_dic[short_name] = "".join(seq_l)  # store previous one

                full_name = line.strip().replace(">", "")
                short_name = full_name.split(" ")[0]
                seq_l = []
            else:  # collect the seq lines
                if len(line) > 8:  # min for fasta file is usually larger than 8
                    seq_line1 = line.strip()
                    seq_l.append(seq_line1)

        fa_dic[short_name] = "".join(seq_l)  # store the last one
    return fa_dic

df_23s = pd.read_csv('23s_sundial_manova.csv')
actual_positives_23s = [2503, 2504, 1911, 1915, 1917, 1939, 2498, 1835, 2552, 2030, 2604, 2251, 746, 745, 747, 1618,
                            2445, 2449, 2069, 1962, 955, 2580, 2457, 2605, 2501]
temp_df = pd.DataFrame(actual_positives_23s,columns=['position'])
fasta=read_fasta_to_dic("23S_rRNA.fasta")
print(1)
# wga_result = pd.read_csv('positive_lr_WGA.csv')
#
# wt_result = pd.merge(wt_result,wga_result[['chrom','position','position_1']],how='left',on=['chrom','position'])
# wt_result = wt_result[wt_result['position_1_y'].isna()]
df_23s['-log10(P-value)'] = -np.log10(df_23s['adj_p'])

visual_col='-log10(P-value)'

df_23s[visual_col] = df_23s[visual_col].abs()

df_23s['position'] = df_23s['position']+1
temp_df = pd.merge(df_23s,temp_df,on='position',how='inner')
temp_df['new_pos']=0
temp_df['new_differ']=0
for idx,item in temp_df.iterrows():
    position = item['position']
    temp = df_23s[(df_23s['position']>=position-4)&(df_23s['position']<=position+4)]
    temp_df.loc[idx,'new_pos'] = temp.loc[temp[visual_col].idxmax()]['position']
    temp_df.loc[idx, 'new_differ'] = temp.loc[temp[visual_col].idxmax()][visual_col]
base_dict={'A':'A',
           'T':'U',
           'C':'C',
           'G':'G',}
# + p9.geom_text(data=temp_df, mapping=p9.aes(x='new_pos', y='new_differ', label='name')) \
#  \
temp_df['name']=temp_df.apply(lambda x:base_dict[fasta[x['chrom']][x['position']-1]]+str(x['position']),axis=1)
df_23s['prediction'] = df_23s.apply(lambda x:'Positive' if x['-log10(P-value)'] > 3 and x['mean_differ'] else 'Negative',axis=1)
pp = p9.ggplot(df_23s, p9.aes(x='position',y=visual_col))\
    + p9.geom_line(color='#2F4858')\
    + p9.geom_point(data=temp_df,mapping=p9.aes(x='new_pos',y='new_differ'),fill='#8B4A48')\
    + p9.geom_text(data=temp_df, mapping=p9.aes(x='new_pos', y='new_differ', label='name')) \
    + p9.theme_bw() \
    + p9.geom_hline(yintercept=0.1)\
    + p9.scale_x_continuous(breaks=[0,500,1000,1500,2000,2500,2904], limits=[0, 2904])\
    + p9.theme(
        figure_size=(8,3),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )
pp.save('mean_differ_23s.pdf',dpi=300)
print(1)
