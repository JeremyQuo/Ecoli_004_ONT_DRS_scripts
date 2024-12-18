import pandas as pd
import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
from statsmodels.stats.multitest import multipletests
# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
df = pd.read_csv('sundial_rRNA_manova.csv')

df['f_position'] = df['chrom'] +':'+ df['position'].astype(str)

df['-log10(P-value)'] = -np.log10(df['adj_p'])

for shift in range(0,1):
    actual_positives_23s = set([2503, 2504, 1911, 1915, 1917, 1939, 2498, 1835, 2552, 2030, 2604, 2251, 746, 745, 747, 1618,
                            2445, 2449, 2069, 1962, 955, 2580, 2457, 2605, 2501])
    actual_positives_16s = set([1407, 516, 1498, 1402, 967, 966, 1516, 1518, 1519, 527, 1207])
    actual_positives_23s_new =set()
    actual_positives_16s_new = set()
    # for item in actual_positives_23s:
    #     for value in range(-shift,shift+1):
    #         actual_positives_23s_new.add(value+item)
    # for item in actual_positives_16s:
    #     for value in range(-shift,shift+1):
    #         actual_positives_16s_new.add(value+item)
    actual_positives_16s = ['J01859.1:' + str(item) for item in actual_positives_16s]
    actual_positives_23s = ['NR_103073.1:' + str(item) for item in actual_positives_23s]

    actual_positives_16s.extend(actual_positives_23s)
    actual_positives = actual_positives_16s
    result_list=[]
    start = 0
    step = 0.01
    stop = 0.4
    mean_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]
    start = 0
    step = 0.001
    stop = 0.04
    std_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]
    start = 0
    step = 0.5
    stop = 15
    dt_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]

    for dwell_cutoff in dt_sequence:
        for mean_cutoff in mean_sequence:

            predicted_positives = df[(df['dwell_differ'].abs() > dwell_cutoff)&(df['mean_differ'].abs() > mean_cutoff)]
            temp = predicted_positives['f_position'].values
            if len(temp) == 0:
                continue
            true_positives = set(temp) & set(actual_positives)

            # 计算精确度（Precision）
            precision = len(true_positives) / len(temp)
            if precision==0:
                continue
            # 计算召回率（Recall）
            recall = len(true_positives) / len(actual_positives)
            F1 = 2 * (precision * recall) / (precision + recall)
            result_list.append([dwell_cutoff,round(mean_cutoff,2),F1,precision,recall,'Dwell time'])
    for dwell_cutoff in std_sequence:
        for mean_cutoff in mean_sequence:

            predicted_positives = df[(df['std_differ'].abs() > dwell_cutoff)&(df['mean_differ'].abs() > mean_cutoff)]
            temp = predicted_positives['f_position'].values
            if len(temp) == 0:
                continue
            true_positives = set(temp) & set(actual_positives)

            # 计算精确度（Precision）
            precision = len(true_positives) / len(temp)
            if precision==0:
                continue
            # 计算召回率（Recall）
            recall = len(true_positives) / len(actual_positives)
            F1 = 2 * (precision * recall) / (precision + recall)
            result_list.append([dwell_cutoff,round(mean_cutoff,2),F1,precision,recall,'STD'])

result_df = pd.DataFrame(result_list,columns=['dwell_differ','mean_differ','F1','precision','recall','Type'])
pp= p9.ggplot(result_df, p9.aes(x='dwell_differ',y='mean_differ',fill='F1'))\
    + p9.geom_raster(p9.aes(fill='F1'), interpolate=True)\
    + p9.theme_bw() \
    + p9.theme(
        figure_size=(9,4.5),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
            )\
    + p9.facet_grid('~ Type',scales='free')

print(pp)
pp.save('co_f1.pdf',dpi=300)
print(1)
# df['status'] = 'Normal'