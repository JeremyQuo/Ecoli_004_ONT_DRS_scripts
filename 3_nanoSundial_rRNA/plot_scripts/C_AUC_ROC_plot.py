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
# -----------------------------------------------------------------------------
# Author: Guo Zhihao
# Created Date:2024/12
# -----------------------------------------------------------------------------
# This script is part of the project hosted at:
# https://github.com/JeremyQuo/Ecoli_004_ONT_DRS_scripts
# -----------------------------------------------------------------------------
def calculate_distance(line,positive_site):
    site_list = positive_site[line['chrom']]
    distance_list = np.array(site_list)-line['position']
    distance_list = np.abs(distance_list)
    return np.min(distance_list)
df = pd.read_csv('sundial_rRNA_manova.csv')

df['-log10(P-value)'] = -np.log10(df['adj_p'])
df['f_position'] = df['chrom']+ ':'+df['position'].astype(str)

predicted_site_set=set()
result_df = pd.DataFrame()

# true positive set
actual_positives_23s = set([2503, 2504, 1911, 1915, 1917, 1939, 2498, 1835, 2552, 2030, 2604, 2251, 746, 745, 747, 1618,
                        2445, 2449, 2069, 1962, 955, 2580, 2457, 2605, 2501])
actual_positives_16s = set([1407, 516, 1498, 1402, 967, 966, 1516, 1518, 1519, 527, 1207])
actual_positives_23s_new =set()
actual_positives_16s_new = set()
actual_positives_16s = ['J01859.1:' + str(item) for item in actual_positives_16s]
actual_positives_23s = ['NR_103073.1:' + str(item) for item in actual_positives_23s]
actual_positives_16s.extend(actual_positives_23s)
actual_positives = actual_positives_16s


result_list=[]
# mean difference roc_curve calculation
start = 0
step = 0.01
stop = 2
mean_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]
for mean_differ in mean_sequence:
    positive_df = df[ (df['mean_differ'].abs() >= mean_differ)]
    predicted_site_set = set(positive_df['f_position'].tolist())
    true_positives = set(predicted_site_set) & set(actual_positives)
    false_positive = set(predicted_site_set) - set(actual_positives)
    TPR = len(true_positives)/len(actual_positives)
    FPR = len(false_positive)/(df.shape[0]-len(actual_positives))

    result_list.append([str(mean_differ), TPR, FPR,'Mean'])

# median difference roc_curve calculation
start = 0
step = 0.01
stop = 2
median_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]
for median_differ in median_sequence:
    positive_df = df[ (df['median_differ'].abs() >= median_differ)]
    predicted_site_set = set(positive_df['f_position'].tolist())
    true_positives = set(predicted_site_set) & set(actual_positives)
    false_positive = set(predicted_site_set) - set(actual_positives)
    TPR = len(true_positives)/len(actual_positives)
    FPR = len(false_positive)/(df.shape[0]-len(actual_positives))
    result_list.append([str(median_differ), TPR, FPR,'Median'])

# std difference roc_curve calculation
start = 0
step = 0.01
stop = 2
std_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]

for std_differ in std_sequence:
    positive_df = df[ (df['std_differ'].abs() >= std_differ)]
    predicted_site_set = set(positive_df['f_position'].tolist())
    true_positives = set(predicted_site_set) & set(actual_positives)
    false_positive = set(predicted_site_set) - set(actual_positives)
    TPR = len(true_positives)/len(actual_positives)
    FPR = len(false_positive)/(df.shape[0]-len(actual_positives))
    result_list.append([str(std_differ), TPR, FPR,'STD'])

# dwell difference roc_curve calculation
start = 0
step = 1
stop = 20
dwell_sequence = [start + step * i for i in range(int((stop - start) / step) + 1)]

for dwell_differ in dwell_sequence:
    positive_df = df[ (df['dwell_differ'].abs() >= dwell_differ)]
    predicted_site_set = set(positive_df['f_position'].tolist())
    true_positives = set(predicted_site_set) & set(actual_positives)
    false_positive = set(predicted_site_set) - set(actual_positives)
    TPR = len(true_positives)/len(actual_positives)
    FPR = len(false_positive)/(df.shape[0]-len(actual_positives))
    result_list.append([str(dwell_differ), TPR, FPR,'Dwell time'])

temp_df = pd.DataFrame(result_list,columns=['shift','TPR','FPR','type'])
temp_df_1 = temp_df[temp_df['type']=='Dwell time']
auc = metrics.auc(temp_df_1['FPR'].values, temp_df_1['TPR'].values)
print(auc)
temp_df_1 = temp_df[temp_df['type']=='Mean']
auc = metrics.auc(temp_df_1['FPR'].values, temp_df_1['TPR'].values)
print(auc)
temp_df_1 = temp_df[temp_df['type']=='Median']
auc = metrics.auc(temp_df_1['FPR'].values, temp_df_1['TPR'].values)
print(auc)
temp_df_1 = temp_df[temp_df['type']=='STD']
auc = metrics.auc(temp_df_1['FPR'].values, temp_df_1['TPR'].values)
print(auc)
result_df = pd.concat([result_df,temp_df],axis=0,ignore_index=True)

pp = p9.ggplot(result_df, p9.aes(x='FPR',y='TPR',color='type'))\
    + p9.geom_line()\
    + p9.scale_color_manual(values={'Mean':'#334574',
                                    'Median':'#935b45',
                                    'Dwell time':'#fb743e',
                                    'STD':'#ffe7d5',
                                    })\
    + p9.theme_bw() \
    + p9.xlim(0,1)\
    + p9.theme(
        figure_size=(4,4.5),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )\

print(pp)
pp.save('roc_score_no_shift.pdf',dpi=300)
