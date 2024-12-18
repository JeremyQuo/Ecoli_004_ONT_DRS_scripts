import pandas as pd
import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
from statsmodels.stats.multitest import multipletests
def calculate_distance(line,positive_site):
    site_list = positive_site[line['chrom']]
    distance_list = np.array(site_list)-line['position']
    distance_list = np.abs(distance_list)
    return np.min(distance_list)



actual_positives_23s = [2503, 2504,1911,1915,1917,1939,2498,1835,2552,2030,2604,2251,746,745,747,1618,2445,2449,2069,1962,955,2580,2457,2605, 2501]
actual_positives_16s = [1407,516,1498,1402,967,966,1516,1518,1519,527,1207]
actual_positives_16s = ['J01859.1:' + str(item) for item in actual_positives_16s]
actual_positives_23s = ['NR_103073.1:' + str(item) for item in actual_positives_23s]

actual_positives_16s.extend(actual_positives_23s)
actual_positives = actual_positives_16s


coverage_list=['50','100','200','all']
result_list=[]
for coverage in coverage_list:
    df = pd.read_csv('xxx/sundial_rRNA_manova_'+coverage+'.csv')
    df['-log10(P-value)'] = -np.log10(df['pvalue'])
    mean_cutoff = 0.18
    shift_list = [0,1,2,3,4,5,6]


    for p_cutoff in range(0,6):
        positive_df = df[
            (df['dwell_differ'].abs() >= 1) & (df['mean_differ'].abs() > mean_cutoff) & (df['-log10(P-value)'] > p_cutoff)]

        positive_df['position'] = positive_df['position'] + 1
        positive_df['f_position'] = positive_df['chrom'] + ':' + positive_df['position'].astype(str)

        predicted_site = positive_df['f_position'].values

        for shift in shift_list:
            predicted_site_set = set()
            for shift_item in range(-shift,shift+1):
                for item in predicted_site:

                    item_list = item.split(':')
                    item_list[1] = str(int(item_list[1]) + shift_item)
                    predicted_site_set.add(':'.join(item_list))
            true_positives = set(predicted_site_set) & set(actual_positives)

            # 计算召回率（Recall）
            precision = len(true_positives) / len(predicted_site)
            precision = np.min([precision,1])
            recall = len(true_positives) / len(actual_positives)
            F1 = 2 * (precision * recall) / (precision + recall)
            result_list.append([shift, precision, recall, F1,str(p_cutoff),coverage])

    vis_col='f1 score'
    result_df = pd.DataFrame(result_list,columns=['shift','precision','recall','f1 score','P-value','coverage'])
    df_melt=pd.melt(result_df,id_vars=['shift','P-value','coverage'],value_vars=['precision','recall','f1 score'])
# result_df =pd.melt(result_df,id_vars=['shift'],value_vars=['precision','recall','f1 score'])
category = pd.api.types.CategoricalDtype(categories=coverage_list, ordered=True)
df_melt['coverage'] = df_melt['coverage'].astype(category)
pp= p9.ggplot(df_melt, p9.aes(x='shift',y='value',color='P-value',group='P-value'))\
    + p9.ylim(0,1)\
    + p9.scale_color_manual(values={'0':'#000000',
                                   '1':'#412728',
                                   '2':'#7f4d3e',
                                   '3':'#b87c4c',
                                   '4':'#e2b659',
                                   '5':'#f9f871',
                                    '6':'#f6f2cb'
                                   })\
    + p9.geom_line(alpha=0.75)\
    + p9.geom_point(alpha=1)\
    + p9.theme_bw() \
    + p9.theme(
        figure_size=(4.5,6),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom'
            )\
    + p9.facet_grid('coverage ~ variable ')

print(pp)
pp.save('pvalue.pdf',dpi=300)
print(1)
