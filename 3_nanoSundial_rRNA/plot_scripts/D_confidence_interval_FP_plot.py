import scipy.stats as stats
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
from sklearn import metrics
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
folder_list=[10,20,50,100,200,500,1000,2000]
folder_list = [str(x) for x in folder_list]

result_df = pd.DataFrame()
result_dict={}
for item in folder_list:
    df_path = 'xxxxx/' + item + '/nanoSundial_manova.csv'
    df = pd.read_csv(df_path)
    df = df[['mean_differ']]
    df['Coverage'] = item
    result_df = pd.concat([result_df, df])
    # 设定均值和标准差
    mean = round(df['mean_differ'].mean(),3)
    std_dev = round(df['mean_differ'].std(),3)
    # print(stats.shapiro(df['mean_differ']))
    # 计算置信区间
    confidence_level = 0.9999
    alpha = 1 - confidence_level

    # 使用norm.ppf函数计算Z分数
    z_score = stats.norm.ppf(1 - alpha / 2)

    # 计算置信区间
    lower_bound = mean - z_score * std_dev
    upper_bound = mean + z_score * std_dev
    # print(str(item)+','+str(upper_bound))
    result_dict[int(item)]=std_dev
print(result_dict)
folder_list.reverse()
category = pd.api.types.CategoricalDtype(categories=folder_list, ordered=True)
result_df['Coverage'] = result_df['Coverage'].astype(category)
pp = p9.ggplot(result_df, p9.aes(x='mean_differ',color='Coverage')) \
    + p9.geom_density()\
    + p9.theme_bw() \
    + p9.xlim(-0.25,0.25)\
    + p9.theme(
    figure_size=(4,4),
    axis_text=p9.element_text(size=12,family='Arial'),
    axis_title=p9.element_text(size=12,family='Arial'),
    panel_grid_minor=p9.element_blank(),
    title=p9.element_text(size=12,family='Arial'),
    strip_background=p9.element_rect(alpha=0),
    strip_text=p9.element_text(size=12,family='Arial'),
    legend_position='bottom'
    )

print(pp)
pp.save('coverage_negative.pdf',dpi=300)
print(1)