import plotnine as p9
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
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
input_path="xxxxx/"

path_dict={
    'new_IVT':"IVT_neg_004",
    'new_WT':"ss&rd_004",
}

df = pd.DataFrame()
for key, value in path_dict.items():
    input_file = input_path + key + "/results/estimated_quality/final_estimated_accuracy.txt"
    tem_df = pd.read_csv(input_file, sep='\t').sample(frac=0.1)
    tem_df['Sample'] = value
    print(tem_df['Q_value'].median())
    print(value)
    tem_df.drop('ID', axis=1, inplace=True)
    df = pd.concat([df, tem_df])

df = pd.melt(df, id_vars=['Sample'], value_vars=['Q_value', 'Read_length'])
# filter
df_group = df.groupby('variable')
df = pd.DataFrame()
for key, value in df_group:
    if key == 'Read_length':
        print(value.shape)
        value = value[value['value'] < 2000]
        print(value.shape)
    else:
        print(value.shape)
        value = value[value['value'] < 30]
        print(value.shape)
    df = pd.concat([df, value], axis=0)

name_list=['ss&rd_004','IVT_neg_004']
color_list = [ '#FF8884', '#2878B5']
name_list.reverse()
color_list.reverse()
category = pd.api.types.CategoricalDtype(categories=name_list, ordered=True)
df['Sample'] = df['Sample'].astype(category)

category = pd.api.types.CategoricalDtype(categories=['Q_value', 'Read_length'], ordered=True)
df['variable'] = df['variable'].astype(category)
df.to_csv('feature.csv')
pp = p9.ggplot(df, p9.aes(x='Sample', y='value', fill='Sample')) \
     + p9.theme_bw() \
     + p9.scale_fill_manual(values=color_list) \
     + p9.labs(y="", x='') \
     + p9.theme(
    figure_size=(6, 3),
    axis_text=p9.element_text(size=12, family='Arial'),
    axis_title=p9.element_text(size=12, family='Arial'),
    panel_grid_minor=p9.element_blank(),
    title=p9.element_text(size=12, family='Arial'),
    strip_background=p9.element_rect(alpha=0),
    strip_text=p9.element_text(size=12, family='Arial'),
    legend_position='none'
) \
     + p9.facet_grid('~ variable', scales='free_x') \
     + p9.geom_violin(width=1.5, style='right',color='none') \
     + p9.geom_boxplot(width=0.1,outlier_shape='') \
     + p9.coord_flip()
print(pp)

# + p9.geom_boxplot(width=0.5 ,position=p9.position_dodge(0.7)) \

pp.save( "estimated_feature.pdf", dpi=300)
