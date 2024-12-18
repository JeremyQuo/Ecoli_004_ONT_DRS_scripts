import pandas as pd
import plotnine as p9
import numpy as np
from matplotlib import pyplot as plt
import math
plt.rcParams['pdf.fonttype'] = 42
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

result_df = pd.DataFrame()
for key,value in path_dict.items():
    input = input_path + key + '/results/observed_quality/final_observed_accuracy.txt'
    temp = pd.read_csv(input,sep='\t').sample(frac=0.1)
    temp['group'] = value
    result_df = pd.concat([result_df,temp],axis=0)
    print(1)

result_df.drop('ID',axis=1,inplace=True)
print(1)


def draw_observed_feature(df):
    df["p_ins"] =  df["Ins"] / (df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
    df["p_del"] = df["Del"] / (df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
    df["p_sub"] = df["Sub"] / (df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
    new_df = df[["p_ins",'p_del','p_sub','Iden',"Acc",'group']]
    new_df.columns = ['Insertion','Deletion','Substitution','Identity',"Accuracy","Samples"]
    # new_df['Q_value'] = new_df['Accuracy'].apply(lambda x:acc2qvalue(x))
    new_df = pd.melt(new_df, id_vars=['Samples'], value_vars=['Insertion','Deletion','Substitution', "Accuracy"])
    # Melt and add sample info
    df_group = new_df.groupby('variable')
    df = pd.DataFrame()
    for key, value in df_group:
        if key == 'Insertion' or key == 'Deletion' or key == 'Substitution':
            print(value.shape)
            value = value[value['value'] < 0.1]
        else:
            value = value[value['value'] > 0.8]
        df = pd.concat([df, value], axis=0)
    category_data=["Accuracy",'Insertion','Deletion','Substitution']
    category = pd.api.types.CategoricalDtype(categories=category_data, ordered=True)
    df['variable'] = df['variable'].astype(category)

    name_list = ['ss&rd_004', 'IVT_neg_004']
    color_list =[ '#FF8884', '#2878B5']
    name_list.reverse()
    color_list.reverse()
    category = pd.api.types.CategoricalDtype(categories=name_list, ordered=True)
    df['Samples'] = df['Samples'].astype(category)
    plot = p9.ggplot(df, p9.aes(x='Samples', y="value",fill='Samples')) \
           + p9.theme_bw() \
           + p9.labs(y="", x='') \
           + p9.scale_fill_manual(values=color_list) \
           + p9.theme(
        figure_size=(10, 3),
        axis_text = p9.element_text(size=12, family="Arial"),
        axis_title = p9.element_text(size=12, family="Arial"),
        panel_grid_minor = p9.element_blank(),
        title = p9.element_text(size=12, family="Arial"),
        legend_text = p9.element_text(size=8, family="Arial"),
        strip_background = p9.element_rect(alpha=0),
        strip_text = p9.element_text(size=12, family="Arial"),
        legend_position='none'
    ) \
           + p9.facet_wrap('~ variable', ncol=5, scales='free_x')
    plot = plot + p9.geom_violin(color='none', width=1, style='right')
    plot = plot + p9.geom_boxplot(outlier_shape='', width=0.1)
    plot = plot + p9.coord_flip()
    # result_df.to_csv("/t1/zhguo/Data/Ecoli_RNA/count_information/Q_value.csv",index=None)
    print(plot)
    # print(1)
    plot.save('observed.pdf',dpi=300)
draw_observed_feature(result_df)