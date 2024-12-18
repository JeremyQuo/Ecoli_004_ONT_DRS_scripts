import pandas as pd
import plotnine as p9
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

gene_list="""rnr
yggX
infC
"""

def plot_tu(input_gene_name):
    # TU_coverage_list.csv is generated from https://github.com/JeremyQuo/ONT_DRS_bacteria_script
    df = pd.read_csv('TU_result/TU_coverage_list.csv')
    df= df[df['TU'].str.contains(input_gene_name)]
    df = df[df['coverage']>=10]
    sorted_df = df.sort_values('coverage', ascending=True).reset_index(drop=True)
    sorted_df['ratio'] = sorted_df['coverage'] / sorted_df['coverage'].sum()

    result_list=[]
    for idx,line in sorted_df.iterrows():
        items = line['TU']
        items = items.split('|')
        for item in items:
            if item == input_gene_name:
                result_list.append([idx,item,str(line['coverage']),line['ratio']])
            else:
                result_list.append([idx, item, str(line['coverage']), None])

    df = pd.DataFrame(result_list)
    df.columns=['idx','gene','coverage','ratio']
    gene_bed = "gene.bed"
    gene_bed = pd.read_csv(gene_bed,sep='\t',header=None)
    gene_bed = gene_bed[[3,1,2,5]]
    gene_bed.columns=['gene','start','end','strand']
    df = pd.merge(df,gene_bed,how='left',on='gene')
    colors=['#FFFFFF','#B95D3A']
    df['ymax'] = df['idx'] + 0.25
    df['ymin'] = df['ymax'] - 0.50

    result_list=[]
    group_df = df.groupby(['idx'])
    max_operon=set()
    for key,item in group_df:
        if item.shape[0] >1:
            result_list.append([key[0], item['start'].min()])
            result_list.append([key[0], item['end'].max()])
        max_operon=max_operon.union(set(item['gene'].values.tolist()))
    max_operon = pd.DataFrame(max_operon)
    max_operon.columns=['gene']
    max_operon = pd.merge(df,max_operon,how='left',on='gene')
    max_operon['idx'] = sorted_df.shape[0]
    max_operon['x_pos'] = (max_operon['start'] + max_operon['end'])/2

    min_x = df['start'].min()
    max_x = df['end'].max()
    result_list.append([sorted_df.shape[0],min_x])
    result_list.append([sorted_df.shape[0], max_x])
    line_df = pd.DataFrame(result_list)
    line_df.columns = ['idx','range']
    x_list= list(sorted_df.index)
    y_list=sorted_df['coverage'].values.tolist()

    plot = p9.ggplot()\
           + p9.geom_line(data=line_df,mapping=p9.aes(x='range',y='idx',group='idx'),size=1)\
           + p9.geom_rect(data=df,color='black',mapping=p9.aes(xmin='start', ymin='ymin',ymax='ymax', xmax='end',fill='ratio')) \
           + p9.geom_label(data=max_operon,mapping=p9.aes(x='x_pos',y='idx',label='gene'))\
           + p9.theme_bw()\
           + p9.xlim(min_x-10,max_x+10)\
           + p9.scale_y_continuous(breaks=x_list, labels=y_list)\
           + p9.labs(x='',y='Detected transcripts units')\
           + p9.scale_fill_gradientn(colors=colors)\
           + p9.theme(
            figure_size=(8, 4),
            axis_text=p9.element_text(size=12, family='Arial'),
            axis_title=p9.element_text(size=12, family='Arial'),
            panel_grid_minor=p9.element_blank(),
            title=p9.element_text(size=12, family='Arial'),
            strip_background=p9.element_rect(alpha=0),
            strip_text=p9.element_text(size=12, family='Arial'),
            legend_position='bottom',
            legend_text=p9.element_text(size=8, family='Arial'),
            legend_title=p9.element_blank(),
            legend_key_size=10
            )

    plot.save(input_gene_name+"_operon.pdf",dpi=300)
gene_list = gene_list.split('\n')
for item in gene_list:
    plot_tu(item)
