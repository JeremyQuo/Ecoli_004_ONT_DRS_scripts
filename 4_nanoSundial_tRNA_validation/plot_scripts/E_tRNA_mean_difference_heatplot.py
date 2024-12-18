import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import plotnine as p9
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
from collections import OrderedDict
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
def reverse_fasta(ref):
    base_dict = {'A': 'T', "T": "A", "C": "G", "G": "C"}
    ref = np.array(list(ref.upper()))
    for index, element in enumerate(ref):
        ref[index] = base_dict[element]
    reverse_fasta = list(reversed(''.join(ref)))
    return ''.join(reverse_fasta)

# tRNA_prediction_all.csv is the result of nanoSundial annotated by gene.bed
df = pd.read_csv('tRNA_prediction_all.csv')
result_df = pd.DataFrame()
df['name'] = df['gene'].apply(lambda x: x[:-1])
df = df[df['name']=='lys']
# df = df[df['gene']!='valU']
gene_df = pd.read_csv('gene.bed')
gene_df=gene_df[['3','ref']]
gene_df.columns=['gene','ref']
result_dict = gene_df.set_index('gene')['ref'].to_dict()

or_cutoff = 0.18
p_cutoff = 3
df.sort_values(by='position', inplace=True)
df['mean_differ']=df['mean_differ'].abs()
base_dict={'A':'A',
           'T':'U',
           'C':'C',
           'G':'G'}
df['status']=df.apply(lambda row: True if (row['-log10(fdr)']>=p_cutoff)&(abs(row['mean_differ']) > or_cutoff) else False, axis=1)
df['ref'] = df.apply(lambda row: result_dict[row['gene']][row['position']], axis=1)
df['ref'] = df['ref'].apply(lambda x: base_dict[x])
# gene_sort=['valV','valW','valT','valX','valY','valZ']
# gene_sort.reverse()
gene_sort=sorted(df['gene'].unique(),reverse=True)
df['gene'] = pd.Categorical(df['gene'], categories=gene_sort, ordered=True)
pp= p9.ggplot(df, p9.aes(x='position',y='gene',fill='mean_differ'))\
    + p9.geom_tile(color='white', size=0.1)\
    + p9.geom_text(p9.aes(label='ref'),size=8)\
    + p9.scale_fill_gradient2(low='#373A40', mid='#EEEEEE', high='#DC5F00', midpoint=0.18)\
    + p9.theme_bw() \
    + p9.labs(x='',y='')\
    + p9.theme(
        figure_size=(10,5),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial'),
        legend_position='bottom',
        legend_title=p9.element_blank(),
            )

print(pp)
pp.save(filename='tRNA_heat_val_mean.pdf',dpi=300)