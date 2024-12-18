import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

gap_len = 100
def calculate_distance(line):
    start = line[7]+gap_len
    end = line[8]-gap_len
    if line[1]>= start and line[1]<= end:
        distance = (line[1]-start)/(end-start) * 100
    elif line[1]<=start:
        distance = line[1]-start
    else:
        distance = line[1] - end +100
    if line[5] == "-":
        distance = 100-distance
    return distance

# positive_mod_U_positive_cds_expand is the filtered sites by gene.bed on CDS region.
# Then expand 100bps both of start and end.
df = pd.read_csv('positive_mod_U_positive_cds_expand.bed', sep='\t',header=None)
df = df[df[12]=='mRNA']
df.drop_duplicates(subset=[1],inplace=True)
df['distance'] = df.apply(lambda row:np.min([abs(row[1]-row[7]-gap_len),abs(row[1]-row[8]+gap_len)]), axis=1)
df.sort_values(by=[1,'distance'], inplace=True)
df.columns=[str(x) for x in list(range(0, 14))]
df.drop_duplicates(keep='first', inplace=True,subset=['1'])
df['new_distance'] = df.apply(lambda row:calculate_distance(row), axis=1)
df['new_distance'] = df['new_distance'].apply(lambda x: x/2 if x<0 else (x-100)/2+100 if x>100 else x)


pp= p9.ggplot(df, p9.aes(x='new_distance'))\
    + p9.geom_density()\
    + p9.labs(x='')\
    + p9.theme_bw() \
    + p9.xlim(-50,150)\
    + p9.theme(
        figure_size=(4.5,3),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        panel_grid_minor=p9.element_blank(),
        title=p9.element_text(size=12,family='Arial'),
        strip_background=p9.element_rect(alpha=0),
        strip_text=p9.element_text(size=12,family='Arial')
            )
print(pp)
pp.save("pseU_density_point.pdf",dpi=300)