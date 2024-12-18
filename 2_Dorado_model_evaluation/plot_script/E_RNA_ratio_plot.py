import pandas as pd
import numpy as np
import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

## positive_mod_U_positive_gene is the filtered sites and annotated the RNA by gene.bed
df = pd.read_csv('positive_mod_U_positive_gene.bed',sep='\t',header=None)

temp = df[12].value_counts().reset_index()
# 数据
sizes = temp['count'].values.tolist()
sizes[0] =sizes[0]-1
sizes[2] =sizes[2]+1
sizes.append(288-df.shape[0])
# 饼图的标签
labels = temp[12].values.tolist()
labels.append('Intergenic region')
# sizes[0]=56
# sizes[1]=8
# 饼图的颜色
colors=[ '#CDFADB','#FFCF96','#FF8080', '#F6FDC3','#D7D7D7']

def autopct_func(pct, allvals):
    absolute = round(pct/100.*sum(allvals))
    return f"{pct:.1f}%\n({absolute:d})"


# 绘制饼图
plt.pie(sizes, labels=labels,colors=colors,pctdistance=0.65,
        autopct=lambda pct: autopct_func(pct, sizes), startangle=90)
# 标题
plt.title("Gene type of pseU sites from Dorado")
plt.savefig('pseU_ratio.pdf', dpi=300)
# 显示图形
plt.show()