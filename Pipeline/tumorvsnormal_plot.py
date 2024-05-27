import matplotlib.pyplot as plt
import pandas as pd





df = pd.read_csv('copynum_per_sample_resum.txt',sep=',')
print(df.columns)

df['pct'] = df['epithelial_and_tumoral']/df['total_epithelial']
print(df)

df.drop(columns=['total_epithelial','epithelial_and_tumoral'  ,'total_tumoral'])
pt=pd.DataFrame()
pt['sample']=df['sample']
pt['tumoral'] =df['pct']*100
pt['epithelials']=100-pt['tumoral']
pt.plot(x='sample',kind='bar', stacked=True,title='Tumoral ratio on epithelials')

print(pt)

plt.ylabel('Proportion of tumoral cells')
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
plt.savefig('Stacked_tumoralcell_barplot.png')
plt.show()