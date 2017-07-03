import pandas as pd
import argparse
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument('--gene1',help='gene_id 1')
parser.add_argument('--gene2',help='gene_id 2')
Gene_ID_1 = parser.parse_args().gene1
Gene_ID_2 = parser.parse_args().gene2

os.system('head -n1 total_Exp_Matrix.txt > head.txt')
os.system('grep -P "^%s$(printf \'\t\')" total_Exp_Matrix.txt > %s.Exp.txt' % (Gene_ID_1, Gene_ID_1))
os.system('grep -P "^%s$(printf \'\t\')" total_Gene_genotype_matrix.txt > %s.Genotype.txt' % (Gene_ID_1, Gene_ID_1))
os.system('cat head.txt %s.Exp.txt %s.Genotype.txt > %s.txt' % (Gene_ID_1, Gene_ID_1, Gene_ID_1))

os.system('grep -P "^%s$(printf \'\t\')" total_Exp_Matrix.txt > %s.Exp.txt' % (Gene_ID_2, Gene_ID_2))
os.system('grep -P "^%s$(printf \'\t\')" total_Gene_genotype_matrix.txt > %s.Genotype.txt' % (Gene_ID_2, Gene_ID_2))
os.system('cat head.txt %s.Exp.txt %s.Genotype.txt > %s.txt' % (Gene_ID_2, Gene_ID_2, Gene_ID_2))

df_Gene1 = pd.read_table(Gene_ID_1+'.txt',index_col=0)
# df_Gene1 = df_Gene1[df_Gene1.index==Gene_ID_1]
df_Gene1 = df_Gene1.T
df_Gene1.columns = ['gene1_exp','gene1_genotype']
df_Gene1['gene1_exp'] = np.log10(df_Gene1['gene1_exp']+0.00000001)


df_Gene2 = pd.read_table(Gene_ID_2+'.txt',index_col=0)
# df_Gene2 = df_Gene2[df_Gene2.index==Gene_ID_2]
df_Gene2 = df_Gene2.T
df_Gene2.columns = ['gene2_exp','gene2_genotype']
df_Gene2['gene2_exp'] = np.log10(df_Gene2['gene2_exp']+0.00000001)


df_merge = pd.concat([df_Gene1, df_Gene2],axis=1)
sns.pairplot(df_merge,hue='gene1_genotype')
plt.show()

os.system('rm %s.*' % (Gene_ID_1))
os.system('rm %s.*' % (Gene_ID_2))
os.system('rm head.txt')

# cd /f/hang/Data/GeneFusion/Database/ICGCdataset/New_Selected_project_output/ALL-Project
# python ../../Pipeline_Current/Statistics/Associate_Two_Genes_considering_genotype.py --gene1 ETS2 --gene2 ERG