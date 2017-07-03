import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
from scipy.stats import ttest_ind
from scipy.stats import ranksums
parser = argparse.ArgumentParser()
parser.add_argument("--gene_id",help="gene_id")

Gene_ID = parser.parse_args().gene_id

os.system("head -n1 total_Exp_Matrix.txt > head.txt")
os.system("grep %s total_Exp_Matrix.txt > %s.Exp.txt" % (Gene_ID, Gene_ID))
os.system("grep %s total_Gene_genotype_matrix.txt > %s.Genotype.txt" % (Gene_ID, Gene_ID))
os.system("cat head.txt %s.Exp.txt %s.Genotype.txt > %s.txt" % (Gene_ID, Gene_ID, Gene_ID))

df_Gene = pd.read_table(Gene_ID+".txt",index_col=0)
df_Gene = df_Gene[df_Gene.index==Gene_ID]
df_Gene = df_Gene.T
df_Gene.columns = ["exp","genotype"]
df_Gene = df_Gene[df_Gene["exp"]>0.00000001]
df_Gene["exp"] = np.log10(df_Gene["exp"])
ttest_output = ttest_ind(df_Gene[df_Gene["genotype"]==12]["exp"], df_Gene[df_Gene["genotype"]==0]["exp"], equal_var = False)
print ttest_output
ranksums_output = ranksums(df_Gene[df_Gene["genotype"]==12]["exp"], df_Gene[df_Gene["genotype"]==0]["exp"])
print ranksums_output
print "Expression comparison:",np.mean(df_Gene[df_Gene["genotype"]==12]["exp"]), np.mean(df_Gene[df_Gene["genotype"]==0]["exp"])
sns.boxplot(x="genotype",y="exp",data=df_Gene)
sns.stripplot(x="genotype", y="exp", data=df_Gene, jitter=True,color="grey")
plt.title(Gene_ID)
plt.show()
os.system("rm %s.*" % (Gene_ID))

# cd /f/hang/Data/GeneFusion/Database/ICGCdataset/New_Selected_project_output/ALL-Project
#python ../../Pipeline_Current/Statistics/associating_genotype_with_expression.py --gene_id SPEN