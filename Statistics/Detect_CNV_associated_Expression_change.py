import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import seaborn as sns
import os
def main():
    work_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Selected_project_output/ALL-Project"
    total_Exp_Matrix = os.path.join(work_dir, "total_Exp_Matrix.txt")
    total_TAD_genotype_matrix = os.path.join(work_dir, "total_TAD_genotype_matrix.txt")
    total_Gene_genotype_matrix = os.path.join(work_dir, "total_Gene_genotype_matrix.txt")
    total_Gene_CNA_matrix = os.path.join(work_dir, "total_Gene_CNA_matrix.txt")
    significant_CNV_gene = os.path.join(work_dir, "Ttest_sig.txt")
    df_Exp = pd.read_table(total_Exp_Matrix, sep="\t", index_col=0)
    df_Gene_genotype = pd.read_table(total_Gene_genotype_matrix, sep="\t", index_col=0)
    

    with open(significant_CNV_gene, 'w') as exp_cnv:
        exp_cnv.write('gene_id\tgenotype\tmean_exp_alt\tmean_exp_ref\tt_statistics\tp_value\fold_change\n')
        for i in range(len(df_Exp)):
            if i %1000==0:
                print "count of genes: %d" % i
            sub_Exp = df_Exp.iloc[i,:]
            sub_genotype = df_Gene_genotype.iloc[i,:]
            df_Gene = pd.DataFrame({"exp":sub_Exp, "genotype":sub_genotype})
            df_Gene = df_Gene[df_Gene["exp"]>0.00000001]
            df_Gene["exp"] = np.log10(df_Gene["exp"])
            gene_id = df_Exp.index[i]
            sub_ref = df_Gene[df_Gene["genotype"]==0]
            for genotype in [-4,-32,-31,-23,-22,-21,-12,-11,11,12,21,22,23,31,32,4]:
                sub_alt = df_Gene[df_Gene["genotype"]==genotype]
                if len(sub_alt.index)<10:
                    continue
                ttest_output=  ttest_ind(sub_alt["exp"],sub_ref["exp"], equal_var = False)
                fold_change = np.mean(sub_alt["exp"])-np.mean(sub_ref["exp"])
                number_of_case = len(sub_alt.index)
                exp_cnv.write("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\n" % (gene_id,genotype,number_of_case, np.mean(sub_alt["exp"]), np.mean(sub_ref["exp"]), ttest_output[0],ttest_output[1],fold_change))
                if ttest_output[1] <0.00001 and fold_change>1:
                    print genotype, gene_id
                    sns.boxplot(x="genotype",y="exp",data=df_Gene)
                    sns.stripplot(x="genotype", y="exp", data=df_Gene, jitter=True,color="grey")
                    plt.title(gene_id)
                    plt.show()
                # if ttest_output[1]<0.000001:


if __name__ == '__main__':
    main()