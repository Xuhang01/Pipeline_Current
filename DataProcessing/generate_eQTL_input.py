import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
# import seaborn as sns

def Preprocess_CNV_Genotype(Exp_matrix_file,CNV_matrix_file, CNV_Genotype_file,eQTL_Exp_file,eQTL_Genotype_file):
    if os.path.exists(eQTL_Exp_file) and os.path.exists(eQTL_Genotype_file):
        print "----%s already exists" % eQTL_Exp_file
        print "----%s already exists" % eQTL_Genotype_file
        return True
    df_Exp = pd.read_table(Exp_matrix_file, sep="\t",index_col=0)
    df_CNV_Genotype = pd.read_table(CNV_Genotype_file, sep="\t",index_col=0)
    df_Gene_CNV = pd.read_table(CNV_matrix_file, sep="\t", index_col=0)
    
    common_specimen_id = [specimen_id for specimen_id in df_Exp if (specimen_id in df_CNV_Genotype and specimen_id in df_Gene_CNV)]
    df_CNV_Genotype = df_CNV_Genotype[common_specimen_id]
    df_Exp = df_Exp[common_specimen_id]
    df_Gene_CNV = df_Gene_CNV[common_specimen_id]

    df_Exp = df_Exp.T; df_Gene_CNV = df_Gene_CNV.T
    common_gene_id = [gene_id for gene_id in df_Gene_CNV if gene_id in df_Exp]
    df_Exp = df_Exp[common_gene_id]
    df_Gene_CNV = df_Gene_CNV[common_gene_id]
    df_Exp = df_Exp.T; df_Gene_CNV = df_Gene_CNV.T
    print "----Number of Specimen: Gene_Exp (%d), Gene_copynumber (%d), CNV_genotype (%d)" % (len(df_Exp.columns), len(df_Gene_CNV.columns), len(df_CNV_Genotype.columns))
    print "----Number of Genes: Gene_copynumber (%d), Gene_Exp (%d)" % (len(df_Exp), len(df_Gene_CNV))

    
    ########## Filtering step
    print "There are %d genes before filtering" % len(df_Exp)
    # (1) revise Expression data with copy number
    Exp_matrix = df_Exp.values
    Exp_matrix = Exp_matrix / df_Gene_CNV.values
    df_Exp = pd.DataFrame(data=Exp_matrix, columns=df_Exp.columns, index=df_Exp.index)
    
    df_Exp = df_Exp.T
    df_Exp = df_Exp[[gene_id for gene_id in df_Exp if [df_Exp[gene_id]>4][0].tolist().count(True) <3]]
    df_Exp = df_Exp.T
    print "There are %d genes after filtering out recurrent CNV" % len(df_Exp)
    # (2) filter out low exp gene. Exp > 0.1
    df_Exp = df_Exp[np.mean(df_Exp.values, axis=1)>0.1]
    print "There are %d genes after filtering out low expression genes" % len(df_Exp)
    # (3) filter out low expression variance gene. below 20th percentile
    Exp_matrix = df_Exp.values
    Exp_matrix = np.log2(Exp_matrix+0.1)
    std_array = np.std(Exp_matrix,axis=1)
    sort_std_array = np.sort(std_array)
    variance_threshold = sort_std_array[len(sort_std_array)/5]
    df_Exp = pd.DataFrame(data=Exp_matrix, columns=df_Exp.columns, index=df_Exp.index)
    df_Exp = df_Exp[np.std(df_Exp.values,axis=1)>variance_threshold]
    print "There are %d genes after filtering out low expression variance" % len(df_Exp)
    
    # (4) reGenotype number of bkpt as 0 for non or 1 for existance
    CNV_Genotype_matrix = df_CNV_Genotype.values
    CNV_Genotype_matrix = [[0 if df_CNV_Genotype[specimen_id][gene_id]==0 else 1 for specimen_id in df_CNV_Genotype.columns] for gene_id in df_CNV_Genotype.index]
    df_CNV_Genotype = pd.DataFrame(data=CNV_Genotype_matrix, columns=df_CNV_Genotype.columns, index=df_CNV_Genotype.index)
    ########## Output
    df_Exp_eQTL = df_Exp
    df_Genotype_eQTL = df_CNV_Genotype

    df_Exp_eQTL.to_csv(eQTL_Exp_file,sep="\t")
    df_Genotype_eQTL.to_csv(eQTL_Genotype_file,sep="\t")
    
    # Simple count of eQTL input
    print "There are %d specimen in Expression file" % len(df_Exp_eQTL.columns)
    print "There are %d genes in Expression file" % len(df_Exp_eQTL)
    print "There are %d specimen in Geontype file" % len(df_Genotype_eQTL.columns)
    print "There are %d genotypes in Genotype file" % len(df_Genotype_eQTL)
    return True
def Preprocess_SV_Genotype(Exp_matrix_file,CNV_matrix_file, SV_Genotype_file,SV_eQTL_Exp_file,SV_eQTL_Genotype_file):
    if os.path.exists(SV_eQTL_Exp_file) and os.path.exists(SV_eQTL_Genotype_file):
        print "----%s already exists" % SV_eQTL_Exp_file
        print "----%s already exists" % SV_eQTL_Genotype_file
        return True
    df_Exp = pd.read_table(Exp_matrix_file, sep="\t",index_col=0)
    df_SV_Genotype = pd.read_table(SV_Genotype_file, sep="\t",index_col=0)
    df_Gene_CNV = pd.read_table(CNV_matrix_file, sep="\t", index_col=0)

    common_specimen_id = [specimen_id for specimen_id in df_Exp if (specimen_id in df_SV_Genotype and specimen_id in df_Gene_CNV)]
    df_SV_Genotype = df_SV_Genotype[common_specimen_id]
    df_Exp = df_Exp[common_specimen_id]
    df_Gene_CNV = df_Gene_CNV[common_specimen_id]

    df_Exp = df_Exp.T; df_Gene_CNV = df_Gene_CNV.T
    common_gene_id = [gene_id for gene_id in df_Exp if gene_id in df_Gene_CNV]
    df_Exp = df_Exp[common_gene_id]
    df_Gene_CNV = df_Gene_CNV[common_gene_id]
    df_Exp = df_Exp.T; df_Gene_CNV = df_Gene_CNV.T
    print "----Number of Specimen: Gene_Exp (%d), Gene_copynumber (%d), CNV_genotype (%d)" % (len(df_Exp.columns), len(df_Gene_CNV.columns), len(df_SV_Genotype.columns))
    print "----Number of Genes: Gene_copynumber (%d), Gene_Exp (%d)" % (len(df_Exp), len(df_Gene_CNV))

    
    ########## Filtering step
    print "There are %d genes before filtering" % len(df_Exp)
    # (1) revise Expression data with copy number
    Exp_matrix = df_Exp.values
    Exp_matrix = Exp_matrix / df_Gene_CNV.values
    df_Exp = pd.DataFrame(data=Exp_matrix, columns=df_Exp.columns, index=df_Exp.index)
    
    df_Exp = df_Exp.T
    df_Exp = df_Exp[[gene_id for gene_id in df_Exp if [df_Exp[gene_id]>4][0].tolist().count(True) <3]]
    df_Exp = df_Exp.T
    print "There are %d genes after filtering out recurrent CNV" % len(df_Exp)
    # (2) filter out low exp gene. Exp > 0.1
    df_Exp = df_Exp[np.mean(df_Exp.values, axis=1)>0.1]
    print "There are %d genes after filtering out low expression genes" % len(df_Exp)
    # (3) filter out low expression variance gene. below 20th percentile
    Exp_matrix = df_Exp.values
    Exp_matrix = np.log2(Exp_matrix+0.1)
    std_array = np.std(Exp_matrix,axis=1)
    sort_std_array = np.sort(std_array)
    variance_threshold = sort_std_array[len(sort_std_array)/5]
    df_Exp = pd.DataFrame(data=Exp_matrix, columns=df_Exp.columns, index=df_Exp.index)
    df_Exp = df_Exp[np.std(df_Exp.values,axis=1)>variance_threshold]
    print "There are %d genes after filtering out low expression variance" % len(df_Exp)
    
    # (4) reGenotype number of bkpt as 0 for non or 1 for existance
    SV_Genotype_matrix = df_SV_Genotype.values
    SV_Genotype_matrix = [[0 if df_SV_Genotype[specimen_id][gene_id]==0 else 1 for specimen_id in df_SV_Genotype.columns] for gene_id in df_SV_Genotype.index]
    df_SV_Genotype = pd.DataFrame(data=SV_Genotype_matrix, columns=df_SV_Genotype.columns, index=df_SV_Genotype.index)
    # (5) remove the snp (tad) with low allele frequency. >5/110
    df_SV_Genotype = df_SV_Genotype[np.sum(df_SV_Genotype.values,axis=1)>5]
    ########## Output
    df_Exp_eQTL = df_Exp
    df_Genotype_eQTL = df_SV_Genotype

    df_Exp_eQTL.to_csv(SV_eQTL_Exp_file,sep="\t")
    df_Genotype_eQTL.to_csv(SV_eQTL_Genotype_file,sep="\t")
    
    # Simple count of eQTL input
    print "There are %d specimen in Expression file" % len(df_Exp_eQTL.columns)
    print "There are %d genes in Expression file" % len(df_Exp_eQTL)
    print "There are %d specimen in Geontype file" % len(df_Genotype_eQTL.columns)
    print "There are %d genotypes in Genotype file" % len(df_Genotype_eQTL)
    return True
def main():
    Exp_matrix_file = "../Temp_file/matrix_file/OV-AU.exp_matrix.tsv"
    CNV_Genotype_file = "../Temp_file/matrix_file/OV-AU.CNV.Genotype.tsv"
    CNV_matrix_file = "../Temp_file/matrix_file/OV-AU.cnv_matrix.tsv"

    SV_Genotype_file = "../Temp_file/matrix_file/OV-AU.SV.Genotype.tsv"
    
    eQTL_Exp_file = "../Temp_file/eQTL/OV-AU.Exp_file.txt"
    eQTL_Genotype_file = "../Temp_file/eQTL/OV-AU.Genotype_file.txt"
    # Preprocess_CNV_Genotype(Exp_matrix_file,CNV_matrix_file, CNV_Genotype_file,eQTL_Exp_file,eQTL_Genotype_file)
    
    SV_eQTL_Genotype_file = "../Temp_file/SV_eQTL/OV-AU.Genotype_file.txt"
    SV_eQTL_Exp_file = "../Temp_file/SV_eQTL/OV-AU.Exp_file.txt"
    Preprocess_SV_Genotype(Exp_matrix_file,CNV_matrix_file, SV_Genotype_file,SV_eQTL_Exp_file,SV_eQTL_Genotype_file)
if __name__ == '__main__':
    main()