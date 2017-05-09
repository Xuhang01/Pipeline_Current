import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import ranksums
import matplotlib.pyplot as plt

def boxplot_show_difference(Exp_file, Genotype_file, cis_eqtl_file):
    df_exp = pd.read_table(Exp_file,index_col=0)
    df_genotype = pd.read_table(Genotype_file,index_col=0)
    df_cis_eQTL = pd.read_table(cis_eqtl_file,sep="\t")
    Specimen_list = df_genotype.columns.values

    Number_of_sign = 5

    for i in range(Number_of_sign):
        gene_id, snp_id = df_cis_eQTL[["gene","SNP"]].iloc[i,:]

        genotype_list = df_genotype[df_genotype.index.values==snp_id].values[0].tolist()
        exp_list = df_exp[df_exp.index.values==gene_id].values[0].tolist()
        df_sub = pd.DataFrame(data=np.array([exp_list,genotype_list]).transpose(), columns=["exp","genotype"],index=Specimen_list)
        
        print snp_id, gene_id
        print np.mean(df_sub[df_sub.genotype==1].exp), np.mean(df_sub[df_sub.genotype==0].exp), np.mean(df_sub[df_sub.genotype==1].exp)-np.mean(df_sub[df_sub.genotype==0].exp)
        print (np.mean(df_sub[df_sub.genotype==1].exp)-np.mean(df_sub[df_sub.genotype==0].exp))/np.std(df_sub[df_sub.genotype==0].exp)
        print ttest_ind(df_sub[df_sub.genotype==1].exp,df_sub[df_sub.genotype==0].exp)
        print ranksums(df_sub[df_sub.genotype==1].exp,df_sub[df_sub.genotype==0].exp)
        plt.boxplot([df_sub[df_sub.genotype==1].exp,df_sub[df_sub.genotype==0].exp], labels=["alt","ref"])
        plt.show()

    return True
    
def Find_special_ciseQTL(cis_eqtl_file, Exp_file, Genotype_file):
    df_exp = pd.read_table(Exp_file,index_col=0)
    df_genotype = pd.read_table(Genotype_file,index_col=0)
    df_cis_eQTL = pd.read_table(cis_eqtl_file,sep="\t")
    Specimen_list = df_genotype.columns.values

    for i in range(1000):
        snp_id, gene_id = df_cis_eQTL[["SNP","gene"]].iloc[i,:]
        sub_exp = df_exp[df_exp.index.values==gene_id].values[0].tolist()
        sub_genotype = df_genotype[df_genotype.index.values==snp_id].values[0].tolist()
        df_sub = pd.DataFrame(data=np.array([sub_exp,sub_genotype]).transpose(), columns=["exp","genotype"],index=Specimen_list)
        

        # condition
        if np.mean(df_sub.exp)>0 and np.mean(df_sub[df_sub.genotype==1].exp)/np.mean(df_sub[df_sub.genotype==0].exp)>2 and np.sum(sub_genotype)>10:
            print snp_id, gene_id
            print np.mean(df_sub[df_sub.genotype==1].exp), np.mean(df_sub[df_sub.genotype==0].exp), np.mean(df_sub[df_sub.genotype==1].exp)/np.mean(df_sub[df_sub.genotype==0].exp)
            print ttest_ind(df_sub[df_sub.genotype==1].exp,df_sub[df_sub.genotype==0].exp)
            print ranksums(df_sub[df_sub.genotype==1].exp,df_sub[df_sub.genotype==0].exp)

            plt.boxplot([df_sub[df_sub.genotype==1].exp,df_sub[df_sub.genotype==0].exp], labels=["alt","ref"])
            plt.show()

def Generate_SV_bed(tad_id,gene_id, SV_file):
    df_SV = pd.read_table(SV_file,sep="\t")
    tad_chrom, tad_start, tad_end = tad_id.split("-"); tad_start= int(tad_start);tad_end=int(tad_end)

    sub_df_SV = df_SV[np.logical_or(df_SV.chr_from==tad_chrom, df_SV.chr_to==tad_chrom)]
    index_from = np.logical_and(sub_df_SV.chr_from_bkpt>tad_start, sub_df_SV.chr_from_bkpt<tad_end)
    index_to = np.logical_and(sub_df_SV.chr_to_bkpt>tad_start, sub_df_SV.chr_to_bkpt<tad_end)
    sub_df_SV = sub_df_SV[np.logical_or(index_from, index_to)]

    print sub_df_SV[["variant_type","chr_from","chr_from_bkpt","chr_to","chr_to_bkpt"]]

def main():
    Exp_file = "../Temp_file/eQTL/reGenotype/OV-AU.Exp_file.txt"
    Genotype_file = "../Temp_file/eQTL/reGenotype/OV-AU.Genotype_file.txt"
    cis_eqtl_file = "../Temp_file/eQTL/reGenotype/cis_eqtl.txt"
    gene_id = "ENSG00000127954"
    snp_id = "7-86720000-87640000"

    boxplot_show_difference(Exp_file, Genotype_file,cis_eqtl_file)

    #SV_file = r"C:\Users\JW Lab user\Desktop\3-22-2017 correlation of gene expression\Data\RawData\structural_somatic_mutation.OV-AU.tsv"
    
    #Generate_SV_bed(snp_id,gene_id, SV_file)

    #Find_special_ciseQTL(cis_eqtl_file, Exp_file, Genotype_file)
if __name__ =="__main__":
    main()