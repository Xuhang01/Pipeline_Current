import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time
import os

Chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

def Generate_BKPT_Genotype_matrix(df_BKPT_TAD, TAD_list, donor_list):
    df_Genotype = pd.DataFrame(data=[[0 for i in donor_list] for j in TAD_list], columns=donor_list, index=TAD_list)
    len_of_BKPT = len(df_BKPT_TAD)
    for i in range(len_of_BKPT):
        donor_id, tad_id  = df_BKPT_TAD.iloc[i,:]
        if tad_id == "noTAD":
            continue
        else:
            df_Genotype[donor_id][tad_id]+=1
    return df_Genotype

def Generate_CNV_Genotype_matrix(TAD_list, CNV_TAD_file,CNV_Genotype_file):
    if os.path.exists(CNV_Genotype_file):
        print "---- ! CNV_Genotype_file already exists"
        return True
    df_CNV_TAD = pd.read_table(CNV_TAD_file,sep="\t")
    donor_list = df_CNV_TAD.icgc_donor_id.unique().tolist()
    df_CNV_start = df_CNV_TAD[["icgc_donor_id","tad_start"]]
    df_CNV_end = df_CNV_TAD[["icgc_donor_id","tad_end"]]
    df_CNV_start_Genotype = Generate_BKPT_Genotype_matrix(df_CNV_start, TAD_list, donor_list)
    df_CNV_end_Genotype = Generate_BKPT_Genotype_matrix(df_CNV_end, TAD_list, donor_list)
    df_CNV_Genotype = pd.DataFrame(data=df_CNV_start_Genotype.values+df_CNV_end_Genotype.values,columns=donor_list, index=TAD_list)
    df_CNV_Genotype.to_csv(CNV_Genotype_file,sep="\t",index_col=0)
    return True
def Generate_SV_Genotype_matrix(TAD_list, SV_TAD_file,SV_Genotype_file):
    if os.path.exists(SV_Genotype_file):
        print "---- ! SV_Genotype_file already exists"
        return True    
    df_SV_TAD = pd.read_table(SV_TAD_file,sep="\t")
    donor_list = df_SV_TAD.icgc_donor_id.unique().tolist()
    df_SV_from = df_SV_TAD[["icgc_donor_id", "tad_from"]]
    df_SV_to = df_SV_TAD[["icgc_donor_id","tad_to"]]
    df_SV_from_Genotype = Generate_BKPT_Genotype_matrix(df_SV_from, TAD_list, donor_list)
    df_SV_to_Genotype = Generate_BKPT_Genotype_matrix(df_SV_to, TAD_list, donor_list)
    df_SV_Genotype = pd.DataFrame(data=df_SV_from_Genotype.values+df_SV_to_Genotype.values, columns=donor_list, index=TAD_list)
    df_SV_Genotype.to_csv(SV_Genotype_file,sep="\t",index_col=0)
    return True

def main():
    from time import time
    TAD_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline_Current\Reference_files\TAD hESC Combined.txt"
    TAD_list = pd.read_table(TAD_file, sep="\t").tad_id.tolist()

    CNV_TAD_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\Temp_file\Add_TAD_file\CNV_TAD_file.OV-AU.txt"
    SV_TAD_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\Temp_file\Add_TAD_file\SV_TAD_file.OV-AU.txt"

    SV_Genotype_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\Genotype\OV-AU.SV.Genotype.tsv"
    CNV_Genotype_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\Genotype\OV-AU.CNV.Genotype.tsv"
    
    # Generate_CNV_Genotype_matrix(TAD_list, CNV_TAD_file, CNV_Genotype_file)
    Generate_SV_Genotype_matrix(TAD_list, SV_TAD_file, SV_Genotype_file)
if __name__ == "__main__":
    main()