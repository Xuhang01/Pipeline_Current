from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time

def Merge_GENE_TAD(df_GENE, df_TAD):
    TAD_info_list=[]
    for i in range(len(df_GENE)):
        chrom, start, end, gene_id = df_GENE.iloc[i,:]
        sub_df_TAD = df_TAD[np.logical_and(df_TAD["chrom"]==chrom,np.logical_and(df_TAD["start"]<start, df_TAD["end"]>end))]
        if len(sub_df_TAD)!=1:
            TAD_info_list.append("noTAD")
            #print chrom, start, end, gene_id
        else:
            tad_id = "-".join([str(each) for each in sub_df_TAD.values[0]])
            TAD_info_list.append(tad_id)
    df_GENE["TAD"]=TAD_info_list
    return df_GENE

def main():
    CNV_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/CNV/copy_number_somatic_mutation.OV-AU.tsv"
    SV_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/SV/structural_somatic_mutation.OV-AU.tsv"
    EXP_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/exp_seq.OV-AU.tsv"
    GENE_file = "Processed_ENSEMBL gene.txt"
    TAD_file = "TAD hESC Combined.txt"
    """
    ######
    df_GENE = pd.read_table(GENE_file,sep="\t")
    df_TAD = pd.read_table(TAD_file,sep="\t")
    df_GENE_TAD = Merge_GENE_TAD(df_GENE, df_TAD)
    df_GENE_TAD.to_csv("ENSEMBL_gene_with_TADinfo.txt",sep="\t",index=False)
    """
    GENE_TAD_file = "ENSEMBL_gene_with_TADinfo.txt"
    df_GT = pd.read_table(GENE_TAD_file)
    #######
    
    df_specimen = pd.read_table("specimen.OV-AU.tsv",sep="\t")
    df_SV = pd.read_table(SV_file,sep="\t")
    df_Exp = pd.read_table(EXP_file, sep="\t")
    SV_specimen = set(df_SV.icgc_specimen_id)
    Exp_specimen = set(df_Exp.icgc_specimen_id)

    print len(SV_specimen), len(Exp_specimen)

    Overlap_specimen = list(SV_specimen.intersection(Exp_specimen)) #110 specimen belongs to 93 donors
    Donor_list = list(df_specimen.icgc_donor_id.unique())
    ###########
    
    sub_df_specimen = df_specimen[[True if df_specimen.icgc_specimen_id[i] in Overlap_specimen else False for i in range(len(df_specimen))]]

    sub_df_specimen.to_csv("subset.specimen.OV-AU.tsv",sep="\t",index=False)
    
    # Pairs_of_Sample = []
    # for i in range(len(df_CoSV)):
    #     for j in range(i+1, len(df_CoSV)):
    #         donor_id_i, specimen_id_i, Count_of_SV_i = df_CoSV.iloc[i,]
    #         donor_id_j, specimen_id_j, Count_of_SV_j = df_CoSV.iloc[j,]

    #         if donor_id_i == donor_id_j:
    #             Pairs_of_Sample.append([donor_id_i, specimen_id_i, specimen_id_j,Count_of_SV_i,Count_of_SV_j])
    # df_Pair_Sample=pd.DataFrame(data=Pairs_of_Sample,columns=["icgc_donor_id","icgc_specimen_id1","icgc_specimen_id2","Count_of_SV1","Count_of_SV2"])
    # #print df_Pair_Sample


    # plt.plot(df_Pair_Sample["Count_of_SV1"],df_Pair_Sample["Count_of_SV2"],"o")

    # x = np.array(range(1000))
    # plt.plot(x,x,"r")
    # plt.plot(x,x*1.3,"g")
    # plt.plot(x,x*0.7,"g")
    # plt.show()
    


    # specimen_type_list =[]

    # for i in range(len(df_CoSV)):
    #     icgc_specimen_id = df_CoSV.iloc[i,1]
    #     specimen_type_list.append(df_specimen[df_specimen.icgc_specimen_id==icgc_specimen_id].specimen_type.tolist()[0])
    # df_CoSV["specimen_type"] = specimen_type_list
    
    # Pairs_of_primary_recurrent=[]
    # for donor_id in df_CoSV.icgc_donor_id.unique():
    #     sub_df_CoSV = df_CoSV[df_CoSV.icgc_donor_id==donor_id]
    #     if len(sub_df_CoSV) !=1:
    #         if "Primary tumour - solid tissue" in sub_df_CoSV.specimen_type.tolist() and "Recurrent tumour - other" in sub_df_CoSV.specimen_type.tolist():
    #             primary_df_CoSV = sub_df_CoSV[sub_df_CoSV.specimen_type=="Primary tumour - solid tissue"]
    #             recurrent_df_CoSV = sub_df_CoSV[sub_df_CoSV.specimen_type=="Recurrent tumour - other"]

    #             Pairs_of_primary_recurrent.append([donor_id, primary_df_CoSV.icgc_specimen_id.tolist()[0], recurrent_df_CoSV.icgc_specimen_id.tolist()[0], primary_df_CoSV.Count_of_SV.tolist()[0], recurrent_df_CoSV.Count_of_SV.tolist()[0]])
    # print Pairs_of_primary_recurrent
    # df_Primary_Recurrent = pd.DataFrame(data=Pairs_of_primary_recurrent, columns=["donor_id","spcimen_primary","specimen_recurrent","SV_primary","SV_recurrent"])
    # df_Primary_Recurrent.to_csv("Primary_Recurrent_comparison_SVcount.txt",index=False,sep="\t")
    # #df_CoSV.to_csv("Cout_of_SV.txt",index=False,sep="\t")
    
    # plt.plot(df_Primary_Recurrent.SV_primary, df_Primary_Recurrent.SV_recurrent, "o")
    # x = np.array(range(1000))
    # plt.plot(x,x,"r")
    # plt.show()



if __name__ == "__main__":
    main()