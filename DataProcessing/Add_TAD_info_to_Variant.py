import pandas as pd
import os
import numpy as np

Chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

def Detect_TAD_of_BKPT_list(df_BKPT, df_TAD_dict):
    TAD_of_BKPT_list = []
    len_of_BKPT_list = len(df_BKPT)
    for i in range(len_of_BKPT_list):
        if i %10000==0:
            print i
        bkpt_chrom, bkpt_position = df_BKPT.iloc[i,:]
        if bkpt_chrom not in Chrom_list:
            TAD_of_BKPT_list.append("noTAD")
        else:
            sub_df_TAD = df_TAD_dict[bkpt_chrom]
            tad_bkpt_index = sub_df_TAD[np.logical_and(sub_df_TAD.start<bkpt_position , sub_df_TAD.end>bkpt_position)].index.tolist()
            if len(tad_bkpt_index)==1:
                TAD_of_BKPT_list.append(sub_df_TAD["tad_id"][tad_bkpt_index[0]])
            elif len(tad_bkpt_index)==0:
                TAD_of_BKPT_list.append("noTAD")
            else:
                print "Error! There must be overlap between TADs"
    return TAD_of_BKPT_list

def Add_TAD_info_to_CNV(CNV_file,df_TAD,CNV_TAD_file):
    if os.path.exists(CNV_TAD_file):
        print "----!CNV_TAD_file already exists"
        return True
    df_CNV = pd.read_table(CNV_file, sep="\t", dtype={"chromosome":"str"})
    print "----There are %d CNVs" % len(df_CNV)
    df_TAD_dict = {}
    for chrom in Chrom_list:
        df_TAD_dict[chrom]=df_TAD[df_TAD.chrom==chrom]
    df_CNV["tad_start"] = Detect_TAD_of_BKPT_list(df_CNV[["chromosome","chromosome_start"]], df_TAD_dict)
    df_CNV["tad_end"] = Detect_TAD_of_BKPT_list(df_CNV[["chromosome","chromosome_end"]], df_TAD_dict)
    df_CNV.to_csv(CNV_TAD_file,sep="\t",index=False)

def Add_TAD_info_to_SV(SV_file, df_TAD, SV_TAD_file):
    if os.path.exists(SV_TAD_file):
        print "----!SV_TAD_file already exists"
        return True
    df_SV = pd.read_table(SV_file,sep="\t",dtype={"chr_from":"str","chr_to":"str"})
    print "----There are %d SVs" % len(df_SV)
    df_TAD_dict = {}
    for chrom in Chrom_list:
        df_TAD_dict[chrom]=df_TAD[df_TAD.chrom==chrom]
    df_SV["tad_from"]= Detect_TAD_of_BKPT_list(df_SV[["chr_from","chr_from_bkpt"]], df_TAD_dict)
    df_SV["tad_to"] = Detect_TAD_of_BKPT_list(df_SV[["chr_to","chr_to_bkpt"]], df_TAD_dict)
    df_SV.to_csv(SV_TAD_file,sep="\t",index=False)

def main():
    Gene_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline_Current\Reference_files\Processed_HNGC.txt"
    TAD_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline_Current\Reference_files\TAD hESC Combined.txt"
    df_TAD= pd.read_table(TAD_file,sep="\t")
    CNV_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\ExampleData\copy_number_somatic_mutation.OV-AU.tsv"
    CNV_TAD_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\Temp_file\Add_TAD_file\CNV_TAD_file.OV-AU.txt"
    SV_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\ExampleData\structural_somatic_mutation.OV-AU.tsv"
    SV_TAD_file = r"C:\Users\JW Lab user\Desktop\SVs\Pipeline\Temp_file\Add_TAD_file\SV_TAD_file.OV-AU.txt"

    Add_TAD_info_to_CNV(CNV_file, df_TAD, CNV_TAD_file)
    Add_TAD_info_to_SV(SV_file, df_TAD, SV_TAD_file)

if __name__ == '__main__':
    main()