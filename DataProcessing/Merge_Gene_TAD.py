import pandas as pd
import numpy as np
import os
def Merge_Gene_TAD(Gene_file, TAD_file, Gene_TAD_file):
    print "start Merge_Gene_TAD:"
    if os.path.exists(Gene_TAD_file):
        print "---- ! Gene_TAD_file already exists"
        return True
    df_Gene = pd.read_table(Gene_file,sep="\t")
    df_TAD = pd.read_table(TAD_file,sep="\t")

    TAD_info_list=[]
    for i in range(len(df_Gene)):
        chrom, start, end, gene_id = df_Gene.iloc[i,:]
        sub_df_TAD = df_TAD[np.logical_and(df_TAD["chrom"]==chrom,np.logical_and(df_TAD["start"]<start, df_TAD["end"]>end))]
        if len(sub_df_TAD)!=1:
            TAD_info_list.append("noTAD")
            #print chrom, start, end, gene_id
        else:
            tad_id = "-".join([str(each) for each in sub_df_TAD.values[0]])
            TAD_info_list.append(tad_id)
    df_Gene["TAD"]=TAD_info_list
    df_Gene.to_csv(Gene_TAD_file,sep="\t",index=False)
    return True

def main():
    
    from time import time
    start = time()
    Gene_file = "../ExampleData/Processed_ENSEMBL gene.txt"
    TAD_file = "../ExampleData/TAD hESC Combined.txt"
    Gene_TAD_file = "../Temp_file/ENSEMBL_gene_with_TADinfo.txt"
    Merge_Gene_TAD(Gene_file, TAD_file, Gene_TAD_file)
    print "Merge_Gene_TAD costs %ds" % (time()-start)

if __name__ == "__main__":
    main()