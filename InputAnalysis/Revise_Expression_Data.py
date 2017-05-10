import os
import pandas as pd
def Revise_Expression_file(Exp_file, revised_Exp_file):

    df_Exp = pd.read_table(Exp_file, sep="\t")
    print "Total Count of expression data: %d" % len(df_Exp)
    Donor_list = df_Exp.icgc_donor_id.unique().tolist()
    print "Count of donor_id: %d" % len(Donor_list)
    Specimen_list = df_Exp.icgc_specimen_id.unique().tolist()
    print "Count of specimen_id: %d " % len(Specimen_list)
    Gene_list =df_Exp.gene_id.unique().tolist()
    print "Count of Gene_list: %d" % len(Gene_list)

    Common_Gene_list =[]
    for gene_id in Gene_list:
        if gene_id != df_Exp.gene_id.iat[0]:
            print "Wrong idea"
        else:
            print "Right idea"
            print len(df_Exp)
            index_0 = df_Exp.index[0]
            index_1 = df_Exp.index[0]
            while df_Exp.iat[index_1,0] == df_Exp.iat[index_0,0]:
                index_1+=1
            df_Exp = df_Exp.drop(range(index_0,index_1))




def main():
    Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/LIRI-JP/sort_exp_seq.LIRI-JP.tsv"
    revised_Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/LIRI-JP/exp_seq.LIRI-JP.tsv"

    Revise_Expression_file(Exp_file, revised_Exp_file)

if __name__ == '__main__':
    main()