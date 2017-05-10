import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time
import os
from scipy.sparse import coo_matrix 

Chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]


def Generate_Exp_matrix(Add_Symbol_Exp_file, df_Gene, Exp_matrix_file):
    # Add_Symbol_Exp_file is the directly download Exp_seq file from ICGC database
    # Exp_matrix_file is the output of this function
    # df_Gene is a dataframe of gene, columns contain ["chrom", "start","end","gene_symbol"]
    if os.path.exists(Exp_matrix_file):
        print "----!Exp_matrix_file already exists"
        return True
    with open(Add_Symbol_Exp_file) as exp:
        Symbol_list = []
        current_donor_id = ""
        for eachline in exp:
            if eachline.startswith("icgc"):
                continue
            else:
                eachline_info = eachline.strip().split("\t")
                gene_symbol = eachline_info[-1]
                donor_id = eachline_info[0]
                if current_donor_id == "":
                    Symbol_list.append(gene_symbol)
                    print gene_symbol
                    current_donor_id = donor_id
                elif current_donor_id != donor_id:
                    print gene_symbol
                    break
                else:
                    Symbol_list.append(gene_symbol)

    count_gene = len(Symbol_list)
    print "----There are %d genes studied in this project" % count_gene
    with open(Add_Symbol_Exp_file) as exp:
        donor_list = []
        read_count_list = []
        line_count = 0
        for eachline in exp:
            if eachline.startswith("icgc"):
                continue
            else:
                line_count+=1
                eachline_info = eachline.strip().split("\t")
                donor_id = eachline_info[0]
                normalized_read_count = float(eachline_info[8])

                if donor_list == []:
                    donor_list.append(donor_id)
                elif line_count > count_gene:
                    if donor_id == donor_list[-1]:
                        continue
                    else:
                        donor_list.append(donor_id)
                        line_count=1
                read_count_list.append(normalized_read_count)


    row = np.array([[i for i in range(len(Symbol_list))] for j in range(len(donor_list))]).reshape(1,-1)[0]
    col = np.array([[j for i in range(len(Symbol_list))] for j in range(len(donor_list))]).reshape(1,-1)[0]
    data = read_count_list
    
    print len(Symbol_list), len(donor_list)
    print len(row), len(col), len(data)

    Exp_matrix = coo_matrix((data,(row, col))).toarray()
    df_Exp_matrix = pd.DataFrame(data=Exp_matrix, columns=donor_list, index=Symbol_list)
    # Gene_Position_Order = []
    # Reference_Symbol_list = df_Gene["gene_symbol"].tolist()
    # for i in range(len(Symbol_list)):
    #     gene_symbol = Symbol_list[i]
    #     Gene_Position_Order.append(-1 if gene_symbol not in Reference_Symbol_list else Reference_Symbol_list.index(gene_symbol)) 
    # df_Exp_matrix["order"]=Gene_Position_Order
    # df_Exp_matrix = df_Exp_matrix[df_Exp_matrix.order!=-1].sort(["order"])
    # df_Exp_matrix = df_Exp_matrix.drop("order",axis=1)
    df_Exp_matrix.to_csv(Exp_matrix_file,sep="\t")
    return True

def Genreate_Gene_CNA_matrix(CNV_file, df_Gene,CNV_matrix_file):
    # CNV_file is the directly download Exp_seq file from ICGC database
    # CNV_matrix_file is the output of this function
    # df_Gene is a dataframe of gene, columns contain ["chrom", "start","end","gene_symbol"]
    if os.path.exists(CNV_matrix_file):
        print "---- ! CNV_matrix_file already exists"
        return True
    df_CNV = pd.read_table(CNV_file, sep="\t", dtype={"chromosome":"str"})
    donor_list = df_CNV.icgc_donor_id.unique().tolist()
    CNV_matrix = []
    # CNV dict is used to simplify the program and save time
    CNV_dict={}  
    for chrom in Chrom_list:
        CNV_dict[chrom]=df_CNV[df_CNV.chromosome==chrom]
    for i in range(len(df_Gene)):
        if i %1000 ==0:
            print i
        gene_symbol, gene_chrom, gene_start, gene_end = df_Gene.iloc[i,:]
        sub_df_CNV = CNV_dict[gene_chrom]; 
        sub_df_CNV = sub_df_CNV[np.logical_and(sub_df_CNV.chromosome_start<gene_start,sub_df_CNV.chromosome_end>gene_end)]
        sub_donor_list = sub_df_CNV.icgc_donor_id.unique().tolist()
        copy_number_list = sub_df_CNV.copy_number.tolist()
        CNV_matrix.append([2 if donor_id not in sub_donor_list else copy_number_list[sub_donor_list.index(donor_id)] for donor_id in donor_list])
    df_CNV_matrix = pd.DataFrame(data=CNV_matrix,columns=donor_list,index=df_Gene["gene_symbol"].tolist())
    df_CNV_matrix.to_csv(CNV_matrix_file, sep="\t")
    return True


def main():
    from time import time
    Gene_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Pipeline_Current/Reference_files/Processed_HGNC.txt"
    TAD_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Pipeline_Current/Reference_files/TAD hESC Combined.txt"
    CNV_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/CNV/copy_number_somatic_mutation.PACA-CA.tsv"
    Add_Symbol_Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/PACA-CA/exp_seq.Symbol.PACA-CA.tsv"

    Gene_CNA_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/PACA-CA/PACA-CA.cnv_matrix.tsv"
    Exp_matrix_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/PACA-CA/PACA-CA.exp_matrix.tsv"
    
    df_Gene = pd.read_table(Gene_file)[["Approved Symbol","chrom","start","end"]].rename(index=str, columns={"Approved Symbol": "gene_symbol"})
    # Generate_Exp_matrix(Add_Symbol_Exp_file, df_Gene, Exp_matrix_file)
    Genreate_Gene_CNA_matrix(CNV_file,df_Gene, Gene_CNA_file)


if __name__ == "__main__":
    main()