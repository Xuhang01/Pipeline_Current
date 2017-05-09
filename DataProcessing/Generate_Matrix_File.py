import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time
import os

Chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]


def Generate_Exp_matrix(Exp_file,Gene_file,Exp_matrix_file):
    if os.path.exists(Exp_matrix_file):
        print "---- ! Exp_matrix_file already exists"
        return True
    df_Exp = pd.read_table(Exp_file, sep="\t")
    specimen_list = df_Exp.icgc_specimen_id.unique().tolist()
    Gene_list = df_Exp[df_Exp.icgc_specimen_id==specimen_list[0]].gene_id.tolist()
    Exp_matrix = []
    for specimen_id in specimen_list:
        Exp_matrix.append(df_Exp[df_Exp.icgc_specimen_id==specimen_id].normalized_read_count.tolist())
    Exp_matrix = np.transpose(np.array(Exp_matrix))
    df_Exp_matrix = pd.DataFrame(data=Exp_matrix, columns=specimen_list, index=Gene_list)
    df_Gene = pd.read_table(Gene_file, sep="\t")
    Gene_list = df_Exp_matrix.index.values
    Gene_Position_Order = []
    Ensembl_Gene_List = df_Gene.gene_id.tolist()
    for i in range(len(Gene_list)):
        gene_id = Gene_list[i]
        Gene_Position_Order.append(-1 if gene_id not in Ensembl_Gene_List else Ensembl_Gene_List.index(gene_id)) 
    df_Exp_matrix["order"]=Gene_Position_Order
    df_Exp_matrix = df_Exp_matrix[df_Exp_matrix.order!=-1].sort(["order"])
    df_Exp_matrix = df_Exp_matrix.drop("order",axis=1)
    df_Exp_matrix.to_csv(Exp_matrix_file,sep="\t")
    return True
def Generate_Gene_CNV_matrix(CNV_file, Gene_file,CNV_matrix_file):
    if os.path.exists(CNV_matrix_file):
        print "---- ! CNV_matrix_file already exists"
        return True
    df_CNV = pd.read_table(CNV_file, sep="\t")
    df_Gene = pd.read_table(Gene_file,sep="\t")
    specimen_list = df_CNV.icgc_specimen_id.unique().tolist()
    CNV_matrix = []
    # CNV dict is used to simplify the program and save time
    CNV_dict={}  
    for chrom in Chrom_list:
        CNV_dict[chrom]=df_CNV[df_CNV.chromosome==chrom]
    df_Gene 
    for i in range(len(df_Gene)):
        gene_chrom, gene_start, gene_end, gene_id = df_Gene.iloc[i,:]
        sub_df_CNV = CNV_dict[gene_chrom]; 
        sub_df_CNV = sub_df_CNV[np.logical_and(sub_df_CNV.chromosome_start<gene_start,sub_df_CNV.chromosome_end>gene_end)]
        sub_specimen_list = sub_df_CNV.icgc_specimen_id.unique().tolist()
        copy_number_list = sub_df_CNV.copy_number.tolist()
        CNV_matrix.append([2 if specimen_id not in sub_specimen_list else copy_number_list[sub_specimen_list.index(specimen_id)] for specimen_id in specimen_list])
    df_CNV_matrix = pd.DataFrame(data=CNV_matrix,columns=specimen_list,index=df_Gene.gene_id.tolist())
    df_CNV_matrix.to_csv(CNV_matrix_file, sep="\t")
    return True
def Generate_CNV_Genotype_matrix(CNV_file,TAD_file,CNV_TAD_file,CNV_Genotype_file):
    if os.path.exists(CNV_Genotype_file):
        print "---- ! CNV_Genotype_file already exists"
        return True
    df_TAD = pd.read_table(TAD_file)
    TAD_list = ["-".join([str(each) for each in df_TAD.iloc[i,:]]) for i in range(len(df_TAD))]
    print "----There are %d TAD" % (len(TAD_list))
    # 1. Add TAD file to each SV
    if not os.path.exists(CNV_TAD_file):
        df_CNV = pd.read_table(CNV_file, sep="\t")
        print "----There are %d CNVs" % len(df_CNV)
        TAD_of_CNV_list = [[],[]]
        TAD_dict = {}
        for chrom in Chrom_list:
            TAD_dict[chrom]=df_TAD[df_TAD.chrom==chrom]
        
        df_CNV_position = df_CNV[["chromosome","chromosome_start","chromosome_end"]]
        len_of_CNV = len(df_CNV)
        for j in range(len_of_CNV):
            cnv_chrom, cnv_start, cnv_end = df_CNV_position.iloc[j,:]
            # variant_type, copy_number = df_CNV[["mutation_type","copy_number"]].iloc[j,:]
            sub_tad= TAD_dict[cnv_chrom]

            tad_start_index = sub_tad[np.logical_and(sub_tad.start<cnv_start, sub_tad.end >cnv_start)].index.values
            tad_end_index = sub_tad[np.logical_and(sub_tad.start<cnv_end,sub_tad.end>cnv_end)].index.values
            if len(tad_start_index)==1:
                TAD_of_CNV_list[0].append(TAD_list[tad_start_index[0]])
            else:
                TAD_of_CNV_list[0].append("noTAD")
            if len(tad_end_index) == 1:
                TAD_of_CNV_list[1].append(TAD_list[tad_end_index[0]])
            else:
                TAD_of_CNV_list[1].append("noTAD")
        df_CNV["tad_start"] = TAD_of_CNV_list[0]
        df_CNV["tad_end"] = TAD_of_CNV_list[1]
        df_CNV_TAD = df_CNV
        df_CNV_TAD.to_csv(CNV_TAD_file,sep="\t",index=False)
    else:
        df_CNV_TAD = pd.read_table(CNV_TAD_file,sep="\t")

    # 2. Generate Genotype matrix of CNV
    specimen_list = df_CNV_TAD.icgc_specimen_id.unique().tolist()
    # filtering step

    CNV_Genotype_matrix = [[0 for specimen_id in specimen_list] for tad_id in TAD_list]
    df_CNV_TAD_part = df_CNV_TAD[["tad_start","tad_end"]]
    icgc_specimen_id_list = df_CNV_TAD["icgc_specimen_id"].tolist()
    for i in range(len(df_CNV_TAD)):
        tad_start, tad_end = df_CNV_TAD_part.iloc[i,:]
        specimen_id = icgc_specimen_id_list[i]

        if tad_start == "noTAD" or tad_end == "noTAD" or tad_start == tad_end:
            continue
        else:
            tad_start_index = TAD_list.index(tad_start)
            tad_end_index = TAD_list.index(tad_end)
            specimen_index = specimen_list.index(specimen_id)
            CNV_Genotype_matrix[tad_start_index][specimen_index] +=1
            CNV_Genotype_matrix[tad_end_index][specimen_index] +=1
    df_CNV_Genotype = pd.DataFrame(data=CNV_Genotype_matrix, columns=specimen_list,index=TAD_list)
    df_CNV_Genotype.to_csv(CNV_Genotype_file,sep="\t",index_col=0)
    return True
def Generate_SV_Genotype_matrix(SV_file,TAD_file,SV_TAD_file,SV_Genotype_file):
    if os.path.exists(SV_Genotype_file):
        print "---- ! SV_Genotype_file already exists"
        return True    
    df_TAD = pd.read_table(TAD_file)
    TAD_list =["-".join([str(each) for each in df_TAD.iloc[i,:]]) for i in range(len(df_TAD))]
    print "----There are %d TADs" % len(df_TAD)  
    # 1. Add TAD file to each SV
    if not os.path.exists(SV_TAD_file):   
        df_SV = pd.read_table(SV_file,sep="\t")
        print "----There are %d SVs" % len(df_SV)
        TAD_of_SV_list = [[],[]]
        TAD_dict = {}
        for chrom in Chrom_list:
            TAD_dict[chrom]=df_TAD[df_TAD.chrom==chrom]
        for j in range(len(df_SV)):
            sv_chrom_from, sv_bkpt_from, sv_chrom_to, sv_bkpt_to = df_SV[["chr_from","chr_from_bkpt","chr_to","chr_to_bkpt"]].iloc[j,:]
            # variant_type = df_SV["variant_type"][j]
            sub_tad_from = TAD_dict[sv_chrom_from]
            tad_from_index = sub_tad_from[np.logical_and(sub_tad_from.start<sv_bkpt_from, sub_tad_from.end>sv_bkpt_from)].index.values
            sub_tad_to = TAD_dict[sv_chrom_to]
            tad_to_index = sub_tad_to[np.logical_and(sub_tad_to.start<sv_bkpt_to, sub_tad_to.end>sv_bkpt_to)].index.values
            if len(tad_from_index)==1:
                TAD_of_SV_list[0].append(TAD_list[tad_from_index[0]])
            else:
                TAD_of_SV_list[0].append("noTAD")
            if len(tad_to_index)==1:
                TAD_of_SV_list[1].append(TAD_list[tad_to_index[0]])
            else:
                TAD_of_SV_list[1].append("noTAD")
        df_SV["tad_from"]= TAD_of_SV_list[0]
        df_SV["tad_to"] = TAD_of_SV_list[1]
        df_SV_TAD = df_SV
        df_SV_TAD.to_csv(SV_TAD_file,sep="\t",index=False)
    
    else:
        df_SV_TAD = pd.read_table(SV_TAD_file,sep="\t")
    
    # 2. Generate Genotype matrix of SV
    specimen_list = df_SV_TAD.icgc_specimen_id.unique().tolist()
    # filtering step
    print "----number of SVs before filtering :", len(df_SV_TAD)
    df_SV_TAD = df_SV_TAD[df_SV_TAD.tad_from != df_SV_TAD.tad_to]
    print "----number of SVs after filtering :", len(df_SV_TAD)

    SV_Genotype_matrix = [[0 for specimen_id in specimen_list] for tad_id in TAD_list]
    for i in range(len(df_SV_TAD)):
        tad_from, tad_to = df_SV_TAD[["tad_from","tad_to"]].iloc[i,:]
        specimen_id = df_SV_TAD["icgc_specimen_id"].iat[i]

        if tad_from == "noTAD" or tad_to == "noTAD":
            continue
        else:
            tad_from_index = TAD_list.index(tad_from)
            tad_to_index = TAD_list.index(tad_to)
            specimen_index = specimen_list.index(specimen_id)
            SV_Genotype_matrix[tad_from_index][specimen_index] +=1
            SV_Genotype_matrix[tad_to_index][specimen_index] +=1
    df_SV_Genotype = pd.DataFrame(data=SV_Genotype_matrix, columns=specimen_list,index=TAD_list)
    df_SV_Genotype.to_csv(SV_Genotype_file,sep="\t",index_col=0)
    return True
def Prepare_matrix_file(SV_file, CNV_file, Exp_file, SV_TAD_file,CNV_TAD_file,SV_Genotype_file, CNV_Genotype_file, CNV_matrix_file, Exp_matrix_file):
    Gene_file = "../ExampleData/Processed_ENSEMBL gene.txt"
    TAD_file = "../ExampleData/TAD hESC Combined.txt"

    print "start Generate_Exp_matrix:"
    Generate_Exp_matrix(Exp_file,Gene_file,Exp_matrix_file)
    print "start Generate_Gene_CNV_matrix:"
    Generate_Gene_CNV_matrix(CNV_file, Gene_file,CNV_matrix_file)
    print "start Generate_CNV_Genotype_matrix:"
    Generate_CNV_Genotype_matrix(CNV_file,TAD_file,CNV_TAD_file,CNV_Genotype_file)
    print "start Generate_SV_Genotype_matrix:"
    Generate_SV_Genotype_matrix(SV_file,TAD_file,SV_TAD_file,SV_Genotype_file)
    return True
def main():
    from time import time
    Gene_file = "../ExampleData/Processed_ENSEMBL gene.txt"
    TAD_file = "../ExampleData/TAD hESC Combined.txt"
    CNV_file = "../ExampleData/copy_number_somatic_mutation.OV-AU.tsv"
    SV_file = "../ExampleData/structural_somatic_mutation.OV-AU.tsv"
    Exp_file = "../ExampleData/exp_seq.OV-AU.tsv"

    SV_TAD_file = "../Temp_file/Add_TAD_file/SV_TAD_file.OV-AU.txt"
    CNV_TAD_file = "../Temp_file/Add_TAD_file/CNV_TAD_file.OV-AU.txt"

    CNV_matrix_file = "../Temp_file/matrix_file/OV-AU.cnv_matrix.tsv"
    Exp_matrix_file = "../Temp_file/matrix_file/OV-AU.exp_matrix.tsv"
    SV_Genotype_file = "../Temp_file/matrix_file/OV-AU.SV.Genotype.tsv"
    CNV_Genotype_file = "../Temp_file/matrix_file/OV-AU.CNV.Genotype.tsv"
    ##### Generate SV genotype matrix
    # print "Start generating SV matrix:"
    # start = time()
    # Generate_SV_Genotype_matrix(SV_file,TAD_file,SV_TAD_file,SV_Genotype_file)
    # print "Generate_SV_Genotype_matrix costs %ds" % (time()-start)
    ##### Generate CNV Genotype matrix
    # print "Start generating CNV genotype matrix:"
    # start = time()
    # Generate_CNV_Genotype_matrix(CNV_file,TAD_file,CNV_TAD_file,CNV_Genotype_file)
    # print "Generate_CNV_Genotype_matrix costs %ds" % (time()-start)
    ##### Generate CNV matrix
    # print "Start Generating CNV matrix:"
    # start = time()
    # Generate_Gene_CNV_matrix(CNV_file,Gene_file,CNV_matrix_file)
    # print "Generate_Gene_CNV_matrix costs %ds" % (time()-start)
    ##### Generate Expression matrix
    # print "start Generating Expression matrix"
    # start = time()
    # Generate_Exp_matrix(Exp_file, Gene_file,Exp_matrix_file)
    # print "Generate_Exp_matrix costs %d s" % (time()-start)
    Prepare_matrix_file(SV_file, CNV_file, Exp_file, SV_TAD_file,CNV_TAD_file,SV_Genotype_file, CNV_Genotype_file, CNV_matrix_file, Exp_matrix_file)

if __name__ == "__main__":
    main()