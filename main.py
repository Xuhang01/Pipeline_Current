import os
import re
import numpy as np

from DataProcessing.Remove_Duplicate_Expression import *
from DataProcessing.Generate_Matrix_File import *
from DataProcessing.Add_TAD_info_to_Variant import *
from DataProcessing.Add_Gene_Symbol import *
from DataProcessing.Generate_Genotype_file import *
from DataProcessing.generate_eQTL_input import *
from OutputAnalysis.Analyze_eQTL_output import boxplot_show_difference


def main():
    # Project_ID = "OV-AU"
    # Exp_Gene_ID_Type = "Ensembl Gene ID"
    # Project_ID = "BRCA-US"
    # Exp_Gene_ID_Type = "Approved Symbol"
    # Project_ID = "PACA-AU"
    # Exp_Gene_ID_Type = "Ensembl Gene ID"
    # Project_ID = "PAEN-AU"
    # Exp_Gene_ID_Type = "Ensembl Gene ID"
    # Project_ID = "PACA-CA"
    # Exp_Gene_ID_Type = "Ensembl Gene ID"
    Project_ID = "LIRI-JP"
    Exp_Gene_ID_Type = "Approved Symbol"

    # filenames uesed in this Pipeline
    Reference_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Pipeline_Current/Reference_files"
    Gene_file = os.path.join(Reference_dir, "Processed_HGNC.txt")
    TAD_file = os.path.join(Reference_dir, "TAD hESC Combined.txt")
    CNV_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/CNV"
    SV_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/SV"
    Exp_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression"
    Temp_file_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/"+Project_ID+"/"
    Code_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Pipeline_Current"

    # #### Remove Duplicate genes and samples
    # print "Remove Duplicate genes and samples:"
    # Exp_file = os.path.join(Exp_dir, re.sub("Project_ID", Project_ID, "exp_seq.Project_ID.tsv"))
    # remove_duplicate_exp_file = os.path.join(Temp_file_dir, re.sub("Project_ID",Project_ID,"remove_duplicate_exp_seq.Project_ID.tsv"))
    # Remove_Duplicate(Exp_file, remove_duplicate_exp_file)
    #### Add Gene_Symbol
    print "Add Gene Symbol:"
    Add_Symbol_Exp_file = os.path.join(Temp_file_dir, re.sub("Project_ID", Project_ID,"exp_seq.Symbol.Project_ID.tsv"))
    df_Gene = pd.read_table(Gene_file,sep="\t")
    Add_Gene_Symbol(remove_duplicate_exp_file, Add_Symbol_Exp_file, df_Gene, Exp_Gene_ID_Type)

    #### Generate Matrix files
    print "Generate Matrix files:"
    CNV_file = os.path.join(CNV_dir, re.sub("Project_ID", Project_ID, "copy_number_somatic_mutation.Project_ID.tsv"))
    Gene_CNA_file = os.path.join(Temp_file_dir, "Gene_CNA.tsv")
    Exp_matrix_file = os.path.join(Temp_file_dir, "Exp_matrix.tsv")

    df_Gene = pd.read_table(Gene_file)[["Approved Symbol","chrom","start","end"]].rename(index=str, columns={"Approved Symbol": "gene_symbol"})
    print "----Generate Exp matrix:"
    Generate_Exp_matrix(Add_Symbol_Exp_file, df_Gene, Exp_matrix_file)
    print "----Generate Gene CNA file:"
    Genreate_Gene_CNA_matrix(CNV_file,df_Gene, Gene_CNA_file)
    
    ########## CNV_eQTL analysis
    #### Add TAD info to variants
    print "Add TAD info to variants:"
    df_TAD= pd.read_table(TAD_file,sep="\t")
    CNV_TAD_file = os.path.join(Temp_file_dir, "CNV_TAD_file.txt")
    Add_TAD_info_to_CNV(CNV_file,df_TAD,CNV_TAD_file)

    #### Generate Genotype files
    print "Generate Genotype files:"
    TAD_list = pd.read_table(TAD_file, sep="\t").tad_id.tolist()
    CNV_Genotype_file = os.path.join(Temp_file_dir, "CNV.Genotype.tsv")
    Generate_CNV_Genotype_matrix(TAD_list, CNV_TAD_file, CNV_Genotype_file)

    ####  Generate eQTL input
    print "Generate eQTL input:"
    CNV_eQTL_dir = os.path.join(Temp_file_dir, "CNV_eQTL/")
    CNV_eQTL_Exp_file = os.path.join(CNV_eQTL_dir, "Exp_file.txt")
    CNV_eQTL_Genotype_file = os.path.join(CNV_eQTL_dir, "Genotype_file.txt")
    Preprocess_CNV_Genotype(Exp_matrix_file,Gene_CNA_file,CNV_Genotype_file,CNV_eQTL_Exp_file,CNV_eQTL_Genotype_file)
    #### Matrix eQTL analysis
    print "Matrix eQTL analysis:"
    if os.path.exists(os.path.join(CNV_eQTL_dir,"cis_eqtl.txt")):
        print "CNV cis_eqtl already exists"
    else:
        os.system("Rscript %s/matrix_eQTL.r %s %s" % (os.path.join(Code_dir,"Matrix_eQTL"),CNV_eQTL_dir,Reference_dir))
    os.system("Rscript %s/QQ-plot.r %s %s" % (os.path.join(Code_dir,"Matrix_eQTL"),CNV_eQTL_dir,Reference_dir))
    

    #### Analyze eQTL output
    print "Analyze eQTL output"
    print "----Analyze CNV cis_eQTL:"
    
    CNV_eQTL_Exp_file = os.path.join(CNV_eQTL_dir, "Exp_file.txt")
    CNV_eQTL_Genotype_file = os.path.join(CNV_eQTL_dir, "Genotype_file.txt")
    CNV_cis_eqtl_file = os.path.join(CNV_eQTL_dir, "cis_eqtl.txt")
    boxplot_show_difference(CNV_eQTL_Exp_file, CNV_eQTL_Genotype_file,CNV_cis_eqtl_file)

    ######### SV_eQTL analysis
    #### Add TAD info to variants
    SV_file = os.path.join(SV_dir, re.sub("Project_ID", Project_ID, "structural_somatic_mutation.Project_ID.tsv"))
    SV_TAD_file = os.path.join(Temp_file_dir, "SV_TAD_file.txt")
    Add_TAD_info_to_SV(SV_file,df_TAD,SV_TAD_file)
    #### Generate Genotype file
    SV_Genotype_file = os.path.join(Temp_file_dir, "SV.Genotype.tsv")
    Generate_SV_Genotype_matrix(TAD_list, SV_TAD_file, SV_Genotype_file)
    ####  Generate eQTL input
    SV_eQTL_dir = os.path.join(Temp_file_dir, "SV_eQTL/")
    SV_eQTL_Exp_file = os.path.join(SV_eQTL_dir, "Exp_file.txt")
    SV_eQTL_Genotype_file = os.path.join(SV_eQTL_dir, "Genotype_file.txt")
    Preprocess_SV_Genotype(Exp_matrix_file,Gene_CNA_file, SV_Genotype_file,SV_eQTL_Exp_file,SV_eQTL_Genotype_file)    
    #### Matrix eQTL analysis
    if os.path.exists(os.path.join(SV_eQTL_dir,"cis_eqtl.txt")):
        print "SV cis_eqtl already exists"
    else:
        os.system("Rscript %s/matrix_eQTL.r %s %s" % (os.path.join(Code_dir,"Matrix_eQTL"),SV_eQTL_dir,Reference_dir))
    os.system("Rscript %s/QQ-plot.r %s %s" % (os.path.join(Code_dir,"Matrix_eQTL"),SV_eQTL_dir,Reference_dir))
    #### Analyze eQTL output
    print "----Analyze SV cis_eQTL:"
    SV_eQTL_Exp_file = os.path.join(SV_eQTL_dir, "Exp_file.txt")
    SV_eQTL_Genotype_file = os.path.join(SV_eQTL_dir, "Genotype_file.txt")
    SV_cis_eqtl_file = os.path.join(SV_eQTL_dir, "cis_eqtl.txt")
    boxplot_show_difference(SV_eQTL_Exp_file, SV_eQTL_Genotype_file,SV_cis_eqtl_file)

if __name__ == '__main__':
    main()
