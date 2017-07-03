import os
import pandas as pd
from time import time

from DataProcessing.Generate_Matrix_File import Generate_Exp_matrix
from DataProcessing.derive_genotype import Generate_Matrix_file_From_bed
# from OutputAnalysis.Analyze_eQTL_output import boxplot_show_difference

# *********Constants**********#
# filenames uesed in this Pipeline
Reference_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Pipeline_Current/Reference_files"
Gene_file = os.path.join(Reference_dir, "gene_loc.txt")
TAD_file = os.path.join(Reference_dir, "TAD_hESC_combined.txt")
Gene_TAD_intersect_file = os.path.join(Reference_dir, "TAD_Gene.intersect.txt")
CNV_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/CNV"
SV_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/SV"
Exp_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression"
Sample_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Specimen"
output_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Selected_project_output/"

Code_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Pipeline_Current"


def Pipeline_CNV_Genotyping(Project_ID):
    Exp_file = os.path.join(Exp_dir, "exp_seq."+Project_ID+".tsv")
    CNV_file = os.path.join(CNV_dir, "copy_number_somatic_mutation."+Project_ID+".tsv")
    Sample_file = os.path.join(Sample_dir, "specimen."+Project_ID+".tsv") 
    work_dir = os.path.join(output_dir, Project_ID)
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    ###### Check if needed files are available
    if not os.path.exists(Exp_file):
        print "----Exp file is not exist----"
    if not os.path.exists(CNV_file):
        print "----CNV file is not exist----"

    ###### Generate Expression matrix file
    print "*********Generate Expression matrix file: Gene_Symbol_list * Donor_list:"
    exp_matrix_file = os.path.join(work_dir,Project_ID+".Exp_matrix.tsv")
    Generate_Exp_matrix(Exp_file,exp_matrix_file,20531) # will output number of patient on screen

    ###### Processing CNV
    print "*********processing CNV:"
    processed_CNV_file = os.path.join(work_dir, Project_ID+".processed_CNV.txt")
    if os.path.exists(processed_CNV_file):
        print "----processed_CNV_file already exists----"
    else:
        os.system("bash ./DataProcessing/Process_CNV.sh %s %s %s"  % (CNV_file, processed_CNV_file, 0.6))

    ###### bedtools intersect TAD_file, gene_file and CNV_file
    print "*********bedtools intersect TAD_file, gene_file and CNV_file:"
    Intersect_file = os.path.join(work_dir, Project_ID+".TAD_GENE_CNV.intersect.txt")
    if os.path.exists(Intersect_file):
        print "----Intersect_file already exists----"
    else:
        os.system("bash ./DataProcessing/generate_cnv_genotype.sh %s %s %s %s" % (TAD_file, Gene_file, processed_CNV_file, Intersect_file))
    
    ###### Generate CNV genotype
    print "*********Generate CNV genotype:"
    tad_genotype_matrix = os.path.join(work_dir, Project_ID+".TAD_genotype_matrix.txt")
    gene_genotype_matrix = os.path.join(work_dir, Project_ID+".Gene_genotype_matrix.txt")
    gene_cna_matrix = os.path.join(work_dir, Project_ID+".Gene_CNA_matrix.txt")
    if os.path.exists(tad_genotype_matrix):
        print "----tad_genotype_matrix already exists----"
    else:
        Generate_Matrix_file_From_bed(Sample_file, TAD_file, Gene_file, Intersect_file, tad_genotype_matrix, gene_genotype_matrix, gene_cna_matrix)

    # ####  Merge_donor_id
    print "*********Merge_donor_id:"
    RemainDonor_dir = os.path.join(work_dir, "RemainDonor/")
    if not os.path.exists(RemainDonor_dir):
        os.mkdir(RemainDonor_dir)
    remain_Exp_matrix = os.path.join(RemainDonor_dir, "Exp_matrix.txt")
    remain_tad_genotype_matrix = os.path.join(RemainDonor_dir, "TAD_genotype_matrix.txt")
    remain_gene_genotype_matrix = os.path.join(RemainDonor_dir, "Gene_genotype_matrix.txt")
    remain_gene_cna_matrix = os.path.join(RemainDonor_dir, "Gene_CNA_matrix.txt")
    
    if os.path.exists(remain_Exp_matrix):
        print "----remaining matrix files already exist----"
    else:
        df_Exp = pd.read_table(exp_matrix_file,sep="\t",index_col=0)
        df_TAD_Genotype = pd.read_table(tad_genotype_matrix,sep="\t",index_col=0)
        df_Gene_Genotype = pd.read_table(gene_genotype_matrix,sep="\t",index_col=0)
        df_Gene_CNA = pd.read_table(gene_cna_matrix,sep="\t",index_col=0)
        common_gene = [gene_id for gene_id in df_Exp.index if gene_id in df_Gene_Genotype.index]
        df_Exp = df_Exp.T[common_gene].T
        df_Gene_Genotype = df_Gene_Genotype.T[common_gene].T
        df_Gene_CNA = df_Gene_CNA.T[common_gene].T
        common_donor = [donor_id for donor_id in df_Exp if donor_id in df_TAD_Genotype]
        df_Exp[common_donor].to_csv(remain_Exp_matrix,sep="\t")
        df_TAD_Genotype[common_donor].to_csv(remain_tad_genotype_matrix,sep="\t")
        df_Gene_Genotype[common_donor].to_csv(remain_gene_genotype_matrix,sep="\t")
        df_Gene_CNA[common_donor].to_csv(remain_gene_cna_matrix,sep="\t")

    # ##### Matrix eQTL analysis
    # print "Matrix eQTL analysis:"
    # if os.path.exists(os.path.join(CNV_eQTL_dir,"cis_eqtl.txt")):
    #     print "CNV cis_eqtl already exists"
    # else:
    #     os.system("Rscript %s/matrix_eQTL.r %s %s" % (os.path.join(Code_dir,"Matrix_eQTL"),CNV_eQTL_dir,Reference_dir))
    # os.system("Rscript %s/QQ-plot.r %s %s" % (os.path.join(Code_dir,"Matrix_eQTL"),CNV_eQTL_dir,Reference_dir))


def main():
    project_code_list = ["BLCA-US","COAD-US","BRCA-US","CESC-US","COAD-US","GBM-US","KIRC-US","KIRP-US","LAML-US",
    "LIHC-US","LUAD-US","LUSC-US","OV-US","PAAD-US","PRAD-US","READ-US","SKCM-US","STAD-US","THCA-US","UCEC-US"]
    # project_code_list = ["BLCA-US"]
    print ">>>>>>>>>>>>Generate matrix files with CNV and exp files"
    for project_code in project_code_list:
        print "****************************"+project_code+"****************************"
        start = time()
        Pipeline_CNV_Genotyping(project_code)
        print "total cost: %d seconds" % (time()-start)
    # merge Exp_matrix, TAD_genotype_matrix, Gene_Genotype_matrix, Gene_CNA_matrix
    print ">>>>>>>>>>>>merge Exp_matrix, TAD_genotype_matrix, Gene_Genotype_matrix, Gene_CNA_matrix"
    total_data_dir = os.path.join(output_dir,"ALL-Project")
    total_Exp_matrix = os.path.join(total_data_dir, "total_Exp_Matrix.txt")
    total_TAD_Genotype_matrix = os.path.join(total_data_dir, "total_TAD_genotype_matrix.txt")
    total_Gene_Genotype_matrix = os.path.join(total_data_dir, "total_Gene_genotype_matrix.txt")
    total_Gene_CNA_matrix = os.path.join(total_data_dir, "total_Gene_CNA_matrix.txt")


    df_Exp_matrix_list = [pd.read_table(os.path.join(output_dir, Project_ID, "RemainDonor", "Exp_matrix.txt"),index_col=0) for Project_ID in project_code_list]
    df_TAD_Genotype_matrix_list = [pd.read_table(os.path.join(output_dir, Project_ID, "RemainDonor","TAD_genotype_matrix.txt"),index_col=0) for Project_ID in project_code_list]
    df_Gene_Genotype_matrix_list = [pd.read_table(os.path.join(output_dir, Project_ID, "RemainDonor","Gene_genotype_matrix.txt"),index_col=0) for Project_ID in project_code_list]
    df_Gene_CNA_matrix_list = [pd.read_table(os.path.join(output_dir, Project_ID, "RemainDonor","Gene_CNA_matrix.txt"),index_col=0) for Project_ID in project_code_list]
    
    df_total_Exp_matrix = pd.concat(df_Exp_matrix_list,axis=1)
    df_total_TAD_Genotype_matrix = pd.concat(df_TAD_Genotype_matrix_list,axis=1)
    df_total_Gene_Genotype_matrix = pd.concat(df_Gene_Genotype_matrix_list,axis=1)
    df_total_Gene_CNA_matrix = pd.concat(df_Gene_CNA_matrix_list,axis=1)

    df_total_Exp_matrix.to_csv(total_Exp_matrix,sep="\t")
    df_total_TAD_Genotype_matrix.to_csv(total_TAD_Genotype_matrix,sep="\t")
    df_total_Gene_Genotype_matrix.to_csv(total_Gene_Genotype_matrix,sep="\t")
    df_total_Gene_CNA_matrix.to_csv(total_Gene_CNA_matrix,sep="\t")



if __name__ == '__main__':
    main()