import os

from Code.Merge_Gene_TAD import Merge_Gene_TAD
from Code.Generate_Matrix_File import Prepare_matrix_file
from Code.generate_eQTL_input import Preprocess_CNV_Genotype, Preprocess_SV_Genotype
from Code.Analyze_eQTL_output import boxplot_show_difference


def main():
    #### 1. Merge Gene and TAD file
    Gene_file = "./ExampleData/Processed_ENSEMBL gene.txt"
    TAD_file = "./ExampleData/TAD hESC Combined.txt"
    Gene_TAD_file = "./Temp_file/ENSEMBL_gene_with_TADinfo.txt"
    Merge_Gene_TAD(Gene_file, TAD_file, Gene_TAD_file)

    #### 2. Generate Matrix files
    print "2. Generate Matrix files"
    CNV_file = "./ExampleData/copy_number_somatic_mutation.OV-AU.tsv"
    SV_file = "./ExampleData/structural_somatic_mutation.OV-AU.tsv"
    Exp_file = "./ExampleData/exp_seq.OV-AU.tsv"

    SV_TAD_file = "./Temp_file/Add_TAD_file/SV_TAD_file.OV-AU.txt"
    CNV_TAD_file = "./Temp_file/Add_TAD_file/CNV_TAD_file.OV-AU.txt"

    CNV_matrix_file = "./Temp_file/matrix_file/OV-AU.cnv_matrix.tsv"
    Exp_matrix_file = "./Temp_file/matrix_file/OV-AU.exp_matrix.tsv"
    SV_Genotype_file = "./Temp_file/matrix_file/OV-AU.SV.Genotype.tsv"
    CNV_Genotype_file = "./Temp_file/matrix_file/OV-AU.CNV.Genotype.tsv"
    Prepare_matrix_file(SV_file, CNV_file, Exp_file, SV_TAD_file,CNV_TAD_file,SV_Genotype_file, CNV_Genotype_file, CNV_matrix_file, Exp_matrix_file)
    
    #### 3. Generate eQTL input
    print "3. Generate eQTL input"
    CNV_eQTL_dir = "./Temp_file/CNV_eQTL/"
    CNV_eQTL_Exp_file = os.path.join(CNV_eQTL_dir, "OV-AU.Exp_file.txt")
    CNV_eQTL_Genotype_file = os.path.join(CNV_eQTL_dir, "OV-AU.Genotype_file.txt")
    SV_eQTL_dir = "./Temp_file/SV_eQTL/"
    SV_eQTL_Exp_file = os.path.join(SV_eQTL_dir, "OV-AU.Exp_file.txt")
    SV_eQTL_Genotype_file = os.path.join(SV_eQTL_dir, "OV-AU.Genotype_file.txt")
    Preprocess_CNV_Genotype(Exp_matrix_file,CNV_matrix_file,CNV_Genotype_file,CNV_eQTL_Exp_file,CNV_eQTL_Genotype_file)
    Preprocess_SV_Genotype(Exp_matrix_file,CNV_matrix_file, SV_Genotype_file,SV_eQTL_Exp_file,SV_eQTL_Genotype_file)

    #### 4. Matrix eQTL analysis
    print "4. Matrix eQTL analysis"
    os.system("Rscript ./Code/matrix_eQTL.r %s" % CNV_eQTL_dir)
    os.system("Rscript ./Code/QQ-plot.r %s" % CNV_eQTL_dir)

    os.system("Rscript ./Code/matrix_eQTL.r %s" % SV_eQTL_dir)
    os.system("Rscript ./Code/QQ-plot.r %s" % SV_eQTL_dir)

    #### 5. Analyze eQTL output
    print "5. Analyze eQTL output"
    print "----Analyze CNV cis_eQTL:"
    
    CNV_eQTL_Exp_file = os.path.join(CNV_eQTL_dir, "OV-AU.Exp_file.txt")
    CNV_eQTL_Genotype_file = os.path.join(CNV_eQTL_dir, "OV-AU.Genotype_file.txt")
    CNV_cis_eqtl_file = os.path.join(CNV_eQTL_dir, "cis_eqtl.txt")
    boxplot_show_difference(CNV_eQTL_Exp_file, CNV_eQTL_Genotype_file,CNV_cis_eqtl_file)
    print "----Analyze SV cis_eQTL:"
    SV_eQTL_Exp_file = os.path.join(SV_eQTL_dir, "OV-AU.Exp_file.txt")
    SV_eQTL_Genotype_file = os.path.join(SV_eQTL_dir, "OV-AU.Genotype_file.txt")
    SV_cis_eqtl_file = os.path.join(SV_eQTL_dir, "cis_eqtl.txt")
    boxplot_show_difference(SV_eQTL_Exp_file, SV_eQTL_Genotype_file,SV_cis_eqtl_file)
if __name__ == '__main__':
    main()