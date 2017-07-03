# All are necessary constants
Chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

Chrom_size = {"1":249250621,"2":243199373,"3":198022430,"4":191154276,"5":180915260,"6":171115067,"7":159138663,
"X":155270560,"8":146364022,"9":141213431,"10":135534747,"11":135006516,"12":133851895,"13":115169878,
"14":107349540,"15":102531392,"16":90354753,"17":81195210,"18":78077248,"20":63025520,"Y":59373566,"19":59128983,
"22":51304566,"21":48129895,"M":16571}

ICGC_filename_dict = {
	"SSM":"simple_somatic_mutation.open.project_code.tsv.gz",
    "CNSM":"copy_number_somatic_mutation.project_code.tsv.gz",
    "StSM":"structural_somatic_mutation.project_code.tsv.gz",
    "SGV":"",
    "EXP-A":"exp_array.project_code.tsv.gz",
    "EXP-S":"exp_seq.project_code.tsv.gz",
    "METH-A":"meth_array.project_code.tsv.gz",
    "METH-S":"meth_seq.project_code.tsv.gz",
    "miRNA-S":"mirna_seq.project_code.tsv.gz",
    "PEXP":"protein_expression.project_code.tsv.gz",
    "JCN":"splice_variant.project_code.tsv.gz"
}

Reference_file_list = {
    "Gene_file":"./Reference_files/Processed_HNGC.txt"
    "TAD_file": "./Reference_files/TAD hESC Combined.txt"
}

project_code_list = ["BLCA-US","COAD-US","BRCA-US","CESC-US","COAD-US","GBM-US","KIRC-US","KIRP-US","LAML-US",
    "LIHC-US","LUAD-US","LUSC-US","OV-US","PAAD-US","PRAD-US","READ-US","SKCM-US","STAD-US","THCA-US","UCEC-US"]
