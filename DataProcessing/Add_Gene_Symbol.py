import pandas as pd
import os
def Add_Gene_Symbol(Exp_file, Add_Symbol_Exp_file, df_Gene, gene_id_type):
    if os.path.exists(Add_Symbol_Exp_file):
        print "----!Add_Symbol_Exp_file already exists" 
        return True
        
    output_symbol_list = []
    first_gene_id = ""
    Gene_ID_list = df_Gene[gene_id_type].tolist()
    Gene_Symbol_list = df_Gene["Approved Symbol"].tolist()
    Gene_Symbol_dict = {}
    for i in range(len(Gene_Symbol_list)):
        gene_id = Gene_ID_list[i]
        gene_symbol = Gene_Symbol_list[i]
        Gene_Symbol_dict[gene_id] = gene_symbol
    with open(Exp_file) as exp:
        for eachline in exp:
            if eachline.startswith("icgc"):
                continue
            gene_id = eachline.split("\t")[7]  
            if first_gene_id == gene_id:
                break
            elif first_gene_id == "":
                first_gene_id = gene_id  
            else:
                pass
            if gene_id in Gene_Symbol_dict:
                output_symbol_list.append(Gene_Symbol_dict[gene_id])
            else:
                output_symbol_list.append("noGeneID")
    gene_count = len(output_symbol_list)
    with open(Add_Symbol_Exp_file, "w") as add_exp:
        line_count=-1
        with open(Exp_file) as exp:
            for eachline in exp:
                if eachline.startswith("icgc"):
                    add_exp.write(eachline.strip()+"\tgene_symbol"+"\n")
                else:
                    line_count+=1
                    index_of_gene = line_count%gene_count
                    add_exp.write(eachline.strip()+"\t"+output_symbol_list[index_of_gene]+"\n")
    return True

def main():
    Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/exp_seq.OV-AU.tsv"
    Add_Symbol_Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/OV-AU/exp_seq.Symbol.OV-AU.tsv"
    df_Gene = pd.read_table("../Reference_files/Processed_HGNC.txt",sep="\t")
    Add_Gene_Symbol(Exp_file, Add_Symbol_Exp_file, df_Gene, "Ensembl Gene ID")
if __name__ == '__main__':
    main()