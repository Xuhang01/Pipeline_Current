import os

def Remove_Duplicate(Exp_file, remove_duplicate_exp_file):
    # 1. Detect the number of genes
    if os.path.exists(remove_duplicate_exp_file):
        print "----!remove_duplicate_exp_file already exists:"
        return True
    Gene_list = []
    with open(Exp_file) as exp:
        for eachline in exp:
            if eachline.startswith("icgc"):
                donor_id = ""
            else:
                if donor_id=="":
                    donor_id = eachline.split("\t")[0]
                elif donor_id != eachline.split("\t")[0]:
                    break
                else:
                    pass
                gene_id = eachline.split("\t")[7]
                Gene_list.append(gene_id)
    set_Gene_list = set(Gene_list)

    Duplicate_Gene_list = []
    print len(Gene_list), len(set_Gene_list)
    if len(Gene_list) == len(set_Gene_list):
        pass
    else:
        for gene_id in Gene_list:
            if Gene_list.count(gene_id)!= 1:
                Duplicate_Gene_list.append(gene_id)
        Duplicate_Gene_list = set(Duplicate_Gene_list)
    print len(Duplicate_Gene_list)
    len_of_remain_gene = len(set_Gene_list)-len(Duplicate_Gene_list)
    with open(remove_duplicate_exp_file,"w") as rexp:
        donor_id = ""
        with open(Exp_file) as exp:
            for eachline in exp:

                if eachline.startswith("icgc"):
                    rexp.write(eachline)
                else:
                    gene_id = eachline.split("\t")[7]
                    current_donor_id = eachline.split("\t")[0]
                    if gene_id in Duplicate_Gene_list:
                        pass
                    elif donor_id != current_donor_id:
                        line_count=1
                        donor_id = current_donor_id
                        rexp.write(eachline)
                    elif donor_id == current_donor_id and line_count < len_of_remain_gene:
                        rexp.write(eachline)
                        line_count+=1
                    else:
                        line_count+=1

    return True


def main():
    Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/exp_seq.BRCA-US.tsv"
    remove_duplicate_exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/remove_duplicate_exp_seq.BRCA-US.tsv"
    Remove_Duplicate(Exp_file, remove_duplicate_exp_file)

    Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/exp_seq.CLLE-ES.tsv"
    remove_duplicate_exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/remove_duplicate_exp_seq.CLLE-ES.tsv"
    Remove_Duplicate(Exp_file, remove_duplicate_exp_file)

if __name__ == '__main__':
    main()