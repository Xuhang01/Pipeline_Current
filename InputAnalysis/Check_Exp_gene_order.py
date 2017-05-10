
def check_order(Exp_file):

    # Check if the exp of each sample contain same gene
    with open(Exp_file) as exp:
        Gene_List_1 = []
        Gene_List_2 = []
        donor_id_1 = ""
        donor_id_2 = ""
        for eachline in exp:
            if eachline.startswith("icgc"):
                continue
            donor_id = eachline.strip().split("\t")[0]
            gene_id = eachline.strip().split("\t")[7]
            if donor_id_1 == "":
                donor_id_1 = donor_id
                Gene_List_1.append(gene_id)
            elif donor_id_1 == donor_id:
                Gene_List_1.append(gene_id)
            elif donor_id_2 == "":
                donor_id_2 = donor_id
                Gene_List_2.append(gene_id)
            elif donor_id_2 == donor_id:
                Gene_List_2.append(gene_id)
            else:
                break
    print "len of gene list in two samples:",len(Gene_List_1), len(Gene_List_2)

    # Check if the order of gene in each sample is the same
    print "the order of gene in each sample is the same: ", Gene_List_2==Gene_List_1

    # Check if the gene list in each sample is the same
    with open(Exp_file) as exp:
        Gene_List = []
        donor_id_current = ""
        for eachline in exp:
            if eachline.startswith("icgc"):
                continue
            donor_id = eachline.strip().split("\t")[2]
            gene_id = eachline.strip().split("\t")[7]
            if donor_id_current == "":
                donor_id_current = donor_id
                Gene_List.append(gene_id)
            elif donor_id_current == donor_id:
                Gene_List.append(gene_id)
            else:
                if len(Gene_List) != len(Gene_List_1):
                    print donor_id_current, len(Gene_List)
                Gene_List=[gene_id]
                donor_id_current = donor_id
    
def main():
    Exp_file = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Expression/exp_seq.BRCA-US.tsv"
    check_order(Exp_file)


if __name__ == '__main__':
    main()