import numpy as np
import pandas as pd


def Generate_Matrix_file_From_bed(sample_file, tad_file, gene_file, intersect_file, tad_genotype_matrix, gene_genotype_matrix, gene_cna_matrix):
    donor_list = pd.read_table(sample_file)["icgc_donor_id"].unique().tolist()
    tad_list = pd.read_table(tad_file,header=None).iloc[:,3].tolist()
    gene_list = pd.read_table(gene_file, header=None).iloc[:,3].tolist()
    df_TAD_matrix = pd.DataFrame(data= np.zeros((len(tad_list),len(donor_list)),dtype=np.int),columns=donor_list, index=tad_list)
    df_GENE_matrix = pd.DataFrame(data= np.zeros((len(gene_list),len(donor_list)),dtype=np.int),columns=donor_list, index=gene_list)
    df_GENE_cna = pd.DataFrame(data= np.zeros((len(gene_list),len(donor_list)),dtype=np.float),columns=donor_list, index=gene_list)

    with open(intersect_file) as intersect:
        i=0
        for eachline in intersect:
            i+=1
            if i %100000 ==0:
                print "Number of records in intersect file:", i
            tad_id,segment_mean, donor_id, specimen_id, sample_id, gene_id, tad_genotype, gene_genotype = eachline.strip().split("\t")
            df_GENE_cna.set_value(gene_id, donor_id, float(segment_mean))
            if float(segment_mean)<0:
                df_TAD_matrix.set_value(tad_id, donor_id, -int(tad_genotype))
                df_GENE_matrix.set_value(gene_id, donor_id, -int(gene_genotype))
            else:
                df_TAD_matrix.set_value(tad_id, donor_id, int(tad_genotype))
                df_GENE_matrix.set_value(gene_id, donor_id, int(gene_genotype))

    df_TAD_matrix.to_csv(tad_genotype_matrix,sep="\t")
    df_GENE_matrix.to_csv(gene_genotype_matrix,sep="\t")
    df_GENE_cna.to_csv(gene_cna_matrix,sep="\t")

def main():
    pass

if __name__ == '__main__':
    main()