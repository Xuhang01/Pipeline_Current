import pandas as pd
import os
import argparse
Chrom_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']


def Generate_Exp_matrix(Exp_file,Exp_matrix_file,gene_number):
    # Exp_file is the directly download Exp_seq file from ICGC database
    # Exp_matrix_file is the output of this function
    # df_Gene is a dataframe of gene, columns contain ['chrom', 'start','end','gene_symbol']
    if os.path.exists(Exp_matrix_file):
        print '----!Exp_matrix_file already exists'
        return True
    
    df_Exp = pd.read_csv(Exp_file,sep='\t',iterator=True)
    loop = True
    chunk_size = gene_number
    i=1
    while loop:
        i+=1
        if i%100==1:
            print 'number of patients:', i
        try:
            chunk = df_Exp.get_chunk(chunk_size)
            chunk = chunk.sort_values(by=['gene_id'])
            chunk.index = chunk['gene_id']
            chunk = chunk[chunk['gene_id']!='?']
            chunk = chunk[chunk['gene_id']!='SLC35E2']
            # if len(chunk['icgc_donor_id'].unique())!=1:
            #     print '--improper icgc_donor_id', chunk['icgc_donor_id'].unique()
            #     break
            icgc_donor_id = chunk['icgc_donor_id'].iat[0]
            chunk[icgc_donor_id]=chunk['normalized_read_count']
            if 'df_Exp_Matrix' not in dir():
                try: 
                    df_Exp_Matrix = chunk[[icgc_donor_id]]
                except ValueError:
                    print '--invalid icgc_donor_id:', icgc_donor_id
                    pass
                
            else:
                try: 
                    df_Exp_Matrix[icgc_donor_id]=chunk[icgc_donor_id]
                except ValueError:
                    print '--invalid icgc_donor_id:', icgc_donor_id
                    pass

        except  StopIteration:
            loop = False
            print 'Coverting Exp file to matrix file has been finished'
    df_Exp_Matrix.to_csv(Exp_matrix_file, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_file',help='expression file as input')
    parser.add_argument('--output_file', help='expression matrix file as output')
    arg = parser.parse_args()

    Exp_file = arg.exp_file
    Exp_matrix_file = arg.output_file

    Generate_Exp_matrix(Exp_file, Exp_matrix_file, 20531)

if __name__ == '__main__':
    main()
