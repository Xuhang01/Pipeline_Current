import pandas as pd
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description=__file__.split(".")[0])
    parser.add_argument('--Exp_matrix_file',help='Exp_matrix_file')
    parser.add_argument('--z_score_matrix_file', help='z_score_matrix_file')
    args = parser.parse_args()
    z_score_matrix_file = args.z_score_matrix_file
    Exp_matrix_file = args.Exp_matrix_file
    
    df_Exp = pd.read_csv(Exp_matrix_file, sep="\t", index_col = 0)
    print "Total number of Genes: %d" % len(df_Exp)
    print "Total number of samples: %d" % len(df_Exp.columns)
    
    df_Exp = np.log10(df_Exp+0.00000001)
    
    Sr_mean = np.mean(df_Exp,axis=1)
    Sr_std = np.std(df_Exp, axis=1)

    df_Exp_score = ((df_Exp.T-Sr_mean)/Sr_std).T
    df_Exp_score.to_csv(z_score_matrix_file,sep="\t")

if __name__ == '__main__':
    main()