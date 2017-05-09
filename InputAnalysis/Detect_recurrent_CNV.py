import pandas as pd

def main():
    CNV_file = "../ExampleData/copy_number_somatic_mutation.OV-AU.tsv"
    CNV_mutation_type=['copy neutral LOH','amp LOH','gain','loss']

    df_CNV = pd.read_table(CNV_file,sep="\t")
    df_CNV = df_CNV[["icgc_specimen_id","mutation_type","chromosome","chromosome_start","chromosome_end"]]
    df_CNV = df_CNV.sort(["chromosome","chromosome_start","chromosome_end"])
    
    print df_CNV.head()
    
    
    
if __name__ == '__main__':
	main()