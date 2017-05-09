import pandas as pd
import numpy as np

def Compare_SV_with_CNV_common(df_CNV, df_SV):
    for i in range(len(df_CNV)):
        if df_CNV.chromosome_start.iat[i]>df_CNV.chromosome_end.iat[i]:
            print "Error! The chromosome_start if CNV is bigger than chromosome_end"
    COMMON =0

    df_SV_position = df_SV[["chr_from","chr_from_bkpt","chr_to","chr_to_bkpt"]]
    CNV_dict = {}; CNV_specimen_list = df_CNV.icgc_specimen_id.unique().tolist(); CNV_chrom_list = df_CNV.chromosome.unique().tolist()
    for chrom in CNV_chrom_list:
        CNV_dict[chrom] = {}
        for specimen_id in CNV_specimen_list:
            CNV_dict[chrom][specimen_id]=df_CNV[np.logical_and(df_CNV.chromosome==chrom, df_CNV.icgc_specimen_id==specimen_id)]
    Mutation_type_pair = []
    for i in range(len(df_SV)):
        sv_chr_from, sv_chr_from_bkpt, sv_chr_to, sv_chr_to_bkpt = df_SV_position.iloc[i,:]
        if sv_chr_from_bkpt > sv_chr_to_bkpt:
            sv_chr_from_bkpt, sv_chr_to_bkpt = sv_chr_to_bkpt, sv_chr_from_bkpt
        if sv_chr_from != sv_chr_to:
            #print "Error! The chr_from and chr_to of SV are different!"
            continue

        sv_specimen = df_SV.icgc_specimen_id.iat[i]
        sub_df_CNV = CNV_dict[sv_chr_from][sv_specimen]
        common_variant = sub_df_CNV[np.logical_and(np.absolute(sub_df_CNV.chromosome_start-sv_chr_from_bkpt)<(sv_chr_to_bkpt-sv_chr_from_bkpt)/2,np.absolute(sub_df_CNV.chromosome_end-sv_chr_to_bkpt)<(sv_chr_to_bkpt-sv_chr_from_bkpt)/2)]

        # common_variant = sub_df_CNV[np.logical_and(np.absolute(sub_df_CNV.chromosome_start-sv_chr_from_bkpt)<10000,np.absolute(sub_df_CNV.chromosome_end-sv_chr_to_bkpt)<10000)]
        if len(common_variant) >0:
            if len(common_variant)>1:
                print "Error! more than one CNV is common:",len(common_variant)
            if common_variant.chromosome_end.iat[0]-common_variant.chromosome_start.iat[0] < (sv_chr_to_bkpt-sv_chr_from_bkpt)/2 or common_variant.chromosome_end.iat[0]-common_variant.chromosome_start.iat[0] >2*(sv_chr_to_bkpt-sv_chr_from_bkpt):
                continue
            if i%100==0:
                print i, sv_specimen, sv_chr_from, sv_chr_from_bkpt, sv_chr_to, sv_chr_to_bkpt, common_variant[["chromosome","chromosome_start","chromosome_end"]].iloc[0,:].tolist()
            COMMON+=1
            Mutation_type_pair.append([df_SV.variant_type.iat[i], common_variant.mutation_type.iat[0]])
    df_Muation_type_pair = pd.DataFrame(Mutation_type_pair,columns=["SV_type","CNV_type"])
    df_Muation_type_pair.to_csv("Mutation_type_pair.txt",sep="\t",index=False)
    return COMMON

def Compare_SV_with_CNV(CNV_file,SV_file):
    df_CNV = pd.read_table(CNV_file,sep="\t")
    df_SV = pd.read_table(SV_file, sep="\t")

    # 1. Donor and Specimen
    CNV_donor_list = df_CNV.icgc_donor_id.unique().tolist()
    SV_donor_list = df_SV.icgc_donor_id.unique().tolist()
    common_donor_list = [donor_id for donor_id in CNV_donor_list if donor_id in SV_donor_list]
    print "CNV_donor: %d, SV_donor: %d, common_donor: %d" % (len(CNV_donor_list), len(SV_donor_list), len(common_donor_list))

    CNV_specimen_list = df_CNV.icgc_specimen_id.unique().tolist()
    SV_specimen_list = df_SV.icgc_specimen_id.unique().tolist()
    common_specimen_list = [specimen_id for specimen_id in CNV_specimen_list if specimen_id in SV_specimen_list]
    print "CNV_specimen: %d, SV_specimen: %d, common_specimen: %d" % (len(CNV_specimen_list), len(SV_specimen_list), len(common_specimen_list))
    
    # 2. Derive mutation type
    CNV_mutation_type = df_CNV.mutation_type.unique()
    SV_variant_type = df_SV.variant_type.unique()
    print "CNV_mutation_type includes:",CNV_mutation_type
    print "SV_variant_type includes:",SV_variant_type
    
    ##### 3.0 The overlap between CNV and SV (all kinds of variant)
    df_SV = df_SV[df_SV.chr_from==df_SV.chr_to]
    print "Count of CNV: %d; Count of SV: %d" % (len(df_CNV),len(df_SV))
    count_common_variant = Compare_SV_with_CNV_common(df_CNV, df_SV)
    print "Number of common variant: %d" % count_common_variant

    # ##### 3.1 The overlap between CNV and SV (deletion)
    # df_CNV_loss = df_CNV[df_CNV.mutation_type=="loss"]
    # df_SV_deletion = df_SV[df_SV.variant_type=="deletion"]
    # print "Count of CNV_loss: %d; Count of SV_deletion:%d" % (len(df_CNV_loss), len(df_SV_deletion))
    # count_common_deletion =Compare_SV_with_CNV_common(df_CNV_loss,df_SV_deletion)
    # print "Number of common deletion: %d" % count_common_deletion

    # ##### 3.2 The overlap between CNV and SV (duplication)
    # df_CNV_gain = df_CNV[np.logical_or(df_CNV.mutation_type=="gain", df_CNV.mutation_type=="amp LOH")]
    # df_SV_duplication = df_SV[df_SV.variant_type=="tandem duplication"]
    # print "Count of CNV_gain: %d; Count of SV_duplication:%d" % (len(df_CNV_gain), len(df_SV_duplication))
    # count_common_duplication = Compare_SV_with_CNV_common(df_CNV_gain,df_SV_duplication)
    # print "Number of common duplication: %d" % count_common_duplication

    ##### 3.3 The overlap between CNV and SV (copy neutral)
    # df_CNV_copy_neutral = df_CNV[df_CNV.mutation_type=="copy neutral LOH"]
    # sv_inversion_type = ["inversion","fold-back inversion", "intrachromosomal rearrangement with inverted orientation","intrachromosomal rearrangement with non-inverted orientation"]
    # df_SV_inversion = df_SV[np.array([True if (variant_type in sv_inversion_type) else False for variant_type in df_SV.variant_type])]
    # print "Count of CNV_copy_neutral: %d; Count of SV inversion: %d" % (len(df_CNV_copy_neutral),len(df_SV_inversion))
    # count_common_copy_neutral = Compare_SV_with_CNV_common(df_CNV_copy_neutral,df_SV_inversion)
    # print "number of common copy neutral variants: %d" % count_common_copy_neutral
    

def main():
    CNV_file = "../ExampleData/copy_number_somatic_mutation.OV-AU.tsv"
    SV_file = "../ExampleData/structural_somatic_mutation.OV-AU.tsv"
    Compare_SV_with_CNV(CNV_file,SV_file)

    CNV_mutation_type=['copy neutral LOH','amp LOH','gain','loss']
    SV_variant_type = ['fold-back inversion','interchromosomal rearrangement - unknown type',
                       'intrachromosomal rearrangement with non-inverted orientation',
                       'deletion',
                       'intrachromosomal rearrangement with inverted orientation',
                       'tandem duplication',
                       'inversion']
    COMMON={}
    df_Muation_type_pair = pd.read_table("Mutation_type_pair.txt",sep="\t")
    CNV_mutation_type_list = df_Muation_type_pair.CNV_type.tolist()
    SV_variant_type_list = df_Muation_type_pair.SV_type.tolist()
    for svt in SV_variant_type:
        for cmt in CNV_mutation_type:
            COMMON["----".join([svt,cmt])]=0
    for i in range(len(df_Muation_type_pair)):
        COMMON["----".join([SV_variant_type_list[i],CNV_mutation_type_list[i]])]+=1
    for key in COMMON:
        print "\t".join(key.split("----")+[str(COMMON[key])])

if __name__ == '__main__':
    main()