import pandas as pd
import os

# def 
from constants import *

def main():
    df_All_ICGC = pd.read_table("/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/How_to_download_all/ICGC_all_data.txt")
    for i in range(len(df_All_ICGC)):
    	if df_All_ICGC.CNSM.iloc[i] !="--" and df_All_ICGC["EXP-S"].iloc[i] !="--":
    		print df_All_ICGC.Code.iloc[i]

if __name__ == '__main__':
	main()