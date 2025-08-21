from glob import glob
import pandas as pd
import os.path
import re

#set input directory
in_dir=""
#create empty df for combined 
freqs=[]


#for file in directory 
for file in glob(os.path.join(in_dir, "*MethylFrequency-filt.tsv")):
    
    #get the basename
    basename= os.path.basename(file)
    #get sample id by substituting using regex
    sample= re.sub(r'_NanoMethPhase_HP[12]_MethylFrequency-filt\.tsv$', '', basename)
    sample= sample.split('_', 1)[1]
    
    #read file
    df= pd.read_csv(file, sep="\t")

    #if file is empty still include the ID but put Na in all columns
    if df.empty:
        df_new = pd.DataFrame([{
            'chromosome': pd.NA,
            'start': pd.NA,
            'end': pd.NA,
            'strand': pd.NA,
            'NumOfAllCalls': pd.NA,
            'NumOfModCalls': pd.NA,
            'MethylFreq': pd.NA,
            'ID': sample,
            'Haplotype': 1 if "_HP1_" in basename else 2
        }])
    else:
         #recalualte the methylation frequency to get single frequency for the DNA strand
        NumOfAllCalls= df["NumOfAllCalls"].sum()
        NumOfModCalls= df["NumOfModCalls"].sum()
        MethylFreq= NumOfModCalls/NumOfAllCalls
        #Add a strand column so we know  
        strands_present = df['strand'].unique()
        strand_info = ",".join(strands_present)

        #make new df with calculated MethylFreq
        df_new= pd.DataFrame([{
            'chromosome': df["chromosome"].iloc[0],
            'start': "207823675",
            'end': "207823676",
            'strand': ".",
            'NumOfAllCalls': NumOfAllCalls,
            'NumOfModCalls': NumOfModCalls,
            'MethylFreq': NumOfModCalls/NumOfAllCalls,
            'ID': sample,
            'Haplotype': 1 if "_HP1_" in basename else 2,
            'strand_in_MethylFreq_Calc': strand_info
            }])
        
    freqs.append(df_new)
#combine all to one 
combi=pd.concat(freqs, ignore_index=True)
combi=combi.sort_values('ID', ascending=True)
#adding ages to df from metadat from ONEil paper
ages= pd.read_excel("normal_samples_dups_filtered.xlsx")
combi = combi.merge(ages[['normal_lib_workflow_name', 'consent_age']], left_on='ID', right_on='normal_lib_workflow_name', how='left')
combi.drop(columns=['normal_lib_workflow_name'], inplace=True)

#get midpoint ages
def middle_age(age_range):
    #make into string
    age_range=str(age_range)
    #find the numbers in the string
    num= re.findall(r'\d+', age_range)

    top= int(num[0]) 
    bottom= int(num[1])
    mid= (top +bottom) / 2
    return mid
combi['Midpoint_age']= combi['consent_age'].apply(middle_age)
combi.to_csv("POG_LRS_MethylFreq_HP1-2_cg10501210_F+R.tsv", sep='\t')

    
