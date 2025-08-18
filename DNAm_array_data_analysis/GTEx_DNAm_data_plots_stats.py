import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from scipy.stats import pearsonr

#get files
#cg10501210 filter beta values for all samples
file_1 = pd.read_csv("cg10501210_GTEx_data_Oliva_etal.txt")
#metadata - sample IDs
file_2 = pd.read_csv("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",  sep="\t")
#metadata - sample ages
file_3 = pd.read_csv("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", sep='\t')

#file 1 transpose
file_1= file_1.T.reset_index()

#make the first row the header 
file_1.columns = file_1.iloc[0]
#drop the repeatedheader from first row
file_1 = file_1[1:]

#rename cg10501210 to beta values
file_1.rename(columns={"cg10501210": "beta_value"}, inplace=True)

#put new column name
file_1.rename(columns={"Unnamed: 0": "sample_ID"}, inplace=True)

#add donor prefix
file_1["donor_prefix"] = file_1["sample_ID"].str.split("-").str[:2].str.join("-")


#merge the tissue types- details
file_1= file_1.merge(file_2[["SAMPID", "SMTS", "SMTSD"]], left_on="sample_ID", right_on="SAMPID", how="left")

#merge age from file 3
file_1= file_1.merge(file_3[["SUBJID", "AGE"]], left_on="donor_prefix", right_on="SUBJID", how="left")

#drop unwanted columns
file_1.drop(columns=["donor_prefix", "SAMPID"], inplace=True)

#name columns
file_1=file_1[["sample_ID", "SUBJID", "SMTS", "SMTSD","AGE", "beta_value"]]

#age midpoint function
def middle_age(age_range):
    #make into string
    age_range=str(age_range)
    #find the numbers in the string
    num= re.findall(r'\d+', age_range)

    top= int(num[0]) 
    bottom= int(num[1])
    mid= (top +bottom) / 2
    return mid

#calculate midpoint for pearson r
file_1['Midpoint_age']= file_1['AGE'].apply(middle_age)
file_1["beta_value"] = pd.to_numeric(file_1["beta_value"])

#remove heart sample - has only 1 sample
file_1= file_1[file_1["SMTS"] != "Heart"]
#save
file_1.to_csv("GSE213478_metadata_GTEx_DNAm_samples.csv", index=False)



##PLOTTING ##
#colour tissues individually
tissues= file_1["SMTS"].unique()
age_order = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89"]
colours = plt.get_cmap("tab20", len(tissues))
#prep for plot axes
n_rows = (len(tissues)) // 3
n_cols=3
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 5))
axes = axes.flatten()
plt.subplots_adjust(hspace=0.5, wspace=0.3)

#violin plot for age ranges 
for i,j in enumerate(tissues):
    tissue_df= file_1[file_1["SMTS"]== j]

    ax = axes[i]

    sns.violinplot(x="Midpoint_age",
                   y="beta_value",
                   data=tissue_df,
                   inner="box",
                   cut=0,
                   density_norm="width",
                   color=colours(i / len(tissues)),
                   ax=ax)
    
    if len(tissue_df) > 1:
        r, p = pearsonr(tissue_df['Midpoint_age'], tissue_df['beta_value'])
    else:
        r, p = float('nan'), float('nan')
    
    ax.set_title(f'Methylation of cg10501210 in {j}\nr={r:.2f}, p={p:.2e}', fontweight='bold')
    ax.set_xlabel('Age Group (Years)')
    ax.set_ylabel('DNAm')

plt.tight_layout()
plt.savefig("GTEx_methyl_violin_plot.png", dpi=300)
plt.close()



#reset axes
n_rows = (len(tissues) + 2) // 3
n_cols = 3
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 5))
axes = axes.flatten()
#regplot and pearson r
for i,j in enumerate(tissues):
    tissue_df= file_1[file_1["SMTS"]== j]

    ax = axes[i]

    sns.regplot(x="Midpoint_age",
                   y="beta_value",
                   data=tissue_df,
                   color=colours(i),
                   ax=ax)
    
    if len(tissue_df) > 1:
        r, p = pearsonr(tissue_df['Midpoint_age'], tissue_df['beta_value'])
    else:
        r, p = float('nan'), float('nan')
    
    ax.set_title(f'Methylation of cg10501210 in {j}\nr={r:.2f}, p={p:.2e}', fontweight='bold')
    ax.set_xlabel('Age (Midpoint Years)')
    ax.set_ylabel('DNAm')

plt.tight_layout()
plt.savefig("GTEx_methyl_regplot_plot.png", dpi=300)
plt.close()

#create csv for pearsons results
results=[]
for tissue, pearr in file_1.groupby("SMTS"):
    r,p = pearsonr(pearr["Midpoint_age"], pearr["beta_value"])
    results.append({"Tissue": tissue, "Pearsons r": r, "p_value": p})
                    
pearsr_res = pd.DataFrame(results)
pearsr_res.to_csv("GTEx_methyl_pearson_r_results.csv")