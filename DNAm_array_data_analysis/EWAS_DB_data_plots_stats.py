import pandas as pd 
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import math

#set input dir
in_dir= ""
#loop through all file in dir
for file in os.listdir(in_dir):
    if file.endswith(".csv"):
        filepath= os.path.join(in_dir, file)
        try:
            print(f"processing {file}")

            df=pd.read_csv(filepath)

            #change male and female columns to sex 
            df_new= df.melt(id_vars="Age (y)", value_vars=["Female", "Male"],
                            var_name="Sex", value_name="Beta")

            #get rid of nas
            df_new= df_new.dropna(subset=["Beta"])

            #sort df
            sorted_df= df_new.sort_values(by="Age (y)")

            #get basename for new file name
            base= os.path.splitext(file)[0]
            out_name= base +"_sorted.csv"
            out_path= os.path.join(in_dir,out_name)

            #add tissue column 
            sorted_df["tissue"]= base
            #save new file
            sorted_df.to_csv(out_path, index=False)


            print(f"done and dusted for {out_name}")
        except Exception as e:
            print(f"error with {file}")
            print(f"{type(e).__name__}: {e}\n")


#get all sorted csv in dir
csv_files = glob.glob("/*_sorted.csv")  

#read each one and pull into one big one
df_list = [pd.read_csv(f) for f in csv_files]
merged_df = pd.concat(df_list)

#save
merged_df.to_csv("EWAS_ATLAS_TISSUES_COMBINED_23_04_2025.csv", index=False, header=True)


## PLOTTING ##
tissues = merged_df["tissue"].unique()
colour = plt.get_cmap("tab20", len(tissues))
n_cols = 3
n_rows = math.ceil(len(tissues) / n_cols)


#set plot axes and size
fig, axes =plt.subplots(n_rows, n_cols, figsize=(16, n_rows * 4))
axes = axes.flatten()

#plot all tissues together
for i, j in enumerate(tissues):
    tissue_df = merged_df[merged_df["tissue"] == j]

    ax = axes[i]

    sns.regplot(
        x="Age (y)",
        y="Beta",
        data=tissue_df,
        color=colour(i),
        ax=ax
    )
    #stats
    if len(tissue_df) > 1:
        r, p = pearsonr(tissue_df['Age (y)'], tissue_df['Beta'])
    else:
        r, p = float('nan'), float('nan')
    
    plt.subplots_adjust(hspace=0.5, wspace=0.3)
    ax.set_title(f'Methylation of cg10501210 in \n{j}\nr={r:.2f}, p={p:.2e}', fontweight='bold')
    ax.set_xlabel('Age')
    ax.set_ylabel('DNAm')
#get rid of empty plots 
for plot in range(i +1, len(axes)):
    fig.delaxes(axes[plot])

plt.savefig("EWAS_cg10501210_me+Age_23_04_2025_results.png", dpi= 300, bbox_inches= 'tight')

#create csv for pearsons results
try:
    results=[]
    for tissue, pearr in merged_df.groupby("tissue"):
        r,p = pearsonr(pearr["Age (y)"], pearr["Beta"])
        results.append({"Tissue": tissue, "Pearsons r": r, "p_value": p})
                        
    pearsr_res = pd.DataFrame(results)
    pearsr_res= pearsr_res.sort_values(by="Pearsons r", ascending=False)
    pearsr_res.to_csv("EWAS_pearson_r_23_04_2025_results.csv")
    print("csv made and saved to folder")
except Exception as e:
    print(f"{type(e).__name__}: {e}\n")
    print(f'{pearsr_res} not saved')