.libPaths("")
library(dplyr)

#input from bash script
args <- commandArgs(trailingOnly = TRUE)
wd<-args[1]
rs_id<-args[2]


#set wd
setwd(wd)

#Read data in 
geno <- read.csv(paste0(rs_id,"_genotype.tsv"), header=TRUE, sep="\t")

phased <- read.csv(paste0(rs_id,"_haplotype_assigned.tsv"), header=TRUE, sep="\t")

meth <- read.csv("methyl_fraction_combined.tsv", header=TRUE, sep="\t")

#Remove ungrounded raeds from methyl data
meth<- meth[meth$Haplotype != "ungrouped", ]


#merge geno and phased first to remove real NA
combined <- geno %>%
  left_join(phased, by ="Sample_ID")
#remove NAs no genotype called from bam
combined <- combined[!is.na(combined$Genotype), ]

#esnure haplotype is character
combined$Haplotype <- as.character(combined$Haplotype)

merged <- combined %>%
  left_join(meth, by = c("Sample_ID", "Haplotype"))

#merge
merged<- combined %>%
  left_join(meth, by=c("Sample_ID", "Haplotype"))
            
#remove NAs for merged
df<- merged[!is.na(merged$fraction_modified),]


#Remove unwanted columns
df<- df[ , -c(18:24)]

##add genotype class column 
#genotype lable column
df <- df %>%
  mutate(Genotype_class = case_when(
    Genotype == "0/0" ~ paste0(REF,REF),
    Genotype == "0/1" ~ paste0(REF,ALT),
    Genotype == "1/1" ~ paste0(ALT,ALT),
    TRUE ~ NA_character_
  ))

#assign homo ref to allele colunm dont change if allele is already hap determined
#only 00 is not hp det but same allele on each hp same for 1/1
df <- df %>%
  mutate(Allele = case_when(
    is.na(Allele) & Genotype == "0/0" ~ REF,
    is.na(Allele) & Genotype == "1/1" ~ ALT,
    TRUE ~ Allele
  ))                            

#save to file
write.table(df, "combined_methylation_genotype.tsv", sep="\t", quote=FALSE, row.names=FALSE)
