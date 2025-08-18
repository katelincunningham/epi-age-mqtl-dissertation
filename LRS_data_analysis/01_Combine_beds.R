.libPaths("")
library(dplyr)

#input from bash script
args <- commandArgs(trailingOnly = TRUE)
bed_dir= args[1]
out_dir=args[2]

#gdt all bed files from bed_dir
bed_list<- list.files(path= bed_dir, pattern = "\\.bed$", full.names = TRUE)

#save empty file
empty_beds <- file.path(out_dir, "empty_bed_files.txt")
empty_bed_files <- c()

#get file and extract info fucntion
to_file_w_ID <- function(filex) {
  if (file.info(filex)$size == 0) {
    empty_bed_files <<- c(empty_bed_files, filex) #<<- to assign variable globally not just in function
    return(NULL) #skip the empty files -they cause problems in tsv creation 
  }
  filename <- basename(filex)
  ID_HP <- strsplit(tools::file_path_sans_ext(filename), "_")[[1]]
  bed<- read.table(filex, header = FALSE, sep="\t", stringsAsFactors = FALSE)
  bed$Sample_ID <- ID_HP[1] #first part is ID
  bed$Haplotype <- ID_HP[2] #second is HP 
  return(bed)
}

#now do this for all bed files and conbine them
cbn_bed <- lapply(bed_list, to_file_w_ID) %>%
  bind_rows()

#order and name columns 
colnames(cbn_bed) <- c('chrom', 'start', 'end', 'modified_base_code', 
    'score', 'strand', 'start_2','end_2', 'colour','N_valid_cov', 
    'fraction_modified','Nmod', 'Ncanonical', 'Nother_mod','Ndelete',
    'Nfail', 'Ndiff','Nnocall', 'Sample_ID', 'Haplotype')

#add column for average methylation for genotype plot
cbn_bed <- cbn_bed %>%
  group_by(Sample_ID) %>%
  mutate(AVG_METHYLATION= mean(fraction_modified, na.rm= TRUE)) %>%
  ungroup()


#save table
write.table(cbn_bed, 
            file = file.path(out_dir, "methyl_fraction_combined.tsv"),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


