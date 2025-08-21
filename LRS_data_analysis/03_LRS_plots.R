#libraies
.libPaths("")

library(ggplot2)
library(dplyr)
library(patchwork)


#input from bash script
args <- commandArgs(trailingOnly = TRUE)
wd= args[1]
data= args[2]
rs_id=args[3]
cpg=args[4]


###Data Preparation ###
#setwd
setwd(wd)

#read methylation data
df<-read.csv(data, sep='\t')

##------------------ FILTER AND SUBSET FOR PLOTS------------------------- ##
## GENOTYPE AND DNA METHYLATION SUBSET##
#data table has two rows per sample(2 each HP3)
#grab only genotype from one of the sample rows
df_sub<-  df %>%
  filter(Haplotype =="1") 

#set up order of genotypes for plot
ref_allele <- df_sub %>%
  filter(Genotype %in% c("0/0", "0/1", "1/0"), REF != ".", !is.na(REF)) %>%
  pull(REF) %>% unique() %>% .[1]

alt_allele <- df_sub %>%
  filter(Genotype %in% c("1/1", "0/1", "1/0"), ALT != ".", !is.na(ALT)) %>%
  pull(ALT) %>% unique() %>% .[1]


#define order 00, 11, 01
allele_order<- c(
  paste0(ref_allele,ref_allele),
  paste0(ref_allele,alt_allele),
  paste0(alt_allele,alt_allele))
df_sub$Genotype_class <- factor(df_sub$Genotype_class, levels = allele_order)

#set colours for plot
colours <- setNames(
  c("#4DAF4A","#FF7F00", "#984EA3"),
  allele_order)

##------------------ PLOTS------------------------- ##

p1 <- ggplot(df_sub, aes(x = Genotype_class, y = AVG_METHYLATION, fill = Genotype_class)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(
    title = paste("Methylation of", cpg, "by","\n", rs_id, "Genotype"),
    x = "Genotype",
    y = "DNAm (%)",
    fill = "Genotype") +
  scale_fill_manual(values = colours) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(filename=paste0(rs_id,"_GT-plot.png"),
       plot = last_plot(),
       device = "png",
       path=wd,
       width=7,
       height=5,
       dpi=300)

## anova or t-tets for significnce -save results to csv
if (length(unique(df_sub$Genotype_class)) == 3) {
  anova_res1 <- aov(AVG_METHYLATION ~ Genotype_class, data = df_sub)
  anova_export <- summary(anova_res1)[[1]]
  write.csv(anova_export, file = file.path(wd, paste0(rs_id, "_methylation_anova_results.csv")), row.names = TRUE)
} else {
  # T test as now only 2 genotypes
  t_result <- t.test(AVG_METHYLATION ~ Genotype_class, data = df_sub)
  
  # Create table from t-test results and save in csv
  t_table <- data.frame(
    tstat = t_result$statistic,
    DF = t_result$parameter,
    pval = t_result$p.value,
    mean_G1 = t_result$estimate[1],
    mean_G2 = t_result$estimate[2],
    confi_Lower = t_result$conf.int[1],
    confi_Upper = t_result$conf.int[2]
  )
  
  write.csv(t_table, file = file.path(wd, paste0(rs_id, "_methylation_t_test.csv")), row.names = FALSE)
}

##PLOT WITH COVERAGE 10##
#create subset over overall cpg calls to 10 
SUB10<- df_sub[df_sub$score >= 10,]

#Plot
sub10p<- ggplot(SUB10, aes(x=Genotype_class, y=AVG_METHYLATION, fill = Genotype_class)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values = colours) +
  labs(
    title = paste("Methylation of", cpg, "by","\n", rs_id, "Genotype"),
    x="Genotype",
    y="DNAm (%)",
    fill= "Genotype") +
  theme(plot.title= element_text(size=12, hjust=0.5))

ggsave(filename=paste0(rs_id,"_GT-plot-coverage10.png"),
       plot = last_plot(),
       device = "png",
       path=wd,
       width=7,
       height=5,
       dpi=300)

if (length(unique(SUB10$Genotype_class)) == 3) {
  anova_res1 <- aov(AVG_METHYLATION ~ Genotype_class, data = SUB10)
  anova_export <- summary(anova_res1)[[1]]
  write.csv(anova_export, file = file.path(wd, paste0(rs_id,"_methylation_anova_resultsSUB10.csv")), row.names = TRUE)
} else {
  # T test as now only 2 genotypes
  t_result <- t.test(AVG_METHYLATION ~ Genotype_class, data = SUB10)
  
  # Create table from t-test results and save in csv
  t_table <- data.frame(
    tstat = t_result$statistic,
    DF = t_result$parameter,
    pval = t_result$p.value,
    mean_G1 = t_result$estimate[1],
    mean_G2 = t_result$estimate[2],
    confi_Lower = t_result$conf.int[1],
    confi_Upper = t_result$conf.int[2]
  )
  
  write.csv(t_table, file = file.path(wd,paste0(rs_id,"_methylation_t_test_SUB10.csv")), row.names = FALSE)
}


##PLOT WITH COVERAGE 12##
#create subset over overall cpg calls to 12
SUB12<- df_sub[df_sub$score >= 12,]

#Plot
sub12p<- ggplot(SUB12, aes(x=Genotype_class, y=AVG_METHYLATION, fill = Genotype_class)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(
    values = colours) +
  labs(
    title = paste("Methylation of", cpg, "by","\n", rs_id, "Genotype"),
    x="Genotype",
    y="DNAm (%)",
    fill="Genotype") +
    theme(plot.title = element_text(hjust=0.5))

ggsave(filename=paste0(rs_id,"_GT-plot-coverage12.png"),
       plot = last_plot(),
       device = "png",
       path=wd,
       width=7,
       height=5,
       dpi=300)

if (length(unique(SUB12$Genotype_class)) == 3) {
  anova_res1 <- aov(AVG_METHYLATION ~ Genotype_class, data = SUB12)
  anova_export <- summary(anova_res1)[[1]]
  write.csv(anova_export, file = file.path(wd, paste0(rs_id,"_methylation_anova_resultsSUB12.csv")), row.names = TRUE)
} else {
  # T test as now only 2 genotypes
  t_result <- t.test(AVG_METHYLATION ~ Genotype_class, data = SUB12)
  
  # Create table from t-test results and save in csv
  t_table <- data.frame(
    tstat = t_result$statistic,
    DF = t_result$parameter,
    pval = t_result$p.value,
    mean_G1 = t_result$estimate[1],
    mean_G2 = t_result$estimate[2],
    confi_Lower = t_result$conf.int[1],
    confi_Upper = t_result$conf.int[2]
  )
  
  write.csv(t_table, file = file.path(wd,paste0(rs_id,"_methylation_t_test_SUB12.csv")), row.names = FALSE)
}


##PLOT WITH COVERAGE 15##
#create subset over overall cpg calls to 15
SUB15<- df_sub[df_sub$score >= 15,]

#Plot
sub15p<- ggplot(SUB15, aes(x=Genotype_class, y=AVG_METHYLATION, fill = Genotype_class)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(
    values = colours) +
  labs(
    title = paste("Methylation of", cpg, "by","\n", rs_id, "Genotype"),
    x="Genotype",
    y="DNAm (%)",
    fill="Genotype")+
    theme(plot.title = element_text(hjust=0.5))

ggsave(filename=paste0(rs_id,"_GT-plot-coverage15.png"),
       plot = last_plot(),
       device = "png",
       path=wd,
       width=7,
       height=5,
       dpi=300)

if (length(unique(SUB15$Genotype_class)) == 3) {
  anova_res1 <- aov(AVG_METHYLATION ~ Genotype_class, data = SUB15)
  anova_export <- summary(anova_res1)[[1]]
  write.csv(anova_export, file = file.path(wd, paste0(rs_id,"_methylation_anova_resultsSUB15.csv")), row.names = TRUE)
} else {
  # T test as now only 2 genotypes
  t_result <- t.test(AVG_METHYLATION ~ Genotype_class, data = SUB15)
  
  # Create table from t-test results and save in csv
  t_table <- data.frame(
    tstat = t_result$statistic,
    DF = t_result$parameter,
    pval = t_result$p.value,
    mean_G1 = t_result$estimate[1],
    mean_G2 = t_result$estimate[2],
    confi_Lower = t_result$conf.int[1],
    confi_Upper = t_result$conf.int[2]
  )
  
  write.csv(t_table, file = file.path(wd,paste0(rs_id,"_methylation_t_test_SUB15.csv")), row.names = FALSE)
}

##PLOT WITH COVERAGE 20##
#create subset over overall cpg calls to 15
SUB20<- df_sub[df_sub$score >= 20,]

#Plot
sub20p<- ggplot(SUB20, aes(x=Genotype_class, y=AVG_METHYLATION, fill = Genotype_class)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(
    values = colours) +
  labs(
    title = paste("Methylation of", cpg, "by","\n", rs_id, "Genotype"),
    x="Genotype",
    y="DNAm (%)",
    fill="Genotype") +
    theme(plot.title = element_text(hjust=0.5))

ggsave(filename=paste0(rs_id,"_GT-plot-coverage20.png"),
       plot = last_plot(),
       device = "png",
       path=wd,
       width=7,
       height=5,
       dpi=300)


if (length(unique(SUB20$Genotype_class)) == 3) {
  anova_res1 <- aov(AVG_METHYLATION ~ Genotype_class, data = SUB20)
  anova_export <- summary(anova_res1)[[1]]
  write.csv(anova_export, file = file.path(wd, paste0(rs_id,"_methylation_anova_resultsSUB20.csv")), row.names = TRUE)
} else {
  # T test as now only 2 genotypes
  t_result <- t.test(AVG_METHYLATION ~ Genotype_class, data = SUB20)
  
  # Create table from t-test results and save in csv
  t_table <- data.frame(
    tstat = t_result$statistic,
    DF = t_result$parameter,
    pval = t_result$p.value,
    mean_G1 = t_result$estimate[1],
    mean_G2 = t_result$estimate[2],
    confi_Lower = t_result$conf.int[1],
    confi_Upper = t_result$conf.int[2]
  )
  
  write.csv(t_table, file = file.path(wd,paste0(rs_id,"_methylation_t_test_SUB20.csv")), row.names = FALSE)
}

## coverage test panel ##
###PLOT VARYING COVERAGE PLOTS ON GRID###
#Add titles to individual plots
plot1 <- sub10p + ggtitle("Coverage 10")
plot2 <- sub12p + ggtitle("Coverage 12")
plot3 <- sub15p + ggtitle("Coverage 15")


combined_plot <- plot1 + plot2 + plot3 + 
  plot_annotation(title = paste("Methylation of", cpg, "by","\n", rs_id, "Genotype"),
                  theme=theme(
                    plot.title = element_text(
                      hjust=0.5, 
                      size=12, 
                      face="bold", 
                      margin=margin(b=20)))) +
  plot_layout(guides = "collect")

ggsave(filename=paste0(rs_id,"_GT-plot-COMBINED.png"),
       plot = last_plot(),
       device = png,
       path=wd,
       width=10,
       height=5,
       dpi=300)


##------------------ FILTER AND SUBSET FOR PLOTS------------------------- ##
## HAPLOTYPE AND DNA METHYLATION SUBSET##
#get rid of any samples where allele coulod not be assigned to HP
df <- df %>% 
  filter(!is.na(Allele))

#set order for HP specifc allele plot
HP_specific_order<- c(
  paste0(ref_allele),
  paste0(alt_allele))
df$Allele <- factor(df$Allele, levels = HP_specific_order)



#set colours for HP
HPcolours <- setNames(
  c("#E41A1C", "#377EB8"),
  HP_specific_order)

##PLOT HAPLOTYPE AND METHYLATION#
hp<- ggplot(df, aes(x=Allele, y=fraction_modified, fill = Allele)) +
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_fill_manual(values= HPcolours) +
  labs(
    title = paste(rs_id,"Allele-Specific Methylation of",cpg ),
    x="Allele",
    y="DNAm (%)"
  )+ 
  theme(plot.title = element_text(hjust = 0.5))

#save
ggsave(filename=paste0(rs_id,"_all_spec_DNAm.png"),
       plot = last_plot(),
       device = png,
       path=wd,
       width=7,
       height=5,
       dpi=300)

#T test as 2 haplotypes
t_result <- t.test(fraction_modified ~ Allele, data = df)
#create table from res and save in csv
t_table_hp <- data.frame(
  tstat = t_result$statistic,
  DF = t_result$parameter,
  pval = t_result$p.value,
  mean_G1 = t_result$estimate[1],
  mean_G2 = t_result$estimate[2],
  confi_Lower = t_result$conf.int[1],
  confi_Upper = t_result$conf.int[2]
)
write.csv(t_table_hp, file=file.path(wd, paste0(rs_id,"_all_spec_ttest_.csv")), row.names = FALSE)

#PLOT WITH COVERAGE 10
SUB10_GA<- df[df$score >= 10,]
#ensure data is paired haplotypes 
SUB10_GA <- SUB10_GA %>%
  group_by(Sample_ID) %>%
  filter(all(c("1", "2") %in% Haplotype)) %>%
  ungroup()

SUB10_GA$Allele <- factor(SUB10_GA$Allele, levels = HP_specific_order)
hp2<- ggplot(SUB10_GA, aes(x=Allele, y=fraction_modified, fill = Allele)) +
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_fill_manual(values= HPcolours) +
  labs(
    title = paste(rs_id, "Allele-Specific Effect on","\n",cpg, "Methylation"),
    x="Allele",
    y="DNAm (%)"
  ) + theme(plot.title = element_text(size=12))

#save
ggsave(filename=paste0(rs_id,"_all_spec_DNAm_sub10.png"),
       plot = last_plot(),
       device = png,
       path=wd,
       width=7,
       height=5,
       dpi=300)
#T test as 2 haplotypes
t_resulthp1 <- t.test(fraction_modified ~ Allele, data = SUB10_GA)
#create table from res and save in csv
t_table_hp10 <- data.frame(
  tstat = t_resulthp1$statistic,
  DF = t_resulthp1$parameter,
  pval = t_resulthp1$p.value,
  mean_G1 = t_resulthp1$estimate[1],
  mean_G2 = t_resulthp1$estimate[2],
  confi_Lower = t_resulthp1$conf.int[1],
  confi_Upper = t_resulthp1$conf.int[2]
)
write.csv(t_table_hp10, file=file.path(wd, paste0(rs_id,"_all_spec_ttest_sub10.csv")), row.names = FALSE)

#PLOT WITH COVERAGE 12
SUB12_GA<- df[df$score >= 12,]
#ensure data is paired haplotypes 
SUB12_GA <- SUB12_GA %>%
  group_by(Sample_ID) %>%
  filter(all(c("1", "2") %in% Haplotype)) %>%
  ungroup()

SUB12_GA$Allele <- factor(SUB12_GA$Allele, levels = HP_specific_order)
hp3<- ggplot(SUB12_GA, aes(x=Allele, y=fraction_modified, fill = Allele)) +
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_fill_manual(values= HPcolours) +
  labs(
    title = paste(rs_id, "Allele-Specific Effect on","\n",cpg, "Methylation"),
    x="Allele",
    y="DNAm (%)"
  ) + theme(plot.title = element_text(hjust=0.5))

ggsave(filename=paste0(rs_id,"_all_spec_DNAm-COVERAGE12.png"),
       plot = last_plot(),
       device = png,
       path=wd,
       width=7,
       height=5,
       dpi=300)

#T test as 2 haplotypes
t_resulthp1 <- t.test(fraction_modified ~ Allele, data = SUB12_GA)
#create table from res and save in csv
t_table_hp12 <- data.frame(
  tstat = t_resulthp1$statistic,
  DF = t_resulthp1$parameter,
  pval = t_resulthp1$p.value,
  mean_G1 = t_resulthp1$estimate[1],
  mean_G2 = t_resulthp1$estimate[2],
  confi_Lower = t_resulthp1$conf.int[1],
  confi_Upper = t_resulthp1$conf.int[2]
)
write.csv(t_table_hp12, file=file.path(wd, paste0(rs_id,"_all_spec_ttest_sub12.csv")), row.names = FALSE)

#PLOT WITH COVERAGE 15
SUB15_GA<- df[df$score >= 15,]
#ensure data is paired haplotypes 
SUB15_GA <- SUB15_GA %>%
  group_by(Sample_ID) %>%
  filter(all(c("1", "2") %in% Haplotype)) %>%
  ungroup()

SUB15_GA$Allele <- factor(SUB15_GA$Allele, levels = HP_specific_order)
hp4<- ggplot(SUB15_GA, aes(x=Allele, y=fraction_modified, fill = Allele)) +
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_fill_manual(values= HPcolours) +
  labs(
    title = paste(rs_id, "Allele-Specific Effect on","\n",cpg, "Methylation"),
    x="Allele",
    y="DNAm (%)"
  ) + theme(plot.title = element_text(hjust=0.5))

ggsave(filename=paste0(rs_id,"_all_spec_DNAm-COVERAGE15.png"),
       plot = last_plot(),
       device = png,
       path=wd,
       width=7,
       height=5,
       dpi=300)

#T test as 2 haplotypes
t_resulthp1 <- t.test(fraction_modified ~ Allele, data = SUB15_GA)
#create table from res and save in csv
t_table_hp15<- data.frame(
  tstat = t_resulthp1$statistic,
  DF = t_resulthp1$parameter,
  pval = t_resulthp1$p.value,
  mean_G1 = t_resulthp1$estimate[1],
  mean_G2 = t_resulthp1$estimate[2],
  confi_Lower = t_resulthp1$conf.int[1],
  confi_Upper = t_resulthp1$conf.int[2]
)
write.csv(t_table_hp15, file=file.path(wd, paste0(rs_id,"_all_spec_ttest_sub15.csv")), row.names = FALSE)


#GRIB PLOT FOR HAPLOTYPE
hptitle <- hp2 + ggtitle("Coverage 10")
hp2title <- hp3 + ggtitle("Coverage 12")
hp3title<- hp4 + ggtitle("coverage 15")

combined_plot_HP <- hptitle + hp2title + hp3title +
  plot_annotation(
    title = paste(rs_id, "Allele-Specific Effect on","\n",cpg, "Methylation"),
    tag_levels ='A',
    theme = theme(
      plot.title = element_text(
        hjust = 0.5, 
        size = 12, 
        face = "bold", 
        margin = margin(b = 20)
      ))) +
  plot_layout(guides = "collect")

ggsave(filename=paste0(rs_id,"_all_spec-plot-COMBINED.png"),
       plot = last_plot(),
       device = png,
       path=wd,
       width=15,
       height=5,
       dpi=300)

sub10p <- sub10p +ggtitle(paste("Methylation of", cpg, "by","\n", rs_id, "Genotype")) + labs(tag = "C")
hp2 <- hp2 + ggtitle(paste(rs_id, "Allele-Specific Effect on","\n",cpg, "Methylation")) + labs(tag = "D")


## panel a 
panela<- sub10p + hp2 +
  plot_annotation(
    #tag_levels = 'A',
    title = "1KG ONT Data- 100 Plus Dataset")

ggsave(filename=paste0(rs_id, "_GT_all_spec-plot-COMBINED_cov10.png"),
       plot = panela,
       device = png,
       width=10,
       height=5,
       dpi=300)


sub12p <- sub12p +ggtitle(paste("Methylation of", cpg, "by","\n", rs_id, "Genotype")) + labs(tag = "C")
hp3 <- hp3 + ggtitle(paste(rs_id, "Allele-Specific Effect on","\n",cpg, "Methylation")) + labs(tag = "D")


## panel a 
panela<- sub12p + hp3 +
  plot_annotation(
    #tag_levels = 'A',
    title = "1KG ONT Data- 100 Plus Dataset")

ggsave(filename=paste0(rs_id, "_GT_all_spec-plot-COMBINED_cov12.png"),
       plot = panela,
       device = png,
       width=10,
       height=5,
       dpi=300)

#data used to plot genotype and DNAm 
write.csv(SUB10, file = file.path(wd, paste0(rs_id, "_SUB10.csv")), row.names = FALSE)
write.csv(SUB12, file = file.path(wd, paste0(rs_id, "_SUB12.csv")), row.names = FALSE)
write.csv(SUB15, file = file.path(wd, paste0(rs_id, "_SUB15.csv")), row.names = FALSE)
write.csv(SUB20, file = file.path(wd, paste0(rs_id, "_SUB20.csv")), row.names = FALSE)

#data used for allele specific DNAm
write.csv(SUB10_GA, file = file.path(wd, paste0(rs_id, "_SUB10_ALLELE_SPEC.csv")), row.names = FALSE)
write.csv(SUB12_GA, file = file.path(wd, paste0(rs_id, "_SUB12_ALLELE_SPEC.csv")), row.names = FALSE)
write.csv(SUB15_GA, file = file.path(wd, paste0(rs_id, "_SUB15_ALLELE_SPEC.csv")), row.names = FALSE)


