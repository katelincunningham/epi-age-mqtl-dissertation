library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

#read in data
data<- read.csv('POG_LRS_MethylFreq_HP1-2_cg10501210_F+R.tsv',sep = "\t")


#turn Haplotype into a categorical
data$Haplotype <- as.factor(data$Haplotype)

## FILTER DATA ##
#remove NAs
data <- na.omit(data) 
df_clean <- data %>%
  filter(Haplotype %in% c(1, 2)) %>% 
  group_by(ID) %>%
  filter(n_distinct(Haplotype) == 2) %>% #only haplotype 1 and 2 inlcuded 
  ungroup() %>%
  #only include samples that had readings for forward and reverese strands in recalulation
  filter(strand_in_MethylFreq_Calc %in% c("-,+", "+,-")) 

## AGE AND METHYLATION ##
plot1<- ggplot(df_clean, aes(x=Midpoint_age, y=MethylFreq, colour = Haplotype))+ 
  geom_point()+
  geom_smooth(method=lm) +
  stat_cor(aes(color = Haplotype), method = "pearson",
           label.x.npc = 0.75,
           label.y.npc = "top")+
  labs(title="cg10501210 Methylation with Age",
       x ="Age midpoint", y = "Methylation Frequency") +
  theme(plot.title = element_text(hjust = 0.6))



## Removing matched skin sample ##
#removing matched skin sample to evaluate different on the results
df_clean<- df_clean %>% filter(ID != "F48645")

Aplot<- ggplot(df_clean, aes(x=Midpoint_age, y=MethylFreq, colour = factor(Haplotype)))+ 
  scale_color_brewer(palette = "Set1", name= "Haplotype")+
  geom_point()+
  geom_smooth(method=lm) +
  stat_cor(aes(color = Haplotype), method = "pearson",
           label.x.npc = 0.75,
           label.y.npc = "top")+
  labs(title="cg10501210 Methylation with Age",
       x ="Age midpoint (Years)", y = "DNAm") +
  theme(plot.title = element_text(hjust = 0.6))

ggsave("Age+me_POG_LRS_dataset.png", 
       plot= Aplot, dpi=300, width = 7, height = 5)

## HAPLOTYPE  AND AGE ##
#valid sample that have both haplotypes more than 3 calls valid_3 used here 
HPplot<- ggplot(df_clean, aes(x=Haplotype, y=MethylFreq, fill = Haplotype)) +
  geom_boxplot() +
  #add mean points+ 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_color_brewer(palette = "Set1") +
  labs(title = "cg10501210 Methylation by Haplotype",
       x="Haplotype",
       y="DNAm") +
  theme(plot.title = element_text(hjust = 0.6))


ggsave("HP+me_POG_LRS_dataset.png", 
       plot= HPplot, dpi=300, width = 7, height = 5)


## T-TEST HAPLOTYPE + METHYLATION##
df_clean$Haplotype <- as.factor(df_clean$Haplotype)
t.test(MethylFreq ~ Haplotype, data = df_clean )

#PLOT PANNEL ##
both_plots<- (Aplot / HPplot)+
  plot_annotation(tag_levels = 'A')

ggsave("pannel_age+HP+me_POG_LRS_dataset.png",
       plot=both_plots,
       dpi=300,width = 8, height = 10)