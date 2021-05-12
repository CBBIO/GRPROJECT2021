
# =============================================
#Sectors_stats.R
#Developed by Ildefonso Cases, Bioinformatics Unit Head, CABD-CSIC. 2020.
#Modified by A. Rojas, CBBio, CABD-CSIC. 2021.
# =============================================

# It calculates the statistical overlap of residues deemed as "sectors" by SCA5 and PysCA using different input alignment.
#requires a summarised table containing wether a residue is a sector or not by a different method. 
#Developed by Ildefonso Cases, Bioinformatics Unit Head, CABD-CSIC. 2020.
#Modified by A. Rojas, CBBio, CABD-CSIC. 2021.

#Requires
#install.packages("writexl")
library(tidyverse)
library(readxl)
library(xlsx)

#Getting the sectors defined by GR01_C2, and data processing

data<-read_xlsx("Inputs/Sector_Analyses/sectors.xlsx") #Check your path here.
data_long<-pivot_longer(data,cols=pySCA:Sca5_onPFAMrp55,names_to="ALING",values_to="sector")
data_long %>% separate(sector,into=c('S1','S2','S3'),sep="\\|",remove = F,fill = "right") -> data_long_s3

data_s1 <- data_long_s3 %>% select(-S2,-S3) %>% rename(Sector=S1)
data_s2 <- data_long_s3 %>% select(-S1,-S3) %>% rename(Sector=S2)
data_s3 <- data_long_s3 %>% select(-S1,-S2) %>% rename(Sector=S3)

data_long_nored<-bind_rows(data_s1,data_s2,data_s3) %>% filter(!is.na(Sector)) %>% filter(Sector!="-") %>%  unique()


res_in_sector <- data_long_nored  %>% select(1:2) %>% unique()
n_res_in_sector<-nrow(res_in_sector)

data_long_nored %>% select(-sector) %>% unite("Sector",ALING,Sector) %>% group_by(Sector) %>% summarise(n=n()) -> res_per_sector
data_long_nored %>% select(-sector) %>% unite("Sector",ALING,Sector) %>% left_join(data_long_nored %>% select(-sector) %>% unite("Sector",ALING,Sector),by=c('res_num',"aa")) %>% group_by(Sector.x,Sector.y) %>% summarise(n=n()) -> overlap_table

overlap_table %>% left_join(res_per_sector,by=c(Sector.x='Sector')) %>% left_join(res_per_sector,by=c(Sector.y='Sector')) %>% rename(overlap=n.x,n_s1=n.y,n_s2=n) -> overlap_table

#calculating  pvals and hypergeometric test
to_pval <- overlap_table %>%
  ungroup %>%
  mutate(total=n_res_in_sector) %>%
  select(-Sector.x,-Sector.y) %>%
  mutate(q=overlap-1,n=total-n_s1) %>%
  rename(m=n_s1,k=n_s2) %>%
  select(-overlap,-total) %>%
  select(q,m,n,k)

pvals<-pmap_dbl(to_pval,phyper,lower.tail=F)

#creating the final table with pvals

overlap_table <- bind_cols(overlap_table,p.val=pvals)
overlap_table <- overlap_table %>% filter(Sector.x!=Sector.y) %>% mutate(adj.p.val=p.adjust(p.val))
write.table(overlap_table, file = "Outputs/Sector_Analyses/sectors_overlap.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)


