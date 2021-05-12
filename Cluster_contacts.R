
# ============================================================================================
# Cluster_contacts.R 
# Developed by Ildefonso Cases and Ana Rojas. 2021. Bioinformatics Unit, CABD-CSIC
# ============================================================================================

# This will cluster different novel dimmers based on common PISA contacts using by mapping positions to a multiple sequence alignment.
# Requires a previous mapping of contacts to a multiple sequence alignment file.
# ============================================================================================


library(tidyverse)

#read the data
data <- read_tsv("Inputs/Contacts_clustering/Clustering/PISAtoMSA.txt")

#get the names
data_wd <- data %>%
  mutate(present=1) %>%
  pivot_wider(names_from = GR,values_from = present )  %>%
  mutate_at(2:41,~if_else(is.na(.),0,1))

my_matrix<-as.matrix(data_wd[,2:41])
rownames(my_matrix)<-data_wd$PosALN

#process and extract data

comunes<-function(x,y){sum(my_matrix[,x] & my_matrix[,y])}
no_comunes_y_no_comunes<- function(x,y){sum(my_matrix[,x] ==1 | my_matrix[,y]==1)}
max_contacts<- function(x,y){max(sum(my_matrix[,x] ==1) ,sum( my_matrix[,y]==1))}
min_contacts<- function(x,y){min(sum(my_matrix[,x] ==1) ,sum( my_matrix[,y]==1))}

comunes_mat<-outer(colnames(my_matrix),colnames(my_matrix),FUN=function(a,b){mapply(comunes,a,b)})
com_y_no_com_mat<-outer(colnames(my_matrix),colnames(my_matrix),FUN=function(a,b){mapply(no_comunes_y_no_comunes,a,b)})
max_mat<-outer(colnames(my_matrix),colnames(my_matrix),FUN=function(a,b){mapply(max_contacts,a,b)})
min_mat<-outer(colnames(my_matrix),colnames(my_matrix),FUN=function(a,b){mapply(min_contacts,a,b)})

#calculate metrics

jaccad_mat<-1-comunes_mat/com_y_no_com_mat
max_dist_mat<-1-comunes_mat/max_mat
min_dist_mat<-1-comunes_mat/min_mat

plotContacts<-function(this_matrix,title){
  rownames(this_matrix)<-colnames(my_matrix)
  hc<-hclust(as.dist(this_matrix))
  pdf(paste0(title,".pdf"), paper = "a4")
  plot(hc, xlab="", main=title,hang=-1,cex=.6,sub='')
  pheatmap(my_matrix,
           cluster_rows=F,
           clustering_distance_cols =  as.dist(this_matrix),
          scale = 'none',
         # breaks=c(0,0.5,1),
          border_color = 'green',
          color=c('grey90','tomato'),
          keep.dendro=T,
          
          fontsize_row   = 2,
          fontsize_col =  5,
          cellwidth = 5,
          cellheight = 1.8,
          legend = FALSE,
          #treeheight_col=60,
          #margins=c(10,10),
          cutree_cols = 6
           )
  dev.off()
}

#Generate PDF files using different metrics

plotContacts(jaccad_mat,"Jaccard")
plotContacts(max_dist_mat,"Max distance")
plotContacts(min_dist_mat,"Min distance")

