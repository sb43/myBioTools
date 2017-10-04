# my functions
library(edgeR)
library(gplots)
library(goseq)
library(geneplotter)
library(limma)
# library for clustering using nb
library(MBCluster.Seq)
library(ShrinkBayes)
library(VennDiagram)
library(lessR)
library(RSQLite)
library(Rsamtools)
library(pathview)
library(org.Sc.sgd.db)
# for RNA-seq GSEA analysis
library(gage)


get_RPKM <- function(edgeR_object)
                        {
                        # this function needds edgeR object to calculate RPKM values
                        edgeR_object_rpkm <- (1e+09 * edgeR_object$counts/(expandAsMatrix(edgeR_object$samples$lib.size, dim(edgeR_object)) * expandAsMatrix(edgeR_object$gene$Length, dim(edgeR_object )))) ;
                        return(edgeR_object_rpkm) }


get_TPM <- function(edgeR_object,read_length)
                       {

                                rl <- read_length
                                T <- (edgeR_object$counts*(expandAsMatrix( read_length,dim(edgeR_object) )) / expandAsMatrix( edgeR_object$gene$Length, dim(edgeR_object) ))
                                edgeR_object_tpm <-  ( ( 1e+06 * edgeR_object$counts * rl )/( expandAsMatrix(colSums(T), dim(edgeR_object)) * expandAsMatrix(edgeR_object$gene$Length, dim(edgeR_object) ) )) ;
                        return(edgeR_object_tpm)
                        }


# function to intersect/union lists functions from library lessR ?? check ...
get_union <- function(lists_for_union) Reduce('union', lists_for_union)
get_intersect <- function(lists_for_intersect) Reduce('intersect', lists_for_intersect)

# rank these genes (function adapted from Mick)...
rank_expression <- function(x) {
  x <- as.vector(x)
  r <- rank(x)
  p <- paste(r, sep="", collapse="")
return(p)
}


#
draw_tree <- function(d,heatmap_title)
{

#d <- d[,11:15]
min <-apply(d, 1, min)
max <- apply(d, 1, max)
# order by desccending order...
ordered_data <- d[order(max - min, decreasing=TRUE),]
# calculate distance
calc_dist <- as.dist(1 - cor(t(ordered_data)))
# do clustering
image_title=paste(heatmap_title,"#data_points","[",dim(d)[1],"]",sep="")
hierarchical_clustering <- hclust(calc_dist, "average")
# draw heatmap...

          heatmap.2(as.matrix(ordered_data),
          Rowv=as.dendrogram(hierarchical_clustering),
          Colv=FALSE,
          cexRow=0.1,
          cexCol=1,
          dendrogram="row",
          scale="row",
          trace="none",
          density.info="none",
          key=TRUE,
          col=greenred.colors(80),
          margins=c(8,7),
          cex=0.2,
            main=paste(image_title,strain,sep=":")
          )

}

## function modified from VennDiagram package

# genes should in the form of single column
plot_venn5<- function(A,B,C,D,E,venn_title)
{
venn.plot <- venn.diagram(
x = list(
T1 = A,
T2 = B,
T3 = C,
T4 = D,
T5 = E
),
filename = venn_title,
col = "black",
fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
alpha = 0.50,
cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
cat.cex = 1.5,
cat.fontface = "bold",
margin = 0.05
);

}

plot_venn4<- function(A,B,C,D,venn_title)
{
  venn.plot <- venn.diagram(
  x = list(
  T1 = A,
  T2 = B,
  T3 = C,
  T4 = D
  ),
  filename = venn_title,
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white",
  "white", "white", "white", "white", "darkblue", "white",
  "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif",
  #rotation.degree = 270,
  margin = 0.2
  );

}


plot_venn3<- function(A,B,C,venn_title)
{
venn.plot <- venn.diagram(
x = list(
T1 = A,
T2 = B,
T3 = C
),
euler.d = TRUE,
filename = venn_title,
fill = c("darkblue", "darkgreen", "orange"),
cat.col = c("darkblue", "darkgreen", "orange"),
cex = 2.5,
cat.cex = 1.5,
cat.pos = 0
);
}


##### MBclust dendrogram

do_MBClustering<- function(dat,strain_name,num_of_clusters)
{

Counts <- dat
GeneID=rownames(Counts)
Treatment=c(colnames(Counts))
#Normalizer=NULL, uses log(Q2) by default
mydata=RNASeq.Data(Counts,Normalize=NULL,Treatment,GeneID)
c0=KmeansPlus.RNASeq(mydata,nK=num_of_clusters)$centers
# poisson model is recommended by authors for data without biological  replicates .....
cls=Cluster.RNASeq(data=mydata, model="poisson", centers=c0,method="EM")$cluster
tr=Hybrid.Tree(data=mydata, cluste=cls, model="poisson")


pdf(paste(strain_name,"MBclustering",".pdf",sep=""), height=7, width=14,paper="USr" )

plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=paste(strain_name,"HybridTree",sep=":"),colorful=TRUE)


par_val <- round(num_of_clusters/4,0)
par(mfrow=c(4,par_val), ps = 10, cex.main = 0.8)
#par(oma=c(2,2,2,4))
par(oma=c(0,0,0,0) )
par(mar=c(4,4.5,2,1))
#par(ps = 12, cex.main = 0.8)

#genes in each cluster

#pdf(paste(strain,"MBclustering",".pdf",sep=""))

 for (i in (1:num_of_clusters))
 {
 gene_cluster <- mydata$logFC[cls==i,]
 #filter genes with FC >=1.5
 
 # check if there are values after filtering...
       if(!is.null(row.names(gene_cluster)))
       {
          logFC_filtered_cluster <- gene_cluster[rowSums(abs(gene_cluster) >= 0.001) >0, ]
           #plotlines(gene_cluster, first.column.origin=FALSE, xlab="Timepoint", ylab="logFC", col=rainbow(7), lwd=1, main=paste(sample_name,"cluster:",i,sep=""))
           plotlines(logFC_filtered_cluster, first.column.origin=FALSE, xlab="Timepoint", ylab="logFC", col=rainbow(7), lwd=1, 
           main=paste(strain_name,"_C:",i,"_genes(",length(row.names(logFC_filtered_cluster)),")",sep=""), cex=2)
           # write 
           #write.table(logFC_filtered_cluster,file=paste(strain,"_cluster_",i,"_logFC.tsv",sep=""),sep="\t")
           #write.table(raw_data[row.names(logFC_filtered_cluster), ],file=paste(strain,"_cluster_",i,"_counts.tsv",sep=""),sep="\t")
           write.table(rownames(logFC_filtered_cluster),file=paste(strain_name,"_cluster_gene_list",i,".tsv",sep=""),row.names=F, quote=F, col.names=F )
           
       
       }
       
 }

dev.off()
 }



 # create line graph for a genelist
 
plot_graph <- function(file_name, data_all_rpkm, db_file)
{
print(paste(file_name,db_file))

# timepoints
#for (tp in 1:5)
#{ 
#file_name <- paste(strain_name[x],"_GFOLD_downregulated_genes_T",tp,sep="")

#my_genes<-read.table(paste(file_name,".tsv",sep=""), head=F, fill=T, stringsAsFactors=FALSE,sep="\t")
my_genes<-read.table(paste(file_name,".txt",sep=""), head=F, fill=T, stringsAsFactors=FALSE,sep="\t")
strain_name <-c("LEB1_vs_LEB3" ,"LEB1_vs_LEB2", "LEB2_vs_LEB3") 
# create sql query list
gene_list<-paste(my_genes,sep="")
gene_list <- gsub("c","", gene_list)
drv <- dbDriver("SQLite")
#db_file <- "M:/ingenza/kegg/yeast_kegg.db"
con <- dbConnect(drv, dbname=db_file )
res <- dbSendQuery(con, paste("select distinct ensid,gene from genes where gene in",gene_list," or ensid in",gene_list,"order by gene",sep=""))
#res <- dbSendQuery(con, paste("select ensid,gene from genes where ensid in",gene_list,sep=""))
data <- fetch(res)
my_genes <- data
# end of SQL query...

# added GLA gene...
genes <- c(my_genes[,1],"GLA_CAS")
annotations <- c(my_genes[,2],"Glucoamylase")
#plotlines(log(as.matrix(data_all_rpkm[genes, ]), first.column.origin=FALSE, xlab="Timepoint", ylab="logFC", col=rainbow(7), lwd=1)

pdf(file=paste(file_name,"_tpm.pdf",sep=""),height=7, width=14,paper="USr")


## ############## draw barcharts################
if(FALSE){     # comment if barchart is required
par(mfrow=c(3,4), ps = 12, cex.main = 0.7)
#layout(matrix(c(1:36), 4, 9, byrow = TRUE))

#par(mfrow=c(4,))

new_gene <- ""
 j<-0;
 for (i in 1:length(genes))
 {
 print(genes[i])
   
if(!is.na(data_all_rpkm[genes[i],1 ]))
  {
       j <- j+1
       # plot in next page after these many graphs
        if(j %% 12 == 0)
         {
            par(new=T)
            par(mfrow=c(3,4), ps = 12, cex.main = 0.7)
         }

          #
          fpkm_vector <- unname(unlist(data_all_rpkm[genes[i],]))
          barplot(fpkm_vector,  main=paste(genes[i],annotations[i],sep=":"),
          ylab="TPM_all",
          xlab="Timepoints",
          names.arg=c(rep(c("T1","T2","T3","T4","T5"),3)),
          ylim=c(0,round(max(data_all_rpkm[genes[i], ]))),
          space=c(0.0,0.8),
          col = c(rep("lightblue",5),rep("mistyrose",5),rep("lavender",5)),
          las=2
           )
          legend("topleft", c("LEB1","LEB2","LEB3"), cex=0.7, bty="n", fill=c("lightblue", "mistyrose", "lavender"), xjust=0, y.intersp=1)
        new_gene[[length(new_gene)+1]] <-genes[i]
          
           #legend(locator(1),c("LEB1","LEB2","LEB3"),fill=c("lightblue", "mistyrose", "lavender"))
          #barplot(as.matrix(data_uc_rpkm[genes[i], ]), xlab="Timepoints", ylab="FPKM_uc", main=genes[i], ylim=c(0,round(max(data_all_rpkm[genes[i], ]))), col = gray.colors(1),las=2)
         # legend("topleft", c("LEB1","LEB2","LEB3"), cex=0.6, bty="n", fill=c("lightblue", "mistyrose", "lavender"))
   }
 }

 } # if FASLE
########## END of draw barchart##############
 
 par(mfrow=c(3,4), ps = 12, cex.main = 0.7)
 
 #### Line graph...##################
 
 new_gene <- ""
 j<-0;
line_col= c("red","green","black")
for (i in 1:length(genes))
{
if(!(genes[i] %in% row.names(data_all_rpkm))) {next} 
if(!is.na(round(max(data_all_rpkm[genes[i], ]))) < 10) {next}
print(genes[i])
   
if(!is.na(data_all_rpkm[genes[i],1 ]))
  {
       
       j <- j+1
       # plot in next page after these many graphs
        if(j %% 12 == 0)
         {
            par(new=T)
            par(mfrow=c(3,4), ps = 12, cex.main = 0.7)
         }

          # timepoints
          t1<-1
          t2<-5
                    
          fpkm_vector <- unname(unlist(data_all_rpkm[genes[i],]))
          # plot  line for first strain
          plot(1:5, fpkm_vector[1:5],  main=paste(genes[i],annotations[i],sep=":"),
          ylab="TPM_all",
          xlab="Timepoints",
          ylim=c(0,round(max(data_all_rpkm[genes[i], ]))),
          col =line_col[1] ,
          type="l"
           )
          #plot line for next two strains....
          for (a in 2:length(strain_name))
          {
            t1<-t2+1 
            t2<-5*a
            lines(1:5,fpkm_vector[t1:t2],col =line_col[a])
            
          }
          #legend("topleft", c("LEB1","LEB2","LEB3"), cex=0.7, bty="n", fill=c("red","green","black"), xjust=0, y.intersp=1)
          legend("topleft", c("LEB1","LEB2","LEB3"), cex=0.7, bty="n", fill=c("red","green","black"), xjust=0, y.intersp=1)
          
          new_gene[[length(new_gene)+1]] <-genes[i]
          
           #legend(locator(1),c("LEB1","LEB2","LEB3"),fill=c("lightblue", "mistyrose", "lavender"))
          #barplot(as.matrix(data_uc_rpkm[genes[i], ]), xlab="Timepoints", ylab="FPKM_uc", main=genes[i], ylim=c(0,round(max(data_all_rpkm[genes[i], ]))), col = gray.colors(1),las=2)
         # legend("topleft", c("LEB1","LEB2","LEB3"), cex=0.6, bty="n", fill=c("lightblue", "mistyrose", "lavender"))
   }
 }

} # line graph function...

# analyse_pathways

pathway_enrichment <- function(genes,gene_length_bias,edgeR_n_Gfold_fc,method_name,tm,strain,kegg.gs,tpm_raw,db_file,pathway_name)
{

pwf=nullp(genes,"sacCer3","sgdGene",bias.data=gene_length_bias, plot.fit=FALSE)
KEGG=goseq(pwf,"sacCer3","sgdGene",gene2cat=path2gene,test.cats="KEGG")

# use uncorrected p-val
#enriched.KEGG=KEGG$category[(KEGG$over_represented_pvalue) < .01]
# corrected p-val can be used more strigent....
enriched.KEGG=KEGG$category[p.adjust(KEGG$over_represented_pvalue,method="BH") < 0.05]

# plot expression on pathway...


for (i in enriched.KEGG)
{
  
  #pathwayList_temp <-addToList(paste(tm,"_",method_name,sep=""),pathway_name[i,])
  #pathwayList<- c(pathwayList,pathwayList_temp)
 kegg_path_temp<-KEGG[KEGG$category == paste(i),]
 pathway<- gsub("_"," ",pathway_name[i,])
 
  write.table(paste(tm,kegg_path_temp[,1],pathway,kegg_path_temp[,2],kegg_path_temp[,4],kegg_path_temp[,5], method_name,sep="\t"),file=paste("../pathway_table_",method_name,"_",strain,".tsv",sep=""),append=TRUE,row.names=F, quote=F, col.names=F )
# source("N:/Shriram/project_yeast/results/my_functions.r") 

  pv.out.list <- sapply(i, function(pid) pathview(gene.data = edgeR_n_Gfold_fc[,1:5], pathway.id = pid, res=300, species = "sce", gene.idtype="KEGG", out.suffix =paste(pathway_name[i,],method_name,tm,strain,sep="_"), same.layer = F, kegg.native = T, node.sum="median"))
  kegg.gs[[grep(i,names(kegg.gs))]]
  #write list of genes in enriched KEGG pathway
  write.table(kegg.gs[[grep(i,names(kegg.gs))]],file=paste(pathway_name[i,],".txt",sep=""), row.names=F, quote=F,col.names=F)
  plot_graph(pathway_name[i,],tpm_raw,db_file)
  dev.off()
}

}

addToList <- function(name, value) {
 pathwayList[[name]] <- value
 return (pathwayList)
}

