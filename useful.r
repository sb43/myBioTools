
# function to plot VAF and depth 
plot_secondary_axis <- function(d,file_path){
# example input for this function
 # V1  V2    V3
 #chr1:104271 123 0.110
  
pdf(file_path,paper='a4r')

par(mar = c(5,5,2,5))
with(d, plot(seq(1,nrow(d)), V3, type="h", col="red3", 
             ylab="VAF", xlab="position",
             ylim=c(min(d[,3]),max(d[,3]))))

#printf <- function(...)print(sprintf(...))

par(new = T)
with(d, plot(seq(1,nrow(d)), V2, pch=16, axes=F, xlab=NA, ylab=NA, cex=0.8))
axis(side = 4)
mtext(side = 4, line = 3, 'Depth')
legend("topleft",
       legend=c(paste0("VAF - Mean(",round(mean(d[,3]),digit=2),")"),paste0("Depth-Mean(",round(mean(d[,2]),digit=2),")")),
       lty=c(1,0), pch=c(NA, 16), col=c("red3", "black"))

dev.off()

}

#================================= plot trinucleotide ==================
plot_trinuc<-function(myfile,outfile_name) {

oesophagus_out<-read.table(myfile,header=T,sep="\t")
library("ggplot2")
library("dplyr")
#library("tidyverse")
library(extrafont)
library(readr)  
###
#pdf(file = paste0(outfile_name,".pdf"))
 
oesophagus_out2 <- filter(oesophagus_out, MutationType != "Total")

#create desired order
oesophagus_out2$MutationType <- factor(oesophagus_out2$MutationType, levels = c("A[C>A]A",
                                                                                "A[C>A]C",
                                                                                "A[C>A]G",
                                                                                "A[C>A]T",
                                                                                "C[C>A]A",
                                                                                "C[C>A]C",
                                                                                "C[C>A]G",
                                                                                "C[C>A]T",
                                                                                "G[C>A]A",
                                                                                "G[C>A]C",
                                                                                "G[C>A]G",
                                                                                "G[C>A]T",
                                                                                "T[C>A]A",
                                                                                "T[C>A]C",
                                                                                "T[C>A]G",
                                                                                "T[C>A]T",
                                                                                "A[C>G]A",
                                                                                "A[C>G]C",
                                                                                "A[C>G]G",
                                                                                "A[C>G]T",
                                                                                "C[C>G]A",
                                                                                "C[C>G]C",
                                                                                "C[C>G]G",
                                                                                "C[C>G]T",
                                                                                "G[C>G]A",
                                                                                "G[C>G]C",
                                                                                "G[C>G]G",
                                                                                "G[C>G]T",
                                                                                "T[C>G]A",
                                                                                "T[C>G]C",
                                                                                "T[C>G]G",
                                                                                "T[C>G]T",
                                                                                "A[C>T]A",
                                                                                "A[C>T]C",
                                                                                "A[C>T]G",
                                                                                "A[C>T]T",
                                                                                "C[C>T]A",
                                                                                "C[C>T]C",
                                                                                "C[C>T]G",
                                                                                "C[C>T]T",
                                                                                "G[C>T]A",
                                                                                "G[C>T]C",
                                                                                "G[C>T]G",
                                                                                "G[C>T]T",
                                                                                "T[C>T]A",
                                                                                "T[C>T]C",
                                                                                "T[C>T]G",
                                                                                "T[C>T]T",
                                                                                "A[T>A]A",
                                                                                "A[T>A]C",
                                                                                "A[T>A]G",
                                                                                "A[T>A]T",
                                                                                "C[T>A]A",
                                                                                "C[T>A]C",
                                                                                "C[T>A]G",
                                                                                "C[T>A]T",
                                                                                "G[T>A]A",
                                                                                "G[T>A]C",
                                                                                "G[T>A]G",
                                                                                "G[T>A]T",
                                                                                "T[T>A]A",
                                                                                "T[T>A]C",
                                                                                "T[T>A]G",
                                                                                "T[T>A]T",
                                                                                "A[T>C]A",
                                                                                "A[T>C]C",
                                                                                "A[T>C]G",
                                                                                "A[T>C]T",
                                                                                "C[T>C]A",
                                                                                "C[T>C]C",
                                                                                "C[T>C]G",
                                                                                "C[T>C]T",
                                                                                "G[T>C]A",
                                                                                "G[T>C]C",
                                                                                "G[T>C]G",
                                                                                "G[T>C]T",
                                                                                "T[T>C]A",
                                                                                "T[T>C]C",
                                                                                "T[T>C]G",
                                                                                "T[T>C]T",
                                                                                "A[T>G]A",
                                                                                "A[T>G]C",
                                                                                "A[T>G]G",
                                                                                "A[T>G]T",
                                                                                "C[T>G]A",
                                                                                "C[T>G]C",
                                                                                "C[T>G]G",
                                                                                "C[T>G]T",
                                                                                "G[T>G]A",
                                                                                "G[T>G]C",
                                                                                "G[T>G]G",
                                                                                "G[T>G]T",
                                                                                "T[T>G]A",
                                                                                "T[T>G]C",
                                                                                "T[T>G]G",
                                                                                "T[T>G]T"))
ggplot()
ggplot(oesophagus_out2, aes(x = MutationType, y = Total, fill = MutationType)) +
  geom_bar(stat = "identity") +
  #facet_wrap(~ data,ncol = 2, scales = "free_y") +
  scale_fill_manual(values=c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue",
                             "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue",
                             "black", "black", "black", "black", "black", "black", "black", "black",
                             "black", "black", "black", "black", "black", "black", "black", "black",
                             "red", "red", "red", "red", "red", "red", "red", "red",
                             "red", "red", "red", "red", "red", "red", "red", "red",
                             "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey",
                             "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey",
                             "dark green", "dark green", "dark green", "dark green", "dark green", "dark green", "dark green", "dark green",
                             "dark green", "dark green", "dark green", "dark green", "dark green", "dark green", "dark green", "dark green",
                             "pink", "pink", "pink", "pink", "pink", "pink", "pink", "pink",
                             "pink", "pink", "pink", "pink", "pink", "pink", "pink", "pink")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.text=element_text(size=6,family="Arial")) +
  ggtitle(outfile_name) +
  xlab("") +
  ylab("")
#ggsave(paste0(outfile_name,".svg"))
# working version wihtout following line
#ggsave(filename=paste0(outfile_name,".pdf"), width=14, height=10, dpi=1200)
#dev.off()
}

plot_data <- function(infile,col1,col2='SAMPLE') {
  #pdf(paste0("~/Desktop/",col2,"_",col1,".pdf"))
  library("ggplot2")
  library("dplyr")
  #library("tidyverse")
  library(extrafont)
  library(readr)  
  
  library(ggplot2)
  df<-read_delim(infile, delim="\t", escape_double = FALSE, trim_ws = TRUE)
  
  # filtering if required
  if(col1 == "CODING") {
    df<-df[df$CODING == 1,]
  }
  
  # violin plot
  print(ggplot(df, aes_string(x=col2, y=col1)) + geom_violin())
  # density plot
  
  print(ggplot(df, aes_string(col1, colour=col2)) + geom_density() + theme(legend.position = c(0.7, 0.3)) )
  # point graph
  print(ggplot(df, aes_string(x=col2, y=col1)) + geom_point())
  # histogram
  
  print(ggplot(df, aes_string(col1) ) + geom_histogram(binwidth = 0.01, colour="black", fill="grey") + facet_wrap(as.formula(paste('~',col2)), ncol = 4,, scales = "free_y"))
  #dev.off()
  
}
