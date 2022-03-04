setwd("K:/Collaborations/Erika_Urdaneta/PTex/20160608_HEK_PTex_LFQ_triplicate/MQ_1.5.1.2/Analyses")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(VennDiagram)
library(gplots)
library(data.table)


#load proteinGroups table

PG <- read.csv("../txt/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE)

#filter out contaminants, reverse and only identified by site
PG <- PG[PG$Reverse != "+",]
#PG <- PG[PG$Potential.contaminant != "+",]         #don't remove contaminants yet
PG <- PG[PG$Only.identified.by.site != "+",]

intensities <- c("Intensity.Input_NoCl_1",
              "Intensity.Input_NoCl_2",
              "Intensity.Input_NoCl_3",
              "Intensity.Input_wCl_1",
              "Intensity.Input_wCl_2",
              "Intensity.Input_wCl_3",
              "Intensity.PTex_NoCl_1",
              "Intensity.PTex_NoCl_2",
              "Intensity.PTex_NoCl_3",
              "Intensity.PTex_wCl_1",
              "Intensity.PTex_wCl_2",
              "Intensity.PTex_wCl_3" )

LFQ <- c("LFQ.intensity.Input_NoCl_1",
      "LFQ.intensity.Input_NoCl_2",
      "LFQ.intensity.Input_NoCl_3",
      "LFQ.intensity.Input_wCl_1",
      "LFQ.intensity.Input_wCl_2",
      "LFQ.intensity.Input_wCl_3",
      "LFQ.intensity.PTex_NoCl_1",
      "LFQ.intensity.PTex_NoCl_2",
      "LFQ.intensity.PTex_NoCl_3",
      "LFQ.intensity.PTex_wCl_1",
      "LFQ.intensity.PTex_wCl_2",
      "LFQ.intensity.PTex_wCl_3")


IBAQ<-c("iBAQ.Input_NoCl_1",
        "iBAQ.Input_NoCl_2",
        "iBAQ.Input_NoCl_3",
        "iBAQ.Input_wCl_1",
        "iBAQ.Input_wCl_2",
        "iBAQ.Input_wCl_3",
        "iBAQ.PTex_NoCl_1",
        "iBAQ.PTex_NoCl_2",
        "iBAQ.PTex_NoCl_3",
        "iBAQ.PTex_wCl_1",
        "iBAQ.PTex_wCl_2",
        "iBAQ.PTex_wCl_3")

#transform intensities and ratios to Log2 (or 10 if you prefer)
PG[c(intensities, IBAQ,LFQ)] = log2(PG[c(intensities,IBAQ,LFQ)])
# change Inf values for na
is.na(PG[c(intensities,IBAQ,LFQ)]) <- sapply(PG[c(intensities,IBAQ,LFQ)], is.infinite)



# how many proteins were identified (have intensity values not NA) in each L-H group pair
for (i in 1:length(IBAQ))
  print(nrow(PG[!is.na(PG[IBAQ[i]]),]))

rm(i)

#nrow(intersect(PG[!is.na(PG[IBAQ[1:3]]),], PG[!is.na(PG[IBAQ[4:6]]),]))
#nrow(intersect(PG[!is.na(PG[IBAQ[7:9]]),], PG[!is.na(PG[IBAQ[10:12]]),]))



## boxplot 
library(reshape2)

#Must remove duplicaed gene names caused by isoforms, if want to show them individually.

#select genes you want to see patterns
proteins = c("CON__P00761")           #Trypsin
#select corresponding colors
cl = c(rgb(0.85,0,0))

for (k in c("IBAQ", "LFQ")) {
  
  all <- melt(PG,id.vars=c("Gene.names", "Majority.protein.IDs"), measure.vars=eval(parse(text = k)))
  
  tiff(paste("Boxplot_",k,".tiff", sep = ""), width = 800, height = 1000, pointsize = 25)
  
  p1 <- ggplot(all, aes(variable,value), col(Gene.names)) +     # Plot boxplot
    geom_boxplot(width = 1) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(face = "bold", color = "black",
                                      angle=90, 
                                      vjust=0.5,
                                      hjust = 1,
                                      size=25),
          axis.title.y = element_text(face="bold",
                                      size=40,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(color = "black", size = 30)) +
    ylab(paste(k," (Log2)", sep = "")) +
    coord_cartesian(ylim = c(12, 40)) +
    scale_y_continuous(breaks=seq(12, 40, 4)) +
    scale_x_discrete(labels= eval(parse(text = k)))
  
  for (j in 1:length(proteins)) {     #add single proteins data points
    vec = c()
    for (i in c(1:length(eval(parse(text = k)))))
      vec = c(vec, PG[PG$Majority.protein.IDs == proteins[j], eval(parse(text = k))[i]])
    p1 <- p1 + geom_point(data = data.frame(vec), x = 1:length(eval(parse(text = k))), y = vec, col = cl[j],size = 4)
  }
  for (j in c(1:length(proteins))) {     #add single proteins data lines
    vec = c()
    for (i in c(1:length(eval(parse(text = k)))))
      vec = c(vec, PG[PG$Majority.protein.IDs == proteins[j], eval(parse(text = k))[i]])
    p1 <- p1 + geom_path(data = data.frame(vec), x = 1:length(eval(parse(text = k))), y = vec, col = cl[j], size = 2, linetype = 1)
  }
  
  plot(p1)
  
  graphics.off()
   
}

rm(cl,i,j,k,p1,proteins,vec)





#### Trypsin normalization ####
#Calculate Trypsin normalization factor = replicate means Input - PTex
trypsin = PG[PG$Majority.protein.IDs == "CON__P00761",]

trp_mean_Inp_NoCl <- mean(trypsin[,LFQ[1]],trypsin[,LFQ[2]],trypsin[,LFQ[3]])
trp_mean_Inp_wCl <- mean(trypsin[,LFQ[4]],trypsin[,LFQ[5]],trypsin[,LFQ[6]])
trp_mean_PTex_NoCl <- mean(trypsin[,LFQ[7]],trypsin[,LFQ[8]],trypsin[,LFQ[9]])
trp_mean_PTex_wCl <- mean(trypsin[,LFQ[10]],trypsin[,LFQ[11]],trypsin[,LFQ[12]])

nf_NoCl <- trp_mean_PTex_NoCl - trp_mean_Inp_NoCl
nf_wCl <- trp_mean_PTex_wCl - trp_mean_Inp_wCl



#Subtract normalization factor from PTex columns
PG <- PG[PG$Potential.contaminant != "+",]
norm_PG <- PG
norm_PG[LFQ[7:9]] <- norm_PG[LFQ[7:9]] - nf_NoCl
norm_PG[LFQ[10:12]] <- norm_PG[LFQ[10:12]] - nf_wCl

















#### Calculate mean intensities ####
norm_PG$LFQ.mean_Input_NoCl <- rowMeans(norm_PG[LFQ[1:3]], na.rm = T)
norm_PG$LFQ.mean_Input_wCl <- rowMeans(norm_PG[LFQ[4:6]], na.rm = T)
norm_PG$LFQ.mean_PTex_NoCl <- rowMeans(norm_PG[LFQ[7:9]], na.rm = T)
norm_PG$LFQ.mean_PTex_wCl <- rowMeans(norm_PG[LFQ[10:12]], na.rm = T)






#### Fold change calculation ####

norm_PG$FC.Input_1 <- norm_PG$LFQ.intensity.Input_wCl_1 - norm_PG$LFQ.intensity.Input_NoCl_1
norm_PG$FC.Input_2 <- norm_PG$LFQ.intensity.Input_wCl_2 - norm_PG$LFQ.intensity.Input_NoCl_2
norm_PG$FC.Input_3 <- norm_PG$LFQ.intensity.Input_wCl_3 - norm_PG$LFQ.intensity.Input_NoCl_3

norm_PG$FC.PTex_1 <- norm_PG$LFQ.intensity.PTex_wCl_1 - norm_PG$LFQ.intensity.PTex_NoCl_1
norm_PG$FC.PTex_2 <- norm_PG$LFQ.intensity.PTex_wCl_2 - norm_PG$LFQ.intensity.PTex_NoCl_2
norm_PG$FC.PTex_3 <- norm_PG$LFQ.intensity.PTex_wCl_3 - norm_PG$LFQ.intensity.PTex_NoCl_3

FC <- c("FC.Input_1", "FC.Input_2", "FC.Input_3", "FC.PTex_1", "FC.PTex_2", "FC.PTex_3")




#Calculate the mean FC
norm_PG$FC.mean_Input <- rowMeans(select(norm_PG, contains("FC.Input_")), na.rm = T)
norm_PG$FC.mean_PTex <- rowMeans(select(norm_PG, contains("FC.PTex_")), na.rm = T)







#### Filter for those identified in at least 2 replicates (2 fold change values) ####
sub_PG <- select(norm_PG, contains("FC.Input_"))
x <- apply (sub_PG, 1, function(x) sum(!is.na(x)) > 1)
norm_PG <- norm_PG[x,]

sub_PG <- select(norm_PG, contains("FC.PTex_"))
x <- apply (sub_PG, 1, function(x) sum(!is.na(x)) > 1)
norm_PG <- norm_PG[x,]

rm(sub_PG,x)













#### pvalue and FDR from the three FCs ####
pval =  eBayes(lmFit(norm_PG[FC[1:3]]))
norm_PG$pval.Input <- pval$p.value  
norm_PG$padj.Input <-  p.adjust(norm_PG$pval.Input, method="BH")

pval =  eBayes(lmFit(norm_PG[FC[4:6]]))
norm_PG$pval.PTex <- pval$p.value  
norm_PG$padj.PTex <-  p.adjust(norm_PG$pval.PTex, method="BH")









##Vulcano plot
#Input
input <- ggplot(data = norm_PG, aes(x = FC.mean_Input, y = -log10(pval.Input))) +
         geom_point() +
  
  coord_cartesian(xlim = c(-4, 12), ylim = c(0, 8)) +
  scale_y_continuous(breaks=seq(0, 8, 1)) + 
  scale_x_continuous(breaks=seq(-4, 12, 2)) + 
  
  ggtitle("Input") +
  xlab("wCL/NoCL Fold Change (log2)") +
  ylab("pValue (-log10)") +
  theme_bw() +
  theme(plot.title = element_text(face="bold",
                                  size = 35,
                                  hjust = 0.5,
                                  vjust = 0.4),
        axis.title.x = element_text(face="bold",
                                    size=30,
                                    hjust = 0.5,
                                    vjust = 0.4),
        axis.text.x  = element_text(face = "bold", color = "black",
                                    angle=0, 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=25),
        axis.title.y = element_text(face="bold",
                                    size=30,
                                    hjust = 0.5,
                                    vjust = 1.5),
        axis.text.y  = element_text(color = "black", size = 25),
        panel.grid=element_blank())

PTex <- ggplot(data = norm_PG, aes(x = FC.mean_PTex, y = -log10(pval.PTex))) +
  geom_point() +

  coord_cartesian(xlim = c(-4, 12), ylim = c(0, 8)) +
  scale_y_continuous(breaks=seq(0, 8, 1)) + 
  scale_x_continuous(breaks=seq(-4, 12, 2)) + 
  
  ggtitle("PTex") +
  xlab("wCL/NoCL Fold Change (log2)") +
  ylab("pValue (-log10)") +
  theme_bw() +
  theme(plot.title = element_text(face="bold",
                                  size = 35,
                                  hjust = 0.5,
                                  vjust = 0.4),
        axis.title.x = element_text(face="bold",
                                    size=30,
                                    hjust = 0.5,
                                    vjust = 0.4),
        axis.text.x  = element_text(face = "bold", color = "black",
                                    angle=0, 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=25),
        axis.title.y = element_text(face="bold",
                                    size=30,
                                    hjust = 0.5,
                                    vjust = 1.5),
        axis.text.y  = element_text(color = "black", size = 25),
        panel.grid=element_blank())


png("VulcanoPlots.png", width = 1000, height = 500, pointsize = 25)
grid.arrange(input, PTex, ncol = 2)
dev.off()
rm(input, PTex)




