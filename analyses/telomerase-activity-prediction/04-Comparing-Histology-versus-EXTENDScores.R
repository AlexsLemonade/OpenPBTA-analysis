suppressPackageStartupMessages({
  library(ggpubr)
  library(grid)
  library(forcats) ### for fct_reorder()
  library(optparse)
  library(cowplot)
})
stringsAsFactors=FALSE

############################################################### Combining Histology with EXTEND Scores (Figure 3) ###############################################################################


root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "data", "pbta-histologies.tsv")   ### Variable representing clinical data 
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data


				 
# set up the command line options
option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Plot output (.pdf)")
)


# parse the command line arguments - this becomes a list
# where each element is named by the long flag option
opt <- parse_args(OptionParser(option_list = option_list))

PBTA_EXTEND_HistologyCompPlot <- opt$output




PTBA_Histology = readr::read_tsv(A, guess_max = 10000)    ## Reading the clinical data

colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'   ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose

PTBA_GE_Standard_TMScores = read.table(B,sep='\t',head=T)  ## Reading Stranded FPKM telomerase scores

PTBA_GE_Standard_Histology = merge(PTBA_Histology,PTBA_GE_Standard_TMScores,by='SampleID')   ### Merging Clinical data with the Telomerase scores

Stranded_Histology = PTBA_GE_Standard_Histology
# If there are any "Other" samples, remove them 
if (any(Stranded_Histology$short_histology == "Other")) {
  Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology == "Other"),]   ### Removing the tumors with catagory labelled as "Others"
}
Frequency = data.frame(table(Stranded_Histology$short_histology))  ### Counting the number of cases for all histologies to avoid less number further
colnames(Frequency)=c('Variables','Freq')
Frequency = Frequency[which(Frequency$Freq == 1),]
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology %in% Frequency$Variables),]     ### Removing the tumors with only one case in histologies


pdf(PBTA_EXTEND_HistologyCompPlot)

## Globally set the theme in one step, so it gets applied to both plots
theme_set(theme_classic() +
          theme(plot.title = element_text(size=10, face="bold"),axis.text.x=element_text(angle=50,size=6,vjust=1,hjust=1),axis.text.y=element_text(size=7), axis.title.x = element_text(size=0), axis.title.y = element_text(size=8),
          legend.position = "top",
          legend.key.size= unit(0.3,"cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.title = element_text(size=7),
          legend.text =element_text(size=6)
        )
)


P1 = ggplot(Stranded_Histology, aes(x=fct_reorder(short_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="pink",alpha=0.5)+ geom_jitter(alpha= 0.5,shape=16, cex=0.3)+  ylab("EXTEND Scores")+ ggtitle("Tumor Short Histology")
P2 = ggplot(Stranded_Histology, aes(x=fct_reorder(broad_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="pink",alpha=0.5)+ geom_jitter(alpha=0.5,shape=16, cex=0.3)+ ylab("EXTEND Scores")+ ggtitle("Tumor Broad Histology")

plot_grid(P1,P2, nrow = 2,labels = "AUTO", label_size = 12, scale=c(0.95,0.95),  align = "v")

dev.off()


