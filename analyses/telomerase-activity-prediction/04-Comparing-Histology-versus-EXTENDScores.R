library('ggpubr')
stringsAsFactors=FALSE
library(grid)
library(forcats) ### for fct_reorder()
library(optparse)
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




PTBA_Histology = read.table(A,sep='\t',head=T)    ## Reading the clinical data

colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'   ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose

PTBA_GE_Standard_TMScores = read.table(B,sep='\t',head=T)  ## Reading Stranded FPKM telomerase scores

PTBA_GE_Standard_Histology = merge(PTBA_Histology,PTBA_GE_Standard_TMScores,by='SampleID')   ### Merging Clinical data with the Telomerase scores

Stranded_Histology = PTBA_GE_Standard_Histology
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology == "Other"),]   ### Removing the tumors with catagory labelled as "Others"
Frequency = data.frame(table(Stranded_Histology$short_histology))  ### Counting the number of cases for all histologies to avoid less number further
colnames(Frequency)=c('Variables','Freq')
Frequency = Frequency[which(Frequency$Freq == 1),]
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology %in% Frequency$Variables),]     ### Removing the tumors with only one case in histologies 


pdf(PBTA_EXTEND_HistologyCompPlot)

## Globally set the theme in one step, so it gets applied to both plots
theme_set(theme_classic() + 
          theme(axis.text.x=element_text(angle=50,size=7,vjust=1,hjust=1),
          legend.position = "top",
          legend.key.size= unit(0.3,"cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.title = element_text(size=7),
          legend.text =element_text(size=6)
        )
)


P1 = ggplot(Stranded_Histology, aes(x=fct_reorder(short_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="cyan3")
P2 = ggplot(Stranded_Histology, aes(x=fct_reorder(broad_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="cyan3")

grid.newpage()
# Create layout : nrow = 2, ncol = 1
pushViewport(viewport(layout = grid.layout(nrow =2, ncol = 1)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =c("black",6),font.ytickslab =c("black",7),font.x = 7,font.y=7,font.legend=6,xlab="Tumor Histology(short)",ylab="EXTEND Scores"),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =c("black",6),font.ytickslab =c("black",7),font.x = 7,font.y=7,font.legend=6,xlab="Tumor Histology(broad)",ylab="EXTEND Scores"),vp = define_region(row = 2, col = 1))


dev.off()


