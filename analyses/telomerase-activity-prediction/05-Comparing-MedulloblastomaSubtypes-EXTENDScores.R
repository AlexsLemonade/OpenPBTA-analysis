library('ggpubr')
stringsAsFactors=FALSE
library(grid)
library(forcats) ### for fct_reorder()
library(optparse)
############################################################### Comparing EXTEND Scores of Medulloblastoma molecular subtypes (Figure 4)  ###############################################################################


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

PBTA_EXTEND_MedulloSubtypes <- opt$output




PTBA_Histology = read.table(A,sep='\t',head=T)    ## Reading the clinical data

colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'   ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose

PTBA_GE_Standard_TMScores = read.table(B,sep='\t',head=T)  ## Reading Stranded FPKM telomerase scores

PTBA_GE_Standard_Histology = merge(PTBA_Histology,PTBA_GE_Standard_TMScores,by='SampleID')   ### Merging Clinical data with the Telomerase scores

PTBA_GE_Standard_Histology = PTBA_GE_Standard_Histology[which(PTBA_GE_Standard_Histology$short_histology == "Medulloblastoma"),]   ### Select tumors with catagory labelled as "Medulloblastoma"

my_comparisons = list(c("Group3","Group4"),c("Group3","SHH"),c("Group3","WNT"),c("SHH","WNT"))


pdf(PBTA_EXTEND_MedulloSubtypes)

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


P1 = ggplot(PTBA_GE_Standard_Histology, aes(x=fct_reorder(molecular_subtype,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.2,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="pink")+ geom_jitter(shape=16, width = 0.2)+stat_compare_means(comparisons = my_comparisons, method= "t.test",size=5)

grid.newpage()
# Create layout : nrow = 1, ncol = 1
pushViewport(viewport(layout = grid.layout(nrow =1, ncol = 1)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =c("black",14),font.ytickslab =c("black",14),font.x = 12,font.y=12,font.legend=6,xlab="Molecular Subgroups Medulloblastoma",ylab="EXTEND Scores"),vp = define_region(row = 1, col = 1))


dev.off()


