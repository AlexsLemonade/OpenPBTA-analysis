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




PTBA_Histology = readr::read_tsv(A, guess_max = 10000)    ## Reading the clinical data

colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'   ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose

PTBA_GE_Standard_TMScores = read.table(B,sep='\t',head=T)  ## Reading Stranded FPKM telomerase scores

PTBA_GE_Standard_Histology = merge(PTBA_Histology,PTBA_GE_Standard_TMScores,by='SampleID')   ### Merging Clinical data with the Telomerase scores

PTBA_GE_Standard_Histology = PTBA_GE_Standard_Histology[which(PTBA_GE_Standard_Histology$short_histology == "Medulloblastoma"),]   ### Select tumors with catagory labelled as "Medulloblastoma"

stat.test <- data.frame(compare_means(
  NormEXTENDScores ~ molecular_subtype, data = PTBA_GE_Standard_Histology,
  method = "t.test"
)) 

combinations = nrow(stat.test)

statistics = stat.test%>%
  mutate(y.position = seq(1,by=0.04,length.out=combinations))


pdf(PBTA_EXTEND_MedulloSubtypes,height = 3, width = 3)

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


P1 = ggplot(PTBA_GE_Standard_Histology, aes(x=fct_reorder(molecular_subtype,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(width=0.5,size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="pink")+ geom_jitter(shape=16, width = 0.1,size=0.4)+stat_pvalue_manual(
    data = statistics, label = "p.adj",size=1.5,
    xmin = "group1", xmax = "group2",tip.length = 0.006,
    y.position = "y.position"
    )


grid.newpage()
# Create layout : nrow = 2, ncol = 1
pushViewport(viewport(layout = grid.layout(nrow =1, ncol = 1)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =c("black",4),font.ytickslab =c("black",5),font.x = 5,font.y=5,font.legend=6,xlab="Molecular Subgroups Medulloblastoma",ylab="EXTEND Scores"),vp = define_region(row = 1, col = 1))


dev.off()


