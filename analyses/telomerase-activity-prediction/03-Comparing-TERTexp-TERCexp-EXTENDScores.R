suppressPackageStartupMessages({
  library(ggpubr)
  library(grid)
  library(optparse)
})
stringsAsFactors=FALSE

###############################################  Comparing TERT and TERC expression with EXTEND Scores (Figure 2)  ################################################################################################




root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

A = file.path(root_dir, "data", "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds")   ### Variable representing Expression data for PolyA FPKM

B = file.path(root_dir, "data", "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")  ### Variable representing Expression data for Stranded FPKM

C = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_FPKM.txt")      ### Variable representing Telomerase activity-prediction for PolyA_FPKM data

D = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data



# set up the command line options
option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Plot output (.pdf)")
)


# parse the command line arguments - this becomes a list
# where each element is named by the long flag option
opt <- parse_args(OptionParser(option_list = option_list))

PBTA_EXTEND_TERTCompPlot <- opt$output

PTBA_GE_PolyA = readRDS(A)  ### reading PolyA FPKM data
TelComp = c('TERT','TERC')
PTBA_GE_PolyA_GE = PTBA_GE_PolyA[which(rownames(PTBA_GE_PolyA) %in% TelComp),]    ### Extracting TERT and TERC genes expression from the PolyA FPKM data
PTBA_GE_PolyA_GE = t(PTBA_GE_PolyA_GE)
SampleID = rownames(PTBA_GE_PolyA_GE)
PTBA_GE_PolyA_GE = data.frame(SampleID,PTBA_GE_PolyA_GE)
row.names(PTBA_GE_PolyA_GE)=NULL
PTBA_GE_PolyA_TMScores2 = read.table(C,sep='\t',head=T)   ### Reading Telomerase scores for PolyA FPKM data from results directory
PTBA_GE_TM_PolyA = merge(PTBA_GE_PolyA_TMScores2,PTBA_GE_PolyA_GE,by='SampleID')     ### Merging TERT and TERC expressions with Telomerase scores




PTBA_GE_Standard = readRDS(B)  ### reading Stranded FPKM data
PTBA_GE_Standard_GE = PTBA_GE_Standard[which(rownames(PTBA_GE_Standard) %in% TelComp),] ### Extracting TERT and TERC genes expression from the PolyA FPKM data
PTBA_GE_Standard_GE = t(PTBA_GE_Standard_GE)
SampleID = rownames(PTBA_GE_Standard_GE)
PTBA_GE_Standard_GE = data.frame(SampleID,PTBA_GE_Standard_GE)
row.names(PTBA_GE_Standard_GE)=NULL
PTBA_GE_Stranded_TMScores2 = read.table(D,sep='\t',head=T) ### Reading Telomerase scores for Stranded FPKM data from results directory
PTBA_GE_TM_Stranded = merge(PTBA_GE_Stranded_TMScores2,PTBA_GE_Standard_GE,by='SampleID') ### Merging TERT and TERC expressions with Telomerase scores



pdf(PBTA_EXTEND_TERTCompPlot)


P1 = ggscatter(PTBA_GE_TM_PolyA, x = "NormEXTENDScores", y = "TERT",color = "red",size = 2,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)


P2 = ggscatter(PTBA_GE_TM_PolyA, x = "NormEXTENDScores", y = "TERC",color = "red",size = 1.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)

P3 = ggscatter(PTBA_GE_TM_Stranded, x = "NormEXTENDScores", y = "TERT",color = "red",size = 1.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)


P4 = ggscatter(PTBA_GE_TM_Stranded, x = "NormEXTENDScores", y = "TERC",color = "red",size = 1.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores PolyA_FPKM",ylab="TERT Expression"),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores PolyA_FPKM",ylab="TERC Expression"),vp = define_region(row = 1, col = 2))
print(ggpar(P3,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores Stranded_FPKM",ylab="TERT Expression"),vp = define_region(row = 2, col = 1))
print(ggpar(P4,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores Stranded_FPKM",ylab="TERC Expression"),vp = define_region(row = 2, col = 2))

dev.off()

