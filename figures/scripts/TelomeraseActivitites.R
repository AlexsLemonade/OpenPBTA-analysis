library('ggpubr')
stringsAsFactors=FALSE
library(stringr)
library(gridBase)
library(gridGraphics)
library(optparse)
library(forcats) ### for fct_reorder()

########################################## Figure A #########################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_counts.txt")   ### Variable representing Telomerase activity-prediction for Polya_counts
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_counts.txt")  ### Variable representing Telomerase activity-prediction for Stranded_counts
				 
C = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_FPKM.txt")      ### Variable representing Telomerase activity-prediction for PolyA_FPKM data
				 
D = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data


				 


# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")
		 
				 				 
PBTA_GEComparisonPlot <- file.path(output_dir, "CompareExtend.png")
PBTA_EXTEND_HistologyCompPlot <- file.path(output_dir, "ExtendAcrossHistologies.png")
PBTA_EXTEND_MedulloSubtypes <- file.path(output_dir, "Extend_MedulloSubtypes.png")
telomerase_png <- file.path(output_dir, "Telomerase_Activities.png")


				 
PTBA_GE_PolyA_TMScores1 = read.table(A,sep='\t',head=T)
colnames(PTBA_GE_PolyA_TMScores1)[colnames(PTBA_GE_PolyA_TMScores1)=="NormEXTENDScores"]='NormEXTENDScores_PolyACounts'

PTBA_GE_Standard_TMScores1 = read.table(B,sep='\t',head=T)
colnames(PTBA_GE_Standard_TMScores1)[colnames(PTBA_GE_Standard_TMScores1)=="NormEXTENDScores"]='NormEXTENDScores_StrandedCounts'

PTBA_GE_PolyA_TMScores2 = read.table(C,sep='\t',head=T)
colnames(PTBA_GE_PolyA_TMScores2)[colnames(PTBA_GE_PolyA_TMScores2)=="NormEXTENDScores"]='NormEXTENDScores_PolyA_FPKM'

PTBA_GE_Standard_TMScores2 = read.table(D,sep='\t',head=T)
colnames(PTBA_GE_Standard_TMScores2)[colnames(PTBA_GE_Standard_TMScores2)=="NormEXTENDScores"]='NormEXTENDScores_Stranded_FPKM'

PTBA_GE_PolyA_TMScores = merge(PTBA_GE_PolyA_TMScores1,PTBA_GE_PolyA_TMScores2,by='SampleID')
PTBA_GE_Standard_TMScores = merge(PTBA_GE_Standard_TMScores1,PTBA_GE_Standard_TMScores2,by='SampleID')



png(PBTA_GEComparisonPlot,width = 4, height = 2, units = "in", res = 1600)


P1 = ggscatter(PTBA_GE_PolyA_TMScores, x = "NormEXTENDScores_PolyACounts", y = "NormEXTENDScores_PolyA_FPKM",color = "red",size = 0.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=1)

		   
P2 = ggscatter(PTBA_GE_Standard_TMScores, x = "NormEXTENDScores_StrandedCounts", y = "NormEXTENDScores_Stranded_FPKM",color = "red",size = 0.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=1)
		   
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =4,font.ytickslab =4,font.x = 4,font.y=4,xlab="EXTEND_PolyA_Counts",ylab="EXTEND_PolyA_FPKM"),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =4,font.ytickslab =4,font.x = 4,font.y=4,xlab="EXTEND_Stranded_Counts",ylab="EXTEND_Stranded_FPKM"),vp = define_region(row = 1, col = 2))

		   
dev.off()

########################################## Figure B #########################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "data", "pbta-histologies.tsv")   ### Variable representing clinical data 
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data



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


png(PBTA_EXTEND_HistologyCompPlot,width = 6, height = 8, units = "in", res = 1500)

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


P1 = ggplot(Stranded_Histology, aes(x=fct_reorder(short_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="pink",alpha=0.7)+ geom_jitter(alpha= 0.7,shape=16, cex=0.3)
P2 = ggplot(Stranded_Histology, aes(x=fct_reorder(broad_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="pink",alpha=0.7)+ geom_jitter(alpha=0.7,shape=16, cex=0.3)

##plot_grid(P1,P2, nrow = 2,labels = "AUTO", label_size = 12, scale=c(0.95,0.95),  align = "v")

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 0,font.y=7,ylab="EXTEND Scores",title = "Tumor Short Histology",font.title = 7),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 0,font.y=7,ylab="EXTEND Scores",title = "Tumor Broad Histology",font.title = 7),vp = define_region(row = 2, col = 1))

		
dev.off()


########################################## Figure C #########################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "data", "pbta-histologies.tsv")   ### Variable representing clinical data 
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data


PTBA_Histology = read.table(A,sep='\t',head=T)    ## Reading the clinical data

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


png(PBTA_EXTEND_MedulloSubtypes,width = 2, height = 2, units = "in", res = 1500)

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

############# Multipanel Figure



telomerase_figure <- multi_panel_figure(columns = 4,
                                            rows = 1,
                                            width = 1200,
                                            height = 300,
                                            panel_label_type = "none")

telomerase_figure <- fill_panel(telomerase_figure,
                                    PBTA_GEComparisonPlot,
                                    col = 1:2,
                                    scaling = "fit")

telomerase_figure <- fill_panel(telomerase_figure,
                                    PBTA_EXTEND_HistologyCompPlot,
                                    col = 3,
                                    scaling = "fit")

telomerase_figure <- fill_panel(telomerase_figure,
                                    PBTA_EXTEND_MedulloSubtypes,
                                    col = 4,
                                    scaling = "fit")

save_multi_panel_figure(telomerase_figure, telomerase_png)



#### Removing temporary pngs


file.remove(c(PBTA_GEComparisonPlot,
              PBTA_EXTEND_HistologyCompPlot,
              PBTA_EXTEND_MedulloSubtypes))