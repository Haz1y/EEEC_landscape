#--------------------------------------------------------------
# filename : Oncoprint.R
# Date : 2024-03-02
# contributor : Zhe Hu, MD
# function: Proteogenomic analysis of early-onset endometrioid endometrial carcinoma
# R version: R/4.1.2
#--------------------------------------------------------------

library(maftools)
library(tidyverse)
library(ggsci)
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggsignif)
library(discover)

ctnnb1_input <- read.table("onco_matrix.txt",sep="\t",header=T)

ctnnb1_input[ctnnb1_input==0] <- ""
mygene_new <- as.matrix(ctnnb1_input)
dim(mygene_new) <- c(ncol(ctnnb1_input)*nrow(ctnnb1_input),1)
unique(mygene_new)
ctnnb1_input[ctnnb1_input=="Missense_Mutation"]= "SNV"
ctnnb1_input[ctnnb1_input=="Multi_Hit"]= "SNV"
ctnnb1_input[ctnnb1_input=="Frame_Shift_Ins"]= "InDel"
ctnnb1_input[ctnnb1_input=="In_Frame_Del"]= "InDel"
ctnnb1_input[ctnnb1_input=="Frame_Shift_Del"]= "InDel"


col = c("InDel" = "#3973B7",  "SNV" = "#AE3033")
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w, h*0.66, 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    InDel = function(x, y, w, h) {
        grid.rect(x, y, w, h*0.66, 
            gp = gpar(fill = col["InDel"], col = NA))
    },
    # small green
    SNV = function(x, y, w, h) {
        grid.rect(x, y, w,h*0.66, 
            gp = gpar(fill = col["SNV"], col = NA))
    }
)

pdf("CTNNB1_mutate.pdf",width = 14,height = 6)
oncoPrint(ctnnb1_input %>% select(-FD_BY_025_C,-FD_BY_053_C),
    alter_fun = alter_fun, col = col,
         top_annotation = NULL,show_column_names=T)
dev.off()


ctnnb1_input = apply(ctnnb1_input,2,function(x) as.numeric(as.character(x)))
ctnnb1_input = as.data.frame(ctnnb1_input)
rownames(ctnnb1_input) = c("CTNNB1","KRAS","SOX17")                    

ctnnb1_events = discover.matrix(ctnnb1_input)
groupwise.discover.test(ctnnb1_events,method = "exclusivity")


UCEC_IA_WES_laml = read.maf(maf = UCEC_IA_WES_maf,clinicalData = UCEC_IA_WES_clin,isTCGA = F)
oncoplot(UCEC_IA_WES_laml,writeMatrix=T,genes=c(mutsigcv),
         removeNonMutated = F,
         clinicalFeatures = "age_type",
         sortByAnnotation = T)

CPTAC_maf =read_tsv("./data/ucec_cptac_2020/data_mutations.txt",skip=0)
CPTAC_laml = read.maf(maf = CPTAC_maf)
CPTAC_WES_titv = titv(CPTAC_laml)

TCGA_maf =read_tsv("./data/ucec_tcga_pan_can_atlas_2018/data_mutations.txt",skip=0)
TCGA_laml = read.maf(maf = TCGA_maf)
TCGA_WES_titv = titv(TCGA_laml)


TCGA_titv = TCGA_WES_titv$fraction.contribution
TCGA_titv$Cohort = "TCGA"
CPTAC_titv = CPTAC_WES_titv$fraction.contribution
CPTAC_titv$Cohort = "CPTAC"

TJFD_titv = UCEC_IA_WES_titv$fraction.contribution
TJFD_titv$Cohort = "TJFD"

final_titv = Reduce(rbind,list(TCGA_titv,CPTAC_titv,TJFD_titv))

final_titv$Cohort = factor(final_titv$Cohort ,
                          levels = c("TCGA","CPTAC","TJFD"))

ggplot(final_titv, aes(x=Cohort,y=`C>T`))+
  geom_violin(aes(fill = Cohort,color=Cohort), trim = FALSE)+
  geom_boxplot(aes(fill=Cohort),notch = F,width=0.3)+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
    ggsci::scale_fill_npg()+
    ggsci::scale_color_npg()+
  #annotate("rect", xmin = 0, xmax =3.5,  ymin =-5 , ymax = 105, alpha = 0.1,fill="#FAE3AD") +
  theme(axis.text=element_text(colour='black',size=10))+
  scale_y_continuous(expand = c(0, 0), limit = c(-5, 110))+
  geom_signif(
    comparisons = list(c("TCGA", "CPTAC")),
    map_signif_level = FALSE,
    y_position = 90, tip_length = 0, vjust = 0.2)+
  geom_signif(
    comparisons = list(c("CPTAC", "TJFD")),
    map_signif_level = FALSE,
    y_position = 95, tip_length = 0, vjust = 0.2)+
  geom_signif(
    comparisons = list(c("TCGA", "TJFD")),
    map_signif_level = FALSE,
    y_position = 100, tip_length = 0, vjust = 0.2
  )


load("../../02.UCEC/07.between_cohort/data//3.5.Compare_Endo_Mut&Clin.RDa")
cgc_list = read_tsv("../../../../Project/database/COSMIC/Census_allTue Apr 26 13_43_25 2022.tsv")
UCEC_clin = read.xlsx("../../02.UCEC_figure/0.raw_data/sample_upload_220615_3.xlsx",sheet=4)
UCEC_clin$Tumor_Sample_Barcode = paste0(UCEC_clin$patient,"_C")
total_maf = read_tsv("../../../../DNAseq/UCEC_un_annovar/10_maf_merge/merge.total.unannovar.maf")
##使用SMG q<0.05
mutsigcv <- read.table("../../02.UCEC/07.between_cohort/mutsig_results/total/output/TJFD_onlyC_no_AEH.sig_genes.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F)
mutsigcv <- rownames(mutsigcv[mutsigcv$q < 0.01,])
AEH_C_maf =total_maf %>% filter(!endsWith(Tumor_Sample_Barcode,"_LN"))
length(unique(AEH_C_maf$Tumor_Sample_Barcode))


UCEC_clin$age_type = ifelse(UCEC_clin$age_at_diagnosis <=40, "Early","Old")
UCEC_IA_WES_clin = UCEC_clin %>% filter(clinical_stage == "IA") %>% filter(Tumor_Sample_Barcode %in% AEH_C_maf$Tumor_Sample_Barcode)
UCEC_IA_WES_maf = total_maf %>% filter(Tumor_Sample_Barcode %in% UCEC_IA_WES_clin$Tumor_Sample_Barcode)
IA_Early_clin = UCEC_IA_WES_clin %>% filter(age_at_diagnosis<=40)
IA_Old_clin = UCEC_IA_WES_clin %>% filter(age_at_diagnosis>40)
IA_Early_maf = UCEC_IA_WES_maf %>% filter(Tumor_Sample_Barcode %in% IA_Early_clin$Tumor_Sample_Barcode)
IA_Old_maf = UCEC_IA_WES_maf %>% filter(Tumor_Sample_Barcode %in% IA_Old_clin$Tumor_Sample_Barcode)


UCEC_IA_WES_laml = read.maf(maf = UCEC_IA_WES_maf,clinicalData = UCEC_IA_WES_clin,isTCGA = F)
oncoplot(UCEC_IA_WES_laml,writeMatrix=T,genes=c(mutsigcv),
         removeNonMutated = F,
         clinicalFeatures = "age_type",
         sortByAnnotation = T)


tmb = tmb(UCEC_IA_WES_laml)


input <- read.table("onco_matrix.txt",sep="\t",header=T)
input[input==0] <- ""
mygene_new <- as.matrix(input)
dim(mygene_new) <- c(ncol(input)*nrow(input),1)
unique(mygene_new)
my_age_type <- function(X){
    ifelse(X>61,"(61-74]",
           ifelse(X > 40,"(40-61]","(19-40]"))
}
UCEC_clin$age_three_subtype = my_age_type(UCEC_clin$age_at_diagnosis)


t_cell <- read.table("../TcellExTRECT//package//t_cell_final.csv",sep=",",header=T)
t_cell <- t_cell[!duplicated(t_cell$sample),]
t_cell <- as.data.frame(t_cell)
tumor_t_cell <- t_cell[ t_cell$sample %in% UCEC_clin$Tumor_Sample_Barcode,]
rownames(tumor_t_cell) <- tumor_t_cell$sample
tumor_t_cell <- tumor_t_cell[ colnames(input),]


## purity
load("../../02.UCEC/purity//final_purity.RDa")
## fga wgd
load("../../02.UCEC/Fig2/data/common_Clin.RDa")
final_purity$Sample = rownames(final_purity)
common_fga$Sample = rownames(common_fga)
common_clin = UCEC_clin %>% 
    filter(Tumor_Sample_Barcode %in% colnames(input))
rownames(common_clin) = common_clin$Tumor_Sample_Barcode
common_clin = common_clin[colnames(input),]


common_clin = merge(common_clin,tumor_t_cell[,c("sample","TCRA.tcell.fraction")],
     by.x="Tumor_Sample_Barcode",by.y="sample",all.x=T)
common_clin = merge(common_clin,final_purity[,c("Sample","median_normalize_purity")],
     by.x="Tumor_Sample_Barcode",by.y="Sample",all.x=T)
common_clin = merge(common_clin,common_fga[,c("Sample","FGA","GAIN","LOSS","IS_WGD")],
     by.x="Tumor_Sample_Barcode",by.y="Sample",all.x=T)
common_clin = merge(common_clin,tmb[,c("Tumor_Sample_Barcode","total")],
     by.x="Tumor_Sample_Barcode",by.y="Tumor_Sample_Barcode",all.x=T)


common_clin = as.data.frame(common_clin)
rownames(common_clin) = common_clin$Tumor_Sample_Barcode
common_clin = common_clin[colnames(input),]


mycol <- colorRampPalette(brewer.pal(11,'Spectral'))(12)[12:1]
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#dcddde", col = NA))
  },  
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "black", col = NA)) 
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#219f94", col = NA)) 
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#1d6ed0", col = NA)) 
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#2a520c", col = NA)) 
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#ff6565", col = NA)) 
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#e4ed2e", col = NA)) 
  },  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#ae441f", col = NA)) 
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#d39655", col = NA)) 
  }
)

col = c("Missense_Mutation" = "#219f94", 
        "Nonsense_Mutation" = "#1d6ed0", 
        "Frame_Shift_Ins" = "#2a520c", 
        "Frame_Shift_Del" = "#ff6565", 
        "In_Frame_Ins" = "#e4ed2e", 
        "In_Frame_Del" = "#ae441f", 
        "Splice_Site" ="#d39655",
        "Multi_Hit" = "black")


##绘制BY的SMG
pdf("maf_oncoprint_sort_tmb_width.pdf",width=10,height=11.5)
# 绘图
m1 <- oncoPrint(input, get_type = function(x) x,
           alter_fun = alter_fun, col = col,
          show_pct = T, #show pct in left
          column_title = "",
          #bottom_annotation = my_annotation,
          show_heatmap_legend=T,
          column_title_gp = gpar(fontsize = 6),
          row_names_gp = gpar(fontsize = 6),
          row_names_side = "left",
          pct_side = "right",
          column_names_gp = gpar(fontsize = 6),
           column_split = factor(common_clin$age_type,levels=c("Early","Old")),
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = list(direction = "horizontal",nrow=2, title = ""),
          column_order=rownames(common_clin)[order(common_clin$total,decreasing=T)],
           top_annotation = HeatmapAnnotation(Nonsyn_Count = anno_barplot(common_clin$total,
                                                                          border = F,gp = gpar(fill = "#1a78b4",col="#1a78b4"),
                                                                          annotation_name="Nonsyn Mutations",bar_width = 0.7,
                                                                          height = unit(2.5, "cm")),
                                             TCGA_Subtype=common_clin$Subtype,
                                              Figo=common_clin$clinical_stage, 
                                             Grade = common_clin$clinical_grade,
                                             Age = common_clin$age_three_subtype, 
                                            T_cell_frac = common_clin$TCRA.tcell.fraction,
                                              consensus_normal_purity = common_clin$median_normalize_purity, 
                                              fga = common_clin$FGA,
                                              gain = common_clin$GAIN,
                                              loss = common_clin$LOSS,
                                              wgd = common_clin$IS_WGD,
                                              annotation_label = c("Nonsyn Mutations","TCGA Subtype","Figo Stage","Figo Grade","Age", "T cell fraction",
                                                                   "Consensus Purity","FGA","Gain","Loss","WGD"),
                                             col = list(TCGA_Subtype = c("AEH"="#01a05c","CN high" = "#d82249","CN low" = "#e5e3c9","MSI_hypermutated" = "#9adcff","POLE_ultramutated" ="#925e9f"),
                                                        Grade = c("G1"="#4fd3c4","G2"="#70b2d2","G3"="#d82148"),
                                                        Figo = c("AEH" = "#0072b5","IA" = "#7876b1","IB-II" = "#ffdc91","III-IV" = "#ee4c97"),
                                                        Age = c("(19-40]" = "#b0e2fe","(40-61]" = "#33cd2e","(61-74]" = "#1c8c1c"),
                                                       wgd = c("FALSE"="#bebebe","TRUE"="#ffdc91")),
                                            show_annotation_name = TRUE,
                                            annotation_name_side = "left",
                                            annotation_name_gp = gpar(fontsize = 6)),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = T)))
draw(m1,merge_legends = FALSE,heatmap_legend_side = "bottom")
dev.off()

common_clin %>% 
ggplot( aes(x=age_type , y=TMB))+
stat_boxplot(geom ="errorbar", width = .3, size=.7,  coef = 1) +
  geom_boxplot(width = .55, outlier.shape = NA,  coef = 1) +
  geom_point(position = position_jitter(seed = 2021, width = .3),
             aes(color = age_type), alpha = .7, size = 6) +
  scale_color_manual(values =  c("#4CA7BE", "#D94B35"))+
 #scale_x_discrete(position = "top") +
#facet_wrap(~ TCGA_Subtype,  nrow = 1, strip.position = "right") 
  theme_bw()+
  theme(strip.text = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 1,color="black"),
        axis.ticks.x.top = element_line(color = "black", size = .5),
        axis.ticks.length = unit(-0.2, "cm"),
        axis.ticks.y = element_line(color = "black", size = .5),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin=unit(c(0.5,0.1,0.5,0.1),"cm")) +
       #geom_rect(xmin=1, xmax=2, ymin=-1.7, ymax=-2.3, color="black", fill="white", size = .5)+
       geom_signif(comparisons = list(c("Early","Old")),
              map_signif_level=T,vjust=0.5,color="black",
              textsize=5,test=t.test,step_increase=0.1)+
        labs(title = "TMB of different Age Group in Total")




ggplot(common_clin,aes(age_type,BMI))+
  ggdist::stat_halfeye(mapping = aes(fill=age_type),width = .34, .width = 0, justification = -.4, point_colour = NA,alpha=0.9) + 
  geom_jitter(mapping = aes(color=age_type),width = .12, alpha = .7,size=4.5)+
  geom_boxplot(width = .17, outlier.shape = NA,fill=NA,size=0.5)+
  labs(title = "BMI of different age group")+
  scale_y_continuous(expand = c(0.02,0))+
  scale_x_discrete(expand = c(0,0.3))+
  scale_fill_manual(values = c("Early"="#4DA8BE","Old"="#DA4C35"))+
  scale_color_manual(values = c("Early"="#4DA8BE","Old"="#DA4C35"))+

scale_y_continuous(breaks=c(20,25,30,35,40))+
 geom_signif(comparisons = list(c("Early","Old")),
              map_signif_level=T,vjust=0.5,color="black",
              textsize=5,test=wilcox.test,step_increase=0.1)+
  theme_classic()+
  theme(
    axis.text.x.bottom = element_text(angle = 0,hjust = 0.5,size = 16,color = "black"),
    axis.text.y.left = element_text(size = 16,color = "black"),
    axis.title.y.left = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.ticks.length=unit(.2, "cm"),
    plot.title = element_text(size = 18,hjust = 0.5),
    legend.position = "none",
  )