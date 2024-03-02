#--------------------------------------------------------------
# filename : 3D-Sig_Explorer.R
# Date : 2024-03-02
# contributor : Zhe Hu, MD
# function: Proteogenomic analysis of early-onset endometrioid endometrial carcinoma
# R version: R/4.1.2
#--------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsci)
library(PResiduals)

load("mt_sig.Rda")
sig_count = mt_sig$Exposure
sig_count = t(sig_count)
sig_count = as.data.frame(sig_count)
sig_count$Tumor_Sample_Barcode = rownames(sig_count)
UCEC_IA_clin = UCEC_clin %>% filter(clinical_stage== "IA")
common_sig_IA_count = sig_count[ UCEC_IA_clin$Tumor_Sample_Barcode,]
common_sig_IA_count_clin = merge(common_sig_IA_count,UCEC_IA_clin,by="Tumor_Sample_Barcode")
common_sig_IA_pct = sig_pct[ UCEC_IA_clin$Tumor_Sample_Barcode,]
common_sig_IA_pct_clin = merge(common_sig_IA_pct,UCEC_IA_clin,by="Tumor_Sample_Barcode")
Early_common_sig_IA_pct_clin = common_sig_IA_pct_clin %>%filter(age_at_diagnosis<40)
SBS_BY_compare_plot = common_sig_IA_pct_clin[,c("Tumor_Sample_Barcode","age_at_diagnosis",paste0("Sig",1:14))]
SBS_BY_compare_plot$age_group = ifelse(SBS_BY_compare_plot$age_at_diagnosis<=40,"Early","Old")
rownames(SBS_BY_compare_plot) = SBS_BY_compare_plot$Tumor_Sample_Barcode

SBS_BY_compare_plot_long$sig_num = as.double(gsub("Sig","",SBS_BY_compare_plot_long$signature))/2
SBS_BY_compare_plot_long %>% 
ggplot( aes(x=sig_num , y=Exposure.norm,fill=age_group))+
    #geom_stripped_cols(odd='white',even='grey90')+

    stat_summary(fun.min = function(x){quantile(x)[2]},fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='#919191',width=0.01,size=0.5,
               position = position_dodge(width = 0.3))+
    geom_point(size=3,stat = 'summary',fun=median,position = position_dodge(width = 0.3),aes(color = age_group))+
    theme_bw()+
    theme(legend.title = element_text(face="bold"),
      legend.position = c(.85, .65),
      plot.title = element_text(hjust=0.5,face="bold", size=15),
        strip.text = element_blank(),
        panel.grid = element_blank())+
    stat_compare_means(
      hide.ns = F,
      method = "t.test",label ='p.signif'
     )+
    scale_x_continuous(limits=c(0.4,7.1),breaks=c(0.5, 1, 1.5, 2,
                                                2.5, 3, 3.5, 4, 
                                                4.5, 5, 5.5, 6, 
                                                6.5, 7),
                      labels = paste0("Sig",1:14))+
    scale_color_npg()


ggscatter(common_sig_IA_pct_clin, x = "age_at_diagnosis", y = "Sig6", 
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#3fa3ae", fill = "#deeded"),
          shape = 7,size= 1,color ="#585958",
          cor.coef = TRUE, cor.method = "spearman")


ggscatter(common_sig_IA_count_clin, x = "age_at_diagnosis", y = "Sig4", 
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#3fa3ae", fill = "#deeded"),
          shape = 7,size= 1,color ="#585958",
          cor.coef = TRUE, cor.method = "spearman")


common_Clin_IA = common_Clin %>% filter(stage_subtype == "IA")
partial_Spearman(TMB_Nonsyno|age_at_diagnosis ~ TCGA_Subtype,data=common_Clin_IA)

common_Clin_IA = common_Clin %>% filter(stage_subtype == "IA")
partial_Spearman(TMB_Nonsyno|age_at_diagnosis ~ stage_subtype,
                  fit.x="poisson", fit.y="lm.emp",data=common_Clin)

partial_Spearman(TMB_Nonsyno|age_at_diagnosis ~ stage_subtype,
                  data=common_Clin)



SBS_BY_compare_plot_long %>% 
ggplot( aes(x=sig_num , y=Exposure.norm,fill=age_group))+
    #geom_stripped_cols(odd='white',even='grey90')+

    stat_summary(fun.min = function(x){quantile(x)[2]},fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='#919191',width=0.01,size=0.5,
               position = position_dodge(width = 0.3))+
    geom_point(size=3,stat = 'summary',fun=median,position = position_dodge(width = 0.3),aes(color = age_group))+
    theme_bw()



sig_pct = mt_sig$Exposure.norm
IA_excel_clin = UCEC_clin %>% filter(Tumor_Sample_Barcode %in% colnames(sig_pct)) %>% 
                            filter(clinical_stage == "IA")
IA_excel_clin$age = ifelse(IA_excel_clin$age_at_diagnosis<=40,"Early","Old")
IA_excel_clin$obesity = ifelse(IA_excel_clin$BMI>=28,"Obesity","Normal")
IA_age_order = IA_excel_clin %>% arrange(age_at_diagnosis) 

IA_maf = total_maf %>% filter(Tumor_Sample_Barcode %in% IA_age_order$Tumor_Sample_Barcode)
ctnnb1_maf =  IA_maf %>% filter(Hugo_Symbol == "CTNNB1") %>% filter(Variant_Classification == "Missense_Mutation")
siglec10_maf =  IA_maf %>% filter(Hugo_Symbol == "SIGLEC10") %>% filter(Variant_Classification == "Missense_Mutation")

IA_age_order$ctnnb1_mut = ifelse(IA_age_order$Tumor_Sample_Barcode  %in% ctnnb1_maf$Tumor_Sample_Barcode,
                                "Mut","WT")
IA_age_order$siglec10_mut = ifelse(IA_age_order$Tumor_Sample_Barcode  %in% siglec10_maf$Tumor_Sample_Barcode,
                                "Mut","WT")
IA_age_order$ctnnb1_mut_binomial = ifelse(IA_age_order$Tumor_Sample_Barcode  %in% ctnnb1_maf$Tumor_Sample_Barcode,
                                1,0)
IA_age_order$siglec10_mut_binomial = ifelse(IA_age_order$Tumor_Sample_Barcode  %in% siglec10_maf$Tumor_Sample_Barcode,
                                1,0)
sig_pct_input = t(sig_pct)
sig_pct_input = as.data.frame(sig_pct_input)
sig_pct_input$Sample = rownames(sig_pct_input)

IA_age_order_sig = cbind(IA_age_order,sig_pct_input[ IA_age_order$Tumor_Sample_Barcode,])

my_sig_name = function(X){names(which.max(X[paste0("Sig",1:14)]))}
my_sig_value = function(X){max(X[paste0("Sig",1:14)])}

for (i in 1:nrow(IA_age_order_sig)){
    IA_age_order_sig$domin_sig[i] = my_sig_name(IA_age_order_sig[i,])
    IA_age_order_sig$domin_sig_pct[i] = my_sig_value(IA_age_order_sig[i,])

    }


IA_maf_sig_age = merge(IA_maf,IA_age_order_sig[,c("Tumor_Sample_Barcode","age","age_at_diagnosis","domin_sig","domin_sig_pct")],
      by="Tumor_Sample_Barcode")
IA_maf_sig_age_non = IA_maf_sig_age %>% filter(Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation",
                                                                            "Frame_Shift_Del","Frame_Shift_Ins",
                                                                            "In_Frame_Del","In_Frame_Ins","Splice_Site"))
IA_early_maf_sig_age_non = IA_maf_sig_age_non %>% filter(age =="Early")
IA_early_maf_sig_age_non$sample_count =69
IA_old_maf_sig_age_non = IA_maf_sig_age_non %>% filter(age =="Old")
IA_old_maf_sig_age_non$sample_count = 39
IA_total_maf_sig_age_non = rbind(IA_early_maf_sig_age_non,IA_old_maf_sig_age_non)

my_short_table = IA_total_maf_sig_age_non %>% filter(Hugo_Symbol == "CTNNB1") %>% 
select(HGVSp_Short,domin_sig,age,Exon_Number) %>% group_by(Exon_Number,HGVSp_Short) %>% summarise(n()) %>% arrange(-`n()`)
hotspot = as.data.frame(my_short_table)[1:9,2]


my_plot = IA_total_maf_sig_age_non %>% filter(Hugo_Symbol == "CTNNB1") %>% 
select(HGVSp_Short,domin_sig,age,Exon_Number)%>% 
filter(age =="Early") 
my_plot[!my_plot$HGVSp_Short %in% hotspot,"HGVSp_Short"] = "Others"
my_plot$HGVSp_Short = factor(my_plot$HGVSp_Short,levels =c(hotspot,"Others"))
my_plot = my_plot %>% arrange(HGVSp_Short)


FigS5A = ggplot(my_plot)+
  geom_col(aes(x = 1,fill = Exon_Number,y=1))+
ggsci::scale_fill_npg()+
  new_scale('fill') +
  geom_col(width = 0.1,aes(x = 1.55,fill =HGVSp_Short,y=1))+
  scale_fill_brewer(palette = 'Paired') +
  coord_polar(theta = 'y') +
theme_void()+ scale_x_continuous(expand = expansion(add = c(1.2,0)))+scale_y_reverse()
my_short_table = IA_total_maf_sig_age_non %>% filter(Hugo_Symbol == "CTNNB1") %>% 
select(HGVSp_Short,domin_sig,age,Exon_Number) %>% group_by(Exon_Number,HGVSp_Short) %>% summarise(n()) %>% arrange(-`n()`)
hotspot = as.data.frame(my_short_table)[1:9,2]

my_plot = IA_total_maf_sig_age_non %>% filter(Hugo_Symbol == "CTNNB1") %>% 
select(HGVSp_Short,domin_sig,age,Exon_Number)%>% 
filter(age =="Old") 
my_plot[!my_plot$HGVSp_Short %in% hotspot,"HGVSp_Short"] = "Others"
my_plot$HGVSp_Short = factor(my_plot$HGVSp_Short,levels =c(hotspot,"Others"))
my_plot = my_plot %>% arrange(HGVSp_Short)
FigS5B = ggplot(my_plot)+
  geom_col(aes(x = 1,fill = Exon_Number,y=1))+
ggsci::scale_fill_npg()+
  new_scale('fill') +
  geom_col(width = 0.1,aes(x = 1.55,fill =HGVSp_Short,y=1))+
  scale_fill_manual(values = c("#1F78B4","#6A3D9A")) +
  coord_polar(theta = 'y') +
theme_void()+ scale_x_continuous(expand = expansion(add = c(1.2,0)))+scale_y_reverse()



denovo_sig_relative_matrix = mt_sig$Exposure.norm
denovo_sig_pattern_matrix = mt_sig$Signature.norm

my_pms = function(X){
    sample_exposure             = denovo_sig_relative_matrix[,X[1,1]]
    signature_motif_pct         = as.double(denovo_sig_pattern_matrix[X[1,2],])
    sample_signature_motif_pct  = sample_exposure * signature_motif_pct
    sample_pms                  = sample_signature_motif_pct/ sum(sample_signature_motif_pct)
    X$dominate_signature        = names(which.max(sample_pms))
    X$dominate_pms              = max(sample_pms)
    X = cbind(X,as.data.frame(t(sample_pms)))
    return(X)
    }

CTNNB1_IA_sig_age_maf = IA_maf_sig_age_non %>% filter(Hugo_Symbol == "CTNNB1")
CTNNB1_laml <- read_maf(maf = CTNNB1_IA_sig_age_maf)
CTNNB1_mats <- CTNNB1_mt_tally <- sig_tally(
  CTNNB1_laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = FALSE,
  mode = "ALL"
)
CTNNB1_SBS96 = CTNNB1_mt_tally$SBS_96
CTNNB1_SBS96 = as.data.frame(CTNNB1_SBS96)
CTNNB1_SBS96$Sample = rownames(CTNNB1_SBS96)
CTNNB1_SBS96_long_m<-melt(CTNNB1_SBS96, id.vars = c("Sample"), 
                  measure.vars = colnames(CTNNB1_mt_tally$SBS_96),
                  variable.name = c('Motif'),
                  value.name = 'Count')
CTNNB1_SBS96_long_m_nozero = CTNNB1_SBS96_long_m %>% filter(Count > 0)

CTNNB1_SBS96_pms = data.frame()

for (i in 1:nrow(CTNNB1_SBS96_long_m_nozero)){
    CTNNB1_SBS96_pms = rbind(CTNNB1_SBS96_pms,
                             my_pms(CTNNB1_SBS96_long_m_nozero[i,]))

    }
CTNNB1_SBS96_pms_age = merge(CTNNB1_SBS96_pms,IA_age_order[,c("Tumor_Sample_Barcode","age_at_diagnosis")]
                             ,by.x="Sample",by.y="Tumor_Sample_Barcode")
CTNNB1_SBS96_pms_age = CTNNB1_SBS96_pms_age[,c(1,19,2:18)]
CTNNB1_SBS96_pms_age_mut = merge(CTNNB1_SBS96_pms_age,ctnnb1_maf[,c("Tumor_Sample_Barcode","HGVSc","HGVSp_Short")],
      by.x="Sample",by.y="Tumor_Sample_Barcode")

CTNNB1_SBS96_pms_age_mut_count = CTNNB1_SBS96_pms_age_mut %>% group_by(HGVSp_Short) %>% summarise(n()) %>% arrange(-`n()`)
CTNNB1_SBS96_pms_age_mut_count_2 = CTNNB1_SBS96_pms_age_mut_count %>% filter(`n()` <=2)
CTNNB1_SBS96_pms_age_mut$HGVSp_Short_merge = CTNNB1_SBS96_pms_age_mut$HGVSp_Short
CTNNB1_SBS96_pms_age_mut[CTNNB1_SBS96_pms_age_mut$HGVSp_Short_merge %in% CTNNB1_SBS96_pms_age_mut_count_2$HGVSp_Short,"HGVSp_Short_merge"]="Others"
plot_hotspot = c("p.D32Y","p.G34V","p.G34R","p.S33C","p.S37C","p.T41I","p.D32G","p.S37F")
CTNNB1_SBS96_pms_age_mut$HGVSp_Short_plot = CTNNB1_SBS96_pms_age_mut$HGVSp_Short
CTNNB1_SBS96_pms_age_mut[!(CTNNB1_SBS96_pms_age_mut$HGVSp_Short %in% plot_hotspot),"HGVSp_Short_plot"]="Others"
CTNNB1_SBS96_pms_age_mut$HGVSp_Short_plot = factor(CTNNB1_SBS96_pms_age_mut$HGVSp_Short_plot,
                                                  levels=c("p.D32Y","p.G34V",
                                                          "p.G34R","p.S33C",
                                                           "p.S37C","p.T41I",
                                                           "p.D32G","p.S37F","Others"))


SIGLEC10_maf = IA_maf %>% filter(Hugo_Symbol == "SIGLEC10") %>%filter(Variant_Classification %in% c( "Missense_Mutation","Nonsense_Mutation"))
SIGLEC10_laml <- read_maf(maf = SIGLEC10_maf)
SIGLEC10_mats <- SIGLEC10_mt_tally <- sig_tally(
 SIGLEC10_laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = TRUE,
  mode = "ALL"
)

SIGELC10_SBS96 = SIGLEC10_mt_tally$SBS_96
SIGELC10_SBS96 = as.data.frame(SIGELC10_SBS96)
SIGELC10_SBS96$Sample = rownames(SIGELC10_SBS96)
SIGELC10_SBS96_long_m <- melt(SIGELC10_SBS96, id.vars = c("Sample"), 
                  measure.vars = colnames(SIGLEC10_mt_tally$SBS_96),
                  variable.name = c('Motif'),
                  value.name = 'Count')
SIGELC10_SBS96_long_m_nozero = SIGELC10_SBS96_long_m %>% filter(Count > 0)
SIGELC10_SBS96_pms = data.frame()

for (i in 1:nrow(SIGELC10_SBS96_long_m_nozero)){
    SIGELC10_SBS96_pms = rbind(SIGELC10_SBS96_pms,
                             my_pms(SIGELC10_SBS96_long_m_nozero[i,]))

    }

SIGELC10_SBS96_pms_age = merge(SIGELC10_SBS96_pms,IA_age_order[,c("Tumor_Sample_Barcode","age_at_diagnosis")],
                          by.x="Sample",by.y="Tumor_Sample_Barcode")
SIGELC10_SBS96_pms_age = SIGELC10_SBS96_pms_age[,c(1,19,2:18)]
SIGELC10_SBS96_pms_age_mut = merge(SIGELC10_SBS96_pms_age,SIGLEC10_maf[,c("Tumor_Sample_Barcode","HGVSc","HGVSp_Short")],
      by.x="Sample",by.y="Tumor_Sample_Barcode")


APC_maf = IA_maf %>% filter(Hugo_Symbol == "APC") %>%filter(Variant_Classification %in% c( "Missense_Mutation","Nonsense_Mutation"))
APC_laml <- read_maf(maf = APC_maf)
APC_mats <- APC_mt_tally <- sig_tally(
 APC_laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = TRUE,
  mode = "ALL"
)
APC_SBS96 = APC_mt_tally$SBS_96
APC_SBS96 = as.data.frame(APC_SBS96)
APC_SBS96$Sample = rownames(APC_SBS96)
APC_SBS96_long_m<-melt(APC_SBS96, id.vars = c("Sample"), 
                  measure.vars = colnames(APC_mt_tally$SBS_96),
                  variable.name = c('Motif'),
                  value.name = 'Count')
APC_SBS96_long_m_nozero = APC_SBS96_long_m %>% filter(Count > 0)
APC_SBS96_pms = data.frame()

for (i in 1:nrow(APC_SBS96_long_m_nozero)){
    APC_SBS96_pms = rbind(APC_SBS96_pms,
                             my_pms(APC_SBS96_long_m_nozero[i,]))

    }
APC_SBS96_pms_age = merge(APC_SBS96_pms,IA_age_order[,c("Tumor_Sample_Barcode","age_at_diagnosis")],
                          by.x="Sample",by.y="Tumor_Sample_Barcode")
APC_SBS96_pms_age = APC_SBS96_pms_age[,c(1,19,2:18)]
APC_SBS96_pms_age_mut = merge(APC_SBS96_pms_age,APC_maf[,c("Tumor_Sample_Barcode","HGVSc","HGVSp_Short")],
      by.x="Sample",by.y="Tumor_Sample_Barcode")

my_comparisions = list(c("CTNNB1","APC"),
                      c("CTNNB1","SIGLEC10"))
Figure_my = ggboxplot(
  total_gene_plot, x = "Hugo_Symbol", y = "Sig6",
   color = "Hugo_Symbol",
add = "jitter" )+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
scale_color_npg()+
geom_signif(comparisons = my_comparisions,
              test = "wilcox.test",
                map_signif_level = T)
