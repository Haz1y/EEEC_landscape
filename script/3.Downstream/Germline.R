#--------------------------------------------------------------
# filename : Germline.R
# Date : 2024-03-02
# contributor : Zhe Hu, MD
# function: Proteogenomic analysis of early-onset endometrioid endometrial carcinoma
# R version: R/4.1.2
#--------------------------------------------------------------


library(tidyverse)
load("sample_upload_0615_clin_maf.RDa")
load("../02.UCEC/Fig2/data/common_Clin.RDa")
germ = read_tsv("../02.UCEC/data//germ_final.txt")


#normal_germ_IA_anno  = normal_germ_IA_anno[,c(1:2,19:22,18,3:17)]
MSI_gene = c("MSH2","MSH3","MSH6","MLH1","MLH3","PMS2")
HRD_gene = c("RAD50","BRCA1","BRCA2","NF1","RB1","PALB2","SMAD4")
wnt_gene = c("APC","ATM","MUTYH","NBN","RAD51C","RAD51D","SMAD4","STK11","TP53","CDK4","CDKN2A","CHEK2","BARD1")
late_count  = normal_germ_IA_anno %>% filter(age_at_diagnosis > 40) %>% group_by(Gene.refGene) %>% summarise(n())
late_count = as.data.frame(late_count)
late_count = late_count[order(-late_count[,2]),]
early_count  = normal_germ_IA_anno %>% filter(age_at_diagnosis <= 40) %>% group_by(Gene.refGene) %>% summarise(n())
early_count = as.data.frame(early_count)
early_count = early_count[order(-early_count[,2]),]
colnames(early_count)[2] = "early_count"
colnames(late_count)[2] = "late_count"
total_count = merge(early_count,late_count, by = "Gene.refGene",all= TRUE)


age_germ_plot = ggplot(sig_count_age_IA_order,aes(x = order_patient,
                y=age_at_diagnosis,color = germ))+
  geom_point(aes(y= age_at_diagnosis_na),size=4)+
  geom_segment(aes(x=reorder(patient,age_at_diagnosis),xend=reorder(patient,age_at_diagnosis),
                   y=20,yend=age_at_diagnosis))+
  scale_color_manual(values = c('#eeeeee','#e48784',"#f8e496",'#eefed4',
                                '#a6dfd3','#e28f66','#f6e4da',"#c2ddd6",
                                '#f2bdd2','#f6e4da'))+
  theme_bw()+
  labs(x="",
      y = "Age  At  Diagnosis",
      color='Germline PTV')+
  geom_hline(yintercept = 40,linetype="dashed")+
  scale_y_continuous(limits = c(20, 70))+
theme(legend.position = c(0.10,0.73),
  panel.grid =element_blank()) +   ## 删去网格线
  theme(axis.text.x = element_blank()) +   ## 删去刻度标签
  theme(axis.ticks.x = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank()) +   ## 删去外层边框
  theme(axis.line.y = element_line(size=0.5, colour = "black"))+## 再加上坐标轴（无刻度、无标签）
  geom_hline(yintercept = 20)
    #guides(color=guide_legend(ncol = 2,order = 0))
age_germ_plot




path_total_gene = c("PTEN","MSH2","MSH3","MSH6","MLH1","MLH3","PMS2","BRCA1","BRCA2",
               "RAD50","PALB2","ATM","FANCA","BAX","SMAD4","NF1","EP300","CTNNB1,CHEK2")
UCEC_gene = normal_germ_IA %>% filter(Gene.refGene %in% total_gene)
sig_count_age_IA$germ = "None"
rownames(sig_count_age_IA)= sig_count_age_IA$patient
for (i in 1:nrow(UCEC_gene)){
    sig_count_age_IA[UCEC_gene$patient[i],"germ"] = UCEC_gene$Gene.refGene[i]
    } 
sig_count_age_IA$germ = factor(sig_count_age_IA$germ,levels = c("None","PTEN","MSH2","PMS2",
                                                                "BRCA2","PALB2","RAD50","NF1",
                                                                "ATM","FANCA"
                                                               ))
sig_count_age_IA_order = sig_count_age_IA%>%
  mutate(order_patient = fct_reorder(patient,age_at_diagnosis))
sig_count_age_IA_order$age_at_diagnosis_na  = sig_count_age_IA_order$age_at_diagnosis
sig_count_age_IA_order[ sig_count_age_IA_order$germ =="None","age_at_diagnosis_na"] = NA


Germ_lodes <- to_lodes_form(sig_count_age_IA_order_input_node_add[,1:ncol(sig_count_age_IA_order_input_node_add)],
                           axes = 1:ncol(sig_count_age_IA_order_input_node_add),
                           id = "Sample")
dim(Germ_lodes)

Germ_lodes_order$stratum = factor(Germ_lodes_order$stratum,levels =c("PTEN","PALB2","ATM","FANCA","BRCA2","RAD50",
                                                        "MSH2","PMS2","NF1",
                                                        "CN low","MSI_hypermutated") )

ggplot(Germ_lodes_order,
       aes(x = x, stratum = stratum, alluvium = Sample,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/4) +
  geom_stratum(alpha = 1,width = 0.14) + 
  geom_text(stat = "stratum", size = 3,color="white") + 
  

  scale_fill_manual(values = mycol_2) +

  xlab("") + ylab("") +
  theme_bw() + 
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
  ggtitle("")+
  guides(fill = FALSE) 

