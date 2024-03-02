#--------------------------------------------------------------
# filename : Proteomics.R
# Date : 2024-03-02
# contributor : Zhe Hu, MD
# function: Proteogenomic analysis of early-onset endometrioid endometrial carcinoma
# R version: R/4.1.2
#--------------------------------------------------------------

library(GseaVis)
library(clusterProfiler)
library(tidyverse)
library(ComplexHeatmap)


hallmark_gmt = "/data/public/database/homo_sapiens/msigdb/h.all.v7.5.1.symbols.gmt"
hallmark = read.gmt(hallmark_gmt)
hallmark$term = gsub("HALLMARK_","",hallmark$term)
hallmark.list = hallmark %>% split(.$term) %>% lapply("[[",2)

gsea_hall_pre <- GSEA(gene_fc,
                  minGSSize = 1,
                  maxGSSize = 10000,
                  #nperm = 100,
                  #eps = 0.0001,
                pvalueCutoff =1,
                  pAdjustMethod ="BH",
             TERM2GENE = hallmark)
head(gsea_hall_pre)

gseaNb(object = gsea_hall_pre,
       geneSetID = 'HALLMARK_FATTY_ACID_METABOLISM',
       arrowType = 'open',
       geneCol = 'black',
      addPval = T)



Figure.6 = ggplot(df_sort_input_five %>% filter(type == "nosign"),aes(logFC, logP))+
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(data = df_sort_input_five %>% filter(type == "nosign"), aes(fill=type,color=type),shape=21,stroke=0.2,alpha = 1,size =7)+
  geom_point(data = df_sort_input_five %>% filter(type %in%  c("up","down")), aes(fill=type,color=type),shape=21,stroke=0.1,alpha = 1,size =7)+
  geom_point(data = df_sort_input_five %>%filter(type %in%  c("more_up","more_down")), aes(fill=type,color=type),shape=21,stroke=0.2,alpha = 1,size =7)+
  scale_fill_manual(values=c('more_down' = '#5DAEC7',
                             'down' = '#ABD3DE',
                             'more_up' = '#DA4C35',
                             'up'='#F0B7AE',
                             'nosign'='#d2dae2'))+
  scale_color_manual(values=c('more_down' = 'black',
                             'down' = '#767777',
                             'more_up' = 'black',
                             'up'='#767777',
                             'nosign'='#d2dae2'))+
 xlim(-5,5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
       #legend.position = c(0.01,0.7),
       legend.justification = c(0,1),
       #text=element_text(family="Arial"
       )+
  guides(col = guide_colorbar(title = "-Log10_q-value"),
        size = "none")+
  geom_text_repel(data= df_sort_input_five,aes(label = label_new),
                  size = 5,box.padding = unit(1, "lines"),
                  point.padding = unit(1, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)+
    xlab("Log2FC")+
    ylab("-Log10(FDR q-value)")+
    labs(title="Early-Onset IA Versus Old-Onset un-KNN")
Figure.6



#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable to be summariezed
# groupnames : vector of column names to be used as grouping variables
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
IHC_Score_Insen_sum <- data_summary(IHC_Score_long_Insen, varname="IHC_Score",groupnames="Gene")
IHC_Score_Sen_sum <- data_summary(IHC_Score_long_Sen, varname="IHC_Score",groupnames="Gene")
IHC_Score_Area = read_xlsx("./data/BY_IHC_Score.xlsx")
colnames(IHC_Score_Area)[1:2] = c("Tumor_Sample_Barcode","CR_Months")
IHC_Score_Area_long = melt(IHC_Score_Area,
                  id.vars=c("Tumor_Sample_Barcode","CR_Months","INTEN_BI"),
                  variable.name="Gene",
                  value.name="IHC_Score"
                  ) %>%
                arrange(Tumor_Sample_Barcode)

IHC_Score_Area_long_Insen = IHC_Score_Area_long %>% filter(INTEN_BI == 1)
IHC_Score_Area_long_Sen   = IHC_Score_Area_long %>% filter(INTEN_BI == 0)
IHC_Score_Area_long_Insen$Gene_Group = paste0(IHC_Score_Area_long_Insen$Gene,"_Insen")
IHC_Score_Area_long_Sen$Gene_Group = paste0(IHC_Score_Area_long_Sen$Gene,"_Sen")
IHC_Score_Area_long_Group = rbind(IHC_Score_Area_long_Insen,IHC_Score_Area_long_Sen)
IHC_Score_Area_long_Group$Gene_Group = factor(IHC_Score_Area_long_Group$Gene_Group, levels = c("EEF1E_Insen","EEF1E_Sen",
                                                                                    "ILVBL_Insen","ILVBL_Sen",
                                                                                    "SRPK1_Insen","SRPK1_Sen",
                                                                                    "NUDT5_Insen","NUDT5_Sen"))
IHC_Score_Area_Insen_sum <- data_summary(IHC_Score_Area_long_Insen, varname="IHC_Score",groupnames="Gene")
IHC_Score_Area_Sen_sum <- data_summary(IHC_Score_Area_long_Sen, varname="IHC_Score",groupnames="Gene")
IHC_Score_Area_Insen_sum$INTEN_BI = 1
IHC_Score_Area_Sen_sum$INTEN_BI   = 0
IHC_Score_Area_Insen_sum$Gene_Group = paste0(IHC_Score_Area_Insen_sum$Gene,"_Insen")
IHC_Score_Area_Sen_sum$Gene_Group = paste0(IHC_Score_Area_Sen_sum$Gene,"_Sen")
IHC_Score_Area_sum_Group = rbind(IHC_Score_Area_Insen_sum,IHC_Score_Area_Sen_sum)
IHC_Score_Area_sum_Group$Gene_Group = factor(IHC_Score_Area_sum_Group$Gene_Group, levels = c("EEF1E_Insen","EEF1E_Sen",
                                                                                    "ILVBL_Insen","ILVBL_Sen",
                                                                                    "SRPK1_Insen","SRPK1_Sen",
                                                                                    "NUDT5_Insen","NUDT5_Sen"))

Figure.8 = ggplot(IHC_Score_Area_long_Group,aes(Gene_Group, IHC_Score))+
  geom_bar(IHC_Score_Area_sum_Group,mapping = aes(x=Gene_Group, y=IHC_Score, fill=Gene_Group),
           stat="identity",width=.8,alpha = 0.4,position=position_dodge(1.6)) +
  geom_errorbar(IHC_Score_Area_sum_Group,mapping = aes(x=Gene_Group, y=IHC_Score,ymin=IHC_Score-sd, ymax=IHC_Score+sd,fill=Gene_Group), 
                width=.2,position=position_dodge(1.2))+
  geom_point(width=0.15,position = 'jitter',aes(fill=Gene_Group),shape=21,size=3.5,color="black",stroke=0.6)+
  geom_signif(IHC_Score_Area_long_Group,mapping = aes(x=Gene_Group, y=IHC_Score),
              comparisons = list(c("EEF1E_Insen","EEF1E_Sen"),c("ILVBL_Insen","ILVBL_Sen"),
                                c("SRPK1_Insen","SRPK1_Sen"),c("NUDT5_Insen","NUDT5_Sen")),
              test = "wilcox.test",
              step_increase = 0,
              tip_length = 0,
              textsize = 6,
              size = 1,
              map_signif_level = F)+
  scale_fill_manual(values=c('#73bbaf','#73bbaf',
                             '#d15e67','#d15e67',
                             '#f5bcaa','#f5bcaa',
                             '#6c43a6','#6c43a6'))+
  scale_color_manual(values=c('#73bbaf','#73bbaf',
                             '#d15e67','#d15e67',
                             '#f5bcaa','#f5bcaa',
                             '#6c43a6','#6c43a6'))+
  scale_y_continuous(expand = c(0,0),limits = c(-0.3,20))+
  labs(x = "", y = "IHC Score")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=14))
Figure.8


## barplot
dat_plot <- data.frame(id = row.names(gsva_hall_diff),
                       t = gsva_hall_diff$t)
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))

dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
library(ggthemes)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, Early-Onset versus Old-Onset') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    )
p



 ggplot(pro_matrix_easy_raw)+
  geom_col(aes(x = 1,fill = unique_pep_type,y=1))+

  scale_fill_manual(values = c('#A6CEE3','#1F78B4','#B2DF8A','#CAB2D6','#6A3D9A'))+
  coord_polar(theta = 'y') +
theme_void()+ scale_x_continuous(expand = expansion(add = c(1.2,0)))+scale_y_reverse()



options(repr.plot.width=12, repr.plot.height=5)


ggplot(pro_matrix_easy_raw %>% filter(Unique.peptides >= 2) , aes(x=as.integer(peptides_new))) +
  geom_histogram(aes(y=..count..),
                color= "black", fill="#4DBBD5",cex=0.4,
                alpha=1,boundary =1.5,
                bins= 20)+
geom_text(aes(label=as.character(round(..count..,2))),
          stat="bin",bins=20,boundary=1.5,vjust=-0.5)+
geom_point(data = pro_matrix_easy_raw_ggplot, color = "#3C5488",size = 5, 
           aes(x=peptides_new,y = cul_pct*1500))+
  scale_x_continuous(limits = c(1,20),
                     breaks = c(2,5,10,15,19),
                     labels = c(2,5,10,15,">19"))+
theme_classic()+
xlab("# of peptides per protein")+
scale_y_continuous(sec.axis = sec_axis(~./1500, name = "percent"),
                  limits = c(0,1500))+
theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
              plot.caption = element_text(size = 12, face = "italic"),
              axis.text = element_text(size=12),  
              axis.title = element_text(size=14, face="bold"))


theme_niwot <- function(){
  theme_bw() +
    theme(
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}




ggplot(Total_pro_input[,!duplicated(names(Total_pro_input))] , aes(x = Sample_type, y = EPCAM)) +
  geom_violin(aes(fill = Sample_type, colour = Sample_type), alpha = 0.5,width =0.8) +
  geom_boxplot(aes(colour = Sample_type), width = 0.2) +
  theme_niwot()+
  ggsci::scale_color_npg()+
  ggsci::scale_fill_npg()+
  ggpubr::stat_compare_means()
