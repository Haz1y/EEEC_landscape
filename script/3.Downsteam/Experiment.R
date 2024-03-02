#--------------------------------------------------------------
# filename : Experiment.R
# Date : 2024-03-02
# contributor : Zhe Hu, MD
# function: Proteogenomic analysis of early-onset endometrioid endometrial carcinoma
# R version: R/4.1.2
#--------------------------------------------------------------


library(ggpubr)
library(plyr)
library(tidyverse)
library(ggbeeswarm)


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


PDO_count_Sum =  PDO_count %>% pivot_longer(cols = -Treat,values_to = "count",names_to = "sample") %>% mutate(Group = gsub('[0-9]','',sample))%>%
mutate(Group_Treat = paste0(Treat,"_",Group)) %>%data_summary(groupnames= "Group_Treat",varname= "count")
PDO_count_long = PDO_count %>% pivot_longer(cols = -Treat,values_to = "count",names_to = "sample") %>% mutate(Group = gsub('[0-9]','',sample))%>%
mutate(Group_Treat = paste0(Treat,"_",Group))
PDO_count_Sum

PDO_count_Sum$Group_Treat = factor(PDO_count_Sum$Group_Treat, levels = c("CON_NC","CON_WT","CON_MUT",
                                                                        "PRG_NC","PRG_WT","PRG_MUT",
                                                                        "MPA_NC","MPA_WT","MPA_MUT"))
PDO_count_long$Group_Treat = factor(PDO_count_long$Group_Treat, levels =  c("CON_NC","CON_WT","CON_MUT",
                                                                        "PRG_NC","PRG_WT","PRG_MUT",
                                                                           "MPA_NC","MPA_WT","MPA_MUT"))

PDO_count_Sum$Group = c("MUT","NC","WT",
                       "MUT","NC","WT",
                       "MUT","NC","WT")
PDO_count_Sum$Treat = c("CON","CON","CON",
                       "MPA","MPA","MPA",
                       "PRG","PRG","PRG")

PDO_count_Sum$Group = factor(PDO_count_Sum$Group,levels =c("NC","WT","MUT"))
PDO_count_long$Group = factor(PDO_count_long$Group,levels =c("NC","WT","MUT"))

PDO_count_long = PDO_count_long%>% arrange(Group_Treat)
PDO_count_Sum = PDO_count_Sum%>% arrange(Group_Treat)

    
anova = aov(count~Group_Treat, PDO_count_long )
summary(anova)
tukey <- TukeyHSD(anova)
tukey <- as.data.frame(tukey$trt)
tukey$pair <- rownames(tukey)
TukeyHSD(anova)


ggplot(PDO_count_long,aes(x_pos, count))+
  geom_bar(PDO_count_Sum,mapping = aes(x=x_pos, y=count, fill=Group),
           stat="identity",width=.45,alpha = 0.4,position=position_dodge(1.6))+
  geom_errorbar(PDO_count_Sum,mapping = aes(x=x_pos, y=count,ymin=count-sd, ymax=count+sd,fill=Group), 
                width=.2,position=position_dodge(1.2))+
  geom_quasirandom(aes(fill=Group),shape=21,width = 0.2,size=3.5,stroke=0.6)+
scale_x_continuous(breaks=c(1,1.5,2,3,3.5,4,5,5.5,6),
                      labels = c("NC","WT","MT","NC","WT","MT","NC","WT","MT"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,160),breaks=c(0,50,100,150))+
  labs(x = "", y = "PDO Counts")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.ticks.length=unit(.25, "cm"),
        axis.text.x = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=14))+
  scale_fill_manual(values=c('#75579F','#F29562','#C4241F'))+
  scale_color_manual(values=c('#75579F','#F29562','#C4241F'))

mice_curve = read_excel("../9.PDO/data_0909.xlsx",sheet=5)
mice_curve = mice_curve[1:7,c(2:12,23:32)]
colnames(mice_curve)[2:6] = paste0("NC_Treat",1:5)
colnames(mice_curve)[7:11] = paste0("NC_Ctrl",1:5)
colnames(mice_curve)[12:16] = paste0("MT_Treat",1:5)
colnames(mice_curve)[17:21] = paste0("MT_Ctrl",1:5)
mice_curve = as.data.frame(mice_curve)
mice_curve[,2] = as.double(mice_curve[,2])


mice_curve %>% pivot_longer(cols = -Days,names_to = "Sample",values_to = "size") %>% 
    mutate(Deatil_Group = gsub('[0-9]','',Sample)) %>% 
    mutate(Deatil_Group = factor(Deatil_Group,levels = c("NC_Ctrl","NC_Treat","MT_Ctrl","MT_Treat"))) 

mice_curve %>% pivot_longer(cols = -Days,names_to = "Sample",values_to = "size") %>% 
    mutate(Deatil_Group = gsub('[0-9]','',Sample)) %>% 
    mutate(Group = ifelse(startsWith(Deatil_Group,"NC"),"NC","MT")) %>%
    mutate(Treat = ifelse(endsWith(Deatil_Group,"Treat"),"Treat","Ctrl")) %>%
    mutate(Deatil_Group = factor(Deatil_Group,levels = c("NC_Ctrl","NC_Treat","MT_Ctrl","MT_Treat"))) %>% 
ggplot(aes(x = Days,y =size,
             group=Deatil_Group,color=Deatil_Group)) + 
  stat_summary(geom = 'line',fun='mean',cex=2)+
  stat_summary(fun.y="mean",geom="line",cex=2) +  
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=1,cex=1)+
  stat_summary(geom = 'point',fun='mean',aes(fill=Deatil_Group),
               size=5,pch=21,color='black')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black",size = 1),
        axis.ticks=element_line(color="black",size=0.8 ),
        axis.ticks.length=unit(.5, "cm"),
        axis.text = element_text(color="black",size=14))+
  theme_classic(base_size = 15)+
  labs(title = "", y="Tumor Size (mm^2)", x = "Day")+
  scale_x_continuous(limits = c(0,22),expand = c(0,0))+
  ggsci::scale_color_npg()+
  ggsci::scale_fill_npg()+
  coord_fixed(1/50)

