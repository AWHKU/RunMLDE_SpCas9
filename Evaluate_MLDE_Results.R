library(ggplot2)
library(tidyverse)
library(ggpmisc)
library(ggseqlogo)
library(readxl)
library(Cairo)
library(stringdist)
### Set plotting format
theme_publication<-theme_classic()+
  theme(text = element_text(size=12, family = "Arial"),
        panel.background = element_rect(colour="black", size=0.5))

Cas9_xl<-read_excel("./41592_2019_473_MOESM4_ESM_CombiSEAL_supp_Fig2.xlsx")
WT="RQKETQKR"
non_target_order=c("---", "-R-", "-H-", "--A", "-RA", "-HA", "-AA", "A--", "AR-", "AH-", "A-A", "ARA", "AHA", "AAA")
target_order=c( "-----", "A----", "-A---", "AA---", "----A", "A---A", "-A--A", "AA--A", "--QGP", 
                "A-QGP", "-AQGP", "AAQGP", "--VRE", "A-VRE", "-AVRE", "AAVRE", "--AWE", "A-AWE", "-AAWE" ,"AAAWE",
                "--MVA", "A-MVA", "-AMVA", "AAMVA", "--KSA", "A-KSA", "-AKSA", "AAKSA", "--RK-", "A-RK-", "-ARK-",  "AARK-", 
                "--CRE", "A-CRE" , "-ACRE", "AACRE", "--QW-", "A-QW-", "-AQW-", "AAQW-", "--LGA", "A-LGA", "-ALGA", "AALGA",
                "--WDE", "A-WDE" , "-AWDE" , "AAWDE", "--HL-", "A-HL-", "-AHL-", "AAHL-", "--VWA", "A-VWA", "-AVWA", "AAVWA",
                "--RRA" , "A-RRA", "-ARRA", "AARRA", "--GDE", "A-GDE", "-AGDE", "AAGDE" , "--MRA", "A-MRA",  "-AMRA", "AAMRA")

#sg5_Fit<-read_csv("/media/achu/新增磁碟區/2021_07_24_MLDE_SpCas9/sg5ON_Fitness.csv")
#sg8_Fit<-read_csv("/media/achu/新增磁碟區/2021_07_24_MLDE_SpCas9/sg8ON_Fitness.csv")

### prepare empirial data
### define true positive - > 70% of WT activity
Cas9_xl<-read_excel("/media/achu/Data/Cas9_MLDE/41592_2019_473_MOESM4_ESM_CombiSEAL_supp_Fig2.xlsx")
WT="RQKETQKR"
colnames(Cas9_xl)<-Cas9_xl[4,]
Cas9_xl<-Cas9_xl[5:956, ]
Cas9_xl[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")]<-sapply(Cas9_xl[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")], as.numeric)
sg5_Fit<-Cas9_xl %>% mutate(AACombo=paste0(`661`, `695`, `848`, `923`, `924`, `926`, `1003`, `1060`, sep="")) %>% select(AACombo, `RFPsg5 ON`)
colnames(sg5_Fit)<-c("AACombo", "Fitness")
sg8_Fit<-Cas9_xl %>% mutate(AACombo=paste0(`661`, `695`, `848`, `923`, `924`, `926`, `1003`, `1060`, sep="")) %>% select(AACombo, `RFPsg8 ON`)
colnames(sg8_Fit)<-c("AACombo", "Fitness")
WT_fit <-sg5_Fit %>% filter(AACombo==WT) %>% select(Fitness) %>% unlist()
sg5_Fit %>% ggplot(aes(x=Fitness)) +geom_histogram(fill=NA, colour="black") +
    geom_vline(xintercept = WT_fit,linetype = "dashed")+
  labs(title="sg5")+theme_publication
sg5_Fit_WT70<-sg5_Fit %>% filter(Fitness >= WT_fit*0.7) %>% select(AACombo) %>% unlist()
sg5_Fit_TN<-sg5_Fit %>% filter(Fitness < WT_fit*0.7) %>% select(AACombo) %>% unlist()
WT_fit <-sg8_Fit %>% filter(AACombo==WT) %>% select(Fitness) %>% unlist()
sg8_Fit %>% ggplot(aes(x=Fitness)) +geom_histogram(fill=NA, colour="black") +
  geom_vline(xintercept = WT_fit,linetype = "dashed")+
  labs(title="sg8")+theme_publication
sg8_Fit_WT70<-sg8_Fit %>% filter(Fitness >= WT_fit*0.7) %>% select(AACombo) %>% unlist()
sg8_Fit_TN<-sg8_Fit %>% filter(Fitness < WT_fit*0.7) %>% select(AACombo) %>% unlist()
WT_fit <-sg8_Fit %>% filter(AACombo==WT) %>% select(Fitness) %>% unlist()
### calculate fold of enrichment 400*Iprediction/N
sg5_top5<-sg5_Fit %>% mutate(rank=percent_rank(Fitness)) %>% filter(rank>=0.95) %>% select(AACombo) %>% unlist()
sg8_top5<-sg8_Fit %>% mutate(rank=percent_rank(Fitness)) %>% filter(rank>=0.95) %>% select(AACombo) %>% unlist()




### survey empirial fitness input
fitness_files<-list.files(pattern = ".csv", path="./SpCas9_diverse_variants")
setwd("./SpCas9_diverse_variants")
AACombo<-Cas9_xl[, 2:9] 
AACombo<-AACombo %>% mutate(AACombo=paste0(`661`,`695`,`848`,`923`,`924`,`926`,`1003`,`1060`))
colnames(AACombo)<-c('R661', 'Q695', 'K848', 'E923', 'T924', 'Q926', 'K1003', 'R1060', "AACombo")
variant_no<-c()
positives<-c()
var_dist<-c()

for (i in 1:length(fitness_files)){
  t<-read_csv(fitness_files[i])
  variant_no<-c(variant_no, nrow(t))
  ### calculate dist occur most frequently ### 
  mat_dist<-stringdistmatrix(t$AACombo, t$AACombo)
  var_dist<-c(var_dist, max(which.max(table(mat_dist))))
  if (grepl("sg5",fitness_files[i] )){
    positives<-c(positives, length(intersect(t$AACombo, sg5_Fit_WT70)))}
  if (grepl("sg8",fitness_files[i] )){
    positives<-c(positives, length(intersect(t$AACombo,sg8_Fit_WT70)))}
  t<-left_join(AACombo, t, by="AACombo")
  t<-t %>% mutate(presence=ifelse(is.na(Fitness), "no", "yes"))
  t$R661[t$R661=="R"]<- "-"
  t$Q695[t$Q695=="Q"]<-"-"
  t$K848[t$K848=="K"]<-"-"
  t$E923[t$E923=="E"]<-"-"
  t$T924[t$T924=="T"]<-"-"
  t$Q926[t$Q926=="Q"]<-"-"
  t$K1003[t$K1003=="K"]<-"-"
  t$R1060[t$R1060=="R"]<-"-"
  t<-t %>% mutate(non_target_strand=paste0(K848,K1003, R1060)) %>% 
    mutate(target_strand=paste0(R661, Q695, E923, T924, Q926)) 
  t$target_strand<-factor(t$target_strand, rev(target_order))
  t$non_target_strand<-factor(t$non_target_strand, non_target_order)
  plot1<-t %>% ggplot(aes(y=target_strand, x=non_target_strand, fill=presence)) +geom_tile(colour="black")+
  scale_fill_manual(values = c("yes"="tomato", "no"="black"))
  ggsave(plot1, filename = gsub("csv", "png", fitness_files[i]), units = "mm", width=120, height = 240, type="cairo")
  
}
fitness_files<-tibble(filename=fitness_files, variant_no, var_dist, positives)
fitness_files$pos_percent<-fitness_files$positives/length(sg5_Fit_WT70)
fitness_files<-fitness_files %>% mutate(datatype=ifelse(grepl("max", filename), "diverse", "random"))
fitness_files<-fitness_files %>% mutate(sgRNA=substring(filename, 8, 12))
pos_tab<-fitness_files %>%group_by(sgRNA, variant_no, datatype) %>% summarise(avg_pos=mean(positives)) %>%
  spread(datatype, avg_pos)


fitness_files %>% group_by( variant_no, dataype) %>% summarise(avg_pos=mean(positives))
### output empirical data info
write_csv(fitness_files, file = "SpCas9_diverse_variants_stats.csv")
summarize_MLDE_models<-function(wkdir, embedding, para){
  setwd(wkdir)
  samplenames<-list.files(pattern = "sg5ON")
  cors<-c()
  rsquare<-c()
  true_pos<-c()
  precision<-c()
  sensitivity<-c()
  specificity<-c()
  NDCG<-c()
  top_variants<-tibble(sgRNA=c(), samplenames=c(), t5=c())
  for (i in (1:length(samplenames))){
    t<-read.csv(paste0("./", samplenames[i],"/","PredictedFitness.csv"))
    t$PredictedFitness<-t$PredictedFitness*(4.273637+2.089742) -2.089742
    t<-full_join(sg5_Fit, t, by="AACombo") 
    c<-cor.test(t$Fitness, t$PredictedFitness)
    reg<-summary(lm(PredictedFitness~Fitness, data = t))
    cors<-c(cors, c$estimate)
    rsquare<-c(rsquare, reg$adj.r.squared)
    top5<-t  %>% mutate(rank=percent_rank(PredictedFitness)) %>% filter(rank>=0.95) %>% select(AACombo) %>% unlist()
    WT_fit<-t  %>% filter(AACombo==WT) %>% select(PredictedFitness) %>% unlist() 
    p<-t %>% ggplot(aes(x=PredictedFitness)) + geom_histogram()+geom_vline(xintercept = WT_fit*0.7)+labs(title=samplenames[i])
    ggsave(p, filename = paste0("./", samplenames[i],"/","PredictedFitness_hist.png"), width = 4, height = 4)
    ###plot heat map
    p<-t %>% select(AACombo, Fitness, PredictedFitness) %>%
      gather(type, value, 2:3) %>% mutate(AA=AACombo) %>% separate(AA, c('p0','R661', 'Q695', 'K848', 'E923', 'T924', 'Q926', 'K1003', 'R1060'), sep="")
    p<-p[, c('AACombo', 'R661', 'Q695', 'K848', 'E923', 'T924', 'Q926', 'K1003', 'R1060', 'type', 'value')]
    p$R661[p$R661=="R"]<- "-"
    p$Q695[p$Q695=="Q"]<-"-"
    p$K848[p$K848=="K"]<-"-"
    p$E923[p$E923=="E"]<-"-"
    p$T924[p$T924=="T"]<-"-"
    p$Q926[p$Q926=="Q"]<-"-"
    p$K1003[p$K1003=="K"]<-"-"
    p$R1060[p$R1060=="R"]<-"-"
    p<-p %>% unite("non_target_strand", c('K848','K1003', 'R1060'), sep="") %>% 
      unite("target_strand", c('R661', 'Q695', 'E923', 'T924', 'Q926'), sep="") 
    p$target_strand<-factor(p$target_strand, rev(target_order))
    p$non_target_strand<-factor(p$non_target_strand, non_target_order)
    #p<- p%>% mutate(value=value*(4.273637+2.089742) -2.089742)
    WT_fit<-p  %>% filter(AACombo==WT)  %>% select(type, value)
    colnames(WT_fit)<-c("type", "WT_fit")
    p<-left_join(p, WT_fit, by="type")%>% mutate(WT_70=ifelse(value>=0.7*WT_fit, "yes", "no"))
    plot1<-p %>%  ggplot(aes(x=non_target_strand, y=target_strand, fill=value))+
      geom_tile()+ 
      geom_tile(data=subset(p, WT_70=="yes"), colour="black") + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, name="Escore")+
      facet_grid(~type)+labs(title=paste0(samplenames[i], "WT_70"))
    ggsave(plot1, filename=paste0("./", samplenames[i],"/","Fitness_Predicted_fitness_WT_70.png"), units = "mm", height=290, width=200, device ="png", type="cairo")
    p<-p %>% group_by(type) %>% mutate(rank=percent_rank(value)) 
    p<-p%>% mutate(top5=ifelse(rank>=0.95, "yes", "no"))
    plot1<-p %>%  ggplot(aes(x=non_target_strand, y=target_strand, fill=value))+
      geom_tile()+ 
      geom_tile(data=subset(p, top5=="yes"), colour="black") + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, name="Escore")+
      facet_grid(~type)+labs(title=paste0(samplenames[i], "top 5%"))
    ggsave(plot1, filename=paste0("./", samplenames[i],"/","Fitness_Predicted_fitness_top5.png"), units = "mm", height=290, width=200, device ="png", type="cairo")
    WT_fit<-t  %>% filter(AACombo==WT) %>% select(PredictedFitness) %>% unlist() 
    variant_70 <- t %>% filter(PredictedFitness>=WT_fit*0.7) %>% select(AACombo) %>% unlist()
    variant_N70 <- t %>% filter(PredictedFitness<WT_fit*0.7) %>% select(AACombo) %>% unlist()
    ###TP, FP, TN, FN
    TP<-length(intersect(variant_70, sg5_Fit_WT70))
    FP<-length(intersect(sg5_Fit_TN, variant_70))
    TN<-length(intersect(variant_N70, sg5_Fit_TN))
    FN<-length(intersect(sg5_Fit_WT70, variant_N70))
    ### calculate sensitivity/recall - TP/TP+FN 
    sensitivity<-c(sensitivity, TP/(TP+FN)*100)
    ### calculate specificity
    specificity<-c(specificity, TN/(TN+FP)*100)
    ### calculate precision TP/TP+FP
    precision<-c(precision, TP/(TP+FP)*100)
    tp<-length(intersect(sg5_top5, top5))
    true_pos<-c(true_pos, tp)
    top_variants<-bind_rows(top_variants, tibble(sgRNA="sg5", samplenames=samplenames[i], t5=top5))
    ### calculate NDCG
    t<-t[!is.na(t$Fitness), ]
    t<-t %>% mutate(F_rank=rank(Fitness), P_rank=rank(PredictedFitness))
    NDCG<-c(NDCG, sum(t$Fitness/log2(t$P_rank +1))/sum(t$Fitness/log2(t$F_rank +1)))
    }
  NGS_stats<-tibble(embedding=embedding, para=para, samplenames, cors, rsquare, true_pos, sensitivity, specificity, precision, NDCG)
  
  samplenames<-list.files(pattern = "sg8ON")
  cors<-c()
  rsquare<-c()
  true_pos<-c()
  precision<-c()
  sensitivity<-c()
  specificity<-c()
  NDCG<-c()
  top_variants<-tibble(sgRNA=c(), samplenames=c(), t5=c())
  for (i in (1:length(samplenames))){
    t<-read.csv(paste0("./", samplenames[i],"/","PredictedFitness.csv"))
    t$PredictedFitness<-t$PredictedFitness*(2.760273+7.206505) - 7.206505
    t<-full_join(sg8_Fit, t, by="AACombo") 
    c<-cor.test(t$Fitness, t$PredictedFitness)
    reg<-summary(lm(PredictedFitness~Fitness, data = t))
    cors<-c(cors, c$estimate)
    rsquare<-c(rsquare, reg$adj.r.squared)
    top5<-t  %>% mutate(rank=percent_rank(PredictedFitness)) %>% filter(rank>=0.95) %>% select(AACombo) %>% unlist()
    WT_fit<-t  %>% filter(AACombo==WT) %>% select(PredictedFitness) %>% unlist() 
    p<-t %>% ggplot(aes(x=PredictedFitness)) + geom_histogram()+geom_vline(xintercept = WT_fit*0.7)+labs(title=samplenames[i])
    ggsave(p, filename = paste0("./", samplenames[i],"/","PredictedFitness_hist.png"), width = 4, height = 4)
    ###plot heat map
    p<-t %>% select(AACombo, Fitness, PredictedFitness) %>%
      gather(type, value, 2:3) %>% mutate(AA=AACombo) %>% separate(AA, c('p0','R661', 'Q695', 'K848', 'E923', 'T924', 'Q926', 'K1003', 'R1060'), sep="")
    p<-p[, c('AACombo', 'R661', 'Q695', 'K848', 'E923', 'T924', 'Q926', 'K1003', 'R1060', 'type', 'value')]
    p$R661[p$R661=="R"]<- "-"
    p$Q695[p$Q695=="Q"]<-"-"
    p$K848[p$K848=="K"]<-"-"
    p$E923[p$E923=="E"]<-"-"
    p$T924[p$T924=="T"]<-"-"
    p$Q926[p$Q926=="Q"]<-"-"
    p$K1003[p$K1003=="K"]<-"-"
    p$R1060[p$R1060=="R"]<-"-"
    p<-p %>% unite("non_target_strand", c('K848','K1003', 'R1060'), sep="") %>% 
      unite("target_strand", c('R661', 'Q695', 'E923', 'T924', 'Q926'), sep="") 
    p$target_strand<-factor(p$target_strand, rev(target_order))
    p$non_target_strand<-factor(p$non_target_strand, non_target_order)
    #p<- p%>% mutate(value=value*(2.760273+7.206505) - 7.206505)
    WT_fit<-p  %>% filter(AACombo==WT)  %>% select(type, value)
    colnames(WT_fit)<-c("type", "WT_fit")
    p<-left_join(p, WT_fit, by="type")%>% mutate(WT_70=ifelse(value>=0.7*WT_fit, "yes", "no"))
    plot1<-p %>%  ggplot(aes(x=non_target_strand, y=target_strand, fill=value))+
      geom_tile()+ 
      geom_tile(data=subset(p, WT_70=="yes"), colour="black") + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, name="Escore")+
      facet_grid(~type)+labs(title=paste0(samplenames[i], "WT_70"))
    ggsave(plot1, filename=paste0("./", samplenames[i],"/","Fitness_Predicted_fitness_WT_70.png"), units = "mm", height=290, width=200, device ="png", type="cairo")
    p<-p %>% group_by(type) %>% mutate(rank=percent_rank(value)) 
    p<-p%>% mutate(top5=ifelse(rank>=0.95, "yes", "no"))
    plot1<-p %>%  ggplot(aes(x=non_target_strand, y=target_strand, fill=value))+
      geom_tile()+ 
      geom_tile(data=subset(p, top5=="yes"), colour="black") + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, name="Escore")+
      facet_grid(~type)+labs(title=paste0(samplenames[i], "top 5%"))
    ggsave(plot1, filename=paste0("./", samplenames[i],"/","Fitness_Predicted_fitness_top5.png"), units = "mm", height=290, width=200, device ="png", type="cairo")
    WT_fit<-t  %>% filter(AACombo==WT) %>% select(PredictedFitness) %>% unlist() 
    variant_70 <- t %>% filter(PredictedFitness>=WT_fit*0.7) %>% select(AACombo) %>% unlist()
    variant_N70 <- t %>% filter(PredictedFitness<WT_fit*0.7) %>% select(AACombo) %>% unlist()
    ###TP, FP, TN, FN
    TP<-length(intersect(variant_70, sg8_Fit_WT70))
    FP<-length(intersect(sg8_Fit_TN, variant_70))
    TN<-length(intersect(variant_N70, sg8_Fit_TN))
    FN<-length(intersect(sg8_Fit_WT70, variant_N70))
    ### calculate sensitivity/recall - TP/TP+FN 
    sensitivity<-c(sensitivity, TP/(TP+FN)*100)
    ### calculate specificity
    specificity<-c(specificity, TN/(TN+FP)*100)
    ### calculate precision TP/TP+FP
    precision<-c(precision, TP/(TP+FP)*100)
    tp<-length(intersect(sg8_top5, top5))
    true_pos<-c(true_pos, tp)
    top_variants<-bind_rows(top_variants, tibble(sgRNA="sg8", samplenames=samplenames[i], t5=top5))
    ### calculate NDCG
    t<-t[!is.na(t$Fitness), ]
    t<-t %>% mutate(F_rank=rank(Fitness), P_rank=rank(PredictedFitness))
    NDCG<-c(NDCG, sum(t$Fitness/log2(t$P_rank +1))/sum(t$Fitness/log2(t$F_rank +1)))
    }
  NGS_stats<-bind_rows(NGS_stats, 
                       tibble(embedding=embedding, para=para, samplenames, cors, rsquare, true_pos, sensitivity, specificity, precision, NDCG))
  write_csv(NGS_stats, file="NGS_stats.csv", quote=F)
  write_csv(top_variants, file="Top5percent_variants.csv", quote=F)
}

wkdir<-"/media/achu/Data/2021_07_24_MLDE_SpCas9/2021_07_29/Belper/ExecuteMlde_para1"
summarize_MLDE_models(wkdir, "Belper", "para1")

wkdir<-"/media/achu/Data/2021_07_24_MLDE_SpCas9/2021_07_29/Belper/ExecuteMlde_para2"
summarize_MLDE_models(wkdir, "Belper", "para2")

wkdir<-"/media/achu/Data/2021_07_24_MLDE_SpCas9/2021_07_29/Georgiev/ExecuteMlde_para1"
summarize_MLDE_models(wkdir, "Georgiev", "para1")

wkdir<-"/media/achu/Data/2021_07_24_MLDE_SpCas9/2021_07_29/Georgiev/ExecuteMlde_para2"
summarize_MLDE_models(wkdir, "Georgiev", "para2")
### analyse Stats
stats_tab<-bind_rows(
  read_csv("/path/to/MLDE/results/Belper/ExecuteMlde_para1/NGS_stats.csv") %>% mutate(embedding="Belper", para="para1"),
  read_csv("/path/to/MLDE/results/Belper/ExecuteMlde_para2/NGS_stats.csv") %>% mutate(embedding="Belper", para="para2"),
  read_csv("/path/to/MLDE/results/Georgiev/ExecuteMlde_para1/NGS_stats.csv")%>% mutate(embedding="Georgiev", para="para1"),
  read_csv("/path/to/MLDE/results/Georgiev/ExecuteMlde_para2/NGS_stats.csv")%>% mutate(embedding="Georgiev", para="para2")
)
stats_tab<-stats_tab %>% mutate(sgRNA=substring(samplenames, 8, 12))
stats_tab<-stats_tab %>% mutate(dataset=ifelse(grepl("max", samplenames), 
                                               substring(samplenames, 22, nchar(samplenames)-5),
                                               substring(samplenames, 22, 23)))
### enrichment fold 400*true_pos/952
stats_tab <- stats_tab %>% mutate(enrichment_fold=400*true_pos/952)


### compare diverse vs random
stats_tab<-stats_tab %>% mutate(datatype=ifelse(grepl("maxnei", samplenames),"diverse", "random"))
fitness_files<-read_csv("./SpCas9_diverse_variants_stats.csv")
fitness_files$filename<-gsub(".csv", "", fitness_files$filename)
colnames(fitness_files)[1]<-c("samplenames")
stats_tab<-left_join(stats_tab, fitness_files, by = "samplenames")
### evaluate performance metrics
stats_tab %>% ggplot(aes(x=datatype, y=specificity))+ geom_boxplot()+facet_grid(interaction(embedding,para)~sgRNA)+theme_publication
stats_tab %>% ggplot(aes(x=datatype, y=precision))+ geom_boxplot()+facet_grid(interaction(embedding,para)~sgRNA)+theme_publication
stats_tab %>% ggplot(aes(x=datatype, y=sensitivity))+ geom_boxplot()+facet_grid(interaction(embedding,para)~sgRNA)+theme_publication
stats_tab %>% ggplot(aes(x=datatype, y=enrichment_fold))+ geom_boxplot()+facet_grid(interaction(embedding,para)~sgRNA)+theme_publication
stats_tab %>% ggplot(aes(x=datatype, y=NDCG))+ geom_boxplot()+facet_grid(interaction(embedding,para)~sgRNA)+theme_publication

### Figure S1
stats_tab %>%  
  select(sgRNA, para, embedding, datatype, precision, sensitivity, specificity, enrichment_fold) %>%
  gather(measure, value, 5:8) %>%ggplot(aes(x=interaction(embedding, para), y=value, group=datatype, colour=datatype))+
  geom_point(position=position_dodge(width=0.7))+facet_grid(measure~sgRNA, scales = "free_y")+
  theme_publication+theme(axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1))

stats_tab %>%  filter(variant_no>60) %>%
  select(sgRNA, para, embedding, datatype, precision, sensitivity, specificity, enrichment_fold, NDCG) %>%
  gather(measure, value, 5:9) %>%
  ggplot(aes(x=interaction( datatype, sgRNA), group=interaction(sgRNA, datatype), y=value, group=datatype, colour=datatype))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.7))+facet_grid(measure~interaction(embedding, para), scales = "free_y")+
  theme_publication+theme(axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1))

### Fig S1. 
stats_tab %>%  filter(sgRNA=="sg5ON", var_no>60) %>% 
  select(sgRNA, para, embedding, datatype, precision, sensitivity, specificity) %>%
  gather(measure, value, 5:7) %>%ggplot(aes(x=datatype, y=value, group=datatype))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.7), alpha=0.7)+
  facet_grid(measure~interaction(embedding, para))+
  theme_publication+theme(axis.title.x = element_blank())+
  labs(y="percentage")

### sg8 only
stats_tab %>%  filter(sgRNA=="sg8ON") %>% 
  select(sgRNA, para, embedding, datatype, precision, sensitivity, specificity, enrichment_fold) %>%
  gather(measure, value, 5:8) %>%ggplot(aes(x=datatype, y=value, group=datatype))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.7), alpha=0.7)+facet_grid(measure~interaction(embedding, para))+
  theme_publication+theme(axis.title.x = element_blank())+
  labs(y="percentage")

### Figure 1B
stats_tab %>% ggplot(aes(x=variant_no, y=specificity, group=interaction(variant_no, datatype), colour=datatype, group=datatype)) +
  geom_point(position=position_jitterdodge(dodge.width = 10, jitter.width=1, jitter.height = 1), alpha=0.3)+
  facet_grid(sgRNA~interaction(embedding, para))+
  labs(x="number of input datapoint")

stats_tab %>% ggplot(aes(x=variant_no, y=precision, group=interaction(variant_no, datatype), colour=datatype, group=datatype)) +
  geom_point(position=position_jitterdodge(dodge.width = 10, jitter.width=1, jitter.height = 1), alpha=0.3)+
  facet_grid(sgRNA~interaction(embedding, para))+
  labs(x="number of input datapoint")


stats_tab %>% ggplot(aes(x=variant_no, y=sensitivity, group=interaction(variant_no, datatype), colour=datatype, group=datatype)) +
  geom_point(position=position_jitterdodge(dodge.width = 10, jitter.width=1, jitter.height = 1), alpha=0.3)+
  facet_grid(sgRNA~interaction(embedding, para))+
  labs(x="number of input datapoint")

stats_tab %>% ggplot(aes(x=variant_no, y=true_pos/65 *100, group=interaction(variant_no, datatype), colour=datatype)) +
  geom_point(position=position_jitterdodge(dodge.width = 10, jitter.width=1, jitter.height = 1), alpha=0.3)+
  facet_grid(sgRNA~interaction(embedding, para))+
  labs(x="number of input datapoint", y="overlap with top 5% in empirical data")


stats_tab %>% group_by(embedding, para, sgRNA) %>% summarize(mean(specificity), mean(precision), mean(sensitivity))

stats_tab %>% filter(sgRNA=="sg5ON", para=="para1", embedding=="Belper", datatype=="random") %>% 
  select(datatype, embedding, para, variant_no, precision, sensitivity, specificity) %>%
  gather(measure, value, 5:7) %>% ggplot(aes(x=variant_no, y=value, group=measure)) +
  geom_point(position=position_jitter(width=1, height = 1), alpha=0.3)+
  facet_grid(measure~.)+
  labs(x="number of input datapoint", y="percentage")+
  theme_publication+
  theme(strip.background = element_blank())


stats_tab %>% filter(sgRNA=="sg5ON", para=="para1", embedding=="Belper", datatype=="random") %>% 
  select(datatype, embedding, para, variant_no, precision, sensitivity, enrichment_fold) %>%
  gather(measure, value, 5:7) %>% ggplot(aes(x=variant_no, y=value, group=measure)) +
  geom_point(position=position_jitter(width=1, height = 1), alpha=0.3)+
  facet_grid(measure~., scale="free_y")+
  labs(x="number of input datapoint", y="percentage")+
  theme_publication+
  theme(strip.background = element_blank())
  
### make summary of sg5 and sg8 figures ###
### overlap between sg5 vs sg8 ###
### var_no<-c(1130, 325, 455, 65, 146, 365, 510, 73)
setwd("/media/achu/Data/2021_07_24_MLDE_SpCas9/2021_07_29/Belper/ExecuteMlde_para1")
focus_files<-list.files(pattern="rep")
fitnesses<-tibble()
for (i in focus_files){
  fitnesses<-bind_rows(fitnesses,
                       read_csv(paste0(i, "/PredictedFitness.csv")) %>% mutate(samplenames=i))
}
fitnesses<-fitnesses %>% mutate(datatype=ifelse(grepl("max", samplenames), "diverse", "random"))
fitnesses<-fitnesses %>% mutate(sgRNA=substring(samplenames, 8, 12))
fitnesses<-left_join(fitnesses%>% ungroup(), stats_tab %>% select(variant_no,samplenames)%>% ungroup(), by="samplenames") 
fitnesses<-fitnesses[!duplicated(fitnesses), ]
fitnesses<-fitnesses %>% mutate(value=ifelse(sgRNA=="sg5ON", 
                                PredictedFitness*(4.273637+2.089742) -2.089742,
                                PredictedFitness*(2.760273+7.206505) - 7.206505))


### Figure S1a
## highlight WT_70
WT_fit<-fitnesses %>% filter(AACombo==WT)%>% select(samplenames, value)
colnames(WT_fit)<-c("samplenames", "WT")
fitnesses<-left_join(fitnesses, WT_fit, by="samplenames")
fitnesses<-fitnesses %>% mutate(WT_70=ifelse(value >= 0.7*WT, "yes", "no"))
### sg5_10
fitnesses$common_top<-c("no")
sg5_5_top<-fitnesses %>% filter(variant_no==33, sgRNA=="sg5ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==33&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_5_top)[table(sg5_5_top)>=1])] <-"yes"


sg5_10_top<-fitnesses %>% filter(variant_no==65, sgRNA=="sg5ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==65&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_10_top)[table(sg5_10_top)>=1])] <-"yes"
sg5_10_top<-fitnesses %>% filter(variant_no==65, sgRNA=="sg5ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==65&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_10_top)[table(sg5_10_top)>=1])] <-"yes"
sg5_20_top<-fitnesses %>% filter(variant_no==130, sgRNA=="sg5ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==130&fitnesses$sgRNA=="sg5ON"&fitnesses$datatype=="diverse"&fitnesses$AACombo %in% names(table(sg5_20_top)[table(sg5_20_top)>=1])] <-"yes"
sg5_20_top<-fitnesses %>% filter(variant_no==130, sgRNA=="sg5ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==130&fitnesses$sgRNA=="sg5ON"&fitnesses$datatype=="random"&fitnesses$AACombo %in% names(table(sg5_20_top)[table(sg5_20_top)>=1])] <-"yes"
sg5_50_top<-fitnesses %>% filter(variant_no==325, sgRNA=="sg5ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==325&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_50_top)[table(sg5_50_top)>=1])] <-"yes"
sg5_50_top<-fitnesses %>% filter(variant_no==325, sgRNA=="sg5ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==325&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_50_top)[table(sg5_50_top)>=1])] <-"yes"
sg5_70_top<-fitnesses %>% filter(variant_no==455, sgRNA=="sg5ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==455&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_70_top)[table(sg5_70_top)>=1])] <-"yes"
sg5_70_top<-fitnesses %>% filter(variant_no==455, sgRNA=="sg5ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==455&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg5ON"&fitnesses$AACombo %in% names(table(sg5_70_top)[table(sg5_70_top)>=1])] <-"yes"

sg8_10_top<-fitnesses %>% filter(variant_no==73, sgRNA=="sg8ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==73&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_10_top)[table(sg8_10_top)>=1])] <-"yes"
sg8_10_top<-fitnesses %>% filter(variant_no==73, sgRNA=="sg8ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==73&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_10_top)[table(sg8_10_top)>=1])] <-"yes"
sg8_20_top<-fitnesses %>% filter(variant_no==146, sgRNA=="sg8ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==146&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_20_top)[table(sg8_20_top)>=1])] <-"yes"
sg8_20_top<-fitnesses %>% filter(variant_no==146, sgRNA=="sg8ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==146&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_20_top)[table(sg8_20_top)>=1])] <-"yes"
sg8_50_top<-fitnesses %>% filter(variant_no==365, sgRNA=="sg8ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==365&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_50_top)[table(sg8_50_top)>=1])] <-"yes"
sg8_50_top<-fitnesses %>% filter(variant_no==365, sgRNA=="sg8ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==365&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_50_top)[table(sg8_50_top)>=1])] <-"yes"
sg8_70_top<-fitnesses %>% filter(variant_no==510, sgRNA=="sg8ON", WT_70=="yes", datatype=="diverse") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==510&fitnesses$datatype=="diverse"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_70_top)[table(sg8_70_top)>=1])] <-"yes"
sg8_70_top<-fitnesses %>% filter(variant_no==510, sgRNA=="sg8ON", WT_70=="yes", datatype=="random") %>% select(AACombo) %>% unlist()
fitnesses$common_top[fitnesses$variant_no==510&fitnesses$datatype=="random"&fitnesses$sgRNA=="sg8ON"&fitnesses$AACombo %in% names(table(sg8_70_top)[table(sg8_70_top)>=1])] <-"yes"

### check what are at 10%
fitnesses %>% filter(variant_no == 65, datatype=="random")  %>% mutate(rank=percent_rank(value)) %>%
  filter(rank >=0.95 &WT_70=="yes") %>% group_by(samplenames) %>% summarize(n())
### how many are top5 $ in TP?
fitnesses %>% filter(variant_no == 65, datatype=="random")  %>% mutate(rank=percent_rank(value)) %>%
  filter(rank >=0.95 & AACombo%in%sg5_Fit_WT70 ) %>% group_by(samplenames) %>% summarize(n())

fitnesses %>% filter(variant_no == 65, datatype=="random")  %>% mutate(rank=percent_rank(value)) %>%
  filter(WT_70=="yes") %>% group_by(samplenames) %>% summarize(n())


fitnesses %>% filter(sgRNA=="sg5ON", datatype=="random")  %>% mutate(rank=percent_rank(value)) %>%
  filter(value >=0.95 & AACombo %in% sg5_Fit_WT70) %>% group_by(samplenames) %>% summarise(counts=n()) 

fitnesses %>% filter(sgRNA=="sg8ON", datatype=="random")  %>% mutate(rank=percent_rank(value)) %>%
  filter(value >=0.95 & AACombo %in% sg8_Fit_WT70) %>% group_by(samplenames) %>% summarise(counts=n())

t<- fitnesses %>% filter(datatype=="random") %>% group_by(variant_no, sgRNA, AACombo, common_top)%>%
  summarise(PredictedFitness=mean(PredictedFitness))

t<- bind_rows(t, 
sg5_Fit %>% filter(!is.na(Fitness))%>% mutate(variant_no=650, sgRNA="sg5ON", common_top=ifelse(AACombo %in%sg5_Fit_WT70, "yes", "no"), PredictedFitness=Fitness) %>% select(variant_no, sgRNA, AACombo, common_top, PredictedFitness),
sg8_Fit %>% filter(!is.na(Fitness))%>% mutate(variant_no=729, sgRNA="sg8ON", common_top=ifelse(AACombo %in%sg8_Fit_WT70, "yes", "no"), PredictedFitness=Fitness) %>% select(variant_no, sgRNA, AACombo, common_top, PredictedFitness)
)



p<-t %>%  mutate(AA=AACombo) %>% separate(AA, c('p0','R661', 'Q695', 'K848', 'E923', 'T924', 'Q926', 'K1003', 'R1060'), sep="")
p$R661[p$R661=="R"]<- "-"
p$Q695[p$Q695=="Q"]<-"-"
p$K848[p$K848=="K"]<-"-"
p$E923[p$E923=="E"]<-"-"
p$T924[p$T924=="T"]<-"-"
p$Q926[p$Q926=="Q"]<-"-"
p$K1003[p$K1003=="K"]<-"-"
p$R1060[p$R1060=="R"]<-"-"
p<-p %>% unite("non_target_strand", c('K848','K1003', 'R1060'), sep="") %>% 
  unite("target_strand", c('R661', 'Q695', 'E923', 'T924', 'Q926'), sep="") 
p$target_strand<-factor(p$target_strand, rev(target_order))
p$non_target_strand<-factor(p$non_target_strand, non_target_order)
p %>%  ggplot(aes(x=non_target_strand, y=target_strand, fill=common_top))+
  geom_tile()+ 
  geom_tile(data=subset(p, common_top=="yes"), colour="black") + 
  scale_fill_manual(values = c("yes"="tomato", "no"="black"), na.value="grey")+
  facet_wrap(sgRNA~variant_no, nrow=2)+theme_publication+
  theme(legend.position = "none",
        panel.background = element_rect(fill="grey"),
        axis.text = element_blank(),
        axis.ticks = element_blank())

p %>% filter(sgRNA=="sg5ON") %>%
  ggplot(aes(x=non_target_strand, y=target_strand, fill=common_top))+
  geom_tile()+ 
  geom_tile(data=subset(p, common_top=="yes"&sgRNA=="sg5ON"), colour="black") + 
  scale_fill_manual(values = c("yes"="tomato", "no"="black"), na.value="grey")+
  facet_wrap(~variant_no, nrow=1)+theme_publication+
  theme(legend.position = "none",
        panel.background = element_rect(fill="grey"),
        axis.text = element_blank(),
        axis.ticks = element_blank())

