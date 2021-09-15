### may have to run again as the min-max normalisation of sp is wrong ###
library(ggplot2)
library(tidyverse)
library(ggpmisc)
library(ggseqlogo)
library(readxl)
library(Cairo)
library(stringdist)
setwd(".")
muts<-read_excel("./41592_2019_473_MOESM4_ESM_CombiSEAL_supp_Fig2.xlsx")
colnames(muts)<-muts[4,]
muts<-muts[5:956, ]
muts[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")]<-sapply(muts[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")], as.numeric)
sg5_Fit<-muts %>% unite(AACombo , c(`661`,`695`,`848`,`923`,`924`,`926`,`1003`,`1060`), sep = "") 
sg5_Fit<-sg5_Fit%>% dplyr::select(AACombo, `RFPsg5 ON`)
colnames(sg5_Fit)<-c("AACombo", "Fitness")
sg8_Fit<-muts %>% unite(AACombo , c(`661`,`695`,`848`,`923`,`924`,`926`,`1003`,`1060`), sep = "") %>% dplyr::select(AACombo, `RFPsg8 ON`)
colnames(sg8_Fit)<-c("AACombo", "Fitness")
### now do min max normalisation (= (x-min)/(max-min))
sg5_Fit <-sg5_Fit[!is.na(sg5_Fit$Fitness), ]
### min max normalise
sg5_Fit$Fitness<-(sg5_Fit$Fitness-min(sg5_Fit$Fitness))/(max(sg5_Fit$Fitness)-min(sg5_Fit$Fitness))
sg5_Fit %>% ggplot(aes(x=Fitness)) + geom_freqpoly()
sg8_Fit <-sg8_Fit[!is.na(sg8_Fit$Fitness), ]
### min max nomalise
sg8_Fit$Fitness<-(sg8_Fit$Fitness-min(sg8_Fit$Fitness))/(max(sg8_Fit$Fitness)-min(sg8_Fit$Fitness))
sg8_Fit %>% ggplot(aes(x=Fitness)) + geom_freqpoly()
WT="RQKETQKR"

write_csv(sg5_Fit, "sg5ON_Fitness.csv")
write_csv(sg8_Fit, "sg8ON_Fitness.csv")
### make diverse variants
### prepare output file
setwd("./SpCas9_diverse_variants")
AACombo<-muts %>% unite(AACombo , c(`661`,`695`,`848`,`923`,`924`,`926`,`1003`,`1060`), sep = "") %>% select(AACombo) %>% unlist()
### randomized data 5%
write_csv(sample_n(sg5_Fit, 33), file = paste0("SpCas9_sg5ON_fitness_33", as.character(var_no), "rep1", ".csv"))
write_csv(sample_n(sg5_Fit, 33), file = paste0("SpCas9_sg5ON_fitness_33", as.character(var_no), "rep2", ".csv"))
write_csv(sample_n(sg5_Fit, 33), file = paste0("SpCas9_sg5ON_fitness_33", as.character(var_no), "rep3", ".csv"))

### diverse datasets
K<-sg5_Fit$AACombo
K_mat<-stringdistmatrix(K, K)
mean(apply(K_mat, 2, function(x) which.max(table(x))))
### 6.313846
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### sg5 variants have average 6 1-mismatch neighbors
K<-sample(K, 65)
K_mat<-stringdistmatrix(K, K)
mean(apply(K_mat, 2, function(x) which.max(table(x))))
### 5.6
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### 0.4
mean(colSums(K_mat==2))
### 1.5

### try to boost diversity by retricting 1-hd and 2-hd neighbor
### 10% = 65
K<-sample(sg5_Fit$AACombo, 100)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=1& colSums(K_mat==2)<=1]
while(length(K[K %in% sg5_Fit$AACombo])<65){
  K<-c(K, sample(sg5_Fit$AACombo, 5))
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1]
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==1)<=1&colSums(K_mat==2)<=1 ]
}
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
K<-sample(K, size=65)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg5_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg5ON_fitness_maxneigbor_65_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg5_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg5_Fit[sg5_Fit$AACombo %in% c(WT, r$AACombo), ]
}
write_csv(r, file = paste0("SpCas9_sg5ON_fitness_65", as.character(var_no), reps[d], ".csv"))

### 20% = 130
K<-sg5_Fit$AACombo
K<-sample(K, 130)
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
### 1.05
mean(colSums(K_mat==2))
### 3.9
mean(colSums(K_mat==3))
### 10

K<-sample(sg5_Fit$AACombo, 100)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=1& colSums(K_mat==2)<=3]
while(length(K[K %in% sg5_Fit$AACombo])<130){
  K<-c(K, sample(sg5_Fit$AACombo, 10))
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1]
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==1)<=3&colSums(K_mat==2)<=5&colSums(K_mat==3)<=14 ]
}
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
mean(colSums(K_mat==3))
mean(colSums(K_mat==4))
K<-sample(K, size=130)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg5_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg5ON_fitness_maxneigbor_130_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg5_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg5_Fit[sg5_Fit$AACombo %in% c(WT, r$AACombo), ]
}
write_csv(r, file = paste0("SpCas9_sg5ON_fitness_130", as.character(var_no), reps[d], ".csv"))

### 50% =  325
K<-sg5_Fit$AACombo
K<-sample(K, 325)
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
### 2.79
mean(colSums(K_mat==2))
### 9.3
mean(colSums(K_mat==3))
### 25.6
K<-sample(sg5_Fit$AACombo, 325)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=5& colSums(K_mat==2)<=12]
while(length(K)<325){
  K<-c(K, sample(sg5_Fit$AACombo, 10))
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1]
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==1)<=5&colSums(K_mat==2)<=14 ]
}
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
K<-sample(K, size=325)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg5_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg5ON_fitness_maxneigbor_325_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg5_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg5_Fit[sg5_Fit$AACombo %in% c(WT, r$AACombo), ]
}
write_csv(r, file = paste0("SpCas9_sg5ON_fitness_325", as.character(var_no), reps[d], ".csv"))

### 70% = 455
K<-sg5_Fit$AACombo
K<-sample(K, 455)
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==0))
mean(colSums(K_mat==1))
### 3.96
max(colSums(K_mat==1))
mean(colSums(K_mat==2))
max(colSums(K_mat==2))
### 13.37
mean(colSums(K_mat==3))
### 35
K<-sample(sg5_Fit$AACombo, 450)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=5]
while(length(K)<455){
  K<-sample(sg5_Fit$AACombo, 465)
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1&colSums(K_mat==1)<=6&colSums(K_mat==2)<=20]
}


K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
K<-sample(K, size=455)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg5_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg5ON_fitness_maxneigbor_455_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg5_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg5_Fit[sg5_Fit$AACombo %in% c(WT, r$AACombo), ]
}
R<-r$AACombo
r_mat<-stringdistmatrix(r$AACombo, r$AACombo)
R<-R[colSums(r_mat==0)<=1]
r_mat<-stringdistmatrix(R, R)
mean(colSums(K_mat==1))
### 3.96
max(colSums(K_mat==1))
mean(colSums(K_mat==2))
max(colSums(K_mat==2))

write_csv(r, file = paste0("SpCas9_sg5ON_fitness_455", as.character(var_no), reps[d], ".csv"))
### 455 cannot make the dataset anymore diverse than random sampling


### sg5 - 65, 130, 325, 455
### sg8 - 73, 146, 365, 510
### 10%
K<-sg8_Fit$AACombo
K_mat<-stringdistmatrix(K, K)
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### sg5 variants have average 6 1-mismatch neighbors
K<-sample(K, 73)
K_mat<-stringdistmatrix(K, K)
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### 0.5
mean(colSums(K_mat==2))
### 2

### try to boost diversity by retricting 1-hd and 2-hd neighbor
K<-sample(sg8_Fit$AACombo, 100)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=1& colSums(K_mat==2)<=1]
while(length(K[K %in% sg8_Fit$AACombo])<73){
  K<-c(K, sample(sg5_Fit$AACombo, 5))
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1]
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==1)<=1&colSums(K_mat==2)<=2 ]
}
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
K<-sample(K, size=73)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg8_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg8ON_fitness_maxneigbor_73_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg8_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg8_Fit[sg8_Fit$AACombo %in% c(WT, r$AACombo), ]
}
R<-r$AACombo
r_mat<-stringdistmatrix(r$AACombo, r$AACombo)
R<-R[colSums(r_mat==0)<=1]
r_mat<-stringdistmatrix(R, R)
mean(colSums(r_mat==1))
max(colSums(r_mat==1))
mean(colSums(r_mat==2))
max(colSums(r_mat==2))
write_csv(r, file = paste0("SpCas9_sg8ON_fitness_73", as.character(var_no), reps[d], ".csv"))


### 20% = 146
K<-sample(sg8_Fit$AACombo, 146)
K_mat<-stringdistmatrix(K, K)
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### 1.2
mean(colSums(K_mat==2))
### 4.6
max(colSums(K_mat==2))
K<-sample(sg8_Fit$AACombo, 146)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=2& colSums(K_mat==2)<=5]
while(length(K[K %in% sg8_Fit$AACombo])<146){
  K<-c(K, sample(sg5_Fit$AACombo, 20))
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1]
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==1)<=2&colSums(K_mat==2)<=6 ]
}
K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
max(colSums(K_mat==1))
max(colSums(K_mat==2))
K<-sample(K, size=146)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg8_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg8ON_fitness_maxneigbor_146_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg8_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg8_Fit[sg8_Fit$AACombo %in% c(WT, r$AACombo), ]
}
R<-r$AACombo
r_mat<-stringdistmatrix(r$AACombo, r$AACombo)
R<-R[colSums(r_mat==0)<=1]
r_mat<-stringdistmatrix(R, R)
mean(colSums(r_mat==1))
max(colSums(r_mat==1))
mean(colSums(r_mat==2))
max(colSums(r_mat==2))
write_csv(r, file = paste0("SpCas9_sg8ON_fitness_146", as.character(var_no), reps[d], ".csv"))

### 50% = 365
K<-sample(sg8_Fit$AACombo,365)
K_mat<-stringdistmatrix(K, K)
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### 2.9
mean(colSums(K_mat==2))
### 9.9
max(colSums(K_mat==1))
max(colSums(K_mat==2))
K<-sample(sg8_Fit$AACombo, 365)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=3& colSums(K_mat==2)<=8]

while(length(K)<365){
  K<-sample(sg8_Fit$AACombo, 365)
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1&colSums(K_mat==1)<=6&colSums(K_mat==2)<=17]
}

K_mat<-stringdistmatrix(K, K)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
max(colSums(K_mat==1))
max(colSums(K_mat==2))
K<-sample(K, size=365)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg8_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg8ON_fitness_maxneigbor_365_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg8_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg8_Fit[sg8_Fit$AACombo %in% c(WT, r$AACombo), ]
}
R<-r$AACombo
r_mat<-stringdistmatrix(r$AACombo, r$AACombo)
R<-R[colSums(r_mat==0)<=1]
r_mat<-stringdistmatrix(R, R)
mean(colSums(r_mat==1))
max(colSums(r_mat==1))
mean(colSums(r_mat==2))
max(colSums(r_mat==2))
write_csv(r, file = paste0("SpCas9_sg8ON_fitness_365", as.character(var_no), reps[d], ".csv"))


### 70% = 510
K<-sample(sg8_Fit$AACombo,510)
K_mat<-stringdistmatrix(K, K)
### count number of neighbor of distance ==1 for random samples
mean(colSums(K_mat==1))
### 4
mean(colSums(K_mat==2))
### 14
max(colSums(K_mat==1))
max(colSums(K_mat==2))
K<-sample(sg8_Fit$AACombo, 510)
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==0)<=1]
K_mat<-stringdistmatrix(K, K)
K<-K[colSums(K_mat==1)<=6& colSums(K_mat==2)<=21]

while(length(K)<510){
  K<-sample(sg8_Fit$AACombo, 510)
  K_mat<-stringdistmatrix(K, K)
  K<-K[colSums(K_mat==0)<=1&colSums(K_mat==1)<=7& colSums(K_mat==2)<=22]
}

K_mat<-stringdistmatrix(K, K)
K<-sample(K, size=510)
if (!WT %in%K){
  K<-K[1:length(K)-1]
  K<-c(WT, K)
}
reps<-c("rep1", "rep2", "rep3")
d=3
write_csv(sg8_Fit %>% filter(AACombo %in%K), file=paste0("SpCas9_sg8ON_fitness_maxneigbor_510_", reps[d], ".csv"))
var_no<-length(K)
r<-sample_n(sg8_Fit, var_no)
if (!WT %in%r$AACombo){
  r<-r[1:nrow(r)-1, ]
  r<-sg8_Fit[sg8_Fit$AACombo %in% c(WT, r$AACombo), ]
}
R<-r$AACombo
r_mat<-stringdistmatrix(r$AACombo, r$AACombo)
R<-R[colSums(r_mat==0)<=1]
r_mat<-stringdistmatrix(R, R)
mean(colSums(K_mat==1))
mean(colSums(K_mat==2))
max(colSums(K_mat==1))
max(colSums(K_mat==2))

mean(colSums(r_mat==1))
max(colSums(r_mat==1))
mean(colSums(r_mat==2))
max(colSums(r_mat==2))
write_csv(r, file = paste0("SpCas9_sg8ON_fitness_510", as.character(var_no), reps[d], ".csv"))

### Double check variants number and make summary table
setwd("./210729_SpCas9_diverse_variants")
fitness_files<-list.files()
K_mat<-stringdistmatrix(sg5_Fit$AACombo, sg5_Fit$AACombo)
variants_info<-as.tibble(table(K_mat)) %>% mutate(samplename="sg5_Fit")
K_mat<-stringdistmatrix(sg8_Fit$AACombo, sg8_Fit$AACombo)
variants_info<-bind_rows(variants_info,
  as.tibble(table(K_mat)) %>% mutate(samplename="sg8_Fit"))

for (i in 1:length(fitness_files)){
  t<-read_csv(fitness_files[i])
  K_mat<-stringdistmatrix(t$AACombo, t$AACombo)
  variants_info<-bind_rows(variants_info,
                           as.tibble(table(K_mat)) %>% mutate(samplename=fitness_files[i]))
}

variants_info<-variants_info %>% spread(K_mat, n)
variants_info<-variants_info %>% mutate(samplesize=`0`, reps=substring(samplename, nchar(samplename) -7, nchar(samplename)-4))
