setwd("/media/achu/Data/Cas9_MLDE")
library(tidyverse)
library(readxl)

muts<-read_xlsx("./41592_2019_473_MOESM4_ESM_CombiSEAL_supp_Fig2.xlsx")
colnames(muts)<-muts[4,]
muts<-muts[5:956, ]
muts[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")]<-sapply(muts[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")], as.numeric)
#muts<-muts[!is.na(muts$"RFPsg5 ON") & !is.na(muts$"RFPsg8 ON"), ]
library(seqinr)
Cas9<-read.fasta(file="./spCas9_uniprot.fasta", 
                 seqtype = c("AA"), as.string = T, seqonly = T)

# header example >GB1_V39A-D40A-E42L
#MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGAAGLWTYDDATKTFTVTE

### make MDLE input fasta and fitness info .csv
output_fasta<-list()
variant_seq<-toString(Cas9[[1]][1])
for (i in 1:nrow(muts)){
  fasta_header<-c("Cas9_", "R661", muts[i, "661"], 
                  "-Q695", muts[i, "695"], 
                  "-K848", muts[i, "848"], 
                  "-E923", muts[i, "923"], 
                  "-T924", muts[i, "924"], 
                  "-Q926", muts[i, "926"], 
                  "-K1003", muts[i, "1003"], 
                  "-R1060", muts[i, "1060"])
  fasta_header<-paste(fasta_header, sep="", collapse = "")
  fasta_variant<-c(substring(variant_seq, 1, 660), muts[i, "661"],
                   substring(variant_seq, 662, 694), muts[i, "695"],
                   substring(variant_seq, 696, 847), muts[i, "848"],
                   substring(variant_seq, 849, 922), muts[i, "923"], muts[i, "924"],
                   substring(variant_seq, 925, 925), muts[i, "926"],
                   substring(variant_seq, 927, 1002), muts[i, "1003"],
                   substring(variant_seq, 1004, 1059), muts[i, "1060"],
                   substring(variant_seq, 1061, 1368))
  fasta_variant<-paste(fasta_variant, sep = "", collapse = "")
  output_fasta[[fasta_header]]<-fasta_variant
}

write.fasta(sequences = output_fasta, names = names(output_fasta), file.out = "Cas9_MLDE_input.fasta", nbchar=1500)

