## Miss Sex all species Plots & Stats ###



#### To make a clean R environment
rm(list=ls())
cat("\014")  

# Packages and functions
#---------------------------------------

dir=dirname(rstudioapi::getActiveDocumentContext()$path)

library(ggplot2)  #ggplot
library(forcats)
library(lme4)
library(ggfortify)
library(cluster)
library(effects)
library(magrittr)
library(lmerTest)
library(ggpubr)
library(factoextra)
library(tidyverse)
library(reshape2)
library(scales)
library(data.table)
library(vegan)
library(devtools)
library(ggConvexHull)
library(ggrepel)       

#####

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=T,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=T) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


# importing data
#---------------------------------------

###############################################################  
###############################################################
###############################
###############################  Small RNA-seq
###############################
###############################################################
###############################################################

###############################################################  
###############################################################
#########
######### Analyse of samples of high quality only (selected sampled, ss, by mgx)
#########
###############################################################  
###############################################################


SSdataFish=read.csv2(paste(dir, "/Analyse_Selection/Selected samples.csv", sep=""), dec=",", header=TRUE, na.strings="na")

######### Analyse Females vs Males all stades and species together

# NMDS On all miRNAs

# Charger le fichier ou il y a tous le nom des  genes:

SSrawdatadesq2MvsF=read.table(paste(dir, "/Analyse_Selection/F_vs_M_species_and_stade_independant/F_vs_M/Geffroy_F_vs_M_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                              header=T, dec=".", sep="\t", fill=T)
rownames(SSrawdatadesq2MvsF)=SSrawdatadesq2MvsF[,1]
SSrawdatadesq2MvsF=SSrawdatadesq2MvsF[,-1]
SSrawdatadesq2MvsFnmds=as.data.frame(t(SSrawdatadesq2MvsF))

SSrawdatadesq2MvsFnmds$ID=rownames(SSrawdatadesq2MvsFnmds)
SSrawdatadesq2MvsFnmds$ID=colsplit(SSrawdatadesq2MvsFnmds$ID, "\\_", c("ID", "Sexe"))[,1]
row.names(SSrawdatadesq2MvsFnmds)=SSrawdatadesq2MvsFnmds$ID


SSdataFish=read.csv2(paste(dir, "/Analyse_Selection/Selected samples.csv", sep=""), dec=",", header=TRUE, na.strings="na")
row.names(SSdataFish)=SSdataFish$ID
all(rownames(SSdataFish) %in% rownames(SSrawdatadesq2MvsFnmds))
SSdataFish <- SSdataFish[rownames(SSrawdatadesq2MvsFnmds),]
all(rownames(SSdataFish) == rownames(SSrawdatadesq2MvsFnmds))

SSrawdatadesq2MvsFnmds1=merge(SSdataFish, SSrawdatadesq2MvsFnmds, by=0, all=TRUE)
SSrawdatadesq2MvsFnmds1$Species[which(SSrawdatadesq2MvsFnmds1$Species=="Jack")]="Blue runner"
SSrawdatadesq2MvsFnmds1$Species[which(SSrawdatadesq2MvsFnmds1$Species=="Sea_bass")]="Sea bass"
SSrawdatadesq2MvsFnmds1$Species[which(SSrawdatadesq2MvsFnmds1$Species=="Red_drum")]="Red drum"
#SSrawdatadesq2MvsFv7PCA=subset(SSrawdatadesq2MvsFv7, select=-c(Row.names, ID.y, ID.x, Sex, Species, Species_VF, Stade))

my_sample_col=SSrawdatadesq2MvsFnmds1[,c("Species", "Sex", "Stade")]
com = SSrawdatadesq2MvsFnmds1[,7:ncol(SSrawdatadesq2MvsFnmds1)]
com=subset(com, select=-c(ID.y))


m_com = as.matrix(com)

set.seed(123)
nmds1 = metaMDS(m_com, distance = "bray")
nmds1

en = envfit(nmds1, my_sample_col, permutations = 999, na.rm = TRUE)
en

plot(nmds1)
plot(en)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds1)$sites)

#add columns to data frame 
data.scores$Sex = SSrawdatadesq2MvsFnmds1$Sex
data.scores$Species = SSrawdatadesq2MvsFnmds1$Species

en_coord_cat = as.data.frame(scores(en, "factors"))
NameAni=c("Blue runner", "Red drum", "Sea bass", "Turbot")
en_coord_cat=en_coord_cat[-c(5:8),]
en_coord_cat=cbind(en_coord_cat, NameAni)

head(data.scores)


ggsexspecies = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Sex, shape=Species), size = 4, alpha = 10) + 
  scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))  + 
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2, shape=NameAni), 
             size = 6, alpha = 0.6, colour = "black") + theme_bw() 

ggsexspecies

ggsave(ggsexspecies, file="Figure 1.pdf", width = 7, height=6, dpi=300)


#########
######### Analyse Females vs Males only immatures all species together
#########


# Charger le fichier ou il y a tous le nom des  genes:

SSrawdatadesq2MvsFImm=read.table(paste(dir, "/Analyse_Selection/Immature/F_vs_M/Geffroy_F_vs_M_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                              header=T, dec=".", sep="\t", fill=T)
rownames(SSrawdatadesq2MvsFImm)=SSrawdatadesq2MvsFImm[,1]

# Keep only significative ones (UNcorrected p-value)

SSsignidesq2MvsFImm=read.csv2(paste(dir, "/Analyse_Selection/Immature/F_vs_M/Geffroy_DESeq2_F_vs_M_ALLp-value.csv", sep=""), dec=",", header=TRUE, row.names = 1)
tokeepsigniImm2=SSsignidesq2MvsFImm[SSsignidesq2MvsFImm$pvalue<0.05,]

all(rownames(SSrawdatadesq2MvsFImm) %in% rownames(tokeepsigniImm2))
SSsignidesq2MvsFImmv0 <- SSrawdatadesq2MvsFImm[rownames(tokeepsigniImm2),]
SSsignidesq2MvsFImmv0=SSsignidesq2MvsFImmv0[,-1]
SSrawdatadesq2MvsFv1Imm=as.data.frame(t(SSsignidesq2MvsFImmv0))

SSrawdatadesq2MvsFv1Imm$ID=rownames(SSrawdatadesq2MvsFv1Imm)
SSrawdatadesq2MvsFv1Imm$ID=colsplit(SSrawdatadesq2MvsFv1Imm$ID, "\\_", c("ID", "Sexe"))[,1]
row.names(SSrawdatadesq2MvsFv1Imm)=SSrawdatadesq2MvsFv1Imm$ID

row.names(SSdataFish)=SSdataFish$ID
all(rownames(SSdataFish) %in% rownames(SSrawdatadesq2MvsFv1Imm))
SSdataFish <- SSdataFish[rownames(SSrawdatadesq2MvsFv1Imm),]
all(rownames(SSdataFish) == rownames(SSrawdatadesq2MvsFv1Imm))

SSrawdatadesq2MvsFv2Imm=merge(SSdataFish, SSrawdatadesq2MvsFv1Imm, by=0, all=TRUE)
SSrawdatadesq2MvsFv2ImmPCA=subset(SSrawdatadesq2MvsFv2Imm, select=-c(Row.names, ID.y, ID.x, Sex, Species, Species_VF, Stade))

names(SSrawdatadesq2MvsFv2ImmPCA)[names(SSrawdatadesq2MvsFv2ImmPCA) == "REV_miR-18a-3p"] <- "miR-18a-3p"

SSrawdatadesq2MvsFv2ImmPCAlog=log(SSrawdatadesq2MvsFv2ImmPCA+1)
summary(SSrawdatadesq2MvsFv2ImmPCAlog)

my.col.var<-rep(c("red","black"),times=c(3,16))
my.col.var1<- c("red"  , "red"  , "red" ,  "black", "black", "black" ,"black" ,"black", "black" ,"black", "black", "black" ,"black" ,"black", "black" ,"black" ,"black" ,"black" ,"black")


SSrawdatadesq2MvsFv2Imm$Species[which(SSrawdatadesq2MvsFv2Imm$Species=="Sea_bass")]="Sea bass"
SSrawdatadesq2MvsFv2Imm$Species[which(SSrawdatadesq2MvsFv2Imm$Species=="Red_drum")]="Red drum"


pImmMvsFemalesigni=autoplot(princomp(SSrawdatadesq2MvsFv2ImmPCAlog, cor = TRUE,scores = TRUE), 
                            data = SSrawdatadesq2MvsFv2Imm, colour = 'Sex',size = 3, frame = F, shape = "Species", loadings = FALSE, axes = c(1, 2),
                            loadings.colour = my.col.var1,loadings.label = TRUE, loadings.label.repel= TRUE, loadings.label.colour=my.col.var1) + scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()+
geom_convexhull(aes(fill= SSrawdatadesq2MvsFv2Imm$Sex), alpha = 0.3)+scale_fill_manual(values=c('#EFC000FF',"darkgoldenrod4"))+
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2), fill=FALSE)
pImmMvsFemalesigni

ggsave(pImmMvsFemalesigni, file="Figure 2a.pdf", width = 8, height=8, dpi=300)

#########
######### Analyse Females vs Males only matures all species together
#########

SSdataFish=read.csv2(paste(dir, "/Analyse_Selection/Selected samples.csv", sep=""), dec=",", header=TRUE, na.strings="na")

SSrawdatadesq2MvsFMature=read.table(paste(dir, "/Analyse_Selection/Mature/F_vs_M/Geffroy_F_vs_M_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                                 header=T, dec=".", sep="\t", fill=T)
rownames(SSrawdatadesq2MvsFMature)=SSrawdatadesq2MvsFMature[,1]
SSrawdatadesq2MvsFMature=SSrawdatadesq2MvsFMature[,-1]

# Keep only significative ones (UNcorrected p-value)

SSsignidesq2MvsFMature=read.csv2(paste(dir, "/Analyse_Selection/Mature/F_vs_M/Geffroy_DESeq2_F_vs_M_ALLp-value.csv", sep=""), dec=",", header=TRUE, row.names = 1)
tokeepsigniMature=SSsignidesq2MvsFMature[SSsignidesq2MvsFMature$pvalue<0.05,]

all(rownames(SSrawdatadesq2MvsFMature) %in% rownames(tokeepsigniMature))
SSsignidesq2MvsFMaturev0 <- SSrawdatadesq2MvsFMature[rownames(tokeepsigniMature),]
#SSsignidesq2MvsFMaturev0=SSsignidesq2MvsFMaturev0[,-1]
SSrawdatadesq2MvsFv1Mature=as.data.frame(t(SSsignidesq2MvsFMaturev0))

SSrawdatadesq2MvsFv1Mature$ID=rownames(SSrawdatadesq2MvsFv1Mature)
SSrawdatadesq2MvsFv1Mature$ID=colsplit(SSrawdatadesq2MvsFv1Mature$ID, "\\_", c("ID", "Sexe"))[,1]
row.names(SSrawdatadesq2MvsFv1Mature)=SSrawdatadesq2MvsFv1Mature$ID

row.names(SSdataFish)=SSdataFish$ID
all(rownames(SSdataFish) %in% rownames(SSrawdatadesq2MvsFv1Mature))
SSdataFish1 <- SSdataFish[rownames(SSrawdatadesq2MvsFv1Mature),]
all(rownames(SSdataFish1) == rownames(SSrawdatadesq2MvsFv1Mature))

SSrawdatadesq2MvsFv2Mature=merge(SSdataFish1, SSrawdatadesq2MvsFv1Mature, by=0, all=TRUE)
SSrawdatadesq2MvsFv2MaturePCA=subset(SSrawdatadesq2MvsFv2Mature, select=-c(Row.names, ID.y, ID.x, Sex, Species, Species_VF, Stade))

SSrawdatadesq2MvsFv2MaturePCAlog=log(SSrawdatadesq2MvsFv2MaturePCA+1)
summary(SSrawdatadesq2MvsFv2MaturePCAlog)

SSrawdatadesq2MvsFv2Mature$Species[which(SSrawdatadesq2MvsFv2Mature$Species=="Sea_bass")]="Sea bass"
SSrawdatadesq2MvsFv2Mature$Species[which(SSrawdatadesq2MvsFv2Mature$Species=="Jack")]="Blue runner"

pMatureMvsFemalesigni=autoplot(princomp(SSrawdatadesq2MvsFv2MaturePCAlog, cor = TRUE,scores = TRUE), 
                               data = SSrawdatadesq2MvsFv2Mature, colour = 'Sex',size = 3, frame = F, shape = "Species", loadings = FALSE, axes = c(1, 2),
                               loadings.colour = 'black',loadings.label = TRUE, loadings.label.repel= TRUE, loadings.label.colour="black") + 
  scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()+ scale_x_reverse()+
  geom_convexhull(aes(fill= SSrawdatadesq2MvsFv2Mature$Sex), alpha = 0.3)+scale_fill_manual(values=c('#EFC000FF',"darkgoldenrod4"))+
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2), fill=FALSE) 

pMatureMvsFemalesigni

ggsave(pMatureMvsFemalesigni, file="Figure 2b.pdf", width = 8, height=8, dpi=300)


Figure2=ggarrange(pImmMvsFemalesigni, pMatureMvsFemalesigni,
                   labels = c("A", "B"),
                   ncol = 1, nrow = 2)
Figure2


ggsave(Figure2, file="Figure 2.pdf", width = 7, height=9, dpi=300)

#########
######### Analyse Females vs Males per species independent of stage 
#########


SSdataFish=read.csv2(paste(dir, "/Analyse_Selection/Selected samples.csv", sep=""), dec=",", header=TRUE, na.strings="na")

########## Blue Runner

SSrawdatadesq2MvsFJacks=read.table(paste(dir, "/Analyse_Selection/F_vs_M_by_species_stade_independant/Carangues/F_vs_M/Geffroy_F_vs_M_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                                    header=T, dec=".", sep="\t", fill=T)
rownames(SSrawdatadesq2MvsFJacks)=SSrawdatadesq2MvsFJacks[,1]

# Keep only significative ones (UNcorrected p-value)

SSsignidesq2MvsFJacks=read.csv2(paste(dir, "/Analyse_Selection/F_vs_M_by_species_stade_independant/Carangues/F_vs_M/Geffroy_DESeq2_F_vs_M_ALLp-value.csv", sep=""), dec=",", header=TRUE, row.names = 1)
tokeepsigniJacks=SSsignidesq2MvsFJacks[SSsignidesq2MvsFJacks$padj<0.05,]

all(rownames(SSrawdatadesq2MvsFJacks) %in% rownames(tokeepsigniJacks))
SSsignidesq2MvsFJacksv0 <- SSrawdatadesq2MvsFJacks[rownames(tokeepsigniJacks),]
SSsignidesq2MvsFJacksv0=SSsignidesq2MvsFJacksv0[,-1]
SSrawdatadesq2MvsFv1Jacks=as.data.frame(t(SSsignidesq2MvsFJacksv0))

SSrawdatadesq2MvsFv1Jacks$ID=rownames(SSrawdatadesq2MvsFv1Jacks)
SSrawdatadesq2MvsFv1Jacks$ID=colsplit(SSrawdatadesq2MvsFv1Jacks$ID, "\\_", c("ID", "Sexe"))[,1]
row.names(SSrawdatadesq2MvsFv1Jacks)=SSrawdatadesq2MvsFv1Jacks$ID

row.names(SSdataFish)=SSdataFish$ID
all(rownames(SSdataFish) %in% rownames(SSrawdatadesq2MvsFv1Jacks))
SSdataFish1 <- SSdataFish[rownames(SSrawdatadesq2MvsFv1Jacks),]
all(rownames(SSdataFish1) == rownames(SSrawdatadesq2MvsFv1Jacks))

SSrawdatadesq2MvsFv2Jacks=merge(SSdataFish1, SSrawdatadesq2MvsFv1Jacks, by=0, all=TRUE)
SSrawdatadesq2MvsFv2Jacks=SSrawdatadesq2MvsFv2Jacks[,-c(10:22)]
SSrawdatadesq2MvsFv2JacksPCA=subset(SSrawdatadesq2MvsFv2Jacks, select=-c(Row.names, ID.y, ID.x, Sex, Species, Species_VF, Stade))

SSrawdatadesq2MvsFv2JacksPCAlog=log(SSrawdatadesq2MvsFv2JacksPCA+1)
summary(SSrawdatadesq2MvsFv2JacksPCAlog)


pJacksMvsFemalesigni=autoplot(princomp(SSrawdatadesq2MvsFv2JacksPCAlog, cor = TRUE,scores = TRUE), 
                               data = SSrawdatadesq2MvsFv2Jacks, colour = 'Sex',size = 3, frame = F, loadings = FALSE, axes = c(1, 2),
                               loadings.colour = 'black',loadings.label = TRUE, loadings.label.vjust = -0.4, loadings.label.colour="black") + 
  scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()

pJacksMvsFemalesigni
ggsave(pJacksMvsFemalesigni, file="JacksMvsFsigniUncorrectedp-value.pdf", width = 8, height=8, dpi=300)



######## Red Drum

SSrawdatadesq2MvsFReddrum=read.table(paste(dir, "/Analyse_Selection/F_vs_M_by_species_stade_independant/Ombrine/F_vs_M/Geffroy_F_vs_M_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                                   header=T, dec=".", sep="\t", fill=T)
rownames(SSrawdatadesq2MvsFReddrum)=SSrawdatadesq2MvsFReddrum[,1]

# Keep only significative ones (UNcorrected p-value)

SSsignidesq2MvsFReddrum=read.csv2(paste(dir, "/Analyse_Selection/F_vs_M_by_species_stade_independant/Ombrine/F_vs_M/Geffroy_DESeq2_F_vs_M_ALLp-value.csv", sep=""), dec=",", header=TRUE, row.names = 1)
tokeepsigniReddrum=SSsignidesq2MvsFReddrum[SSsignidesq2MvsFReddrum$pvalue<0.05,]
tokeepsigniReddrum=na.omit(tokeepsigniReddrum)


all(rownames(SSrawdatadesq2MvsFReddrum) %in% rownames(tokeepsigniReddrum))
SSsignidesq2MvsFReddrumv0 <- SSrawdatadesq2MvsFReddrum[rownames(tokeepsigniReddrum),]
SSsignidesq2MvsFReddrumv0=SSsignidesq2MvsFReddrumv0[,-1]
SSrawdatadesq2MvsFv1Reddrum=as.data.frame(t(SSsignidesq2MvsFReddrumv0))

SSrawdatadesq2MvsFv1Reddrum$ID=rownames(SSrawdatadesq2MvsFv1Reddrum)
SSrawdatadesq2MvsFv1Reddrum$ID=colsplit(SSrawdatadesq2MvsFv1Reddrum$ID, "\\_", c("ID", "Sexe"))[,1]
row.names(SSrawdatadesq2MvsFv1Reddrum)=SSrawdatadesq2MvsFv1Reddrum$ID

row.names(SSdataFish)=SSdataFish$ID
all(rownames(SSdataFish) %in% rownames(SSrawdatadesq2MvsFv1Reddrum))
SSdataFish1 <- SSdataFish[rownames(SSrawdatadesq2MvsFv1Reddrum),]
all(rownames(SSdataFish1) == rownames(SSrawdatadesq2MvsFv1Reddrum))

SSrawdatadesq2MvsFv2Reddrum=merge(SSdataFish1, SSrawdatadesq2MvsFv1Reddrum, by=0, all=TRUE)
SSrawdatadesq2MvsFv2ReddrumPCA=subset(SSrawdatadesq2MvsFv2Reddrum, select=-c(Row.names, ID.y, ID.x, Sex, Species, Species_VF, Stade))

SSrawdatadesq2MvsFv2ReddrumPCAlog=log(SSrawdatadesq2MvsFv2ReddrumPCA+1)
summary(SSrawdatadesq2MvsFv2ReddrumPCAlog)


pReddrumMvsFemalesigni=autoplot(princomp(SSrawdatadesq2MvsFv2ReddrumPCAlog, cor = TRUE,scores = TRUE), 
                              data = SSrawdatadesq2MvsFv2Reddrum, colour = 'Sex',size = 3, frame = F, loadings = FALSE, axes = c(1, 2),
                              loadings.colour = 'black',loadings.label = TRUE, loadings.label.vjust = -0.4, loadings.label.colour="black") + 
  scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()

pReddrumMvsFemalesigni
ggsave(pReddrumMvsFemalesigni, file="ReddrumMvsFsigniUncorrectedp-value.pdf", width = 8, height=8, dpi=300)

######## Turbot

SSrawdatadesq2MvsFTurbot=read.table(paste(dir, "/Analyse_Selection/F_vs_M_by_species_stade_independant/Turbot/F_vs_M/Geffroy_F_vs_M_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                                     header=T, dec=".", sep="\t", fill=T)
rownames(SSrawdatadesq2MvsFTurbot)=SSrawdatadesq2MvsFTurbot[,1]

# Keep only significative ones (UNcorrected p-value)

SSsignidesq2MvsFTurbot=read.csv2(paste(dir, "/Analyse_Selection/F_vs_M_by_species_stade_independant/Turbot/F_vs_M/Geffroy_DESeq2_F_vs_M_ALLp-value.csv", sep=""), dec=",", header=TRUE, row.names = 1)
tokeepsigniTurbot=SSsignidesq2MvsFTurbot[SSsignidesq2MvsFTurbot$pvalue<0.05,]
tokeepsigniTurbot=na.omit(tokeepsigniTurbot)


all(rownames(SSrawdatadesq2MvsFTurbot) %in% rownames(tokeepsigniTurbot))
SSsignidesq2MvsFTurbotv0 <- SSrawdatadesq2MvsFTurbot[rownames(tokeepsigniTurbot),]
SSsignidesq2MvsFTurbotv0=SSsignidesq2MvsFTurbotv0[,-1]
SSrawdatadesq2MvsFv1Turbot=as.data.frame(t(SSsignidesq2MvsFTurbotv0))

SSrawdatadesq2MvsFv1Turbot$ID=rownames(SSrawdatadesq2MvsFv1Turbot)
SSrawdatadesq2MvsFv1Turbot$ID=colsplit(SSrawdatadesq2MvsFv1Turbot$ID, "\\_", c("ID", "Sexe"))[,1]
row.names(SSrawdatadesq2MvsFv1Turbot)=SSrawdatadesq2MvsFv1Turbot$ID

row.names(SSdataFish)=SSdataFish$ID
all(rownames(SSdataFish) %in% rownames(SSrawdatadesq2MvsFv1Turbot))
SSdataFish1 <- SSdataFish[rownames(SSrawdatadesq2MvsFv1Turbot),]
all(rownames(SSdataFish1) == rownames(SSrawdatadesq2MvsFv1Turbot))

SSrawdatadesq2MvsFv2Turbot=merge(SSdataFish1, SSrawdatadesq2MvsFv1Turbot, by=0, all=TRUE)
SSrawdatadesq2MvsFv2TurbotPCA=subset(SSrawdatadesq2MvsFv2Turbot, select=-c(Row.names, ID.y, ID.x, Sex, Species, Species_VF, Stade))

SSrawdatadesq2MvsFv2TurbotPCAlog=log(SSrawdatadesq2MvsFv2TurbotPCA+1)
summary(SSrawdatadesq2MvsFv2TurbotPCAlog)


pTurbotMvsFemalesigni=autoplot(princomp(SSrawdatadesq2MvsFv2TurbotPCAlog, cor = TRUE,scores = TRUE), 
                                data = SSrawdatadesq2MvsFv2Turbot, colour = 'Sex',size = 3, frame = F, loadings = FALSE, axes = c(1, 2),
                                loadings.colour = 'black',loadings.label = TRUE, loadings.label.vjust = -0.4, loadings.label.colour="black") + 
  scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()

pTurbotMvsFemalesigni
ggsave(pTurbotMvsFemalesigni, file="TurbotMvsFsigniUncorrectedp-value.pdf", width = 8, height=8, dpi=300)


###############################################################  
###############################################################
#########
######### Analyse of Seabream stages
#########
###############################################################  
###############################################################


SBdataFishSexchange=read.csv2(paste(dir, "/Daurade_Gac/All Info Seabream.csv", sep=""), dec=",", header=TRUE, na.strings="na")

######### Analyse Females vs Males all stades and species together

# Charger le fichier ou il y a tous le nom des  genes:

SBrawdatadesqSC=read.table(paste(dir, "/Daurade_Gac/ALL/Geffroy_Daurade_Gac_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                              header=T, dec=".", sep="\t", fill=T)
rownames(SBrawdatadesqSC)=SBrawdatadesqSC[,1]

# Comparing PM vs GF

SBsignidesq2PMvsGF=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/P_M_vs_G_F/miSS_Daurade_Gac_miRNA_DESeq2_P_M_vs_G_F_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2PMvsGF=SBsignidesq2PMvsGF[SBsignidesq2PMvsGF$padj<0.005,]

all(rownames(SBrawdatadesqSC) %in% rownames(SBsignidesq2PMvsGF))
SBsignidesq2PMvsGFv0 <- SBrawdatadesqSC[rownames(SBsignidesq2PMvsGF),]
SBsignidesq2PMvsGFv0=SBsignidesq2PMvsGFv0[,-1]
SBsignidesq2PMvsGFv1=as.data.frame(t(SBsignidesq2PMvsGFv0))

row.names(SBdataFishSexchange)=SBdataFishSexchange$Name.sequencing
SBdataFishSexchange=SBdataFishSexchange[,c(4,5,7)]
all(rownames(SBdataFishSexchange) %in% rownames(SBsignidesq2PMvsGFv1))
SBdataFishSexchangok <- SBdataFishSexchange[rownames(SBsignidesq2PMvsGFv1),]
all(rownames(SBdataFishSexchangok) == rownames(SBsignidesq2PMvsGFv1))

SBsignidesq2PMvsGFv2=merge(SBsignidesq2PMvsGFv1, SBdataFishSexchangok, by=0, all=TRUE)

SBsignidesq2PMvsGFv3=SBsignidesq2PMvsGFv2[which(SBsignidesq2PMvsGFv2$sex.name.Sequencing %in%  c("P_M","G_F")),]
SBsignidesq2PMvsGFv3PCA=subset(SBsignidesq2PMvsGFv3, select=-c(sex.name.Sequencing, Poids, Taille, Row.names))

SBsignidesq2PMvsGFv3PCAlog=log(SBsignidesq2PMvsGFv3PCA+1)
summary(SBsignidesq2PMvsGFv3PCAlog, getOption("max.print"))


pAllPMvsGF=autoplot(princomp(SBsignidesq2PMvsGFv3PCAlog, cor = TRUE,scores = TRUE), 
                            data = SBsignidesq2PMvsGFv3, colour = 'sex.name.Sequencing',size = 3, frame = F, loadings = FALSE, axes = c(1, 2),
                            loadings.colour = 'black',loadings.label = TRUE, loadings.label.vjust = -0.4, loadings.label.colour="black") + 
  scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()

pAllPMvsGF
ggsave(pAllPMvsGF, file="AllMvsFsigni.pdf", width = 8, height=8, dpi=300)



######### Comparing all stages (all minus M_MPP_vs_M_I; M_MP_vs_M_I)

SBdataFishSexchange=read.csv2(paste(dir, "/Daurade_Gac/All Info Seabream.csv", sep=""), dec=",", header=TRUE, na.strings="na")

######### Analyse Females vs Males all stades and species together

# Charger le fichier ou il y a tous le nom des  genes:

SBrawdatadesqSC=read.table(paste(dir, "/Daurade_Gac/ALL/Geffroy_Daurade_Gac_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                           header=T, dec=".", sep="\t", fill=T)
rownames(SBrawdatadesqSC)=SBrawdatadesqSC[,1]

SBsignidesq2IvsGF=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/M_I_vs_G_F/miSS_Daurade_Gac_miRNA_DESeq2_M_I_vs_G_F_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2PMvsGF=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/P_M_vs_G_F/miSS_Daurade_Gac_miRNA_DESeq2_P_M_vs_G_F_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2PMvsMP=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/P_M_vs_M_MP/miSS_Daurade_Gac_miRNA_DESeq2_P_M_vs_M_MP_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2PMvsI=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/P_M_vs_M_I/miSS_Daurade_Gac_miRNA_DESeq2_P_M_vs_M_I_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2PMvsMPP=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/P_M_vs_M_MPP/miSS_Daurade_Gac_miRNA_DESeq2_P_M_vs_M_MPP_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2MPPvsGF=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/M_MPP_vs_G_F/miSS_Daurade_Gac_miRNA_DESeq2_M_MPP_vs_G_F_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2MPvsMPP=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/M_MP_vs_M_MPP/miSS_Daurade_Gac_miRNA_DESeq2_M_MP_vs_M_MPP_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)
SBsignidesq2MPvsGF=read.csv2(paste(dir, "/Daurade_Gac/comparaisons/M_MP_vs_G_F/miSS_Daurade_Gac_miRNA_DESeq2_M_MP_vs_G_F_p-value0_05.csv", sep=""), dec=",", header=TRUE, row.names = 1)

SBsignidesq2IvsGF$ID=rownames(SBsignidesq2IvsGF)
SBsignidesq2PMvsGF$ID=rownames(SBsignidesq2PMvsGF)
SBsignidesq2PMvsMP$ID=rownames(SBsignidesq2PMvsMP)
SBsignidesq2PMvsI$ID=rownames(SBsignidesq2PMvsI)
SBsignidesq2PMvsMPP$ID=rownames(SBsignidesq2PMvsMPP)
SBsignidesq2MPPvsGF$ID=rownames(SBsignidesq2MPPvsGF)
SBsignidesq2MPvsMPP$ID=rownames(SBsignidesq2MPvsMPP)
SBsignidesq2MPvsGF$ID=rownames(SBsignidesq2MPvsGF)


Allsignif=bind_rows(SBsignidesq2IvsGF,SBsignidesq2PMvsGF, SBsignidesq2PMvsMP,SBsignidesq2PMvsI,SBsignidesq2PMvsMPP,SBsignidesq2MPPvsGF,SBsignidesq2MPvsMPP, SBsignidesq2MPvsGF)
Allsignif<- Allsignif[!duplicated(Allsignif$ID), ]
rownames(Allsignif)=Allsignif$ID

write.csv2(Allsignif, "allsigniseabream.csv",row.names = T)

Allsignif=Allsignif[Allsignif$padj<0.01,]

all(rownames(SBrawdatadesqSC) %in% rownames(Allsignif))
SBsignidesq2ALLv0 <- SBrawdatadesqSC[rownames(Allsignif),]
SBsignidesq2ALLv0=SBsignidesq2ALLv0[,-1]
SBsignidesq2Allv1=as.data.frame(t(SBsignidesq2ALLv0))

row.names(SBdataFishSexchange)=SBdataFishSexchange$Name.sequencing
SBdataFishSexchange=SBdataFishSexchange[,c(4,5,7)]
all(rownames(SBdataFishSexchange) %in% rownames(SBsignidesq2Allv1))
SBdataFishSexchangok <- SBdataFishSexchange[rownames(SBsignidesq2Allv1),]
all(rownames(SBdataFishSexchangok) == rownames(SBsignidesq2Allv1))

SBsignidesq2ALLv2=merge(SBsignidesq2Allv1, SBdataFishSexchangok, by=0, all=TRUE)
SBsignidesq2ALLv3PCA=subset(SBsignidesq2ALLv2, select=-c(sex.name.Sequencing, Poids, Taille, Row.names))

SBsignidesq2ALLv3PCAlog=log(SBsignidesq2ALLv3PCA+1)
summary(SBsignidesq2ALLv3PCAlog, getOption("max.print"))


pAll=autoplot(princomp(SBsignidesq2ALLv3PCAlog, cor = TRUE,scores = TRUE), 
                    data = SBsignidesq2ALLv2, colour = 'sex.name.Sequencing',size = 3, frame = F, loadings = FALSE, axes = c(1, 2),
                    loadings.colour = 'black',loadings.label = TRUE, loadings.label.vjust = -0.4, loadings.label.colour="black") #+ 
  #scale_colour_manual(values=c('#EFC000FF',"darkgoldenrod4"))+theme_bw()

pAll
ggsave(pAllMvsFemalesigni, file="AllMvsFsigni.pdf", width = 8, height=8, dpi=300)

SBsignidesq2ALLv2$sex = factor(SBsignidesq2ALLv2$sex.name.Sequencing, levels=c('P_M','M_MPP','M_MP','M_I','G_F'))

pcasexDaurade=princomp(SBsignidesq2ALLv3PCAlog, cor = TRUE,scores = TRUE)
summary(pcasexDaurade)
p1=fviz_pca_ind(pcasexDaurade, 
                   col.ind = SBsignidesq2ALLv2$sex, geom="point",pointsize = 3,
                   palette = c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"),
                   addEllipses = TRUE, # Concentration ellipses
                   ellipse.type = "confidence",loadings.label.vjust = -0.4,repel = TRUE,
                   axes = c(1, 2) 
)
pEllipsesSexDaurade=p1+ theme_bw()
pEllipsesSexDaurade

ggsave(pEllipsesSexDaurade, file="SeabreamPCAellipsesallstages.pdf", width = 8, height=8, dpi=300)


# Finding common miRNAs detected in gonchoristic and hermaphrodites : i.e. related to sex.


miRNAseabream=as.data.frame(sub("*gac-", "", Allsignif$ID))
rownames(miRNAseabream)=miRNAseabream[,1]
miRNAseabream <- data.frame(do.call('rbind', strsplit(as.character(rownames(miRNAseabream)),'_',fixed=TRUE))) # seabream
miRNAseabream<- miRNAseabream[!duplicated(miRNAseabream$X2), ]
rownames(miRNAseabream)=miRNAseabream$X2

tokeepsigniALL # Analyse Females vs Males all stades and species together
tokeepsigniImm2 # Analyse Females vs Males all species together Immatures
tokeepsigniMature # nalyse Females vs Males all species together matures

data_common2 <- inner_join(miRNAseabream, tokeepsigniALL) 

Test1=miRNAseabream[which(rownames(miRNAseabream) %in% rownames(tokeepsigniALL)),] 
Test2=miRNAseabream[which(rownames(miRNAseabream) %in% rownames(tokeepsigniImm2)),] # miR-125b-2/3-5p in MPPvsGF and miR-21a-3p
Test3=miRNAseabream[which(rownames(miRNAseabream) %in% rownames(tokeepsigniMature)),] # miR-125b-2/3-5p in MPPvsGF (more in females in seabream and matures) and miR-21a-3p (more in males in matures and seabream)

#### NMDS Daurade

######### Comparing all stages (all minus M_MPP_vs_M_I; M_MP_vs_M_I)

SBdataFishSexchange=read.csv2(paste(dir, "/Daurade_Gac/All Info Seabream.csv", sep=""), dec=",", header=TRUE, na.strings="na")

######### Analyse Females vs Males all stades and species together

# Charger le fichier ou il y a tous le nom des  genes:

SBrawdatadesqSC=read.table(paste(dir, "/Daurade_Gac/ALL/Geffroy_Daurade_Gac_DESeq2_Normalized_Counts.txt", sep=""),             #    with condition details
                           header=T, dec=".", sep="\t", fill=T)
rownames(SBrawdatadesqSC)=SBrawdatadesqSC[,1]
SBrawdatadesqSC=SBrawdatadesqSC[,-1]
SBrawdatadesqSC=as.data.frame(t(SBrawdatadesqSC))


row.names(SBdataFishSexchange)=SBdataFishSexchange$Name.sequencing
SBdataFishSexchange=SBdataFishSexchange[,c(4,5,7)]
all(rownames(SBdataFishSexchange) %in% rownames(SBrawdatadesqSC))
SBdataFishSexchangok <- SBdataFishSexchange[rownames(SBrawdatadesqSC),]
all(rownames(SBdataFishSexchangok) == rownames(SBrawdatadesqSC))

SBALLv2=merge(SBrawdatadesqSC, SBdataFishSexchangok, by=0, all=TRUE)


my_sample_colSB=as.data.frame(SBdataFishSexchangok[,3])

com1 = SBrawdatadesqSC
m_com1 = as.matrix(com1)

set.seed(123)
nmds2 = metaMDS(m_com1, distance = "bray")
nmds2

en = envfit(nmds2, my_sample_colSB, permutations = 999, na.rm = TRUE)
en

plot(nmds2)
plot(en)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds2)$sites)

#add columns to data frame 
data.scores$Stages = SBALLv2$sex.name.Sequencing


en_coord_cat = as.data.frame(scores(en, "factors"))
NameStage=c("Female Mature", "Intersex 2", "Intersex 1", "Male Mature", "Immature Male")
en_coord_cat1=cbind(en_coord_cat, NameStage)


head(data.scores)

data.scores$Stages[which(data.scores$Stages=="G_F")]="Female Mature"
data.scores$Stages[which(data.scores$Stages=="M_I")]="Intersex 2"
data.scores$Stages[which(data.scores$Stages=="M_MP")]="Intersex 1"
data.scores$Stages[which(data.scores$Stages=="M_MPP")]="Male Mature"
data.scores$Stages[which(data.scores$Stages=="P_M")]="Immature Male"
data.scores$Stages = factor(data.scores$Stages, levels=c("Immature Male","Male Mature", "Intersex 1", "Intersex 2", "Female Mature"))


ggsexSBreamAllmiRNAS = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Stages), size = 4, alpha = 10) + 
  scale_colour_manual(values=c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"))  + 
  geom_point(data = en_coord_cat1, aes(x = NMDS1, y = NMDS2, colour = NameStage), 
             size = 6) + theme_bw()

ggsexSBreamAllmiRNAS

ggsave(ggsexSBreamAllmiRNAS, file="Figure 3B.pdf", width = 10, height=8, dpi=300)

#### NMDS Daurade only significant miRNAs

SBsignidesq2ALLv2=merge(SBsignidesq2Allv1, SBdataFishSexchangok, by=0, all=TRUE)
SBsignidesq2ALLv3PCA=subset(SBsignidesq2ALLv2, select=-c(sex.name.Sequencing, Poids, Taille, Row.names))

my_sample_colSB2=as.data.frame(SBsignidesq2ALLv2[,37])

com2 = SBsignidesq2ALLv3PCA
m_com2 = as.matrix(com2)

set.seed(123)
nmds3 = metaMDS(m_com2, distance = "bray")
nmds3

en1 = envfit(nmds3, my_sample_colSB2, permutations = 999, na.rm = TRUE)
en1

plot(nmds3)
plot(en1)

#extract NMDS scores (x and y coordinates)
data.scores2 = as.data.frame(scores(nmds3)$sites)

#add columns to data frame 
data.scores2$Stages = SBsignidesq2ALLv2$sex.name.Sequencing


en_coord_cat2 = as.data.frame(scores(en1, "factors"))
NameStage=c("Female Mature", "Intersex 2","Intersex 1","Male Mature","Immature Male")
en_coord_cat3=cbind(en_coord_cat2, NameStage)



head(data.scores2)

data.scores2$Stages[which(data.scores2$Stages=="G_F")]="Female Mature"
data.scores2$Stages[which(data.scores2$Stages=="M_I")]="Intersex 2"
data.scores2$Stages[which(data.scores2$Stages=="M_MP")]="Intersex 1"
data.scores2$Stages[which(data.scores2$Stages=="M_MPP")]="Male Mature"
data.scores2$Stages[which(data.scores2$Stages=="P_M")]="Immature Male"
data.scores2$Stages = factor(data.scores2$Stages, levels=c("Immature Male","Male Mature", "Intersex 1", "Intersex 2", "Female Mature"))


#Test with ellipses
# ggsexSBreamSignimiRNAS = ggplot(data = data.scores2, aes(x = NMDS1, y = NMDS2)) + 
#   geom_point(data = data.scores2, aes(colour = Stages), size = 4, alpha = 10) + stat_ellipse(data = data.scores2, aes(colour = Stages, fill=Stages), alpha = 0.2, geom = "polygon")+
#   scale_colour_manual(values=c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1")) + 
#   scale_fill_manual(values=c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"))+
#   geom_point(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2, colour = NameStage), 
#              size = 6) + theme_bw()
# 
# ggsexSBreamSignimiRNAS

ggsexSBreamSignimiRNAS = ggplot(data = data.scores2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores2, aes(colour = Stages), size = 4, alpha = 10) +
  scale_colour_manual(values=c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1")) + 
  geom_point(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2, colour = NameStage), 
             size = 6) + theme_bw()

ggsexSBreamSignimiRNAS


ggsexSBreamAllmiRNAS
ggsexSBreamSignimiRNAS


FigureSbreamALL2=ggarrange(ggsexSBreamAllmiRNAS, ggsexSBreamSignimiRNAS, 
                          labels = c("B", "C"),
                          ncol = 2, nrow=1,
                          common.legend = TRUE)
FigureSbreamALL2

ggsave(FigureSbreamALL2, file="Figure 3.pdf", width = 6, height=4, dpi=300)


########## miR-125b-2/3-5p

SBdataFishSexchange=read.csv2(paste(dir, "/Daurade_Gac/All Info Seabream.csv", sep=""), dec=",", header=TRUE, na.strings="na")
row.names(SBdataFishSexchange)=SBdataFishSexchange$Name.sequencing
SBmir125b5p=read.csv2(paste(dir, "/Daurade_Gac/Geffroy_Daurade_Gac_DESeq2_Normalized_Counts_miR-125b_5p.csv", sep=""), dec=",", header=TRUE, na.strings="na")
rownames(SBmir125b5p)=SBmir125b5p[,1]

all(rownames(SBdataFishSexchange) %in% rownames(SBmir125b5p))
SBdataFishSexchangok <- SBdataFishSexchange[rownames(SBmir125b5p),]
all(rownames(SBdataFishSexchangok) == rownames(SBmir125b5p))

SBmir125info=merge(SBmir125b5p, SBdataFishSexchangok, by=0, all=TRUE)
SBmir125info$Stages = factor(SBmir125info$sex.name.Sequencing, levels=c('P_M','M_MPP','M_MP','M_I','G_F'))
SSrawdatadesq2MvsFnmds1$"miR-125b-2/3-5p"
names(SBmir125info)[names(SBmir125info) == "miR.125b.2.3.5p"] <- "gene"


my_sample_col1=SSrawdatadesq2MvsFnmds1[,c("Sex", "Species", "Stade", "miR-125b-2/3-5p")]
my_sample_col2=SBmir125info[,c("Origine", "Stages", "gene")]

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

names(my_sample_col1)[names(my_sample_col1) == "miR-125b-2/3-5p"] <- "gene"
my_sample_col1$Stade[which(my_sample_col1$Stade=="Immatures")]="Immature"

########## miR-125b-2/3-5p

my_sample_col2$Names=as.character(my_sample_col2$Stages)

my_sample_col2$Names[which(my_sample_col2$Names=="G_F")]="Female Mature"
my_sample_col2$Names[which(my_sample_col2$Names=="M_I")]="Intersex 2"
my_sample_col2$Names[which(my_sample_col2$Names=="M_MP")]="Intersex 1"
my_sample_col2$Names[which(my_sample_col2$Names=="M_MPP")]="Male Mature"
my_sample_col2$Names[which(my_sample_col2$Names=="P_M")]="Immature Male"
my_sample_col2$Names = factor(my_sample_col2$Names, levels=c("Immature Male","Male Mature", "Intersex 1", "Intersex 2", "Female Mature"))


Dorp <- ggplot(my_sample_col2, aes(x=Names, y=gene)) +
      geom_jitter(aes(color = Names), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"))

Dorp= Dorp+ theme_pubr()
Dorp=Dorp+ labs(title = "Seabream") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                                           axis.title.y = element_text(color = "black", size = 13, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none")+
  stat_summary(fun.data=data_summary, color=c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"), geom = "errorbar", width=.6, size=2)+
  stat_summary(fun.data=data_summary,  geom = "point", color="black", size=3, shape= 5, fill="black") + scale_y_continuous(name="number of reads of miR-125b-2/3-5p")

Dorp


#my_sample_col1=my_sample_col1[which(my_sample_col1$Species!="Red_drum"),] 
Sea_bass=my_sample_col1[which(my_sample_col1$Species=="Sea bass"),] 
Turbot=my_sample_col1[which(my_sample_col1$Species=="Turbot"),] 
Jack=my_sample_col1[which(my_sample_col1$Species=="Blue runner"),] 
Red_drum=my_sample_col1[which(my_sample_col1$Species=="Red drum"),] 


mir125 <- summarySE(my_sample_col1, measurevar="gene", groupvars=c("Species","Sex"))
head(mir125)

### By species
# Seabass

Seabassmir125 <- summarySE(Sea_bass, measurevar="gene", groupvars=c("Stade","Sex"))
head(Seabassmir125)

SBp1 <- ggplot(Sea_bass, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

SBp1= SBp1+ theme_pubr()
SBp1=SBp1+ labs(title = "Sea bass") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                         axis.title.x = element_blank(),  axis.title.y = element_text(color = "black", size = 13, hjust = 0.5), legend.position = "none")+
geom_errorbar(
  aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,
  data = Seabassmir125) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4",'#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Seabassmir125)+
  scale_y_continuous(name="number of reads of miR-125b-2/3-5p")
SBp1


# Turbots

Turbotmir125 <- summarySE(Turbot, measurevar="gene", groupvars=c("Stade","Sex"))
head(Turbotmir125)

Tp1 <- ggplot(Turbot, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

Tp1= Tp1+ theme_pubr()
Tp1=Tp1+ labs(title = "Turbot") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                                   axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")+
  geom_errorbar(
    aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,
    data = Turbotmir125) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4",'#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Turbotmir125)
Tp1


# Jacks

Jackmir125 <- summarySE(Jack, measurevar="gene", groupvars=c("Stade","Sex"))
head(Jackmir125)

Jp1 <- ggplot(Jack, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

Jp1= Jp1+ theme_pubr()
Jp1=Jp1+ labs(title = "Blue runner") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                                  axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")
Jp1=Jp1 +  geom_errorbar(
    aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,data = Jackmir125) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Jackmir125)
Jp1

# Red_drum

Red_drum125 <- summarySE(Red_drum, measurevar="gene", groupvars=c("Stade","Sex"))
head(Red_drum125)

RD1 <- ggplot(Red_drum, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

RD1= RD1+ theme_pubr()
RD1=RD1+ labs(title = "Red Drum") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                                       axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")
RD1=RD1 +  geom_errorbar(
  aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,data = Red_drum125) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Red_drum125)
RD1





AllmiR125b=ggarrange(Dorp,                                           
          ggarrange(SBp1,Tp1,Jp1, RD1, ncol = 4), 
          nrow = 2, heights = c(1, 1)) 
AllmiR125b

ggsave(AllmiR125b, file="Figure 5.pdf", width = 12, height=8, dpi=300)



########## miR-21a-3p

SBdataFishSexchange=read.csv2(paste(dir, "/Daurade_Gac/All Info Seabream.csv", sep=""), dec=",", header=TRUE, na.strings="na")
row.names(SBdataFishSexchange)=SBdataFishSexchange$Name.sequencing
SBmir21a3p=read.csv2(paste(dir, "/Daurade_Gac/Geffroy_Daurade_Gac_DESeq2_Normalized_Counts_miR-21a-3p.csv", sep=""), dec=",", header=TRUE, na.strings="na")
rownames(SBmir21a3p)=SBmir21a3p[,1]

all(rownames(SBdataFishSexchange) %in% rownames(SBmir21a3p))
SBdataFishSexchangok <- SBdataFishSexchange[rownames(SBmir21a3p),]
all(rownames(SBdataFishSexchangok) == rownames(SBmir21a3p))

SBmir23info=merge(SBmir21a3p, SBdataFishSexchangok, by=0, all=TRUE)
SBmir23info$Stages = factor(SBmir23info$sex.name.Sequencing, levels=c('P_M','M_MPP','M_MP','M_I','G_F'))
SSrawdatadesq2MvsFnmds1$"miR-21a-3p"
names(SBmir23info)[names(SBmir23info) == "miR.21a.3p"] <- "gene"


my_sample_col1=SSrawdatadesq2MvsFnmds1[,c("Sex", "Species", "Stade", "miR-21a-3p")]
my_sample_col1$Stade[which(my_sample_col1$Stade=="Immatures")]="Immature"
my_sample_col2=SBmir23info[,c("Origine", "Stages", "gene")]

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

names(my_sample_col1)[names(my_sample_col1) == "miR-21a-3p"] <- "gene"

########## miR-21a-3p

my_sample_col2$Names=as.character(my_sample_col2$Stages)

my_sample_col2$Names[which(my_sample_col2$Names=="G_F")]="Female Mature"
my_sample_col2$Names[which(my_sample_col2$Names=="M_I")]="Intersex 2"
my_sample_col2$Names[which(my_sample_col2$Names=="M_MP")]="Intersex 1"
my_sample_col2$Names[which(my_sample_col2$Names=="M_MPP")]="Male Mature"
my_sample_col2$Names[which(my_sample_col2$Names=="P_M")]="Immature Male"
my_sample_col2$Names = factor(my_sample_col2$Names, levels=c("Immature Male","Male Mature", "Intersex 1", "Intersex 2", "Female Mature"))


Dorp <- ggplot(my_sample_col2, aes(x=Names, y=gene)) +
  geom_jitter(aes(color = Names), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"))

Dorp= Dorp+ theme_pubr()
Dorp=Dorp+ labs(title = "Seabream") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                           axis.title.y = element_text(color = "black", size = 13, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none")+
  stat_summary(fun.data=data_summary, color=c('gold4',"gold3","goldenrod", "gold","lightgoldenrod1"), geom = "errorbar", width=.6, size=2)+
  stat_summary(fun.data=data_summary,  geom = "point", color="black", size=3, shape= 5, fill="black")+ scale_y_continuous(name="number of reads of miR-21a-3p")

Dorp


my_sample_col1=my_sample_col1[which(my_sample_col1$Species!="Red drum"),] 
Sea_bass=my_sample_col1[which(my_sample_col1$Species=="Sea bass"),] 
Turbot=my_sample_col1[which(my_sample_col1$Species=="Turbot"),] 
Jack=my_sample_col1[which(my_sample_col1$Species=="Blue runner"),] 

mir21 <- summarySE(my_sample_col1, measurevar="gene", groupvars=c("Species","Sex"))
head(mir21)

### By species
# Seabass

Seabassmir21 <- summarySE(Sea_bass, measurevar="gene", groupvars=c("Stade","Sex"))
head(Seabassmir21)

SBp1 <- ggplot(Sea_bass, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

SBp1= SBp1+ theme_pubr()
SBp1=SBp1+ labs(title = "Sea bass") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                           axis.title.y = element_text(color = "black", size = 13, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none")+
  geom_errorbar(
    aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,
    data = Seabassmir21) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4",'#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Seabassmir21)+
  scale_y_continuous(name="number of reads of miR-21a-3p")
SBp1


# Turbots

Turbotmir21 <- summarySE(Turbot, measurevar="gene", groupvars=c("Stade","Sex"))
head(Turbotmir21)

Tp1 <- ggplot(Turbot, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

Tp1= Tp1+ theme_pubr()
Tp1=Tp1+ labs(title = "Turbot") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                                    axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")+
  geom_errorbar(
    aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,
    data = Turbotmir21) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4",'#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Turbotmir21)
Tp1


# Jacks

Jackmir21 <- summarySE(Jack, measurevar="gene", groupvars=c("Stade","Sex"))
head(Jackmir21)

Jp1 <- ggplot(Jack, aes(x=Stade, y=gene, color="Sex")) +
  geom_jitter(aes(color = Sex), size = 4, position=position_jitter(0.2))+
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))

Jp1= Jp1+ theme_pubr()
Jp1=Jp1+ labs(title = "Blue runner") +theme(plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "italic"), 
                                                  axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")
Jp1=Jp1 +  geom_errorbar(
  aes(ymin = gene-se, ymax = gene+se, color=Sex), width=0.2, size=2,data = Jackmir21) +
  scale_colour_manual(values =  c('#EFC000FF',"darkgoldenrod4"))+ geom_point(color="black", size=3, shape= 5, fill="black", data = Jackmir21)
Jp1



AllmiR21a=ggarrange(Dorp,                                                 # First row with scatter plot
                    ggarrange(SBp1,Tp1,Jp1, ncol = 3), # Second row with box and dot plots
                    nrow = 2, heights = c(1.7, 2)) 
AllmiR21a

ggsave(AllmiR21a, file="Figure 4.pdf", width = 9, height=8, dpi=300)

