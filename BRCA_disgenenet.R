library(Kendall)
library(dplyr)

timestart <- Sys.time()
cat("Start time is: ",timestart,"\n")

pval_diff <- read.csv("BRCAdata_filter.csv",header=T,check.names=FALSE) ## Differentially expressed dataset
pval_diff_miRNA <- read.csv("miRNAdata_filter.csv",header=T,check.names=FALSE) ## miRNA Differentially expressed dataset
key_genes <- read.csv("key_gene_filter.csv",header=T,check.names=FALSE)
TF_miRNA_association <- read.csv("TF_miRNA_association.csv",header=T,check.names=FALSE)
miRNA_gene_es <- read.csv("mirna_gene_noTF.csv",header=T,check.names=FALSE)


pval_diff[1:5,1:5]
pval_diff_miRNA[1:5,1:5]

dim(pval_diff)
dim(miRNA_gene_es)
dim(TF_miRNA_association)
n1 <- 80        ## control group
n2 <- 461        ## case group


TFsexpress <- pval_diff[,which(colnames(pval_diff)%in%colnames(pval_diff[,2:176]))]  ## TFs expression dataset

dim(TFsexpress) 
TFsexpress[1:5,1:5]

full_data <- pval_diff[,-1]  ## full TF-GENE dataset
full_data_miRNA <- pval_diff_miRNA[,-1]   ## full miRNA dataset
full_data_miRNA_gene_es <- miRNA_gene_es[,-1]

dim(full_data_miRNA_gene_es)

exp_new <- full_data[1:461,]   ## case group tf/gene
ctl_new <- full_data[462:541,]  ## control group tf/gene
exp_new_miRNA <- full_data_miRNA[1:461,]   ## case group miRNA
ctl_new_miRNA <- full_data_miRNA[462:541,]   ## control group miRNA
exp_new_miRNA_gene_es <- full_data_miRNA_gene_es[1:461,]
ctl_new_miRNA_gene_es <- full_data_miRNA_gene_es[462:541,]


comm_ctl_exp=intersect(colnames(exp_new),colnames(ctl_new))
comm_ctl_exp_miRNA=intersect(colnames(exp_new_miRNA),colnames(ctl_new_miRNA))
ctl_new1=ctl_new[,which(colnames(ctl_new)%in%comm_ctl_exp)]
exp_new1=exp_new[,which(colnames(exp_new)%in%comm_ctl_exp)]
ctl_new1_miRNA=ctl_new_miRNA[,which(colnames(ctl_new_miRNA)%in%comm_ctl_exp_miRNA)]
exp_new1_miRNA=exp_new_miRNA[,which(colnames(exp_new_miRNA)%in%comm_ctl_exp_miRNA)]
comm_ctl_exp_miRNA_gene_es=intersect(colnames(exp_new_miRNA_gene_es),colnames(ctl_new_miRNA_gene_es))
ctl_new1_miRNA_gene_es=ctl_new_miRNA_gene_es[,which(colnames(ctl_new_miRNA_gene_es)%in%comm_ctl_exp_miRNA_gene_es)]
exp_new1_miRNA_gene_es=exp_new_miRNA_gene_es[,which(colnames(exp_new_miRNA_gene_es)%in%comm_ctl_exp_miRNA_gene_es)]

dim(exp_new1)
dim(ctl_new1)
dim(exp_new1_miRNA)
dim(ctl_new1_miRNA)
dim(exp_new1_miRNA_gene_es)
dim(ctl_new1_miRNA_gene_es)

TFs <- colnames(ctl_new1)[colnames(ctl_new1)%in%colnames(TFsexpress)]  ## TFs

use_TFs <- TFs
only_genes <- colnames(ctl_new1)[!colnames(ctl_new1)%in%TFs]  ## genes not TFs
only_miRNAs <- colnames(ctl_new1_miRNA)[!colnames(ctl_new1_miRNA)%in%TFs]  ## miRNA



## calculating the differential connectivity
s <- NULL

top_genes <- c(1,key_genes[,2])
##key <- matrix(top_genes,1,length(top_genes))
key <- diag(top_genes)

##partnumber = c(8,7,6,5,4,3,2,1)
for(alp in 1:8){
for(bet in 1:8){
alpha=alp/10
beta=bet/10
if(alpha+beta<1){
gamma=1-alpha-beta


for(l in 1:length(use_TFs)){
sub <- c(use_TFs[l],only_genes,only_miRNAs)



pop_exp_new <- exp_new1[,colnames(exp_new1)%in%sub]
pop_ctl_new <- ctl_new1[,colnames(ctl_new1)%in%sub]
pop_exp_new_miRNA <- exp_new1_miRNA[,colnames(exp_new1_miRNA)%in%sub]
pop_ctl_new_miRNA <- ctl_new1_miRNA[,colnames(ctl_new1_miRNA)%in%sub]

## next change
new_TF <- as.character(TF_miRNA_association$TF)
miRNA_rebyTF <- filter(TF_miRNA_association,new_TF == use_TFs[l])
use_miRNAs <-miRNA_rebyTF[,2]

sum_s_miRNA_Rreby_TF <- 0

if(length(use_miRNAs)>0){
s_miRNA_Rreby_TF <- NULL

for(m in 1:length(use_miRNAs)){
sub2 <- c(use_miRNAs[m],only_genes)
pop_exp_new_miRNA_gene_es <- exp_new1_miRNA_gene_es[,colnames(exp_new1_miRNA_gene_es)%in%sub2]
pop_ctl_new_miRNA_gene_es <- ctl_new1_miRNA_gene_es[,colnames(ctl_new1_miRNA_gene_es)%in%sub2]
s1_miRNA_rebyTF=cor(pop_exp_new_miRNA_gene_es[,use_miRNAs[m]],pop_exp_new_miRNA_gene_es,method="pearson")
s2_miRNA_Rreby_TF=cor(pop_ctl_new_miRNA_gene_es[,use_miRNAs[m]],pop_ctl_new_miRNA_gene_es,method="pearson")
s_miRNA_Rreby_TF[m] <- sum(abs(s1_miRNA_rebyTF-s2_miRNA_Rreby_TF))/length(only_genes)
}
sum_s_miRNA_Rreby_TF <- sum(s_miRNA_Rreby_TF)/length(use_miRNAs)
}




##

s1=cor(pop_exp_new[,use_TFs[l]],pop_exp_new,method="pearson")
s2=cor(pop_ctl_new[,use_TFs[l]],pop_ctl_new,method="pearson")
s3=cor(pop_exp_new_miRNA[,use_TFs[l]],pop_exp_new_miRNA,method="pearson")
s4=cor(pop_ctl_new_miRNA[,use_TFs[l]],pop_ctl_new_miRNA,method="pearson")

##s[l] <- alpha*(sum(abs(s1-s2)%*%key)/length(only_genes))+beta*(sum(abs(s3-s4))/length(only_miRNAs))+gamma*sum_s_miRNA_Rreby_TF
s[l] <- alpha*(sum(abs(s1-s2)%*%key)/length(only_genes))+beta*(sum(abs(s3-s4))/length(only_miRNAs))+gamma*sum_s_miRNA_Rreby_TF

}

DEgenes_pop <- data.frame(use_TFs,s)
TfTestStat_pop=DEgenes_pop[with(DEgenes_pop, order(-DEgenes_pop$s)),]
                           ## Differential connectivity scores for TFs
filename<-paste("DEindividualtestpval1_",alpha,"_",beta,".csv")
write.csv(TfTestStat_pop,filename)


#########################################################################
## Step2b)


PulledTFexpressTrans=TFsexpress[,colnames(TFsexpress)%in%TFs]  ## Expression set for TFs only
PulledTFexpressTrans[1:5,1:5]
dim(PulledTFexpressTrans)

temp <- PulledTFexpressTrans[, !duplicated(colnames(PulledTFexpressTrans))]

cormatTF=cor(temp,method="pearson")  ## Correlation among TFs

#####################################################################
## Step 3 & 4

#TfTestStat_pop=read.csv("DEindividualtestpval1.csv",header=T,check.names=FALSE)[,-1]

## Calculating the Kendall's statistic
stat_tau=NULL
for(t in 1:length(TfTestStat_pop[,1])){
TranInterested=as.character(TfTestStat_pop[t,1])

orderCorrTF=sort(abs(cormatTF[TranInterested,]),decreasing=TRUE)

d1=data.frame(gene=TfTestStat_pop[,1],TfTestStat_pop)
d2=data.frame(orderCorrTF)
d22=data.frame(d2,gene=rownames(d2))
data=merge(d22, d1, by="gene")

test_tau=Kendall(data[,2], data[,4])
stat_tau[t]=test_tau$tau
}
allTFsStat_pop=data.frame(TransFacs=TfTestStat_pop[,1],stat_tau)  ## Kendall's statistic value for all TFs
max_tau=max(allTFsStat_pop[,2])
allTFsStat_pop1 <- allTFsStat_pop[with(allTFsStat_pop,order(-allTFsStat_pop[,2])),]
                                      ## Ordered list of TFs based on their test statistic values.



filename<-paste("allTFsStat_pop1_",alpha,"_",beta,".csv")
write.csv(allTFsStat_pop1,filename)

print(" Now ! It is time for Bootstrap!!")







############################  P-Value calculation #######################################################

###########################################################################################################
##############################################################################################################
## Sampling (p-value calculation)

###############################################################################
full_data <- rbind(exp_new1,ctl_new1)
full_data_miRNA <- rbind(exp_new1_miRNA,ctl_new1_miRNA)   ## full miRNA dataset
full_data_miRNA_gene_es <-rbind(exp_new1_miRNA_gene_es,ctl_new1_miRNA_gene_es)
all_genes <- unique(colnames(exp_new1))
all_miRNAs <- unique(colnames(exp_new1_miRNA))
only_genes <- all_genes[!all_genes%in%TFs]
only_miRNAs <- all_miRNAs[!all_miRNAs%in%TFs]  ## miRNA



p_tau <- NULL
sum <- 0
recount <- 100
p_value <- 0

sam_top_genes <- c(1,key_genes[,2])
sam_key <- diag(sam_top_genes)

for(k in 1:recount){ ## change   (No. of iterations)9
n <- 541


x <- sample(1:541,size=n,replace=TRUE)
sam_exp=full_data[x[1:461],]
sam_ctl=full_data[x[462:541],]                ## Bootstrap samples
sam_exp_miRNA=full_data_miRNA[x[1:461],]
sam_ctl_miRNA=full_data_miRNA[x[462:541],]
sam_exp_miRNA_gene_es=full_data_miRNA_gene_es[x[1:461],]
sam_ctl_miRNA_gene_es=full_data_miRNA_gene_es[x[462:541],]


## Bootstrap based differential connectivity
s <- NULL
for(l in 1:length(use_TFs)){  
sub <- c(use_TFs[l],only_genes,only_miRNAs)



sam_exp_new <- sam_exp[,colnames(sam_exp)%in%sub]
sam_ctl_new <- sam_ctl[,colnames(sam_ctl)%in%sub]
sam_exp_new_miRNA <- sam_exp_miRNA[,colnames(sam_exp_miRNA)%in%sub]
sam_ctl_new_miRNA <- sam_ctl_miRNA[,colnames(sam_ctl_miRNA)%in%sub]


## next change
new_TF <- as.character(TF_miRNA_association$TF)
miRNA_rebyTF <- filter(TF_miRNA_association,new_TF == use_TFs[l])
use_miRNAs <-miRNA_rebyTF[,2]

sum_s_miRNA_Rreby_TF <- 0

if(length(use_miRNAs)>0){
s_miRNA_Rreby_TF <- NULL

for(m in 1:length(use_miRNAs)){
sub2 <- c(use_miRNAs[m],only_genes)
sam_exp_new_miRNA_gene_es <- sam_exp_miRNA_gene_es[,colnames(sam_exp_miRNA_gene_es)%in%sub2]
sam_ctl_new_miRNA_gene_es <- sam_ctl_miRNA_gene_es[,colnames(sam_ctl_miRNA_gene_es)%in%sub2]
s1_miRNA_rebyTF=cor(sam_exp_new_miRNA_gene_es[,use_miRNAs[m]],sam_exp_new_miRNA_gene_es,method="pearson")
s2_miRNA_Rreby_TF=cor(sam_ctl_new_miRNA_gene_es[,use_miRNAs[m]],sam_ctl_new_miRNA_gene_es,method="pearson")
s_miRNA_Rreby_TF[m] <- sum(abs(s1_miRNA_rebyTF-s2_miRNA_Rreby_TF))/length(only_genes)
}
sum_s_miRNA_Rreby_TF <- sum(s_miRNA_Rreby_TF)/length(use_miRNAs)
}


s1=cor(sam_exp_new[,use_TFs[l]],sam_exp_new,method="pearson")
s2=cor(sam_ctl_new[,use_TFs[l]],sam_ctl_new,method="pearson")
s3=cor(sam_exp_new_miRNA[,use_TFs[l]],sam_exp_new_miRNA,method="pearson")
s4=cor(sam_ctl_new_miRNA[,use_TFs[l]],sam_ctl_new_miRNA,method="pearson")   

##s[l] <- alpha*(sum(abs(s1-s2)%*%sam_key)/length(only_genes))+beta*(sum(abs(s3-s4))/length(only_miRNAs))+gamma*sum_s_miRNA_Rreby_TF
s[l] <- alpha*(sum(abs(s1-s2)%*%sam_key)/length(only_genes))+beta*(sum(abs(s3-s4))/length(only_miRNAs))+gamma*sum_s_miRNA_Rreby_TF

}
DEgenes_sam <- data.frame(use_TFs,s) 
TfTestStat=DEgenes_sam[with(DEgenes_sam, order(-DEgenes_sam$s)),]


data_req=rbind(sam_ctl,sam_exp)
TFexpress=data_req[,which(colnames(data_req)%in%TFs)]

cormatTF=cor(TFexpress,method="pearson")

## Bootstrap based Kendall's test statistic values
stat=NULL
for(t in 1:length(TFs)){ 
TranInterested=as.character(TfTestStat[t,1])

orderCorrTF=sort(abs(cormatTF[TranInterested,]),decreasing=TRUE)

d1=data.frame(gene=TfTestStat[,1],TfTestStat)
d2=data.frame(orderCorrTF)
d22=data.frame(d2,gene=rownames(d2))
data=merge(d22, d1, by="gene")

test_tau=Kendall(data[,2], data[,4])
stat[t]=test_tau$tau
}
allTFsStat=data.frame(TransFacs=TfTestStat[,1],stat) 
Stat_obs_tau=max(allTFsStat[,2])
p_tau[k] <- Stat_obs_tau

status <- p_tau  ## Bootstrap based test statistic values
if(p_tau[k]>max_tau){
sum=sum+1;
}

cat("This is the times of: ",k,"\n")

}

filename<-paste("statusvec1_",alpha,"_",beta,".csv")
write.csv(status,filename)

p_value=sum/recount

filename<-paste("p_value_",alpha,"_",beta,".csv")
write.csv(p_value,filename)


cat("The p_value is : ",p_value,"\n")

timeend <- Sys.time()
cat("End time is: ",timeend,"\n")

cat("The runtime is : ",timeend-timestart,"\n")



}
}
}
