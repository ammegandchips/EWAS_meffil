###################################################################################################################
###################################################################################################################
#####                                Full code for any  EWAS analysis                                         #####
#####                  Compiled 12/06/2015 by Gemma Sharp  - UPDATED FOR VERSION 3 20/05/2016                 #####
#####                               UPDATED FOR MEFFIL 19/05/2017                                             #####
###################################################################################################################
###################################################################################################################

# pull in args
args	<- 	commandArgs(trailingOnly = TRUE)
Trait 	<- 	toString(args[1]) #Phenotypic exposure or outcome of interest
CellData <-toString(args[2]) #Which cell counts should we use? (houseman, gse68456, gervinandlyle, andrews-and-bakulski, houseman_eos)
CellAdj   <-      toString(args[3]) #Cell adjusted? (noCells, Cells)
Phenofile       <-      toString(args[4]) #Path to file containing all phenotype information (must be a dta stata version 12)
BorM    <-      toString(args[5]) #Betas or M-values (B or M)
TP <-  toString(args[6]) #Time point (cord or F7 or 15up or antenatal or FOM)
Covariates <-toString(args[7]) #list of covariates (eg: m_age,mum_uni,matsm,parity i.e. commas but no spaces or quotation marks)
WD<- toString(args[8]) #Working directory (eg /panfs/panasas01/sscm/gs8094/EWAS/example_project)

print(Trait)
print(CellData)
print(CellAdj)
print(Phenofile)
print(BorM)
print(TP)
print(Covariates)
print(WD)

#set working directory
setwd(WD)

#load packages
library(foreign) #to read stata file
library(meffil)

#Load description of samples
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/samplesheet/data.Robj")
samplesheet<-subset(samplesheet, time_point==TP)
if(TP !="antenatal" & TP !="FOM"){
qletB<-samplesheet$ALN[which(samplesheet$QLET=="B")] #find alns for multiple pregnancies
samplesheet<-samplesheet[-which(samplesheet$ALN %in% qletB),] #remove multiple pregnancies
}

#load the methylation data
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj")
meth <- norm.beta.random[,samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(norm.beta.random)

#load detection P-values (used to filter all probes with a high detection P-value)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/detection_p_values/data.Robj")
pvals <- detp[,samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(detp)

#load annotation data
annotation <- meffil.get.features("450k")

#Filter meth data (remove sex chromosomes and SNPs and probes with high detection P-values)
pvalue_over_0.05 <- pvals > 0.05
count_over_0.05 <- rowSums(sign(pvalue_over_0.05))
Probes_to_exclude_Pvalue <- rownames(pvals)[which(count_over_0.05 > ncol(pvals)*0.05)]
XY <- as.character(annotation$name[which(annotation$chromosome %in% c("chrX", "chrY"))])
SNPs.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])
annotation<- annotation[-which(annotation$name %in% c(XY,SNPs.and.controls,Probes_to_exclude_Pvalue)),]
meth <- subset(meth, row.names(meth) %in% annotation$name)
paste("There are now ",nrow(meth), " probes")
paste(length(XY),"were removed because they were XY")
paste(length(SNPs.and.controls),"were removed because they were SNPs/controls")
paste(length(Probes_to_exclude_Pvalue),"were removed because they had a high detection P-value")
rm(XY, SNPs.and.controls, pvals, count_over_0.05, pvalue_over_0.05, Probes_to_exclude_Pvalue)

#Load phenotype data (this should be stored in your working directory)
Pheno<-read.dta(paste0(Phenofile,".dta"))

#Load cell-counts
if(TP=="cord"){
cells<-read.table(paste0("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/derived/cellcounts/cord/",CellData,"/data.txt"),header=T)
}else{
if(CellData=="houseman_eos"){
load("/panfs/panasas01/sscm/gs8094/Common_files/aries-detailed-cell-counts-20150409.rda")
cells<-detailed.cell.counts[[TP]]
}else{
cells<-read.table("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/derived/cellcounts/houseman/data.txt", header=TRUE)
}}

#Add Sample_Name to Pheno (assuming Pheno contains aln)
Pheno<-merge(Pheno,samplesheet[,c("ALN","Sample_Name")],by.x="aln",by.y="ALN")

#Prepare phenotype data
Covs<-strsplit(Covariates,split=",")[[1]]
Pheno<-na.omit(Pheno[,c("Sample_Name",Trait,Covs)])

#Merge Pheno with cell-counts
colnames(cells)[1]<-"Sample_Name"
Pheno<-merge(Pheno,cells,by.x="Sample_Name",by.y="Sample_Name")

#Match meth to Pheno
meth<-meth[,na.omit(match(Pheno$Sample_Name,colnames(meth)))]
Pheno<-Pheno[match(colnames(meth),Pheno$Sample_Name),]

ifelse(all(Pheno$Sample_Name==colnames(meth)), "meth and phenotype data successfully matched :) ","Data not matched :(")

Pheno<-droplevels(Pheno) # get rid of any empty factor levels 

#Little summary of the phenotype data at this point
paste("There are ",nrow(Pheno)," people in this analysis")
"Here's a summary of the phenotype data:"
summary(Pheno)

# Convert to M-values if necessary
if(BorM=="M"){
	meth <- log2(meth/(1-meth))
}

#Include cell counts in the EWAS model?
if(CellAdj=="noCells"){
Pheno<-Pheno[,-which(colnames(Pheno) %in% colnames(cells)[-1])]
}

#Run EWAS using meffil

obj <- meffil.ewas(meth, variable=Pheno[,2], covariates=Pheno[,-(1:2)], winsorize.pct = NA  ,most.variable = min(nrow(meth), 20000), outlier.iqr.factor=3, verbose=TRUE)
ewas_res<-data.frame(
	probeID=rownames(obj$analyses$none$table),
	coef.none=obj$analyses$none$table$coefficient,
	se.none=obj$analyses$none$table$coefficient.se,
	p.none=obj$analyses$none$table$p.value,
	coef.all=obj$analyses$all$table$coefficient,
	se.all=obj$analyses$all$table$coefficient.se,
	p.all=obj$analyses$all$table$p.value,
	coef.sva=obj$analyses$sva$table$coefficient,
	se.sva=obj$analyses$sva$table$coefficient.se,
	p.sva=obj$analyses$sva$table$p.value,
	coef.isva=obj$analyses$isva$table$coefficient,
	se.isva=obj$analyses$isva$table$coefficient.se,
	p.isva=obj$analyses$isva$table$p.value
	)

ewas.parameters<-meffil.ewas.parameters(sig.threshold=1e-5, max.plots=5, model="isva")
ewas.summary <- meffil.ewas.summary(obj, meth, parameters=ewas.parameters)   
savefile <- paste("ewas_results/",Trait,TP,Covariates,CellAdj,Sys.Date(),".html", sep = "_")
meffil.ewas.report(ewas.summary, output.file=savefile,author="gemma sharp", study="alspac")

#Create N_for_probe column
ewas_res$original_n<-rowSums(!is.na(meth))
outliers_n<-data.frame(table(c(rownames(obj$too.hi),rownames(obj$too.lo))))
colnames(outliers_n)<-c("probeID","n_outliers")
ewas_res<-merge(ewas_res,outliers_n,by="probeID",all.x=TRUE)
ewas_res$n_outliers<-ewas_res$n_outliers*-1
ewas_res$final_n<-rowSums(ewas_res[,c("original_n","n_outliers")],na.rm=TRUE)

#Create N_cases column if necessary
if(class(obj$variable)=="factor"|any(as.numeric(Pheno[,2])!=0&as.numeric(Pheno[,2])!=1)==FALSE|class(Pheno[,2])=="character"){
print("Phenotype of interest is binary")
outliers_cases<-as.data.frame(rbind(obj$too.hi,obj$too.lo))
outliers_cases<-table(outliers_cases$row, outliers_cases$col)
outliers_cases<-data.frame(probeID=rownames(meth)[as.numeric(rownames(outliers_cases))],
	n_outliers_cases=rowSums(outliers_cases[,as.factor(obj$variable)==levels(obj$variable)[2]]),
	n_outliers_controls=rowSums(outliers_cases[,as.factor(obj$variable)==levels(obj$variable)[1]]))
ewas_res<-merge(ewas_res,outliers_cases,by="probeID",all=TRUE)
}

# Now we can include more information in our EWAS results file:

#Load the Naeem list of problematic probes
Naeem<-read.csv("/panfs/panasas01/sscm/gs8094/Common_files/naeem_list.csv")
ewas_res$OnNaeem<-ifelse(ewas_res$probeID %in% Naeem$EXCLUDE_PROBES,"yes","no")

# Adjustment for multiple testing
ewas_res$fdr.none<-p.adjust(ewas_res$p.none, method="fdr") 
ewas_res$bonferroni.none<-p.adjust(ewas_res$p.none, method="bonferroni")
ewas_res$fdr.all<-p.adjust(ewas_res$p.all, method="fdr") 
ewas_res$bonferroni.all<-p.adjust(ewas_res$p.all, method="bonferroni")
ewas_res$fdr.sva<-p.adjust(ewas_res$p.sva, method="fdr") 
ewas_res$bonferroni.sva<-p.adjust(ewas_res$p.sva, method="bonferroni")
ewas_res$fdr.isva<-p.adjust(ewas_res$p.isva, method="fdr") 
ewas_res$bonferroni.isva<-p.adjust(ewas_res$p.isva, method="bonferroni")

##Annotate and sort
#Add information about EWAS
ewas_res$Trait<-Trait
ewas_res$Covariates<-paste0(ls(obj$covariates),collapse=", ")
ewas_res$nSVs<-ncol(obj$analyses$sva$design)
ewas_res$nISVs<-ncol(obj$analyses$isva$design)
ewas_res$TP<-TP
ewas_res$BorM<-BorM
ewas_res$CellData<-CellData
ewas_res$CellAdj<-CellAdj

#Add Lambda (a measure of test statistic inflation)
Lambda<-function(P){
chisq <- qchisq(1-P,1)
median(chisq,na.rm=T)/qchisq(0.5,1)
}

ewas_res$lambda.none<-Lambda(ewas_res$p.none)
ewas_res$lambda.all<-Lambda(ewas_res$p.all)
ewas_res$lambda.sva<-Lambda(ewas_res$p.sva)
ewas_res$lambda.isva<-Lambda(ewas_res$p.isva)

#Add annotation information about probes and sort by P-value
ewas_res<-merge(ewas_res,annotation,by.x="probeID", by.y="name",all.x=TRUE)
ewas_res<-ewas_res[order(ewas_res$p.isva),]

#Save as an Rdata file
savefile <- paste("ewas_results/",Trait,TP,Covariates,CellAdj,Sys.Date(),".Rdata", sep = "_")
save(ewas_res,file=savefile)
