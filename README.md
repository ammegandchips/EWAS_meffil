# This isn't finished yet! I'm hoping to finish it over the next couple of days. Come back then :)
My EWAS pipeline using meffil
Written by Gemma Sharp on the 19th of May 2017

This is step-by-step description of how I would run an EWAS using ARIES (ALSPAC) data. I'm not saying this is the BEST way to run an EWAS, but it's the best way I've found (and I've spent a fair bit of time thinking about what the best way might be).
If you find a better way, then please let me know!

# Set up
In my home directory on Bluecrystal 3, I have a folder called `Common_files`. In this folder, I store two files that I think I'll need for multiple EWAS:
`naeem_list.csv` : a list of probes that have been identified as potentially problematic (on a SNP, cross-hybridising, etc) according to Naeem et al.
`aries-detailed-cell-counts-20150409.rda` : a list of dataframes (one for each ARIES time point) with cell counts estimated using the Houseman methods but where Granulocytes are split into Eosinophils and Neutrophils.
`meffil_EWAS_script.r` : my EWAS R script (more on this below). This should be generic (i.e. you can use the same script to run any ARIES EWAS, providing you don't want to do anything too fancy).

I've uploaded these three files to this github repository, but please note that `aries-detailed-cell-counts-20150409.rda`is password protected and only relevant when you want to adjust for 'detailed' cell counts (i.e. including eosinophils and neutrophils). You can email me for the password.

I then have another folder called `EWAS` with subfolders for different projects. For example, I have projects called `pace_mat_bmi`, `pace_mat_alcohol`, `alspac_mat_depression` `cleft_subtypes`. Within each project folder I have further subfolders:

`submission_scripts` : all the .sh files I use to 'call' the R EWAS script. These .sh files are unique for each EWAS because they list the arguments that are read into the R script at the start of the EWAS.
`ewas_results` : This is where I save all the outputs from the EWAS(s).

# . sh files
For each EWAS I run, I make one of these .sh files and then submit it as a job to Bluecrystal using:
`qsub EWAS/project_name/submission_scripts/model_1.sh`
So if I wanted to run 10 EWAS (e.g. 10 different models within the same project), I would make 10 .sh files and then submit them as separate jobs. There must be an easier way, but I don't know what it is and to be honest this isn't too arduous.
I've uploaded an example of a .sh file to this repository (`model_1.sh`) but here's a better annotated version:

`
#!/bin/bash
#
#
#PBS -l nodes=1:ppn=1,walltime=12:00:00


WORK_DIR="/panfs/panasas01/sscm/gs8094/EWAS/Maternal_depression"
module add languages/R-3.0.2

cd $WORK_DIR

R CMD BATCH --no-save --no-restore '--args Trait CellData CellAdj Phenofile Method Removal BorM TP PACE Covariates Crude_or_Adj WD' /panfs/panasas01/sscm/gs8094/EWAS/GestationalAge/ALSPAC/GA_ALSPAC.r /panfs/panasas01/sscm/gs8094/EWAS/Maternal_depression/Model1.out
`

# The EWAS R script
The EWAS runs in R

# pull in args
args	<- 	commandArgs(trailingOnly = TRUE)
Trait 	<- 	toString(args[1]) #Phenotypic exposure or outcome of interest
CellData <-toString(args[2]) #Which cell counts should we use? (houseman, gse68456, gervinandlyle, andrews-and-bakulski, houseman_eos)
CellAdj   <-      toString(args[3]) #Cell adjusted? (noCells, Cells)
Phenofile       <-      toString(args[4]) #Path to file containing all phenotype information (must be a dta stata version 12)
BorM    <-      toString(args[5]) #Betas or M-values (B or M)
TP <-  toString(args[6]) #Time point (cord or F7 or 15up or antenatal or FOM)
Covariates <-toString(args[7]) #list of covariates (eg: m_age,mum_uni,matsm,parity i.e. commas but no spaces or quotation marks)
WD<- toString(args[8]) #Working directory (eg /panfs/panasas01/sscm/gs8094/EWAS/)

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
ewas.summary <- meffil.ewas.summary(ewas.results, meth, parameters=ewas.parameters)   
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
