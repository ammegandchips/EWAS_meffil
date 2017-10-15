# My EWAS pipeline using meffil

This is step-by-step description of how I would run an EWAS using ARIES (ALSPAC) data and the R package meffil. I'm not saying this is the BEST way to run an EWAS, but it's the best way I've found (and I've spent a fair bit of time thinking about what the best way might be).
If you find a better way, then please let me know!

# What type of EWAS will this run?

Features of the EWAS performed through this pipeline:
* linear regression 
* methylation is the outcome, your trait of interest is the exposure
* four models are actually run for each input file, these are:

   1) **none:** no covariates
   2) **all:** all covariates
   3) **sva:** all covariates + surrogate variables generated using sva (surrogate variable analysis)
   4) **isva:** all covariates + surrogate variables generated using isva (independent surrogate variable analysis)
   
* you'll get the results for all four in your EWAS results

If you want to run a more advanced EWAS (e.g., methylation as the exposure, logistic regression, etc), this script isn't the one! I'm hoping to write a similar tutorial for more advanced EWAS when I get a chance. For now, you can email me and I'll share what code I have.

# Set up
In my home directory on Bluecrystal 3, I have a folder called `Common_files`. In this folder, I store three files that I think I'll need for multiple EWAS:

* `naeem_list.csv` : a list of probes that have been identified as potentially problematic (on a SNP, cross-hybridising, etc) according to Naeem et al.

* `aries-detailed-cell-counts-20150409.rda` : a list of dataframes (one for each ARIES time point) with cell counts estimated using the Houseman methods but where Granulocytes are split into Eosinophils and Neutrophils.

* `meffil_EWAS_script.r` : my EWAS R script (more on this below). This should be generic (i.e. you can use the same script to run any ARIES EWAS, providing you don't want to do anything too fancy).


I've uploaded these three files to this github repository, but please note that `aries-detailed-cell-counts-20150409.rda`is password protected and only relevant when you want to adjust for 'detailed' cell counts (i.e. including eosinophils and neutrophils). You can email me for the password.

I then have another folder called `EWAS` with subfolders for different projects. For example, I have projects called `pace_mat_bmi`, `pace_mat_alcohol`, `alspac_mat_depression` `cleft_subtypes`. Within each project folder I have further subfolders:

* `submission_scripts` : all the .sh files I use to 'call' the R EWAS script. These .sh files are unique for each EWAS because they list the arguments that are read into the R script at the start of the EWAS.

* `ewas_results` : This is where I save all the outputs from the EWAS(s).

I also store my phenotype file (i.e. my file containing the variable of interest and any covariates) in this project folder.

# . sh files

For each EWAS I run, I make one of these .sh files.

So if I wanted to run 10 EWAS (e.g. 10 different models within the same project), I would make 10 .sh files and then submit them as separate jobs. There must be an easier way, but I don't know what it is and to be honest this isn't too arduous.
I've uploaded an example of a .sh file to this repository (`example.sh`) but here's a walk-through:

At the top of the file, I have this little bit of code. I'm not totally sure what the first three lines are doing, but they're important (at least the first one is!). The fourth line is where I set the number of nodes and the time I think it's going to take to run. 12 hours is an over-estimation, but this means the job gets sent to the "very short" queue and is unlikely to be killed before it has finished running, so it probably makes sense to leave all this as it is:

```
#!/bin/bash
#
#
#PBS -l nodes=1:ppn=1,walltime=12:00:00
```

Next, I set the working directory, which should be the project folder I mentioned above:

```
WORK_DIR="/panfs/panasas01/sscm/gs8094/EWAS/example_project"
```

Next, I tell it that I want to install R, specifying which version. You can get a list of versions on bluecrystal by typing `module avail` (R is listed under languages). In this example, I'm using R-3.0.2, but you can use a more up-to-date version if you want. (At the time of writing, the most recent version available on bluecrystal is R-3.4.1-ATLAS.)

```
module add languages/R-3.0.2
```

The next line just changes the directory to the working directory I specified in line 5:

```
cd $WORK_DIR
```

Everything I've entered so far is likely to be the same for each EWAS I run within the same project. The next line is where I have to change things for each EWAS I want to run.

The first bit of this line (`R CMD BATCH --no-save --no-restore`) is just saying that I want to run this in R. This bit stays the same for all EWAS. 

The next bit is where we sent the 'arguments', which are basically just inputs to the EWAS R script:

```
R CMD BATCH --no-save --no-restore '--args Trait CellData CellAdj Phenofile Method Removal BorM TP PACE Covariates Crude_or_Adj WD' /panfs/panasas01/sscm/gs8094/Common_files/meffil_EWAS_script.r /panfs/panasas01/sscm/gs8094/EWAS/example_project/example.out
```

`--args ` stays the same, but `Trait CellData CellAdj Phenofile Method Removal BorM TP PACE Covariates Crude_or_Adj WD` are the bits I change for each EWAS I run.

* **Trait:** The trait of interest, exactly as it appears in the Phenofile (i.e. the exact variable name)
* **CellData:** Which cell counts should we use? Options are houseman or houseman_eos for whole blood or one of the cord blood references (gse68456, gervinandlyle, andrews-and-bakulski)
* **CellAdj:** Do we want to adjust for estimated cells? Options are Cells or noCells
* **Phenofile:** Path to the file containing all the phenotype information, which must be a .dta file from stata version 12 (saveold in STATA). If this is a hassle (which it often is), you could just change the line where this file is read in the R script (e.g. instead of read.dta use read.csv and save Phenofile as a .csv file). 
* **BorM:** Do we want to run our EWAS on beta values (scale of 0 to 1) or M-values (logit transformation of beta values)? You can read about the pros and cons of each in a paper by Du et al. However, I always go for beta values because the resulting coefficient is more intuitive. Options for this argument are B or M.
* **TP:** What ARIES time-point are we interested in? Options are cord, F7, 15up, antenatal or FOM
* **Covariates:** List of covariates we want to adjust for in this model. All these variables need to be in Phenofile. Separate the list here by commas, but don't include spaces or quotation marks (for example: mat_age,mat_ses,mat_smoke,parity).
* **WD:** The directory for this EWAS project (i.e. the project folder)  

The final bits of this line are just the file path to the R script and the file path to an `.out` file, which will show the output from R and is basically just there so you can check how your EWAS is doing and troubleshoot any problems if it doesn't finish (i.e. you can see how far it got through the script and any error messages R might spit out).

The path to the EWAS script should stay the same (i.e. it should be in your Common_files folder). The filename of the out file will change depending on what EWAS you're running in your model. Here, it's called `continuous_cell_adj.out`. I like to call my .out files the same thing as the corresponding .sh file, because I find it easier to track the input to the output. Note that these `.out` files are saved in the `ewas_results` folder within the project folder.

So an example of this line, having set the arguments and output files would be:

```
R CMD BATCH --no-save --no-restore '--args continuous_bmi houseman Cells pace_mat_bmi_pheno B F7 mat_age,mat_ses,mat_smoke,parity /panfs/panasas01/sscm/gs8094/EWAS/' /panfs/panasas01/sscm/gs8094/Common_files/meffil_EWAS_script.r /panfs/panasas01/sscm/gs8094/EWAS/pace_mat_bmi/ewas_results/continuous_cell_adj.out
```

And that's the .sh file all sorted :)

Make sure you have an .sh file for each EWAS model you want to run. For example, your next .sh file might be exactly the same as the one described above, but now you don't want to adjust for cell counts. In that case, I'd change the `CellAdj` argument to `noCells` and the output file to `continuous_cell_unadj.out`. I'd then call the new .sh file `continuous_cell_unadj.sh`.

Every .sh file should be saved in the submission_scripts folder.

Then to run the EWAS, you just submit the .sh file to Bluecrystal using:

```
qsub EWAS/project_name/submission_scripts/example.sh
```

# The EWAS R script
The actual EWAS all happens in R, but this is automated by the .sh script, so I shouldn't actually need to open R myself (unless my EWAS doesn't run properly and I want to troubleshoot why that is).

The R script is called `meffil_EWAS_script.r` and is stored in the Common_files folder.

Here's a walkthrough of what it's doing:

First of all, I want to pull in the arguments that I set in my .sh file, so the top of the R script looks like this:

```
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
```

Then I want to check that they've been set properly:

```
print(Trait)
print(CellData)
print(CellAdj)
print(Phenofile)
print(BorM)
print(TP)
print(Covariates)
print(WD)
```

Next I set the working directory (i.e. the EWAS project folder):

```
#set working directory
setwd(WD)
```

Load necessary R packages:

```
#load packages
library(foreign) #to read stata file
library(meffil) #to run EWAS
```

The description of all the ARIES samples is stored in an R object in the IEU directory on bluecrystal. I load that using:

```
#Load description of samples
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/samplesheet/data.Robj")
```

I don't want ALL the ARIES samples, I just want the samples for one time point (if I want to look at multiple time points, I either run multiple EWAS looking at each time point cross-sectionally, or I would have to write a more advanced script if I wanted to do something fancy like longitudinal modelling).

This line just extracts the info on the samples specific to the time point of interest:

```
samplesheet<-subset(samplesheet, time_point==TP)
```

If we're interested in any of the 3 time points in childhood, we'll want to make sure we have independent participants, so it makes sense to get rid of any twins. Since this is all wrapped up in `if()`, you don't need to worry about changing this if you're interested in mothers - mothers of twins won't be excluded.

```
if(TP !="antenatal" & TP !="FOM"){
qletB<-samplesheet$ALN[which(samplesheet$QLET=="B")] #find alns for multiple pregnancies
samplesheet<-samplesheet[-which(samplesheet$ALN %in% qletB),] #remove multiple pregnancies
}
```

The ARIES methylation data is also stored on bluecrystal. It takes a little while to load:

```
#load the methylation data
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj")
```

Next I get rid of any samples in the ARIES methylation data that correspond to time points I'm not interested in:

```
meth <- norm.beta.random[,samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
```

Just a bit of housekeeping (I changed norm.beta.random to meth because it's shorter and neater :) so we can get rid of norm.beta.random now):

```
rm(norm.beta.random)
```

Next, I load the "detection P-values", which give an idea of how well each probe and sample has performed. Then I select just those p-values that correspond to the time point I'm interested in.

```
#load detection P-values (used to filter all probes with a high detection P-value)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/detection_p_values/data.Robj")
pvals <- detp[,samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(detp)
```

Next, I load the annotation data, which for ARIES is for the 450k array. This annotation data includes information like mapped gene, genomic location and relation to CpG island for each CpG on the 450k array.

```
#load annotation data
annotation <- meffil.get.features("450k")
```

Now I want to filter out certain probes from the methylation data.

First I'll get rid of anything with a detection P-value over 0.05 for over 5% of samples (this threshold is up to you, 0.05 is commonly used, but there's some evidence that 1E-10 might be more appropriate), then remove anything on the X and Y chromosomes (if you don't want to do this, just delete these lines), then remove SNPs and control probes included on the array for quality control purposes. The final few lines here just let me know how many probes I have left and how many were removed for each reason:

```
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
```

Next, I load my phenotype data, which is saved as a .dta file in the EWAS project folder. As I mentioned above somewhere, it needs to be saved as an old (version 12) STATA file if you use this approach. Alternatively, you could just save it as something like a .csv file and then update this line to `Pheno <- read.csv(paste0(Phenofile,".csv"),stringsAsFactors=FALSE)`.

```
Load phenotype data (this should be stored in your working directory)
Pheno<-read.dta(paste0(Phenofile,".dta"))
```

Next I load my cell count data, i.e. the estimated cell proportions for each sample. These will be different depending on the reference dataset and time point used, so these lines just make sure we're reading and selecting the right cell counts (based on the CellData and TP arguments).

```
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
```

In the DNA methylation data, each sample has a "sentrix ID", which describes where the sample appeared on the Illumina chip (row and column number). In our ALSPAC phenotype data, individuals have an 'aln' ID. So I have to add the sentrix ID (called 'Sample_Name') to my phenotype data so that I can then match my phenotype data to the methylation data.

```
#Add Sample_Name to Pheno (assuming Pheno contains aln)
Pheno<-merge(Pheno,samplesheet[,c("ALN","Sample_Name")],by.x="aln",by.y="ALN")
```

Next, I want to get rid of any phenotype variables I don't need and get rid of any individuals with NAs (missing data). If I wanted to impute missing data (instead of doing a complete case analysis) I'd have to use a fancy imputation approach (talk to Harriet Mills and/or Kate Tilling).

```
#Prepare phenotype data
Covs<-strsplit(Covariates,split=",")[[1]]
Pheno<-na.omit(Pheno[,c("Sample_Name",Trait,Covs)])
```

Now I add estimated cell counts to my phenotype data:

```
#Merge Pheno with cell-counts
colnames(cells)[1]<-"Sample_Name"
Pheno<-merge(Pheno,cells,by.x="Sample_Name",by.y="Sample_Name")
```

And then match my methylation data to my phenotype data. This just means I get the individuals in my phenotype data in the same order as the individuals in the methylation data, so the samples match up. The last line prints a message to let me know whether the matching was successful.

```
#Match meth to Pheno
meth<-meth[,na.omit(match(Pheno$Sample_Name,colnames(meth)))]
Pheno<-Pheno[match(colnames(meth),Pheno$Sample_Name),]
ifelse(all(Pheno$Sample_Name==colnames(meth)), "meth and phenotype data successfully matched :) ","Data not matched :(")
```

Next, I get rid of any "empty" factor levels. This just means that, if any of my phenotype variables were being treated by R as factor variables, I would be getting rid of any levels that didn't apply to anyone in my dataset. E.g. if everyone in my dataset had either 1 or 0 previous children, but the variable for number of previous children included levels 0, 1 and 2 I would be dropping 2 as an empty level. Don't worry too much if you don't understand this. It's just an easy step that can stop later functions from failing.

```
Pheno<-droplevels(Pheno) # get rid of any empty factor levels 
```

At this point, I like to give myself a little summary of the data to make sure I've roughly got the number of people I was expecting 

```
#Little summary of the phenotype data at this point
paste("There are ",nrow(Pheno)," people in this analysis")
"Here's a summary of the phenotype data:"
summary(Pheno)
```

If I want to convert beta values to M-values (set using the BorM argument), then this bit will do that:

```
# Convert to M-values if necessary
if(BorM=="M"){
	meth <- log2(meth/(1-meth))
}
```

If I don't want to include cell counts in the EWAS model (set using the CellAdj argument), this bit will get rid of cell counts from my phenotype data:

```
#Include cell counts in the EWAS model?
if(CellAdj=="noCells"){
Pheno<-Pheno[,-which(colnames(Pheno) %in% colnames(cells)[-1])]
}
```

Now I've finally finished setting up my data for the EWAS :)
It's time to actually run the EWAS using the meffil package.
This is achieved using a single function called meffil.ewas. I feed it the following arguments:

* methylation data matrix (meth)
* my trait of interest [i.e. the second column of Pheno] (variable=Pheno[,2])
* my covariates [i.e. everything in Pheno that isn't the ID or my trait of interest] (covariates=Pheno[,-(1:2)]) 
* do I want to remove outliers by winsorizing? No (winsorize.pct = NA)
* how many CpGs do we want to base the surrogate variable analysis on? All the CpGs, or if there are >20,000, the 20,000 most variable [most.variable = min(nrow(meth), 20000)] 
* do I want to remove outliers using the Tukey Method? Yes, and my threshod is the IQR multiplied by 3 (outlier.iqr.factor=3)
* do I want to know how meffil is getting on with my EWAS? Yes [verbose=TRUE]

```
#Run EWAS using meffil
obj <- meffil.ewas(meth, variable=Pheno[,2], covariates=Pheno[,-(1:2)], winsorize.pct = NA  ,most.variable = min(nrow(meth), 20000), outlier.iqr.factor=3, verbose=TRUE)
```

This should take a while to run, because it's basically running four EWAS at once (none, all, sva, isva), each one involving 450,000 regression analyses on hundreds of people.

From the EWAS output, I now select just the information I'm interested in, which is probe ID plus the coefficients, standard errors and p-values for each of the four EWAS. I save this in a dataframe called 'ewas_res'.

```
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
```

This next section makes a HTML report summarizing the results of the EWAS. Here I select just to look at the isva model, but you could change this. The HTML fine is saved in the `ewas_results` folder:

```
ewas.parameters<-meffil.ewas.parameters(sig.threshold=1e-5, max.plots=5, model="isva")
ewas.summary <- meffil.ewas.summary(ewas.results, meth, parameters=ewas.parameters)   
savefile <- paste("ewas_results/",Trait,TP,Covariates,CellAdj,Sys.Date(),".html", sep = "_")
meffil.ewas.report(ewas.summary, output.file=savefile,author="gemma sharp", study="alspac")
```

Now I add a bit more information to EWAS results dataframe. First I create a column describing how many samples were included in the analysis. This will vary by probe because I removed outliers using the IQR method. `original_n` is the number of samples before removal of outliers. `n_outliers` and `final_n` are self-explanatory :) 

```
#Create N_for_probe column
ewas_res$original_n<-rowSums(!is.na(meth))
outliers_n<-data.frame(table(c(rownames(obj$too.hi),rownames(obj$too.lo))))
colnames(outliers_n)<-c("probeID","n_outliers")
ewas_res<-merge(ewas_res,outliers_n,by="probeID",all.x=TRUE)
ewas_res$n_outliers<-ewas_res$n_outliers*-1
ewas_res$final_n<-rowSums(ewas_res[,c("original_n","n_outliers")],na.rm=TRUE)
```

If my trait of interest is binary (e.g. case/control), I will also want to know how many cases were included in the EWAS for each probe. This section firstly identifies if the trait of interest is binary and then (if it is) it will add more columns to the EWAS results dataframe to described the number of outliers removed within the cases and the final number of cases and controls.

```
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
```

Next, I add a column indicating whether probes are on the Naeem list of possibly problematic probes or not

```
#Load the Naeem list of problematic probes
Naeem<-read.csv("/panfs/panasas01/sscm/gs8094/Common_files/naeem_list.csv")
ewas_res$OnNaeem<-ifelse(ewas_res$probeID %in% Naeem$EXCLUDE_PROBES,"yes","no")
```

I adjust P-values from each EWAS (none, all, sva, isva) for multiple testing using the FDR and the Bonferroni approaches

```
# Adjustment for multiple testing
ewas_res$fdr.none<-p.adjust(ewas_res$p.none, method="fdr") 
ewas_res$bonferroni.none<-p.adjust(ewas_res$p.none, method="bonferroni")
ewas_res$fdr.all<-p.adjust(ewas_res$p.all, method="fdr") 
ewas_res$bonferroni.all<-p.adjust(ewas_res$p.all, method="bonferroni")
ewas_res$fdr.sva<-p.adjust(ewas_res$p.sva, method="fdr") 
ewas_res$bonferroni.sva<-p.adjust(ewas_res$p.sva, method="bonferroni")
ewas_res$fdr.isva<-p.adjust(ewas_res$p.isva, method="fdr") 
ewas_res$bonferroni.isva<-p.adjust(ewas_res$p.isva, method="bonferroni")
```

Next I add more information about my EWAS design. This won't vary by probe, but it is useful to have this information in the dataframe in case I forget what the inputs were (very possible).

```
##Annotate and sort
#Add information about EWAS
ewas_res$Trait<-Trait
ewas_res$Covariates<-paste0(ls(obj$covariates),collapse=", ")
ewas_res$nSVs<-ncol(obj$analyses$sva$design) #Number 
ewas_res$nISVs<-ncol(obj$analyses$isva$design)
ewas_res$TP<-TP
ewas_res$BorM<-BorM
ewas_res$CellData<-CellData
ewas_res$CellAdj<-CellAdj
```

I add the lambda (a measure of test statitic inflation), which again, won't vary by probe, but it's useful to have this information available quickly. Interpretation of lambda in EWAS is difficult and there is no consensus on its use, but if I get a very high Lambda (e.g. >1.5) I would worry about residual confounding and think carefully about my models... But this really depends on the phenotype of interest and its proposed biological effect, so don't worry about it too much.

```
#Add Lambda (a measure of test statistic inflation)
Lambda<-function(P){
chisq <- qchisq(1-P,1)
median(chisq,na.rm=T)/qchisq(0.5,1)
}

ewas_res$lambda.none<-Lambda(ewas_res$p.none)
ewas_res$lambda.all<-Lambda(ewas_res$p.all)
ewas_res$lambda.sva<-Lambda(ewas_res$p.sva)
ewas_res$lambda.isva<-Lambda(ewas_res$p.isva)
```

Next I add information from Illumina about each probe. This annotation data includes information like mapped gene, genomic location and relation to CpG island for each CpG on the 450k array.

```
#Add annotation information about probes
ewas_res<-merge(ewas_res,annotation,by.x="probeID", by.y="name",all.x=TRUE)
```

Finally, I sort by P-value (from the isva model) and save the EWAS results as a .Rdata file in the 'ewas_res' folder:

```
ewas_res<-ewas_res[order(ewas_res$p.isva),]
#Save as an Rdata file
savefile <- paste("ewas_results/",Trait,TP,Covariates,CellAdj,Sys.Date(),".Rdata", sep = "_")
save(ewas_res,file=savefile)
```

